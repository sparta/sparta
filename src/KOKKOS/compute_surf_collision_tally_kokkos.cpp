/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "compute_surf_collision_tally_kokkos.h"
#include <type_traits>
#include "particle_kokkos.h"
#include "surf_kokkos.h"
#include "update.h"
#include "domain.h"
#include "mixture.h"
#include "memory_kokkos.h"
#include "sparta_masks.h"
#include "error.h"

using namespace SPARTA_NS;

#define DELTA 4096

/* ---------------------------------------------------------------------- */

ComputeSurfCollisionTallyKokkos::ComputeSurfCollisionTallyKokkos(SPARTA *sparta,
                                                                int narg, char **arg) :
  ComputeSurfCollisionTally(sparta, narg, arg)
{
  kokkos_flag = 1;
  ntally_mark = 0;

  d_ntally = DAT::t_int_scalar("surf/collision/tally/kk:ntally");
  h_ntally = HAT::t_int_scalar("surf/collision/tally/kk:ntally_mirror");

  // flatten the value list once; it never changes after construction

  DAT::tdual_int_1d k_which("surf/collision/tally/kk:which",nvalue);
  for (int m = 0; m < nvalue; m++) k_which.view_host()[m] = which[m];
  k_which.modify_host();
  k_which.sync_device();
  d_which = k_which.view_device();

  maxtally = DELTA;
  maxtally_host = 0;
  MemKK::realloc_kokkos(k_array_tally,"surf/collision/tally/kk:array_tally",
                        maxtally,nvalue);
  d_array_tally = k_array_tally.view_device();
}

/* ---------------------------------------------------------------------- */

ComputeSurfCollisionTallyKokkos::ComputeSurfCollisionTallyKokkos(SPARTA *sparta) :
  ComputeSurfCollisionTally(sparta)
{
  copy = 1;
  ntally_mark = 0;
}

/* ---------------------------------------------------------------------- */

ComputeSurfCollisionTallyKokkos::~ComputeSurfCollisionTallyKokkos()
{
  if (copy || copymode) return;

  // the host base class frees array_tally, which this class never allocated

  memory->destroy(array_tally);
  array_tally = NULL;
  maxtally = 0;
}

/* ---------------------------------------------------------------------- */

void ComputeSurfCollisionTallyKokkos::clear()
{
  ntally = 0;
}

/* ----------------------------------------------------------------------
   called by UpdateKokkos before the move kernel
------------------------------------------------------------------------- */

void ComputeSurfCollisionTallyKokkos::pre_surf_tally()
{
  SurfKokkos* surf_kk = (SurfKokkos*) surf;
  surf_kk->sync(Device,LINE_MASK|TRI_MASK);
  d_lines = surf_kk->k_lines.view_device();
  d_tris = surf_kk->k_tris.view_device();

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,SPECIES_MASK);
  d_s2g = particle_kk->k_species2group.view_device();

  dim = domain->dimension;
  dt = update->dt;

  Kokkos::deep_copy(d_ntally,0);
}

/* ----------------------------------------------------------------------
   called by UpdateKokkos after the move kernel
   bring the row count and the rows themselves to the host, where
     dump tally and Compute::tallyinfo() read them
------------------------------------------------------------------------- */

void ComputeSurfCollisionTallyKokkos::post_surf_tally()
{
  Kokkos::deep_copy(h_ntally,d_ntally);
  ntally = h_ntally();

  // an overflowed attempt is discarded and repeated by UpdateKokkos, so do
  //   not publish its partial rows.  the device counter kept climbing past
  //   the end of the buffer, so ntally is not a row count here -- leaving it
  //   in the host base class would make dump tally read that many rows out
  //   of an array that never held them

  if (ntally > (int) d_array_tally.extent(0)) { ntally = 0; return; }

  // the host base class hands out array_tally, row-indexed to ntally
  //   (compute_surf_collision_tally.cpp:172).  Publish only the rows this
  //   step produced: k_array_tally is sized to maxtally, a high-water mark,
  //   so sync_host() moved the whole buffer and the scalar loop then touched
  //   every element of it again.  The host half of the dual view has no
  //   other reader, so copy the device rows straight into array_tally --
  //   d_array_tally is LayoutRight and Memory::create hands back one
  //   contiguous row-major block with the same nvalue row stride, so an
  //   unmanaged wrap of the leading ntally rows lines up exactly.  Keep the
  //   allocation on its own high-water mark rather than free/malloc-ing it
  //   every step.

  if (ntally) {
    static_assert(std::is_same<F_FLOAT,double>::value,
                  "array_tally is double**, so F_FLOAT must be double "
                  "for the unmanaged wrap below to alias it");
    if (ntally > maxtally_host) {
      memory->destroy(array_tally);
      maxtally_host = MAX(ntally,maxtally);
      memory->create(array_tally,maxtally_host,nvalue,
                     "surf/collision/tally/kk:array_tally_host");
    }
    Kokkos::View<F_FLOAT**,Kokkos::LayoutRight,Kokkos::HostSpace,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      h_rows(array_tally[0],ntally,nvalue);
    Kokkos::deep_copy(h_rows,
                      Kokkos::subview(d_array_tally,
                                      Kokkos::make_pair(0,ntally),
                                      Kokkos::ALL()));
  }
}

/* ----------------------------------------------------------------------
   grow the device row buffer to hold at least N rows
   called by UpdateKokkos when an attempt overflowed, before repeating it
------------------------------------------------------------------------- */

void ComputeSurfCollisionTallyKokkos::grow_tally_kokkos(int n)
{
  if (n <= (int) d_array_tally.extent(0)) return;

  // grow past what the failed attempt asked for, so a step whose collision
  //   count is still climbing does not repeat the move again and again

  maxtally = MAX(n + DELTA, (int)(1.5*n));

  // resize, not realloc: MemKK::realloc_kokkos drops the old allocation and
  //   allocates WithoutInitializing, so every existing row becomes garbage.
  //   The rows below ntally_mark were written by earlier migration iterations
  //   of this same step and must survive -- rewind_ntally() only takes the
  //   count back to the mark, and the re-run appends from there, so a
  //   discarded prefix would be published as uninitialized tally output.
  //   DualView::resize preserves contents, as memoryKK->grow_kokkos does.
  //
  //   No deck reaches that combination today: maxtally is a high-water mark
  //   that persists across steps, so by the time a step has more than one
  //   tallying migration iteration the buffer already fits the whole step and
  //   never overflows again.  Probing in.surf.collision.tally at 4 ranks with
  //   DELTA forced to 8 and then to no headroom at all gave 13 grows, every
  //   one with ntally_mark == 0, against 40 nonzero marks elsewhere in the
  //   same run.  This is a latent hazard closed on principle, not a
  //   reproduced failure.

  k_array_tally.resize(maxtally,nvalue);
  d_array_tally = k_array_tally.view_device();
}

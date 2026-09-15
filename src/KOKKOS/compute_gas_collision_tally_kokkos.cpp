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

#include "compute_gas_collision_tally_kokkos.h"
#include <type_traits>
#include "particle_kokkos.h"
#include "grid_kokkos.h"
#include "update.h"
#include "domain.h"
#include "mixture.h"
#include "memory_kokkos.h"
#include "sparta_masks.h"
#include "error.h"

using namespace SPARTA_NS;

#define DELTA 4096

/* ---------------------------------------------------------------------- */

ComputeGasCollisionTallyKokkos::ComputeGasCollisionTallyKokkos(SPARTA *sparta,
                                                                int narg, char **arg) :
  ComputeGasCollisionTally(sparta, narg, arg)
{
  kokkos_flag = 1;
  ntally_mark = 0;

  d_ntally = DAT::t_int_scalar("gas/collision/tally/kk:ntally");
  h_ntally = HAT::t_int_scalar("gas/collision/tally/kk:ntally_mirror");

  // flatten the value list once; it never changes after construction

  DAT::tdual_int_1d k_which("gas/collision/tally/kk:which",nvalue);
  for (int m = 0; m < nvalue; m++) k_which.view_host()[m] = which[m];
  k_which.modify_host();
  k_which.sync_device();
  d_which = k_which.view_device();

  maxtally = DELTA;
  maxtally_host = 0;
  MemKK::realloc_kokkos(k_array_tally,"gas/collision/tally/kk:array_tally",
                        maxtally,nvalue);
  d_array_tally = k_array_tally.view_device();
}

/* ---------------------------------------------------------------------- */

ComputeGasCollisionTallyKokkos::ComputeGasCollisionTallyKokkos(SPARTA *sparta) :
  ComputeGasCollisionTally(sparta)
{
  copy = 1;
  ntally_mark = 0;
}

/* ---------------------------------------------------------------------- */

ComputeGasCollisionTallyKokkos::~ComputeGasCollisionTallyKokkos()
{
  if (copy || copymode) return;

  // the host base class frees array_tally, which this class never allocated

  memory->destroy(array_tally);
  array_tally = NULL;
  maxtally = 0;
}

/* ---------------------------------------------------------------------- */

void ComputeGasCollisionTallyKokkos::clear()
{
  ntally = 0;
}

/* ----------------------------------------------------------------------
   called by CollideVSSKokkos before the collision kernel
------------------------------------------------------------------------- */

void ComputeGasCollisionTallyKokkos::pre_gas_tally()
{
  GridKokkos* grid_kk = (GridKokkos*) grid;
  grid_kk->sync(Device,CELL_MASK|CINFO_MASK);
  d_cells = grid_kk->k_cells.view_device();
  d_cinfo = grid_kk->k_cinfo.view_device();

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,SPECIES_MASK);
  d_s2g = particle_kk->k_species2group.view_device();

  Kokkos::deep_copy(d_ntally,0);
}

/* ----------------------------------------------------------------------
   called by CollideVSSKokkos after the collision kernel
   bring the row count and the rows themselves to the host, where
     dump tally and Compute::tallyinfo() read them
------------------------------------------------------------------------- */

void ComputeGasCollisionTallyKokkos::post_gas_tally()
{
  Kokkos::deep_copy(h_ntally,d_ntally);
  ntally = h_ntally();

  // an overflowed attempt is discarded and repeated by CollideVSSKokkos, so do
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
                     "gas/collision/tally/kk:array_tally_host");
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
   called by CollideVSSKokkos when an attempt overflowed, before repeating it
------------------------------------------------------------------------- */

void ComputeGasCollisionTallyKokkos::grow_tally_kokkos(int n)
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
  //   discarded prefix is published as uninitialized tally output.
  //   DualView::resize preserves contents, as memoryKK->grow_kokkos does.

  k_array_tally.resize(maxtally,nvalue);
  d_array_tally = k_array_tally.view_device();
}

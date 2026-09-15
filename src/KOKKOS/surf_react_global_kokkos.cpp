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

#include "math.h"
#include "surf_react_global_kokkos.h"
#include "input.h"
#include "update.h"
#include "comm.h"
#include "error.h"
#include "random_knuth.h"
#include "particle_kokkos.h"
#include "sparta_masks.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

SurfReactGlobalKokkos::SurfReactGlobalKokkos(SPARTA *sparta, int narg, char **arg) :
  SurfReactGlobal(sparta,narg,arg),
  rand_pool(12345 + comm->me
#ifdef SPARTA_KOKKOS_EXACT
            , sparta
#endif
            )
{
  kokkosable = 1;

  k_nsingle = DAT::tdual_int_scalar("surf_react:nsingle");
  d_nsingle = k_nsingle.view_device();
  h_nsingle = k_nsingle.view_host();

  k_tally_single = DAT::tdual_bigint_1d("surf_react:tally_single",nlist);
  d_tally_single = k_tally_single.view_device();
  h_tally_single = k_tally_single.view_host();

  random_backup = NULL;

#ifdef SPARTA_KOKKOS_EXACT
  // allocate on the real class instance: backup() is only ever invoked on a
  //  KKCopy of this class, and a pointer stored there is dropped by the next
  //  copy() of the original over it

  random_backup = new RanKnuth(12345 + comm->me);
#endif
}

SurfReactGlobalKokkos::SurfReactGlobalKokkos(SPARTA *sparta) :
  SurfReactGlobal(sparta),
  rand_pool(12345 // seed will be copied over
#ifdef SPARTA_KOKKOS_EXACT
            , sparta
#endif
            )
{
  copy = 1;
}

/* ---------------------------------------------------------------------- */

SurfReactGlobalKokkos::~SurfReactGlobalKokkos()
{
  if (copy) return;

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.destroy();
  if (random_backup)
    delete random_backup;
#endif
}

/* ---------------------------------------------------------------------- */

void SurfReactGlobalKokkos::init()
{
  SurfReact::init();

  Kokkos::deep_copy(d_nsingle,0);
  Kokkos::deep_copy(d_tally_single,0);

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.init(random);
#endif
}

/* ---------------------------------------------------------------------- */

void SurfReactGlobalKokkos::tally_reset()
{
  SurfReact::tally_reset();

  Kokkos::deep_copy(d_nsingle,0);
  Kokkos::deep_copy(d_tally_single,0);
}

/* ---------------------------------------------------------------------- */

void SurfReactGlobalKokkos::tally_update()
{
  Kokkos::deep_copy(h_nsingle,d_nsingle);
  Kokkos::deep_copy(h_tally_single,d_tally_single);

  // also store the per-step counts on the host: compute_vector() reports
  //   nsingle and tally_single_all alongside the running totals

  nsingle = h_nsingle();
  ntotal += h_nsingle();
  for (int i = 0; i < nlist; i++) {
    tally_single[i] = h_tally_single[i];
    tally_total[i] += h_tally_single[i];
  }
}

/* ---------------------------------------------------------------------- */

void SurfReactGlobalKokkos::pre_react()
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,PARTICLE_MASK);
  d_particles = particle_kk->k_particles.view_device();
}

/* ---------------------------------------------------------------------- */

void SurfReactGlobalKokkos::post_react()
{
  d_particles = {};
}

/* ---------------------------------------------------------------------- */

void SurfReactGlobalKokkos::backup()
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  d_particles = particle_kk->k_particles.view_device();

  if (!d_nsingle_backup.data())
    d_nsingle_backup = DAT::t_int_scalar(
      Kokkos::view_alloc("surf_react:nsingle_backup",Kokkos::WithoutInitializing));
  Kokkos::deep_copy(d_nsingle_backup,d_nsingle);

  if (d_tally_single_backup.extent(0) != d_tally_single.extent(0))
    d_tally_single_backup = DAT::t_bigint_1d(
      Kokkos::view_alloc("surf_react:tally_single_backup",Kokkos::WithoutInitializing),
      d_tally_single.extent(0));
  Kokkos::deep_copy(d_tally_single_backup,d_tally_single);

#ifdef SPARTA_KOKKOS_EXACT
  if (!random_backup)
    random_backup = new RanKnuth(12345 + comm->me);
  memcpy(random_backup,random,sizeof(RanKnuth));
#endif
}

/* ---------------------------------------------------------------------- */

void SurfReactGlobalKokkos::restore()
{
  Kokkos::deep_copy(d_nsingle,d_nsingle_backup);
  Kokkos::deep_copy(d_tally_single,d_tally_single_backup);

#ifdef SPARTA_KOKKOS_EXACT
  memcpy(random,random_backup,sizeof(RanKnuth));
#endif
}

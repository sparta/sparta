/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
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

  d_scalars = t_int_3("surf_react:scalars");
  d_nsingle = Kokkos::subview(d_scalars,0);
  d_tally_single = Kokkos::subview(d_scalars,std::make_pair(1,3));

  h_scalars = t_host_int_3("surf_react:scalars_mirror");
  h_nsingle = Kokkos::subview(h_scalars,0);
  h_tally_single = Kokkos::subview(h_scalars,std::make_pair(1,3));

  random_backup = NULL;
}

SurfReactGlobalKokkos::SurfReactGlobalKokkos(SPARTA *sparta) :
  SurfReactGlobal(sparta),
  rand_pool(12345 // seed will be copied over
#ifdef SPARTA_KOKKOS_EXACT
            , sparta
#endif
            )
{
  random = NULL;
  random_backup = NULL;

  id = NULL;
  style = NULL;
  tally_single = NULL;
  tally_total = NULL;
  tally_single_all = NULL;
  tally_total_all = NULL;
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

  Kokkos::deep_copy(d_scalars,0);

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.init(random);
#endif
}

/* ---------------------------------------------------------------------- */

void SurfReactGlobalKokkos::tally_reset()
{
  SurfReact::tally_reset();

  Kokkos::deep_copy(d_scalars,0);
}

/* ---------------------------------------------------------------------- */

void SurfReactGlobalKokkos::tally_update()
{
  Kokkos::deep_copy(h_scalars,d_scalars);
  ntotal += h_nsingle();
  for (int i = 0; i < nlist; i++) tally_total[i] += h_tally_single[i];
}

/* ---------------------------------------------------------------------- */

void SurfReactGlobalKokkos::pre_react()
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,PARTICLE_MASK);
  d_particles = particle_kk->k_particles.d_view;
}

/* ---------------------------------------------------------------------- */

void SurfReactGlobalKokkos::backup()
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  d_particles = particle_kk->k_particles.d_view;

#ifdef SPARTA_KOKKOS_EXACT
  if (!random_backup)
    random_backup = new RanKnuth(12345 + comm->me);
  memcpy(random_backup,random,sizeof(RanKnuth));
#endif
}

/* ---------------------------------------------------------------------- */

void SurfReactGlobalKokkos::restore()
{
  Kokkos::deep_copy(d_scalars,0);

#ifdef SPARTA_KOKKOS_EXACT
  memcpy(random,random_backup,sizeof(RanKnuth));
#endif
}

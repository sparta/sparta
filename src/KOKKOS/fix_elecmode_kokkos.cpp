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
#include "fix_elecmode_kokkos.h"
#include "update.h"
#include "particle.h"
#include "collide.h"
#include "comm.h"
#include "memory_kokkos.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "math_const.h"
#include "error.h"
#include "sparta_masks.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{INT,DOUBLE};                      // several files
enum{NONE,DISCRETE,SMOOTH};            // several files

/* ---------------------------------------------------------------------- */

FixElecmodeKokkos::FixElecmodeKokkos(SPARTA *sparta, int narg, char **arg) :
  FixElecmode(sparta, narg, arg),
  rand_pool(12345 + comm->me
#ifdef SPARTA_KOKKOS_EXACT
            , sparta
#endif
            )
{
  kokkos_flag = 1;
  execution_space = Device;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.init(random);
#endif
}

/* ---------------------------------------------------------------------- */

FixElecmodeKokkos::FixElecmodeKokkos(SPARTA *sparta) :
  FixElecmode(sparta),
  rand_pool(12345 // seed doesn't matter since it will just be copied over
#ifdef SPARTA_KOKKOS_EXACT
            , sparta
#endif
            )
{
  random = NULL;
  id = NULL;
  style = NULL;
}

/* ---------------------------------------------------------------------- */

FixElecmodeKokkos::~FixElecmodeKokkos()
{
  if (copy) return;

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.destroy();
#endif
}

/* ---------------------------------------------------------------------- */

void FixElecmodeKokkos::pre_update_custom_kokkos()
{
  boltz = update->boltz;

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,PARTICLE_MASK|SPECIES_MASK|CUSTOM_MASK);
  d_particles = particle_kk->k_particles.d_view;
  d_species = particle_kk->k_species.d_view;
  d_nelecstates = particle_kk->d_nelecstates;
  auto h_ewhich = particle_kk->k_ewhich.h_view;
  auto k_edvec = particle_kk->k_edvec;
  auto k_eivec = particle_kk->k_eivec;
  d_eelec = k_edvec.h_view[h_ewhich[index_eelec]].k_view.d_view;
  d_elecstate = k_eivec.h_view[h_ewhich[index_elecstate]].k_view.d_view;
  if (particle->maxlocal > (int)d_cumulative_probabilities.extent(0))
    MemKK::realloc_kokkos(d_cumulative_probabilities,"cumulative_probabilities",particle->maxlocal,particle->maxelecstate);

  elecstyle = NONE;
  if (collide) elecstyle = collide->elecstyle;
}

/* ----------------------------------------------------------------------
   called when a particle with index is created
    or when temperature dependent properties need to be updated
   populate an electronic state and set eelec
------------------------------------------------------------------------- */

void FixElecmodeKokkos::update_custom(int index, double temp_thermal,
                                     double temp_rot, double temp_vib, double temp_elec,
                                     double *vstream)
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Host,PARTICLE_MASK|SPECIES_MASK|CUSTOM_MASK);
  FixElecmode::update_custom(index, temp_thermal, temp_rot, temp_vib, temp_elec, vstream);
  particle_kk->modify(Host,PARTICLE_MASK|CUSTOM_MASK);
}

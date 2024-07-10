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
#include "stdlib.h"
#include "string.h"
#include "fix_ambipolar_kokkos.h"
#include "update.h"
#include "particle.h"
#include "comm.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "error.h"
#include "sparta_masks.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{INT,DOUBLE};                      // several files

/* ---------------------------------------------------------------------- */

FixAmbipolarKokkos::FixAmbipolarKokkos(SPARTA *sparta, int narg, char **arg) :
  FixAmbipolar(sparta, narg, arg),
  rand_pool(12345 + comm->me
#ifdef SPARTA_KOKKOS_EXACT
            , sparta
#endif
            )
{
  kokkos_flag = 1;

  // random = RNG for electron velocity creation

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.init(random);
#endif

  d_ions = DAT::t_int_1d("fix_ambipolar:ions",maxion);
  auto h_ions = Kokkos::create_mirror_view(d_ions);
  for (int i = 0; i < maxion; i++) h_ions[i] = ions[i];
  Kokkos::deep_copy(d_ions,h_ions);
}

/* ---------------------------------------------------------------------- */

FixAmbipolarKokkos::FixAmbipolarKokkos(SPARTA *sparta) :
  FixAmbipolar(sparta),
  rand_pool(12345 // seed doesn't matter since it will just be copied over
#ifdef SPARTA_KOKKOS_EXACT
            , sparta
#endif
            )
{
  ions = NULL;
  random = NULL;
  id = NULL;
  style = NULL;
}

/* ---------------------------------------------------------------------- */

FixAmbipolarKokkos::~FixAmbipolarKokkos()
{
  if (copy) return;

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.destroy();
#endif
}

/* ---------------------------------------------------------------------- */

void FixAmbipolarKokkos::pre_update_custom_kokkos()
{
  boltz = update->boltz;
  temp_thermal = update->temp_thermal;
  vstream[0] = update->vstream[0];
  vstream[1] = update->vstream[1];
  vstream[2] = update->vstream[2];

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,PARTICLE_MASK|SPECIES_MASK|CUSTOM_MASK);
  d_particles = particle_kk->k_particles.d_view;
  d_species = particle_kk->k_species.d_view;
  auto h_ewhich = particle_kk->k_ewhich.h_view;
  auto k_eivec = particle_kk->k_eivec;
  auto k_edarray = particle_kk->k_edarray;
  d_ionambi = k_eivec.h_view[h_ewhich[ionindex]].k_view.d_view;
  d_velambi = k_edarray.h_view[h_ewhich[velindex]].k_view.d_view;
}

/* ----------------------------------------------------------------------
   called when a particle with index is created
   creation used temp_thermal and vstream to set particle velocity
   if an ion, set ionambi and velambi for particle
------------------------------------------------------------------------- */

void FixAmbipolarKokkos::update_custom(int index, double temp_thermal,
                                       double temp_rot, double temp_vib,
                                       double *vstream)
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Host,PARTICLE_MASK|SPECIES_MASK|CUSTOM_MASK);
  FixAmbipolar::update_custom(index, temp_thermal, temp_rot, temp_vib, vstream);
  particle_kk->modify(Host,CUSTOM_MASK);
}

/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
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

FixAmbipolarKokkos::~FixAmbipolarKokkos()
{
  if (copymode) return;

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.destroy();
#endif
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

/* ----------------------------------------------------------------------
   called when a surface reaction occurs
   iorig = particle I before reaction
   I,J = indices of two particles after reaction
         either can be -1, meaning particle does not exist
------------------------------------------------------------------------- */

void FixAmbipolarKokkos::surf_react(Particle::OnePart *iorig, int &i, int &j)
{
  // not yet supported by Kokkos package
}

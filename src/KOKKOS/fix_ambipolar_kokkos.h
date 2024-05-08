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

#ifdef FIX_CLASS

FixStyle(ambipolar/kk,FixAmbipolarKokkos)

#else

#ifndef SPARTA_FIX_AMBIPOLAR_KOKKOS_H
#define SPARTA_FIX_AMBIPOLAR_KOKKOS_H

#include "fix_ambipolar.h"
#include "particle_kokkos.h"
#include "kokkos_type.h"
#include "math_const.h"

namespace SPARTA_NS {

class FixAmbipolarKokkos : public FixAmbipolar {
 public:
  DAT::t_int_1d d_ions;                  // 1 if a particle species is an ionx

  FixAmbipolarKokkos(class SPARTA *, int, char **);
  FixAmbipolarKokkos(class SPARTA *);
  ~FixAmbipolarKokkos();
  void pre_update_custom_kokkos();
  void update_custom(int, double, double, double, double *);

 private:

  double boltz;
  double temp_thermal;
  double vstream[3];

  t_particle_1d d_particles;
  t_species_1d d_species;

  DAT::t_int_1d d_ionambi;
  DAT::t_float_2d_lr d_velambi;

#ifndef SPARTA_KOKKOS_EXACT
  Kokkos::Random_XorShift64_Pool<DeviceType> rand_pool;
  typedef typename Kokkos::Random_XorShift64_Pool<DeviceType>::generator_type rand_type;

  //Kokkos::Random_XorShift1024_Pool<DeviceType> rand_pool;
  //typedef typename Kokkos::Random_XorShift1024_Pool<DeviceType>::generator_type rand_type;
#else
  RandPoolWrap rand_pool;
  typedef RandWrap rand_type;
#endif

 public:

  /* ----------------------------------------------------------------------
     called when a particle with index is created
     creation used temp_thermal and vstream to set particle velocity
     if an ion, set ionambi and velambi for particle
  ------------------------------------------------------------------------- */

  KOKKOS_INLINE_FUNCTION
  void update_custom_kokkos(int index, double temp_thermal,
                            double, double,
                            const double *vstream) const
  {
    // if species is not ambipolar ion, set ionambi off and return

    int ispecies = d_particles[index].ispecies;

    if (d_ions[ispecies] == 0) {
      d_ionambi[index] = 0;
      return;
    }

    // set velocity of electron
    // based on electron mass, thermal temperature, and streaming velocity

    d_ionambi[index] = 1;

    const double vscale = sqrt(2.0 * boltz * temp_thermal /
                               d_species[especies].mass);

    rand_type rand_gen = rand_pool.get_state();

    const double vn = vscale * sqrt(-log(rand_gen.drand()));
    const double vr = vscale * sqrt(-log(rand_gen.drand()));
    const double theta1 = MathConst::MY_2PI * rand_gen.drand();
    const double theta2 = MathConst::MY_2PI * rand_gen.drand();

    d_velambi(index,0) = vstream[0] + vn*cos(theta1);
    d_velambi(index,1) = vstream[1] + vr*cos(theta2);
    d_velambi(index,2) = vstream[2] + vr*sin(theta2);

    rand_pool.free_state(rand_gen);
  }

  /* ----------------------------------------------------------------------
     called when a surface reaction occurs
     iorig = particle I before reaction
     I,J = indices of two particles after reaction
           either can be -1, meaning particle does not exist
  ------------------------------------------------------------------------- */

  KOKKOS_INLINE_FUNCTION
  void surf_react_kokkos(Particle::OnePart *iorig, int &i, int &j) const
  {
    int ispecies = iorig->ispecies;

    // recombination reaction, just return

    if (i < 0) return;

    // exchange reaction
    // if ion -> non-ion, unset ionambi flag

    if (j < 0) {
      if (d_ions[ispecies] == 0) return;
      if (d_ions[d_particles[i].ispecies] == 1) return;
      d_ionambi[i] = 0;
    }

    // dissociation reaction
    // if non-ion -> ion + electron, create an ambipolar ion
    // use global temp_thermal and vstream for electron creation
    // set j = -1 to delete electron that was just created by caller

    else {
      if (d_ions[ispecies] == 1) return;
      if (d_ions[d_particles[i].ispecies] == 0) return;
      if (d_particles[j].ispecies != especies) return;
      update_custom_kokkos(i,temp_thermal,temp_thermal,
                           temp_thermal,vstream);
      j = -1;
    }
  }

};

}

#endif
#endif

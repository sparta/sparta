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

FixStyle(elecmode/kk,FixElecmodeKokkos)

#else

#ifndef SPARTA_FIX_ELECMODE_KOKKOS_H
#define SPARTA_FIX_ELECMODE_KOKKOS_H

#include "fix_elecmode.h"
#include "kokkos_type.h"
#include "particle_kokkos.h"
#include "Kokkos_Random.hpp"
#include "rand_pool_wrap.h"

namespace SPARTA_NS {

class FixElecmodeKokkos : public FixElecmode {
 public:
  FixElecmodeKokkos(class SPARTA *, int, char **);
  FixElecmodeKokkos(class SPARTA *sparta);
  ~FixElecmodeKokkos();
  void pre_update_custom_kokkos();
  void update_custom(int, double, double, double, double, double *);

  KOKKOS_INLINE_FUNCTION
  void update_custom_kokkos(int, double, double, double, double, const double *) const;

 private:
  int boltz;

#ifndef SPARTA_KOKKOS_EXACT
  Kokkos::Random_XorShift64_Pool<DeviceType> rand_pool;
  typedef typename Kokkos::Random_XorShift64_Pool<DeviceType>::generator_type rand_type;

  //Kokkos::Random_XorShift1024_Pool<DeviceType> rand_pool;
  //typedef typename Kokkos::Random_XorShift1024_Pool<DeviceType>::generator_type rand_type;
#else
  RandPoolWrap rand_pool;
  typedef RandWrap rand_type;
#endif

  t_particle_1d d_particles;
  t_species_1d d_species;

  DAT::t_float_1d d_eelec;
  DAT::t_int_1d d_elecstate;
};

/* ----------------------------------------------------------------------
   called when a particle with index is created
    or when temperature dependent properties need to be updated
   populate an electronic state and set eelec
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixElecmodeKokkos::update_custom_kokkos(int index, double temp_thermal,
                                            double temp_rot, double temp_vib,
                                            double temp_elec, const double *) const
{
  int isp = d_particles[index].ispecies;
  int nstate = d_species[isp].nelecstate;

  // no states, just return

  if (nstate == 0) return;

  d_elecstate[index] = 0; // Need to update somehow or remove
  d_eelec[index] = particle->eelec(isp,temp_elec,random);
}

}

#endif
#endif

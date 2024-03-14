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
  int boltz,elecstyle;

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
  DAT::t_int_1d d_nelecstates;
  t_elecstate_2d d_elecstates;
  DAT::t_float_2d d_cumulative_probabilities;

  KOKKOS_INLINE_FUNCTION
  void electronic_distribution_func(int, int, double) const;

  KOKKOS_INLINE_FUNCTION
  double ielec(int, int, double, rand_type &) const;
};

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixElecmodeKokkos::electronic_distribution_func(int index, int isp, double temp_elec) const
{
  auto &d_distribution = d_cumulative_probabilities;
  double partition_function = 0.0;
  const int nelecstate = d_nelecstates[isp];

  for (int i = 0; i < nelecstate; ++i) {
    // Calculate boltzmann fractions
    d_distribution(index,i) = d_elecstates(isp,i).degen*exp(-d_elecstates(isp,i).temp/temp_elec);
    // Calculate partition function
    partition_function += d_distribution(index,i);
  }

  for (int i = 0; i < nelecstate; ++i)
    d_distribution(index,i) /= partition_function;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
double FixElecmodeKokkos::ielec(int index, int isp, double temp_elec, rand_type &erandom) const
{
  enum{NONE,DISCRETE,SMOOTH};            // several files

  int ielec = 0;

  if (elecstyle == DISCRETE) {
    int nelecstate = d_nelecstates[isp];
    if (!nelecstate) return 0.0;

    electronic_distribution_func(index, isp, temp_elec);

    for (int i = 1; i < nelecstate; ++i)
      d_cumulative_probabilities(index,i) += d_cumulative_probabilities(index,i-1);

    double ran = erandom.drand();
    ielec = 0;
    while (ran > d_cumulative_probabilities(index,ielec))
      ++ielec;
  }
  return ielec;
}

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
  int nstate = d_nelecstates[isp];

  // no states, just return

  if (nstate == 0) return;

  rand_type rand_gen = rand_pool.get_state();

  d_elecstate[index] = 0;
  d_eelec[index] = ielec(index,isp,temp_elec,rand_gen);

  rand_pool.free_state(rand_gen);
}

}

#endif
#endif

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

FixStyle(vibmode/kk,FixVibmodeKokkos)

#else

#ifndef SPARTA_FIX_VIBMODE_KOKKOS_H
#define SPARTA_FIX_VIBMODE_KOKKOS_H

#include "fix_vibmode.h"
#include "kokkos_type.h"
#include "particle_kokkos.h"
#include "Kokkos_Random.hpp"
#include "rand_pool_wrap.h"

namespace SPARTA_NS {

class FixVibmodeKokkos : public FixVibmode {
 public:
  FixVibmodeKokkos(class SPARTA *, int, char **);
  FixVibmodeKokkos(class SPARTA *);
  ~FixVibmodeKokkos();
  void pre_update_custom_kokkos();
  void update_custom(int, double, double, double, double *);

  KOKKOS_INLINE_FUNCTION
  void update_custom_kokkos(int, double, double, double, const double *) const;

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

  DAT::t_int_2d_lr d_vibmode;
};

/* ----------------------------------------------------------------------
   called when a particle with index is created
    or when temperature dependent properties need to be updated
   populate all vibrational modes and set evib = sum of mode energies
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixVibmodeKokkos::update_custom_kokkos(int index, double temp_thermal,
                                            double temp_rot, double temp_vib,
                                            const double *) const
{
  int isp = d_particles[index].ispecies;
  int nmode = d_species[isp].nvibmode;

  // no modes, just return

  if (nmode == 0) return;

  // single mode, evib already set by Particle::evib()
  // just convert evib back to mode level

  if (nmode == 1) {
    d_vibmode(index,0) = static_cast<int>
      (d_particles[index].evib / boltz /
       d_species[isp].vibtemp[0]);
    return;
  }

  // loop over modes and populate each
  // accumlate new total evib

  double evib = 0.0;

  rand_type rand_gen = rand_pool.get_state();

  for (int imode = 0; imode < nmode; imode++) {
    const int ivib = static_cast<int> (-log(rand_gen.drand()) * temp_vib /
                                       d_species[isp].vibtemp[imode]);
    d_vibmode(index,imode) = ivib;
    evib += ivib * boltz * d_species[isp].vibtemp[imode];
  }

  d_particles[index].evib = evib;

  rand_pool.free_state(rand_gen);
}

}

#endif
#endif

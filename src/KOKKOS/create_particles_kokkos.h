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

#ifdef COMMAND_CLASS

CommandStyle(create_particles/kk,CreateParticlesKokkos)

#else

#ifndef SPARTA_CREATE_PARTICLES_KOKKOS_H
#define SPARTA_CREATE_PARTICLES_KOKKOS_H

#include "create_particles.h"
#include "particle_kokkos.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "Kokkos_Random.hpp"
#include "rand_pool_wrap.h"

namespace SPARTA_NS {

  class CreateParticlesKokkos : public CreateParticles {

  public:
    CreateParticlesKokkos(class SPARTA *);
    ~CreateParticlesKokkos();

    void create_local(bigint);
    void create_local_twopass(bigint np) { create_local(np); };

#ifndef SPARTA_KOKKOS_EXACT
    Kokkos::Random_XorShift64_Pool<DeviceType> rand_pool;
    typedef typename Kokkos::Random_XorShift64_Pool<DeviceType>::generator_type rand_type;

    //Kokkos::Random_XorShift1024_Pool<DeviceType> rand_pool;
    //typedef typename Kokkos::Random_XorShift1024_Pool<DeviceType>::generator_type rand_type;
#else
    RandPoolWrap rand_pool;
    typedef RandWrap rand_type;
#endif
  };
}

#endif
#endif

/* ERROR/WARNING messages:
*/


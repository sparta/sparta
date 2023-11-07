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

#ifdef SURF_REACT_CLASS

SurfReactStyle(global/kk,SurfReactGlobalKokkos)

#else

#ifndef SPARTA_SURF_REACT_GLOBAL_KOKKOS_H
#define SPARTA_SURF_REACT_GLOBAL_KOKKOS_H

#include "surf_react_global.h"
#include "kokkos_type.h"
#include "rand_pool_wrap.h"
#include "Kokkos_Random.hpp"
#include "particle_kokkos.h"

namespace SPARTA_NS {

class SurfReactGlobalKokkos : public SurfReactGlobal {
 public:
  SurfReactGlobalKokkos(class SPARTA *, int, char **);
  SurfReactGlobalKokkos(class SPARTA *);
  ~SurfReactGlobalKokkos();
  void init();
  void tally_reset();
  void tally_update();

  void pre_react();
  void backup();
  void restore();

#ifndef SPARTA_KOKKOS_EXACT
  Kokkos::Random_XorShift64_Pool<DeviceType> rand_pool;
  typedef typename Kokkos::Random_XorShift64_Pool<DeviceType>::generator_type rand_type;

  //Kokkos::Random_XorShift1024_Pool<DeviceType> rand_pool;
  //typedef typename Kokkos::Random_XorShift1024_Pool<DeviceType>::generator_type rand_type;
#else
  RandPoolWrap rand_pool;
  typedef RandWrap rand_type;
#endif

  RanKnuth* random_backup;

  typedef Kokkos::DualView<int[3], DeviceType::array_layout, DeviceType> tdual_int_3;
  typedef tdual_int_3::t_dev t_int_3;
  typedef tdual_int_3::t_host t_host_int_3;
  t_int_3 d_scalars;
  t_host_int_3 h_scalars;

  DAT::t_int_scalar d_nsingle;
  DAT::t_int_1d d_tally_single;

  HAT::t_int_scalar h_nsingle;
  HAT::t_int_1d h_tally_single;

  t_particle_1d d_particles;

 public:

  /* ----------------------------------------------------------------------
     select surface reaction to perform for particle with ptr IP on surface
     return which reaction 1 (destroy), 2 (create), 0 = no reaction
     if create, add particle and return ptr JP
  ------------------------------------------------------------------------- */

  template<int ATOMIC_REDUCTION>
  KOKKOS_INLINE_FUNCTION
  int react_kokkos(Particle::OnePart *&ip, int, const double *,
                   Particle::OnePart *&jp, int &,
                   const DAT::t_int_scalar &d_retry, const DAT::t_int_scalar &d_nlocal) const
  {
    rand_type rand_gen = rand_pool.get_state();
    double r = rand_gen.drand();

    // perform destroy reaction

    if (r < prob_destroy) {
      if (ATOMIC_REDUCTION == 0) {
        d_nsingle()++;
        d_tally_single(0)++;
      } else {
        Kokkos::atomic_increment(&d_nsingle());
        Kokkos::atomic_increment(&d_tally_single(0));
      }
      ip = NULL;
      rand_pool.free_state(rand_gen);
      return 1;
    }

    // perform create reaction
    // clone 1st particle to create 2nd particle
    // if add_particle performs a realloc must retry

    if (r < prob_destroy+prob_create) {
      if (ATOMIC_REDUCTION == 0) {
        d_nsingle()++;
        d_tally_single(1)++;
      } else {
        Kokkos::atomic_increment(&d_nsingle());
        Kokkos::atomic_increment(&d_tally_single(1));
      }
      double x[3],v[3];
      int id = MAXSMALLINT*rand_gen.drand();
      memcpy(x,ip->x,3*sizeof(double));
      memcpy(v,ip->v,3*sizeof(double));

      int index;
      if (ATOMIC_REDUCTION == 0) {
        index = d_nlocal();
        d_nlocal()++;
      } else
        index = Kokkos::atomic_fetch_add(&d_nlocal(),1);

      int reallocflag = ParticleKokkos::add_particle_kokkos(d_particles,index,id,ip->ispecies,ip->icell,x,v,0.0,0.0);
      if (reallocflag) {
        d_retry() = 1;
        rand_pool.free_state(rand_gen);
        return 0;
      }
      jp = &d_particles[index];
      rand_pool.free_state(rand_gen);
      return 2;
    }

    // no reaction

    rand_pool.free_state(rand_gen);
    return 0;
  }

};

}

#endif
#endif

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

SurfReactStyle(prob/kk,SurfReactProbKokkos)

#else

#ifndef SPARTA_SURF_REACT_PROB_KOKKOS_H
#define SPARTA_SURF_REACT_PROB_KOKKOS_H

#include "surf_react_prob.h"
#include "kokkos_type.h"
#include "rand_pool_wrap.h"
#include "Kokkos_Random.hpp"
#include "particle_kokkos.h"

namespace SPARTA_NS {

enum{DISSOCIATION,EXCHANGE,RECOMBINATION};

class SurfReactProbKokkos : public SurfReactProb {
 public:
  SurfReactProbKokkos(class SPARTA *, int, char **);
  SurfReactProbKokkos(class SPARTA *);
  ~SurfReactProbKokkos();
  void init();
  void tally_reset();
  void tally_update();

  void pre_react();
  void backup();
  void restore();

 private:
  DAT::t_int_1d d_reactions_n;       // # of reactions in list
  DAT::t_int_2d d_list;

  DAT::t_int_1d d_type;
  DAT::t_int_2d d_reactants;
  DAT::t_int_2d d_products;
  DAT::t_float_2d d_coeffs;

  void init_reactions();

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

  DAT::t_int_1d d_scalars;
  HAT::t_int_1d h_scalars;

  DAT::t_int_scalar d_nsingle;
  DAT::t_int_1d d_tally_single;

  HAT::t_int_scalar h_nsingle;
  HAT::t_int_1d h_tally_single;

  t_particle_1d d_particles;

 public:

  /* ----------------------------------------------------------------------
     select surface reaction to perform for particle with ptr IP on surface
     return which reaction 1 to N, 0 = no reaction
     if dissociation, add particle and return ptr JP
  ------------------------------------------------------------------------- */

  template<int ATOMIC_REDUCTION>
  KOKKOS_INLINE_FUNCTION
  int react_kokkos(Particle::OnePart *&ip, int, const double *,
                   Particle::OnePart *&jp, int &,
                   const DAT::t_int_scalar &d_retry, const DAT::t_int_scalar &d_nlocal) const
  {
    int n = d_reactions_n[ip->ispecies];
    if (n == 0) return 0;

    // probablity to compare to reaction probability

    double react_prob = 0.0;
    rand_type rand_gen = rand_pool.get_state();
    double random_prob = rand_gen.drand();

    // loop over possible reactions for this species
    // if dissociation performs a realloc:
    //   make copy of x,v with new species
    //   rot/vib energies will be reset by SurfCollide
    //   repoint ip to new particles data struct if reallocated

    for (int i = 0; i < n; i++) {
      int j = d_list(ip->ispecies,i);
      react_prob += d_coeffs(j,0);

      if (react_prob > random_prob) {
        if (ATOMIC_REDUCTION == 0) {
          d_nsingle()++;
          d_tally_single(j)++;
        } else {
          Kokkos::atomic_increment(&d_nsingle());
          Kokkos::atomic_increment(&d_tally_single(j));
        }
        switch (d_type(j)) {
        case DISSOCIATION:
          {
            double x[3],v[3];
            ip->ispecies = d_products(j,0);
            int id = MAXSMALLINT*rand_gen.drand();
            memcpy(x,ip->x,3*sizeof(double));
            memcpy(v,ip->v,3*sizeof(double));

            int index;
            if (ATOMIC_REDUCTION == 0) {
              index = d_nlocal();
              d_nlocal()++;
            } else
              index = Kokkos::atomic_fetch_add(&d_nlocal(),1);

            int reallocflag = ParticleKokkos::add_particle_kokkos(d_particles,index,id,d_products(j,1),ip->icell,x,v,0.0,0.0);
            if (reallocflag) {
              d_retry() = 1;
              rand_pool.free_state(rand_gen);
              return 0;
            }
            jp = &d_particles[index];
            rand_pool.free_state(rand_gen);
            return (j + 1);
          }
        case EXCHANGE:
          {
            ip->ispecies = d_products(j,0);
            rand_pool.free_state(rand_gen);
            return (j + 1);
          }
        case RECOMBINATION:
          {
            ip = NULL;
            rand_pool.free_state(rand_gen);
            return (j + 1);
          }
        }
      }
    }

    // no reaction

    rand_pool.free_state(rand_gen);
    return 0;
  }

};

}

#endif
#endif

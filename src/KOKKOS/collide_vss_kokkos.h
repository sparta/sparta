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

#ifdef COLLIDE_CLASS

CollideStyle(vss/kk,CollideVSSKokkos)

#else

#ifndef SPARTA_COLLIDE_VSS_KOKKOS_H
#define SPARTA_COLLIDE_VSS_KOKKOS_H

#include "collide_vss.h"
#include "collide_vss_kokkos.h"
#include "particle_kokkos.h"
#include "grid_kokkos.h"
#include "kokkos_type.h"
#include "Kokkos_Random.hpp"
#include "rand_pool_wrap.h"
#include "kokkos_copy.h"

namespace SPARTA_NS {

struct s_COLLIDE_REDUCE {
  int nattempt_one,ncollide_one;
  KOKKOS_INLINE_FUNCTION
  s_COLLIDE_REDUCE() {
    nattempt_one = 0;
    ncollide_one = 0;
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(const s_COLLIDE_REDUCE &rhs) {
    nattempt_one += rhs.nattempt_one;
    ncollide_one += rhs.ncollide_one;
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(const volatile s_COLLIDE_REDUCE &rhs) volatile {
    nattempt_one += rhs.nattempt_one;
    ncollide_one += rhs.ncollide_one;
  }
};
typedef struct s_COLLIDE_REDUCE COLLIDE_REDUCE;

struct TagCollideResetVremax{};
struct TagCollideZeroNN{};

template < int NEARCP, int ATOMIC_REDUCTION >
struct TagCollideCollisionsOne{};

class CollideVSSKokkos : public CollideVSS {
 public:
  typedef ArrayTypes<DeviceType> AT;
  typedef COLLIDE_REDUCE value_type;

  CollideVSSKokkos(class SPARTA *, int, char **);
  ~CollideVSSKokkos();
  void init();
  void reset_vremax();
  void collisions();
  void sync(ExecutionSpace, unsigned int);
  void modify(ExecutionSpace, unsigned int);

#ifndef SPARTA_KOKKOS_EXACT
  Kokkos::Random_XorShift64_Pool<DeviceType> rand_pool;
  typedef typename Kokkos::Random_XorShift64_Pool<DeviceType>::generator_type rand_type;

  //Kokkos::Random_XorShift1024_Pool<DeviceType> rand_pool;
  //typedef typename Kokkos::Random_XorShift1024_Pool<DeviceType>::generator_type rand_type;
#else
  RandPoolWrap rand_pool;
  typedef RandWrap rand_type;
#endif

  KOKKOS_INLINE_FUNCTION
  double attempt_collision_kokkos(int, int, double, rand_type &) const;
  KOKKOS_INLINE_FUNCTION
  int test_collision_kokkos(int, int, int, Particle::OnePart *, Particle::OnePart *, struct State &, rand_type &) const;
  KOKKOS_INLINE_FUNCTION
  void setup_collision_kokkos(Particle::OnePart *, Particle::OnePart *, struct State &, struct State &) const;
  KOKKOS_INLINE_FUNCTION
  int perform_collision_kokkos(Particle::OnePart *&, Particle::OnePart *&, 
                        Particle::OnePart *&, struct State &, struct State &, rand_type &) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagCollideResetVremax, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagCollideZeroNN, const int&) const;

  template < int NEARCP, int ATOMIC_REDUCTION >
  KOKKOS_INLINE_FUNCTION
  void operator()(TagCollideCollisionsOne< NEARCP, ATOMIC_REDUCTION >, const int&) const;

  template < int NEARCP, int ATOMIC_REDUCTION >
  KOKKOS_INLINE_FUNCTION
  void operator()(TagCollideCollisionsOne< NEARCP, ATOMIC_REDUCTION >, const int&, COLLIDE_REDUCE&) const;

 private:
  int pack_grid_one(int, char *, int);
  int unpack_grid_one(int, char *);
  void compress_grid();
  void adapt_grid();
  void grow_percell(int);

  KKCopy<GridKokkos> grid_kk_copy;

  t_particle_1d d_particles;
  t_species_1d d_species;

  DAT::tdual_int_1d k_next;
  typename AT::t_int_1d d_next;

  DAT::tdual_float_2d k_vremax_initial;
  typename AT::t_float_2d d_vremax_initial;
  DAT::tdual_float_3d k_vremax;
  typename AT::t_float_3d d_vremax;
  DAT::tdual_float_3d k_remain;
  typename AT::t_float_3d d_remain;

  DAT::tdual_int_scalar k_nattempt_one;
  typename AT::t_int_scalar d_nattempt_one;
  HAT::t_int_scalar h_nattempt_one;

  DAT::tdual_int_scalar k_ncollide_one;
  typename AT::t_int_scalar d_ncollide_one;
  HAT::t_int_scalar h_ncollide_one;

  DAT::tdual_int_scalar k_error_flag;
  typename AT::t_int_scalar d_error_flag;
  HAT::t_int_scalar h_error_flag;

  typename AT::t_int_2d d_nn_last_partner; 

  template < int NEARCP > void collisions_one(COLLIDE_REDUCE&);

  // VSS specific

  DAT::tdual_float_2d k_prefactor;
  typename AT::t_float_2d d_prefactor;

  typedef Kokkos::
    DualView<Params*, Kokkos::LayoutRight, DeviceType> tdual_params_1d;
  typedef tdual_params_1d::t_dev t_params_1d;
  tdual_params_1d k_params;
  t_params_1d d_params;


  double dt,fnum,boltz;

  KOKKOS_INLINE_FUNCTION
  void SCATTER_TwoBodyScattering(Particle::OnePart *, 
                                 Particle::OnePart *,
                                 struct State &, struct State &, rand_type &) const;
  KOKKOS_INLINE_FUNCTION
  void EEXCHANGE_NonReactingEDisposal(Particle::OnePart *, 
                                      Particle::OnePart *,
                                      struct State &, struct State &, rand_type &) const;

  KOKKOS_INLINE_FUNCTION
  double sample_bl(rand_type &, double, double) const;
  KOKKOS_INLINE_FUNCTION
  double rotrel (int, double) const;
  KOKKOS_INLINE_FUNCTION
  double vibrel (int, double) const;

  KOKKOS_INLINE_FUNCTION
  int find_nn(rand_type &, int, int, int) const;
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/

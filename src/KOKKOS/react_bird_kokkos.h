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

#ifndef SPARTA_REACT_BIRD_KOKKOS_H
#define SPARTA_REACT_BIRD_KOKKOS_H

#include "react_bird.h"
#include "kokkos_type.h"
#include "Kokkos_Random.hpp"
#include "rand_pool_wrap.h"

#define MAXREACTANT 2
#define MAXPRODUCT 3
#define MAXCOEFF 7               // 5 in file, extra for pre-computation

namespace SPARTA_NS {

class ReactBirdKokkos : public ReactBird {
 public:
  ReactBirdKokkos(class SPARTA *, int, char **);
  ReactBirdKokkos(class SPARTA* sparta) : ReactBird(sparta),
                                          rand_pool(1
#ifdef SPARTA_KOKKOS_EXACT
                                          , sparta
#endif
                                          ) {random_backup = NULL;};
  virtual ~ReactBirdKokkos();
  virtual void init();
  virtual int attempt(Particle::OnePart *, Particle::OnePart *,
                      double, double, double, double &, int &) = 0;
  double extract_tally(int);
  void backup();
  void restore();

  // tallies for reactions

  DAT::tdual_bigint_1d k_tally_reactions;
  DAT::t_bigint_1d d_tally_reactions;
  DAT::tdual_bigint_1d k_tally_reactions_all;
  DAT::t_bigint_1d d_tally_reactions_backup;

  struct OneReactionKokkos {
    int active;                    // 1 if reaction is active
    int initflag;                  // 1 if reaction params have been init
    int type;                      // reaction type = DISSOCIATION, etc
    int style;                     // reaction style = ARRHENIUS, etc
    int ncoeff;                    // # of numerical coeffs
    int nreactant,nproduct;        // # of reactants and products
    int d_reactants[MAXREACTANT],d_products[MAXPRODUCT];      // species indices of reactants/products
    double d_coeff[MAXCOEFF];                 // numerical coeffs for reaction
  };

  // all reactions a pair of reactant species is part of

  struct ReactionIJKokkos {
    DAT::t_int_1d d_list;       // N-length list of rlist indices
                     //   for reactions defined for this IJ pair,
                     //   just a ptr into sub-section of long list_ij vector
                     //   for all pairs
    DAT::t_int_1d d_sp2recomb;  // Nspecies-length list of rlist indices
                     //   for recomb reactions defined for this IJ pair,
                     //   one index for all 3rd particle species,
                     //   just a ptr into sub-section of long sp2recomb_ij
                     //   vector for all pairs which have recomb reactions
    int n;           // # of reactions in list
  };

 protected:

  typedef Kokkos::
    DualView<OneReactionKokkos*, DeviceType::array_layout, DeviceType> tdual_reaction_1d;
  typedef tdual_reaction_1d::t_dev t_reaction_1d;
  typedef tdual_reaction_1d::t_host t_host_reaction_1d;

  t_reaction_1d d_rlist;              // list of all reactions read from file

  typedef Kokkos::
    DualView<ReactionIJKokkos**, DeviceType::array_layout, DeviceType> tdual_reactionIJ_2d;
  typedef tdual_reactionIJ_2d::t_dev t_reactionIJ_2d;
  typedef tdual_reactionIJ_2d::t_host t_host_reactionIJ_2d;

  tdual_reactionIJ_2d k_reactions;     // reaction info for all IJ pairs of species
  t_reactionIJ_2d d_reactions;     // reaction info for all IJ pairs of species

  RanKnuth* random_backup;

 public:
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

/* ERROR/WARNING messages:

*/

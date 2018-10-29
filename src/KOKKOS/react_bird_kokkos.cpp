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
#include "string.h"
#include "stdlib.h"
#include "react_bird_kokkos.h"
#include "input.h"
#include "collide.h"
#include "update.h"
#include "particle.h"
#include "comm.h"
#include "modify.h"
#include "fix.h"
#include "fix_ambipolar.h"
#include "random_park.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{DISSOCIATION,EXCHANGE,IONIZATION,RECOMBINATION};  // other react files
enum{ARRHENIUS,QUANTUM};                               // other react files

#define MAXLINE 1024
#define DELTALIST 16

/* ---------------------------------------------------------------------- */

ReactBirdKokkos::ReactBirdKokkos(SPARTA *sparta, int narg, char **arg) :
  ReactBird(sparta, narg, arg),
  rand_pool(12345 + comm->me
#ifdef SPARTA_KOKKOS_EXACT
            , sparta
#endif
            )
{

}

/* ---------------------------------------------------------------------- */

ReactBirdKokkos::~ReactBirdKokkos()
{
  if (copy) return;

  // deallocate views of views in serial to prevent race conditions in external tools

  //for (int i = 0; i < maxlist; i++) {
  //  d_rlist(i).d_reactants = DAT::t_int_1d();
  //  d_rlist(i).d_products = DAT::t_int_1d();
  //  d_rlist(i).d_coeff = DAT::t_int_1d();
  //}

  int nspecies = k_reactions.h_view.extent(0);
  for (int i = 0; i < nspecies; i++) {
    for (int j = 0; j < nspecies; j++) {
      if (d_reactions.data()) {
        k_reactions.h_view(i,j).d_list = DAT::t_int_1d();
        k_reactions.h_view(i,j).d_sp2recomb = DAT::t_int_1d();
      }
    }
  }

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.destroy();
#endif
}

/* ---------------------------------------------------------------------- */

void ReactBirdKokkos::init() 
{
  ReactBird::init();

  // copy data into device views

  d_rlist = t_reaction_1d("react/bird:rlist",maxlist);
  auto h_rlist = Kokkos::create_mirror_view(d_rlist);
  for (int i = 0; i < maxlist; i++) {
    h_rlist[i].active = rlist[i].active;
    h_rlist[i].initflag = rlist[i].initflag;
    h_rlist[i].type = rlist[i].type;
    h_rlist[i].style = rlist[i].style;
    h_rlist[i].ncoeff = rlist[i].ncoeff;
    h_rlist[i].nreactant = rlist[i].nreactant;
    h_rlist[i].nproduct = rlist[i].nproduct;
    for (int j = 0; j < MAXREACTANT; j++)
      h_rlist[i].d_reactants[j] = rlist[i].reactants[j];
    for (int j = 0; j < MAXPRODUCT; j++)
      h_rlist[i].d_products[j] = rlist[i].products[j];
    for (int j = 0; j < MAXCOEFF; j++)
      h_rlist[i].d_coeff[j] = rlist[i].coeff[j];
  }
  Kokkos::deep_copy(d_rlist,h_rlist);

  int nspecies = particle->nspecies;
  // Doesn't work with deep_copy instead of k_reactions DualView, not sure why
  k_reactions = tdual_reactionIJ_2d("react/bird:reactions",nspecies,nspecies);
  for (int i = 0; i < nspecies; i++) {
    for (int j = 0; j < nspecies; j++) {
      const int n = reactions[i][j].n;
      k_reactions.h_view(i,j).n = n;

      k_reactions.h_view(i,j).d_list = DAT::t_int_1d("react/bird:list",n);
      auto h_list = Kokkos::create_mirror_view(k_reactions.h_view(i,j).d_list);
      for (int k = 0; k < n; k++)
        h_list(k) = reactions[i][j].list[k];
      Kokkos::deep_copy(k_reactions.h_view(i,j).d_list,h_list);

      if (!recombflag) continue;

      k_reactions.h_view(i,j).d_sp2recomb = DAT::t_int_1d("react/bird:sp2recomb",nspecies);
      auto h_sp2recomb = Kokkos::create_mirror_view(k_reactions.h_view(i,j).d_sp2recomb);
      for (int k = 0; k < nspecies; k++)
        h_sp2recomb(k) = reactions[i][j].sp2recomb[k];
      Kokkos::deep_copy(k_reactions.h_view(i,j).d_sp2recomb,h_sp2recomb);
    }
  }
  k_reactions.modify<SPAHostType>();
  k_reactions.sync<DeviceType>();
  d_reactions = k_reactions.d_view;

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.init(random);
#endif
}

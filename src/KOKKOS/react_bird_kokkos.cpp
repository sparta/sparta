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
#include "random_knuth.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "kokkos.h"
#include "memory_kokkos.h"

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
  delete [] tally_reactions;
  delete [] tally_reactions_all;
  memoryKK->create_kokkos(k_tally_reactions,tally_reactions,nlist,"react_bird:tally_reactions");
  memoryKK->create_kokkos(k_tally_reactions_all,tally_reactions_all,nlist,"react_bird:tally_reactions_all");
  d_tally_reactions = k_tally_reactions.d_view;

  random_backup = NULL;
}

/* ---------------------------------------------------------------------- */

ReactBirdKokkos::~ReactBirdKokkos()
{
  if (copy) return;

  tally_reactions = NULL;
  tally_reactions_all = NULL;

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.destroy();
  if (random_backup)
    delete random_backup;
#endif

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

      if (!recombflag || !reactions[i][j].sp2recomb) continue;

      k_reactions.h_view(i,j).d_sp2recomb = DAT::t_int_1d("react/bird:sp2recomb",nspecies);
      auto h_sp2recomb = Kokkos::create_mirror_view(k_reactions.h_view(i,j).d_sp2recomb);
      for (int k = 0; k < nspecies; k++)
        h_sp2recomb(k) = reactions[i][j].sp2recomb[k];
      Kokkos::deep_copy(k_reactions.h_view(i,j).d_sp2recomb,h_sp2recomb);
    }
  }
  k_reactions.modify_host();
  k_reactions.sync_device();
  d_reactions = k_reactions.d_view;

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.init(random);
#endif
}

/* ----------------------------------------------------------------------
   return tally associated with a reaction
------------------------------------------------------------------------- */

double ReactBirdKokkos::extract_tally(int m)
{
  if (!tally_flag) {
    tally_flag = 1;

    if (sparta->kokkos->gpu_aware_flag) {
      MPI_Allreduce(d_tally_reactions.data(),k_tally_reactions_all.d_view.data(),nlist,
                    MPI_SPARTA_BIGINT,MPI_SUM,world);
      k_tally_reactions_all.modify_device();
      k_tally_reactions_all.sync_host();
    } else {
      k_tally_reactions.modify_device();
      k_tally_reactions.sync_host();
      MPI_Allreduce(k_tally_reactions.h_view.data(),k_tally_reactions_all.h_view.data(),nlist,
                    MPI_SPARTA_BIGINT,MPI_SUM,world);
    }

  }

  return 1.0*tally_reactions_all[m];
};

/* ---------------------------------------------------------------------- */

void ReactBirdKokkos::backup()
{
  d_tally_reactions_backup = decltype(d_tally_reactions)(Kokkos::view_alloc("react_bird:tally_reactions_backup",Kokkos::WithoutInitializing),d_tally_reactions.extent(0));

  Kokkos::deep_copy(d_tally_reactions_backup,d_tally_reactions);

#ifdef SPARTA_KOKKOS_EXACT
  if (!random_backup)
    random_backup = new RanKnuth(12345 + comm->me);
  memcpy(random_backup,random,sizeof(RanKnuth));
#endif

}

/* ---------------------------------------------------------------------- */

void ReactBirdKokkos::restore()
{
  Kokkos::deep_copy(d_tally_reactions,d_tally_reactions_backup);

#ifdef SPARTA_KOKKOS_EXACT
  memcpy(random,random_backup,sizeof(RanKnuth));
#endif

  d_tally_reactions_backup = decltype(d_tally_reactions_backup)();
}

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
#include "surf_react_prob_kokkos.h"
#include "input.h"
#include "update.h"
#include "comm.h"
#include "error.h"
#include "particle_kokkos.h"
#include "sparta_masks.h"

using namespace SPARTA_NS;

#define MAXREACTANT 1
#define MAXPRODUCT 2
#define MAXCOEFF 2


/* ---------------------------------------------------------------------- */

SurfReactProbKokkos::SurfReactProbKokkos(SPARTA *sparta, int narg, char **arg) :
  SurfReactProb(sparta,narg,arg),
  rand_pool(12345 + comm->me
#ifdef SPARTA_KOKKOS_EXACT
            , sparta
#endif
            )
{
  kokkosable = 1;

  d_scalars = DAT::t_int_1d("surf_react:scalars",nlist+1);
  d_nsingle = Kokkos::subview(d_scalars,0);
  d_tally_single = Kokkos::subview(d_scalars,std::make_pair(1,nlist+1));

  h_scalars = HAT::t_int_1d("surf_react:scalars_mirror",nlist+1);
  h_nsingle = Kokkos::subview(h_scalars,0);
  h_tally_single = Kokkos::subview(h_scalars,std::make_pair(1,nlist+1));

  random_backup = NULL;
}

SurfReactProbKokkos::SurfReactProbKokkos(SPARTA *sparta) :
  SurfReactProb(sparta),
  rand_pool(12345 // seed will be copied over
#ifdef SPARTA_KOKKOS_EXACT
            , sparta
#endif
            )
{
  random = NULL;
  random_backup = NULL;
  
  id = NULL;
  style = NULL;
  tally_single = NULL;
  tally_total = NULL;
  tally_single_all = NULL;
  tally_total_all = NULL;

  rlist = NULL;
  reactions = NULL;
  indices = NULL;
}

/* ---------------------------------------------------------------------- */

SurfReactProbKokkos::~SurfReactProbKokkos()
{
 if (copy) return;

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.destroy();
  if (random_backup)
    delete random_backup;
#endif
}

/* ---------------------------------------------------------------------- */

void SurfReactProbKokkos::init()
{
  SurfReactProb::init();

  Kokkos::deep_copy(d_scalars,0);

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.init(random);
#endif
}

/* ---------------------------------------------------------------------- */

void SurfReactProbKokkos::tally_reset()
{
  SurfReact::tally_reset();

  Kokkos::deep_copy(d_scalars,0);
}

/* ---------------------------------------------------------------------- */

void SurfReactProbKokkos::tally_update()
{
  Kokkos::deep_copy(h_scalars,d_scalars);
  ntotal += h_nsingle();
  for (int i = 0; i < nlist; i++) tally_total[i] += h_tally_single[i];
}

/* ---------------------------------------------------------------------- */

void SurfReactProbKokkos::init_reactions()
{
  SurfReactProb::init_reactions();

  int nspecies = particle->nspecies;

  int nmax = 0;

  d_reactions_n = DAT::t_int_1d("surf_react_prob:offsets",nspecies);
  auto h_reactions_n = Kokkos::create_mirror_view(d_reactions_n);

  for (int i = 0; i < nspecies; i++) {
    int n = reactions[i].n;
    h_reactions_n(i) = n;
    nmax = MAX(nmax,n);
  }

  d_list = DAT::t_int_2d("surf_react_prob:list",nspecies,nmax);
  auto h_list = Kokkos::create_mirror_view(d_list);

  for (int i = 0; i < nspecies; i++)
    for (int j = 0; j < reactions[i].n; j++)
      h_list(i,j) = reactions[i].list[j];

  d_type = DAT::t_int_1d("surf_react_prob:type",maxlist_prob);
  d_reactants = DAT::t_int_2d("surf_react_prob:reactants",maxlist_prob,MAXREACTANT);
  d_products = DAT::t_int_2d("surf_react_prob:products",maxlist_prob,MAXPRODUCT);
  d_coeffs = DAT::t_float_2d("surf_react_prob:coeffs",maxlist_prob,MAXCOEFF);

  auto h_type = Kokkos::create_mirror_view(d_type);
  auto h_reactants = Kokkos::create_mirror_view(d_reactants);
  auto h_products = Kokkos::create_mirror_view(d_products);
  auto h_coeffs = Kokkos::create_mirror_view(d_coeffs);

  for (int i = 0; i < nlist_prob; i++) {
    OneReaction *r = &rlist[i];

    h_type(i) = r->type;

    for (int j = 0; j < MAXREACTANT; j++)
      h_reactants(i,j) = r->reactants[j];

    for (int j = 0; j < MAXPRODUCT; j++)
      h_products(i,j) = r->products[j];

    for (int j = 0; j < MAXCOEFF; j++)
      h_coeffs(i,j) = r->coeff[j];
  }

  Kokkos::deep_copy(d_reactions_n,h_reactions_n);
  Kokkos::deep_copy(d_list,h_list);
  Kokkos::deep_copy(d_type,h_type);
  Kokkos::deep_copy(d_reactants,h_reactants);
  Kokkos::deep_copy(d_products,h_products);
  Kokkos::deep_copy(d_coeffs,h_coeffs);
}

/* ---------------------------------------------------------------------- */

void SurfReactProbKokkos::pre_react()
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,PARTICLE_MASK);
  d_particles = particle_kk->k_particles.d_view;
}

/* ---------------------------------------------------------------------- */

void SurfReactProbKokkos::backup()
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  d_particles = particle_kk->k_particles.d_view;

#ifdef SPARTA_KOKKOS_EXACT
  if (!random_backup)
    random_backup = new RanKnuth(12345 + comm->me);
  memcpy(random_backup,random,sizeof(RanKnuth));
#endif
}

/* ---------------------------------------------------------------------- */

void SurfReactProbKokkos::restore()
{
  Kokkos::deep_copy(d_scalars,0);

#ifdef SPARTA_KOKKOS_EXACT
  memcpy(random,random_backup,sizeof(RanKnuth));
#endif
}

/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "surf_collide_adiabatic_kokkos.h"
#include "surf_kokkos.h"
#include "input.h"
#include "variable.h"
#include "particle.h"
#include "domain.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "math_const.h"
#include "math_extra.h"
#include "error.h"
#include "particle_kokkos.h"
#include "sparta_masks.h"
#include "collide.h"

using namespace SPARTA_NS;
using namespace MathConst;

#define VAL_1(X) X
#define VAL_2(X) VAL_1(X), VAL_1(X)

/* ---------------------------------------------------------------------- */

SurfCollideAdiabaticKokkos::SurfCollideAdiabaticKokkos(SPARTA *sparta, int narg, char **arg) :
  SurfCollideAdiabatic(sparta, narg, arg),
  fix_ambi_kk_copy(sparta),
  fix_vibmode_kk_copy(sparta),
#ifdef SPARTA_KOKKOS_FIXED_LISTS
  sr_kk_global_copy{VAL_2(KKCopy<SurfReactGlobalKokkos>(sparta))},
  sr_kk_prob_copy{VAL_2(KKCopy<SurfReactProbKokkos>(sparta))},
  sr_kk_adsorb_copy{VAL_2(KKCopy<SurfReactAdsorbKokkos>(sparta))},
#endif
  rand_pool(12345 + comm->me
#ifdef SPARTA_KOKKOS_EXACT
            , sparta
#endif
           )
{
  kokkosable = 1;

  random_backup = NULL;

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.init(random);

  // allocate on the real class instance: backup() is only ever invoked on a
  //  KKCopy of this class, and a pointer stored there is dropped by the next
  //  copy() of the original over it

  random_backup = new RanKnuth(12345 + comm->me);
#endif

  // use 1D view for scalars to reduce GPU memory operations

  d_scalars = t_bigint_2("surf_collide_adiabatic:scalars");
  d_nsingle = Kokkos::subview(d_scalars,0);
  d_nreact_one = Kokkos::subview(d_scalars,1);

  h_scalars = t_host_bigint_2("surf_collide_adiabatic:scalars_mirror");
  h_nsingle = Kokkos::subview(h_scalars,0);
  h_nreact_one = Kokkos::subview(h_scalars,1);
}

SurfCollideAdiabaticKokkos::SurfCollideAdiabaticKokkos(SPARTA *sparta) :
  SurfCollideAdiabatic(sparta),
  fix_ambi_kk_copy(sparta),
  fix_vibmode_kk_copy(sparta),
#ifdef SPARTA_KOKKOS_FIXED_LISTS
  sr_kk_global_copy{VAL_2(KKCopy<SurfReactGlobalKokkos>(sparta))},
  sr_kk_prob_copy{VAL_2(KKCopy<SurfReactProbKokkos>(sparta))},
  sr_kk_adsorb_copy{VAL_2(KKCopy<SurfReactAdsorbKokkos>(sparta))},
#endif
  rand_pool(12345 // seed doesn't matter since it will just be copied over
#ifdef SPARTA_KOKKOS_EXACT
            , sparta
#endif
           )
{
  copy = 1;
}

/* ---------------------------------------------------------------------- */

SurfCollideAdiabaticKokkos::~SurfCollideAdiabaticKokkos()
{
  if (copy) return;

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.destroy();
  if (random_backup)
    delete random_backup;
#endif
}

/* ---------------------------------------------------------------------- */

void SurfCollideAdiabaticKokkos::init()
{
  SurfCollideAdiabatic::init();

  // scan the fix list directly rather than testing modify->n_update_custom
  //  first: SPARTA::init() runs surf->init() before modify->init(), so that
  //  count still holds its value from the previous run (0 on the first one)

  ambi_flag = vibmode_flag = 0;
  afix_kk = NULL;
  vfix_kk = NULL;

  for (int ifix = 0; ifix < modify->nfix; ifix++) {
    if (strcmp(modify->fix[ifix]->style,"ambipolar") == 0) {
      ambi_flag = 1;
      FixAmbipolar *afix = (FixAmbipolar *) modify->fix[ifix];
      if (!afix->kokkos_flag)
        error->all(FLERR,"Must use fix ambipolar/kk when Kokkos is enabled");
      afix_kk = (FixAmbipolarKokkos*)afix;
    } else if (strcmp(modify->fix[ifix]->style,"vibmode") == 0) {
      vibmode_flag = 1;
      FixVibmode *vfix = (FixVibmode *) modify->fix[ifix];
      if (!vfix->kokkos_flag)
        error->all(FLERR,"Must use fix vibmode/kk when Kokkos is enabled");
      vfix_kk = (FixVibmodeKokkos*)vfix;
    }
  }
}

/* ---------------------------------------------------------------------- */

void SurfCollideAdiabaticKokkos::pre_collide()
{
  if (ambi_flag) {
    afix_kk->pre_update_custom_kokkos();
    fix_ambi_kk_copy.copy(afix_kk);
  }

  if (vibmode_flag) {
    vfix_kk->pre_update_custom_kokkos();
    fix_vibmode_kk_copy.copy(vfix_kk);
  }

#ifdef SPARTA_KOKKOS_FIXED_LISTS
  if (surf->nsr > KOKKOS_MAX_TOT_SURF_REACT)
    error->all(FLERR,"Kokkos currently supports a limited number of surface reaction methods");
#else

  // the buffers must be sized before anything is blitted into them.  surf->nsr
  //   bounds the count of every individual style, so sizing all of them to it
  //   needs no counting pass, and the loop below still runs pre_react() in
  //   surf react list order

  sr_idx_resize(k_sr_type_list,d_sr_type_list,surf->nsr);
  sr_idx_resize(k_sr_map,d_sr_map,surf->nsr);
  sr_buf_resize<SurfReactGlobalKokkos>(k_sr_global,d_sr_global,surf->nsr);
  sr_buf_resize<SurfReactProbKokkos>(k_sr_prob,d_sr_prob,surf->nsr);
  sr_buf_resize<SurfReactAdsorbKokkos>(k_sr_adsorb,d_sr_adsorb,surf->nsr);
#endif

  if (surf->nsr > 0) {
    int nglob,nprob,nadsorb;
    nglob = nprob = nadsorb = 0;
    for (int n = 0; n < surf->nsr; n++) {
      if (!surf->sr[n]->kokkosable)
        error->all(FLERR,"Must use Kokkos-enabled surface reaction method with Kokkos");
      if (strcmp(surf->sr[n]->style,"global") == 0) {
#ifdef SPARTA_KOKKOS_FIXED_LISTS
        if (nglob >= KOKKOS_MAX_SURF_REACT_PER_TYPE)
          error->all(FLERR,"Kokkos currently supports two instances of each surface reaction method");
        sr_kk_global_copy[nglob].copy((SurfReactGlobalKokkos*)(surf->sr[n]));
#else
        sr_buf_blit(k_sr_global,nglob,(SurfReactGlobalKokkos*)(surf->sr[n]));
#endif
        KK_SR_H_GLOBAL(nglob).pre_react();
        KK_SR_H_TYPE(n) = 0;
        KK_SR_H_MAP(n) = nglob;
        nglob++;
      } else if (strcmp(surf->sr[n]->style,"prob") == 0) {
#ifdef SPARTA_KOKKOS_FIXED_LISTS
        if (nprob >= KOKKOS_MAX_SURF_REACT_PER_TYPE)
          error->all(FLERR,"Kokkos currently supports two instances of each surface reaction method");
        sr_kk_prob_copy[nprob].copy((SurfReactProbKokkos*)(surf->sr[n]));
#else
        sr_buf_blit(k_sr_prob,nprob,(SurfReactProbKokkos*)(surf->sr[n]));
#endif
        KK_SR_H_PROB(nprob).pre_react();
        KK_SR_H_TYPE(n) = 1;
        KK_SR_H_MAP(n) = nprob;
        nprob++;
      } else if (strcmp(surf->sr[n]->style,"adsorb") == 0) {
#ifdef SPARTA_KOKKOS_FIXED_LISTS
        if (nadsorb >= KOKKOS_MAX_SURF_REACT_PER_TYPE)
          error->all(FLERR,"Kokkos currently supports two instances of each surface reaction method");
        sr_kk_adsorb_copy[nadsorb].copy((SurfReactAdsorbKokkos*)(surf->sr[n]));
#else
        sr_buf_blit(k_sr_adsorb,nadsorb,(SurfReactAdsorbKokkos*)(surf->sr[n]));
#endif
        KK_SR_H_ADSORB(nadsorb).pre_react();
        KK_SR_H_TYPE(n) = 2;
        KK_SR_H_MAP(n) = nadsorb;
        nadsorb++;
      } else {
        error->all(FLERR,"Unknown Kokkos surface reaction method");
      }
    }

#ifndef SPARTA_KOKKOS_FIXED_LISTS

    // the models were blitted into the host image of the buffers and their
    //   pre_react() ran there; push the result to the device

    sr_buf_sync(k_sr_global,d_sr_global);
    sr_buf_sync(k_sr_prob,d_sr_prob);
    sr_buf_sync(k_sr_adsorb,d_sr_adsorb);
    sr_idx_sync(k_sr_type_list,d_sr_type_list);
    sr_idx_sync(k_sr_map,d_sr_map);
#endif

  }

  if (random == NULL) {
    // initialize RNG

    random = new RanKnuth(update->ranmaster->uniform());
    double seed = update->ranmaster->uniform();
    random->reset(seed,comm->me,100);

#ifdef SPARTA_KOKKOS_EXACT
    rand_pool.init(random);
#endif
  }

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,PARTICLE_MASK|SPECIES_MASK);
  d_particles = particle_kk->k_particles.view_device();
  d_species = particle_kk->k_species.view_device();

  Kokkos::deep_copy(d_scalars,0);
}

/* ---------------------------------------------------------------------- */

void SurfCollideAdiabaticKokkos::post_collide()
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  if (ambi_flag || vibmode_flag) particle_kk->modify(Device,CUSTOM_MASK);

  Kokkos::deep_copy(h_scalars,d_scalars);

  int m = surf->find_collide(id);
  auto sc = surf->sc[m]; // can't modify the copy directly, use the original
  sc->nsingle += h_nsingle();
  surf->nreact_one += h_nreact_one();

  d_particles = {};
  d_species = {};

  // pre_collide() runs on this KKCopy and has each active surf react model
  //  retain a reference to the particle list.  Release it before the next
  //  copy() blits over the member: a blit does not release, so the reference
  //  would be orphaned and its allocation never freed

  for (int n = 0; n < surf->nsr; n++) {
    if (KK_SR_H_TYPE(n) == 0) KK_SR_H_GLOBAL(KK_SR_H_MAP(n)).post_react();
    else if (KK_SR_H_TYPE(n) == 1) KK_SR_H_PROB(KK_SR_H_MAP(n)).post_react();
    else KK_SR_H_ADSORB(KK_SR_H_MAP(n)).post_react();
  }
}

/* ---------------------------------------------------------------------- */

void SurfCollideAdiabaticKokkos::backup()
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  d_particles = particle_kk->k_particles.view_device();

  if (surf->nsr > 0) {
    int nglob,nprob,nadsorb;
    nglob = nprob = nadsorb = 0;
    for (int n = 0; n < surf->nsr; n++) {
      if (strcmp(surf->sr[n]->style,"global") == 0) {
        KK_SR_H_GLOBAL(nglob).backup();
        nglob++;
      } else if (strcmp(surf->sr[n]->style,"prob") == 0) {
        KK_SR_H_PROB(nprob).backup();
        nprob++;
      } else if (strcmp(surf->sr[n]->style,"adsorb") == 0) {
        KK_SR_H_ADSORB(nadsorb).backup();
        nadsorb++;
      }
    }

#ifndef SPARTA_KOKKOS_FIXED_LISTS

    // backup() rewrites members of each model -- d_particles above all, which
    //   a grow reallocates -- so the device image of the buffers is stale
    //   until it is pushed again

    sr_buf_sync(k_sr_global,d_sr_global);
    sr_buf_sync(k_sr_prob,d_sr_prob);
    sr_buf_sync(k_sr_adsorb,d_sr_adsorb);
#endif
  }

#ifdef SPARTA_KOKKOS_EXACT
  if (!random_backup)
    random_backup = new RanKnuth(12345 + comm->me);
  memcpy(random_backup,random,sizeof(RanKnuth));
#endif
}

/* ---------------------------------------------------------------------- */

void SurfCollideAdiabaticKokkos::restore()
{
  if (surf->nsr > 0) {
    int nglob,nprob,nadsorb;
    nglob = nprob = nadsorb = 0;
    for (int n = 0; n < surf->nsr; n++) {
      if (strcmp(surf->sr[n]->style,"global") == 0) {
        KK_SR_H_GLOBAL(nglob).restore();
        nglob++;
      } else if (strcmp(surf->sr[n]->style,"prob") == 0) {
        KK_SR_H_PROB(nprob).restore();
        nprob++;
      } else if (strcmp(surf->sr[n]->style,"adsorb") == 0) {
        KK_SR_H_ADSORB(nadsorb).restore();
        nadsorb++;
      }
    }

#ifndef SPARTA_KOKKOS_FIXED_LISTS

    // restore() writes no member of a model today, only deep_copies into Views
    //   the model already holds, but the buffers are pushed for the same reason
    //   backup() pushes them: the device image must never trail the host one

    sr_buf_sync(k_sr_global,d_sr_global);
    sr_buf_sync(k_sr_prob,d_sr_prob);
    sr_buf_sync(k_sr_adsorb,d_sr_adsorb);
#endif
  }

  Kokkos::deep_copy(d_scalars,0);

#ifdef SPARTA_KOKKOS_EXACT
  memcpy(random,random_backup,sizeof(RanKnuth));
#endif
}

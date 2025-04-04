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
#include "surf_collide_diffuse_kokkos.h"
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

enum{INT,DOUBLE};                        // several files
enum{NUMERIC,CUSTOM,VARIABLE,VAREQUAL,VARSURF};   // surf_collide classes

#define VAL_1(X) X
#define VAL_2(X) VAL_1(X), VAL_1(X)

/* ---------------------------------------------------------------------- */

SurfCollideDiffuseKokkos::SurfCollideDiffuseKokkos(SPARTA *sparta, int narg, char **arg) :
  SurfCollideDiffuse(sparta, narg, arg),
  fix_ambi_kk_copy(sparta),
  fix_vibmode_kk_copy(sparta),
  sr_kk_global_copy{VAL_2(KKCopy<SurfReactGlobalKokkos>(sparta))},
  sr_kk_prob_copy{VAL_2(KKCopy<SurfReactProbKokkos>(sparta))},
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
#endif

  // use 1D view for scalars to reduce GPU memory operations

  d_scalars = t_int_2("surf_collide_diffuse:scalars");
  d_nsingle = Kokkos::subview(d_scalars,0);
  d_nreact_one = Kokkos::subview(d_scalars,1);

  h_scalars = t_host_int_2("surf_collide_diffuse:scalars_mirror");
  h_nsingle = Kokkos::subview(h_scalars,0);
  h_nreact_one = Kokkos::subview(h_scalars,1);
}

SurfCollideDiffuseKokkos::SurfCollideDiffuseKokkos(SPARTA *sparta) :
  SurfCollideDiffuse(sparta),
  fix_ambi_kk_copy(sparta),
  fix_vibmode_kk_copy(sparta),
  sr_kk_global_copy{VAL_2(KKCopy<SurfReactGlobalKokkos>(sparta))},
  sr_kk_prob_copy{VAL_2(KKCopy<SurfReactProbKokkos>(sparta))},
  rand_pool(12345 // seed doesn't matter since it will just be copied over
#ifdef SPARTA_KOKKOS_EXACT
            , sparta
#endif
           )
{
  copy = 1;
}

/* ---------------------------------------------------------------------- */

SurfCollideDiffuseKokkos::~SurfCollideDiffuseKokkos()
{
  if (uncopy) {
    fix_ambi_kk_copy.uncopy();
    fix_vibmode_kk_copy.uncopy();

    for (int i = 0; i < KOKKOS_MAX_SURF_REACT_PER_TYPE; i++) {
      sr_kk_global_copy[i].uncopy();
      sr_kk_prob_copy[i].uncopy();
    }
  }

  if (copy) return;

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.destroy();
  if (random_backup)
    delete random_backup;
#endif
}

/* ---------------------------------------------------------------------- */

void SurfCollideDiffuseKokkos::init()
{
  SurfCollideDiffuse::init();

  ambi_flag = vibmode_flag = 0;
  if (modify->n_update_custom) {
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
}

/* ----------------------------------------------------------------------
   recalculate Tsurf values which are dynamic
   called by Update::setup() and Update::run()
---------------------------------------------------------------------- */

void SurfCollideDiffuseKokkos::dynamic()
{
  // VAREQUAL mode
  // equal-style variable sets single tsurf value for all surfs

  if (tmode == VAREQUAL) {

    // only evaluate variable if timestep is multiple of tfreq

    if (update->ntimestep % tfreq) return;
    tsurf = input->variable->compute_equal(tindex_var);
    if (tsurf <= 0.0) error->all(FLERR,"Surf_collide tsurf <= 0.0");

  // VARSURF mode
  // surf-style variable sets new tsurf values for all surfs
  // particle/surf collisions access t_persurf for local+ghost values

  } else if (tmode == VARSURF) {

    // only evaluate variable if timestep is multiple of tfreq

    int spreadflag = 0;
    if (update->ntimestep % tfreq == 0) {
      if (n_owned != surf->nown) {
        memory->destroy(t_owned);
        n_owned = surf->nown;
        memory->create(t_owned,n_owned,"surfcollide:t_owned");
      }

      input->variable->compute_surf(tindex_var,t_owned,1,0);
      spreadflag = 1;
    }

    // spread t_owned values to t_localghost values via spread_own2local()
    // if just re-computed variable OR surfs are
    //   distributed and load balance/adaptation took place on previous step

    if (spreadflag ||
        (surf->distributed && surf->localghost_changed_step == update->ntimestep-1)) {
      if (n_localghost != surf->nlocal + surf->nghost) {
        memory->destroy(t_localghost);
        n_localghost = surf->nlocal + surf->nghost;
        memory->create(t_localghost,n_localghost,"surfcollide:t_localghost");
      }

      surf->spread_own2local(1,DOUBLE,t_owned,t_localghost);
      t_persurf = t_localghost;

      auto h_t_persurf = HAT::t_float_1d(t_persurf,n_localghost);
      d_t_persurf = Kokkos::create_mirror_view_and_copy(SPADeviceType(),h_t_persurf);
    }

  // CUSTOM mode
  // ensure access to custom per-surf vec for tsurf values for all surfs
  // particle/surf collisions access t_persurf for local+ghost values

  } else if (tmode == CUSTOM) {
    SurfKokkos* surf_kk = (SurfKokkos*) surf;
    auto h_edvec_local = surf_kk->k_edvec_local.h_view;

    // spread owned values to local+ghost values via spread_custom()
    // estatus == 1 means owned values already spread to local+ghost values
    // if estatus == 0: owned values are new OR
    //   surfs are distributed and load balance/adaptation took place

    if (surf->estatus[tindex_custom] == 0) surf->spread_custom(tindex_custom);

    h_edvec_local[tindex_custom].k_view.sync_device();
    d_t_persurf = h_edvec_local[tindex_custom].k_view.d_view;
  }
}

/* ---------------------------------------------------------------------- */

void SurfCollideDiffuseKokkos::pre_collide()
{
  if (ambi_flag) {
    afix_kk->pre_update_custom_kokkos();
    fix_ambi_kk_copy.copy(afix_kk);
  }

  if (vibmode_flag) {
    vfix_kk->pre_update_custom_kokkos();
    fix_vibmode_kk_copy.copy(vfix_kk);
  }

  if (surf->nsr > KOKKOS_MAX_TOT_SURF_REACT)
    error->all(FLERR,"Kokkos currently supports two instances of each surface reaction method");

  if (surf->nsr > 0) {
    int nglob,nprob;
    nglob = nprob = 0;
    for (int n = 0; n < surf->nsr; n++) {
      if (!surf->sr[n]->kokkosable)
        error->all(FLERR,"Must use Kokkos-enabled surface reaction method with Kokkos");
      if (strcmp(surf->sr[n]->style,"global") == 0) {
        sr_kk_global_copy[nglob].copy((SurfReactGlobalKokkos*)(surf->sr[n]));
        sr_kk_global_copy[nglob].obj.pre_react();
        sr_type_list[n] = 0;
        sr_map[n] = nprob;
        nglob++;
      } else if (strcmp(surf->sr[n]->style,"prob") == 0) {
        sr_kk_prob_copy[nprob].copy((SurfReactProbKokkos*)(surf->sr[n]));
        sr_kk_prob_copy[nprob].obj.pre_react();
        sr_type_list[n] = 1;
        sr_map[n] = nprob;
        nprob++;
      } else {
        error->all(FLERR,"Unknown Kokkos surface reaction method");
      }
    }

    if (nglob > KOKKOS_MAX_SURF_REACT_PER_TYPE || nprob > KOKKOS_MAX_SURF_REACT_PER_TYPE)
      error->all(FLERR,"Kokkos currently supports two instances of each surface reaction method");
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
  d_particles = particle_kk->k_particles.d_view;
  d_species = particle_kk->k_species.d_view;
  boltz = update->boltz;

  rotstyle = NONE;
  if (Pointers::collide) rotstyle = Pointers::collide->rotstyle;
  vibstyle = NONE;
  if (Pointers::collide) vibstyle = Pointers::collide->vibstyle;

  Kokkos::deep_copy(d_scalars,0);
}

/* ---------------------------------------------------------------------- */

void SurfCollideDiffuseKokkos::post_collide()
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  if (ambi_flag || vibmode_flag) particle_kk->modify(Device,CUSTOM_MASK);

  Kokkos::deep_copy(h_scalars,d_scalars);

  int m = surf->find_collide(id);
  auto sc = surf->sc[m]; // can't modify the copy directly, use the original
  sc->nsingle += h_nsingle();
  surf->nreact_one += h_nreact_one();

  d_particles = {};
}

/* ---------------------------------------------------------------------- */

void SurfCollideDiffuseKokkos::backup()
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  d_particles = particle_kk->k_particles.d_view;

  if (surf->nsr > 0) {
    int nglob,nprob;
    nglob = nprob = 0;
    for (int n = 0; n < surf->nsr; n++) {
      if (strcmp(surf->sr[n]->style,"global") == 0) {
        sr_kk_global_copy[nglob].obj.backup();
        nglob++;
      } else if (strcmp(surf->sr[n]->style,"prob") == 0) {
        sr_kk_prob_copy[nprob].obj.backup();
        nprob++;
      }
    }
  }

#ifdef SPARTA_KOKKOS_EXACT
  if (!random_backup)
    random_backup = new RanKnuth(12345 + comm->me);
  memcpy(random_backup,random,sizeof(RanKnuth));
#endif
}

/* ---------------------------------------------------------------------- */

void SurfCollideDiffuseKokkos::restore()
{
  if (surf->nsr > 0) {
    int nglob,nprob;
    nglob = nprob = 0;
    for (int n = 0; n < surf->nsr; n++) {
      if (strcmp(surf->sr[n]->style,"global") == 0) {
        sr_kk_global_copy[nglob].obj.restore();
        nglob++;
      } else if (strcmp(surf->sr[n]->style,"prob") == 0) {
        sr_kk_prob_copy[nprob].obj.restore();
        nprob++;
      }
    }
  }

  Kokkos::deep_copy(d_scalars,0);

#ifdef SPARTA_KOKKOS_EXACT
  memcpy(random,random_backup,sizeof(RanKnuth));
#endif
}

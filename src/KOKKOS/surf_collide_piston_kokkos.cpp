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

/* ----------------------------------------------------------------------
   Contributing author: Alan Stagg (Sandia)
------------------------------------------------------------------------- */

#include "surf_collide_piston_kokkos.h"
#include "fix.h"
#include "modify.h"
#include "error.h"
#include "sparta_masks.h"
#include "surf.h"

using namespace SPARTA_NS;

#define VAL_1(X) X
#define VAL_2(X) VAL_1(X), VAL_1(X)

/* ---------------------------------------------------------------------- */

SurfCollidePistonKokkos::SurfCollidePistonKokkos(SPARTA *sparta, int narg, char **arg) :
  SurfCollidePiston(sparta, narg, arg),
  fix_ambi_kk_copy(sparta),
  fix_vibmode_kk_copy(sparta),
  sr_kk_global_copy{VAL_2(KKCopy<SurfReactGlobalKokkos>(sparta))},
  sr_kk_prob_copy{VAL_2(KKCopy<SurfReactProbKokkos>(sparta))}
{
  kokkosable = 1;

  // use 1D view for scalars to reduce GPU memory operations

  d_scalars = t_int_2("surf_collide_piston:scalars");
  d_nsingle = Kokkos::subview(d_scalars,0);
  d_nreact_one = Kokkos::subview(d_scalars,1);

  h_scalars = t_host_int_2("surf_collide_piston:scalars_mirror");
  h_nsingle = Kokkos::subview(h_scalars,0);
  h_nreact_one = Kokkos::subview(h_scalars,1);
}

/* ---------------------------------------------------------------------- */

SurfCollidePistonKokkos::SurfCollidePistonKokkos(SPARTA *sparta) :
  SurfCollidePiston(sparta),
  fix_ambi_kk_copy(sparta),
  fix_vibmode_kk_copy(sparta),
  sr_kk_global_copy{VAL_2(KKCopy<SurfReactGlobalKokkos>(sparta))},
  sr_kk_prob_copy{VAL_2(KKCopy<SurfReactProbKokkos>(sparta))}
{
  id = NULL;
  style = NULL;
}

/* ---------------------------------------------------------------------- */

SurfCollidePistonKokkos::~SurfCollidePistonKokkos()
{
  if (copy) return;

  fix_ambi_kk_copy.uncopy(1);
  fix_vibmode_kk_copy.uncopy(1);

  for (int i = 0; i < KOKKOS_MAX_SURF_REACT_PER_TYPE; i++) {
    sr_kk_global_copy[i].uncopy();
    sr_kk_prob_copy[i].uncopy();
  }
}

/* ---------------------------------------------------------------------- */

void SurfCollidePistonKokkos::init()
{
  SurfCollidePiston::init();

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

/* ---------------------------------------------------------------------- */

void SurfCollidePistonKokkos::pre_collide()
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

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,PARTICLE_MASK|SPECIES_MASK);
  d_particles = particle_kk->k_particles.d_view;

  Kokkos::deep_copy(d_scalars,0);
}

/* ---------------------------------------------------------------------- */

void SurfCollidePistonKokkos::post_collide()
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  if (ambi_flag || vibmode_flag) particle_kk->modify(Device,CUSTOM_MASK);

  Kokkos::deep_copy(h_scalars,d_scalars);

  int m = surf->find_collide(id);
  auto sc = surf->sc[m]; // can't modify the copy directly, use the original
  sc->nsingle += h_nsingle();
  surf->nreact_one += h_nreact_one();
}

/* ---------------------------------------------------------------------- */

void SurfCollidePistonKokkos::backup()
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
}

/* ---------------------------------------------------------------------- */

void SurfCollidePistonKokkos::restore()
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
}

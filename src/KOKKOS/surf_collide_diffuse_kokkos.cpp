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
#include "stdlib.h"
#include "string.h"
#include "surf_collide_diffuse_kokkos.h"
#include "surf_kokkos.h"
#include "surf_react.h"
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

enum{INT,DOUBLE};                      // several files

/* ---------------------------------------------------------------------- */

SurfCollideDiffuseKokkos::SurfCollideDiffuseKokkos(SPARTA *sparta, int narg, char **arg) :
  SurfCollideDiffuse(sparta, narg, arg),
  fix_ambi_kk_copy(sparta),
  fix_vibmode_kk_copy(sparta),
  rand_pool(12345 + comm->me
#ifdef SPARTA_KOKKOS_EXACT
            , sparta
#endif
           )
{
#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.init(random);
#endif

  k_nsingle = DAT::tdual_int_scalar("SurfCollide:nsingle");
  d_nsingle = k_nsingle.d_view;
  h_nsingle = k_nsingle.h_view;
}

SurfCollideDiffuseKokkos::SurfCollideDiffuseKokkos(SPARTA *sparta) :
  SurfCollideDiffuse(sparta),
  fix_ambi_kk_copy(sparta),
  fix_vibmode_kk_copy(sparta),
  rand_pool(12345 // seed doesn't matter since it will just be copied over
#ifdef SPARTA_KOKKOS_EXACT
            , sparta
#endif
           )
{
  tstr = NULL;
  random = NULL;
  id = NULL;
  style = NULL;
}

/* ---------------------------------------------------------------------- */

SurfCollideDiffuseKokkos::~SurfCollideDiffuseKokkos()
{
  if (copy) return;

  fix_ambi_kk_copy.uncopy(1);
  fix_vibmode_kk_copy.uncopy(1);

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.destroy();
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

  if (tmode == CUSTOM) {
    SurfKokkos* surf_kk = (SurfKokkos*) surf;
    surf_kk->sync(Device,SURF_CUSTOM_MASK);

    int tindex = surf->find_custom(tstr);
    auto h_ewhich = surf_kk->k_ewhich.h_view;
    auto h_edvec = surf_kk->k_edvec.h_view;
    d_tvector = h_edvec[h_ewhich[tindex]].k_view.d_view;
  }

  rotstyle = NONE;
  if (Pointers::collide) rotstyle = Pointers::collide->rotstyle;
  vibstyle = NONE;
  if (Pointers::collide) vibstyle = Pointers::collide->vibstyle;
}

/* ---------------------------------------------------------------------- */

void SurfCollideDiffuseKokkos::post_collide()
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  if (ambi_flag || vibmode_flag) particle_kk->modify(Device,CUSTOM_MASK);
  d_particles = decltype(d_particles)();
}

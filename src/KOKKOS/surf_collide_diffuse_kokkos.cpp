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
#include "stdlib.h"
#include "string.h"
#include "surf_collide_diffuse_kokkos.h"
#include "surf.h"
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

/* ---------------------------------------------------------------------- */

SurfCollideDiffuseKokkos::SurfCollideDiffuseKokkos(SPARTA *sparta, int narg, char **arg) :
  SurfCollideDiffuse(sparta, narg, arg),
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
  rand_pool(12345 // seed doesn't matter since it will just be copied over
#ifdef SPARTA_KOKKOS_EXACT
            , sparta
#endif
           )
{
  // ID and style
  // ID must be all alphanumeric chars or underscores

  int narg = 4;
  const char* arg[] = {"sc_kk_diffuse_copy","diffuse","300.0","1.0"};

  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  for (int i = 0; i < n-1; i++)
    if (!isalnum(id[i]) && id[i] != '_')
      error->all(FLERR,"Surf_collide ID must be alphanumeric or "
                 "underscore characters");

  dynamicflag = 1;
  n = strlen(arg[1]) + 1;
  style = new char[n];
  strcpy(style,arg[1]);

  vector_flag = 1;
  size_vector = 2;

  nsingle = ntotal = 0;

  copy = 0;

  if (narg < 4) error->all(FLERR,"Illegal surf_collide diffuse command");

  allowreact = 1;

  tstr = NULL;

  if (strstr(arg[2],"v_") == arg[2]) {
    int n = strlen(&arg[2][2]) + 1;
    tstr = new char[n];
    strcpy(tstr,&arg[2][2]);
  } else {
    twall = input->numeric(FLERR,(char*)arg[2]);
    if (twall <= 0.0) error->all(FLERR,"Surf_collide diffuse temp <= 0.0");
  }

  acc = input->numeric(FLERR,arg[3]);
  if (acc < 0.0 || acc > 1.0)
    error->all(FLERR,"Illegal surf_collide diffuse command");

  // optional args

  tflag = rflag = 0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"translate") == 0) {
      if (iarg+4 > narg)
        error->all(FLERR,"Illegal surf_collide diffuse command");
      tflag = 1;
      vx = atof(arg[iarg+1]);
      vy = atof(arg[iarg+2]);
      vz = atof(arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"rotate") == 0) {
      if (iarg+7 > narg)
        error->all(FLERR,"Illegal surf_collide diffuse command");
      rflag = 1;
      px = atof(arg[iarg+1]);
      py = atof(arg[iarg+2]);
      pz = atof(arg[iarg+3]);
      wx = atof(arg[iarg+4]);
      wy = atof(arg[iarg+5]);
      wz = atof(arg[iarg+6]);
      if (domain->dimension == 2 && pz != 0.0)
        error->all(FLERR,"Surf_collide diffuse rotation invalid for 2d");
      if (domain->dimension == 2 && (wx != 0.0 || wy != 0.0))
        error->all(FLERR,"Surf_collide diffuse rotation invalid for 2d");
      iarg += 7;
    } else error->all(FLERR,"Illegal surf_collide diffuse command");
  }

  if (tflag && rflag) error->all(FLERR,"Illegal surf_collide diffuse command");
  if (tflag || rflag) trflag = 1;
  else trflag = 0;

  random = NULL;

  k_nsingle = DAT::tdual_int_scalar("SurfCollide:nsingle");
  d_nsingle = k_nsingle.d_view;
  h_nsingle = k_nsingle.h_view;

  allowreact = 0;
}

/* ---------------------------------------------------------------------- */

SurfCollideDiffuseKokkos::~SurfCollideDiffuseKokkos()
{
  if (copy) return;

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.destroy();
#endif
}

/* ---------------------------------------------------------------------- */

void SurfCollideDiffuseKokkos::pre_collide()
{
  if (modify->n_update_custom)
    error->all(FLERR,"Cannot yet use surf_collide diffuse/kk"
               "with fix vibmode or fix ambipolar");

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
  particle_kk->sync(Device,SPECIES_MASK);
  d_species = particle_kk->k_species.d_view;
  boltz = update->boltz;

  rotstyle = NONE;
  if (Pointers::collide) rotstyle = Pointers::collide->rotstyle;
  vibstyle = NONE;
  if (Pointers::collide) vibstyle = Pointers::collide->vibstyle;
}

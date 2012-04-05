/* ----------------------------------------------------------------------
   DSMC - Sandia parallel DSMC code
   www.sandia.gov/~sjplimp/dsmc.html
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2011) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level DSMC directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "surf_collide_diffuse.h"
#include "math_extra.h"
#include "input.h"
#include "variable.h"
#include "update.h"
#include "comm.h"
#include "random_mars.h"
#include "random_park.h"
#include "math_const.h"
#include "error.h"

using namespace DSMC_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

SurfCollideDiffuse::SurfCollideDiffuse(DSMC *dsmc, int narg, char **arg) :
  SurfCollide(dsmc, narg, arg)
{
  if (narg != 4) error->all(FLERR,"Illegal surf_collide diffuse command");

  tstr = NULL;

  acccoeff = atof(arg[2]);
  if (acccoeff < 0.0 || acccoeff > 1.0) 
    error->all(FLERR,"Illegal surf_collide diffuse command");

  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    tstr = new char[n];
    strcpy(tstr,&arg[3][2]);
  } else {
    twall = atof(arg[3]);
    if (twall < 0.0) error->all(FLERR,"Illegal surf_collide diffuse command");
  }

  random = NULL;
}

/* ---------------------------------------------------------------------- */

SurfCollideDiffuse::~SurfCollideDiffuse()
{
  delete [] tstr;
  delete random;
}

/* ---------------------------------------------------------------------- */

void SurfCollideDiffuse::init()
{
  // initialize RNG

  if (random == NULL) {
    random = new RanPark(update->ranmaster->uniform());
    double seed = update->ranmaster->uniform();
    random->reset(seed,comm->me,100);
  }

  // check variable

  if (tstr) {
    tvar = input->variable->find(tstr);
    if (tvar < 0) 
      error->all(FLERR,"Surf_collide diffuse variable name does not exist");
    if (!input->variable->equal_style(tvar))
      error->all(FLERR,"Surf_collide diffuse variable is invalid style");
  }
}

/* ----------------------------------------------------------------------
   set current surface temperature
------------------------------------------------------------------------- */

void SurfCollideDiffuse::dynamic()
{
  twall = input->variable->compute_equal(tvar);
}

/* ----------------------------------------------------------------------
   particle collision with surface
   p = particle with current x = collision pt, current v = incident v
   norm = surface normal unit vector
------------------------------------------------------------------------- */

void SurfCollideDiffuse::collide(Particle::OnePart *p, double *norm)
{
  // specular reflection
  // reflect incident v around norm

  if (random->uniform() > acccoeff) {
    MathExtra::reflect3(p->v,norm);

  // diffuse reflection
  // vrm = most probable speed of species isp, eqns (4.1) and (4.7)
  // vparallel = velocity component along surface normal, eqn (12.3)
  // vperp12 = 2 velocity components perpendicular to surface normal

  } else {
    Particle::Species *species = particle->species;
    double vrm = sqrt(2.0*update->boltz*twall/species[p->ispecies].mass);
    double vparallel = vrm * sqrt(-log(random->uniform()));
    double theta = MY_2PI * random->uniform();
    double vperp1 = vrm*sqrt(-log(random->uniform()))*sin(theta);
    theta = MY_2PI * random->uniform();
    double vperp2 = vrm*sqrt(-log(random->uniform()))*sin(theta);

    /*
      erot(isp);
      evib(isp);
    */
  }
}

void SurfCollideDiffuse2::collide(Particle::OnePart *p, double *norm)
{
  // specular reflection
  // reflect incident v around norm

  if (random->uniform() > acccoeff) {
    MathExtra::reflect3(p->v,norm);

  // diffuse reflection
  // vrm = most probable speed of species isp, eqns (4.1) and (4.7)
  // vparallel = velocity component along surface normal, eqn (12.3)
  // vperp12 = 2 velocity components perpendicular to surface normal

  } else {
    Particle::Species *species = particle->species;
    double vrm = sqrt(2.0*update->boltz*twall/species[p->ispecies].mass);
    double vparallel = vrm * sqrt(-log(random->uniform()));
    double theta = MY_2PI * random->uniform();
    double up = vrm*sqrt(-log(random->uniform()))
    double vperp1 = up * sin(theta);
    double vperp2 = up * cos(theta);

    v[0] = vparallel*norm[0] + vperp1*AAA[0] + vperp2*BBB[0];
    v[1] = vparallel*norm[1] + vperp1*AAA[1] + vperp2*BBB[1];
    v[2] = vparallel*norm[2] + vperp1*AAA[2] + vperp2*BBB[2];

    /*
      erot(isp);
      evib(isp);
    */
  }
}


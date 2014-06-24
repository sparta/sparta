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
#include "surf_collide_diffuse.h"
#include "math_extra.h"
#include "input.h"
#include "variable.h"
#include "particle.h"
#include "domain.h"
#include "update.h"
#include "comm.h"
#include "random_mars.h"
#include "random_park.h"
#include "math_const.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

SurfCollideDiffuse::SurfCollideDiffuse(SPARTA *sparta, int narg, char **arg) :
  SurfCollide(sparta, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal surf_collide diffuse command");

  tstr = NULL;

  if (strstr(arg[2],"v_") == arg[2]) {
    int n = strlen(&arg[2][2]) + 1;
    tstr = new char[n];
    strcpy(tstr,&arg[2][2]);
  } else {
    twall = atof(arg[2]);
    if (twall < 0.0) error->all(FLERR,"Illegal surf_collide diffuse command");
  }

  acc = atof(arg[3]);
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
   particle collision with surface
   p = particle with current x = collision pt, current v = incident v
   norm = surface normal unit vector
   reset p->v to post-collision outward velocity
     with optional translation and rotation
   reset erot and evib of particle to random new values
------------------------------------------------------------------------- */

void SurfCollideDiffuse::collide(Particle::OnePart *p, double *norm)
{
  // specular reflection
  // reflect incident v around norm

  if (random->uniform() > acc) {
    MathExtra::reflect3(p->v,norm);

  // diffuse reflection
  // vrm = most probable speed of species, eqns (4.1) and (4.7)
  // vperp = velocity component perpendicular to surface along norm, eqn (12.3)
  // vtan12 = 2 velocity components tangential to surface
  // tangent1 = component of particle v tangential to surface,
  //   check if tangent1 = 0 (normal collision), set randomly
  // tangent2 = norm x tangent1 = orthogonal tangential direction
  // tangent12 are both unit vectors

  } else {
    double tangent1[3],tangent2[3];
    Particle::Species *species = particle->species;
    int ispecies = p->ispecies;

    double vrm = sqrt(2.0*update->boltz*twall/species[ispecies].mass);
    double vperp = vrm * sqrt(-log(random->uniform()));

    double theta = MY_2PI * random->uniform();
    double vtangent = vrm * sqrt(-log(random->uniform()));
    double vtan1 = vtangent * sin(theta);
    double vtan2 = vtangent * cos(theta);

    double *v = p->v;
    double dot = MathExtra::dot3(v,norm);
    tangent1[0] = v[0] - dot*norm[0];
    tangent1[1] = v[1] - dot*norm[1];
    tangent1[2] = v[2] - dot*norm[2];
    if (MathExtra::lensq3(tangent1) == 0.0) {
      tangent2[0] = random->uniform();
      tangent2[1] = random->uniform();
      tangent2[2] = random->uniform();
      MathExtra::cross3(norm,tangent2,tangent1);
    }

    MathExtra::norm3(tangent1);
    MathExtra::cross3(norm,tangent1,tangent2);

    // add in translation or rotation term if specified

    if (trflag) {
      double vxdelta,vydelta,vzdelta;
      if (tflag) {
        vxdelta = vx; vydelta = vy; vzdelta = vz;
      } else {
        double *x = p->x;
        vxdelta = wy*(x[2]-pz) - wz*(x[1]-py);
        vydelta = wz*(x[0]-px) - wx*(x[2]-pz);
        vzdelta = wx*(x[1]-py) - wy*(x[0]-px);
      }
      v[0] = vperp*norm[0] + vtan1*tangent1[0] + vtan2*tangent2[0] + vxdelta;
      v[1] = vperp*norm[1] + vtan1*tangent1[1] + vtan2*tangent2[1] + vydelta;
      v[2] = vperp*norm[2] + vtan1*tangent1[2] + vtan2*tangent2[2] + vzdelta;

    // no translation or rotation

    } else {
      v[0] = vperp*norm[0] + vtan1*tangent1[0] + vtan2*tangent2[0];
      v[1] = vperp*norm[1] + vtan1*tangent1[1] + vtan2*tangent2[1];
      v[2] = vperp*norm[2] + vtan1*tangent1[2] + vtan2*tangent2[2];
    }

    p->erot = particle->erot(ispecies,random);
    p->evib = particle->evib(ispecies,random);
  }
}

/* ----------------------------------------------------------------------
   set current surface temperature
------------------------------------------------------------------------- */

void SurfCollideDiffuse::dynamic()
{
  twall = input->variable->compute_equal(tvar);
}

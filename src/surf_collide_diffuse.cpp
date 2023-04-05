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
#include "surf_collide_diffuse.h"
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

using namespace SPARTA_NS;
using namespace MathConst;

enum{INT,DOUBLE};                      // several files
enum{NUMERIC,VARIABLE,CUSTOM};

/* ---------------------------------------------------------------------- */

SurfCollideDiffuse::SurfCollideDiffuse(SPARTA *sparta, int narg, char **arg) :
  SurfCollide(sparta, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal surf_collide diffuse command");

  tstr = NULL;

  if (strstr(arg[2],"v_") == arg[2]) {
    dynamicflag = 1;
    tmode = VARIABLE;
    int n = strlen(&arg[2][2]) + 1;
    tstr = new char[n];
    strcpy(tstr,&arg[2][2]);
  } else if (strstr(arg[2],"s_") == arg[2]) {
    tmode = CUSTOM;
    int n = strlen(&arg[2][2]) + 1;
    tstr = new char[n];
    strcpy(tstr,&arg[2][2]);
  } else {
    tmode = NUMERIC;
    twall = input->numeric(FLERR,arg[2]);
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

      if (domain->dimension == 2) {
        if (pz != 0.0)
          error->all(FLERR,"Surf_collide diffuse rotation invalid for 2d");
        if (!domain->axisymmetric && (wx != 0.0 || wy != 0.0))
          error->all(FLERR,"Surf_collide diffuse rotation invalid for 2d");
        if (domain->axisymmetric && (wy != 0.0 || wz != 0.0))
          error->all(FLERR,
                     "Surf_collide diffuse rotation invalid for 2d axisymmetric");
      }

      iarg += 7;
    } else error->all(FLERR,"Illegal surf_collide diffuse command");
  }

  if (tflag && rflag) error->all(FLERR,"Illegal surf_collide diffuse command");
  if (tflag || rflag) trflag = 1;
  else trflag = 0;

  vstream[0] = vstream[1] = vstream[2] = 0.0;

  // initialize RNG

  random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,comm->me,100);
}

/* ---------------------------------------------------------------------- */

SurfCollideDiffuse::~SurfCollideDiffuse()
{
  if (copy) return;

  delete [] tstr;
  delete random;
}

/* ---------------------------------------------------------------------- */

void SurfCollideDiffuse::init()
{
  SurfCollide::init();

  // check variable and custom surf vector

  if (tmode == VARIABLE) {
    tvar = input->variable->find(tstr);
    if (tvar < 0)
      error->all(FLERR,"Surf_collide diffuse variable name does not exist");
    if (!input->variable->equal_style(tvar))
      error->all(FLERR,"Surf_collide diffuse variable is invalid style");
  } else if (tmode == CUSTOM) {
    int tindex = surf->find_custom(tstr);
    if (tindex < 0)
      error->all(FLERR,"Surf_collide diffuse could not find "
                 "custom per-surf vector");
    if (surf->etype[tindex] != DOUBLE || surf->esize[tindex] != 0)
      error->all(FLERR,"Surf_collide diffuse custom per-surf vector in invalid");
    tvector = surf->edvec[surf->ewhich[tindex]];
  }
}

/* ----------------------------------------------------------------------
   particle collision with surface with optional chemistry
   ip = particle with current x = collision pt, current v = incident v
   isurf = index of surface element
   norm = surface normal unit vector
   isr = index of reaction model if >= 0, -1 for no chemistry
   ip = reset to NULL if destroyed by chemistry
   return jp = new particle if created by chemistry
   return reaction = index of reaction (1 to N) that took place, 0 = no reaction
   resets particle(s) to post-collision outward velocity
------------------------------------------------------------------------- */

Particle::OnePart *SurfCollideDiffuse::
collide(Particle::OnePart *&ip, double &,
        int isurf, double *norm, int isr, int &reaction)
{
  nsingle++;

  // if surface chemistry defined, attempt reaction
  // reaction = 1 to N for which reaction took place, 0 for none
  // velreset = 1 if reaction reset post-collision velocity, else 0

  Particle::OnePart iorig;
  Particle::OnePart *jp = NULL;
  reaction = 0;
  int velreset = 0;

  if (isr >= 0) {
    if (modify->n_surf_react) memcpy(&iorig,ip,sizeof(Particle::OnePart));
    reaction = surf->sr[isr]->react(ip,isurf,norm,jp,velreset);
    if (reaction) surf->nreact_one++;
  }

  // diffuse reflection for each particle
  // only if SurfReact did not already reset velocities
  // also both partiticles need to trigger any fixes
  //   to update per-particle properties which depend on
  //   temperature of the particle, e.g. fix vibmode and fix ambipolar

  if (tmode == CUSTOM) twall = tvector[isurf];

  if (ip) {
    if (!velreset) diffuse(ip,norm);
    if (modify->n_update_custom) {
      int i = ip - particle->particles;
      modify->update_custom(i,twall,twall,twall,vstream);
    }
  }
  if (jp) {
    if (!velreset) diffuse(jp,norm);
    if (modify->n_update_custom) {
      int j = jp - particle->particles;
      modify->update_custom(j,twall,twall,twall,vstream);
    }
  }

  // call any fixes with a surf_react() method
  // they may reset j to -1, e.g. fix ambipolar
  //   in which case newly created j is deleted

  if (reaction && modify->n_surf_react) {
    int i = -1;
    if (ip) i = ip - particle->particles;
    int j = -1;
    if (jp) j = jp - particle->particles;
    modify->surf_react(&iorig,i,j);
    if (jp && j < 0) {
      jp = NULL;
      particle->nlocal--;
    }
  }

  return jp;
}

/* ----------------------------------------------------------------------
   diffusive particle collision with surface
   p = particle with current x = collision pt, current v = incident v
   norm = surface normal unit vector
   resets particle(s) to post-collision outward velocity
------------------------------------------------------------------------- */

void SurfCollideDiffuse::diffuse(Particle::OnePart *p, double *norm)
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

    double vrm = sqrt(2.0*update->boltz * twall / species[ispecies].mass);
    double vperp = vrm * sqrt(-log(random->uniform()));

    double theta = MY_2PI * random->uniform();
    double vtangent = vrm * sqrt(-log(random->uniform()));
    double vtan1 = vtangent * sin(theta);
    double vtan2 = vtangent * cos(theta);

    double *v = p->v;
    double dot = MathExtra::dot3(v,norm);

    double beta_un,normalized_distbn_fn;

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

    // add in translation or rotation vector if specified
    // only keep portion of vector tangential to surface element

    if (trflag) {
      double vxdelta,vydelta,vzdelta;
      if (tflag) {
        vxdelta = vx; vydelta = vy; vzdelta = vz;
        double dot = vxdelta*norm[0] + vydelta*norm[1] + vzdelta*norm[2];

        if (fabs(dot) > 0.001) {
          dot /= vrm;
          do {
            do {
              beta_un = (6.0*random->uniform() - 3.0);
            } while (beta_un + dot < 0.0);
            normalized_distbn_fn = 2.0 * (beta_un + dot) /
              (dot + sqrt(dot*dot + 2.0)) *
              exp(0.5 + (0.5*dot)*(dot-sqrt(dot*dot + 2.0)) -
                  beta_un*beta_un);
          } while (normalized_distbn_fn < random->uniform());
          vperp = beta_un*vrm;
        }

      } else {
        double *x = p->x;
        vxdelta = wy*(x[2]-pz) - wz*(x[1]-py);
        vydelta = wz*(x[0]-px) - wx*(x[2]-pz);
        vzdelta = wx*(x[1]-py) - wy*(x[0]-px);
        double dot = vxdelta*norm[0] + vydelta*norm[1] + vzdelta*norm[2];
        vxdelta -= dot*norm[0];
        vydelta -= dot*norm[1];
        vzdelta -= dot*norm[2];
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

    // initialize rot/vib energy

    p->erot = particle->erot(ispecies,twall,random);
    p->evib = particle->evib(ispecies,twall,random);
  }
}

/* ----------------------------------------------------------------------
   wrapper on diffuse() method to perform collision for a single particle
   pass in 2 coefficients to match command-line args for style diffuse
   flags, coeffs can be NULL
   called by SurfReactAdsorb
------------------------------------------------------------------------- */

void SurfCollideDiffuse::wrapper(Particle::OnePart *p, double *norm,
                                 int *flags, double *coeffs)
{
  if (coeffs) {
    twall = coeffs[0];
    acc = coeffs[1];
  }

  diffuse(p,norm);
}

/* ----------------------------------------------------------------------
   return flags and coeffs for this SurfCollide instance to caller
------------------------------------------------------------------------- */

void SurfCollideDiffuse::flags_and_coeffs(int *flags, double *coeffs)
{
  if (tmode == CUSTOM)
    error->all(FLERR,"Surf_collide diffuse with custom per-surf Twall "
               "does not support external caller");

  coeffs[0] = twall;
  coeffs[1] = acc;
}

/* ----------------------------------------------------------------------
   set current surface temperature from equal-style variable
------------------------------------------------------------------------- */

void SurfCollideDiffuse::dynamic()
{
  twall = input->variable->compute_equal(tvar);
  if (twall <= 0.0) error->all(FLERR,"Surf_collide diffuse temp <= 0.0");
}

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

/* ----------------------------------------------------------------------
   Contributing author: Krishnan Gopalan (NASA Ames)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "surf_collide_supra.h"
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
#include "random_park.h"
#include "math_const.h"
#include "math_extra.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

SurfCollideSupra::SurfCollideSupra(SPARTA *sparta, int narg, char **arg) :
  SurfCollide(sparta, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal surf_collide supra command");

  tstr = NULL;

  if (strstr(arg[2],"v_") == arg[2]) {
    int n = strlen(&arg[2][2]) + 1;
    tstr = new char[n];
    strcpy(tstr,&arg[2][2]);
  } else {
    twall = atof(arg[2]);
    if (twall < 0.0) error->all(FLERR,"Illegal surf_collide supra command");
  }

  e = atof(arg[3]);
  if (e < 0.0 || e > 1.0) error->all(FLERR,"Illegal surf_collide supra command");

  // optional args

  tflag = rflag = 0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"translate") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal surf_collide supra command");
      tflag = 1;
      vx = atof(arg[iarg+1]);
      vy = atof(arg[iarg+2]);
      vz = atof(arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"rotate") == 0) {
      if (iarg+7 > narg) error->all(FLERR,"Illegal surf_collide supra command");
      rflag = 1;
      px = atof(arg[iarg+1]);
      py = atof(arg[iarg+2]);
      pz = atof(arg[iarg+3]);
      wx = atof(arg[iarg+4]);
      wy = atof(arg[iarg+5]);
      wz = atof(arg[iarg+6]);
      if (domain->dimension == 2 && pz != 0.0) 
        error->all(FLERR,"Surf_collide supra rotation invalid for 2d");
      if (domain->dimension == 2 && (wx != 0.0 || wy != 0.0))
        error->all(FLERR,"Surf_collide supra rotation invalid for 2d");
      iarg += 7;
    } else error->all(FLERR,"Illegal surf_collide supra command");
  }

  if (tflag && rflag) error->all(FLERR,"Illegal surf_collide supra command");
  if (tflag || rflag) trflag = 1;
  else trflag = 0;

  vstream[0] = vstream[1] = vstream[2] = 0.0;

  // initialize RNG
  
  random = new RanPark(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,comm->me,100);
}

/* ---------------------------------------------------------------------- */

SurfCollideSupra::~SurfCollideSupra()
{
  delete [] tstr;
  delete random;
}

/* ---------------------------------------------------------------------- */

void SurfCollideSupra::init()
{
  SurfCollide::init();

  // check variable

  if (tstr) {
    tvar = input->variable->find(tstr);
    if (tvar < 0) 
      error->all(FLERR,"Surf_collide supra variable name does not exist");
    if (!input->variable->equal_style(tvar))
      error->all(FLERR,"Surf_collide supra variable is invalid style");
  }
}

/* ----------------------------------------------------------------------
   particle collision with surface with optional chemistry
   ip = particle with current x = collision pt, current v = incident v
   norm = surface normal unit vector
   ip = set to NULL if destroyed by chemsitry
   return jp = new particle if created by chemistry
   resets particle(s) to post-collision outward velocity
------------------------------------------------------------------------- */

Particle::OnePart *SurfCollideSupra::
collide(Particle::OnePart *&ip, double *norm, double &, int isr)
{
  nsingle++;

  // if surface chemistry defined, attempt reaction
  // reaction = 1 if reaction took place

  Particle::OnePart iorig;
  Particle::OnePart *jp = NULL;
  int reaction = 0, ad_react = 0;

  if (isr >= 0) {
    if (modify->n_surf_react) memcpy(&iorig,ip,sizeof(Particle::OnePart));
    reaction = surf->sr[isr]->react(ip,norm,jp);
    if (reaction) surf->nreact_one++;
  }

  // supra reflection for each particle
  // if new particle J created, also need to trigger any fixes

  if (ip) supra(ip,norm);
  if (jp) {
    supra(jp,norm);
    if (modify->n_add_particle) {
      int j = jp - particle->particles;
      modify->add_particle(j,twall,twall,twall,vstream);
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
   supra reflection
   vrm = most probable speed of species, eqns (4.1) and (4.7)
   vperp = velocity component perpendicular to surface along norm, eqn (12.3)
   vtan12 = 2 velocity components tangential to surface
   tangent1 = component of particle v tangential to surface,
     check if tangent1 = 0 (normal collision), set randomly
   tangent2 = norm x tangent1 = orthogonal tangential direction
   tangent12 are both unit vectors
------------------------------------------------------------------------- */

void SurfCollideSupra::supra(Particle::OnePart *p, double *norm)
{
  double tangent1[3],tangent2[3];
  Particle::Species *species = particle->species;
  int ispecies = p->ispecies;
  double beta_un,normalized_distbn_fn;
    
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
    
  double v_mag = MathExtra::lensq3(v);
  double tan1 = MathExtra::dot3(v,tangent1);
  double tan2 = MathExtra::dot3(v,tangent2);
    
  double theta_i = acos(dot/v_mag);
  double psi_i = acos(cos(theta_i)*cos(theta_i));
  double phi_i = atan2(tan2,tan1);
    
  double P = 0;
  double phi_f,psi_f,cos_beta;
    
  while (random->uniform() > P) {
    phi_f = MY_2PI *  random->uniform();
    psi_f = acos(random->uniform());
    cos_beta = cos(psi_f)*cos(psi_i) + sin(psi_f)*sin(psi_i)*cos(phi_f - phi_i);
    P = (1-e)/(1-e*cos_beta);
  } 
    
  double theta_f = acos(sqrt(cos(psi_f)));
  double vrm = sqrt(2.0*update->boltz * twall / species[ispecies].mass);
    
  // SUPRA model normal velocity

  double vperp = vrm * cos(theta_f);
    
  // SUPRA model tangential velocities

  double vtan1 = vrm * sin(theta_f) * cos(phi_f);
  double vtan2 = vrm * sin(theta_f) * sin(phi_f);

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
            beta_un = (6.0*random->gaussian() - 3.0);
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

  p->erot = particle->erot(ispecies,twall,random);
  p->evib = particle->evib(ispecies,twall,random);
}

/* ----------------------------------------------------------------------
   set current surface temperature
------------------------------------------------------------------------- */

void SurfCollideSupra::dynamic()
{
  twall = input->variable->compute_equal(tvar);
}

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
#include "surf_collide_td.h"
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

/* ---------------------------------------------------------------------- */

SurfCollideTD::SurfCollideTD(SPARTA *sparta, int narg, char **arg) :
  SurfCollide(sparta, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal surf_collide td command");

  tstr = NULL;

  if (strstr(arg[2],"v_") == arg[2]) {
    int n = strlen(&arg[2][2]) + 1;
    tstr = new char[n];
    strcpy(tstr,&arg[2][2]);
  } else {
    twall = atof(arg[2]);
    if (twall < 0.0) error->all(FLERR,"Illegal surf_collide td command");
  }

  // optional args

  barrier_flag = initen_flag = bond_flag = 0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"barrier") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal surf_collide td command");
      barrier_flag = 1;
      barrier_val = atof(arg[iarg+1]);
      if (barrier_val < 0.0)
        error->all(FLERR,"Illegal surf_collide td barrier value");
      iarg += 2;
    } else if (strcmp(arg[iarg],"initenergy") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal surf_collide td command");
      initen_flag = 1;
      initen_trans = atof(arg[iarg+1]);
      initen_rot = atof(arg[iarg+2]);
      initen_vib = atof(arg[iarg+3]);
      if (initen_trans < 0.0 || initen_rot < 0.0 || initen_vib < 0.0)
        error->all(FLERR,"Illegal surf_collide td initenergy value");
      iarg += 4;
    } else if (strcmp(arg[iarg],"bond") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal surf_collide td command");
      bond_flag = 1;
      bond_trans = atof(arg[iarg+1]);
      bond_rot = atof(arg[iarg+2]);
      bond_vib = atof(arg[iarg+3]);
      if (bond_trans < 0.0 || bond_rot < 0.0 || bond_vib < 0.0)
        error->all(FLERR,"Illegal surf_collide td bond value");
      iarg += 4;
    } else error->all(FLERR,"Illegal surf_collide td command");
  }

  vstream[0] = vstream[1] = vstream[2] = 0.0;

  // initialize RNG

  random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,comm->me,100);
}

/* ---------------------------------------------------------------------- */

SurfCollideTD::~SurfCollideTD()
{
  delete [] tstr;
  delete random;
}

/* ---------------------------------------------------------------------- */

void SurfCollideTD::init()
{
  SurfCollide::init();

  // check variable

  if (tstr) {
    tvar = input->variable->find(tstr);
    if (tvar < 0)
      error->all(FLERR,"Surf_collide td variable name does not exist");
    if (!input->variable->equal_style(tvar))
      error->all(FLERR,"Surf_collide td variable is invalid style");
  }
}

/* ----------------------------------------------------------------------
   particle collision with surface with optional chemistry
   ip = particle with current x = collision pt, current v = incident v
   norm = surface normal unit vector
   ip = set to NULL if destroyed by chemsitry
   return jp = new particle if created by chemistry
   return reaction = index of reaction (1 to N) that took place, 0 = no reaction
   resets particle(s) to post-collision outward velocity
------------------------------------------------------------------------- */

Particle::OnePart *SurfCollideTD::
collide(Particle::OnePart *&ip, double *norm, double &, int isr, int &reaction)
{
  nsingle++;

  // if surface chemistry defined, attempt reaction
  // reaction = 1 if reaction took place

  Particle::OnePart iorig;
  Particle::OnePart *jp = NULL;

  if (isr >= 0) {
    if (modify->n_surf_react) memcpy(&iorig,ip,sizeof(Particle::OnePart));
    reaction = surf->sr[isr]->react(ip,norm,jp);
    if (reaction) surf->nreact_one++;
  }

  // td reflection for each particle
  // particle I needs to trigger any fixes to update per-particle
  //  properties which depend on the temperature of the particle
  //  (e.g. fix vibmode and fix ambipolar)
  // if new particle J created, also need to trigger any fixes

  if (ip) {
    td(ip,norm);
    if (modify->n_update_custom) {
      int i = ip - particle->particles;
      modify->update_custom(i,twall,twall,twall,vstream);
    }
  }
  if (jp) {
    td(jp,norm);
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
   TD reflection
   vrm = most probable speed of species, eqns (4.1) and (4.7)
   vperp = velocity component perpendicular to surface along norm, eqn (12.3)
   vtan12 = 2 velocity components tangential to surface
   tangent1 = component of particle v tangential to surface,
     check if tangent1 = 0 (normal collision), set randomly
   tangent2 = norm x tangent1 = orthogonal tangential direction
   tangent12 are both unit vectors
------------------------------------------------------------------------- */

void SurfCollideTD::td(Particle::OnePart *p, double *norm)
{
  double tangent1[3],tangent2[3];
  Particle::Species *species = particle->species;
  int ispecies = p->ispecies;
  double boltz = update->boltz;

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

  double mass = species[ispecies].mass;
  double E_i = 0.5 * mass * MathExtra::lensq3(v);

  double E_t = boltz*twall;
  if (bond_flag) E_t += boltz*bond_trans;
  if (initen_flag) E_t += E_i*initen_trans;

  double E_n = E_t;
  if (barrier_flag) E_n += boltz*barrier_val;

  double vrm_n = sqrt(2.0*E_n / mass);
  double vrm_t = sqrt(2.0*E_t / mass);
  double vperp = vrm_n * sqrt(-log(random->uniform()));

  double theta = MY_2PI * random->uniform();
  double vtangent = vrm_t * sqrt(-log(random->uniform()));
  double vtan1 = vtangent * sin(theta);
  double vtan2 = vtangent * cos(theta);

  v[0] = vperp*norm[0] + vtan1*tangent1[0] + vtan2*tangent2[0];
  v[1] = vperp*norm[1] + vtan1*tangent1[1] + vtan2*tangent2[1];
  v[2] = vperp*norm[2] + vtan1*tangent1[2] + vtan2*tangent2[2];

  double twall_rot = twall;
  double twall_vib = twall;

  if (bond_flag) {
    twall_rot += bond_rot;
    twall_vib += bond_vib;
  }

  if (initen_flag) {
    twall_rot += E_i*initen_rot/boltz;
    twall_vib += E_i*initen_vib/boltz;
  }

  p->erot = particle->erot(ispecies,twall_rot,random);
  p->evib = particle->evib(ispecies,twall_vib,random);
}

/* ----------------------------------------------------------------------
   set current surface temperature
------------------------------------------------------------------------- */

void SurfCollideTD::dynamic()
{
  twall = input->variable->compute_equal(tvar);
}

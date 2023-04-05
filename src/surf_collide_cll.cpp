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
   Contributing author: Krishnan Gopalan (NASA Ames)
------------------------------------------------------------------------- */

#include "collide.h"
#include "mpi.h"
#include "stdio.h"
#include "sparta.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "surf_collide_cll.h"
#include "surf.h"
#include "surf_react.h"
#include "input.h"
#include "variable.h"
#include "particle.h"
#include "collide.h"
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

enum{NONE,DISCRETE,SMOOTH};
enum{NUMERIC,VARIABLE,CUSTOM};

/* ---------------------------------------------------------------------- */

SurfCollideCLL::SurfCollideCLL(SPARTA *sparta, int narg, char **arg) :
  SurfCollide(sparta, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal surf_collide cll command");

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
    if (twall <= 0.0) error->all(FLERR,"Surf_collide cll temp <= 0.0");
  }

  acc_n = atof(arg[3]);
  acc_t = atof(arg[4]);
  acc_rot = atof(arg[5]);
  acc_vib = atof(arg[6]);

  if (acc_n < 0.0 || acc_n > 1.0 || acc_t < 0.0 || acc_t > 1.0 ||
      acc_rot < 0.0 || acc_rot > 1.0 || acc_vib < 0.0 || acc_vib > 1.0)
    error->all(FLERR,"Surf_collide cll accommodation coeffs "
               "must be >= 0 and <= 1");

  // optional args

  eccen = 0.0;
  pflag = 0;
  tflag = rflag = 0;

  int iarg = 7;
  while (iarg < narg) {

    if (strcmp(arg[iarg],"partial") == 0) {
        if (iarg+2 > narg) error->all(FLERR,"Illegal surf_collide cll command");
        if (acc_n != acc_t)
          error->all(FLERR,"Surf_collide cll partial requires acc_n = acc_t");
        pflag = 1;
        eccen = atof(arg[iarg+1]);
        if (eccen < 0.0 || eccen >= 1.0 )
          error->all(FLERR,"Surf_collide cll eccentricity must be >= 0 and <= 1");
        iarg += 2;
    } else if (strcmp(arg[iarg],"translate") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal surf_collide cll command");
      tflag = 1;
      vx = atof(arg[iarg+1]);
      vy = atof(arg[iarg+2]);
      vz = atof(arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"rotate") == 0) {
      if (iarg+7 > narg) error->all(FLERR,"Illegal surf_collide cll command");
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
    } else error->all(FLERR,"Illegal surf_collide cll command");
  }

  if (tflag && rflag) error->all(FLERR,"Illegal surf_collide cll command");
  if (tflag || rflag) trflag = 1;
  else trflag = 0;

  vstream[0] = vstream[1] = vstream[2] = 0.0;

  // initialize RNG

  random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,comm->me,100);
}

/* ---------------------------------------------------------------------- */

SurfCollideCLL::~SurfCollideCLL()
{
  delete [] tstr;
  delete random;
}

/* ---------------------------------------------------------------------- */

void SurfCollideCLL::init()
{
  SurfCollide::init();

  // check variable

  if (tstr) {
    tvar = input->variable->find(tstr);
    if (tvar < 0)
      error->all(FLERR,"Surf_collide cll variable name does not exist");
    if (!input->variable->equal_style(tvar))
      error->all(FLERR,"Surf_collide cll variable is invalid style");
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

Particle::OnePart *SurfCollideCLL::
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

  // CLL reflection for each particle
  // only if SurfReact did not already reset velocities
  // also both partiticles need to trigger any fixes
  //   to update per-particle properties which depend on
  //   temperature of the particle, e.g. fix vibmode and fix ambipolar

  if (tmode == CUSTOM) twall = tvector[isurf];

  if (ip) {
    if (!velreset) cll(ip,norm);
    if (modify->n_update_custom) {
      int i = ip - particle->particles;
      modify->update_custom(i,twall,twall,twall,vstream);
    }
  }
  if (jp) {
    if (!velreset) cll(jp,norm);
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
  cll reflection
  vrm = most probable speed of species, eqns (4.1) and (4.7)
  vperp = velocity component perpendicular to surface along norm, eqn (12.3)
  vtan12 = 2 velocity components tangential to surface
  tangent1 = component of particle v tangential to surface,
  check if tangent1 = 0 (normal collision), set randomly
  tangent2 = norm x tangent1 = orthogonal tangential direction
  tangent12 are both unit vectors
------------------------------------------------------------------------- */

void SurfCollideCLL::cll(Particle::OnePart *p, double *norm)
{
  double tangent1[3],tangent2[3];
  Particle::Species *species = particle->species;
  int ispecies = p->ispecies;
  double beta_un,normalized_distbn_fn;

  double *v = p->v;
  double dot = MathExtra::dot3(v,norm);
  double vrm, vperp, vtan1, vtan2;

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

  double tan1 = MathExtra::dot3(v,tangent1);

  vrm = sqrt(2.0*update->boltz * twall / species[ispecies].mass);

  // CLL model normal velocity

  double r_1 = sqrt(-acc_n*log(random->uniform()));
  double theta_1 = MY_2PI * random->uniform();
  double dot_norm = dot/vrm * sqrt(1-acc_n);
  vperp = vrm * sqrt(r_1*r_1 + dot_norm*dot_norm + 2*r_1*dot_norm*cos(theta_1));

  // CLL model tangential velocities

  double r_2 = sqrt(-acc_t*log(random->uniform()));
  double theta_2 = MY_2PI * random->uniform();
  double vtangent = tan1/vrm * sqrt(1-acc_t);
  vtan1 = vrm * (vtangent + r_2*cos(theta_2));
  vtan2 = vrm * r_2 * sin(theta_2);

  // partial keyword
  // incomplete energy accommodation with partial/fully diffuse scattering
  // adjust the final angle of the particle while keeping
  //   the velocity magnitude or speed according to CLL scattering

  if (pflag) {
    double tan2 = MathExtra::dot3(v,tangent2);
    double theta_i, phi_i, psi_i, theta_f, phi_f, psi_f, cos_beta;

    theta_i = acos(dot/sqrt(MathExtra::lensq3(v)));
    psi_i = acos(dot*dot/MathExtra::lensq3(v));
    phi_i = atan2(tan2,tan1);

    double v_mag = sqrt(vperp*vperp + vtan1*vtan1 + vtan2*vtan2);

    double P = 0;
    while (random->uniform() > P) {
      phi_f = MY_2PI*random->uniform();
      psi_f = acos(1-random->uniform());
      cos_beta =  cos(psi_i)*cos(psi_f) +
        sin(psi_i)*sin(psi_f)*cos(phi_i - phi_f);
      P = (1-eccen)/(1-eccen*cos_beta);
    }

    theta_f = acos(sqrt(cos(psi_f)));

    vperp = v_mag * cos(theta_f);
    vtan1 = v_mag * sin(theta_f) * cos(phi_f);
    vtan2 = v_mag * sin(theta_f) * sin(phi_f);
  }

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
            exp(0.5 + (0.5*dot)*(dot-sqrt(dot*dot + 2.0)) - beta_un*beta_un);
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

  // rotational component

  if (!sparta->collide || sparta->collide->rotstyle == NONE ||
      species[ispecies].rotdof < 2) p->erot = 0.0;

  else {
    double erot_mag = sqrt(p->erot*(1-acc_rot)/(update->boltz*twall));

    double r_rot,cos_theta_rot,A_rot,X_rot;
    if (species[ispecies].rotdof == 2) {
      r_rot = sqrt(-acc_rot*log(random->uniform()));
      cos_theta_rot = cos(MY_2PI*random->uniform());
    }
    else if (species[ispecies].rotdof > 2) {
      A_rot = 0;
      while (A_rot < random->uniform()) {
        X_rot = 4*random->uniform();
        A_rot = 2.71828182845904523536028747*X_rot*X_rot*exp(-X_rot*X_rot);
      }
      r_rot = sqrt(acc_rot)*X_rot;
      cos_theta_rot = 2*random->uniform() - 1;
    }

    p->erot = update->boltz * twall *
      (r_rot*r_rot + erot_mag*erot_mag + 2*r_rot*erot_mag*cos_theta_rot);
    }

  // vibrational component
  // NOTE: check all references to species[]->vibtmp

  int vibdof = species[ispecies].vibdof;
  double r_vib, cos_theta_vib, A_vib, X_vib, evib_mag, evib_val;

  if (!sparta->collide || sparta->collide->vibstyle == NONE || vibdof < 2)
    p->evib = 0.0;

  else if (sparta->collide->vibstyle == DISCRETE && vibdof == 2) {
    double evib_star =
      -log(1 - random->uniform() *
           (1 - exp(-update->boltz*species[ispecies].vibtemp[0])));
    evib_val = p->evib + evib_star;
    evib_mag = sqrt(evib_val*(1-acc_vib)/(update->boltz*twall));
    r_vib = sqrt(-acc_vib*log(random->uniform()));
    cos_theta_vib = cos(MY_2PI*random->uniform());
    evib_val = update->boltz * twall *
      (r_vib*r_vib + evib_mag*evib_mag + 2*r_vib*evib_mag*cos_theta_vib);
    int ivib =  evib_val / (update->boltz*species[ispecies].vibtemp[0]);
    p->evib = ivib * update->boltz * species[ispecies].vibtemp[0];
  }

  else if (sparta->collide->vibstyle == SMOOTH || vibdof >= 2) {
    evib_mag = sqrt(p->evib*(1-acc_vib)/(update->boltz*twall));
    if (vibdof == 2) {
      r_vib = sqrt(-acc_vib*log(random->uniform()));
      cos_theta_vib = cos(MY_2PI*random->uniform());
    } else if (vibdof > 2) {
      A_vib = 0;
      while (A_vib < random->uniform()) {
        X_vib = 4*random->uniform();
        A_vib = 2.71828182845904523536028747*X_vib*X_vib*exp(-X_vib*X_vib);
      }
      r_vib = sqrt(acc_vib)*X_vib;
      cos_theta_vib = 2*random->uniform() - 1;
    }

    p->evib = update->boltz * twall *
      (r_vib*r_vib + evib_mag*evib_mag + 2*r_vib*evib_mag*cos_theta_vib);
  }
}

/* ----------------------------------------------------------------------
   wrapper on cll() method to perform collision for a single particle
   pass in flags/coefficients to match command-line args for style cll
   flags, coeffs can be NULL
   called by SurfReactAdsorb
------------------------------------------------------------------------- */

void SurfCollideCLL::wrapper(Particle::OnePart *p, double *norm,
                             int *flags, double *coeffs)
{
  if (flags) {
    twall = coeffs[0];
    acc_n = coeffs[1];
    acc_t = coeffs[2];
    acc_rot = coeffs[3];
    acc_vib = coeffs[4];

    if (flags[0]) eccen = coeffs[5];
    else eccen = 0.0;
  }

  cll(p,norm);
}

/* ----------------------------------------------------------------------
   return flags and coeffs for this SurfCollide instance to caller
------------------------------------------------------------------------- */

void SurfCollideCLL::flags_and_coeffs(int *flags, double *coeffs)
{
  if (tmode == CUSTOM)
    error->all(FLERR,"Surf_collide cll with custom per-surf Twall "
               "does not support external caller");

  coeffs[0] = twall;
  coeffs[1] = acc_n;
  coeffs[2] = acc_t;
  coeffs[3] = acc_rot;
  coeffs[4] = acc_vib;

  flags[0] = 0;
  if (eccen != 0.0) {
    flags[0] = 1;
    coeffs[5] = eccen;
  }
}

/* ----------------------------------------------------------------------
   set current surface temperature
------------------------------------------------------------------------- */

void SurfCollideCLL::dynamic()
{
  twall = input->variable->compute_equal(tvar);
  if (twall <= 0.0) error->all(FLERR,"Surf_collide cll temp <= 0.0");
}

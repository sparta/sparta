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
   Contributing author: Krishnan Swaminathan Gopalan (NASA Ames)
------------------------------------------------------------------------- */

#include "collide.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "surf_collide_impulsive.h"
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

enum{NONE,DISCRETE,SMOOTH};
enum{NUMERIC,VARIABLE,CUSTOM};

/* ---------------------------------------------------------------------- */

SurfCollideImpulsive::SurfCollideImpulsive(SPARTA *sparta, int narg, char **arg) :
  SurfCollide(sparta, narg, arg)
{
  if (narg < 10) error->all(FLERR,"Illegal surf_collide impulsive command");

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
    if (twall <= 0.0) error->all(FLERR,"Surf_collide impulsive temp <= 0.0");
  }

  softsphere_flag = 0;
  if (strcmp(arg[3],"softsphere") == 0) {
      eng_ratio = atof(arg[4]);
      eff_mass = atof(arg[5]);
      if (eng_ratio > 1.0)
        error->all(FLERR,"Illegal surf_collide impulsive energy ratio");
      if (eff_mass <= 0.0)
        error->all(FLERR,"Illegal surf_collide impulsive effective mass");
      softsphere_flag = 1;
    } else if (strcmp(arg[3],"tempvar") == 0) {
      u0_a = atof(arg[4]);
      u0_b = atof(arg[5]);
    } else error->all(FLERR,"Illegal surf_collide impulsive command");

  var_alpha = atof(arg[6]);
  theta_peak = atof(arg[7]);
  cos_theta_pow = atof(arg[8]);
  cos_phi_pow = atof(arg[9]);

  if (var_alpha < 0.0) error->all(FLERR,"Illegal surf_collide impulsive alpha");
  if (theta_peak < 0.0 || theta_peak > 90.0)
    error->all(FLERR,"Illegal surf_collide impulsive theta peak");
  if (cos_theta_pow < 0.0 || cos_phi_pow < 0.0)
    error->all(FLERR,"Illegal surf_collide impulsive cosine power");

  theta_peak *= MY_PI/180;
  var_alpha_sq = var_alpha * var_alpha;

  // optional args

  step_flag = double_flag = intenergy_flag = 0;
  step_size = 0;
  cos_theta_pow_2 = 0;

  int iarg = 10;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"step") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal surf_collide impulsive command");
      step_flag = 1;
      step_size = atof(arg[iarg+1]);
      if (step_size < 0.0)
        error->all(FLERR,"Illegal surf_collide impulsive step size");
      iarg += 2;
    } else if (strcmp(arg[iarg],"double") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal surf_collide impulsive command");
      double_flag = 1;
      cos_theta_pow_2 = atof(arg[iarg+1]);
      if (cos_theta_pow_2 < 0.0)
        error->all(FLERR,"Illegal surf_collide impulsive cosine power2");
      iarg += 2;
    } else if (strcmp(arg[iarg],"intenergy") == 0) {
      if (iarg+3 > narg)
        error->all(FLERR,"Illegal surf_collide impulsive command");
      intenergy_flag = 1;
      rot_frac = atof(arg[iarg+1]);
      vib_frac = atof(arg[iarg+2]);
      if (rot_frac < 0.0 || rot_frac > 1.0 )
        error->all(FLERR,"Illegal surf_collide impulsive internal energy "
                   "rotational fraction");
      if (vib_frac < 0.0 || vib_frac > 1.0 )
        error->all(FLERR,"Illegal surf_collide impulsive internal energy "
                   "vibrational fraction");
      if (rot_frac + vib_frac > 1.0)
        error->all(FLERR,"Illegal surf_collide impulsive internal energy "
                   "rot-vib fraction sum");
      iarg += 3;
    } else error->all(FLERR,"Illegal surf_collide impulsive command");
  }

  vstream[0] = vstream[1] = vstream[2] = 0.0;

  // initialize RNG

  random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,comm->me,100);
}

/* ---------------------------------------------------------------------- */

SurfCollideImpulsive::~SurfCollideImpulsive()
{
  delete [] tstr;
  delete random;
}

/* ---------------------------------------------------------------------- */

void SurfCollideImpulsive::init()
{
  SurfCollide::init();

  // check variable

  if (tstr) {
    tvar = input->variable->find(tstr);
    if (tvar < 0)
      error->all(FLERR,"Surf_collide impulsive variable name does not exist");
    if (!input->variable->equal_style(tvar))
      error->all(FLERR,"Surf_collide impulsive variable is invalid style");
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

Particle::OnePart *SurfCollideImpulsive::
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

  // impulsive reflection for each particle
  // only if SurfReact did not already reset velocities
  // also both partiticles need to trigger any fixes
  //   to update per-particle properties which depend on
  //   temperature of the particle, e.g. fix vibmode and fix ambipolar

  if (tmode == CUSTOM) twall = tvector[isurf];

  if (ip) {
    if (!velreset) impulsive(ip,norm);
    if (modify->n_update_custom) {
      int i = ip - particle->particles;
      modify->update_custom(i,twall,twall,twall,vstream);
    }
  }
  if (jp) {
    if (!velreset) impulsive(jp,norm);
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
   impulsive reflection
   vrm = most probable speed of species, eqns (4.1) and (4.7)
   vperp = velocity component perpendicular to surface along norm, eqn (12.3)
   vtan12 = 2 velocity components tangential to surface
   tangent1 = component of particle v tangential to surface,
     check if tangent1 = 0 (normal collision), set randomly
   tangent2 = norm x tangent1 = orthogonal tangential direction
   tangent12 are both unit vectors
------------------------------------------------------------------------- */

void SurfCollideImpulsive::impulsive(Particle::OnePart *p, double *norm)
{
  double tangent1[3],tangent2[3];
  Particle::Species *species = particle->species;
  int ispecies = p->ispecies;

  double vperp, vtan1, vtan2;
  double mass = species[ispecies].mass;
  //double vrm = sqrt(2.0*update->boltz * twall / mass);

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

  // compute final polar (theta) and azimuthal (phi) angles

  double tan1 = MathExtra::dot3(v,tangent1);
  double tan2 = MathExtra::dot3(v,tangent2);

  double v_i_mag_sq = MathExtra::lensq3(v);
  double E_i = 0.5 * mass * v_i_mag_sq;
  double theta_i = acos(-dot/sqrt(v_i_mag_sq));
  double phi_i = atan2(tan2,tan1);
  double phi_peak = MY_2PI - phi_i;

  double theta_f, phi_f;
  double P = 0.0;

  // theta_f calculation

  while (random->uniform() > P) {
    theta_f = MY_PI2 * random->uniform();
    P = pow(cos( theta_f - theta_peak ),cos_theta_pow) * sin(theta_f);
    if (double_flag) {
      if (theta_f > theta_peak)
        P = pow(cos( theta_f - theta_peak ),cos_theta_pow_2) * sin(theta_f);
    }

    if (step_flag) {
      double func_step = 0.0;
      double tan_theta = tan(theta_f);
      double cotangent = 1.0/tan_theta;
      if (cotangent > step_size) func_step = 1 - step_size*tan_theta;
      P *= func_step;
    }
  }

  // phi_f calculations

  P = 0.0;
  while (random->uniform() > P) {
    phi_f = phi_peak + MY_PI * (2*random->uniform() - 1);
    P = pow(cos( 0.5*(phi_f - phi_peak) ),cos_phi_pow);
  }

  if (phi_f > MY_PI) phi_f -= MY_2PI;
  else if (phi_f < -MY_PI) phi_f += MY_2PI;

  v_f_avg = 0.0;
  if (softsphere_flag) {
    double mu = species[ispecies].molwt/eff_mass;
    double cos_khi = cos(MY_PI - theta_i - theta_f);
    double sin_khi_sq = 1 - cos_khi*cos_khi;
    double dE, E_f_avg, v_f_mag;

    dE = 2*mu/((mu+1)*(mu+1)) *
      (1 + mu*sin_khi_sq + eng_ratio*(mu+1)/(2*mu) -
       cos_khi*sqrt(1 - mu*mu*sin_khi_sq - eng_ratio*(mu + 1)));
    E_f_avg = E_i * (1 - dE);
    v_f_avg = var_alpha_sq * sqrt(mass/(2*E_f_avg)) *
      (2*E_f_avg/(mass*var_alpha_sq) - 1);
  } else {
    v_f_avg = u0_a*twall + u0_b;
  }

  double v_f_max = 0.5 * (v_f_avg + sqrt(v_f_avg*v_f_avg + 6*var_alpha_sq));
  double f_max = v_f_max*v_f_max*v_f_max *
    exp(-(v_f_max - v_f_avg) * (v_f_max - v_f_avg)/(var_alpha_sq));

  double v_f_mag;
  P = 0.0;
  while (random->uniform() > P) {
    v_f_mag = v_f_max + 3 * var_alpha * ( 2 * random->uniform() - 1 );
    P = v_f_mag*v_f_mag*v_f_mag/(f_max) *
      exp(-(v_f_mag - v_f_avg)*(v_f_mag - v_f_avg)/(var_alpha_sq));
  }

  vperp = v_f_mag * cos(theta_f);
  vtan1 = v_f_mag * sin(theta_f) * cos(phi_f);
  vtan2 = v_f_mag * sin(theta_f) * sin(phi_f);

  v[0] = vperp*norm[0] + vtan1*tangent1[0] + vtan2*tangent2[0];
  v[1] = vperp*norm[1] + vtan1*tangent1[1] + vtan2*tangent2[1];
  v[2] = vperp*norm[2] + vtan1*tangent1[2] + vtan2*tangent2[2];

  //p->erot = particle->erot(ispecies,twall,random);
  //p->evib = particle->evib(ispecies,twall,random);

  if (intenergy_flag) {
    double E_f = 0.5 * mass * v_f_mag * v_f_mag;
    double extra_energy = E_i - E_f;

    // rotational component

    if (!sparta->collide || sparta->collide->rotstyle == NONE ||
        species[ispecies].rotdof < 2) p->erot = 0.0;
    else p->erot += rot_frac*extra_energy;

    // vibrational component

    int vibdof = species[ispecies].vibdof;

    if (!sparta->collide || sparta->collide->vibstyle == NONE || vibdof < 2) {
      p->evib = 0.0;
    } else {
      double *vibtemp = species[ispecies].vibtemp;
      double evib_val = p->evib + vib_frac*extra_energy;

      if (sparta->collide->vibstyle == SMOOTH) p->evib = evib_val;
      if (sparta->collide->vibstyle == DISCRETE && vibdof==2) {
        int ivib =  evib_val / (update->boltz*vibtemp[0]);
        p->evib = ivib * update->boltz * vibtemp[0];
      } else {
        int nvibmode = species[ispecies].nvibmode;
        int *vibdegen = species[ispecies].vibdegen;
        double tot_temp=0.0;
        double evib_sum = 0.0;

        for (int imode=0; imode<nvibmode; imode++)
          tot_temp += vibtemp[imode]*vibdegen[imode];

        for (int imode=0; imode<nvibmode; imode++) {
          int ivib = evib_val / (update->boltz*tot_temp);
          evib_sum += ivib * update->boltz * vibtemp[imode]*vibdegen[imode];
        }

        p->evib = evib_sum;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   wrapper on impulsive() method to perform collision for a single particle
   pass in flags/coefficients to match command-line args for style impulsive
   flags, coeffs can be NULL
   called by SurfReactAdsorb
------------------------------------------------------------------------- */

void SurfCollideImpulsive::wrapper(Particle::OnePart *p, double *norm,
                                   int *flags, double *coeffs)
{
  if (flags) {
    twall = coeffs[0];

    softsphere_flag = flags[0];
    if (softsphere_flag) {
      eng_ratio = coeffs[1];
      eff_mass = coeffs[2];
    } else {
      u0_a = coeffs[1];
      u0_b = coeffs[2];
    }

    var_alpha = coeffs[3];
    theta_peak = coeffs[4];
    cos_theta_pow = coeffs[5];
    cos_phi_pow = coeffs[6];

    step_flag = flags[1];
    double_flag = flags[2];
    intenergy_flag = flags[3];

    int m = 7;

    if (step_flag) {
      step_size = coeffs[m++];
    }
    if (double_flag) {
      cos_theta_pow_2 = coeffs[m++];
    }
    if (intenergy_flag) {
      rot_frac = coeffs[m++];
      vib_frac = coeffs[m++];
    }
  }

  impulsive(p,norm);
}

/* ----------------------------------------------------------------------
   return flags and coeffs for this SurfCollide instance to caller
------------------------------------------------------------------------- */

void SurfCollideImpulsive::flags_and_coeffs(int *flags, double *coeffs)
{
  if (tmode == CUSTOM)
    error->all(FLERR,"Surf_collide impulsive with custom per-surf Twall "
               "does not support external caller");

  coeffs[0] = twall;

  flags[0] = softsphere_flag;
  if (softsphere_flag) {
    coeffs[1] = eng_ratio;
    coeffs[2] = eff_mass;
  } else {
    coeffs[1] = u0_a;
    coeffs[2] = u0_b;
  }

  coeffs[3] = var_alpha;
  coeffs[4] = theta_peak;
  coeffs[5] = cos_theta_pow;
  coeffs[6] = cos_phi_pow;

  flags[1] = step_flag;
  flags[2] = double_flag;
  flags[3] = intenergy_flag;

  int m = 7;

  if (step_flag) {
    coeffs[m++] = step_size;
  }
  if (double_flag) {
    coeffs[m++] = cos_theta_pow_2;
  }
  if (intenergy_flag) {
    coeffs[m++] = rot_frac;
    coeffs[m++] = vib_frac;
  }
}

/* ----------------------------------------------------------------------
   set current surface temperature
------------------------------------------------------------------------- */

void SurfCollideImpulsive::dynamic()
{
  twall = input->variable->compute_equal(tvar);
  if (twall <= 0.0) error->all(FLERR,"Surf_collide impulsive temp <= 0.0");
}

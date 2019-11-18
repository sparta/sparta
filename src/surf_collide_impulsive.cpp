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
#include "random_park.h"
#include "math_const.h"
#include "math_extra.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

SurfCollideImpulsive::SurfCollideImpulsive(SPARTA *sparta, int narg, char **arg) :
  SurfCollide(sparta, narg, arg)
{
  if (narg < 9) error->all(FLERR,"Illegal surf_collide impulsive command");

  tstr = NULL;

  if (strstr(arg[2],"v_") == arg[2]) {
    int n = strlen(&arg[2][2]) + 1;
    tstr = new char[n];
    strcpy(tstr,&arg[2][2]);
  } else {
    twall = atof(arg[2]);
    if (twall < 0.0) error->all(FLERR,"Illegal surf_collide impulsive command");
  }

  eng_ratio = atof(arg[3]);
  eff_mass = atof(arg[4]);
  var_alpha = atof(arg[5]);
  theta_peak = atof(arg[6]);
  cos_theta_pow = atof(arg[7]);
  cos_phi_pow = atof(arg[8]);
  
  if (eng_ratio > 1.0) 
    error->all(FLERR,"Illegal surf_collide impulsive energy ratio");
  if (eff_mass <= 0.0) 
    error->all(FLERR,"Illegal surf_collide impulsive effective mass");
  if (var_alpha < 0.0) 
    error->all(FLERR,"Illegal surf_collide impulsive alpha");
  if (theta_peak < 0.0 || theta_peak > 90.0) 
    error->all(FLERR,"Illegal surf_collide impulsive theta peak");
  if (cos_theta_pow < 0.0 || cos_phi_pow < 0.0) 
    error->all(FLERR,"Illegal surf_collide impulsive cosine power");

  theta_peak *= MY_PI/180; 
  var_alpha_sq = var_alpha * var_alpha;

  // optional args

  step_flag = double_flag = slow_flag = 0;
  tflag = rflag = 0;
  step_size = 0;
  cos_theta_pow_2 = 0;
  slow_frac = slow_alpha = slow_u0_a = slow_u0_b = slow_u0 = slow_barrier = 0;

  int iarg = 9;
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
    } else if (strcmp(arg[iarg],"slow") == 0) {
      if (iarg+6 > narg)
        error->all(FLERR,"Illegal surf_collide impulsive command");
      slow_flag = 1;
      slow_frac = atof(arg[iarg+1]);
      slow_alpha = atof(arg[iarg+2]); 
      slow_u0_a = atof(arg[iarg+3]);
      slow_u0_b = atof(arg[iarg+4]); 
      slow_barrier = atof(arg[iarg+5]);
      if (slow_alpha < 0.0) 
        error->all(FLERR,"Illegal surf_collide impulsive slow alpha");
      iarg += 6;
    } else if (strcmp(arg[iarg],"translate") == 0) {
      if (iarg+4 > narg) 
        error->all(FLERR,"Illegal surf_collide impulsive command");
      tflag = 1;
      vx = atof(arg[iarg+1]);
      vy = atof(arg[iarg+2]);
      vz = atof(arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"rotate") == 0) {
      if (iarg+7 > narg) 
        error->all(FLERR,"Illegal surf_collide impulsive command");
      rflag = 1;
      px = atof(arg[iarg+1]);
      py = atof(arg[iarg+2]);
      pz = atof(arg[iarg+3]);
      wx = atof(arg[iarg+4]);
      wy = atof(arg[iarg+5]);
      wz = atof(arg[iarg+6]);
      if (domain->dimension == 2 && pz != 0.0) 
        error->all(FLERR,"Surf_collide impulsive rotation invalid for 2d");
      if (domain->dimension == 2 && (wx != 0.0 || wy != 0.0))
        error->all(FLERR,"Surf_collide impulsive rotation invalid for 2d");
      iarg += 7;
    } else error->all(FLERR,"Illegal surf_collide impulsive command");
  }

  if (tflag && rflag) error->all(FLERR,"Illegal surf_collide impulsive command");
  if (tflag || rflag) trflag = 1;
  else trflag = 0;
  
  if (slow_flag) {
    slow_u0 = slow_u0_a * twall + slow_u0_b;
    slow_alpha_sq = slow_alpha * slow_alpha;
  }

  vstream[0] = vstream[1] = vstream[2] = 0.0;

  // initialize RNG

  random = new RanPark(update->ranmaster->uniform());
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
   norm = surface normal unit vector
   ip = set to NULL if destroyed by chemsitry
   return jp = new particle if created by chemistry
   resets particle(s) to post-collision outward velocity
------------------------------------------------------------------------- */

Particle::OnePart *SurfCollideImpulsive::
collide(Particle::OnePart *&ip, double *norm, double &, int isr)
{
  nsingle++;

  // if surface chemistry defined, attempt reaction
  // reaction = 1 if reaction took place

  Particle::OnePart iorig;
  Particle::OnePart *jp = NULL;
  int reaction = 0;
  
  if (isr >= 0) {
    if (modify->n_surf_react) memcpy(&iorig,ip,sizeof(Particle::OnePart));    
    reaction = surf->sr[isr]->react(ip,norm,jp);
    if (reaction) surf->nreact_one++;    
  }

  // impulsive reflection for each particle
  // if new particle J created, also need to trigger any fixes

  if (ip) impulsive(ip,norm);
  if (jp) {
    impulsive(jp,norm);
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
  double vrm = sqrt(2.0*update->boltz * twall / mass);
    
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
    
  if (slow_flag && random->uniform() < slow_frac) {
    double vrm_n = sqrt(2.0*update->boltz * (twall + slow_barrier) / mass);
    double vrm_t = vrm;
    
    double slow_vf_max = 0.5 * 
      (slow_u0 + sqrt(slow_u0*slow_u0 + 6*slow_alpha_sq));
    double slow_f_max = slow_vf_max*slow_vf_max*slow_vf_max * 
      exp(-(slow_vf_max - slow_u0)*(slow_vf_max - slow_u0)/(slow_alpha_sq));

    double P = 0, slow_vf_mag;
    while (random->uniform() > P) {
      slow_vf_mag = slow_vf_max + 3 * slow_alpha * ( 2 * random->uniform() - 1 );
      P = slow_vf_mag*slow_vf_mag*slow_vf_mag/(slow_f_max) * 
        exp(-(slow_vf_mag - slow_u0)*(slow_vf_mag - slow_u0)/(slow_alpha_sq));
    }
    
    double slow_phi = MY_2PI * random->uniform();
    double slow_theta = atan2(vrm_t * sqrt(-log(random->uniform())),
                              vrm_n * sqrt(-log(random->uniform())));
    
    vperp = slow_vf_mag * cos(slow_theta);
    vtan1 = slow_vf_mag * sin(slow_theta) * cos(slow_phi);
    vtan2 = slow_vf_mag * sin(slow_theta) * sin(slow_phi);   
    
  } else {
    double mu = species[ispecies].molwt/eff_mass;
    double tan1 = MathExtra::dot3(v,tangent1);
    double tan2 = MathExtra::dot3(v,tangent2);   
    
    double v_i_mag_sq = MathExtra::lensq3(v);
    double E_i = 0.5 * mass * v_i_mag_sq;
    double theta_i = acos(-dot/sqrt(v_i_mag_sq));
    double phi_i = atan2(tan2,tan1);
    double phi_peak = MY_2PI - phi_i;
    
    double theta_f, phi_f;
    double P = 0.0; 
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
    
    P = 0.0; 
    while (random->uniform() > P) {
      phi_f = phi_peak + MY_PI * (2*random->uniform() - 1);
      P = pow(cos( 0.5*(phi_f - phi_peak) ),cos_phi_pow);        
    }
    
    if (phi_f > MY_PI) phi_f -= MY_2PI;
    else if (phi_f < -MY_PI) phi_f += MY_2PI;
    
    double cos_khi = cos(MY_PI - theta_i - theta_f);
    double sin_khi_sq = 1 - cos_khi*cos_khi;
    double E_f_avg, v_f_avg, v_f_mag;
    
    E_f_avg = E_i * (1 - 2*mu/((mu+1)*(mu+1)) * 
                     (1 + mu*sin_khi_sq + eng_ratio*(mu+1)/(2*mu) - 
                      cos_khi*sqrt(1 - mu*mu*sin_khi_sq - eng_ratio*(mu + 1))));
    v_f_avg = var_alpha_sq * sqrt(mass/(2*E_f_avg)) * 
      (2*E_f_avg/(mass*var_alpha_sq) - 1);
    
    double v_f_max = 0.5 * (v_f_avg + sqrt(v_f_avg*v_f_avg + 6*var_alpha_sq));
    double f_max = v_f_max*v_f_max*v_f_max * 
      exp(-(v_f_max - v_f_avg) * (v_f_max - v_f_avg)/(var_alpha_sq));

    P = 0.0;
    while (random->uniform() > P) {
      v_f_mag = v_f_max + 3 * var_alpha * ( 2 * random->uniform() - 1 );
      P = v_f_mag*v_f_mag*v_f_mag/(f_max) * 
        exp(-(v_f_mag-v_f_avg)*(v_f_mag-v_f_avg)/(var_alpha_sq));
    }
    
    vperp = v_f_mag * cos(theta_f);
    vtan1 = v_f_mag * sin(theta_f) * cos(phi_f);
    vtan2 = v_f_mag * sin(theta_f) * sin(phi_f); 
  }

  // add in translation or rotation vector if specified
  // only keep portion of vector tangential to surface element
    
  double beta_un,normalized_distbn_fn;

  if (trflag) {
    double vxdelta,vydelta,vzdelta;
    if (tflag) {
      fprintf(screen,"t flag\n");
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

void SurfCollideImpulsive::dynamic()
{
  twall = input->variable->compute_equal(tvar);
}

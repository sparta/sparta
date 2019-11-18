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
#include "surf_collide_cll_impulsive.h"
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

SurfCollideCLLImpulsive::
SurfCollideCLLImpulsive(SPARTA *sparta, int narg, char **arg) :
  SurfCollide(sparta, narg, arg)
{
  if (narg < 14) 
    error->all(FLERR,"Illegal surf_collide cll/impulsive command");

  tstr = NULL;

  if (strstr(arg[2],"v_") == arg[2]) {
    int n = strlen(&arg[2][2]) + 1;
    tstr = new char[n];
    strcpy(tstr,&arg[2][2]);
  } else {
    twall = atof(arg[2]);
    if (twall < 0.0) 
      error->all(FLERR,"Illegal surf_collide cll/impulsive command");
  }

  eng_thresh = atof(arg[3]);
  eng_ratio = atof(arg[4]);
  eff_mass = atof(arg[5]);
  var_alpha = atof(arg[6]);
  theta_peak = atof(arg[7]);
  cos_theta_pow = atof(arg[8]);
  cos_phi_pow = atof(arg[9]);
  
  acc_n = atof(arg[10]);
  acc_t = atof(arg[11]);
  acc_rot = atof(arg[12]);
  acc_vib = atof(arg[13]);
  
  if (eng_ratio > 1.0) 
    error->all(FLERR,"Surf_collide cll/impulsive "
               "energy ratio must be less than 1");
  if (eff_mass <= 0.0) 
    error->all(FLERR,"Surf_collide cll/impulsive "
               "effective mass must be greater than 0");
  if (var_alpha < 0.0) 
    error->all(FLERR,"Surf_collide cll/impulsive "
               "alpha must be greater than or equal to 0");
  if (theta_peak < 0.0 || theta_peak > 90.0) 
    error->all(FLERR,"Surf_collide cll/impulsive "
               "theta peak must be between 0 and 90");
  if (cos_theta_pow < 0.0 || cos_phi_pow < 0.0) 
    error->all(FLERR,"Surf_collide cll/impulsive "
               "cosine powers must be greater than or equal to 0");
  if (acc_n < 0.0 || acc_n > 1.0 || acc_t < 0.0 || acc_t > 1.0 || 
      acc_rot < 0.0 || acc_rot > 1.0 || acc_vib < 0.0 || acc_vib > 1.0) 
    error->all(FLERR,"Surf_collide cll/impulsive accommodation coeffs "
               "must be >= 0 and <= 1");  
  
  //eng_thresh /= update->J2eV;    // NOTE: what is this setting in update?
  theta_peak *= MY_PI/180; 
  var_alpha_sq = var_alpha * var_alpha;

  // optional args

  step_flag = double_flag = slow_flag = 0;
  pflag = 0;
  eccen = 0;
  tflag = rflag = 0;
  step_size = 0;
  cos_theta_pow_2 = 0;
  slow_frac = slow_alpha = slow_u0_a = slow_u0_b = slow_u0 = slow_barrier = 0;

  int iarg = 14;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"step") == 0) {
      if (iarg+2 > narg) 
        error->all(FLERR,"Illegal surf_collide cll/impulsive command");
      step_flag = 1;
      step_size = atof(arg[iarg+1]);
      if (step_size < 0.0) 
        error->all(FLERR,"Illegal surf_collide cll/impulsive step size");
      iarg += 2;
    } else if (strcmp(arg[iarg],"double") == 0) {
      if (iarg+2 > narg) 
        error->all(FLERR,"Illegal surf_collide cll/impulsive command");
      double_flag = 1;
      cos_theta_pow_2 = atof(arg[iarg+1]);
      if (cos_theta_pow_2 < 0.0) 
        error->all(FLERR,"Illegal surf_collide cll/impulsive cosine power2");
      iarg += 2;
    } else if (strcmp(arg[iarg],"slow") == 0) {
      if (iarg+6 > narg) 
        error->all(FLERR,"Illegal surf_collide cll/impulsive command");
      slow_flag = 1;
      slow_frac = atof(arg[iarg+1]);
      slow_alpha = atof(arg[iarg+2]); 
      slow_u0_a = atof(arg[iarg+3]);
      slow_u0_b = atof(arg[iarg+4]); 
      slow_barrier = atof(arg[iarg+5]);
      if (slow_alpha < 0.0) 
        error->all(FLERR,"Illegal surf_collide cll/impulsive slow alpha");
      iarg += 6;
    } else if (strcmp(arg[iarg],"partial") == 0) {
        if (iarg+2 > narg) 
          error->all(FLERR,"Illegal surf_collide cll/impulsive command");
        if (acc_n != acc_t) 
          error->all(FLERR,"Surf_collide cll/impulsive "
                     "partial requires acc_n = acc_t");
        pflag = 1;
        eccen = atof(arg[iarg+1]);
        if (eccen < 0.0 || eccen >= 1.0 ) 
          error->all(FLERR,"Surf_collide cll/impulsive eccentricity "
                     "must be >= 0 and <= 1");
        iarg += 2;  
    } else if (strcmp(arg[iarg],"translate") == 0) {
      if (iarg+4 > narg) 
        error->all(FLERR,"Illegal surf_collide cll/impulsive command");
      tflag = 1;
      vx = atof(arg[iarg+1]);
      vy = atof(arg[iarg+2]);
      vz = atof(arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"rotate") == 0) {
      if (iarg+7 > narg) 
        error->all(FLERR,"Illegal surf_collide cll/impulsive command");
      rflag = 1;
      px = atof(arg[iarg+1]);
      py = atof(arg[iarg+2]);
      pz = atof(arg[iarg+3]);
      wx = atof(arg[iarg+4]);
      wy = atof(arg[iarg+5]);
      wz = atof(arg[iarg+6]);
      if (domain->dimension == 2 && pz != 0.0) 
        error->all(FLERR,"Surf_collide cll/impulsive rotation invalid for 2d");
      if (domain->dimension == 2 && (wx != 0.0 || wy != 0.0))
        error->all(FLERR,"Surf_collide cll/impulsive rotation invalid for 2d");
      iarg += 7;
    } else error->all(FLERR,"Illegal surf_collide cll/impulsive command");
  }

  if (tflag && rflag) 
    error->all(FLERR,"Illegal surf_collide cll/impulsive command");
  if (tflag || rflag) trflag = 1;
  else trflag = 0;

  vstream[0] = vstream[1] = vstream[2] = 0.0;

  // initialize RNG

  random = new RanPark(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,comm->me,100);
}

/* ---------------------------------------------------------------------- */

SurfCollideCLLImpulsive::~SurfCollideCLLImpulsive()
{
  delete [] tstr;
  delete random;
}

/* ---------------------------------------------------------------------- */

void SurfCollideCLLImpulsive::init()
{
  SurfCollide::init();

  // check variable

  if (tstr) {
    tvar = input->variable->find(tstr);
    if (tvar < 0) 
      error->all(FLERR,"Surf_collide cll/impulsive variable name does not exist");
    if (!input->variable->equal_style(tvar))
      error->all(FLERR,"Surf_collide cll/impulsive variable is invalid style");
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

Particle::OnePart *SurfCollideCLLImpulsive::
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

  // CLL/impulsive reflection for each particle
  // if new particle J created, also need to trigger any fixes

  if (ip) cll_impulsive(ip,norm);
  if (jp) {
    cll_impulsive(jp,norm);
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
   cll/impulsive reflection
   vrm = most probable speed of species, eqns (4.1) and (4.7)
   vperp = velocity component perpendicular to surface along norm, eqn (12.3)
   vtan12 = 2 velocity components tangential to surface
   tangent1 = component of particle v tangential to surface,
     check if tangent1 = 0 (normal collision), set randomly
   tangent2 = norm x tangent1 = orthogonal tangential direction
   tangent12 are both unit vectors
------------------------------------------------------------------------- */

void SurfCollideCLLImpulsive::cll_impulsive(Particle::OnePart *p, double *norm)
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
    
  double tan1 = MathExtra::dot3(v,tangent1);
  double tan2 = MathExtra::dot3(v,tangent2);   
    
  double v_i_mag_sq = MathExtra::lensq3(v);
  double E_i = 0.5 * mass * v_i_mag_sq;
    
  // impulsive

  if (E_i > eng_thresh) {
    
    if (slow_flag && random->uniform() < slow_frac) {
      double vrm_n = sqrt(2.0*update->boltz * (twall + slow_barrier) / mass);
      double vrm_t = vrm;
      double slow_vf_max = 0.5 * 
        (slow_u0 + sqrt(slow_u0*slow_u0 + 6*slow_alpha_sq));
      double slow_f_max = slow_vf_max*slow_vf_max*slow_vf_max * 
        exp(-(slow_vf_max - slow_u0)*(slow_vf_max - slow_u0)/(slow_alpha_sq));
      double P = 0, slow_vf_mag;

      while (random->uniform() > P) {
        slow_vf_mag = slow_vf_max + 3.0 * slow_alpha * 
          (2.0 * random->uniform() - 1.0);
        P = slow_vf_mag*slow_vf_mag*slow_vf_mag/(slow_f_max) * 
          exp(-(slow_vf_mag - slow_u0)*(slow_vf_mag - slow_u0)/(slow_alpha_sq));
      }
      
      double slow_phi = MY_2PI * random->uniform();
      double slow_theta = 
        atan2(vrm_t * sqrt(-log(random->uniform())),vrm_n * 
              sqrt(-log(random->uniform())));
    
      vperp = slow_vf_mag * cos(slow_theta);
      vtan1 = slow_vf_mag * sin(slow_theta) * cos(slow_phi);
      vtan2 = slow_vf_mag * sin(slow_theta) * sin(slow_phi);   
      
    } else {
    
      double mu = species[ispecies].molwt/eff_mass;
      double theta_i = acos(-dot/sqrt(v_i_mag_sq));
      double phi_i = atan2(tan2,tan1);
      double phi_peak = MY_2PI - phi_i;
    
      double theta_f,phi_f;
      double P = 0; 

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
    
      P = 0; 
      while (random->uniform() > P) {        
        phi_f = phi_peak + MY_PI * (2*random->uniform() - 1);
        P = pow(cos( 0.5*(phi_f - phi_peak) ),cos_phi_pow);        
      } 
      
      if (phi_f > MY_PI) phi_f -= MY_2PI;
      else if (phi_f < -MY_PI) phi_f += MY_2PI;
    
      double cos_khi = cos(MY_PI - theta_i - theta_f);
      double sin_khi_sq = 1 - cos_khi*cos_khi;
      double E_f_avg, v_f_avg, v_f_mag;
    
      E_f_avg = E_i * 
        (1.0 - 2.0*mu/((mu+1)*(mu+1)) * 
         (1.0 + mu*sin_khi_sq + eng_ratio*(mu+1)/(2*mu) - 
          cos_khi*sqrt(1 - mu*mu*sin_khi_sq - eng_ratio*(mu + 1.0))));

      v_f_avg = var_alpha_sq * sqrt(mass/(2.0*E_f_avg)) * 
        (2.0*E_f_avg/(mass*var_alpha_sq) - 1.0);

      double v_f_max = 0.5 * (v_f_avg + sqrt(v_f_avg*v_f_avg + 6*var_alpha_sq));
      double f_max = v_f_max*v_f_max*v_f_max * 
        exp(-(v_f_max - v_f_avg)*(v_f_max - v_f_avg)/(var_alpha_sq));

      P = 0;
      while (random->uniform() > P) {
        v_f_mag = v_f_max + 3 * var_alpha * ( 2 * random->uniform() - 1 );
        P = v_f_mag*v_f_mag*v_f_mag/(f_max) * 
          exp(-(v_f_mag-v_f_avg)*(v_f_mag-v_f_avg)/(var_alpha_sq));
      }
    
      vperp = v_f_mag * cos(theta_f);
      vtan1 = v_f_mag * sin(theta_f) * cos(phi_f);
      vtan2 = v_f_mag * sin(theta_f) * sin(phi_f); 
    }
    
  // CLL

  } else {
    
    vrm = sqrt(2.0*update->boltz * twall / species[ispecies].mass);
    
    // CLL model normal velocity

    double r_1 = sqrt(-acc_n*log(random->uniform()));
    double theta_1 = MY_2PI * random->uniform();
    double dot_norm = dot/vrm * sqrt(1-acc_n);
    vperp = vrm * sqrt(r_1*r_1 + dot_norm*dot_norm + 
                       2.0*r_1*dot_norm*cos(theta_1) );
    
    // CLL model tangential velocities
    
    double r_2 = sqrt(-acc_t*log(random->uniform()));
    double theta_2 = MY_2PI * random->uniform();
    double vtangent = tan1/vrm * sqrt(1-acc_t);
    vtan1 = vrm * (vtangent + r_2*cos(theta_2));
    vtan2 = vrm * r_2 * sin(theta_2);
    
    // partial keyword:
    // incomplete energy accommodation with partial/fully diffuse scattering
    // adjust final angle of the particle while keeping 
    // the velocity magnitude or speed according to CLL scattering
    
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
  }

  // add in translation or rotation vector if specified
  // only keep portion of vector tangential to surface element
  
  double beta_un,normalized_distbn_fn;

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

void SurfCollideCLLImpulsive::dynamic()
{
  twall = input->variable->compute_equal(tvar);
}

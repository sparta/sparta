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
#include "string.h"
#include "stdlib.h"
#include "collide_vss.h"
#include "grid.h"
#include "update.h"
#include "particle.h"
#include "mixture.h"
#include "collide.h"
#include "react.h"
#include "comm.h"
#include "random_park.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{NONE,DISCRETE,SMOOTH};            // several files
enum{CONSTANT,VARIABLE};

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

CollideVSS::CollideVSS(SPARTA *sparta, int narg, char **arg) :
  Collide(sparta, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal collide command");

  // optional args

  relaxflag = CONSTANT;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"relax") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal collide command");
      if (strcmp(arg[iarg+1],"constant") == 0) relaxflag = CONSTANT;
      else if (strcmp(arg[iarg+1],"variable") == 0) relaxflag = VARIABLE;
      else error->all(FLERR,"Illegal collide command");
      iarg += 2;
    } else error->all(FLERR,"Illegal collide command");
  }

  // proc 0 reads file to extract params for current species
  // broadcasts params to all procs

  nparams = particle->nspecies;
  params = new Params[nparams];
  if (comm->me == 0) read_param_file(arg[2]);
  MPI_Bcast(params,nparams*sizeof(Params),MPI_BYTE,0,world);

  // check that params were read for all species

  for (int i = 0; i < nparams; i++)
    if (params[i].diam < 0.0) {
      char str[128];
      sprintf(str,"Species %s did not appear in VSS parameter file",
	      particle->species[i].id);
      error->one(FLERR,str);
    }

  // allocate per-species prefactor array

  memory->create(prefactor,nparams,nparams,"collide:prefactor");
}

/* ---------------------------------------------------------------------- */

CollideVSS::~CollideVSS()
{
  delete [] params;
  memory->destroy(prefactor);
}

/* ---------------------------------------------------------------------- */

void CollideVSS::init()
{
  // initially read-in per-species params must match current species list

  if (nparams != particle->nspecies)
    error->all(FLERR,"VSS parameters do not match current species");

  Collide::init();
}

/* ----------------------------------------------------------------------
   estimate a good value for vremax for a group pair in any grid cell
   called by Collide parent in init()
------------------------------------------------------------------------- */

double CollideVSS::vremax_init(int igroup, int jgroup)
{
  // parent has set mixture ptr

  Particle::Species *species = particle->species;
  double *vscale = mixture->vscale;
  int *mix2group = mixture->mix2group;
  int nspecies = particle->nspecies;

  double vrmgroup = 0.0;

  for (int isp = 0; isp < nspecies; isp++) {
    if (mix2group[isp] != igroup) continue;
    for (int jsp = 0; jsp < nspecies; jsp++) {
      if (mix2group[jsp] != jgroup) continue;

      double diam = 0.5 * (params[isp].diam + params[jsp].diam);
      double omega = 0.5 * (params[isp].omega + params[jsp].omega);
      double tref = 0.5 * (params[isp].tref + params[jsp].tref);
      double mr = species[isp].mass * species[jsp].mass /
	(species[isp].mass + species[jsp].mass);
      double cxs = diam*diam*MY_PI;
      prefactor[isp][jsp] = cxs *
	pow(2.0*update->boltz*tref/mr,omega-0.5)/tgamma(2.5-omega);
      double beta = MAX(vscale[isp],vscale[jsp]);
      double vrm = 2.0 * cxs * beta;
      vrmgroup = MAX(vrmgroup,vrm);
    }
  }

  return vrmgroup;
}

/* ---------------------------------------------------------------------- */

double CollideVSS::attempt_collision(int icell, int np, double volume)
{
 double fnum = update->fnum;
 double dt = update->dt;

 double nattempt;

 if (remainflag) {
   nattempt = 0.5 * np * (np-1) *
     vremax[icell][0][0] * dt * fnum / volume + remain[icell][0][0];
   remain[icell][0][0] = nattempt - static_cast<int> (nattempt);
 } else 
   nattempt = 0.5 * np * (np-1) *
     vremax[icell][0][0] * dt * fnum / volume + random->uniform();

  return nattempt;
}

/* ---------------------------------------------------------------------- */

double CollideVSS::attempt_collision(int icell, int igroup, int jgroup, 
				     double volume)
{
 double fnum = update->fnum;
 double dt = update->dt;

 double nattempt;

 double npairs;
 if (igroup == jgroup) npairs = 0.5 * ngroup[igroup] * (ngroup[igroup]-1);
 else npairs = 0.5 * ngroup[igroup] * ngroup[jgroup];

 nattempt = npairs * vremax[icell][igroup][jgroup] * dt * fnum / volume;

 if (remainflag) {
   nattempt += remain[icell][igroup][jgroup];
   remain[icell][igroup][jgroup] = nattempt - static_cast<int> (nattempt);
 } else nattempt += random->uniform();

 return nattempt;
}

/* ----------------------------------------------------------------------
   determine if collision actually occurs
   1 = yes, 0 = no
   update vremax either way
------------------------------------------------------------------------- */

int CollideVSS::test_collision(int icell, int igroup, int jgroup,
			       Particle::OnePart *ip, Particle::OnePart *jp)
{
  Particle::Species *species = particle->species;
  Particle::OnePart *particles = particle->particles;

  double *vi = ip->v;
  double *vj = jp->v;
  int ispecies = ip->ispecies;
  int jspecies = jp->ispecies;
  double du  = vi[0] - vj[0];
  double dv  = vi[1] - vj[1];
  double dw  = vi[2] - vj[2];
  double vr2 = du*du + dv*dv + dw*dw;
  double omega1 = params[ispecies].omega;
  double omega2 = params[jspecies].omega;
  double omega = 0.5 * (omega1+omega2);
  double vro  = pow(vr2,1.0-omega);

  // although the vremax is calcualted for the group,
  // the individual collisions calculated species dependent vre

  double vre = vro*prefactor[ispecies][jspecies];
  vremax[icell][igroup][jgroup] = MAX(vre,vremax[icell][igroup][jgroup]);
  if (vre/vremax[icell][igroup][jgroup] < random->uniform()) return 0;
  precoln.vr2 = vr2;
  return 1;
}

/* ---------------------------------------------------------------------- */

void CollideVSS::setup_collision(Particle::OnePart *ip, Particle::OnePart *jp)
{
  Particle::Species *species = particle->species;
  Particle::OnePart *particles = particle->particles;

  int isp = ip->ispecies;
  int jsp = jp->ispecies;

  precoln.vr = sqrt(precoln.vr2);

  precoln.rotdof_i = species[isp].rotdof;
  precoln.rotdof_j = species[jsp].rotdof;

  precoln.vibdof_i = species[isp].vibdof;
  precoln.vibdof_j = species[jsp].vibdof;

  precoln.ave_rotdof = (species[isp].rotdof + species[jsp].rotdof)/2.;
  precoln.ave_vibdof = (species[isp].vibdof + species[jsp].vibdof)/2.;

  precoln.ave_dof = (precoln.ave_rotdof  + precoln.ave_vibdof)/2.;

  precoln.mass_i = species[isp].mass;
  precoln.mass_j = species[jsp].mass;

  precoln.mr = species[isp].mass * species[jsp].mass /
    (species[isp].mass + species[jsp].mass);

  precoln.etrans = 0.5 * precoln.mr * precoln.vr2;
  precoln.erot = ip->erot + jp->erot;
  precoln.evib = ip->evib + jp->evib;

  precoln.eint   = precoln.erot + precoln.evib;
  precoln.etotal = precoln.etrans + precoln.eint;

  postcoln.etrans = precoln.etrans;
  postcoln.erot = 0.0;
  postcoln.evib = 0.0;
  postcoln.eint = 0.0;
  postcoln.etotal = precoln.etotal;
}

/* ---------------------------------------------------------------------- */

Particle::OnePart *CollideVSS::perform_collision(Particle::OnePart *&ip, 
                                                 Particle::OnePart *&jp)
{
  double x[3],v[3];

  Particle::Species *species = particle->species;
  int reaction,kspecies;

  if (react) { 
    reaction = react->attempt(ip,jp,
                              precoln.etrans,precoln.erot,
                              precoln.evib,postcoln.etotal,kspecies);
  }
  else reaction = 0;
 
  // add a 3rd particle if necessary, index = nlocal-1
  // if add_particle performs a realloc:
  //   make copy of x,v, then repoint ip,jp to new particles data struct

  Particle::OnePart *kp = NULL;

  if (reaction) {
    nreact_one++;
    if (kspecies >= 0) {
      int id = MAXSMALLINT*random->uniform();

      Particle::OnePart *particles = particle->particles;
      memcpy(x,ip->x,3*sizeof(double));
      memcpy(v,ip->v,3*sizeof(double));
      int reallocflag = 
        particle->add_particle(id,kspecies,ip->icell,x,v,0.0,0.0);
      if (reallocflag) {
        ip = particle->particles + (ip - particles);
        jp = particle->particles + (jp - particles);
      }

      kp = &particle->particles[particle->nlocal-1];
      EEXCHANGE_ReactingEDisposal(ip,jp,kp);
      SCATTER_ThreeBodyScattering(ip,jp,kp);

    } else {
      EEXCHANGE_ReactingEDisposal(ip,jp,kp);
      SCATTER_TwoBodyScattering(ip,jp);
    }

  } else { 
    if (precoln.ave_dof > 0.0) EEXCHANGE_NonReactingEDisposal(ip,jp);
    SCATTER_TwoBodyScattering(ip,jp);
  }

  return kp;
}

/* ---------------------------------------------------------------------- */

void CollideVSS::SCATTER_TwoBodyScattering(Particle::OnePart *ip, 
					   Particle::OnePart *jp)
{
  double ua,vb,wc;
  double vrc[3];

  Particle::Species *species = particle->species;
  double *vi = ip->v;
  double *vj = jp->v;
  int isp = ip->ispecies;
  int jsp = jp->ispecies;
  double mass_i = species[isp].mass;
  double mass_j = species[jsp].mass;

  double alpha_r = 2.0 / (params[isp].alpha + params[jsp].alpha);
  double mr = species[isp].mass * species[jsp].mass /
    (species[isp].mass + species[jsp].mass);

  double eps = random->uniform() * 2*MY_PI;
  if (fabs(alpha_r - 1.0) < 0.001) { 
    double vr = sqrt(2.0 * postcoln.etrans / mr);
    double cosX = 2.0*random->uniform() - 1.0;
    double sinX = sqrt(1.0 - cosX*cosX);
    ua = vr*cosX;
    vb = vr*sinX*cos(eps);
    wc = vr*sinX*sin(eps); 
  } else {
    double scale = sqrt((2.0 * postcoln.etrans) / (mr * precoln.vr2));
    double cosX = 2.0*pow(random->uniform(),alpha_r) - 1.0;
    double sinX = sqrt(1.0 - cosX*cosX);
    vrc[0] = vi[0]-vj[0];
    vrc[1] = vi[1]-vj[1];
    vrc[2] = vi[2]-vj[2];
    double d = sqrt(vrc[1]*vrc[1]+vrc[2]*vrc[2]);
    if (d > 1.0e-6) { 
      ua = scale * ( cosX*vrc[0] + sinX*d*sin(eps) );
      vb = scale * ( cosX*vrc[1] + sinX*(precoln.vr*vrc[2]*cos(eps) - 
                                         vrc[0]*vrc[1]*sin(eps))/d );
      wc = scale * ( cosX*vrc[2] - sinX*(precoln.vr*vrc[1]*cos(eps) + 
                                         vrc[0]*vrc[2]*sin(eps))/d );
    } else {
      ua = scale * ( cosX*vrc[0] ); 
      vb = scale * ( sinX*vrc[1]*cos(eps) ); 
      wc = scale * ( sinX*vrc[2]*sin(eps) );
    }
  }
  
  // COM velocity calculated using reactant masses

  double divisor = precoln.mass_i + precoln.mass_j;
  double ucmf = ((precoln.mass_i*vi[0]) + (precoln.mass_j*vj[0])) / divisor;
  double vcmf = ((precoln.mass_i*vi[1]) + (precoln.mass_j*vj[1])) / divisor;
  double wcmf = ((precoln.mass_i*vi[2]) + (precoln.mass_j*vj[2])) / divisor;

  // new velocities for the products

  divisor = mass_i + mass_j;
  vi[0] = ucmf - (mass_j/divisor)*ua;
  vi[1] = vcmf - (mass_j/divisor)*vb;
  vi[2] = wcmf - (mass_j/divisor)*wc;
  vj[0] = ucmf + (mass_i/divisor)*ua;
  vj[1] = vcmf + (mass_i/divisor)*vb;
  vj[2] = wcmf + (mass_i/divisor)*wc;
}

/* ---------------------------------------------------------------------- */

void CollideVSS::EEXCHANGE_NonReactingEDisposal(Particle::OnePart *ip, 
						Particle::OnePart *jp)
{
  double Exp_1,Exp_2,State_prob,Fraction_Rot,Fraction_Vib,E_Dispose;
  int i,rotdof,vibdof,max_level,ivib;

  Particle::OnePart *p;
  Particle::Species *species = particle->species;

  double AdjustFactor = 0.99999999;
  postcoln.erot = 0.0;
  postcoln.evib = 0.0;

  // handle each kind of energy disposal for non-reacting reactants

  if (precoln.ave_dof == 0) {
    ip->erot = 0.0;
    jp->erot = 0.0;
    ip->evib = 0.0;
    jp->evib = 0.0;

  } else {
    E_Dispose = precoln.etrans;

    for (i = 0; i < 2; i++) {
      if (i == 0) p = ip; 
      else p = jp;

      int sp = p->ispecies;
      rotdof = species[sp].rotdof;
      double rotn_phi = species[sp].rotrel; 

      if (rotdof) {
        if (relaxflag == VARIABLE) rotn_phi = rotrel(sp,E_Dispose);
        if (rotn_phi >= random->uniform()) {
          if (rotstyle == NONE) {
            p->erot = 0.0 ; 

          } else if (rotstyle != NONE && rotdof == 2) {
            E_Dispose += p->erot;
            Fraction_Rot = 
              1- pow(random->uniform(),(1/(2.5-params[sp].omega)));
            p->erot = Fraction_Rot * E_Dispose;
            E_Dispose -= p->erot;
          } else {
            E_Dispose += p->erot;
            p->erot = E_Dispose * 
              sample_bl(random,0.5*species[sp].rotdof-1.0,
                        1.5-params[sp].omega);
            E_Dispose -= p->erot;
          }
        }
      }
      postcoln.erot += p->erot;
 
      vibdof = species[sp].vibdof;
      double vibn_phi = species[sp].vibrel; 

      if (vibdof) {
        if (relaxflag == VARIABLE) vibn_phi = vibrel(sp,E_Dispose+p->evib);
        if (vibn_phi >= random->uniform()) {
          if (vibstyle == NONE) {
            p->evib = 0.0; 
          } else if (vibdof == 2 && vibstyle == DISCRETE) {

            E_Dispose += p->evib;
            max_level = static_cast<int>
              (E_Dispose / (update->boltz * species[sp].vibtemp));
            do {
              ivib = static_cast<int> 
                (random->uniform()*(max_level+AdjustFactor));
              p->evib = ivib * update->boltz * species[sp].vibtemp;
              State_prob = pow((1.0 - p->evib / E_Dispose),
                             (1.5 - params[sp].omega));
            } while (State_prob < random->uniform());
            E_Dispose -= p->evib;

          } else if (vibdof == 2 && vibstyle == SMOOTH) {
            E_Dispose += p->evib;
            Fraction_Vib = 
              1.0 - pow(random->uniform(),(1.0/(2.5-params[sp].omega)));
            p->evib= Fraction_Vib * E_Dispose;
            E_Dispose -= p->evib;

          } else if (vibdof > 2) {
            E_Dispose += p->evib;
            p->evib = E_Dispose * 
              sample_bl(random,0.5*species[sp].vibdof-1.0,
                        1.5-params[sp].omega);
            E_Dispose -= p->evib;
          }

        }
      }
      postcoln.evib += p->evib;
    }
  }

  // compute portion of energy left over for scattering

  postcoln.eint = postcoln.erot + postcoln.evib;
  postcoln.etrans = E_Dispose;
}

/* ---------------------------------------------------------------------- */

void CollideVSS::SCATTER_ThreeBodyScattering(Particle::OnePart *ip, 
			  		     Particle::OnePart *jp,
			  		     Particle::OnePart *kp)
{
  double vrc[3],ua,vb,wc;

  Particle::Species *species = particle->species;
  int isp = ip->ispecies;
  int jsp = jp->ispecies;
  int ksp = kp->ispecies;
  double mass_i = species[isp].mass;
  double mass_j = species[jsp].mass;
  double mass_k = species[ksp].mass;
  double mass_ij = mass_i + mass_j;
  double *vi = ip->v;
  double *vj = jp->v;
  double *vk = kp->v;

  double alpha_r = 2.0 / (params[isp].alpha + params[jsp].alpha);
  double mr = mass_ij * mass_k / (mass_ij + mass_k);
  postcoln.eint = ip->erot + jp->erot + ip->evib + jp->evib 
                + kp->erot + kp->evib;

  double cosX = 2.0*pow(random->uniform(), alpha_r) - 1.0;
  double sinX = sqrt(1.0 - cosX*cosX);
  double eps = random->uniform() * 2*MY_PI;

  if (fabs(alpha_r - 1.0) < 0.001) { 
    double vr = sqrt(2*postcoln.etrans/mr);
    ua = vr*cosX;
    vb = vr*sinX*cos(eps);
    wc = vr*sinX*sin(eps);
  } else {
    double scale = sqrt((2.0*postcoln.etrans) / (mr*precoln.vr2));
    vrc[0] = vi[0]-vj[0];
    vrc[1] = vi[1]-vj[1];
    vrc[2] = vi[2]-vj[2];
    double d = sqrt(vrc[1]*vrc[1]+vrc[2]*vrc[2]);
    if (d > 1.E-6 ) { 
      ua = scale * (cosX*vrc[0] + sinX*d*sin(eps));
      vb = scale * (cosX*vrc[1] + sinX*(precoln.vr*vrc[2]*cos(eps) - 
                                        vrc[0]*vrc[1]*sin(eps))/d);
      wc = scale * (cosX*vrc[2] - sinX*(precoln.vr*vrc[1]*cos(eps) + 
                                        vrc[0]*vrc[2]*sin(eps))/d);
    } else {
      ua = scale * cosX*vrc[0]; 
      vb = scale * sinX*vrc[1]*cos(eps); 
      wc = scale * sinX*vrc[2]*sin(eps);
    }
  }

  // COM velocity calculated using reactant masses

  double divisor = precoln.mass_i + precoln.mass_j;
  double ucmf = ((precoln.mass_i*vi[0]) + (precoln.mass_j*vj[0])) / divisor;
  double vcmf = ((precoln.mass_i*vi[1]) + (precoln.mass_j*vj[1])) / divisor;
  double wcmf = ((precoln.mass_i*vi[2]) + (precoln.mass_j*vj[2])) / divisor;

  // new velocities for the products

  divisor = mass_ij + mass_k;
  vi[0] = ucmf - (mass_ij/divisor)*ua;
  vi[1] = vcmf - (mass_ij/divisor)*vb;
  vi[2] = wcmf - (mass_ij/divisor)*wc;
  vk[0] = ucmf + (mass_k/divisor)*ua;
  vk[1] = vcmf + (mass_k/divisor)*vb;
  vk[2] = wcmf + (mass_k/divisor)*wc;
  vj[0] = vi[0];
  vj[1] = vi[1];
  vj[2] = vi[2];
}

/* ---------------------------------------------------------------------- */

void CollideVSS::EEXCHANGE_ReactingEDisposal(Particle::OnePart *ip, 
                                             Particle::OnePart *jp,
                                             Particle::OnePart *kp)
{
  double Exp_1,Exp_2,State_prob,Fraction_Rot,Fraction_Vib;
  int i,numspecies,rotdof,vibdof,max_level,ivib;

  Particle::OnePart *p;
  Particle::Species *species = particle->species;
  double AdjustFactor = 0.99999999;

  if (!kp) {
    ip->erot = 0.0;
    jp->erot = 0.0;
    ip->evib = 0.0;
    jp->evib = 0.0;
    numspecies = 2;
  } else {
    ip->erot = 0.0;
    jp->erot = 0.0;
    kp->erot = 0.0;
    ip->evib = 0.0;
    jp->evib = 0.0;
    kp->evib = 0.0;
    numspecies = 3;
  }

  // handle each kind of energy disposal for non-reacting reactants
  // clean up memory for the products
  
  double E_Dispose = postcoln.etotal;

  for (i = 0; i < numspecies; i++) {
    if (i == 0) p = ip; 
    else if (i == 1) p = jp; 
    else p = kp;

    int sp = p->ispecies;
    rotdof = species[sp].rotdof;

    if (rotdof) {
      if (rotstyle == NONE) {
        p->erot = 0.0 ;
      } else if (rotdof == 2) {
        Fraction_Rot =
          1- pow(random->uniform(),(1/(2.5-params[sp].omega)));
        p->erot = Fraction_Rot * E_Dispose;
        E_Dispose -= p->erot;
        
      } else if (rotdof > 2) {
        p->erot = E_Dispose * 
          sample_bl(random,0.5*species[sp].rotdof-1.0,
                    1.5-params[sp].omega);
        E_Dispose -= p->erot;
      }
    }
    
    vibdof = species[sp].vibdof;

    if (vibdof) {
      if (vibstyle == NONE) {
        p->evib = 0.0;
      } else if (vibdof == 2 && vibstyle == DISCRETE) {
        max_level = static_cast<int> 
          (E_Dispose / (update->boltz * species[sp].vibtemp));
        do {
          ivib = static_cast<int> 
            (random->uniform()*(max_level+AdjustFactor));
          p->evib = (double)
            (ivib * update->boltz * species[sp].vibtemp);
          State_prob = pow((1.0 - p->evib / E_Dispose),
                           (1.5 - params[sp].omega));
        } while (State_prob < random->uniform());
        E_Dispose -= p->evib;
        
      } else if (vibdof == 2 && vibstyle == SMOOTH) {
        Fraction_Vib =
          1.0 - pow(random->uniform(),(1.0 / (2.5-params[sp].omega)));
        p->evib = Fraction_Vib * E_Dispose;
        E_Dispose -= p->evib;
        
      } else if (vibdof > 2) {
        p->evib = E_Dispose * 
          sample_bl(random,0.5*species[sp].vibdof-1.0,
                    1.5-params[sp].omega);
        E_Dispose -= p->evib;
      }
    }
  }
  
  // compute post-collision internal energies
  
  postcoln.erot = ip->erot + jp->erot;
  postcoln.evib = ip->evib + jp->evib;
  
  if (kp) {
    postcoln.erot += kp->erot;
    postcoln.evib += kp->evib;
  }
  
  // compute portion of energy left over for scattering
  
  postcoln.eint = postcoln.erot + postcoln.evib;
  postcoln.etrans = E_Dispose;
}

/* ---------------------------------------------------------------------- */

double CollideVSS::sample_bl(RanPark *random, double Exp_1, double Exp_2)
{
  double Exp_s = Exp_1 + Exp_2;
  double x,y;
  do {
    x = random->uniform();
    y = pow(x*Exp_s/Exp_1, Exp_1)*pow((1.0-x)*Exp_s/Exp_2, Exp_2);
  } while (y < random->uniform());
  return x;
}

/* ----------------------------------------------------------------------
   compute a variable rotational relaxation parameter
------------------------------------------------------------------------- */

double CollideVSS::rotrel(int isp, double Ec)
{
  double Tr = Ec /(update->boltz * (2.5-params[isp].omega));
  double rotphi = (1.0+params[isp].rotc2/sqrt(Tr) + params[isp].rotc3/Tr)
                / params[isp].rotc1; 
  return rotphi;
}

/* ----------------------------------------------------------------------
   compute a variable vibrational relaxation parameter
------------------------------------------------------------------------- */

double CollideVSS::vibrel(int isp, double Ec)
{
  double Tr = Ec /(update->boltz * (3.5-params[isp].omega));
  double omega = params[isp].omega;
  double vibphi = 1.0 / (params[isp].vibc1/pow(Tr,omega) * 
                         exp(params[isp].vibc2/pow(Tr,1.0/3.0)));
  return vibphi;
}

/* ----------------------------------------------------------------------
   read list of species defined in species file
   store info in filespecies and nfilespecies
   only invoked by proc 0
------------------------------------------------------------------------- */

void CollideVSS::read_param_file(char *fname)
{
  FILE *fp = fopen(fname,"r");
  if (fp == NULL) {
    char str[128];
    sprintf(str,"Cannot open VSS parameter file %s",fname);
    error->one(FLERR,str);
  }

  // set all diameters to -1, so can detect if not read

  for (int i = 0; i < nparams; i++) params[i].diam = -1.0;

  // read file line by line
  // skip blank lines or comment lines starting with '#'
  // all other lines must have at least NWORDS 

  int NWORDS = 5;
  if (relaxflag == VARIABLE) NWORDS = 9;
  char **words = new char*[NWORDS];
  char line[MAXLINE],copy[MAXLINE];
  int isp;

  while (fgets(line,MAXLINE,fp)) {
    int pre = strspn(line," \t\n");
    if (pre == strlen(line) || line[pre] == '#') continue;

    strcpy(copy,line);
    int nwords = wordcount(copy);
    if (nwords < NWORDS)
      error->one(FLERR,"Incorrect line format in VSS parameter file");
    wordparse(NWORDS,line,words);

    isp = particle->find_species(words[0]);
    if (isp < 0) continue;

    params[isp].diam = atof(words[1]);
    params[isp].omega = atof(words[2]);
    params[isp].tref = atof(words[3]);
    params[isp].alpha = atof(words[4]);
    if (relaxflag == VARIABLE) {
      params[isp].rotc1 = atof(words[5]);
      params[isp].rotc2 = atof(words[6]);
      params[isp].rotc3 =  (MY_PI+MY_PI2*MY_PI2)*params[isp].rotc2;
      params[isp].rotc2 =  (MY_PI*MY_PIS/2.)*sqrt(params[isp].rotc2);
      params[isp].vibc1 = atof(words[7]);
      params[isp].vibc2 = atof(words[8]);
    }
  }

  delete [] words;
  fclose(fp);
}

/* ----------------------------------------------------------------------
   count whitespace-delimited words in line
------------------------------------------------------------------------- */

int CollideVSS::wordcount(char *line)
{
  int nwords = 0;
  char *word = strtok(line," \t");
  while (word) {
    nwords++;
    word = strtok(NULL," \t");
  }
  return nwords;
}

/* ----------------------------------------------------------------------
   parse first N whitespace-delimited words in line
   store ptr to each word in words
------------------------------------------------------------------------- */

void CollideVSS::wordparse(int n, char *line, char **words)
{
  for (int i = 0; i < n; i++) {
    if (i == 0) words[i] = strtok(line," \t");
    else words[i] = strtok(NULL," \t");
  }
}

/* ----------------------------------------------------------------------
   return a per-species parameter to caller
------------------------------------------------------------------------- */

double CollideVSS::extract(int isp, const char *name)
{
  if (strcmp(name,"diam") == 0) return params[isp].diam;
  else if (strcmp(name,"omega") == 0) return params[isp].omega;
  else if (strcmp(name,"tref") == 0) return params[isp].tref;
  else error->all(FLERR,"Request for unknown parameter from collide");
  return 0.0;
}

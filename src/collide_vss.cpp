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

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

CollideVSS::CollideVSS(SPARTA *sparta, int narg, char **arg) :
  Collide(sparta, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal collision vss command");

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

 if (remainflag) {
   nattempt = 0.5 * ngroup[igroup] * (ngroup[jgroup]-1) *
     vremax[icell][igroup][jgroup] * dt * fnum / volume + 
     remain[icell][igroup][jgroup];
   remain[icell][igroup][jgroup] = nattempt - static_cast<int> (nattempt);
 } else
   nattempt = 0.5 * ngroup[igroup] * (ngroup[jgroup]-1) *
     vremax[icell][igroup][jgroup] * dt * fnum / volume + random->uniform();

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
  precoln.evib = ip->ivib * update->boltz*species[isp].vibtemp + 
    jp->ivib * update->boltz*species[jsp].vibtemp;

  precoln.eint   = precoln.erot + precoln.evib;
  precoln.etotal = precoln.etrans + precoln.eint;

  postcoln.etrans = precoln.etrans;
  postcoln.erot = 0.0;
  postcoln.evib = 0.0;
  postcoln.eint = 0.0;
  postcoln.etotal = precoln.etotal;
}

/* ---------------------------------------------------------------------- */

Particle::OnePart *CollideVSS::perform_collision(Particle::OnePart *ip, 
                                                 Particle::OnePart *jp)
{
  Particle::Species *species = particle->species;
  int reaction,kspecies;

  if (react) 
    reaction = react->attempt(ip,jp,
                              precoln.etrans,precoln.erot,
                              precoln.evib,postcoln.etotal,kspecies);
  else reaction = 0;
 
  // add a 3rd particle if necessary, index = nlocal-1

  Particle::OnePart *kp = NULL;

  if (reaction) {
    nreact_one++;
    if (kspecies >= 0) {
      int id = MAXSMALLINT*random->uniform();
      particle->add_particle(id,kspecies,ip->icell,ip->x,ip->v,0.0,0);
      kp = &particle->particles[particle->nlocal-1];
      double rotdof = precoln.ave_dof+species[kspecies].rotdof;
      if (rotdof > 1.0) EEXCHANGE_ReactingEDisposal(ip,jp,kp);
      SCATTER_ThreeBodyScattering(ip,jp,kp);
    } else {
      if (precoln.ave_dof > 1.0) EEXCHANGE_ReactingEDisposal(ip,jp,kp);
      SCATTER_TwoBodyScattering(ip,jp);
    }
  } else { 
      if (precoln.ave_dof > 1.0) EEXCHANGE_NonReactingEDisposal(ip,jp);
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

  double alpha = 0.5 * (params[isp].alpha + params[jsp].alpha);
  double mr = species[isp].mass * species[jsp].mass /
    (species[isp].mass + species[jsp].mass);

  double eps = random->uniform() * 2*MY_PI;
  double vr = sqrt(2.0 * postcoln.etrans / mr);
  if (abs(alpha) - 1.0 < 0.001) { 
    double cosX = 2.0*random->uniform() - 1.0;
    double sinX = sqrt(1.0 - cosX*cosX);
    ua = vr*cosX;
    vb = vr*sinX*cos(eps);
    wc = vr*sinX*sin(eps); 
  } else {
    double cosX = 2.0*pow(random->uniform(),alpha) - 1.0;
    double sinX = sqrt(1.0 - cosX*cosX);
    vrc[0] = vi[0]-vj[0];
    vrc[1] = vi[1]-vj[1];
    vrc[2] = vi[2]-vj[2];
    double d = sqrt(vrc[1]*vrc[1]+vrc[2]*vrc[2]);
    if (d > 1.E-6 ) { 
      ua = cosX*vrc[0] + sinX*d*sin(eps);
      vb = cosX*vrc[1] + sinX*(vr*vrc[2]*cos(eps) - vrc[0]*vrc[1]*sin(eps))/d;
      wc = cosX*vrc[2] - sinX*(vr*vrc[1]*cos(eps) + vrc[0]*vrc[2]*sin(eps))/d;
    } else {
      ua = cosX*vrc[0]; 
      vb = sinX*vrc[1]*cos(eps); 
      wc = sinX*vrc[2]*sin(eps);
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

  return;
}

/* ---------------------------------------------------------------------- */

void CollideVSS::EEXCHANGE_NonReactingEDisposal(Particle::OnePart *ip, 
						Particle::OnePart *jp)
{
  double Exp_1,Exp_2,State_prob,Fraction_Rot,evib;
  long i,Max_Level,Vib_state;

  Particle::OnePart *p;
  Particle::Species *species = particle->species;

  double AdjustFactor = 0.99999;
  double rotn_phi = 0.5;
  double vibn_phi = 0.5;

  // handle each kind of energy disposal for non-reacting reactants

  if (precoln.ave_dof == 0) {
    ip->erot  = 0.0;
    jp->erot  = 0.0;
    ip->ivib = 0;
    jp->ivib = 0;

  } else {
    double E_Dispose = precoln.etrans;

    for (i = 0; i < 2; i++) {
      if (i == 0) p = ip; 
      else p = jp;
 
      int sp = p->ispecies;

      if (species[sp].vibdof >= 2 && vibn_phi >= random->uniform()) {
        evib = (double) (p->ivib*update->boltz*species[sp].vibtemp);
        E_Dispose += evib;
        Max_Level = (long) (E_Dispose/(update->boltz * species[sp].vibtemp));
        
        do {
          p->ivib = (int) (random->uniform()*(Max_Level+AdjustFactor));
          evib = (double) 
            (p->ivib * update->boltz * species[sp].vibtemp);
          State_prob = pow((1 - evib / E_Dispose),
                           (1.5 - params[sp].omega));
        } while (State_prob < random->uniform());
        
        E_Dispose -= evib;
      }
      
      if (species[sp].rotdof == 2 && rotn_phi >= random->uniform()) {
        E_Dispose += p->erot;
        Fraction_Rot = 
          1- pow(random->uniform(),(1/(2.5-params[sp].omega)));
        p->erot = Fraction_Rot * E_Dispose;
        E_Dispose -= p->erot;
        
      } else if (species[sp].rotdof > 2 && rotn_phi >= random->uniform()) {
        E_Dispose += p->erot;
        Exp_1 = species[sp].rotdof / 2;
        Exp_2 = 2.5 - params[sp].omega;
        do {
          Fraction_Rot = random->uniform();
          State_prob = 
            ((Exp_1+Exp_2-2)/pow(((Exp_1-1)*Fraction_Rot),(Exp_1-1)) *
             ((Exp_1+Exp_2)/(Exp_2-1)*pow((1-Fraction_Rot),(Exp_2-1))));
        } while (State_prob < random->uniform());
        p->erot = Fraction_Rot*E_Dispose;
        E_Dispose -= p->erot;
      }
    }
  }

  // compute post-collision internal energies

  postcoln.erot = ip->erot + jp->erot;
  int isp = ip->ispecies;
  int jsp = jp->ispecies;
  postcoln.evib = ip->ivib * update->boltz*species[isp].vibtemp +
    jp->ivib * update->boltz*species[jsp].vibtemp;
  
  // compute portion of energy left over for scattering

  postcoln.eint = postcoln.erot + postcoln.evib;
  postcoln.etrans = precoln.etotal - postcoln.eint;

 return;
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
  double mass_i = species[isp].mass;
  double mass_j = species[jsp].mass;
  double *vi = ip->v;
  double *vj = jp->v;
  double *vk = kp->v;

  double alpha = 0.5 * (params[isp].alpha + params[jsp].alpha);
  double mr = species[isp].mass * species[jsp].mass /
	     (species[isp].mass + species[jsp].mass);

  postcoln.eint = ip->erot + jp->erot + 
    ip->ivib * update->boltz*species[isp].vibtemp + 
    jp->ivib * update->boltz*species[jsp].vibtemp;

  double cosX = 2.0*pow(random->uniform(), alpha) - 1.0;
  double sinX = sqrt(1.0 - cosX*cosX);
  double eps = random->uniform() * 2*MY_PI;
  double vr = precoln.vr * sqrt(postcoln.etrans/precoln.etrans);

  if (alpha - 1.0 < 0.001) { 
    ua = vr*cosX;
    vb = vr*sinX*cos(eps);
    wc = vr*sinX*sin(eps);
  } else {
    vrc[0] = vi[0]-vj[0];
    vrc[1] = vi[1]-vj[1];
    vrc[2] = vi[2]-vj[2];
    double d = sqrt(vrc[1]*vrc[1]+vrc[2]*vrc[2]);
    if (d > 1.E-6 ) { 
      ua = cosX*vrc[0] + sinX*d*sin(eps);
      vb = cosX*vrc[1] + sinX*(vr*vrc[2]*cos(eps) - vrc[0]*vrc[1]*sin(eps))/d;
      wc = cosX*vrc[2] - sinX*(vr*vrc[1]*cos(eps) + vrc[0]*vrc[2]*sin(eps))/d;
    } else {
      ua = cosX*vrc[0]; 
      vb = sinX*vrc[1]*cos(eps); 
      wc = sinX*vrc[2]*sin(eps);
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
  vk[0] = vj[0];
  vk[1] = vj[1];
  vk[2] = vj[2];

  return;
}

/* ---------------------------------------------------------------------- */

void CollideVSS::EEXCHANGE_ReactingEDisposal(Particle::OnePart *ip, 
                                             Particle::OnePart *jp,
                                             Particle::OnePart *kp)
{
  double Exp_1,Exp_2,State_prob,Fraction_Rot,evib;
  long i,Max_Level,Vib_state;
  long numspecies;

  Particle::OnePart *p;
  Particle::Species *species = particle->species;
  double AdjustFactor = 0.99999;

  // handle each kind of energy disposal for non-reacting reactants

  
  // clean up memory for the products
  
  if (!kp) {
    ip->erot  = 0.0;
    jp->erot  = 0.0;
    ip->ivib = 0;
    jp->ivib = 0;
    numspecies = 2;
  } else {
    ip->erot  = 0.0;
    jp->erot  = 0.0;
    kp->erot  = 0.0;
    ip->ivib = 0;
    jp->ivib = 0;
    kp->ivib = 0;
    numspecies = 3;
  }

  double E_Dispose = postcoln.etotal;

  for (i = 0; i < numspecies; i++) {
    if (i == 0) p = ip; 
    else if (i == 1) p = jp; 
    else p = kp;
    int sp = p->ispecies;
    
    if (species[sp].vibdof >= 2 ) {
      evib = (double) (p->ivib*update->boltz*species[sp].vibtemp);
      E_Dispose += evib;
      Max_Level = (long) (E_Dispose/(update->boltz * species[sp].vibtemp));
      
      do {
        p->ivib = (int) (random->uniform()*(Max_Level+AdjustFactor));
        evib = (double) 
          (p->ivib * update->boltz * species[sp].vibtemp);
        State_prob = pow((1 - evib / E_Dispose),
                         (1.5 - params[sp].omega));
      } while (State_prob < random->uniform());
      
      E_Dispose -= evib;
    }
    
    if ((species[sp].rotdof == 2)) {
      E_Dispose += p->erot;
      Fraction_Rot = 
        1- pow(random->uniform(),(1/(2.5-params[sp].omega)));
      p->erot = Fraction_Rot * E_Dispose;
      E_Dispose -= p->erot;
      
    } else if ((species[sp].rotdof > 2)) {
      E_Dispose += p->erot;
      Exp_1 = species[sp].rotdof / 2;
      Exp_2 = 2.5 - params[sp].omega;
      do {
        Fraction_Rot = random->uniform();
        State_prob = 
          ((Exp_1+Exp_2-2.0)/pow(((Exp_1-1.0)*Fraction_Rot),(Exp_1-1.0)) *
           ((Exp_1+Exp_2)/(Exp_2-1.0)*pow((1.0-Fraction_Rot),(Exp_2-1.0))));
      } while (State_prob < random->uniform());
      p->erot = Fraction_Rot*E_Dispose;
      E_Dispose -= p->erot;
    }
  }
  
  // compute post-collision internal energies
  
  postcoln.erot = ip->erot + jp->erot;
  int isp = ip->ispecies;
  int jsp = jp->ispecies;
  postcoln.evib = ip->ivib * update->boltz*species[isp].vibtemp +
    jp->ivib * update->boltz*species[jsp].vibtemp;
  postcoln.erot = ip->erot + jp->erot;
  
  if (kp) {
    int ksp = kp->ispecies;
    postcoln.evib =+ kp->ivib * update->boltz*species[ksp].vibtemp;
    postcoln.erot =+ kp->erot;
  }
  
  // compute portion of energy left over for scattering
  
  postcoln.eint = postcoln.erot + postcoln.evib;
  postcoln.etrans = postcoln.etotal - postcoln.eint;
  
  return;
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
  // all other lines must have NWORDS 

  int NWORDS = 5;
  char **words = new char*[NWORDS];
  char line[MAXLINE],copy[MAXLINE];
  int isp;

  while (fgets(line,MAXLINE,fp)) {
    int pre = strspn(line," \t\n");
    if (pre == strlen(line) || line[pre] == '#') continue;

    strcpy(copy,line);
    int nwords = wordcount(copy,NULL);
    if (nwords != NWORDS)
      error->one(FLERR,"Incorrect line format in VSS parameter file");
    nwords = wordcount(line,words);

    isp = particle->find_species(words[0]);
    if (isp < 0) continue;

    params[isp].diam = atof(words[1]);
    params[isp].omega = atof(words[2]);
    params[isp].tref = atof(words[3]);
    params[isp].alpha = atof(words[4]);
  }

  delete [] words;
  fclose(fp);
}

/* ----------------------------------------------------------------------
   count whitespace-delimited words in line
   line will be modified, since strtok() inserts NULLs
   if words is non-NULL, store ptr to each word
------------------------------------------------------------------------- */

int CollideVSS::wordcount(char *line, char **words)
{
  int nwords = 0;
  char *word = strtok(line," \t");

  while (word) {
    if (words) words[nwords] = word;
    nwords++;
    word = strtok(NULL," \t");
  }

  return nwords;
}

/* ----------------------------------------------------------------------
   return a per-species parameter to caller
------------------------------------------------------------------------- */

double CollideVSS::extract(int isp, const char *name)
{
  if (strcmp(name,"diam") == 0) return params[isp].diam;
  else if (strcmp(name,"omega") == 0) return params[isp].omega;
  else if (strcmp(name,"tref") == 0) return params[isp].tref;
  else error->all(FLERR,"Request for unknown param from collide");
  return 0.0;
}

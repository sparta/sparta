/* ----------------------------------------------------------------------
   DSMC - Sandia parallel DSMC code
   www.sandia.gov/~sjplimp/dsmc.html
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2011) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level DSMC directory.
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
#include "comm.h"
#include "random_park.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace DSMC_NS;
using namespace MathConst;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

CollideVSS::CollideVSS(DSMC *dsmc, int narg, char **arg) :
  Collide(dsmc, narg, arg)
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

  prefactor = NULL;
  vremax = NULL;
  vrm = NULL;
}

/* ---------------------------------------------------------------------- */

CollideVSS::~CollideVSS()
{
  delete [] params;
  memory->destroy(prefactor);
  memory->destroy(vremax);
  memory->destroy(vrm);
}

/* ---------------------------------------------------------------------- */

void CollideVSS::init()
{
  Collide::init();

  // initially read-in per-species params must match current species list

  if (nparams != particle->nspecies)
    error->all(FLERR,"VSS parameters do not match current species");

  Particle::Species *species = particle->species;
  int nspecies = particle->nspecies;
  int imix = particle->find_mixture(mixID);

  double *vscale = particle->mixture[imix]->vscale;
  int *mix2group = particle->mixture[imix]->mix2group;

  memory->destroy(prefactor);
  memory->destroy(vremax);
  memory->destroy(vrm);

  memory->create(prefactor,nspecies,nspecies,"collide:prefactor");
  memory->create(vremax,grid->nlocal,nspecies,nspecies,"collide:vremax");
  memory->create(vrm,nspecies,nspecies,"collide:vrm");

  // prefactor = static contributions to collision attempt frequencies

  double beta = 1.0;
  double cxs = 1.0;

  for (int isp = 0; isp < nspecies; isp++)
    for (int jsp = 0; jsp < nspecies; jsp++) {
      double diam = 0.5 * (params[isp].diam + params[jsp].diam);
      double omega = 0.5 * (params[isp].omega + params[jsp].omega);
      double tref = 0.5 * (params[isp].tref + params[jsp].tref);
      double mr = species[isp].mass * species[jsp].mass /
	(species[isp].mass + species[jsp].mass);
      cxs = diam*diam*MY_PI;
      prefactor[isp][jsp] = cxs *
	pow(2.0*update->boltz*tref/mr,omega-0.5)/tgamma(2.5-omega);
      //printf(" Prefactor %e %e %e \n", cxs, omega, tgamma(2.5-omega));
      double beta = MIN(vscale[isp],vscale[jsp]);
      double max_thermal_velocity = 3.0/beta;
      vrm[isp][jsp] = cxs * max_thermal_velocity;
    }

  // vremax = max relative velocity factors on per-grid, per-species basis
  // vremax value assignment should be done at a group level.
  // each group should get the maximum vremax of all species involved.
  
  int nglocal = grid->nlocal;

  for (int icell = 0; icell < nglocal; icell++)
    for (int isp = 0; isp < nspecies; isp++)
      for (int jsp = 0; jsp < nspecies; jsp++) {
        int igroup = mix2group[isp];
        int jgroup = mix2group[jsp];
	vremax[icell][igroup][jgroup] = vrm[isp][jsp];
      }
}

/* ---------------------------------------------------------------------- */

double CollideVSS::attempt_collision(int ilocal, int igroup, int jgroup, 
				     double volume)
{
 double fnum = update->fnum;
 double dt = update->dt;

 double nattempt = 0.5 * ngroup[igroup] * (ngroup[jgroup]-1) *
   vremax[ilocal][igroup][jgroup] * dt * fnum / volume + random->uniform();

 // printf(" Attempts = %d %d %d %d\n", icell, ngroup[igroup],ngroup[jgroup],nattempt );
 // if (nattempt > 0) printf(" Attempts = %e %e %e \n", dt, fnum, volume);
 /* no residual number of attempts is kept. The number of attempts is randomly scaled
    up or down. Better formulation for transient flows.
 */

//  printf("%e %e \n",nattempt,vremax[ilocal][igroup][jgroup]);
  return nattempt;
}

/* ----------------------------------------------------------------------
   determine if collision actually occurs
   1 = yes, 0 = no
   update vremax either way
------------------------------------------------------------------------- */

int CollideVSS::test_collision(int ilocal, int igroup, int jgroup,
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
  double vr  = pow(vr2,1-omega);

  // although the vremax is calcualted for the group,
  // the individual collisions calculated species dependent vre

  double vre = vr*prefactor[ispecies][jspecies];

//  printf("INSIDE %e %e \n", vr, prefactor[ispecies][jspecies]);
  // update vremax if new max

  vremax[ilocal][igroup][jgroup] = MAX(vre,vremax[ilocal][igroup][jgroup]);

//  printf("INSIDE %e %e \n", vre, vremax[ilocal][igroup][jgroup]);

  if (vre/vremax[ilocal][igroup][jgroup] < random->uniform()) return 0;
  return 1;
}

/* ---------------------------------------------------------------------- */

void CollideVSS::setup_collision(Particle::OnePart *ip, Particle::OnePart *jp)
{
  Particle::Species *species = particle->species;
  Particle::OnePart *particles = particle->particles;

  int isp = ip->ispecies;
  int jsp = jp->ispecies;

  double *vi = ip->v;
  double *vj = jp->v;

  double du  = vi[0] - vj[0];
  double dv  = vi[1] - vj[1];
  double dw  = vi[2] - vj[2];

  precoln.vr2 = du*du + dv*dv + dw*dw;

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

  postcoln.etrans = 0.0;
  postcoln.erot = 0.0;
  postcoln.evib = 0.0;
  postcoln.eint = 0.0;
  postcoln.etotal = precoln.etotal;
}

/* ---------------------------------------------------------------------- */

Particle::OnePart *CollideVSS::perform_collision(Particle::OnePart *ip, 
						  Particle::OnePart *jp)
{
  EEXCHANGE_NonReactingEDisposal(ip,jp);
  SCATTER_TwoBodyScattering(ip,jp);
  return NULL;
}

/* ---------------------------------------------------------------------- */

void CollideVSS::SCATTER_TwoBodyScattering(Particle::OnePart *ip, 
					   Particle::OnePart *jp)
{
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

  postcoln.eint = ip->erot + jp->erot + 
    ip->ivib * update->boltz*species[isp].vibtemp + 
    jp->ivib * update->boltz*species[jsp].vibtemp;

  double cosX = 2.0*pow(random->uniform(), alpha) - 1.0;
  double sinX = sqrt(1.0 - cosX*cosX);
  double eps = random->uniform() * 2*MY_PI;
  double vr = sqrt(2.0 * postcoln.etrans / mr);

  double ua = vr*cosX;
  double vb = vr*sinX*cos(eps);
  double wc = vr*sinX*sin(eps);

  // centre_of_mass velocity calculated using reactant masses

  double divisor = precoln.mass_i + precoln.mass_j;
  double ucmf = ((precoln.mass_i*vi[0]) + (precoln.mass_j*vj[0])) / divisor;
  double vcmf = ((precoln.mass_i*vi[1]) + (precoln.mass_j*vj[1])) / divisor;
  double wcmf = ((precoln.mass_i*vi[2]) + (precoln.mass_j*vj[2])) / divisor;

  // new velocities for the products

  // printf("In Scatter Pre %e %e %e \n", vi[0], vi[1], vi[2]);
  divisor = mass_i + mass_j;
  vi[0] = ucmf - (mass_j/divisor)*ua;
  vi[1] = vcmf - (mass_j/divisor)*vb;
  vi[2] = wcmf - (mass_j/divisor)*wc;
  vj[0] = ucmf + (mass_i/divisor)*ua;
  vj[1] = vcmf + (mass_i/divisor)*vb;
  vj[2] = wcmf + (mass_i/divisor)*wc;
  // printf("In Scatter aft %e %e %e \n", vi[0], vi[1], vi[2]);

  return;
}

/* ---------------------------------------------------------------------- */

void CollideVSS::EEXCHANGE_NonReactingEDisposal(Particle::OnePart *ip, 
						Particle::OnePart *jp)
{
  Particle::OnePart *kp;
  Particle::Species *species = particle->species;

  static double AdjustFactor = 0.99999;
  static double rotn_phi     = 0.5;
  static double vibn_phi     = 0.5;
  double Exp_1,Exp_2,State_prob,Fraction_Rot,evib;
  long i,Max_Level,Vib_state;

  // handle each kind of energy disposal for non-reacting reactants

  if (precoln.ave_dof == 0) {
    ip->erot  = 0.0;
    jp->erot  = 0.0;
    ip->ivib = 0;
    jp->ivib = 0;

  } else {
    double E_Dispose = precoln.etrans;
    for (i = 0; i < 2; i++) {
     if (i == 0) kp = ip; 
     if (i == 1) kp = jp; 
     int ksp = kp->ispecies;

     if (species[ksp].vibdof >= 2 && (vibn_phi>=random->uniform())) {
       evib = (double) (kp->ivib*update->boltz*species[ksp].vibtemp);
       E_Dispose += evib;
       Max_Level = (long) (E_Dispose/(update->boltz * species[ksp].vibtemp));
       
       do {
	 kp->ivib = (int) (random->uniform()*(Max_Level+AdjustFactor));
	 evib = (double) 
	   (kp->ivib * update->boltz * species[ksp].vibtemp);
	 State_prob = pow((1 - evib / E_Dispose),
			  (1.5 - params[ksp].omega));
       } while (State_prob < random->uniform());
       
       E_Dispose -= evib;
     }
     
     if ((species[ksp].rotdof == 2) && (rotn_phi >= random->uniform())) {
       E_Dispose += kp->erot;
       Fraction_Rot = 
	 1- pow(random->uniform(),(1/(2.5-params[ksp].omega)));
       kp->erot = Fraction_Rot * E_Dispose;
       E_Dispose -= kp->erot;
       
     } else if ((species[ksp].rotdof > 2) && 
		(rotn_phi >= random->uniform())) {
       E_Dispose += kp->erot;
       Exp_1 = species[ksp].rotdof / 2;
       Exp_2 = 2.5 - params[ksp].omega;
       do {
	 Fraction_Rot = random->uniform();
	 State_prob = 
	   ((Exp_1+Exp_2-2)/pow(((Exp_1-1)*Fraction_Rot),(Exp_1-1)) *
	    ((Exp_1+Exp_2)/(Exp_2-1)*pow((1-Fraction_Rot),(Exp_2-1))));
       } while (State_prob < random->uniform());
       kp->erot = Fraction_Rot*E_Dispose;
       E_Dispose -= kp->erot;
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

  /*
  if(postcoln->etranslation<0) {
     fprintf(stderr,
	"NEGATIVE ETRANSLATION IN Energy Disposal! Exiting!!\n");
     printf("NEGATIVE ETRANSLATION IN Energy Disposal! Exiting!!\n");
     printf("ET = %e, EV = %e, ER = %e, ETOT = %e \n", 
	    postcoln->etranslation , 
	    postcoln->erotation , postcoln->evibration,
	    precoln->etotal); 
     exit(0);
  }
  */

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

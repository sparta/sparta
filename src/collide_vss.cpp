/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
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
#include "fix_vibmode.h"
#include "random_knuth.h"
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
  if (nparams == 0)
    error->all(FLERR,"Cannot use collide command with no species defined");

  memory->create(params,nparams,nparams,"collide:params");
  if (comm->me == 0) read_param_file(arg[2]);
  MPI_Bcast(params[0],nparams*nparams*sizeof(Params),MPI_BYTE,0,world);

  // allocate per-species prefactor array

  memory->create(prefactor,nparams,nparams,"collide:prefactor");
}

/* ---------------------------------------------------------------------- */

CollideVSS::~CollideVSS()
{
  if (copymode) return;

  memory->destroy(params);
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

  double *vscale = mixture->vscale;
  int *mix2group = mixture->mix2group;
  int nspecies = particle->nspecies;

  double vrmgroup = 0.0;

  for (int isp = 0; isp < nspecies; isp++) {
    if (mix2group[isp] != igroup) continue;
    for (int jsp = 0; jsp < nspecies; jsp++) {
      if (mix2group[jsp] != jgroup) continue;

      double cxs = params[isp][jsp].diam*params[isp][jsp].diam*MY_PI;
      prefactor[isp][jsp] = cxs * pow(2.0*update->boltz*params[isp][jsp].tref/
        params[isp][jsp].mr,params[isp][jsp].omega-0.5) /
        tgamma(2.5-params[isp][jsp].omega);
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
  } else {
    nattempt = 0.5 * np * (np-1) *
      vremax[icell][0][0] * dt * fnum / volume + random->uniform();
  }
  return nattempt;
}

/* ---------------------------------------------------------------------- */
// ========================================================================
// Modify the number of attempted collision to account for the species
// weight, leading to an artificial increase of the number of 
// trace species numerical particles. Thus, in below equations,
// N from the conventional DSMC is lower than N from SWS methods. 
// Conventional DSMC (sws=0): Ncoll = 1/2 N fnum (N-1)
// SWS (sws=1): Ncoll = 1/2 count_wi fnum (N-1)
// SWSmax (sws=2): Ncoll = 1/2 N fnum wi_max (N-1)
// ========================================================================
double CollideVSS::attempt_collision_SWS(int icell, int np, double volume, double count_wi, double maxwi)
{
  double fnum = update->fnum;
  double dt = update->dt;

  double nattempt;

  int sws = particle->sws;  // SWS
  if (sws==1) {             // SWS
    if (remainflag) {
      nattempt = 0.5 * count_wi * (np-1) *
        vremax[icell][0][0] * dt * fnum / volume + remain[icell][0][0];
      remain[icell][0][0] = nattempt - static_cast<int> (nattempt);
    } else {
      nattempt = 0.5 * count_wi * (np-1) *
        vremax[icell][0][0] * dt * fnum / volume + random->uniform();
    }
  } else if (sws==2) {      // SWS
    if (remainflag) {
      nattempt = 0.5 * np * maxwi * (np-1) *
        vremax[icell][0][0] * dt * fnum / volume + remain[icell][0][0];
      remain[icell][0][0] = nattempt - static_cast<int> (nattempt);
    } else {
      nattempt = 0.5 * np * maxwi * (np-1) *
        vremax[icell][0][0] * dt * fnum / volume + random->uniform();
    }
  } else {
    if (remainflag) {
      nattempt = 0.5 * np * (np-1) *
        vremax[icell][0][0] * dt * fnum / volume + remain[icell][0][0];
      remain[icell][0][0] = nattempt - static_cast<int> (nattempt);
    } else {
      nattempt = 0.5 * np * (np-1) *
        vremax[icell][0][0] * dt * fnum / volume + random->uniform();
    }
  }
  return nattempt;
}

/* ---------------------------------------------------------------------- */

double CollideVSS::attempt_collision(int icell, int igroup, int jgroup,
  double volume)
{
double fnum = update->fnum;
double dt = update->dt;

double nattempt;

// return 2x the value for igroup != jgroup, since no J,I pairing

double npairs;
if (igroup == jgroup) npairs = 0.5 * ngroup[igroup] * (ngroup[igroup]-1);
else npairs = ngroup[igroup] * (ngroup[jgroup]);

//else npairs = 0.5 * ngroup[igroup] * (ngroup[jgroup]);

nattempt = npairs * vremax[icell][igroup][jgroup] * dt * fnum / volume;

if (remainflag) {
nattempt += remain[icell][igroup][jgroup];
remain[icell][igroup][jgroup] = nattempt - static_cast<int> (nattempt);
} else nattempt += random->uniform();

return nattempt;
}

/* ---------------------------------------------------------------------- */
// ========================================================================
// SWS keyword: use n = count_wi * fnum instead of n = np * fnum 
// to obtain the correct number of attempted collisions when using 
// different species weight.
// SWSmax keyword: use n = np * fnum * max(wi)
// ========================================================================
double CollideVSS::attempt_collision_SWS(int icell, int igroup, int jgroup,
                                     double volume)
{
 double fnum = update->fnum;
 double dt = update->dt;

 double nattempt;

 // return 2x the value for igroup != jgroup, since no J,I pairing

 double npairs;

 int sws = particle->sws;  // SWS
 if (igroup == jgroup) {   
  if (sws==1) npairs = 0.5 * count_wi_group[igroup] * (ngroup[igroup]-1);  // SWS
  else if (sws==2) npairs = 0.5 * ngroup[igroup] * maxwigr[igroup] * (ngroup[igroup]-1);  // SWS
  else npairs = 0.5 * ngroup[igroup] * (ngroup[igroup]-1);
 }
 else {
  if (sws==1) npairs = count_wi_group[igroup] * (ngroup[jgroup]);  // SWS
  else if (sws==2) npairs = ngroup[igroup] * maxwigr[igroup] * (ngroup[jgroup]);  // SWS
  else npairs = ngroup[igroup] * (ngroup[jgroup]);
 }

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
  double *vi = ip->v;
  double *vj = jp->v;
  int ispecies = ip->ispecies;
  int jspecies = jp->ispecies;
  double du  = vi[0] - vj[0];
  double dv  = vi[1] - vj[1];
  double dw  = vi[2] - vj[2];
  double vr2 = du*du + dv*dv + dw*dw;
  double vro  = pow(vr2,1.0-params[ispecies][jspecies].omega);

  // although the vremax is calculated for the group,
  // the individual collisions calculated species dependent vre

  double vre = vro*prefactor[ispecies][jspecies];
  vremax[icell][igroup][jgroup] = MAX(vre,vremax[icell][igroup][jgroup]);
  if (vre/vremax[icell][igroup][jgroup] < random->uniform()) return 0;
  precoln.vr2 = vr2;
  return 1;
}

/* ---------------------------------------------------------------------- */
// ========================================================================
// Modify the acceptance-rejection method so that particle pair
// candidates are less likely to be selected if both particles
// have low weights. 
// This modification goes along with: Ncoll = 1/2 N fnum w_imax (N-1),
// N being increased by SWS.
// This method provides the best accuracy for physical particles species 
// collision rates with the cost of a slight increase in the total number of 
// attempted collision.
// Keyword to use this method: SWSmax
// ========================================================================
int CollideVSS::test_collision_SWS(int icell, int igroup, int jgroup,
  Particle::OnePart *ip, Particle::OnePart *jp, double maxwi)
{
double *vi = ip->v;
double *vj = jp->v;
int ispecies = ip->ispecies;
int jspecies = jp->ispecies;
double du  = vi[0] - vj[0];
double dv  = vi[1] - vj[1];
double dw  = vi[2] - vj[2];
double vr2 = du*du + dv*dv + dw*dw;
double vro  = pow(vr2,1.0-params[ispecies][jspecies].omega);

// although the vremax is calculated for the group,
// the individual collisions calculated species dependent vre

double vre = vro*prefactor[ispecies][jspecies];
vremax[icell][igroup][jgroup] = MAX(vre,vremax[icell][igroup][jgroup]);

int sws = particle->sws;  // SWS
if (sws==2) {             // SWS
  Particle::Species *species = particle->species;   // SWS
  double w_ipart = species[ispecies].specwt;   // SWS
  double w_jpart = species[jspecies].specwt;   // SWS
  if ((vre/vremax[icell][igroup][jgroup])*(MAX(w_ipart,w_jpart)/maxwi) < random->uniform()) return 0;  // SWS
} else {
  if (vre/vremax[icell][igroup][jgroup] < random->uniform()) return 0;
}
precoln.vr2 = vr2;
return 1;
}

/* ---------------------------------------------------------------------- */

void CollideVSS::setup_collision(Particle::OnePart *ip, Particle::OnePart *jp)
{
  Particle::Species *species = particle->species;

  int isp = ip->ispecies;
  int jsp = jp->ispecies;

  precoln.vr = sqrt(precoln.vr2);

  precoln.ave_rotdof = 0.5 * (species[isp].rotdof + species[jsp].rotdof);
  precoln.ave_vibdof = 0.5 * (species[isp].vibdof + species[jsp].vibdof);
  precoln.ave_dof = (precoln.ave_rotdof  + precoln.ave_vibdof)/2.;

  double imass = precoln.imass = species[isp].mass;
  double jmass = precoln.jmass = species[jsp].mass;

  precoln.etrans = 0.5 * params[isp][jsp].mr * precoln.vr2;
  precoln.erot = ip->erot + jp->erot;
  precoln.evib = ip->evib + jp->evib;

  precoln.eint   = precoln.erot + precoln.evib;
  precoln.etotal = precoln.etrans + precoln.eint;

  // COM velocity calculated using reactant masses

  double divisor = 1.0 / (imass+jmass);
  double *vi = ip->v;
  double *vj = jp->v;
  precoln.ucmf = ((imass*vi[0])+(jmass*vj[0])) * divisor;
  precoln.vcmf = ((imass*vi[1])+(jmass*vj[1])) * divisor;
  precoln.wcmf = ((imass*vi[2])+(jmass*vj[2])) * divisor;

  postcoln.etrans = precoln.etrans;
  postcoln.erot = 0.0;
  postcoln.evib = 0.0;
  postcoln.eint = 0.0;
  postcoln.etotal = precoln.etotal;
}

/* ---------------------------------------------------------------------- */
// ========================================================================
// Set the weights according to the colliding particle species and 
// determine the maximum weight over all the simulated particles.
// Add the stored energy lost due to differently weighted
// collision to the translation energy of the current collision if
// the colliding particles are major species.
// ========================================================================
void CollideVSS::setup_collision_SWS(Particle::OnePart *ip, Particle::OnePart *jp)
{
  Particle::Species *species = particle->species;

  int isp = ip->ispecies;
  int jsp = jp->ispecies;

  double w_i = species[isp].specwt;  // SWS
  double w_j = species[jsp].specwt;  // SWS
  int nspecies = particle->nspecies; // SWS
  double w_max = 0.0;                // SWS

  for (int i = 0; i < nspecies; i++){  // SWS
    w_max = std::max(species[i].specwt,w_max);
  }

  precoln.vr = sqrt(precoln.vr2);

  precoln.ave_rotdof = 0.5 * (species[isp].rotdof + species[jsp].rotdof);
  precoln.ave_vibdof = 0.5 * (species[isp].vibdof + species[jsp].vibdof);
  precoln.ave_dof = (precoln.ave_rotdof  + precoln.ave_vibdof)/2.;

  double imass = precoln.imass = species[isp].mass;
  double jmass = precoln.jmass = species[jsp].mass;

   if ((w_i==w_max) && (w_j==w_max)){  // SWS
    precoln.etrans = 0.5 * params[isp][jsp].mr * precoln.vr2 + Ewilost;
    Ewilost = 0.0;
  } else {
    precoln.etrans = 0.5 * params[isp][jsp].mr * precoln.vr2;
  }
  
  precoln.erot = ip->erot + jp->erot;
  precoln.evib = ip->evib + jp->evib;

  precoln.eint   = precoln.erot + precoln.evib;
  precoln.etotal = precoln.etrans + precoln.eint;

  // COM velocity calculated using reactant masses

  double divisor = 1.0 / (imass+jmass);
  double *vi = ip->v;
  double *vj = jp->v;
  precoln.ucmf = ((imass*vi[0])+(jmass*vj[0])) * divisor;
  precoln.vcmf = ((imass*vi[1])+(jmass*vj[1])) * divisor;
  precoln.wcmf = ((imass*vi[2])+(jmass*vj[2])) * divisor;

  postcoln.etrans = precoln.etrans;
  postcoln.erot = 0.0;
  postcoln.evib = 0.0;
  postcoln.eint = 0.0;
  postcoln.etotal = precoln.etotal;
}

/* ---------------------------------------------------------------------- */

int CollideVSS::perform_collision(Particle::OnePart *&ip,
                                  Particle::OnePart *&jp,
                                  Particle::OnePart *&kp)
{
  int reactflag,kspecies;
  double x[3],v[3];
  Particle::OnePart *p3;
  
  // if gas-phase chemistry defined, attempt and perform reaction
  // if a 3rd particle is created, its kspecies >= 0 is returned
  // if 2nd particle is removed, its jspecies is set to -1
  
  if (react) 
    reactflag = react->attempt(ip,jp,
                               precoln.etrans,precoln.erot,
                               precoln.evib,postcoln.etotal,kspecies);
  else reactflag = 0;
  
  // repartition energy and perform velocity scattering for I,J,K particles
  // reaction may have changed species of I,J particles
  // J,K particles may have been removed or created by reaction
  
  kp = NULL;
  
  if (reactflag) {
    
    // add 3rd K particle if reaction created it
    // index of new K particle = nlocal-1
    // if add_particle() performs a realloc:
    //   make copy of x,v, then repoint ip,jp to new particles data struct
    //   unless electron
    
    if (kspecies >= 0) {
      int id = MAXSMALLINT*random->uniform();
      
      Particle::OnePart *particles = particle->particles;
      memcpy(x,ip->x,3*sizeof(double));
      memcpy(v,ip->v,3*sizeof(double));
      int ielectron_flag = (ambiflag && ip->ispecies == ambispecies);
      int jelectron_flag = (ambiflag && jp->ispecies == ambispecies);
      int reallocflag =
        particle->add_particle(id,kspecies,ip->icell,x,v,0.0,0.0);
      if (reallocflag) {
        if (!ielectron_flag)
          ip = particle->particles + (ip - particles);
        if (!jelectron_flag)
          jp = particle->particles + (jp - particles);
      }
      
      kp = &particle->particles[particle->nlocal-1];
      EEXCHANGE_ReactingEDisposal(ip,jp,kp);
      SCATTER_ThreeBodyScattering(ip,jp,kp);
    
  // remove 2nd J particle if recombination reaction removed it
  // p3 is 3rd particle participating in energy exchange
  
    } else if (jp->ispecies < 0) {
      double *vi = ip->v;
      double *vj = jp->v;
      
      double divisor = 1.0 / (precoln.imass + precoln.jmass);
      double ucmf = ((precoln.imass*vi[0]) + (precoln.jmass*vj[0])) * divisor;
      double vcmf = ((precoln.imass*vi[1]) + (precoln.jmass*vj[1])) * divisor;
      double wcmf = ((precoln.imass*vi[2]) + (precoln.jmass*vj[2])) * divisor;
      
      vi[0] = ucmf;
      vi[1] = vcmf;
      vi[2] = wcmf;
      
      jp = NULL;
      p3 = react->recomb_part3;
      
      // properly account for 3rd body energy with another call to setup_collision()
      // it needs relative velocity of recombined species and 3rd body
      
      double *vp3 = p3->v;
      double du  = vi[0] - vp3[0];
      double dv  = vi[1] - vp3[1];
      double dw  = vi[2] - vp3[2];
      double vr2 = du*du + dv*dv + dw*dw;
      precoln.vr2 = vr2;
      
      // internal energy of ip particle is already included
      //   in postcoln.etotal returned from react->attempt()
      // but still need to add 3rd body internal energy
      
      double partial_energy =  postcoln.etotal + p3->erot + p3->evib;
      
      ip->erot = 0;
      ip->evib = 0;
      p3->erot = 0;
      p3->evib = 0;
      
      // returned postcoln.etotal will increment only the
      //   relative translational energy between recombined species and 3rd body
      // add back partial_energy to get full total energy
      
      setup_collision(ip,p3);
      postcoln.etotal += partial_energy;
      
      if (precoln.ave_dof > 0.0) EEXCHANGE_ReactingEDisposal(ip,p3,jp);
      SCATTER_TwoBodyScattering(ip,p3);
    
    } else {
      EEXCHANGE_ReactingEDisposal(ip,jp,kp);
      SCATTER_TwoBodyScattering(ip,jp);
    }
  
  } else {
    if (precoln.ave_dof > 0.0) EEXCHANGE_NonReactingEDisposal(ip,jp);
    SCATTER_TwoBodyScattering(ip,jp);
  }
  
  return reactflag;
}

/* ---------------------------------------------------------------------- */

int CollideVSS::perform_collision_SWS(Particle::OnePart *&ip,
                                  Particle::OnePart *&jp,
                                  Particle::OnePart *&kp,
                                  int &n_i,
                                  int &n_j,
                                  int &n_k,
                                  int &n_pre)
{
  int reactflag,kspecies;
  double x[3],v[3];
  Particle::OnePart *p3;
  Particle::OnePart *p_pre;   // SWS
  
  int i;   // SWS
  Particle::Species *species = particle->species;   // SWS

  //Additional parameters for Species weighting scheme
  //All of particles are pre-collision parameters
  int isp = ip->ispecies;
  int jsp = jp->ispecies;
  int ksp = kp->ispecies;
  Particle::OnePart ip_pre = *ip;
  Particle::OnePart jp_pre = *jp;
  Particle::OnePart maxp_pre = *jp;
  // SWS - variables for weighting scheme
  double w_i = species[isp].specwt;
  double w_j = species[jsp].specwt;
  double w_k = species[ksp].specwt;
  double w_min = std::min(w_i, w_j);
  double w_max = std::max(w_i, w_j);
  double phi_i;
  double phi_j;
  double phi_k;

  // if gas-phase chemistry defined, attempt and perform reaction
  // if a 3rd particle is created, its kspecies >= 0 is returned
  // if 2nd particle is removed, its jspecies is set to -1

  if (react)
    reactflag = react->attempt(ip,jp,
                               precoln.etrans,precoln.erot,
                               precoln.evib,postcoln.etotal,kspecies);
  else reactflag = 0;

  // repartition energy and perform velocity scattering for I,J,K particles
  // reaction may have changed species of I,J particles
  // J,K particles may have been removed or created by reaction

  kp = NULL;

  if (reactflag) {
    // compute the number of split particle（0 or 1）
    if (w_i == w_j) n_pre = 0;         // SWS
    else n_pre = ((((w_max-w_min)/w_max)/random->uniform()>1)?1:0);   // SWS
    maxp_pre = ((w_i > w_j) ? ip_pre : jp_pre);      // SWS

    // particle creation of major reactant. This should also be taken account.
    // p_pre is the pointer for a reactant particle that continue to exist
    if ( n_pre == 1) {   // SWS
       int id = MAXSMALLINT*random->uniform();
       Particle::OnePart *particles = particle->particles;
       memcpy(x,maxp_pre.x,3*sizeof(double));
       memcpy(v,maxp_pre.v,3*sizeof(double));
       int reallocflag = 
       particle->add_particle(id,maxp_pre.ispecies,maxp_pre.icell,x,v,maxp_pre.erot,maxp_pre.evib);
       if (reallocflag) {
        ip = particle->particles + (ip - particles);
        jp = particle->particles + (jp - particles);
      }
       p_pre = &particle->particles[particle->nlocal-1];
     }          
    isp = ip->ispecies;  // SWS
    jsp = jp->ispecies;  // SWS
    w_i = species[isp].specwt;  // SWS
    w_j = species[jsp].specwt;  // SWS
    phi_i = w_min/w_i;  // SWS
    phi_j = w_min/w_j;  // SWS
    
    // add 3rd K particle if reaction created it
    // index of new K particle = nlocal-1
    // if add_particle() performs a realloc:
    //   make copy of x,v, then repoint ip,jp to new particles data struct
    //   unless electron

    if (kspecies >= 0) {
      int id = MAXSMALLINT*random->uniform();

      Particle::OnePart *particles = particle->particles;
      memcpy(x,ip->x,3*sizeof(double));
      memcpy(v,ip->v,3*sizeof(double));
      int ielectron_flag = (ambiflag && ip->ispecies == ambispecies);
      int jelectron_flag = (ambiflag && jp->ispecies == ambispecies);
      int reallocflag =
        particle->add_particle(id,kspecies,ip->icell,x,v,0.0,0.0);
      if (reallocflag) {
        if (!ielectron_flag)
          ip = particle->particles + (ip - particles);
        if (!jelectron_flag)
          jp = particle->particles + (jp - particles);
      }

      kp = &particle->particles[particle->nlocal-1];
      
      ksp = kp->ispecies;  // SWS
      w_k = species[ksp].specwt;  // SWS
      phi_k = w_min/w_k;  // SWS
      
      // !! if the reaction is impact ionization, 
      // If the ramdom number [0:1] is lower than the creteria made by two weight,
      // integer will be 1 , else 0
      // we can have number of additional product here
      if (phi_i < 1.0) n_i = (((phi_i)/random->uniform()>1)?1:0);  // SWS
      else if (phi_i >= 1.0) n_i = int(phi_i)+(((phi_i-int(phi_i))/random->uniform()>1)?1:0);  // SWS  
      if (phi_j < 1.0) n_j = (((phi_j)/random->uniform()>1)?1:0);  // SWS
      else if (phi_j >= 1.0) n_j = int(phi_j)+(((phi_j-int(phi_j))/random->uniform()>1)?1:0);  // SWS
      if (phi_k < 1.0) n_k = (((phi_k)/random->uniform()>1)?1:0);  // SWS
      else if (phi_k >= 1.0) n_k = int(phi_k)+(((phi_k-int(phi_k))/random->uniform()>1)?1:0);  // SWS
        
      EEXCHANGE_ReactingEDisposal(ip,jp,kp);
      SCATTER_ThreeBodyScattering(ip,jp,kp);

    // remove 2nd J particle if recombination reaction removed it
    // p3 is 3rd particle participating in energy exchange

    } else if (jp->ispecies < 0) {
      double *vi = ip->v;
      double *vj = jp->v;

      double divisor = 1.0 / (precoln.imass + precoln.jmass);
      double ucmf = ((precoln.imass*vi[0]) + (precoln.jmass*vj[0])) * divisor;
      double vcmf = ((precoln.imass*vi[1]) + (precoln.jmass*vj[1])) * divisor;
      double wcmf = ((precoln.imass*vi[2]) + (precoln.jmass*vj[2])) * divisor;

      vi[0] = ucmf;
      vi[1] = vcmf;
      vi[2] = wcmf;

      jp = NULL;
      p3 = react->recomb_part3;

      // properly account for 3rd body energy with another call to setup_collision()
      // it needs relative velocity of recombined species and 3rd body

      double *vp3 = p3->v;
      double du  = vi[0] - vp3[0];
      double dv  = vi[1] - vp3[1];
      double dw  = vi[2] - vp3[2];
      double vr2 = du*du + dv*dv + dw*dw;
      precoln.vr2 = vr2;

      // internal energy of ip particle is already included
      //   in postcoln.etotal returned from react->attempt()
      // but still need to add 3rd body internal energy

      double partial_energy =  postcoln.etotal + p3->erot + p3->evib;

      ip->erot = 0;
      ip->evib = 0;
      p3->erot = 0;
      p3->evib = 0;

      // returned postcoln.etotal will increment only the
      //   relative translational energy between recombined species and 3rd body
      // add back partial_energy to get full total energy

      setup_collision_SWS(ip,p3);  // SWS
      postcoln.etotal += partial_energy;
      
      if (phi_i < 1.0) n_i = (((phi_i)/random->uniform()>1)?1:0);  // SWS
      else if (phi_i >= 1.0) n_i = int(phi_i)+(((phi_i-int(phi_i))/random->uniform()>1)?1:0);   // SWS
      // because this reaction produce 1 particles, number of j,k particle is 0
      n_j = n_k = 0;    // SWS

      if (precoln.ave_dof > 0.0) EEXCHANGE_ReactingEDisposal(ip,p3,jp);

      // Add reactflag to scattering routine as 
      // splitting-merging should not be used when reaction
      // with differently weighted reactant and/or product occur.
      // Instead, the mass, momentum and energy are conserved
      // through the probability of creation/production of the
      // involved species.
      SCATTER_TwoBodyScattering_SWS(ip,p3,reactflag);  // SWS
    } else {
      // exchange reaction or associative ionization
      // compute number of particle after the reaction
      // !! if the reaction is charge exchange or associative ionization
      // !! reactant and product charge should be same.
      if (phi_i < 1.0) n_i = (((phi_i)/random->uniform()>1)?1:0);  // SWS
      else n_i = int(phi_i)+(((phi_i-int(phi_i))/random->uniform()>1)?1:0);  // SWS
      if (phi_j < 1.0) n_j = (((phi_j)/random->uniform()>1)?1:0);  // SWS
      else n_j = int(phi_j)+(((phi_j-int(phi_j))/random->uniform()>1)?1:0);  // SWS

      // because this reaction produce 2 particles, number of k particle is 0
      n_k = 0;  // SWS
      EEXCHANGE_ReactingEDisposal(ip,jp,kp);
      SCATTER_TwoBodyScattering_SWS(ip,jp,reactflag);  // SWS
    }

  } else {
    // no reaction
    // if reaction is not triggered, particle creation part is skipped
    if (precoln.ave_dof > 0.0) EEXCHANGE_NonReactingEDisposal_SWS(ip,jp);  // SWS
    SCATTER_TwoBodyScattering_SWS(ip,jp,reactflag); // SWS
  }

  return reactflag;
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
  
  double alpha_r = 1.0 / params[isp][jsp].alpha;
  
  double eps = random->uniform() * 2*MY_PI;
  if (fabs(alpha_r - 1.0) < 0.001) {
    double vr = sqrt(2.0 * postcoln.etrans / params[isp][jsp].mr);
    double cosX = 2.0*random->uniform() - 1.0;
    double sinX = sqrt(1.0 - cosX*cosX);
    ua = vr*cosX;
    vb = vr*sinX*cos(eps);
    wc = vr*sinX*sin(eps);
  } else {
    double scale = sqrt((2.0 * postcoln.etrans) / (params[isp][jsp].mr * precoln.vr2));
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
      vb = scale * ( sinX*vrc[0]*cos(eps) );
      wc = scale * ( sinX*vrc[0]*sin(eps) );
    }
  }
  
  // new velocities for the products
  
  double divisor = 1.0 / (mass_i + mass_j);
  vi[0] = precoln.ucmf + (mass_j*divisor)*ua;
  vi[1] = precoln.vcmf + (mass_j*divisor)*vb;
  vi[2] = precoln.wcmf + (mass_j*divisor)*wc;
  vj[0] = precoln.ucmf - (mass_i*divisor)*ua;
  vj[1] = precoln.vcmf - (mass_i*divisor)*vb;
  vj[2] = precoln.wcmf - (mass_i*divisor)*wc;
}
  
/* ---------------------------------------------------------------------- */
// ========================================================================
// Compute the post-collision velocity according using the 
// splitting-merging method.
// Compute and store the energy lost due to non conservation of the energy
// during differently weighted particles collision.
// If reaction occurs, the conservation is ensured by the 
// probobility of creation/deletion of the particles. Thus,
// the splitting-merging method is not used.
// Remark : if equally weighted particles collide phi=1.
// ========================================================================
void CollideVSS::SCATTER_TwoBodyScattering_SWS(Particle::OnePart *ip,
					   Particle::OnePart *jp, int reactflag)
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
  
  // SWS - variables
  double w_i = species[isp].specwt;  
  double w_j = species[jsp].specwt; 
  double phi = 1.0; 

  double vi_pre[3];  
  double vi_post[3]; 
  double vj_pre[3];
  double vj_post[3];

  vi_pre[0]=vi[0];
  vi_pre[1]=vi[1];
  vi_pre[2]=vi[2];

  vj_pre[0]=vj[0];
  vj_pre[1]=vj[1];
  vj_pre[2]=vj[2];

  if ((w_i>0) && (w_j>0)){
    if (w_i>w_j){
      phi = w_j/w_i;
    } else {
      phi = w_i/w_j;
    } 
  }

  double alpha_r = 1.0 / params[isp][jsp].alpha;

  double eps = random->uniform() * 2*MY_PI;
  if (fabs(alpha_r - 1.0) < 0.001) {
    double vr = sqrt(2.0 * postcoln.etrans / params[isp][jsp].mr);
    double cosX = 2.0*random->uniform() - 1.0;
    double sinX = sqrt(1.0 - cosX*cosX);
    ua = vr*cosX;
    vb = vr*sinX*cos(eps);
    wc = vr*sinX*sin(eps);
  } else {
    double scale = sqrt((2.0 * postcoln.etrans) / (params[isp][jsp].mr * precoln.vr2));
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
      vb = scale * ( sinX*vrc[0]*cos(eps) );
      wc = scale * ( sinX*vrc[0]*sin(eps) );
    }
  }

  // new velocities for the products

  double divisor = 1.0 / (mass_i + mass_j);
  vi[0] = precoln.ucmf + (mass_j*divisor)*ua;
  vi[1] = precoln.vcmf + (mass_j*divisor)*vb;
  vi[2] = precoln.wcmf + (mass_j*divisor)*wc;
  vj[0] = precoln.ucmf - (mass_i*divisor)*ua;
  vj[1] = precoln.vcmf - (mass_i*divisor)*vb;
  vj[2] = precoln.wcmf - (mass_i*divisor)*wc;

  if (!(reactflag)){     // SWS
    if (w_i>w_j){
      vj_post[0] = vj[0];
      vj_post[1] = vj[1];
      vj_post[2] = vj[2];
      
      vi_post[0] = (1-phi)*vi_pre[0]+phi*vi[0];
      vi_post[1] = (1-phi)*vi_pre[1]+phi*vi[1];
      vi_post[2] = (1-phi)*vi_pre[2]+phi*vi[2];
  
      Ewilost += w_i*0.5*mass_i*phi*(1-phi)*(
      pow((vi_pre[0]-vi[0]),2.0)+
      pow((vi_pre[1]-vi[1]),2.0)+
      pow((vi_pre[2]-vi[2]),2.0)
      );
  
      vi[0]=vi_post[0];
      vi[1]=vi_post[1];
      vi[2]=vi_post[2];
  
    } else {
      vi_post[0] = vi[0];
      vi_post[1] = vi[1];
      vi_post[2] = vi[2];
      
      vj_post[0] = (1-phi)*vj_pre[0]+phi*vj[0];
      vj_post[1] = (1-phi)*vj_pre[1]+phi*vj[1];
      vj_post[2] = (1-phi)*vj_pre[2]+phi*vj[2];
  
      Ewilost += w_j*0.5*mass_j*phi*(1-phi)*(
      pow((vj_pre[0]-vj[0]),2.0)+
      pow((vj_pre[1]-vj[1]),2.0)+
      pow((vj_pre[2]-vj[2]),2.0)
      );
  
      vj[0]=vj_post[0];
      vj[1]=vj_post[1];
      vj[2]=vj_post[2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void CollideVSS::EEXCHANGE_NonReactingEDisposal(Particle::OnePart *ip,
  Particle::OnePart *jp)
{

  double State_prob,Fraction_Rot,Fraction_Vib,E_Dispose;
  int i,rotdof,vibdof,max_level,ivib;
  
  Particle::OnePart *p;
  Particle::Species *species = particle->species;
  
  double AdjustFactor = 0.99999999;
  postcoln.erot = 0.0;
  postcoln.evib = 0.0;
  double pevib = 0.0;
  
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
        if (relaxflag == VARIABLE) rotn_phi = rotrel(sp,E_Dispose+p->erot);
        if (rotn_phi >= random->uniform()) {
          if (rotstyle == NONE) {
            p->erot = 0.0;
          } else if (rotstyle != NONE && rotdof == 2) {
            E_Dispose += p->erot;
            Fraction_Rot =
              1- pow(random->uniform(),
                    (1/(2.5-params[ip->ispecies][jp->ispecies].omega)));
            p->erot = Fraction_Rot * E_Dispose;
            E_Dispose -= p->erot;
          } else {
            E_Dispose += p->erot;
            p->erot = E_Dispose *
              sample_bl(random,0.5*species[sp].rotdof-1.0,
                        1.5-params[ip->ispecies][jp->ispecies].omega);
            E_Dispose -= p->erot;
          }
        }
      }
      postcoln.erot += p->erot;
      
      vibdof = species[sp].vibdof;
      double vibn_phi = species[sp].vibrel[0];
      
      if (vibdof) {
        if (relaxflag == VARIABLE) vibn_phi = vibrel(sp,E_Dispose+p->evib);
        if (vibn_phi >= random->uniform()) {
          if (vibstyle == NONE) {
            p->evib = 0.0;
      
          } else if (vibdof == 2) {
            if (vibstyle == SMOOTH) {
              E_Dispose += p->evib;
              Fraction_Vib =
              1.0 - pow(random->uniform(),
                       (1.0/(2.5-params[ip->ispecies][jp->ispecies].omega)));
              p->evib= Fraction_Vib * E_Dispose;
              E_Dispose -= p->evib;
      
            } else if (vibstyle == DISCRETE) {
              E_Dispose += p->evib;
              max_level = static_cast<int>
                (E_Dispose / (update->boltz * species[sp].vibtemp[0]));
              do {
                ivib = static_cast<int>
                (random->uniform()*(max_level+AdjustFactor));
                p->evib = ivib * update->boltz * species[sp].vibtemp[0];
                State_prob = pow((1.0 - p->evib / E_Dispose),
                                 (1.5 - params[ip->ispecies][jp->ispecies].omega));
              } while (State_prob < random->uniform());
              E_Dispose -= p->evib;
            }
      
          } else if (vibdof > 2) {
            if (vibstyle == SMOOTH) {
              E_Dispose += p->evib;
              p->evib = E_Dispose *
                sample_bl(random,0.5*species[sp].vibdof-1.0,
                          1.5-params[ip->ispecies][jp->ispecies].omega);
              E_Dispose -= p->evib;
      
            } else if (vibstyle == DISCRETE) {
              p->evib = 0.0;
      
              int nmode = particle->species[sp].nvibmode;
              int **vibmode =
                particle->eiarray[particle->ewhich[index_vibmode]];
              int pindex = p - particle->particles;
      
              for (int imode = 0; imode < nmode; imode++) {
                ivib = vibmode[pindex][imode];
                E_Dispose += ivib * update->boltz *
                  particle->species[sp].vibtemp[imode];
                max_level = static_cast<int>
                  (E_Dispose / (update->boltz * species[sp].vibtemp[imode]));
      
                do {
                  ivib = static_cast<int>
                    (random->uniform()*(max_level+AdjustFactor));
                  pevib = ivib * update->boltz * species[sp].vibtemp[imode];
                  State_prob = pow((1.0 - pevib / E_Dispose),
                                   (1.5 - params[ip->ispecies][jp->ispecies].omega));
                } while (State_prob < random->uniform());
                
                vibmode[pindex][imode] = ivib;
                p->evib += pevib;
                E_Dispose -= pevib;
              }
            }
          } // end of vibstyle/vibdof if
        }
        postcoln.evib += p->evib;
      } // end of vibdof if
    }
  }
  
  // compute portion of energy left over for scattering
  
  postcoln.eint = postcoln.erot + postcoln.evib;
  postcoln.etrans = E_Dispose;
}

/* ---------------------------------------------------------------------- */
// ========================================================================
// Compute the post-collision rotational and vibrational energy
// using the splitting-merging method.
// ========================================================================
void CollideVSS::EEXCHANGE_NonReactingEDisposal_SWS(Particle::OnePart *ip,
                                                Particle::OnePart *jp)
{

  double State_prob,Fraction_Rot,Fraction_Vib,E_Dispose;
  int i,rotdof,vibdof,max_level,ivib;

  Particle::OnePart *p, *p2;
  Particle::Species *species = particle->species;

  double AdjustFactor = 0.99999999;
  postcoln.erot = 0.0;
  postcoln.evib = 0.0;
  double pevib = 0.0;

  // handle each kind of energy disposal for non-reacting reactants

  if (precoln.ave_dof == 0) {
    ip->erot = 0.0;
    jp->erot = 0.0;
    ip->evib = 0.0;
    jp->evib = 0.0;

  } else {
    E_Dispose = precoln.etrans;

    for (i = 0; i < 2; i++) {
      if (i == 0) {
        p = ip;
        p2 = jp;
      }
      else {
        p = jp;
        p2 = ip;
      }  

      // Two different methods are used:
      // 1) Exclusive method:
      // Allow all exchange scenario but the one that imply internal
      // energy exchange of a major particle when colliding with a 
      // minor particle. In this case the energy exchanged is negligible 
      // (proportional to phi). Particles still exchange internal energy 
      // when colliding with identically weighted particles or 
      // more highly weighted particles.
      // 2) Explicit method:
      // Use the splitting merging principle: E''_i=phi*E'_i+(1-phi)*E_i
      // Explicit conservation of the internal energy when a major
      // particle trigger internal energy exchange.
      // NOTE: 
      // The exclusive method is only used for discrete
      // vibrational energy exchange of major particles, as the discrete  
      // energy mode are proportional to integers (ivib) that 
      // are not mergeable without loosing this integer
      // property. All the other exchanges are handled via spitting
      // merging method.

      // SWS - variables
      int psp = p->ispecies;
      int p2sp = p2->ispecies;
      double w_p = species[psp].specwt;
      double w_p2 = species[p2sp].specwt;
      double phi=1.0;
      if ((w_p>0) && (w_p2>0)){ 
        if (w_p>w_p2){
          phi = w_p2/w_p;
        } else {
          phi = w_p/w_p2;
        } 
      }

      int sp = p->ispecies;
      rotdof = species[sp].rotdof;
      double rotn_phi = species[sp].rotrel;

      if (rotdof) {
        if (relaxflag == VARIABLE) rotn_phi = rotrel(sp,E_Dispose+p->erot);
        if (rotn_phi >= random->uniform()) {
          if (w_p>w_p2) {  // SWS
            if (rotstyle == NONE) {
              p->erot = 0.0;
            } else if (rotstyle != NONE && rotdof == 2) {
              double erot_pre = p->erot;
              E_Dispose += p->erot;
              Fraction_Rot =
               1- pow(random->uniform(),
		        (1/(2.5-params[ip->ispecies][jp->ispecies].omega)));
              double erot_post = Fraction_Rot * E_Dispose;
              E_Dispose -= erot_post;
              p->erot = phi * erot_post + (1-phi) * erot_pre;
            } else {
              double erot_pre = p->erot;
              E_Dispose += p->erot;
              double erot_post = E_Dispose *
                sample_bl(random,0.5*species[sp].rotdof-1.0,
                           1.5-params[ip->ispecies][jp->ispecies].omega);
              E_Dispose -= erot_post;
              p->erot = phi * erot_post + (1-phi) * erot_pre; // SWS
            } 
          } else {
            if (rotstyle == NONE) {
              p->erot = 0.0;
            } else if (rotstyle != NONE && rotdof == 2) {
              E_Dispose += p->erot;
              Fraction_Rot =
               1- pow(random->uniform(),
		        (1/(2.5-params[ip->ispecies][jp->ispecies].omega)));
              p->erot = Fraction_Rot * E_Dispose;
              E_Dispose -= p->erot;
            } else {
              E_Dispose += p->erot;
              p->erot = E_Dispose *
                sample_bl(random,0.5*species[sp].rotdof-1.0,
                           1.5-params[ip->ispecies][jp->ispecies].omega);
              E_Dispose -= p->erot;
            }
          }
        }
      }
      postcoln.erot += p->erot;

      vibdof = species[sp].vibdof;
      double vibn_phi = species[sp].vibrel[0];

      if (vibdof) {
        if (relaxflag == VARIABLE) vibn_phi = vibrel(sp,E_Dispose+p->evib);
        if (vibn_phi >= random->uniform()) {
          if (vibstyle == NONE) {
            p->evib = 0.0;

          } else if (vibdof == 2) {
            if (vibstyle == SMOOTH) {
              double e_vib_pre =  p->evib; // SWS
              E_Dispose += p->evib;
              Fraction_Vib =
                1.0 - pow(random->uniform(),
			            (1.0/(2.5-params[ip->ispecies][jp->ispecies].omega)));
              if (w_p>w_p2) {  // SWS
                double e_vib_post = Fraction_Vib * E_Dispose;
                E_Dispose -= e_vib_post;
                p->evib = phi * e_vib_post + (1-phi) * e_vib_pre;  // SWS
              } else {
                p->evib= Fraction_Vib * E_Dispose;
                E_Dispose -= p->evib;
              }

            } else if (vibstyle == DISCRETE) {
              // Added condition (Exclusive mode):
              if (!(w_p>w_p2)) {  // SWS
              E_Dispose += p->evib;
              max_level = static_cast<int>
                (E_Dispose / (update->boltz * species[sp].vibtemp[0]));
              do {
                ivib = static_cast<int>
                  (random->uniform()*(max_level+AdjustFactor));
                p->evib = ivib * update->boltz * species[sp].vibtemp[0];
                State_prob = pow((1.0 - p->evib / E_Dispose),
                                 (1.5 - params[ip->ispecies][jp->ispecies].omega));
              } while (State_prob < random->uniform());
              E_Dispose -= p->evib;
            }
            }

          } else if (vibdof > 2) {
            if (vibstyle == SMOOTH) {
              double e_vib_pre =  p->evib;  // SWS
              E_Dispose += p->evib;
              if (w_p>w_p2) {  // SWS
                double e_vib_post = E_Dispose *
                  sample_bl(random,0.5*species[sp].vibdof-1.0,
                            1.5-params[ip->ispecies][jp->ispecies].omega);
                E_Dispose -= e_vib_post;
                p->evib = phi * e_vib_post + (1-phi) * e_vib_pre;  // SWS
              } else {
                p->evib = E_Dispose *
                  sample_bl(random,0.5*species[sp].vibdof-1.0,
                            1.5-params[ip->ispecies][jp->ispecies].omega);
                E_Dispose -= p->evib;
              }

            } else if (vibstyle == DISCRETE) {
              // Added condition (Exclusive mode):
              if (!(w_p>w_p2)) {  // SWS
              p->evib = 0.0;
              
              int nmode = particle->species[sp].nvibmode;
              int **vibmode =
                particle->eiarray[particle->ewhich[index_vibmode]];
              int pindex = p - particle->particles;
              
              for (int imode = 0; imode < nmode; imode++) {
                ivib = vibmode[pindex][imode];
                E_Dispose += ivib * update->boltz *
                  particle->species[sp].vibtemp[imode];
                max_level = static_cast<int>
                  (E_Dispose / (update->boltz * species[sp].vibtemp[imode]));
              
                do {
                  ivib = static_cast<int>
                    (random->uniform()*(max_level+AdjustFactor));
                  pevib = ivib * update->boltz * species[sp].vibtemp[imode];
                  State_prob = pow((1.0 - pevib / E_Dispose),
                                   (1.5 - params[ip->ispecies][jp->ispecies].omega));
                } while (State_prob < random->uniform());
              
                vibmode[pindex][imode] = ivib;
                p->evib += pevib;
                E_Dispose -= pevib;
              }
              }
            }
          } // end of vibstyle/vibdof if
        }
        postcoln.evib += p->evib;
      } // end of vibdof if
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

  double alpha_r = 1.0 / params[isp][jsp].alpha;
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
      vb = scale * sinX*vrc[0]*cos(eps);
      wc = scale * sinX*vrc[0]*sin(eps);
    }
  }

  // new velocities for the products

  double divisor = 1.0 / (mass_ij + mass_k);
  vi[0] = precoln.ucmf + (mass_k*divisor)*ua;
  vi[1] = precoln.vcmf + (mass_k*divisor)*vb;
  vi[2] = precoln.wcmf + (mass_k*divisor)*wc;
  vk[0] = precoln.ucmf - (mass_ij*divisor)*ua;
  vk[1] = precoln.vcmf - (mass_ij*divisor)*vb;
  vk[2] = precoln.wcmf - (mass_ij*divisor)*wc;
  vj[0] = vi[0];
  vj[1] = vi[1];
  vj[2] = vi[2];
}

/* ---------------------------------------------------------------------- */

void CollideVSS::EEXCHANGE_ReactingEDisposal(Particle::OnePart *ip,
                                             Particle::OnePart *jp,
                                             Particle::OnePart *kp)
{
  double State_prob,Fraction_Rot,Fraction_Vib;
  int i,numspecies,rotdof,vibdof,max_level,ivib;
  double aveomega,pevib;

  Particle::OnePart *p;
  Particle::Species *species = particle->species;
  double AdjustFactor = 0.99999999;

  if (!kp) {
    ip->erot = 0.0;
    jp->erot = 0.0;
    ip->evib = 0.0;
    jp->evib = 0.0;
    numspecies = 2;
    aveomega = params[ip->ispecies][jp->ispecies].omega;
  } else {
    ip->erot = 0.0;
    jp->erot = 0.0;
    kp->erot = 0.0;
    ip->evib = 0.0;
    jp->evib = 0.0;
    kp->evib = 0.0;
    numspecies = 3;
    aveomega = (params[ip->ispecies][ip->ispecies].omega + params[jp->ispecies][jp->ispecies].omega +
                params[kp->ispecies][kp->ispecies].omega)/3;
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
        p->erot = 0.0;
      } else if (rotdof == 2) {
        Fraction_Rot =
          1- pow(random->uniform(),(1/(2.5-aveomega)));
        p->erot = Fraction_Rot * E_Dispose;
        E_Dispose -= p->erot;

      } else if (rotdof > 2) {
        p->erot = E_Dispose *
          sample_bl(random,0.5*species[sp].rotdof-1.0,
                    1.5-aveomega);
        E_Dispose -= p->erot;
      }
    }

    vibdof = species[sp].vibdof;

    if (vibdof) {
      if (vibstyle == NONE) {
        p->evib = 0.0;
      } else if (vibdof == 2 && vibstyle == DISCRETE) {
        max_level = static_cast<int>
          (E_Dispose / (update->boltz * species[sp].vibtemp[0]));
        do {
          ivib = static_cast<int>
            (random->uniform()*(max_level+AdjustFactor));
          p->evib = (double)
            (ivib * update->boltz * species[sp].vibtemp[0]);
          State_prob = pow((1.0 - p->evib / E_Dispose),
                           (1.5 - aveomega));
        } while (State_prob < random->uniform());
        E_Dispose -= p->evib;

      } else if (vibdof == 2 && vibstyle == SMOOTH) {
        Fraction_Vib =
          1.0 - pow(random->uniform(),(1.0 / (2.5-aveomega)));
        p->evib = Fraction_Vib * E_Dispose;
        E_Dispose -= p->evib;

      } else if (vibdof > 2 && vibstyle == SMOOTH) {
          p->evib = E_Dispose *
          sample_bl(random,0.5*species[sp].vibdof-1.0,
                   1.5-aveomega);
          E_Dispose -= p->evib;
      } else if (vibdof > 2 && vibstyle == DISCRETE) {
          p->evib = 0.0;

          int nmode = particle->species[sp].nvibmode;
          int **vibmode = particle->eiarray[particle->ewhich[index_vibmode]];
          int pindex = p - particle->particles;

          for (int imode = 0; imode < nmode; imode++) {
            ivib = vibmode[pindex][imode];
            E_Dispose += ivib * update->boltz *
            particle->species[sp].vibtemp[imode];
            max_level = static_cast<int>
            (E_Dispose / (update->boltz * species[sp].vibtemp[imode]));
            do {
              ivib = static_cast<int>
              (random->uniform()*(max_level+AdjustFactor));
              pevib = ivib * update->boltz * species[sp].vibtemp[imode];
              State_prob = pow((1.0 - pevib / E_Dispose),
                               (1.5 - aveomega));
            } while (State_prob < random->uniform());

            vibmode[pindex][imode] = ivib;
            p->evib += pevib;
            E_Dispose -= pevib;
          }
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

double CollideVSS::sample_bl(RanKnuth *random, double Exp_1, double Exp_2)
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
  // Because we are only relaxing one of the particles in each call, we only
  //  include its DoF, consistent with Bird 2013 (3.32)

  double Tr = Ec /(update->boltz *
                   (2.5-params[isp][isp].omega +
                    particle->species[isp].rotdof/2.0));
  double rotphi = (1.0+params[isp][isp].rotc2/sqrt(Tr) +
                   params[isp][isp].rotc3/Tr) / params[isp][isp].rotc1;
  return rotphi;
}

/* ----------------------------------------------------------------------
   compute a variable vibrational relaxation parameter
------------------------------------------------------------------------- */

double CollideVSS::vibrel(int isp, double Ec)
{
  double Tr = Ec /(update->boltz * (3.5-params[isp][isp].omega));
  double vibphi = 1.0 / (params[isp][isp].vibc1/pow(Tr,params[isp][isp].omega) *
                         exp(params[isp][isp].vibc2/pow(Tr,1.0/3.0)));
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

  // set all species diameters to -1, so can detect if not read
  // set all cross-species parameters to -1 to catch no-reads, as
  // well as user-selected average

  for (int i = 0; i < nparams; i++) {
    params[i][i].diam = -1.0;
    for ( int j = i+1; j<nparams; j++) {
      params[i][j].diam = params[i][j].omega = params[i][j].tref = -1.0;
      params[i][j].alpha = params[i][j].rotc1 = params[i][j].rotc2 = -1.0;
      params[i][j].rotc3 = params[i][j].vibc1 = params[i][j].vibc2 = -1.0;
    }
  }

  // read file line by line
  // skip blank lines or comment lines starting with '#'
  // all other lines must have at least REQWORDS, which depends on VARIABLE flag

  int REQWORDS = 5;
  if (relaxflag == VARIABLE) REQWORDS = 9;
  char **words = new char*[REQWORDS+1]; // one extra word in cross-species lines
  char line[MAXLINE];
  int isp,jsp;

  while (fgets(line,MAXLINE,fp)) {
    int pre = strspn(line," \t\n\r");
    if (pre == strlen(line) || line[pre] == '#') continue;

    int nwords = wordparse(REQWORDS+1,line,words);
    if (nwords < REQWORDS)
      error->one(FLERR,"Incorrect line format in VSS parameter file");

    isp = particle->find_species(words[0]);
    if (isp < 0) continue;

    jsp = particle->find_species(words[1]);

    // if we don't match a species with second word, but it's not a number,
    // skip the line (it involves a species we aren't using)
    if ( jsp < 0 &&  !(atof(words[1]) > 0) ) continue;

    if (jsp < 0 ) {
      params[isp][isp].diam = atof(words[1]);
      params[isp][isp].omega = atof(words[2]);
      params[isp][isp].tref = atof(words[3]);
      params[isp][isp].alpha = atof(words[4]);
      if (relaxflag == VARIABLE) {
        params[isp][isp].rotc1 = atof(words[5]);
        params[isp][isp].rotc2 = atof(words[6]);
        params[isp][isp].rotc3 = (MY_PI+MY_PI2*MY_PI2)*params[isp][isp].rotc2;
        params[isp][isp].rotc2 = (MY_PI*MY_PIS/2.)*sqrt(params[isp][isp].rotc2);
        params[isp][isp].vibc1 = atof(words[7]);
        params[isp][isp].vibc2 = atof(words[8]);
      }
    }else {
      if (nwords < REQWORDS+1)  // one extra word in cross-species lines
        error->one(FLERR,"Incorrect line format in VSS parameter file");
      params[isp][jsp].diam = params[jsp][isp].diam = atof(words[2]);
      params[isp][jsp].omega = params[jsp][isp].omega = atof(words[3]);
      params[isp][jsp].tref = params[jsp][isp].tref = atof(words[4]);
      params[isp][jsp].alpha = params[jsp][isp].alpha = atof(words[5]);
      if (relaxflag == VARIABLE) {
        params[isp][jsp].rotc1 = params[jsp][isp].rotc1 = atof(words[6]);
        params[isp][jsp].rotc2 = atof(words[7]);
        params[isp][jsp].rotc3 = params[jsp][isp].rotc3 =
                        (MY_PI+MY_PI2*MY_PI2)*params[isp][jsp].rotc2;
        if(params[isp][jsp].rotc2 > 0)
                params[isp][jsp].rotc2 = params[jsp][isp].rotc2 =
                                (MY_PI*MY_PIS/2.)*sqrt(params[isp][jsp].rotc2);
        params[isp][jsp].vibc1 = params[jsp][isp].vibc1= atof(words[8]);
        params[isp][jsp].vibc2 = params[jsp][isp].vibc2= atof(words[9]);
      }
    }
  }

  delete [] words;
  fclose(fp);

  // check that params were read for all species
  for (int i = 0; i < nparams; i++) {

    if (params[i][i].diam < 0.0) {
      char str[128];
      sprintf(str,"Species %s did not appear in VSS parameter file",
              particle->species[i].id);
      error->one(FLERR,str);
    }
  }

  for ( int i = 0; i<nparams; i++) {
    params[i][i].mr = particle->species[i].mass / 2;
    for ( int j = i+1; j<nparams; j++) {
      params[i][j].mr = params[j][i].mr = particle->species[i].mass *
        particle->species[j].mass / (particle->species[i].mass + particle->species[j].mass);

      if(params[i][j].diam < 0) params[i][j].diam = params[j][i].diam =
                                  0.5*(params[i][i].diam + params[j][j].diam);
      if(params[i][j].omega < 0) params[i][j].omega = params[j][i].omega =
                                   0.5*(params[i][i].omega + params[j][j].omega);
      if(params[i][j].tref < 0) params[i][j].tref = params[j][i].tref =
                                  0.5*(params[i][i].tref + params[j][j].tref);
      if(params[i][j].alpha < 0) params[i][j].alpha = params[j][i].alpha =
                                   0.5*(params[i][i].alpha + params[j][j].alpha);

      if (relaxflag == VARIABLE) {
        if(params[i][j].rotc1 < 0) params[i][j].rotc1 = params[j][i].rotc1 =
                                     0.5*(params[i][i].rotc1 + params[j][j].rotc1);
        if(params[i][j].rotc2 < 0) params[i][j].rotc2 = params[j][i].rotc2 =
                                     0.5*(params[i][i].rotc2 + params[j][j].rotc2);
        if(params[i][j].rotc3 < 0) params[i][j].rotc3 = params[j][i].rotc3 =
                                     0.5*(params[i][i].rotc3 + params[j][j].rotc3);
        if(params[i][j].vibc1 < 0) params[i][j].vibc1 = params[j][i].vibc1 =
                                     0.5*(params[i][i].vibc1 + params[j][j].vibc1);
        if(params[i][j].vibc2 < 0) params[i][j].vibc2 = params[j][i].vibc2 =
                                     0.5*(params[i][i].vibc2 + params[j][j].vibc2);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   parse up to n=maxwords whitespace-delimited words in line
   store ptr to each word in words and count number of words
------------------------------------------------------------------------- */

int CollideVSS::wordparse(int maxwords, char *line, char **words)
{
  int nwords = 1;
  char * word;

  words[0] = strtok(line," \t\n");
  while ((word = strtok(NULL," \t\n")) != NULL && nwords < maxwords) {
    words[nwords++] = word;
  }
  return nwords;
}

/* ----------------------------------------------------------------------
   return a per-species parameter to caller
------------------------------------------------------------------------- */

double CollideVSS::extract(int isp, int jsp, const char *name)
{
  if (strcmp(name,"diam") == 0) return params[isp][jsp].diam;
  else if (strcmp(name,"omega") == 0) return params[isp][jsp].omega;
  else if (strcmp(name,"tref") == 0) return params[isp][jsp].tref;
  else error->all(FLERR,"Request for unknown parameter from collide");
  return 0.0;
}

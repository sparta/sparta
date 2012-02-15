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
#include "collide_vss.h"
#include "grid.h"
#include "update.h"
#include "random_park.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace DSMC_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

CollideVSS::CollideVSS(DSMC *dsmc, int narg, char **arg) :
  Collide(dsmc, narg, arg)
{
  if (narg != 2) error->all(FLERR,"Illegal collision vss command");

  prefactor = NULL;
  vremax = NULL;
}

/* ---------------------------------------------------------------------- */

CollideVSS::~CollideVSS()
{
  memory->destroy(prefactor);
  memory->destroy(vremax);
}

/* ---------------------------------------------------------------------- */

void CollideVSS::init()
{
  Collide::init();

  /*
  if (nspecies != oldspecies) {
    memory->destroy(prefactor);
    memory->destroy(vremax);
    memory->create(prefactor,nspecies,nspecies,"collide:prefactor");
    memory->create(vremax,grid->nlocal,nspecies,nspecies,"collide:vremax");
  }

  // prefactor = static contributions to collision attempt frequencies

  Particle::Species *species = particle->species;

  double kboltz = 1.0;
  double fratio = 1.0;

  for (int isp = 0; isp < nspecies; isp++)
    for (int jsp = 0; jsp < nspecies; jsp++) {
      double diam = 0.5 * (species[isp].diam + species[jsp].diam);
      double omega = 0.5 * (species[isp].omega + species[jsp].omega);
      double tref = 0.5 * (species[isp].tref + species[jsp].tref);
      double mr = species[isp].mass * species[jsp].mass /
	(species[isp].mass + species[jsp].mass);
      double b = diam*diam / pow(mr,omega);
      prefactor[isp][jsp] = 0.5 * MY_PI * fratio * update->dt * b *
	pow(2.0*(2.0-omega)*kboltz*tref,omega);
    }

  // vremax = max relative velocity factors on per-grid, per-species basis
  
  // NOTE: what is beta and OMEGA

  double beta = 1.0;
  double OMEGA = 1.0;

  vr_indice = 1.0 - 2.0 * OMEGA;
  double max_thermal_velocity = 5.0/beta;
  double vrm = pow(max_thermal_velocity,vr_indice);
   
  int nglocal = grid->nlocal;

  for (int icell = 0; icell < nglocal; icell++)
    for (int isp = 0; isp < nspecies; isp++)
      for (int jsp = 0; jsp < nspecies; jsp++)
	vremax[icell][isp][jsp] = vrm;

  */
}

/* ---------------------------------------------------------------------- */

double CollideVSS::attempt_collision(int icell, 
				     int igroup, int jgroup, double volume)
{
  // dummy attempt frequency

  double attempt = 0.05 * ngroup[igroup] * (ngroup[jgroup]-1);

  /*
  //int nattempt = prefactor[isp][jsp] * nsp[isp] * (nsp[jsp]-1) * 
  //  vremax[icell][isp][jsp] * volume;
  */

  return attempt;
}

/* ----------------------------------------------------------------------
   test if collision actually occurs
   also update vremax
------------------------------------------------------------------------- */

int CollideVSS::test_collision(int icell, int igroup, int jgroup,
			       Particle::OnePart *ip, Particle::OnePart *jp)
{
  if (random->uniform() < 0.5) return 1;
  return 0;

  /*
  double *vi = ip->v;
  double *vj = jp->v;
  double du  = vi[0] - vj[0];
  double dv  = vi[1] - vj[1];
  double dw  = vi[2] - vj[2];
  double vr  = sqrt(du*du + dv*dv + dw*dw);
  double vre = pow(vr,vr_indice);
  
  vremax[icell][isp][jsp] = MAX(vre,vremax[icell][isp][jsp]);
  if (vre/vremax[icell][isp][jsp] < random->uniform()) return 0;
  return 1;
  */
}

/* ---------------------------------------------------------------------- */

void CollideVSS::setup_collision(Particle::OnePart *ip, Particle::OnePart *jp)
{
  /*
  State *precolln  = &cdynamics.precolln;
  State *postcolln = &cdynamics.postcolln;

  precolln->ave_rotn_dof = (SPECIES_Rot_DoF(mola->species_idx) + 
			    SPECIES_Rot_DoF(molb->species_idx)) / 2;
  precolln->ave_vibn_dof = (SPECIES_Vib_DoF(mola->species_idx) +
			    SPECIES_Vib_DoF(molb->species_idx)) / 2;
  precolln->ave_dof      = 
     (precolln->ave_rotn_dof  + precolln->ave_vibn_dof ) / 2;

  precolln->etranslation = 0.5 * precolln->mr * precolln->vrs;
  precolln->erotation    = mola->erotation  + molb->erotation;
  precolln->evibration   = mola->evibration + molb->evibration;
  precolln->einternal    = precolln->erotation + precolln->evibration;
  precolln->etotal       = precolln->etranslation + precolln->einternal;

  postcolln->etranslation = 0.0;
  postcolln->erotation    = 0.0;
  postcolln->evibration   = 0.0;
  postcolln->einternal    = 0.0;
  postcolln->etotal       = precolln->etotal;
  */
}

/* ---------------------------------------------------------------------- */

Particle::OnePart *CollideVSS::perform_collision(Particle::OnePart *ip, 
						  Particle::OnePart *jp)
{
  // dummy collision = just swap velocities

  double tmp;

  double *vi = ip->v;
  double *vj = jp->v;
  tmp = vi[0];
  vi[0] = vj[0];
  vj[0] = tmp;
  tmp = vi[1];
  vi[1] = vj[1];
  vj[1] = tmp;
  tmp = vi[2];
  vi[2] = vj[2];
  vj[2] = tmp;

  return NULL;

  /*
  if (eng_exchange) EEXCHANGE_NonReactingEDisposal(ip,jp);
  SCATTER_TwoBodyScattering(ip,jp);
  return NULL;
  */
}

/* ---------------------------------------------------------------------- */

void CollideVSS::SCATTER_TwoBodyScattering(Particle::OnePart *ip, 
					    Particle::OnePart *jp)
{
  /*
  State *precolln  = &cdynamics.precolln;
  State *postcolln = &cdynamics.postcolln;

  double mass_a    = SPECIES_Mass(mola->species_idx);
  double mass_b    = SPECIES_Mass(molb->species_idx);

  double einternal    = mola->evibration + mola->erotation +
		        molb->evibration + molb->erotation;
  double etranslation = postcolln->etotal - einternal;

  double cosX = 2.0*pow(RANF(), INVALPHA) - 1.0;
  double sinX = sqrt(1.0 - cosX*cosX);
  double eps  = RANF() * TWOPI;
  double vr   = sqrt(2.0 * etranslation / postcolln->mr);

  double ua   = vr*cosX;
  double vb   = vr*sinX*cos(eps);
  double wc   = vr*sinX*sin(eps);

  double du, dv, dw;
  double mr;
  double ke1, ke2;
  double etotal, erotation, evibration;

  // The centre_of_mass velocity is calculated using reactant masses

  double divisor = precolln->mass_a + precolln->mass_b;
  double ucmf =
      ((precolln->mass_a * mola->velocity.u) +
       (precolln->mass_b * molb->velocity.u)) / divisor;
  double vcmf =
      ((precolln->mass_a * mola->velocity.v)+
       (precolln->mass_b * molb->velocity.v)) / divisor;
  double wcmf =
      ((precolln->mass_a * mola->velocity.w)+
       (precolln->mass_b * molb->velocity.w)) / divisor;

  // Compute new velocities for the products

  divisor = mass_a + mass_b;
  mola->velocity.u = ucmf - (mass_b/divisor)*ua;
  mola->velocity.v = vcmf - (mass_b/divisor)*vb;
  mola->velocity.w = wcmf - (mass_b/divisor)*wc;
  molb->velocity.u = ucmf + (mass_a/divisor)*ua;
  molb->velocity.v = vcmf + (mass_a/divisor)*vb;
  molb->velocity.w = wcmf + (mass_a/divisor)*wc;

  // Econservation test

  du = mola->velocity.u - molb->velocity.u;
  dv = mola->velocity.v - molb->velocity.v;
  dw = mola->velocity.w - molb->velocity.w;
  mr = SPECIES_Mr(mola->species_idx, molb->species_idx);

  if (!EQUAL(0.5*mr*(du*du+dv*dv+dw*dw)*1.0E10, etranslation*1.0E10))
    fprintf(stderr,"2-Body Scatter: Violation of EConservation!\n");

  etranslation = 0.5*mr*(du*du+dv*dv+dw*dw);
  erotation    = mola->erotation  + molb->erotation;
  evibration   = mola->evibration + molb->evibration;
  etotal       = (etranslation + erotation + evibration);

  if (!EQUAL(etotal*1.0E10, postcolln->etotal*1.0E10))
    fprintf(stderr,"2-Body Scatter: Violation of EConservation!\n");
  */
}

/* ---------------------------------------------------------------------- */

void CollideVSS::EEXCHANGE_NonReactingEDisposal(Particle::OnePart *ip, 
						 Particle::OnePart *jp)
{
  /*
  static double AdjustFactor = 0.99999;
  static double rotn_phi     = 0.5;

  State *precolln  = &cdynamics.precolln;
  State *postcolln = &cdynamics.postcolln;

  Molecule * mol;
  double E_Dispose,State_prob, Exp_1,Exp_2,Fraction_Rot;
  double rdof, vdof;
  long i, Max_Level,Vib_state;

  // Handle each kind of EDisposal for NonReacting reactants

  switch(cdynamics.state)
  {
     case monatom_monatom:
       mola->erotation  = 0.0;
       molb->erotation  = 0.0;
       mola->evibration = 0.0;
       molb->evibration = 0.0;
       break;

     case monatom_diatom:
     case diatom_diatom:
       E_Dispose = precolln->etranslation;
       for (i = 0; i < 2; i++) {
         if (i == 0) mol = mola; 
         if (i == 1) mol = molb; 
	 rdof = SPECIES_Rot_DoF(mol->species_idx);
	 vdof = SPECIES_Vib_DoF(mol->species_idx);

         if (vdof >= 2 && (vibn_phi>=RANF())) {
	   E_Dispose += mol->evibration;
	   Max_Level  = (long) 
	     (E_Dispose/(BOLTZMANN * SPECIES_Theta(mol->species_idx)));

	   do {
	     Vib_state = (long) (RANF()*(Max_Level+AdjustFactor));
	     mol->evibration = (double) 
 	       (Vib_state*BOLTZMANN *SPECIES_Theta(mol->species_idx));
	     State_prob = pow((1 - mol->evibration/E_Dispose),
			      (1.5-SPECIES_Omega(mol->species_idx)));
	   } while (State_prob < RANF());

	   E_Dispose -= mol->evibration;
         }

         if ((rdof == 2) && (rotn_phi >= RANF())) {
           E_Dispose += mol->erotation;
           Fraction_Rot = 
	     1- pow(RANF(),(1/(2.5-SPECIES_Omega(mol->species_idx))));
	   mol->erotation = Fraction_Rot*E_Dispose;
	   E_Dispose -= mol->erotation;

         } else if ((vdof > 2) && (rotn_phi >= RANF())) {
           E_Dispose += mol->erotation;
	   Exp_1 = rdof / 2;
	   Exp_2 = 2.5 - SPECIES_Omega(mol->species_idx);
	   do {
	     Fraction_Rot = RANF();
	     State_prob   = 
	       ((Exp_1+Exp_2-2)/pow(((Exp_1-1)*Fraction_Rot),(Exp_1-1)) *
	       ((Exp_1+Exp_2)/(Exp_2-1)*pow((1-Fraction_Rot),(Exp_2-1))));
	   } while (State_prob < RANF());
	   mol->erotation = Fraction_Rot*E_Dispose;
	   E_Dispose -= mol->erotation;
         }
       } 
       break;

     default:
       fprintf(stderr,"Nonsensical cdynamics - NonReactingEDisposal");
       printf("Nonsensical cdynamics - NonReactingEDisposal");
       exit(0);
   }

  // Compute the post-collision internal energies

  postcolln->erotation = mola->erotation + molb->erotation;
  postcolln->evibration = mola->evibration + molb->evibration;

  // Finally, compute the proportion of energy left over for scattering

  postcolln->einternal    = postcolln->erotation + postcolln->evibration;
  postcolln->etranslation = precolln->etotal -  postcolln->einternal;

  if(postcolln->etranslation<0) {
     fprintf(stderr,
	"NEGATIVE ETRANSLATION IN Energy Disposal! Exiting!!\n");
     printf("NEGATIVE ETRANSLATION IN Energy Disposal! Exiting!!\n");
     printf("ET = %e, EV = %e, ER = %e, ETOT = %e \n", 
	    postcolln->etranslation , 
	    postcolln->erotation , postcolln->evibration,
	    precolln->etotal); 
     exit(0);
  }
  */
}

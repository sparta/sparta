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
#include "particle.h"
#include "mixture.h"
#include "collide.h"


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

  Particle::Species *species = particle->species;
  int nspecies = particle->nspecies;
  int imix = particle->find_mixture(mixID);


  double *vscale = particle->mixture[imix]->vscale;
  int *mix2group = particle->mixture[imix]->mix2group;


//  if (particle->nspecies != mixture->oldspecies) {
    memory->destroy(prefactor);
    memory->destroy(vremax);
//    memory->destroy(vrm);
    memory->create(prefactor,nspecies,nspecies,"collide:prefactor");
    memory->create(vremax,grid->nlocal,nspecies,nspecies,"collide:vremax");
    memory->create(vrm,nspecies,nspecies,"collide:vrm");
//  }

  // prefactor = static contributions to collision attempt frequencies

  double beta = 1.0;
  double cxs = 1.0;

  for (int isp = 0; isp < nspecies; isp++)
    for (int jsp = 0; jsp < nspecies; jsp++) {
      double diam = 0.5 * (species[isp].diam + species[jsp].diam);
      double omega = 0.5 * (species[isp].omega + species[jsp].omega);
      double tref = 0.5 * (species[isp].tref + species[jsp].tref);
      double mr = species[isp].mass * species[jsp].mass /
	     (species[isp].mass + species[jsp].mass);
             cxs = diam*diam*MY_PI;
             prefactor[isp][jsp] = cxs *
             pow(2.0*update->boltz*tref/mr,omega-0.5)/gamma(2.5-omega);
      double beta=MIN(vscale[isp],vscale[jsp]);
      double max_thermal_velocity = 3.0/beta;
             vrm[isp][jsp] = cxs * max_thermal_velocity;
    }

  // vremax = max relative velocity factors on per-grid, per-species basis
  
  int nglocal = grid->nlocal;

  for (int icell = 0; icell < nglocal; icell++)
    for (int isp = 0; isp < nspecies; isp++)
      for (int jsp = 0; jsp < nspecies; jsp++) {
// vremax value assignment should be done at a group level.
// each group should get the maximum vremax of all species involved.
        int igroup = mix2group[isp];
        int jgroup = mix2group[jsp];
	vremax[icell][igroup][jgroup] = MAX(vremax[icell][igroup][jgroup],vrm[isp][jsp]);
      }
}

/* ---------------------------------------------------------------------- */

double CollideVSS::attempt_collision(int icell, 
				     int igroup, int jgroup, double volume)
{
 double fnum = update->fnum;
 double dt = update->dt;

 int nattempt = 0.5 * ngroup[igroup] * (ngroup[jgroup]-1) *
    vremax[icell][igroup][jgroup] * dt * fnum / volume + random->uniform();
// printf(" Attempts = %d %d %d %d\n", icell, ngroup[igroup],ngroup[jgroup],nattempt );
// if (nattempt > 0) printf(" Attempts = %e %e %e \n", dt, fnum, volume);
/* no residual number of attempts is kept. The number of attempts is randomly scaled
   up or down. Better formulation for transient flows.
*/

  return nattempt;
}

/* ----------------------------------------------------------------------
   test if collision actually occurs
   also update vremax
------------------------------------------------------------------------- */

int CollideVSS::test_collision(int icell, int igroup, int jgroup,
			       Particle::OnePart *ip, Particle::OnePart *jp)
{
  Particle::Species *species = particle->species;
  Particle::OnePart *particles = particle->particles;

  if (random->uniform() < 0.5) return 1;
  return 0;

  double *vi = ip->v;
  double *vj = jp->v;
  int ispecies = ip->ispecies;
  int jspecies = jp->ispecies;
  double du  = vi[0] - vj[0];
  double dv  = vi[1] - vj[1];
  double dw  = vi[2] - vj[2];
  double vr2 = du*du + dv*dv + dw*dw;
  double omega1 = species[ispecies].omega;
  double omega2 = species[jspecies].omega;
  double omega = 0.5 * (omega1+omega2);
  double vr  = pow(vr2,1-omega);
// although the vremax is calcualted for the group the individual collisions
// calculated species dependent vre
  double vre = vr*prefactor[ispecies][jspecies];

// update vremax if needed
  vremax[icell][igroup][jgroup] = MAX(vre,vremax[icell][igroup][jgroup]);
  if (vre/vremax[icell][igroup][jgroup] < random->uniform()) return 0;
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
  precoln.erot  = ip->erot + jp->erot;
  precoln.evib  = ip->ivib * update->boltz*species[isp].vibtemp
                + jp->ivib * update->boltz*species[jsp].vibtemp;

  precoln.eint   = precoln.erot + precoln.evib;
  precoln.etotal = precoln.etrans + precoln.eint;

  postcoln.etrans = 0.0;
  postcoln.erot   = 0.0;
  postcoln.evib   = 0.0;
  postcoln.eint   = 0.0;
  postcoln.etotal = precoln.etotal;
}

/* ---------------------------------------------------------------------- */

Particle::OnePart *CollideVSS::perform_collision(Particle::OnePart *ip, 
						  Particle::OnePart *jp)
{
//  if (eng_exchange) EEXCHANGE_NonReactingEDisposal(ip,jp);
  EEXCHANGE_NonReactingEDisposal(ip,jp);
  SCATTER_TwoBodyScattering(ip,jp);
  return NULL;
}

/* ---------------------------------------------------------------------- */

void CollideVSS::SCATTER_TwoBodyScattering(Particle::OnePart *ip, 
					    Particle::OnePart *jp)
{
//  Particle::OnePart *particles = particle->particles;
  Particle::Species *species = particle->species;

  double *vi = ip->v;
  double *vj = jp->v;
  int isp = ip->ispecies;
  int jsp = jp->ispecies;
  double mass_i = species[isp].mass;
  double mass_j = species[jsp].mass;


  double alpha = 0.5 * (species[isp].alpha + species[jsp].alpha);
  double mr = species[isp].mass * species[jsp].mass /
	     (species[isp].mass + species[jsp].mass);

  postcoln.eint = ip->erot + jp->erot
                + ip->ivib * update->boltz*species[isp].vibtemp
                + jp->ivib * update->boltz*species[jsp].vibtemp;


  double cosX = 2.0*pow(random->uniform(), alpha) - 1.0;
  double sinX = sqrt(1.0 - cosX*cosX);
  double eps  = random->uniform() * 2*MY_PI;
  double vr   = sqrt(2.0 * postcoln.etrans / mr);

  double ua   = vr*cosX;
  double vb   = vr*sinX*cos(eps);
  double wc   = vr*sinX*sin(eps);

  // The centre_of_mass velocity is calculated using reactant masses

  double divisor = precoln.mass_i + precoln.mass_j;

  double ucmf =
      ((precoln.mass_i * vi[0]) +
       (precoln.mass_j * vj[0])) / divisor;
  double vcmf =
      ((precoln.mass_i * vi[1])+
       (precoln.mass_j * vj[1])) / divisor;
  double wcmf =
      ((precoln.mass_i * vi[2])+
       (precoln.mass_j * vj[2])) / divisor;

  // Compute new velocities for the products

//  printf("In Scatter Pre %e %e %e \n", vi[0], vi[1], vi[2]);
  divisor = mass_i + mass_j;
  vi[0] = ucmf - (mass_j/divisor)*ua;
  vi[1] = vcmf - (mass_j/divisor)*vb;
  vi[2] = wcmf - (mass_j/divisor)*wc;
  vj[0] = ucmf + (mass_i/divisor)*ua;
  vj[1] = vcmf + (mass_i/divisor)*vb;
  vj[2] = wcmf + (mass_i/divisor)*wc;

//  printf("In Scatter aft %e %e %e \n", vi[0], vi[1], vi[2]);
 return ;
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


  long i, Max_Level,Vib_state;

  // Handle each kind of EDisposal for NonReacing reactants

  if (precoln.ave_dof == 0)
  {
    ip->erot  = 0.0;
    jp->erot  = 0.0;
    ip->ivib = 0;
    jp->ivib = 0;
  }
  else 
  {
    double E_Dispose = precoln.etrans;
    for (i = 0; i < 2; i++) {
     if (i == 0) kp = ip; 
     if (i == 1) kp = jp; 
         int ksp = kp->ispecies;

         if (species[ksp].vibdof >= 2 && (vibn_phi>=random->uniform())) {
           evib = (double) (kp->ivib*update->boltz*species[ksp].vibtemp);
	   E_Dispose += evib;
	   Max_Level  = (long) 
	     (E_Dispose/(update->boltz * species[ksp].vibtemp));

	   do {
	     kp->ivib = (int) (random->uniform()*(Max_Level+AdjustFactor));
	     evib = (double) 
 	       (kp->ivib * update->boltz * species[ksp].vibtemp);
	     State_prob = pow((1 - evib / E_Dispose),
			      (1.5 - species[ksp].omega));
	   } while (State_prob < random->uniform());

	   E_Dispose -= evib;
         }

         if ((species[ksp].rotdof == 2) && (rotn_phi >= random->uniform())) {
           E_Dispose += kp->erot;
           Fraction_Rot = 
	     1- pow(random->uniform(),(1/(2.5-species[ksp].omega)));
	   kp->erot = Fraction_Rot * E_Dispose;
	   E_Dispose -= kp->erot;

         } else if ((species[ksp].rotdof > 2) && (rotn_phi >= random->uniform())) {
           E_Dispose += kp->erot;
	   Exp_1 = species[ksp].rotdof / 2;
	   Exp_2 = 2.5 - species[ksp].omega;
	   do {
	     Fraction_Rot = random->uniform();
	     State_prob   = 
	       ((Exp_1+Exp_2-2)/pow(((Exp_1-1)*Fraction_Rot),(Exp_1-1)) *
	       ((Exp_1+Exp_2)/(Exp_2-1)*pow((1-Fraction_Rot),(Exp_2-1))));
	   } while (State_prob < random->uniform());
	   kp->erot = Fraction_Rot*E_Dispose;
	   E_Dispose -= kp->erot;
         }
       } 
   }


  // Compute the post-collision internal energies

  postcoln.erot = ip->erot + jp->erot;
  int isp = ip->ispecies;
  int jsp = jp->ispecies;
  postcoln.evib = ip->ivib * update->boltz*species[isp].vibtemp 
                + jp->ivib * update->boltz*species[jsp].vibtemp;

  // Finally, compute the proportion of energy left over for scattering

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

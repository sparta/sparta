/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2012) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "react_qk.h"
#include "update.h"
#include "particle.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{DISSOCIATION,EXCHANGE,IONISATION,RECOMBINATION};

/* ---------------------------------------------------------------------- */

ReactQK::ReactQK(SPARTA *sparta, int narg, char **arg) :
  React(sparta, narg, arg)
{
  if (narg != 2) error->all(FLERR,"Illegal react qk command");
}

/* ---------------------------------------------------------------------- */

void ReactQK::init() 
{
  int nspecies = particle->nspecies;
  reactions = memory->create(reactions,nspecies,nspecies,"react/qk:reactions");
  nreactions = 0;

  /*
  for (i = 0; i < NSPECIES; i++)
  {
     for (j = 0; j < NSPECIES; j++)
     {
        reactions[i][j].formula             = NULL;
	reactions[i][j].outcome             = undefined;
	reactions[i][j].nproducts           = 0;
	reactions[i][j].product_C           = -1;
	reactions[i][j].product_D           = -1;
	reactions[i][j].product_E           = -1;
	reactions[i][j].e_forward           = 0.0;
	reactions[i][j].e_backward          = 0.0;
     }
  }
  */
}

/* ---------------------------------------------------------------------- */

int ReactQK::attempt(Particle::OnePart *ip, Particle::OnePart *jp, 
                     double pre_etrans, double pre_erot,
                     double pre_evib, double &post_etotal, int &kspecies)
{
  double rprobability,pre_etotal,distbn,prior,max_distbn;
  double fe,fc,fmax,lambda_v,cutoff,e_excess;
  double Teff,a,b,c,d,e,sqrta;
  double V_DoF,R_DoF,Omega_avg;
  double Exponent_1,Exponent_2;
  Reaction *reaction;

  Particle::Species *species = particle->species;
  int isp = ip->ispecies;
  int jsp = jp->ispecies;
  double mass_i = species[isp].mass;
  double mass_j = species[jsp].mass;

  double pre_rotdof_i = species[isp].rotdof;
  double pre_rotdof_j = species[jsp].rotdof;
  double pre_vibdof_i = species[isp].vibdof;
  double pre_vibdof_j = species[jsp].vibdof;
  double pre_ave_rotdof = (species[isp].rotdof + species[jsp].rotdof)/2.;
  double pre_ave_vibdof = (species[isp].vibdof + species[jsp].vibdof)/2.;
  double pre_ave_dof = (pre_ave_rotdof  + pre_ave_vibdof)/2.;

  // try both A+B and B+A reaction
  // if two species identical, only try reaction once

  int nattempt = 2;
  if (isp == jsp) nattempt = 1;

  for (int iattempt = 0; iattempt < nattempt; iattempt++) {
    if (iattempt == 0) reaction = &reactions[isp][jsp];
    else reaction = &reactions[jsp][isp];
    if (!reaction->formula) continue;

    pre_etotal = pre_etrans + pre_erot + pre_evib;
    e_excess = pre_etotal - reaction->e_forward;
    
    // crude test to filter out energetically impossible reactions
    
    if (e_excess <= 0.0) continue;

    // compute probability of reaction
    
    switch (reaction->outcome) {
    case IONISATION:
      rprobability = 0.0;
      if (rprobability > 1.0) {
        printf("rprobability > 1.0 in polyatomic test reaction!\n");
        printf("(rprob = %f)\n", rprobability);
      }
      break;

    case DISSOCIATION:
      {
        double ecc = pre_etrans + ip->ivib*update->boltz*species[isp].vibtemp;
        int maxlev = (int) (ecc/update->boltz*species[isp].vibtemp); 
        /*
        int limlev = update->boltz*species[isp].distemp /
          update->boltz*species[isp].vibtemp;
        if (maxlev > limlev) { 
          post_etotal = pre_etotal - 
            update->boltz*species[isp].distemp*update->boltz;
          ip->ispecies = reaction->product_C;
          jp->ispecies = reaction->product_D;
          kp->ispecies = reaction->product_D;
        } 
        */
        break;
      }

    case EXCHANGE:
      {
      // NOTE: what if neither case met?
        int asp,bsp;
        if (species[isp].vibdof > 0) {
          asp = isp;
          bsp = jsp;
        } else if (species[jsp].vibdof > 0) {
          asp = jsp;
          bsp = isp;
        }
        //double omega = 0.5 * (params[asp].omega + params[jsp].omega);

        // pre-collision molecule that splits is of species asp

        int reac_vmode = 1;

        // endothermic reaction

        /*
        if (reaction->e_forward < 0.0 && reac_vmode > 0.0) {
          double ecc = pre_etrans + ap->ivib*update->boltz*species[asp].vibtemp;
          int maxlev = (int) (ecc/update->boltz*species[asp].vibtemp);
          if (ecc > abs(reaction->e_forward)) {
            double prob = 0.0;
            do {
              iv=random->uniform()*(maxlev+0.99999999);
              evib= (double) iv*update->boltz*species[asp].vibtemp;
              if (evib < ecc) prob = pow(1.0-evib/ecc,1.5-omega);
              //PROB is the probability ratio of eqn (5.61)
            } while (random->uniform < prob);

            ilevel = (int) abs(reaction->e_forward/update->boltz/
                               species[asp].vibtemp);
          
            if (iv >= ilevel) { // GT MODEL
              ip->species = reaction->product_C;
              jp->species = reaction->product_D;
              post_etotal = pre_etotal + reaction->e_forward;
              ec = post_etotal;
            }
          }
          // NOTE: any additional vibrational modes must be set to zero
        } // NOTE: else for endothermic?

        // exothermic reaction
        
        //double omega = 0.5 * (params[asp].omega + params[jsp].omega);
        double ecc = pre_etrans + ap->ivib*update->boltz*species[asp].vibtemp;
        
        // the potential post-collision species of the molecule
        int mspec = reaction->product_C; 
        // the potential post-collision species of the atom
        int nspec = reaction->product_D; 
        // potential post-collision energy
        ecc += reaction->e_forward;
        
        do {
          iv=random->uniform()*(maxlev+0.99999999);
          evib= (double) iv*update->boltz*species[mspec].vibtemp;
          if (evib < ecc) prob = pow(1.0-evib/ecc,1.5-omega);
          //PROB is the probability ratio of eqn (5.61)
        } while (random->uniform < prob);
          
        ilevel = (int) abs(reaction->e_forward/update->boltz/
                           species[mspec].vibtemp);

        if (iv >= ilevel) { // GT MODEL
          ip->species = reaction->product_C;
          jp->species = reaction->product_D;
          post_etotal = pre_etotal+reaction->e_forward;
          ec = postcoln_etotal;
        }
        */        
        break;
      }

    default:
      error->one(FLERR,"Unknown outcome in reaction");
      break;
    }
  }

  return 0;
}

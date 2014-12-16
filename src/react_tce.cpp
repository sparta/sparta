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
#include "react_tce.h"
#include "particle.h"
#include "collide.h"
#include "random_park.h"
#include "error.h"

using namespace SPARTA_NS;

enum{DISSOCIATION,EXCHANGE,IONIZATION,RECOMBINATION};   // other files

/* ---------------------------------------------------------------------- */

ReactTCE::ReactTCE(SPARTA *sparta, int narg, char **arg) :
  ReactBird(sparta, narg, arg) {}

/* ---------------------------------------------------------------------- */

void ReactTCE::init()
{
  if (!collide || strcmp(collide->style,"vss") != 0)
    error->all(FLERR,"React tce can only be used with collide vss");

  ReactBird::init();
}

/* ---------------------------------------------------------------------- */

int ReactTCE::attempt(Particle::OnePart *ip, Particle::OnePart *jp, 
                      double pre_etrans, double pre_erot, double pre_evib,
                      double &post_etotal, int &kspecies)
{
  double pre_etotal,ecc,e_excess;
  OneReaction *r;

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
  double pre_ave_dof = 0.5 * (pre_ave_rotdof + pre_ave_vibdof);

  int n = reactions[isp][jsp].n;
              
  if (n == 0) return 0;
  int *list = reactions[isp][jsp].list;

  // probablity to compare to reaction probability

  double react_prob = 0.0;
  double random_prob = random->uniform(); 

  // loop over possible reactions for these 2 species

  for (int i = 0; i < n; i++) {
    r = &rlist[list[i]];

    // ignore energetically impossible reactions

    pre_etotal = pre_etrans + pre_erot + pre_evib;

    ecc = pre_etrans; 
    if (pre_ave_rotdof > 0.1) ecc += pre_erot*r->coeff[0]/pre_ave_rotdof;

    e_excess = ecc - r->coeff[1];
    if (e_excess <= 0.0) continue;
        
    // compute probability of reaction
        
    switch (r->type) {
    case DISSOCIATION:
    case EXCHANGE:
      {
        react_prob += r->coeff[2] * 
          pow(ecc-r->coeff[1],r->coeff[3]) *
          pow(1.0-r->coeff[1]/ecc,r->coeff[5]);
        break;
      }
        
    default:
      error->one(FLERR,"Unknown outcome in reaction");
      break;
    }
      
    // test against random number to see if this reaction occurs

    if (react_prob > random_prob) {
      ip->ispecies = r->products[0];
      jp->ispecies = r->products[1];
      
      post_etotal = pre_etotal + r->coeff[4];
      if (r->nproduct > 2) kspecies = r->products[2];
      else kspecies = -1;
      
      return 1;
    }
  }
  
  return 0;
}

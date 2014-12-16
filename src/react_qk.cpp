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
#include "react_qk.h"
#include "update.h"
#include "particle.h"
#include "collide.h"
#include "random_park.h"
#include "math_const.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{DISSOCIATION,EXCHANGE,IONIZATION,RECOMBINATION};   // other files

/* ---------------------------------------------------------------------- */

ReactQK::ReactQK(SPARTA *sparta, int narg, char **arg) :
  ReactBird(sparta, narg, arg) {}

/* ---------------------------------------------------------------------- */

void ReactQK::init()
{
  if (!collide || strcmp(collide->style,"vss") != 0)
    error->all(FLERR,"React qk can only be used with collide vss");

  ReactBird::init();
}

/* ---------------------------------------------------------------------- */

int ReactQK::attempt(Particle::OnePart *ip, Particle::OnePart *jp, 
                     double pre_etrans, double pre_erot, double pre_evib,
                     double &post_etotal, int &kspecies)
{
  double pre_etotal,ecc,e_excess;
  double reac_prob,prob,evib;
  int iv,ilevel,maxlev,limlev;
  int mspec,aspec;
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

  double iomega = collide->extract(isp,"omega");
  double jomega = collide->extract(jsp,"omega");
  double omega = 0.5 * (iomega+jomega);

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
      {
        ecc = pre_etrans + ip->evib;
        maxlev = static_cast<int> (ecc/(update->boltz*species[isp].vibtemp));
        limlev = static_cast<int> 
          (fabs(r->coeff[1])/(update->boltz*species[isp].vibtemp));
        if (maxlev > limlev) reac_prob = 1.0;
        break; 
      }
    case EXCHANGE:
      {
        if (r->coeff[1] > 0.0 && species[isp].rotdof > 0) {

          // endothermic reaction 

          ecc = pre_etrans + ip->evib;
          maxlev = static_cast<int> (ecc/(update->boltz*species[isp].vibtemp));
          if (ecc > r->coeff[1]) {

            // PROB is the probability ratio of eqn (5.61)

            prob = 0.0;
            do {
              iv =  static_cast<int> (random->uniform()*(maxlev+0.99999999));
              evib = static_cast<double> 
                (iv*update->boltz*species[isp].vibtemp);
              if (evib < ecc) reac_prob = pow(1.0-evib/ecc,1.5-omega);
            } while (random->uniform() < reac_prob);
            
            ilevel = static_cast<int> 
              (abs(fabs(r->coeff[4]))/(update->boltz*species[isp].vibtemp));

            // GT MODEL

            if (iv >= ilevel) reac_prob = 1.0;
          }

        } else if (species[isp].rotdof > 0) {

          ecc = pre_etrans + ip->evib;
        
          // the potential post-collision species of the particle
          mspec = r->products[0];
          // the potential post-collision species of the atom
          aspec = r->products[1];

          if (species[mspec].rotdof < 2.0)  {
            mspec = r->products[1];
            aspec = r->products[0];
          }
              
          // potential post-collision energy

          ecc += r->coeff[4];
          maxlev = static_cast<int> (ecc/(update->boltz*species[isp].vibtemp));
          do {
            iv = random->uniform()*(maxlev+0.99999999);
            evib = static_cast<double> 
              (iv*update->boltz*species[mspec].vibtemp);
            if (evib < ecc) prob = pow(1.0-evib/ecc,1.5-omega);
          } while (random->uniform() < prob);

          ilevel = static_cast<int> 
            (fabs(r->coeff[1]/update->boltz/species[mspec].vibtemp));

          // GT MODEL

          if (iv >= ilevel) reac_prob = 1.0;
        }
        
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

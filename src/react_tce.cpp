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
#include "grid.h"
#include "modify.h"
#include "compute.h"

// DEBUG
#include "update.h"



using namespace SPARTA_NS;

enum{NONE,DISCRETE,SMOOTH};
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
  double pre_etotal,ecc,e_excess,z,Tt;
  int imode;
  OneReaction *r;
   
  Grid::ChildInfo *cinfo = grid->cinfo;
  Particle::OnePart *particles = particle->particles;
  int *next = particle->next;

  Particle::Species *species = particle->species;
  int isp = ip->ispecies;
  int jsp = jp->ispecies;
  
  double pre_ave_rotdof = (species[isp].rotdof + species[jsp].rotdof)/2.0;

  int n = reactions[isp][jsp].n;
  if (n == 0) return 0;
  int *list = reactions[isp][jsp].list;
  double nfrac_p3;

  // probablity to compare to reaction probability

  double react_prob = 0.0;
  double random_prob = random->uniform(); 

  // loop over possible reactions for these 2 species

  for (int i = 0; i < n; i++) {
    r = &rlist[list[i]];

    // ignore energetically impossible reactions

    pre_etotal = pre_etrans + pre_erot + pre_evib;

    // two options for total energy in TCE model
    // 1: rDOF model
    // 2: TCE: Rotation + Vibration
      
    // average DOFs participating in the reaction
             
    if (birdflag) {
       ecc = pre_etrans;
       z = 0.0;
       if (r->type == DISSOCIATION) z = 1.0;
       if (pre_ave_rotdof > 0.1) ecc += pre_erot*z/pre_ave_rotdof;
    }
    else {
       ecc = pre_etotal;
       if (pre_etotal+r->coeff[4] <= 0.0) continue; // Cover cases where coeff[1].neq.coeff[4]
       z = pre_ave_rotdof + (species[isp].vibdof + species[jsp].vibdof)/2.0;
    }

    e_excess = ecc - r->coeff[1];
    if (e_excess <= 0.0) continue;
      
//    if (!birdflag) z += (species[isp].vibdof + species[jsp].vibdof)/2.0;
//    {
//        if (collide->vibstyle == SMOOTH) z += (species[isp].vibdof + species[jsp].vibdof)/2.0;
//        else if (collide->vibstyle == DISCRETE) {
////            Tt = pre_etrans / (update->boltz * (2.5-r->coeff[5]));
//            Tt = modify->compute[1]->vector_grid[0];
//            if (species[isp].vibdof == 2) z += (species[isp].vibtemp[0]/Tt) / (exp(species[isp].vibtemp[0]/Tt)-1);
//            else if (species[isp].vibdof > 2) {
//                imode = 0;
//                while (imode < 4) z += (species[isp].vibtemp[imode]/Tt) / (exp(species[isp].vibtemp[imode]/Tt)-1);
//            }
//            if (species[jsp].vibdof == 2) z += (species[jsp].vibtemp[0]/Tt) / (exp(species[jsp].vibtemp[0]/Tt)-1);
//            else if (species[jsp].vibdof > 2) {
//                imode = 0;
//                while (imode < 4) z += (species[jsp].vibtemp[imode]/Tt) / (exp(species[jsp].vibtemp[imode]/Tt)-1);
//            }
//        }
//    }

        
    // compute probability of reaction
        
    switch (r->type) {
    case DISSOCIATION:
    case IONIZATION:
    case EXCHANGE:
      {
        react_prob += r->coeff[2] * tgamma(z+2.5-r->coeff[5]) / MAX(1.0e-6,tgamma(z+r->coeff[3]+1.5)) *
          pow(ecc-r->coeff[1],r->coeff[0]) *
          pow(1.0-r->coeff[1]/ecc,z + 1.5 - r->coeff[5]);
//        react_prob += r->coeff[2] *
//          pow(ecc-r->coeff[1],r->coeff[3]) *
//          pow(1.0-r->coeff[1]/ecc,r->coeff[5]);
        break;
      }

    case RECOMBINATION:
      {
        // skip if no 3rd particle chosen by Collide::collisions()
        //   this includes effect of boost factor to skip recomb reactions
        // check if this recomb reaction is the same one 
        //   that the 3rd particle species maps to, else skip it
        // this effectively skips all recombinations reactions 
        //   if selected a 3rd particle species that matches none of them
        // scale probability by boost factor to restore correct stats

        if (recomb_species < 0) continue;
        int *sp2recomb = reactions[isp][jsp].sp2recomb;
        if (sp2recomb[recomb_species] != list[i]) continue;
         
        //  In this TCE implemenation recomb_density is taken as the
        //  3rd body number density instead of the total number density
         
        int np3 = 0;
        int icell = ip->icell;
        int np = cinfo[icell].count;
        int ip3 = cinfo[icell].first;
        while (ip3 >= 0) {
          if (particles[ip3].ispecies == react->recomb_species) np3 += 1;
          ip3 = next[ip3];
        }
        if (np3 > 0) nfrac_p3 = np3/double(np);
        else nfrac_p3 = 0.0;

        react_prob += recomb_boost * recomb_density * r->coeff[2] * tgamma(z+2.5-r->coeff[5]) / MAX(1.0e-6,tgamma(z+r->coeff[3]+1.5)) *
          pow(ecc-r->coeff[1],r->coeff[0]) *  // extended to general recombination case with non-zero activation energy
          pow(1.0-r->coeff[1]/ecc,z + 1.5 - r->coeff[5]);
//        react_prob += recomb_boost * recomb_density * r->coeff[2] *
//          pow(ecc-r->coeff[1],r->coeff[3]) *
//          pow(1.0-r->coeff[1]/ecc,r->coeff[5]);
        break;
      }

    default:
      error->one(FLERR,"Unknown outcome in reaction");
      break;
    }
      
    // test against random number to see if this reaction occurs
    // if it does, reset species of I,J and optional K to product species
    // J particle can be destroyed in recombination reaction, set species = -1
    // K particle can be created in a dissociation or ionization reaction,
    //   set its kspecies, parent will create it
    // important NOTE:
    //   does not matter what order I,J reactants are in compared
    //     to order the reactants are listed in the reaction file
    //   for two reasons:
    //   a) list of N possible reactions above includes all reactions
    //      that I,J species are in, regardless of order
    //   b) properties of pre-reaction state, stored in precoln,
    //      as computed by setup_collision(),
    //      and used by perform_collision() after reaction has taken place,
    //      only store combined properties of I,J,
    //      nothing that is I-specific or J-specific

    if (react_prob > random_prob) {
//      fprintf(screen,"%d %d z = %f Tt = %f\n",r->reactants[0],r->reactants[1],z,Tt);
      tally_reactions[list[i]]++;
      ip->ispecies = r->products[0];

      // Previous statment did not destroy the 2nd species (B) if
      //   recombination was specified as A+B->AB+M (which has nproductus=2)
      //   but only for the A+B->AB specication form (which has nproductus=1)

      switch (r->type) {
      case DISSOCIATION:
      case IONIZATION:
      case EXCHANGE:
        {
          jp->ispecies = r->products[1];
          break;
        }
      case RECOMBINATION:
        {
          // always destroy 2nd reactant species

          jp->ispecies = -1;
          break;
        }
      }

      if (r->nproduct > 2) kspecies = r->products[2];
      else kspecies = -1;

      post_etotal = pre_etotal + r->coeff[4];

      return 1;
    }
  }
  
  return 0;
}

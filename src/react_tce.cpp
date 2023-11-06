/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
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
#include "react_tce.h"
#include "particle.h"
#include "collide.h"
#include "update.h"
#include "random_knuth.h"
#include "error.h"

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
  double pre_etotal,ecc,e_excess,z;
  int inmode,jnmode;
  OneReaction *r;

  Particle::Species *species = particle->species;
  int isp = ip->ispecies;
  int jsp = jp->ispecies;
  double ievib = ip->evib;
  double jevib = jp->evib;

  double pre_ave_rotdof = (species[isp].rotdof + species[jsp].rotdof)/2.0;

  int n = reactions[isp][jsp].n;
  if (n == 0) return 0;
  int *list = reactions[isp][jsp].list;

  // probablity to compare to reaction probability

  double react_prob = 0.0;
  double random_prob = random->uniform();
  double zi = 0.0;
  double zj = 0.0;
  int avei = 0;
  int avej = 0;
  double iTvib = 0.0;
  double jTvib = 0.0;

  // loop over possible reactions for these 2 species

  for (int i = 0; i < n; i++) {
    r = &rlist[list[i]];

    // ignore energetically impossible reactions

    pre_etotal = pre_etrans + pre_erot + pre_evib;

    // two options for total energy in TCE model
    // 0: partialEnergy = true: rDOF model
    // 1: partialEnergy = false: TCE: Rotation + Vibration

    // average DOFs participating in the reaction

    if (partialEnergy) {
       ecc = pre_etrans;
       z = r->coeff[0];
       if (pre_ave_rotdof > 0.1) ecc += pre_erot*z/pre_ave_rotdof;
    } else {
       ecc = pre_etotal;
       z = pre_ave_rotdof;
    }

    // Cover cases where coeff[1].neq.coeff[4]
    if (r->coeff[1]>((-1)*r->coeff[4])) e_excess = ecc - r->coeff[1];
    else e_excess = ecc + r->coeff[4];
    if (e_excess <= 0.0) continue;


    if (!partialEnergy) {

       if (collide->vibstyle == SMOOTH) z += (species[isp].vibdof + species[jsp].vibdof)/2.0;
       else if (collide->vibstyle == DISCRETE) {
            inmode = species[isp].nvibmode;
            jnmode = species[jsp].nvibmode;
            //Instantaneous z for diatomic molecules
            if (inmode == 1) {
                avei = static_cast<int>
                        (ievib / (update->boltz * species[isp].vibtemp[0]));
                if (avei > 0) zi = 2.0 * avei * log(1.0 / avei + 1.0);
                else zi = 0.0;
            } else if (inmode > 1) {
                if (ievib < 1e-26 ) zi = 0.0; //Low Energy Cut-Off to prevent nan solutions to newtonTvib
                //Instantaneous T for polyatomic
                else {
                  iTvib = newtonTvib(inmode,ievib,species[isp].vibtemp,3000,1e-4,1000);
                  zi = (2 * ievib)/(update->boltz * iTvib);
                }
            } else zi = 0.0;

            if (jnmode == 1) {
                avej = static_cast<int>
                        (jevib / (update->boltz * species[jsp].vibtemp[0]));
                if (avej > 0) zj = 2.0 * avej * log(1.0 / avej + 1.0);
                else zj = 0.0;
            } else if (jnmode > 1) {
                if (jevib < 1e-26) zj = 0.0;
                else {
                  jTvib = newtonTvib(jnmode,jevib,species[jsp].vibtemp,3000,1e-4,1000);
                  zj = (2 * jevib)/(update->boltz * jTvib);
                }
            } else zj = 0.0;

            if (isnan(zi) || isnan(zj) || zi < 0 || zj < 0) error->one(FLERR,"Root-Finding Error");
            z += 0.5 * (zi+zj);
       }
    }

    // compute probability of reaction

    switch (r->type) {
    case DISSOCIATION:
    case IONIZATION:
    case EXCHANGE:
      {
        react_prob += r->coeff[2] * tgamma(z+2.5-r->coeff[5]) / MAX(1.0e-6,tgamma(z+r->coeff[3]+1.5)) *
          pow(ecc-r->coeff[1],r->coeff[3]-1+r->coeff[5]) *
          pow(1.0-r->coeff[1]/ecc,z+1.5-r->coeff[5]);
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

        react_prob += recomb_boost * recomb_density * r->coeff[2] *
          tgamma(z+2.5-r->coeff[5]) / MAX(1.0e-6,tgamma(z+r->coeff[3]+1.5)) *
          pow(ecc-r->coeff[1],r->coeff[3]-1+r->coeff[5]) *  // extended to general recombination case with non-zero activation energy
          pow(1.0-r->coeff[1]/ecc,z+1.5-r->coeff[5]);
        break;
      }

    if (react_prob < 0) error->warning(FLERR,"Negative reaction probability");
    else if (react_prob > 1) error->warning(FLERR,"Reaction probability greater than 1");

    default:
      error->one(FLERR,"Unknown outcome in reaction");
      break;
    }

    // test against random number to see if this reaction occurs
    // if it does, reset species of I,J and optional K to product species
    // J particle is destroyed in recombination reaction, set species = -1
    // K particle can be created in a dissociation or ionization reaction,
    //   set its kspecies, parent will create it
    // important NOTE:
    //   does not matter what order I,J reactants are in compared
    //     to order the reactants are listed in the reaction file
    //   for two reasons:
    //   a) list of N possible reactions above includes all reactions
    //      that I,J species are in, regardless of order
    //   b) properties of pre-reaction state are stored in precoln:
    //      computed by setup_collision()
    //      used by perform_collision() after reaction has taken place
    //      precoln only stores combined properties of I,J
    //      nothing that is I-specific or J-specific

    if (react_prob > random_prob) {
      tally_reactions[list[i]]++;

      if (!computeChemRates) {
          ip->ispecies = r->products[0];

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
      }

      return 1;
    }
  }

  return 0;
}

/* ---------------------------------------------------------------------- */

double ReactTCE::bird_Evib(int nmode, double Tvib,
                            double vibtemp[],
                            double Evib)
{

  // COMPUTES f FOR NEWTON'S SEARCH METHOD OUTLINED IN "newtonTvib".

  double f = -Evib;
  double kb = 1.38064852e-23;

  for (int i = 0; i < nmode; i++) {
    f += (((kb*vibtemp[i])/(exp(vibtemp[i]/Tvib)-1)));
  }

  return f;

}


/* ---------------------------------------------------------------------- */

double ReactTCE::bird_dEvib(int nmode, double Tvib, double vibtemp[])
{

  // COMPUTES df FOR NEWTON'S SEARCH METHOD

  double df = 0.0;
  double kb = 1.38064852e-23;

  for (int i = 0; i < nmode; i++) {
    df += ((pow(vibtemp[i],2)*kb*exp(vibtemp[i]/Tvib))/(pow(Tvib,2)*pow(exp(vibtemp[i]/Tvib)-1,2)));
  }

  return df;

}

/* ---------------------------------------------------------------------- */

double ReactTCE::newtonTvib(int nmode, double Evib, double vibTemp[],
               double Tvib0,
               double tol,
               int nmax)
{


  // FUNCTION FOR CONVERTING VIBRATIONAL ENERGY TO VIBRATIONAL TEMPERATURE
  // Computes Tvib assuming the vibrational energy levels occupy a simple harmonic oscillator (SHO)
  // spacing.
  // Search for Tvib begins at some initial value "Tvib0" until the search reaches a tolerance level "tol".


  double f;
  double df;
  double Tvib, Tvib_prev;
  double err;
  int i;

  // Uses Newton's method to solve for a vibrational temperature given a
  // distribution of vibrational energy levels.

  // f and df are computed for Newton's search
  f = bird_Evib(nmode,Tvib0,vibTemp,Evib);
  df = bird_dEvib(nmode,Tvib0,vibTemp);

  // Update guess for Tvib and compute error
  Tvib = Tvib0 - (f/df);
  err = fabs(Tvib-Tvib0);

  i=2;

  // Continue to search for Tvib until the error is greater than the tolerance:
  while((err >= tol) && (i <= nmax))
  {

    Tvib_prev = Tvib;

    f = bird_Evib(nmode,Tvib,vibTemp,Evib);
    df = bird_dEvib(nmode,Tvib,vibTemp);

    Tvib = Tvib_prev-(f/df);
    err = fabs(Tvib-Tvib_prev);

    i=i+1;

  }

  return Tvib;

}

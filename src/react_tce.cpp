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
  ReactBird(sparta, narg, arg)
{
  prob_warn_flag = 0;
}

/* ---------------------------------------------------------------------- */

void ReactTCE::init()
{
  if (!collide || strcmp(collide->style,"vss") != 0)
    error->all(FLERR,"React tce can only be used with collide vss");

  prob_warn_flag = 0;

  ReactBird::init();

  // vib_energy micro: validate settings and build the per-reaction
  // microcanonical energy-factor tables

  if (vibEnergyMode == VIB_MICRO) {
    if (partialEnergy)
      error->all(FLERR,"react_modify vib_energy micro requires "
                 "partial_energy no");
    if (collide->vibstyle != DISCRETE)
      error->all(FLERR,"react_modify vib_energy micro requires "
                 "collide_modify vibrate discrete");
    if (elecEnergyMode == ELEC_INCLUDE)
      error->all(FLERR,"react_modify vib_energy micro cannot be combined "
                 "with elec_energy yes");
    build_micro_tables();
  } else free_micro_tables();
}

/* ---------------------------------------------------------------------- */

int ReactTCE::attempt(Particle::OnePart *ip, Particle::OnePart *jp,
                      double pre_etrans, double pre_erot, double pre_evib, double pre_eelec,
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

    pre_etotal = pre_etrans + pre_erot + pre_evib + pre_eelec;

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
       // input Arrhenius rates are equilibrium-calibrated, so the excited-
       // state contribution is already contained in (A,eta,Ea). ELEC_MICRO
       // (the default) keeps eelec in ecc and replaces the TCE energy factor
       // by its microcanonical ladder average (see below), which keeps the
       // equilibrium rate on the input Arrhenius rate. ELEC_EXCLUDE instead
       // drops eelec from ecc entirely; ELEC_INCLUDE adds it with per-state
       // DOF (historical behavior), which raises the equilibrium rate above
       // the input rate. eelec stays in pre_etotal/post_etotal either way so
       // energy is conserved.
       if (elecEnergyMode == ELEC_EXCLUDE) ecc -= pre_eelec;
       z = pre_ave_rotdof;
    }

    // Cover cases where coeff[1].neq.coeff[4]

    if (r->coeff[1]>((-1)*r->coeff[4])) e_excess = ecc - r->coeff[1];
    else e_excess = ecc + r->coeff[4];
    if (e_excess <= 0.0) continue;


    // with vib_energy micro the vibrational (and electronic) ladders enter
    // through the microcanonical energy-factor table instead of z:
    // z stays at the rotational pair average and the instantaneous-DOF
    // blocks below are skipped

    if (!partialEnergy && vibEnergyMode != VIB_MICRO) {

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

      // per-state electronic DoF only participates when electronic energy
      // is included in the reaction energy (react_modify elec_energy yes)

      if (collide->elecstyle == DISCRETE && elecEnergyMode == ELEC_INCLUDE) {
        zi = 0.0;
        if (species[isp].elecdat != NULL) {
          int ielec = particle->eivec[particle->ewhich[collide->index_elecstate]][ip - particle->particles];
          zi = species[isp].elecdat->states[ielec].dof;
        }
        zj = 0.0;
        if (species[jsp].elecdat != NULL) {
          int ielec = particle->eivec[particle->ewhich[collide->index_elecstate]][jp - particle->particles];
          zj = species[jsp].elecdat->states[ielec].dof;
        }
        z += 0.5*(zi + zj);
      }
    }

    // energy-dependent factor of the TCE probability
    // standard form: (ecc-Ea)^(eta-1+omega) * (1-Ea/ecc)^(z+1.5-omega)
    // for elec_energy micro, ecc is the TOTAL collision energy including
    // electronic, and the factor is replaced by its microcanonical average
    // over the pair's electronic ladder at fixed total energy:
    //   sum_p g_p (ecc-eps_p-Ea)_+^(z+eta+0.5)
    //     / sum_p g_p (ecc-eps_p)_+^(z+1.5-omega)
    // (identical to the standard form when neither species has electronic
    // states). This lets electronic energy count toward the barrier while
    // keeping the equilibrium rate on the input Arrhenius rate.

    double efactor;
    if (!partialEnergy && vibEnergyMode == VIB_MICRO &&
        mtab && mtab[list[i]])
      efactor = vib_micro_factor(list[i],ecc);
    else if (!partialEnergy && elecEnergyMode == ELEC_MICRO &&
        collide->elecstyle == DISCRETE &&
        (species[isp].elecdat != NULL || species[jsp].elecdat != NULL))
      efactor = elec_micro_factor(isp,jsp,ecc,z,r);
    else
      efactor = pow(ecc-r->coeff[1],r->coeff[3]-1+r->coeff[5]) *
                pow(1.0-r->coeff[1]/ecc,z+1.5-r->coeff[5]);

    // compute probability of reaction

    switch (r->type) {
    case DISSOCIATION:
    case IONIZATION:
    case EXCHANGE:
      {
        react_prob += r->coeff[2] * tgamma(z+2.5-r->coeff[5]) / MAX(1.0e-6,tgamma(z+r->coeff[3]+1.5)) *
          efactor;
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
          efactor;   // extended to general recombination case with non-zero activation energy
        break;
      }

    default:
      error->one(FLERR,"Unknown outcome in reaction");
      break;
    }

    // sum of reaction probabilities should be < 1 for a valid TCE scheme,
    //   else reaction rates are biased by clipping
    // warn only once per run to avoid flooding output

    if (!prob_warn_flag) {
      if (react_prob < 0.0) {
        prob_warn_flag = 1;
        error->warning(FLERR,"Negative TCE reaction probability, "
                       "check reaction file coefficients "
                       "(further warnings suppressed)");
      } else if (react_prob > 1.0) {
        prob_warn_flag = 1;
        error->warning(FLERR,"TCE reaction probability exceeded 1.0, "
                       "chemistry may be under-resolved, "
                       "consider reducing timestep or fnum "
                       "(further warnings suppressed)");
      }
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

        // return reaction from 1 to N

        return list[i] + 1;
      }

      break;
    }
  }

  // no reaction performed

  return 0;
}

/* ----------------------------------------------------------------------
   energy factor of the TCE probability for react_modify elec_energy micro
   ecc = TOTAL collision energy including the pair's electronic energy
   returns the microcanonical (density-of-states weighted) average of the
   standard TCE energy factor over the two reactants' electronic ladders:
     sum_p g_p (ecc - eps_p - Ea)_+^(z+eta+0.5)
       / sum_p g_p (ecc - eps_p)_+^(z+1.5-omega)
   where p runs over pair states (eps_p = eps_i + eps_j, g_p = g_i g_j).
   For a single pair state at eps = 0 this reduces to
     (ecc-Ea)^(eta-1+omega) * (1-Ea/ecc)^(z+1.5-omega).
   The average is chosen so the equilibrium reaction rate stays on the
   input Arrhenius rate while electronic energy counts toward the barrier.
------------------------------------------------------------------------- */

double ReactTCE::elec_micro_factor(int isp, int jsp, double ecc, double z,
                                   OneReaction *r)
{
  Particle::Species *species = particle->species;
  double boltz = update->boltz;

  static const Particle::ElecState ground = {0.0, 1, 1, 0.0};

  const Particle::ElecState *istates,*jstates;
  int ni,nj;

  if (species[isp].elecdat) {
    istates = species[isp].elecdat->states;
    ni = species[isp].elecdat->nelecstate;
  } else {
    istates = &ground;
    ni = 1;
  }
  if (species[jsp].elecdat) {
    jstates = species[jsp].elecdat->states;
    nj = species[jsp].elecdat->nelecstate;
  } else {
    jstates = &ground;
    nj = 1;
  }

  double ea = r->coeff[1];
  double exp_num = z + r->coeff[3] + 0.5;
  double exp_den = z + 1.5 - r->coeff[5];

  double num = 0.0;
  double den = 0.0;

  for (int i = 0; i < ni; i++) {
    for (int j = 0; j < nj; j++) {
      double x = ecc - boltz*(istates[i].temp + jstates[j].temp);
      if (x <= 0.0) continue;
      double g = istates[i].degen * jstates[j].degen;
      den += g * pow(x,exp_den);
      if (x > ea) num += g * pow(x-ea,exp_num);
    }
  }

  if (den == 0.0) return 0.0;
  return num/den;
}

/* ---------------------------------------------------------------------- */

double ReactTCE::bird_Evib(int nmode, double Tvib,
                            double vibtemp[],
                            double Evib)
{
  // Comutes f for Newton's search method outlined in newtonTvib()

  double f = -Evib;
  double kb = update->boltz;

  for (int i = 0; i < nmode; i++) {
    const double vti = vibtemp[i];
    f += (((kb*vti)/(exp(vti/Tvib)-1)));
  }

  return f;
}

/* ---------------------------------------------------------------------- */

double ReactTCE::bird_dEvib(int nmode, double Tvib, double vibtemp[])
{
  // Comutes df for Newton's search method

  double df = 0.0;
  double kb = update->boltz;

  for (int i = 0; i < nmode; i++) {
    const double vti = vibtemp[i];
    const double vti2 = vti * vti;
    const double Tvib2 = Tvib * Tvib;
    const double k1 = vti/Tvib;
    const double ek1 = exp(k1);
    const double k2 = ek1 - 1.0;
    const double k22 = k2 * k2;
    df += (vti2*kb*ek1)/(Tvib2*k22);
  }

  return df;
}

/* ---------------------------------------------------------------------- */

double ReactTCE::newtonTvib(int nmode, double Evib, double vibTemp[],
               double Tvib0,
               double tol,
               int nmax)
{
  // Function for converting vibrational energy to vibrational temperature
  // Computes Tvib assuming the vibrational energy levels occupy a simple harmonic oscillator (SHO) spacing
  // Search for Tvib begins at some initial value "Tvib0" until the search reaches a tolerance level "tol"

  double f;
  double df;
  double Tvib, Tvib_prev;
  double err;
  int i;

  // Uses Newton's method to solve for a vibrational temperature given a
  // distribution of vibrational energy levels

  // f and df are computed for Newton's search

  f = bird_Evib(nmode,Tvib0,vibTemp,Evib);
  df = bird_dEvib(nmode,Tvib0,vibTemp);

  // Update guess for Tvib and compute error

  Tvib = Tvib0 - (f/df);
  err = fabs(Tvib-Tvib0);

  i = 2;

  // Continue to search for Tvib until the error is greater than the tolerance:

  while((err >= tol) && (i <= nmax))
  {
    Tvib_prev = Tvib;

    f = bird_Evib(nmode,Tvib,vibTemp,Evib);
    df = bird_dEvib(nmode,Tvib,vibTemp);

    Tvib = Tvib_prev-(f/df);
    err = fabs(Tvib-Tvib_prev);

    i++;
  }

  return Tvib;
}

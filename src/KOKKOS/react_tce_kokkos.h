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

#ifdef REACT_CLASS

ReactStyle(tce/kk,ReactTCEKokkos)

#else

#ifndef SPARTA_REACT_TCE_KOKKOS_H
#define SPARTA_REACT_TCE_KOKKOS_H

#include "math.h"
#include "react_bird_kokkos.h"
#include "kokkos_type.h"
#include "update.h"

namespace SPARTA_NS {

class ReactTCEKokkos : public ReactBirdKokkos {
 public:
  ReactTCEKokkos(class SPARTA *, int, char **);
  ReactTCEKokkos(class SPARTA* sparta) : ReactBirdKokkos(sparta) {};
  void init();
  int attempt(Particle::OnePart *, Particle::OnePart *,
              double, double, double, double &, int &) { return 0; }

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
double bird_Evib(const int& nmode, const double& Tvib,
                 const double vibtemp[],
                 const double& Evib) const
{
  // Comutes f for Newton's search method outlined in newtonTvib()

  double f = -Evib;
  const double kb = 1.38064852e-23;

  for (int i = 0; i < nmode; i++) {
    const double vti = vibtemp[i];
    f += (((kb*vti)/(exp(vti/Tvib)-1)));
  }

  return f;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
double bird_dEvib(const int& nmode, const double& Tvib, const double vibtemp[]) const
{
  // Comutes df for Newton's search method

  double df = 0.0;
  const double kb = 1.38064852e-23;

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

KOKKOS_INLINE_FUNCTION
double newtonTvib(const int &nmode, const double& Evib, const double vibTemp[],
               const double &Tvib0,
               const double &tol,
               const int& nmax) const
{
  // Function for converting vibrational energy to vibrational temperature
  // Computes Tvib assuming the vibrational energy levels occupy a simple harmonic oscillator (SHO) spacing
  // Search for Tvib begins at some initial value "Tvib0" until the search reaches a tolerance level "tol"

  // Uses Newton's method to solve for a vibrational temperature given a
  // distribution of vibrational energy levels

  double Tvib_prev;

  // f and df are computed for Newton's search
  double f = bird_Evib(nmode,Tvib0,vibTemp,Evib);
  double df = bird_dEvib(nmode,Tvib0,vibTemp);

  // Update guess for Tvib and compute error
  double Tvib = Tvib0 - (f/df);
  double err = fabs(Tvib-Tvib0);

  int i = 2;

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
/* ---------------------------------------------------------------------- */

enum{NONE,DISCRETE,SMOOTH};
enum{DISSOCIATION,EXCHANGE,IONIZATION,RECOMBINATION};   // other files

KOKKOS_INLINE_FUNCTION
int attempt_kk(Particle::OnePart *ip, Particle::OnePart *jp,
         double pre_etrans, double pre_erot, double pre_evib,
         double &post_etotal, int &kspecies,
         int &recomb_species, double &recomb_density, const t_species_1d_const &d_species) const
{
  OneReactionKokkos *r;

  const int isp = ip->ispecies;
  const int jsp = jp->ispecies;
  const double ievib = ip->evib;
  const double jevib = jp->evib;

  const double pre_ave_rotdof = (d_species[isp].rotdof + d_species[jsp].rotdof)/2.0;

  const int n = d_reactions(isp,jsp).n;
  if (n == 0) return 0;
  auto& d_list = d_reactions(isp,jsp).d_list;

  // probablity to compare to reaction probability

  double react_prob = 0.0;
  rand_type rand_gen = rand_pool.get_state();
  const double random_prob = rand_gen.drand();
  rand_pool.free_state(rand_gen);
  double zi = 0.0;
  double zj = 0.0;
  int avei = 0;
  int avej = 0;
  double iTvib = 0.0;
  double jTvib = 0.0;

  // loop over possible reactions for these 2 species

  for (int i = 0; i < n; i++) {
    r = &d_rlist[d_list[i]];

    // ignore energetically impossible reactions

    const double pre_etotal = pre_etrans + pre_erot + pre_evib;

    // two options for total energy in TCE model
    // 0: partialEnergy = true: rDOF model
    // 1: partialEnergy = false: TCE: Rotation + Vibration

    // average DOFs participating in the reaction

    double ecc,z;
    double e_excess = 0.0;

    if (partialEnergy) {
      ecc = pre_etrans;
      z = r->d_coeff[0];
      if (pre_ave_rotdof > 0.1)
        ecc += pre_erot*z/pre_ave_rotdof;
    } else {
      ecc = pre_etotal;
      z = pre_ave_rotdof;
    }

    // Cover cases where coeff[1].neq.coeff[4]
    if (r->d_coeff[1]>((-1)*r->d_coeff[4])) e_excess = ecc - r->d_coeff[1];
    else e_excess = ecc + r->d_coeff[4];
    if (e_excess <= 0.0) continue;

    if (!partialEnergy) {

       if (vibstyle == SMOOTH) z += (d_species[isp].vibdof + d_species[jsp].vibdof)/2.0;
       else if (vibstyle == DISCRETE) {
            //Instantaneous z for diatomic molecules
            if (d_species[isp].nvibmode == 1) {
                avei = static_cast<int>
                        (ievib / (boltz * d_species[isp].vibtemp[0]));
                if (avei > 0) zi = 2.0 * avei * log(1.0 / avei + 1.0);
                else zi = 0.0;
            } else if (d_species[isp].nvibmode > 1) {
                if (ievib < 1e-26 ) zi = 0.0; //Low Energy Cut-Off to prevent nan solutions to newtonTvib
                //Instantaneous T for polyatomic
                else {
                  iTvib = newtonTvib(d_species[isp].nvibmode,ievib,d_species[isp].vibtemp,3000,1e-4,1000);
                  zi = (2 * ievib)/(boltz * iTvib);
                }
            } else zi = 0.0;

            if (d_species[jsp].nvibmode == 1) {
                avej = static_cast<int>
                        (jevib / (boltz * d_species[jsp].vibtemp[0]));
                if (avej > 0) zj = 2.0 * avej * log(1.0 / avej + 1.0);
                else zj = 0.0;
            } else if (d_species[jsp].nvibmode > 1) {
                if (jevib < 1e-26) zj = 0.0;
                else {
                  jTvib = newtonTvib(d_species[jsp].nvibmode,jevib,d_species[jsp].vibtemp,3000,1e-4,1000);
                  zj = (2 * jevib)/(boltz * jTvib);
                }
            } else zj = 0.0;

            if (isnan(zi) || isnan(zj) || zi < 0 || zj < 0) Kokkos::abort("Root-Finding Error\n");
            z += 0.5 * (zi+zj);
       }
    }

    // compute probability of reaction

    switch (r->type) {
    case DISSOCIATION:
    case IONIZATION:
    case EXCHANGE:
      {
        react_prob += r->d_coeff[2] * tgamma(z+2.5-r->d_coeff[5]) / MAX(1.0e-6,tgamma(z+r->d_coeff[3]+1.5)) *
          pow(ecc-r->d_coeff[1],r->d_coeff[3]-1+r->d_coeff[5]) *
          pow(1.0-r->d_coeff[1]/ecc,z+1.5-r->d_coeff[5]);
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
        auto& d_sp2recomb = d_reactions(isp,jsp).d_sp2recomb;
        if (d_sp2recomb[recomb_species] != d_list[i]) continue;

        react_prob += recomb_boost * recomb_density * r->d_coeff[2] *
          tgamma(z+2.5-r->d_coeff[5]) / MAX(1.0e-6,tgamma(z+r->d_coeff[3]+1.5)) *
          pow(ecc-r->d_coeff[1],r->d_coeff[3]-1+r->d_coeff[5]) *  // extended to general recombination case with non-zero activation energy
          pow(1.0-r->d_coeff[1]/ecc,z+1.5-r->d_coeff[5]);
        break;
      }

      //if (react_prob < 0) error->warning(FLERR,"Negative reaction probability");
      //else if (react_prob > 1) error->warning(FLERR,"Reaction probability greater than 1");

    default:
      //error->one(FLERR,"Unknown outcome in reaction");
      //d_error_flag() = 1;
      Kokkos::abort("ReactTCEKokkos: Unknown outcome in reaction\n");
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
      Kokkos::atomic_increment(&d_tally_reactions[d_list[i]]);
      if (!computeChemRates) {
        ip->ispecies = r->d_products[0];

        switch (r->type) {
        case DISSOCIATION:
        case IONIZATION:
        case EXCHANGE:
          {
            jp->ispecies = r->d_products[1];
            break;
          }
        case RECOMBINATION:
          {
            // always destroy 2nd reactant species

            jp->ispecies = -1;
            break;
          }
        }

        if (r->nproduct > 2) kspecies = r->d_products[2];
        else kspecies = -1;

        post_etotal = pre_etotal + r->d_coeff[4];

        return 1;
      } else {
        return 0;
      }
    }
  }

  return 0;
}

/* ---------------------------------------------------------------------- */

 protected:
  int vibstyle;
  double boltz;

  DAT::tdual_int_scalar k_error_flag;
  DAT::t_int_scalar d_error_flag;
  HAT::t_int_scalar h_error_flag;
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/

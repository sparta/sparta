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

#ifdef REACT_CLASS

ReactStyle(tce/kk,ReactTCEKokkos)

#else

#ifndef SPARTA_REACT_TCE_KOKKOS_H
#define SPARTA_REACT_TCE_KOKKOS_H

#include "react_bird_kokkos.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

class ReactTCEKokkos : public ReactBirdKokkos {
 public:
  ReactTCEKokkos(class SPARTA *, int, char **);
  ReactTCEKokkos(class SPARTA* sparta) : ReactBirdKokkos(sparta) {};
  void init();
  int attempt(Particle::OnePart *, Particle::OnePart *,
              double, double, double, double &, int &) { return 0; }

/* ---------------------------------------------------------------------- */

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

  const double pre_ave_rotdof = (d_species[isp].rotdof + d_species[jsp].rotdof)/2.0;

  const int n = d_reactions(isp,jsp).n;
  if (n == 0) return 0;
  auto& d_list = d_reactions(isp,jsp).d_list;

  // probablity to compare to reaction probability

  double react_prob = 0.0;
  rand_type rand_gen = rand_pool.get_state();
  const double random_prob = rand_gen.drand();
  rand_pool.free_state(rand_gen);

  // loop over possible reactions for these 2 species

  for (int i = 0; i < n; i++) {
    r = &d_rlist[d_list[i]];

    // ignore energetically impossible reactions

    const double pre_etotal = pre_etrans + pre_erot + pre_evib;

    double ecc = pre_etrans;
    if (pre_ave_rotdof > 0.1) ecc += pre_erot*r->d_coeff[0]/pre_ave_rotdof;

    const double e_excess = ecc - r->d_coeff[1];
    if (e_excess <= 0.0) continue;

    // compute probability of reaction

    switch (r->type) {
    case DISSOCIATION:
    case IONIZATION:
    case EXCHANGE:
      {
        react_prob += r->d_coeff[2] *
          pow(ecc-r->d_coeff[1],r->d_coeff[3]) *
          pow(1.0-r->d_coeff[1]/ecc,r->d_coeff[5]);
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
          pow(ecc,r->d_coeff[3]) *
          pow(1.0-r->d_coeff[1]/ecc,r->d_coeff[5]);
        break;
      }

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
      ip->ispecies = r->d_products[0];

      // Previous statment did not destroy the 2nd species (B) if
      //   recombination was specified as A+B->AB+M (which has nproductus=2)
      //   but only for the A+B->AB specication form (which has nproductus=1)

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
    }
  }

  return 0;
}

/* ---------------------------------------------------------------------- */

 protected:
  DAT::tdual_int_scalar k_error_flag;
  DAT::t_int_scalar d_error_flag;
  HAT::t_int_scalar h_error_flag;
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/

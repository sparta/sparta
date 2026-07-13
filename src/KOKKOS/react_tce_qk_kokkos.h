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

ReactStyle(tce/qk/kk,ReactTCEQKKokkos)

#else

#ifndef SPARTA_REACT_TCE_QK_KOKKOS_H
#define SPARTA_REACT_TCE_QK_KOKKOS_H

#include "math.h"
#include "react_bird_kokkos.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

class ReactTCEQKKokkos : public ReactBirdKokkos {
 public:
  ReactTCEQKKokkos(class SPARTA *, int, char **);
  ReactTCEQKKokkos(class SPARTA* sparta) : ReactBirdKokkos(sparta) {copy = 1;}
  void init();
  int attempt(Particle::OnePart *, Particle::OnePart *,
              double, double, double, double &, int &) {return 0;}

  enum{DISSOCIATION,EXCHANGE,IONIZATION,RECOMBINATION};   // reaction types
  enum{ARRHENIUS,QUANTUM};                                // reaction styles

/* ----------------------------------------------------------------------
   hybrid TCE/QK reaction attempt, device version
   mirrors ReactTCEQK::attempt(): per reaction, ARRHENIUS style uses the simple
   TCE probability, QUANTUM style uses the QK model; each evaluated reaction
   draws its own random number (matching the host RNG order)
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
int attempt_kk(Particle::OnePart *ip, Particle::OnePart *jp,
         double pre_etrans, double pre_erot, double pre_evib,
         double &post_etotal, int &kspecies,
         int & /*recomb_species*/, double & /*recomb_density*/,
         const t_species_1d_const &d_species) const
{
  const int isp = ip->ispecies;
  const int jsp = jp->ispecies;

  const int n = d_reactions(isp,jsp).n;
  if (n == 0) return 0;
  auto& d_list = d_reactions(isp,jsp).d_list;

  const double pre_ave_rotdof = (d_species[isp].rotdof + d_species[jsp].rotdof)/2.0;
  const double omega = d_omega(isp,jsp);

  rand_type rand_gen = rand_pool.get_state();

  for (int i = 0; i < n; i++) {
    OneReactionKokkos *r = &d_rlist[d_list[i]];

    const double pre_etotal = pre_etrans + pre_erot + pre_evib;

    // top-level energetic-possibility screen (uses total energy)

    double ecc = pre_etotal;
    if (ecc - r->d_coeff[1] <= 0.0) continue;

    int fired = 0;

    // per-reaction probability + its own random draw (helper semantics)

    const double random_prob = rand_gen.drand();

    double react_prob = 0.0;
    double ecc2 = pre_etrans;
    if (pre_ave_rotdof > 0.1) ecc2 += pre_erot*r->d_coeff[0]/pre_ave_rotdof;
    const double e_excess = ecc2 - r->d_coeff[1];

    if (e_excess > 0.0) {
      if (r->style == ARRHENIUS) {                 // attempt_tce
        switch (r->type) {
        case DISSOCIATION:
        case EXCHANGE:
          react_prob += r->d_coeff[2] *
            pow(ecc2-r->d_coeff[1],r->d_coeff[3]) *
            pow(1.0-r->d_coeff[1]/ecc2,r->d_coeff[5]);
          break;
        default:
          Kokkos::abort("ReactTCEQKKokkos: Unknown outcome in reaction\n");
          break;
        }
        if (react_prob > random_prob) fired = 1;

      } else {                                     // attempt_qk
        const double inverse_kT = 1.0 / (boltz * d_species[isp].vibtemp[0]);
        int iv = 0,ilevel,maxlev,limlev;
        double eccq;
        switch (r->type) {
        case DISSOCIATION:
          {
            eccq = pre_etrans + ip->evib;
            maxlev = static_cast<int> (eccq * inverse_kT);
            limlev = static_cast<int> (fabs(r->d_coeff[1]) * inverse_kT);
            if (maxlev > limlev) react_prob = 1.0;
            break;
          }
        case EXCHANGE:
          {
            if (r->d_coeff[4] < 0.0 && d_species[isp].rotdof > 0) {
              eccq = pre_etrans + ip->evib;
              maxlev = static_cast<int> (eccq * inverse_kT);
              if (eccq > r->d_coeff[1]) {
                do {
                  iv = static_cast<int> (rand_gen.drand()*(maxlev+0.99999999));
                  double evib = static_cast<double> (iv / inverse_kT);
                  if (evib < eccq) react_prob = pow(1.0-evib/eccq,1.5-omega);
                } while (rand_gen.drand() < react_prob);
                ilevel = static_cast<int> (fabs(r->d_coeff[4]) * inverse_kT);
                if (iv >= ilevel) react_prob = 1.0;
              }
            } else if (r->d_coeff[4] > 0.0 && d_species[isp].rotdof > 0) {
              eccq = pre_etrans + ip->evib;
              int mspec = r->d_products[0];
              if (d_species[mspec].rotdof < 2.0) mspec = r->d_products[1];
              eccq += r->d_coeff[4];
              maxlev = static_cast<int> (eccq * inverse_kT);
              double prob = 0.0;
              do {
                iv = rand_gen.drand()*(maxlev+0.99999999);
                double evib = static_cast<double> (iv * boltz*d_species[mspec].vibtemp[0]);
                if (evib < eccq) prob = pow(1.0-evib/eccq,1.5 - r->d_coeff[6]);
              } while (rand_gen.drand() < prob);
              ilevel = static_cast<int> (fabs(r->d_coeff[4]/boltz/d_species[mspec].vibtemp[0]));
              if (iv >= ilevel) react_prob = 1.0;
            }
            break;
          }
        default:
          Kokkos::abort("ReactTCEQKKokkos: Unknown outcome in reaction\n");
          break;
        }
        if (react_prob > random_prob) fired = 1;
      }
    }

    if (fired) {
      Kokkos::atomic_inc(&d_tally_reactions[d_list[i]]);
      ip->ispecies = r->d_products[0];
      jp->ispecies = r->d_products[1];
      post_etotal = pre_etotal + r->d_coeff[4];
      if (r->nproduct > 2) kspecies = r->d_products[2];
      else kspecies = -1;
      rand_pool.free_state(rand_gen);
      return d_list[i] + 1;
    }
  }

  rand_pool.free_state(rand_gen);
  return 0;
}

 protected:
  double boltz;
  DAT::t_float_2d d_omega;     // VSS omega for each species pair
};

}

#endif
#endif

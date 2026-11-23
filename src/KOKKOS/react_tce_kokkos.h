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
  ReactTCEKokkos(class SPARTA* sparta) : ReactBirdKokkos(sparta) {copy = 1;}
  void init();
  void end_of_run();
  int attempt(Particle::OnePart *, Particle::OnePart *,
              double, double, double, double, double &, int &) { return 0; }

/* ---------------------------------------------------------------------- */

enum{NONE,DISCRETE,SMOOTH};
enum{DISSOCIATION,EXCHANGE,IONIZATION,RECOMBINATION};   // other files
enum{PROBBIG_NONE = 1073741824};                  // no d_probbig entry

KOKKOS_INLINE_FUNCTION
int attempt_kk(Particle::OnePart *ip, Particle::OnePart *jp,
         double pre_etrans, double pre_erot, double pre_evib, double pre_eelec,
         double &post_etotal, int &kspecies,
         int &recomb_species, double &recomb_density,
         Particle::OnePart *recomb_part3, const double recomb_p3_eelec,
         const t_species_1d_const &d_species,
         const double tgas_cell) const
{
  OneReactionKokkos *r;

  const int isp = ip->ispecies;
  const int jsp = jp->ispecies;

  const double pre_ave_rotdof = (d_species[isp].rotdof + d_species[jsp].rotdof)/2.0;

  const int n = d_reactions(isp,jsp).n;
  if (n == 0) return 0;
  auto& d_list = d_reactions(isp,jsp).d_list;

  rand_type rand_gen = rand_pool.get_state();
  const double random_prob = rand_gen.drand();
  rand_pool.free_state(rand_gen);

  // pass 1: cumulative channel selection, each channel's probability
  // individually capped at 1 (sigma_R <= sigma_VHS), matching
  // ReactTCE::attempt exactly; no early exit so the channel split is not
  // biased by the ordering of reactions in the input file

  double total_prob = 0.0;
  int isel = -1;

  for (int i = 0; i < n; i++) {
    double p = channel_prob_kk(d_list[i],ip,jp,pre_etrans,pre_erot,
                               pre_evib,pre_eelec,pre_ave_rotdof,
                               recomb_species,recomb_density,recomb_part3,
                               recomb_p3_eelec,d_species,tgas_cell);
    if (p < 0.0) Kokkos::atomic_max(&d_probneg(),1);
    int over = 0;
    if (p > 1.0) { p = 1.0; over = 1; }
    const double before = total_prob;
    total_prob += p;
    // a channel whose own probability was capped, or the channel at which
    // the pair's summed probability crossed 1 (matches ReactTCE::attempt)
    if (over || (before <= 1.0 && total_prob > 1.0))
      Kokkos::atomic_min(&d_probbig(),d_list[i]);
    if (isel < 0 && total_prob > random_prob) isel = i;
  }

  // saturation: the pair reacts with probability 1 and the channel is
  // selected proportionally (P_i/total) by re-walking the deterministic
  // cumulative sum against u*total; this slow second pass runs only for
  // saturated pairs (matches ReactTCE::attempt)

  if (total_prob > 1.0) {
    const double target = random_prob*total_prob;
    double cum = 0.0;
    isel = -1;
    for (int i = 0; i < n; i++) {
      double p = channel_prob_kk(d_list[i],ip,jp,pre_etrans,pre_erot,
                                 pre_evib,pre_eelec,pre_ave_rotdof,
                                 recomb_species,recomb_density,recomb_part3,
                                 recomb_p3_eelec,d_species,tgas_cell);
      if (p > 1.0) p = 1.0;
      cum += p;
      if (p > 0.0) isel = i;               // round-off fallback: last live channel
      if (cum > target) break;
    }
  }

  if (isel < 0) return 0;

  // selected channel: tally it; in compute_chem_rates mode the reaction is
  // tallied but not performed.  Exactly one channel is tallied per
  // reacting collision (matches ReactTCE::attempt)

  r = &d_rlist[d_list[isel]];
  Kokkos::atomic_inc(&d_tally_reactions[d_list[isel]]);

  if (computeChemRates) return 0;

  // perform the reaction: reset species of I,J and optional K to products.
  // J particle is destroyed in a recombination reaction (species = -1); a
  // K particle can be created in a dissociation or ionization reaction
  // (parent creates it from kspecies)

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

  post_etotal = pre_etrans + pre_erot + pre_evib + pre_eelec + r->d_coeff[4];

  // return reaction from 1 to N

  return d_list[isel] + 1;
}

/* ----------------------------------------------------------------------
   TCE probability of reaction rindex for the collision of ip,jp; device
   twin of ReactTCE::channel_prob, kept line-for-line consistent for
   SPARTA_KOKKOS_EXACT parity.  Deterministic: no random numbers drawn.
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
double channel_prob_kk(int rindex, Particle::OnePart *ip,
         Particle::OnePart *jp,
         double pre_etrans, double pre_erot, double pre_evib,
         double pre_eelec, double pre_ave_rotdof,
         int &recomb_species, double &recomb_density,
         Particle::OnePart *recomb_part3, const double recomb_p3_eelec,
         const t_species_1d_const &d_species,
         const double tgas_cell) const
{
  const int isp = ip->ispecies;
  const int jsp = jp->ispecies;

  OneReactionKokkos *r = &d_rlist[rindex];

    // an inactive channel contributes nothing (matches ReactTCE::channel_prob)

    if (!r->active) return 0.0;

    // ignore energetically impossible reactions

    const double pre_etotal = pre_etrans + pre_erot + pre_evib + pre_eelec;

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
      // total energy model, matching ReactTCE::attempt: ecc is the TOTAL
      // collision energy including electronic; z is the continuum
      // internal DOF (rotation, plus vibration when it is a continuous
      // mode); the discrete ladders enter through the tables below
      ecc = pre_etotal;
      z = pre_ave_rotdof;
      if (vibstyle == SMOOTH)
        z += (d_species[isp].vibdof + d_species[jsp].vibdof)/2.0;
    }

    // a reverse (B-style) recombination is fully microcanonical, matching
    // ReactTCE::attempt: the 3rd particle's energy counts toward the
    // barrier, so the pair-energy threshold below does not apply to it; an
    // external-Keq recombination uses the same 3-body table times the
    // residual thermal factor applied below

    int micro3 = 0;
    if (r->reverse && r->type == RECOMBINATION &&
        d_mtab_n[rindex] > 0 && d_mtab_num_flag[rindex]) micro3 = 1;

    // Cover cases where coeff[1].neq.coeff[4]

    if (r->d_coeff[1]>((-1)*r->d_coeff[4])) e_excess = ecc - r->d_coeff[1];
    else e_excess = ecc + r->d_coeff[4];
    if (e_excess <= 0.0 && !micro3) return 0.0;

    // energy-dependent factor of the TCE probability: the tabulated
    // microcanonical average over the reactants' discrete vibrational
    // and electronic ladders, matching ReactBird::vib_micro_factor
    // exactly; reactions without discrete ladders use the standard
    // analytic factor

    double efactor = 0.0;
    if (!micro3) {   // a micro3 reaction uses its own 3-body factor below
      if (!partialEnergy && d_mtab_n[rindex] > 0) {
        const int ir = rindex;
        const int n = d_mtab_n[ir];
        const double x = ecc/d_mtab_du[ir];
        const int k = (int) x;
        if (k >= n-1) efactor = d_mtab(ir,n-1);
        else {
          const double f = x - k;
          efactor = (1.0-f)*d_mtab(ir,k) + f*d_mtab(ir,k+1);
        }
      } else
        efactor = pow(ecc-r->d_coeff[1],r->d_coeff[3]-1+r->d_coeff[5]) *
                  pow(1.0-r->d_coeff[1]/ecc,z+1.5-r->d_coeff[5]);
    }

    // residual thermal correction for an external-Keq reverse reaction,
    // matching ReactTCE::attempt: multiply the detailed-balance table by
    // R(T) = Keq_statmech/Keq_ext at the cell temperature so the thermal
    // average matches the external Keq while the shape stays energy-
    // resolved; other reverse reactions need no temperature (keq_resid = 1)

    double keq_resid = 1.0;
    if (r->reverse && r->keq_flag) {
      if (tgas_cell > 0.0) {
        // clamp to the residual fit window (1000-60000 K), matching
        // ReactTCE::attempt, to avoid an uncontrolled Park extrapolation
        const double tr = tgas_cell < 1000.0 ? 1000.0 :
          (tgas_cell > 60000.0 ? 60000.0 : tgas_cell);
        const double z10 = 10000.0/tr;
        keq_resid = exp(r->keq_resid_coeff[0]/z10 + r->keq_resid_coeff[1] +
                        r->keq_resid_coeff[2]*log(z10) +
                        r->keq_resid_coeff[3]*z10 +
                        r->keq_resid_coeff[4]*z10*z10);
      } else return 0.0;
    }
    double prefactor = r->d_coeff[2] * keq_resid;

    // compute probability of reaction
    // gamma function denominator is non-positive if the temperature
    //   exponent is below the TCE bound, then it is clamped to 1.0e-6
    //   and the reaction rate is incorrect, warned about at init by
    //   ReactBird::check_tce_bounds()

    switch (r->type) {
    case DISSOCIATION:
    case IONIZATION:
    case EXCHANGE:
      {
        return prefactor * tgamma(z+2.5-r->d_coeff[5]) / MAX(1.0e-6,tgamma(z+r->d_coeff[3]+1.5)) *
          efactor;
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

        if (recomb_species < 0) return 0.0;
        auto& d_sp2recomb = d_reactions(isp,jsp).d_sp2recomb;
        if (d_sp2recomb[recomb_species] != rindex) return 0.0;

        // fully microcanonical reverse recombination, matching
        // ReactTCE::attempt exactly (see ReactBird::build_db3_table);
        // no cell temperature is used

        if (micro3) {
          Particle::OnePart *p3 = recomb_part3;
          const int sp3 = p3->ispecies;

          const double mi = d_species[isp].mass;
          const double mj = d_species[jsp].mass;
          const double m3 = d_species[sp3].mass;
          const double divisor = 1.0/(mi + mj);
          double *vi = ip->v;
          double *vj = jp->v;
          double *v3 = p3->v;
          const double du3 = v3[0] - (mi*vi[0] + mj*vj[0])*divisor;
          const double dv3 = v3[1] - (mi*vi[1] + mj*vj[1])*divisor;
          const double dw3 = v3[2] - (mi*vi[2] + mj*vj[2])*divisor;
          const double mu3 = m3*(mi + mj)/(m3 + mi + mj);
          const double eps_t = 0.5*mu3*(du3*du3 + dv3*dv3 + dw3*dw3);

          const double eelec3 = recomb_p3_eelec;

          // continuum density weights of the 3rd particle's energies,
          // matching the flat-measure dimension of the table; skip the
          // attempt at the (measure-zero) singular points

          double c3 = sqrt(eps_t);
          const int rotdof3 = d_species[sp3].rotdof;
          if (rotdof3 > 0 && rotdof3 != 2)
            c3 *= pow(p3->erot,0.5*rotdof3-1.0);
          const int vibdof3 = d_species[sp3].vibdof;
          if (vibstyle == SMOOTH && vibdof3 > 0 && vibdof3 != 2)
            c3 *= pow(p3->evib,0.5*vibdof3-1.0);
          if (!(c3 > 0.0)) return 0.0;

          const int ir = rindex;
          const int ntab = d_mtab_n[ir];
          const double dutab = d_mtab_du[ir];

          double xpair;
          {
            const double x = ecc/dutab;
            const int k = (int) x;
            if (k >= ntab-1) xpair = d_mtab(ir,ntab-1);
            else {
              const double f = x - k;
              xpair = (1.0-f)*d_mtab(ir,k) + f*d_mtab(ir,k+1);
            }
          }
          if (!(xpair > 0.0)) return 0.0;

          const double w = ecc + eps_t + p3->erot + p3->evib + eelec3;
          double numw;
          {
            const double x = w/dutab;
            const int k = (int) x;
            if (k >= ntab-1) numw = d_mtab_num(ir,ntab-1);
            else {
              const double f = x - k;
              numw = (1.0-f)*d_mtab_num(ir,k) + f*d_mtab_num(ir,k+1);
            }
          }

          return keq_resid * recomb_boost * recomb_density * numw / (xpair*c3);
        }

        return recomb_boost * recomb_density * prefactor *
          tgamma(z+2.5-r->d_coeff[5]) / MAX(1.0e-6,tgamma(z+r->d_coeff[3]+1.5)) *
          efactor;   // extended to general recombination case with non-zero activation energy
      }

      //if (react_prob < 0) error->warning(FLERR,"Negative reaction probability");
      //else if (react_prob > 1) error->warning(FLERR,"Reaction probability greater than 1");

    default:
      //error->one(FLERR,"Unknown outcome in reaction");
      //d_error_flag() = 1;
      Kokkos::abort("ReactTCEKokkos: Unknown outcome in reaction\n");
      break;
    }

  return 0.0;
}

/* ---------------------------------------------------------------------- */

 protected:
  int vibstyle;
  int elecstyle;
  double boltz;

  // per-reaction microcanonical energy-factor tables, device mirrors of
  // the ReactBird host tables; n = 0 marks "no table"; d_mtab_num rows
  // are only nonzero for reverse recombinations (3-body detailed
  // balance, see ReactBird::build_db3_table)

  // kept in double whatever SPA_PRECISION is: the entries span many
  // decades (a recombination row holds N_F/V_3 and the pair density of
  // states), and in single precision the small ones underflow to zero

  typedef Kokkos::View<double**, DeviceType::array_layout, DeviceType> t_table_2d;
  typedef Kokkos::View<double*, DeviceType::array_layout, DeviceType> t_table_1d;
  t_table_2d d_mtab;
  t_table_2d d_mtab_num;
  t_table_1d d_mtab_du;
  DAT::t_int_1d d_mtab_n;
  DAT::t_int_1d d_mtab_num_flag;   // 1 if mtab_num row is a real 3-body
                                   //   detailed-balance table (mirrors the
                                   //   host mtab_num[i] != NULL test)

  DAT::tdual_int_scalar k_error_flag;
  DAT::t_int_scalar d_error_flag;
  HAT::t_int_scalar h_error_flag;

  // invalid reaction probabilities recorded during the run, reported by
  //   end_of_run() (match ReactTCE::prob_neg_flag and prob_big_index):
  //   d_probneg = 1 if a negative probability occurred,
  //   d_probbig = lowest rlist index at which a pair's summed probability
  //     exceeded 1, PROBBIG_NONE if none
  DAT::t_int_scalar d_probneg,d_probbig;
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/

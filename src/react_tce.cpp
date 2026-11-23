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
#include <string>
#include "react_tce.h"
#include "particle.h"
#include "collide.h"
#include "update.h"
#include "random_knuth.h"
#include "math_const.h"
#include "comm.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{NONE,DISCRETE,SMOOTH};
enum{DISSOCIATION,EXCHANGE,IONIZATION,RECOMBINATION};   // other files

#define MAXCACHE 32     // channels per species pair whose first-pass
                        //   probabilities attempt() keeps

/* ---------------------------------------------------------------------- */

ReactTCE::ReactTCE(SPARTA *sparta, int narg, char **arg) :
  ReactBird(sparta, narg, arg)
{
  prob_neg_flag = 0;
  prob_big_index = -1;
  tce_style = 1;
}

/* ---------------------------------------------------------------------- */

void ReactTCE::init()
{
  if (!collide || strcmp(collide->style,"vss") != 0)
    error->all(FLERR,"React tce can only be used with collide vss");

  prob_neg_flag = 0;
  prob_big_index = -1;

  ReactBird::init();

  // custom electronic energy of the 3rd particle, needed by the 3-body
  // detailed-balance probability of a reverse recombination

  index_eelec = -1;
  if (collide->elecstyle == DISCRETE)
    index_eelec = particle->find_custom((char *) "eelec");

  // with partial_energy no, the TCE energy factor is the microcanonical
  // average over the reactants' discrete ladders (SHO vibrational levels
  // when vibrate discrete, electronic states when electronic discrete):
  // build the per-reaction lookup tables; reactions whose reactants carry
  // no discrete ladders fall back to the standard analytic factor, which
  // is the same average with an empty ladder

  if (!partialEnergy) build_micro_tables();
  else free_micro_tables();

  // warn if temperature exponent of any reaction is out of bounds
  //   for the TCE reaction probability

  check_tce_bounds();
}

/* ----------------------------------------------------------------------
   warn once, over all procs, if any invalid TCE reaction probability
     occurred during the run
   called by Finish on all procs at the end of each run
------------------------------------------------------------------------- */

void ReactTCE::end_of_run()
{
  report_tce_prob_flags(prob_neg_flag,prob_big_index);
  prob_neg_flag = 0;
  prob_big_index = -1;
}

/* ---------------------------------------------------------------------- */

int ReactTCE::attempt(Particle::OnePart *ip, Particle::OnePart *jp,
                      double pre_etrans, double pre_erot, double pre_evib, double pre_eelec,
                      double &post_etotal, int &kspecies)
{
  OneReaction *r;

  Particle::Species *species = particle->species;
  int isp = ip->ispecies;
  int jsp = jp->ispecies;

  double pre_ave_rotdof = (species[isp].rotdof + species[jsp].rotdof)/2.0;

  int n = reactions[isp][jsp].n;
  if (n == 0) return 0;
  int *list = reactions[isp][jsp].list;

  double random_prob = random->uniform();

  // pass 1: cumulative channel selection, each channel's probability
  // individually capped at 1 (sigma_R <= sigma_VHS: the TCE derivation
  // assumes the reaction cross-section is a fraction of the VHS
  // cross-section, so a computed probability above 1 is outside the model
  // and is limited to the collision rate; Higdon 2018, following Strand &
  // Goldstein).  No early exit on selection: the full channel sum is
  // needed to detect saturation of the PAIR, so that the channel split is
  // not biased by the ordering of reactions in the input file.  The extra
  // cost is only the remainder of this pair's (short) reaction list on
  // the rare collisions that select a reaction; nothing is computed for
  // collisions of pairs with no reactions, and the collision-acceptance
  // machinery ((sigma_T c_r)_max) is untouched.

  double total_prob = 0.0;
  int isel = -1;

  // the capped channel probabilities are kept for the saturation pass
  // below (channel_prob() is deterministic, so this is exact); a pair with
  // more than MAXCACHE channels recomputes the rest there

  double pcache[MAXCACHE];

  for (int i = 0; i < n; i++) {
    double p = channel_prob(list[i],ip,jp,pre_etrans,pre_erot,
                            pre_evib,pre_eelec,pre_ave_rotdof);
    if (p < 0.0) prob_neg_flag = 1;
    int over = 0;
    if (p > 1.0) { p = 1.0; over = 1; }
    if (i < MAXCACHE) pcache[i] = p;
    double before = total_prob;
    total_prob += p;
    // record a channel whose own probability was capped, or the channel
    // at which the pair's summed probability crossed 1 (end_of_run()
    // names the lowest such reaction index)
    if ((over || (before <= 1.0 && total_prob > 1.0)) &&
        (prob_big_index < 0 || list[i] < prob_big_index))
      prob_big_index = list[i];
    if (isel < 0 && total_prob > random_prob) isel = i;
  }

  // saturation: the pair's total reaction probability exceeds 1, so the
  // collision reacts with probability 1 and the channel must be selected
  // proportionally (P_i/total) to remove the ordering bias of the
  // sequential test above.  Re-walk the (deterministic) cumulative sum
  // against u*total, from the probabilities of the first pass.

  if (total_prob > 1.0) {
    double target = random_prob*total_prob;
    double cum = 0.0;
    isel = -1;
    for (int i = 0; i < n; i++) {
      double p = (i < MAXCACHE) ? pcache[i] :
        channel_prob(list[i],ip,jp,pre_etrans,pre_erot,
                     pre_evib,pre_eelec,pre_ave_rotdof);
      if (p > 1.0) p = 1.0;
      cum += p;
      if (p > 0.0) isel = i;               // round-off fallback: last live channel
      if (cum > target) break;
    }
  }

  // the probability is per collision, so it does not depend on timestep
  //   or fnum; only record it here, end_of_run() warns once over all procs

  if (isel < 0) return 0;

  // selected channel: tally it; in compute_chem_rates mode the reaction is
  // tallied but not performed.  Exactly one channel is tallied per reacting
  // collision, so each channel's tally rate is P_i (or P_i/total when the
  // pair is saturated), independent of its position in the reaction file.

  r = &rlist[list[isel]];
  tally_reactions[list[isel]]++;

  if (computeChemRates) return 0;

  // perform the reaction: reset species of I,J and optional K to products.
  // J particle is destroyed in a recombination reaction (species = -1); a
  // K particle can be created in a dissociation or ionization reaction
  // (parent creates it from kspecies).
  // important NOTE:
  //   it does not matter what order the I,J reactants are in compared to
  //   the order listed in the reaction file: the list of possible
  //   reactions includes all reactions the I,J species are in, and the
  //   pre-reaction state (precoln) stores only combined I,J properties

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

  post_etotal = pre_etrans + pre_erot + pre_evib + pre_eelec + r->coeff[4];

  // return reaction from 1 to N

  return list[isel] + 1;
}

/* ----------------------------------------------------------------------
   TCE probability of reaction rindex (index into rlist) for the collision
   of ip,jp with the given pre-collision energies; returns 0 for channels
   that are skipped (energetically impossible, no matching 3rd body, no
   cell temperature for an external-Keq reverse).  Deterministic: no
   random numbers are drawn, so the saturation path of attempt() can
   re-walk the cumulative sum exactly.
------------------------------------------------------------------------- */

double ReactTCE::channel_prob(int rindex, Particle::OnePart *ip,
                              Particle::OnePart *jp,
                              double pre_etrans, double pre_erot,
                              double pre_evib, double pre_eelec,
                              double pre_ave_rotdof)
{
  double ecc,e_excess,z;

  Particle::Species *species = particle->species;
  int isp = ip->ispecies;
  int jsp = jp->ispecies;

  OneReaction *r = &rlist[rindex];

    // an inactive channel (e.g. a reverse reaction that was switched off)
    // must contribute nothing.  ReactBird::init() builds the per-pair
    // lists after all deactivations, so this is only a safeguard (matched
    // in ReactTCEKokkos::channel_prob_kk)

    if (!r->active) return 0.0;

    // ignore energetically impossible reactions

    double pre_etotal = pre_etrans + pre_erot + pre_evib + pre_eelec;

    // two options for total energy in TCE model
    // 0: partialEnergy = true: rDOF model
    // 1: partialEnergy = false: TCE: Rotation + Vibration

    // average DOFs participating in the reaction

    if (partialEnergy) {
       ecc = pre_etrans;
       z = r->coeff[0];
       if (pre_ave_rotdof > 0.1) ecc += pre_erot*z/pre_ave_rotdof;
    } else {
       // total energy model: ecc is the TOTAL collision energy including
       // the pair's electronic energy; z is the continuum internal DOF
       // (rotation, plus vibration when it is a continuous mode); the
       // discrete vibrational and electronic ladders enter the reaction
       // probability through the microcanonical energy-factor tables below
       ecc = pre_etotal;
       z = pre_ave_rotdof;
       if (collide->vibstyle == SMOOTH)
         z += (species[isp].vibdof + species[jsp].vibdof)/2.0;
    }

    // a reverse (B-style) recombination is fully microcanonical: the 3rd
    // particle's energy counts toward the barrier (its probability is
    // resolved in the total available energy), so the pair-energy threshold
    // below does not apply to it.  This holds for an external-Keq
    // recombination too: it uses the same 3-body table times the residual
    // thermal factor R(T) applied below.

    int micro3 = 0;
    if (r->reverse && r->type == RECOMBINATION &&
        mtab && mtab[rindex] && mtab_num[rindex]) micro3 = 1;

    // Cover cases where coeff[1].neq.coeff[4]

    e_excess = tce_excess(r,ecc);
    if (e_excess <= 0.0 && !micro3) return 0.0;

    // energy-dependent factor of the TCE probability:
    // the microcanonical (density-of-states weighted) average of the
    // standard factor over the reactants' discrete vibrational and
    // electronic ladders at fixed total collision energy, tabulated per
    // reaction at init; reactions whose reactants carry no discrete
    // ladders use the standard analytic factor
    //   (ecc-Ea)^(eta-1+omega) * (1-Ea/ecc)^(z+1.5-omega)
    // which is the same average with an empty ladder. The average keeps
    // the equilibrium rate on the input Arrhenius rate while the actual
    // internal states of the colliding pair count toward the barrier.

    double efactor = 0.0;
    if (!micro3)     // a micro3 reaction uses its own 3-body factor below
      efactor = tce_efactor(rindex,ecc,z,!partialEnergy);

    // residual thermal correction for a reverse reaction matched to an
    // external Keq curve fit (react_modify keq_file): its detailed-balance
    // table already reproduces the reverse rate for the statistical-
    // mechanics Keq, so multiplying by R(T) = Keq_statmech/Keq_ext at the
    // local cell temperature React::tgas makes the thermal average match
    // the EXTERNAL Keq while the energy-resolved shape stays microscopically
    // reversible.  If no cell temperature is available (fewer than 2
    // particles in the cell), skip the reaction for this collision.  Other
    // reverse reactions need no temperature (keq_resid = 1).

    double keq_resid = 1.0;
    if (r->reverse && r->keq_flag) {
      if (tgas > 0.0) {
        // clamp to the residual fit window (1000-60000 K, see
        // ReactBird::fit_keq_residual) so a very cold or very hot cell
        // gets the edge correction rather than an uncontrolled Park
        // extrapolation
        double tr = tgas < 1000.0 ? 1000.0 : (tgas > 60000.0 ? 60000.0 : tgas);
        keq_resid = keq_eval(r->keq_resid_coeff,tr);
      } else return 0.0;
    }
    double prefactor = r->coeff[2] * keq_resid;

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
        return tce_prefactor(r,z,prefactor) * efactor;
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
        int *sp2recomb = reactions[isp][jsp].sp2recomb;
        if (sp2recomb[recomb_species] != rindex) return 0.0;

        // fully microcanonical reverse recombination (issue #472): the
        // probability is resolved in the total available energy
        //   w = ecc + eps_t + erot3 + evib3 + eelec3
        // where eps_t is the 3rd particle's translational energy
        // relative to the pair's center of mass; dividing the
        // detailed-balance numerator table by the pair's density of
        // states and the 3rd particle's continuum density weights makes
        // the thermal average reproduce k_f(T)/K_eq(T) at every
        // temperature (see ReactBird::build_db3_table); no cell
        // temperature is used

        if (micro3) {
          Particle::OnePart *p3 = recomb_part3;
          int sp3 = p3->ispecies;

          double mi = species[isp].mass;
          double mj = species[jsp].mass;
          double m3 = species[sp3].mass;
          double divisor = 1.0/(mi + mj);
          double *vi = ip->v;
          double *vj = jp->v;
          double *v3 = p3->v;
          double du3 = v3[0] - (mi*vi[0] + mj*vj[0])*divisor;
          double dv3 = v3[1] - (mi*vi[1] + mj*vj[1])*divisor;
          double dw3 = v3[2] - (mi*vi[2] + mj*vj[2])*divisor;
          double mu3 = m3*(mi + mj)/(m3 + mi + mj);
          double eps_t = 0.5*mu3*(du3*du3 + dv3*dv3 + dw3*dw3);

          double eelec3 = 0.0;
          if (index_eelec >= 0 && species[sp3].elecdat)
            eelec3 = particle->edvec[particle->ewhich[index_eelec]]
              [p3 - particle->particles];

          // continuum density weights of the 3rd particle's energies,
          // matching the flat-measure dimension of the table; skip the
          // attempt at the (measure-zero) singular points

          double c3 = sqrt(eps_t);
          int rotdof3 = species[sp3].rotdof;
          if (rotdof3 > 0 && rotdof3 != 2)
            c3 *= pow(p3->erot,0.5*rotdof3-1.0);
          int vibdof3 = species[sp3].vibdof;
          if (collide->vibstyle == SMOOTH && vibdof3 > 0 && vibdof3 != 2)
            c3 *= pow(p3->evib,0.5*vibdof3-1.0);
          if (!(c3 > 0.0)) return 0.0;

          double xpair = vib_micro_factor(rindex,ecc);
          if (!(xpair > 0.0)) return 0.0;

          double w = ecc + eps_t + p3->erot + p3->evib + eelec3;
          return keq_resid * recomb_boost * recomb_density *
            db3_num_factor(rindex,w) / (xpair*c3);
        }

        return recomb_boost * recomb_density * prefactor *
          tgamma(z+2.5-r->coeff[5]) / MAX(1.0e-6,tgamma(z+r->coeff[3]+1.5)) *
          efactor;   // extended to general recombination case with non-zero activation energy
      }

    default:
      error->one(FLERR,"Unknown outcome in reaction");
      break;
    }

  return 0.0;
}

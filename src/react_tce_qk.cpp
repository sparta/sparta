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
#include "react_tce_qk.h"
#include "update.h"
#include "collide.h"
#include "random_knuth.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{DISSOCIATION,EXCHANGE,IONIZATION,RECOMBINATION};   // other files
enum{ARRHENIUS,QUANTUM};                                // other files

/* ---------------------------------------------------------------------- */

ReactTCEQK::ReactTCEQK(SPARTA *sparta, int narg, char **arg) :
  ReactBird(sparta, narg, arg) {}

/* ---------------------------------------------------------------------- */

void ReactTCEQK::init()
{
  if (!collide || strcmp(collide->style,"vss") != 0)
    error->all(FLERR,"React tce/qk can only be used with collide vss");

  ReactBird::init();

  // do not allow recombination reactions for now

  for (int i = 0; i < nlist; i++)
    if (rlist[i].active && rlist[i].type == RECOMBINATION)
      error->all(FLERR,
                 "React qk does not currently support recombination reactions");
}

/* ---------------------------------------------------------------------- */

int ReactTCEQK::attempt(Particle::OnePart *ip, Particle::OnePart *jp,
                        double pre_etrans, double pre_erot, double pre_evib,
                        double &post_etotal, int &kspecies)
{
  double pre_etotal;
  double e_excess;
  double ecc;
  int reaction;

  int isp = ip->ispecies;
  int jsp = jp->ispecies;

  int n = reactions[isp][jsp].n;

  if (n == 0) return 0;
  int *list = reactions[isp][jsp].list;

  // loop over possible reactions for these 2 species

  for (int i = 0; i < n; i++) {
    OneReaction *r = &rlist[list[i]];

    // ignore energetically impossible reactions

    pre_etotal = pre_etrans + pre_erot + pre_evib;

    ecc = pre_etotal;

    e_excess = ecc - r->coeff[1];
    if (e_excess <= 0.0) continue;

    // compute probability of reaction

    if (r->style == ARRHENIUS)
      reaction = attempt_tce(ip,jp,r,
                             pre_etrans,pre_erot,
                             pre_evib,post_etotal,kspecies);
    else if (r->style == QUANTUM)
      reaction = attempt_qk(ip,jp,r,
                            pre_etrans,pre_erot,
                            pre_evib,post_etotal,kspecies);

    if (reaction) tally_reactions[list[i]]++;
  }

  return 0;
}

/* ---------------------------------------------------------------------- */

int ReactTCEQK::attempt_tce(Particle::OnePart *ip, Particle::OnePart *jp,
                            OneReaction *r,
                            double pre_etrans, double pre_erot, double pre_evib,
                            double &post_etotal, int &kspecies)
{
  Particle::Species *species = particle->species;
  int isp = ip->ispecies;
  int jsp = jp->ispecies;

  double pre_ave_rotdof = (species[isp].rotdof + species[jsp].rotdof)/2.;

  // probablity to compare to reaction probability

  double react_prob = 0.0;
  double random_prob = random->uniform();

  double pre_etotal = pre_etrans + pre_erot + pre_evib;

  double ecc = pre_etrans;
  if (pre_ave_rotdof > 0.1) ecc += pre_erot*r->coeff[0]/pre_ave_rotdof;

  double e_excess = ecc - r->coeff[1];
  if (e_excess > 0.0) return 0;

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

  return 0;
}

/* ---------------------------------------------------------------------- */

int ReactTCEQK::attempt_qk(Particle::OnePart *ip, Particle::OnePart *jp,
                           OneReaction * r,
                           double pre_etrans, double pre_erot, double pre_evib,
                           double &post_etotal, int &kspecies)
{
  double prob,evib,inverse_kT;
  int iv,ilevel,maxlev,limlev;
  int mspec,aspec;

  Particle::Species *species = particle->species;
  int isp = ip->ispecies;
  int jsp = jp->ispecies;

  double pre_ave_rotdof = (species[isp].rotdof + species[jsp].rotdof)/2.0;

  double omega = collide->extract(isp,jsp,"omega");

  // probablity to compare to reaction probability

  double react_prob = 0.0;
  double random_prob = random->uniform();

  double pre_etotal = pre_etrans + pre_erot + pre_evib;

  double ecc = pre_etrans;
  if (pre_ave_rotdof > 0.1) ecc += pre_erot*r->coeff[0]/pre_ave_rotdof;

  double e_excess = ecc - r->coeff[1];
  if (e_excess > 0.0) return 0;

  // compute probability of reaction

  inverse_kT = 1.0 / (update->boltz * species[isp].vibtemp[0]);

  switch (r->type) {
  case DISSOCIATION:
    {
      ecc = pre_etrans + ip->evib;
      maxlev = static_cast<int> (ecc * inverse_kT);
      limlev = static_cast<int> (fabs(r->coeff[1]) * inverse_kT);
      if (maxlev > limlev) react_prob = 1.0;
      break;
    }
  case EXCHANGE:
    {
      if (r->coeff[4] < 0.0 && species[isp].rotdof > 0) {

        // endothermic reaction

        ecc = pre_etrans + ip->evib;
        maxlev = static_cast<int> (ecc * inverse_kT);
        if (ecc > r->coeff[1]) {

          // PROB is the probability ratio of eqn (5.61)

          prob = 0.0;
          do {
            iv =  static_cast<int> (random->uniform()*(maxlev+0.99999999));
            evib = static_cast<double> (iv / inverse_kT);
            if (evib < ecc) react_prob = pow(1.0-evib/ecc,1.5-omega);
          } while (random->uniform() < react_prob);

          ilevel = static_cast<int> (fabs(r->coeff[4]) * inverse_kT);
          if (iv >= ilevel) react_prob = 1.0;
        }

       } else if (r->coeff[4] > 0.0 && species[isp].rotdof > 0) {

        ecc = pre_etrans + ip->evib;

        // mspec = post-collision species of the particle
        // aspec = post-collision species of the atom

        mspec = r->products[0];
        aspec = r->products[1];

        if (species[mspec].rotdof < 2.0)  {
          mspec = r->products[1];
            aspec = r->products[0];
        }

        // potential post-collision energy

        ecc += r->coeff[4];
        maxlev = static_cast<int> (ecc * inverse_kT);
        do {
          iv = random->uniform()*(maxlev+0.99999999);
          evib = static_cast<double>
            (iv * update->boltz*species[mspec].vibtemp[0]);
          if (evib < ecc) prob = pow(1.0-evib/ecc,1.5 - r->coeff[6]);
        } while (random->uniform() < prob);

        ilevel = static_cast<int>
          (fabs(r->coeff[4]/update->boltz/species[mspec].vibtemp[0]));
        if (iv >= ilevel) react_prob = 1.0;
      }

      break;
    }

  default:
    error->one(FLERR,"Unknown outcome in reaction");
    break;
  }

  // test against random number to see if this reaction occurs
  // if it does, reset species of I,J and optional K to product species
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
    ip->ispecies = r->products[0];
    jp->ispecies = r->products[1];

    post_etotal = pre_etotal + r->coeff[4];
    if (r->nproduct > 2) kspecies = r->products[2];
    else kspecies = -1;

    return 1;
  }

  return 0;
}

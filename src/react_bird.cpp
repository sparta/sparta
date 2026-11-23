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

#include <cmath>
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include <string>
#include "react_bird.h"
#include "input.h"
#include "collide.h"
#include "update.h"
#include "particle.h"
#include "comm.h"
#include "modify.h"
#include "fix.h"
#include "fix_ambipolar.h"
#include "random_knuth.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{DISSOCIATION,EXCHANGE,IONIZATION,RECOMBINATION};  // other react files
enum{ARRHENIUS,QUANTUM};                               // other react files
enum{NONE,DISCRETE,SMOOTH};                            // several files

#define MAXREACTANT 2
#define MAXPRODUCT 3
#define MAXCOEFF 7               // 5 in file, extra for pre-computation

#define MAXLINE 1024
#define DELTALIST 16

// detailed-balance table self-check thresholds: a small drift between the
// two calibration temperatures is grid discretization (warn), a gross one
// means the reverse rate cannot reproduce k_f/Keq at all (fatal)

#define DB_DRIFT_WARN 0.02
#define DB_DRIFT_FATAL 0.25

/* ----------------------------------------------------------------------
   return 1 if the two small unordered species-index sets are identical
   used to match a reverse reaction with its forward partner
------------------------------------------------------------------------- */

static int set_match(int *a, int na, int *b, int nb)
{
  if (na != nb) return 0;
  int used[MAXPRODUCT] = {0};
  for (int i = 0; i < na; i++) {
    int found = 0;
    for (int j = 0; j < nb; j++) {
      if (!used[j] && a[i] == b[j]) { used[j] = 1; found = 1; break; }
    }
    if (!found) return 0;
  }
  return 1;
}

/* ---------------------------------------------------------------------- */

ReactBird::ReactBird(SPARTA *sparta, int narg, char **arg) :
  React(sparta, narg, arg)
{
  if (narg != 2) error->all(FLERR,"Illegal react tce or qk command");

  nlist = maxlist = 0;
  rlist = NULL;
  readfile(arg[1]);
  check_duplicate();

  tally_reactions = new bigint[nlist];
  tally_reactions_all = new bigint[nlist];
  tally_flag = 0;
  tce_bounds_checked = 0;
  reverse_off_warned = 0;
  reverse_modes_checked = 0;

  reactions = NULL;
  list_ij = NULL;
  sp2recomb_ij = NULL;

  mtab = NULL;
  mtab_num = NULL;
  mtab_du = NULL;
  mtab_n = NULL;
  mtab_nlist = 0;

  keqfits = NULL;
  nkeqfits = 0;
  generated_flag = 0;
  tce_style = 0;
}

/* ---------------------------------------------------------------------- */

ReactBird::ReactBird(SPARTA *sparta) : React(sparta)
{
  rlist = NULL;
  reactions = NULL;
  list_ij = NULL;
  sp2recomb_ij = NULL;

  mtab = NULL;
  mtab_num = NULL;
  mtab_du = NULL;
  mtab_n = NULL;
  mtab_nlist = 0;

  keqfits = NULL;
  nkeqfits = 0;
  generated_flag = 0;
  tce_style = 0;
  tally_reactions = NULL;
  tally_reactions_all = NULL;
  tce_bounds_checked = 0;
  reverse_off_warned = 0;
  reverse_modes_checked = 0;
}

/* ---------------------------------------------------------------------- */

ReactBird::~ReactBird()
{
  if (copy) return;

  free_micro_tables();

  delete [] tally_reactions;
  delete [] tally_reactions_all;

  if (rlist) {
    for (int i = 0; i < maxlist; i++) {
      for (int j = 0; j < rlist[i].nreactant; j++)
        delete [] rlist[i].id_reactants[j];
      for (int j = 0; j < rlist[i].nproduct; j++)
        delete [] rlist[i].id_products[j];
      delete [] rlist[i].id_reactants;
      delete [] rlist[i].id_products;
      delete [] rlist[i].reactants;
      delete [] rlist[i].products;
      delete [] rlist[i].coeff;
      delete [] rlist[i].id;
    }
  }
  memory->destroy(rlist);

  memory->destroy(reactions);
  memory->destroy(list_ij);
  memory->destroy(sp2recomb_ij);
  memory->sfree(keqfits);
}

/* ---------------------------------------------------------------------- */

void ReactBird::init()
{
  tally_flag = 0;
  for (int i = 0; i < nlist; i++) tally_reactions[i] = 0;

  // convert species IDs to species indices
  // flag reactions as active/inactive depending on whether all species exist
  // mark recombination reactions inactive if recombflag_user = 0,
  //   after resolving their species, and remember which ones so that
  //   generate_reverses() still sees them as providing a reverse

  int *recomb_off = new int[nlist];

  for (int m = 0; m < nlist; m++) {
    recomb_off[m] = 0;
    OneReaction *r = &rlist[m];
    r->active = 1;
    r->keq_flag = 0;
    for (int j = 0; j < 5; j++) r->keq_resid_coeff[j] = 0.0;

    // auto-generated reverse reactions are inert when reverse_auto is off

    if (r->generated && !reverse_auto) {
      r->active = 0;
      continue;
    }

    for (int i = 0; i < r->nreactant; i++) {
      r->reactants[i] = particle->find_species(r->id_reactants[i]);
      if (r->reactants[i] < 0) {
        r->active = 0;
        break;
      }
    }

    for (int i = 0; i < r->nproduct; i++) {
      r->products[i] = particle->find_species(r->id_products[i]);
      if (r->products[i] < 0) {

        // special case: recombination reaction with 2nd product = atom/mol

        if (r->type == RECOMBINATION && i == 1) {
          if (strcmp(r->id_products[i],"atom") == 0) {
            r->products[i] = -1;
            continue;
          } else if (strcmp(r->id_products[i],"mol") == 0) {
            r->products[i] = -2;
            continue;
          }
        }

        r->active = 0;
        break;
      }
    }

    if (r->active && r->type == RECOMBINATION && recombflag_user == 0) {
      r->active = 0;
      recomb_off[m] = 1;
    }
  }

  // detailed-balance reverse rates are implemented only by the tce styles
  //   with the total-energy model (react_modify partial_energy no)
  // otherwise B lines read from a file are switched off with one warning,
  //   so a file that carries a reverse for every forward reaction (e.g.
  //   data/air.tce) still runs forward-only, as it did before B lines
  //   existed; react_modify reverse auto is an explicit request, so an error

  if (!tce_style || partialEnergy) {
    if (reverse_auto) {
      if (!tce_style)
        error->all(FLERR,"React_modify reverse auto requires react tce");
      error->all(FLERR,"React_modify reverse auto requires "
                 "react_modify partial_energy no");
    }
    int noff = 0;
    for (int m = 0; m < nlist; m++)
      if (rlist[m].active && rlist[m].reverse) {
        rlist[m].active = 0;
        noff++;
      }
    // warn on the first run, and again only if the count changes
    if (noff && noff != reverse_off_warned && comm->me == 0) {
      char str[256];
      snprintf(str,sizeof(str),"%d reverse (B-style) reaction(s) "
               "deactivated: detailed-balance reverse rates require "
               "react tce with react_modify partial_energy no",noff);
      error->warning(FLERR,str);
    }
    reverse_off_warned = noff;
  }

  // auto-generate reverse (B-style) reactions for eligible forward
  // reactions (react_modify reverse auto)

  if (reverse_auto) generate_reverses(recomb_off);
  delete [] recomb_off;

  // issue #472: pair each reverse (detailed-balance, B-style) reaction with
  //   its forward partner and seed the backward Arrhenius coefficients.
  // supported pairs:
  //   EXCHANGE B      <-> EXCHANGE forward with reactants/products swapped
  //   RECOMBINATION B <-> DISSOCIATION forward with the same third body:
  //                       B: A + B -> AB + M  pairs with  F: AB + M -> A + B + M
  // following Bird94 sec 6.6 and Boyd & Schwartzentruber sec 7.5.2-7.5.3, with
  //   dHf = F.coeff[4] (signed reaction energy, negative for endothermic F):
  //   Ea_B  = Ea_F + dHf        (barrier seen from the product side)
  //   dHr_B = -dHf              (reverse reaction energy)
  //   b_B, z_B inherit from F   (temperature exponent, effective DOF)
  //   A_B(T) = A_F * q_react,F(T)/q_prod,F(T)
  // the equilibrium-constant exponential cancels against the shifted barrier,
  //   leaving only the partition-function ratio, which ReactTCE::attempt()
  //   applies at run time using the cell temperature React::tgas; for a
  //   recombination pair the ratio has one net translational factor (units of
  //   volume), converting the m^3/s dissociation prefactor into the m^6/s
  //   recombination prefactor.
  // here we store the RAW forward A into B.coeff[2] so the standard TCE
  //   transform below yields B's geometric prefactor G_B * A_F.
  // guarded by initflag so the seeding happens exactly once (like the transform)

  for (int m = 0; m < nlist; m++) {
    OneReaction *b = &rlist[m];
    if (!b->reverse) continue;

    // skip inactive/already-initialized entries BEFORE validating the
    // reaction type: a reaction file may legitimately contain B lines for
    // species absent from this run, and those are deactivated at read time
    // like every other reaction rather than being an error

    if (!b->active || b->initflag) continue;

    if (b->type != EXCHANGE && b->type != RECOMBINATION) {
      print_reaction(b);
      error->all(FLERR,"Reverse (B-style) reactions support exchange "
                 "reactions (including dissociative recombination "
                 "AB+ + e -> A + B written as an E-style line, paired "
                 "with associative ionization) and recombination "
                 "reactions paired with dissociation; the reverse of "
                 "electron-impact ionization must be supplied explicitly");
    }
    if (b->type == RECOMBINATION && b->products[1] < 0) {
      print_reaction(b);
      error->all(FLERR,"Reverse (B-style) recombination requires an "
                 "explicit third-body species (not atom/mol), so it can "
                 "pair with the per-partner forward dissociation reaction");
    }

    int found = -1;
    if (b->type == EXCHANGE) {
      // an EXCHANGE B line also pairs with a 2-product IONIZATION
      // forward (associative ionization A + B -> AB+ + e): its reverse
      // is dissociative recombination AB+ + e -> A + B, a 2 -> 2
      // reaction whose Saha-type equilibrium constant comes from the
      // same partition functions (the free electron contributes spin
      // degeneracy 2, see partition_function)

      for (int f = 0; f < nlist; f++) {
        OneReaction *r = &rlist[f];
        if (f == m || !r->active || r->reverse) continue;
        // only an Arrhenius forward can seed a detailed-balance reverse:
        // a QK reaction's coefficients are not (A, Ea, eta) and would be
        // silently misread as such (cf. generate_reverses)
        if (r->style != ARRHENIUS) continue;
        if (r->type != EXCHANGE &&
            !(r->type == IONIZATION && r->nproduct == 2)) continue;
        if (set_match(r->reactants,r->nreactant,b->products,b->nproduct) &&
            set_match(r->products,r->nproduct,b->reactants,b->nreactant)) {
          found = f;
          break;
        }
      }
    } else {                    // RECOMBINATION paired with DISSOCIATION
      for (int f = 0; f < nlist; f++) {
        OneReaction *r = &rlist[f];
        if (f == m || !r->active || r->reverse) continue;
        if (r->style != ARRHENIUS) continue;   // see note above
        if (r->type != DISSOCIATION) continue;
        // F: AB + M -> A + B + M   vs   B: A + B -> AB + M
        if (r->reactants[0] != b->products[0]) continue;   // AB
        if (r->reactants[1] != b->products[1]) continue;   // same third body
        if (r->products[2] != r->reactants[1]) continue;   // F third body sane
        int ab_pair[2] = {r->products[0],r->products[1]};
        if (!set_match(ab_pair,2,b->reactants,b->nreactant)) continue;
        found = f;
        break;
      }
    }
    if (found < 0) {
      print_reaction(b);
      error->all(FLERR,
                 "No forward partner found for reverse (B-style) reaction");
    }

    OneReaction *f = &rlist[found];
    b->reverse_partner = found;

    b->coeff[0] = f->coeff[0];                // effective internal DOF
    b->coeff[1] = f->coeff[1] + f->coeff[4];  // Ea_B = Ea_F + dHf
    // raw A_F: the forward's coeff[2] already holds its TCE-transformed
    // prefactor if it was initialized in an earlier run (a reverse
    // generated or activated after that run), so use the saved raw value

    double araw_f = f->initflag ? f->a_raw : f->coeff[2];
    b->coeff[2] = araw_f;                     // raw A_F (scaled at run time)
    b->coeff[3] = 0.0;                        // T dependence handled by the
    b->reverse_bf = f->coeff[3];              //   microcanonical detailed-
    b->reverse_A = araw_f;                    //   balance tables (see
                                              //   build_db_table and
                                              //   build_db3_table), which
                                              //   also keep the backward
                                              //   probability integrable for
                                              //   any forward b_F (a seeded
                                              //   b_B <= -(z+3/2) would make
                                              //   the standard TCE energy
                                              //   factor non-integrable, e.g.
                                              //   for atom + atom
                                              //   recombination); with an
                                              //   external Keq fit the T^b_F
                                              //   factor is applied at the
                                              //   cell temperature instead
    b->coeff[4] = -f->coeff[4];               // reverse reaction energy
    if (b->coeff[1] < 0.0) b->coeff[1] = 0.0; // clamp small negative barrier
  }

  // read equilibrium-constant curve fits (react_modify keq_file) and
  // assign them to the reverse reactions whose forward partner they
  // describe; a matched reaction evaluates k_b = k_f/Keq_fit at the cell
  // temperature instead of using the internal partition functions

  read_keq_file();
  assign_keq_fits();

  // modify Arrhenius coefficients for TCE model
  // C1,C2 Bird 94, p 127
  // initflag logic insures only done once per reaction

  Particle::Species *species = particle->species;

  for (int m = 0; m < nlist; m++) {
    OneReaction *r = &rlist[m];
    if (!r->active) continue;
    if (r->initflag) continue;
    r->initflag = 1;

    int isp = r->reactants[0];
    int jsp = r->reactants[1];

    // symmetry parameter

    double epsilon = 1.0;
    if (isp == jsp) epsilon = 2.0;

    double diam = collide->extract(isp,jsp,"diam");
    double omega = collide->extract(isp,jsp,"omega");
    double tref = collide->extract(isp,jsp,"tref");

    // double pre_ave_vibdof = (species[isp].vibdof + species[jsp].vibdof)/2.0;
    double mr = species[isp].mass * species[jsp].mass /
        (species[isp].mass + species[jsp].mass);
    double sigma = MY_PI*diam*diam;

    // read effective internal DOFs participating in the reaction

    //double z = r->coeff[0];

    // add additional coeff for effective DOF

    r->a_raw = r->coeff[2];
    double c1 = MY_PIS*epsilon*r->coeff[2]/(2.0*sigma) *
      sqrt(mr/(2.0*update->boltz*tref)) *
      pow(tref,1.0-omega)/pow(update->boltz,r->coeff[3]-1.0+omega);
    double c2 = r->coeff[3] - 1.0 + omega;

    r->coeff[2] = c1;
    r->coeff[5] = omega;

    // add additional coeff for post-collision effective omega
    // mspec = post-collision species of the particle
    // aspec = post-collision species of the atom

    double momega,aomega;

    if (r->nproduct > 1) {
      int mspec = r->products[0];
      int aspec = r->products[1];

      if (species[mspec].rotdof < 2.0) {
        mspec = r->products[1];
        aspec = r->products[0];
      }

      int ncount = 0;
      if (mspec >= 0) {
        momega = collide->extract(mspec,mspec,"omega");
        ncount++;
      } else momega = 0.0;
      if (aspec >= 0) {
        aomega = collide->extract(aspec,aspec,"omega");
        ncount++;
      } else aomega = 0.0;

      r->coeff[6] = (momega+aomega) / ncount;

    } else {
      int mspec = r->products[0];
      momega = collide->extract(mspec,mspec,"omega");
      r->coeff[6] = momega;
    }
  }

  // a reverse reaction derived from a forward reaction whose temperature
  //   exponent is below the TCE bound (see check_tce_bounds()) has no
  //   consistent rate: the forward probability is clamped and does not
  //   reproduce k_f, so the detailed-balance table cannot be calibrated
  //   against it; switch the reverse off, check_tce_bounds() warns about
  //   the forward reaction

  for (int m = 0; m < nlist; m++) {
    OneReaction *b = &rlist[m];
    if (!b->active || !b->reverse || b->reverse_partner < 0) continue;
    OneReaction *f = &rlist[b->reverse_partner];
    if (f->active && f->coeff[3] > tce_eta_bound(f)) continue;
    b->active = 0;
  }

  // the equilibrium constant of every active reverse reaction comes from
  //   partition functions: a rotating species with no rotfile data would
  //   silently get q_rot = 1, a constant error in the reverse rate that the
  //   table drift check cannot detect when it cancels in T, so require it

  for (int m = 0; m < nlist; m++) {
    OneReaction *b = &rlist[m];
    if (!b->active || !b->reverse) continue;
    for (int i = 0; i < b->nreactant + b->nproduct; i++) {
      int isp = (i < b->nreactant) ? b->reactants[i] :
        b->products[i-b->nreactant];
      if (isp < 0) continue;
      Particle::Species *sp = &particle->species[isp];
      if (sp->rotdof == 0) continue;
      int ok = 0;
      if (sp->rotdof == 2 && sp->nrottemp >= 1 && sp->rottemp[0] > 0.0)
        ok = 1;
      else if (sp->rotdof == 3 && sp->nrottemp == 3 &&
               sp->rottemp[0] > 0.0 && sp->rottemp[1] > 0.0 &&
               sp->rottemp[2] > 0.0) ok = 1;
      if (!ok) {
        char str[MAXLINE+256];
        snprintf(str,sizeof(str),"Reverse reaction %s needs the rotational "
                 "temperature(s) of species %s for its equilibrium "
                 "constant: add a rotfile to the species command",
                 b->id,sp->id);
        error->all(FLERR,str);
      }
    }
  }

  // set reverse_active if any active reverse reaction uses an external
  // Keq curve fit: only the residual thermal correction R(T) needs the
  // per-cell temperature at run time; the energy-resolved shape of every
  // reverse reaction comes from its temperature-free detailed-balance table

  reverse_active = 0;
  for (int m = 0; m < nlist; m++)
    if (rlist[m].active && rlist[m].reverse && rlist[m].keq_flag)
      reverse_active = 1;

  // set reverse_recomb_active if any active reverse recombination uses a
  // 3-body detailed-balance table (with or without an external Keq fit):
  // the KOKKOS collision loop fetches the 3rd particle's electronic energy
  // for the reaction probability only when this is set

  reverse_recomb_active = 0;
  for (int m = 0; m < nlist; m++)
    if (rlist[m].active && rlist[m].reverse &&
        rlist[m].type == RECOMBINATION &&
        rlist[m].reverse_partner >= 0)
      reverse_recomb_active = 1;

  // the recombination boost tests recombination in 1/rboost of the
  // collisions with rboost times the probability.  The derived (B-style)
  // recombination probability is unbounded for soft collisions and is
  // capped at 1 at run time, so the boost moves rate into the cap and
  // lowers the realized recombination rate (by about 25% at 1e24 m^-3 and
  // 6000 K with the default 1000).  The default is kept; warn

  if (reverse_recomb_active && recomb_boost > 1.0 && comm->me == 0)
    error->warning(FLERR,"Reverse (B-style) recombination reactions with "
                   "react_modify rboost > 1: the boosted probability is "
                   "capped at 1, which lowers the recombination rate below "
                   "its detailed-balance value; use react_modify rboost 1");

  // count possible active reactions for each species pair
  // include J,I reactions in I,J list and vice versa
  // this allows collision pair I,J to be in either order in Collide

  memory->destroy(reactions);
  int nspecies = particle->nspecies;
  reactions = memory->create(reactions,nspecies,nspecies,
                             "react/bird:reactions");

  for (int i = 0; i < nspecies; i++)
    for (int j = 0; j < nspecies; j++)
      reactions[i][j].n = 0;

  int n = 0;
  for (int m = 0; m < nlist; m++) {
    OneReaction *r = &rlist[m];
    if (!r->active) continue;
    int i = r->reactants[0];
    int j = r->reactants[1];
    reactions[i][j].n++;
    n++;
    if (i == j) continue;
    reactions[j][i].n++;
    n++;
  }

  // allocate list_IJ = contiguous list of reactions for each IJ pair

  memory->destroy(list_ij);
  memory->create(list_ij,n,"react/bird:list_ij");

  // reactions[i][j].list = pointer into full list_ij vector

  int offset = 0;
  for (int i = 0; i < nspecies; i++)
    for (int j = 0; j < nspecies; j++) {
      reactions[i][j].list = &list_ij[offset];
      offset += reactions[i][j].n;
    }

  // reactions[i][j].list = indices of reactions for each species pair
  // include J,I reactions in I,J list and vice versa

  for (int i = 0; i < nspecies; i++)
    for (int j = 0; j < nspecies; j++)
      reactions[i][j].n = 0;

  for (int m = 0; m < nlist; m++) {
    OneReaction *r = &rlist[m];
    if (!r->active) continue;
    int i = r->reactants[0];
    int j = r->reactants[1];
    reactions[i][j].list[reactions[i][j].n++] = m;
    if (i == j) continue;
    reactions[j][i].list[reactions[j][i].n++] = m;
  }

  // set recombflag = 0/1 if any recombination reactions are defined & active
  // check for user disabling them is at top of this method

  recombflag = 0;
  for (int m = 0; m < nlist; m++) {
    if (!rlist[m].active) continue;
    if (rlist[m].type == RECOMBINATION) recombflag = 1;
  }

  if (!recombflag) return;

  // count how many IJ pairs have a recombination reaction
  // allocate sp2recomb_ij = contiguous list of reactions
  //   for all species for each IJ pair that has a recombination reaction

  OneReaction *r;

  int nij = 0;
  for (int i = 0; i < nspecies; i++)
    for (int j = 0; j < nspecies; j++) {
      int n = reactions[i][j].n;
      int *list = reactions[i][j].list;
      for (int m = 0; m < n; m++) {
        r = &rlist[list[m]];
        if (r->type == RECOMBINATION) {
          nij++;
          break;
        }
      }
    }

  memory->destroy(sp2recomb_ij);
  memory->create(sp2recomb_ij,nij*nspecies,"react/bird:sp2recomb_ij");

  // reactions[i][j].sp2recomb = pointer into full sp2recomb_ij vector

  offset = 0;
  for (int i = 0; i < nspecies; i++)
    for (int j = 0; j < nspecies; j++) {
      int n = reactions[i][j].n;
      int *list = reactions[i][j].list;
      if (offset < nij*nspecies)
        reactions[i][j].sp2recomb = &sp2recomb_ij[offset];
      else
        reactions[i][j].sp2recomb = NULL; // Needed for Kokkos
      for (int m = 0; m < n; m++) {
        r = &rlist[list[m]];
        if (r->type == RECOMBINATION) {
          offset += nspecies;
          break;
        }
      }
    }

  // loop over species K for each IJ pair
  // if the IJ pair has any recombination reactions,
  //   then fill in its reactions[i][j].sp2recomb entries,
  //   which are indices into rlist of specific recomb reaction
  //   for each 3rd particle species
  // if IJ pair has no recombination reactions, then do NOT set sp2recomb vec
  // matching recombination reaction is one most specific to species K
  // 4 levels of specificity from most to least, in 4 inner loops
  //   explicit K, K = atom/mol, any K, no match at all (sp2recomb = -1)

  int m;

  for (int i = 0; i < nspecies; i++)
    for (int j = 0; j < nspecies; j++) {
      int n = reactions[i][j].n;
      int *list = reactions[i][j].list;
      for (int k = 0; k < nspecies; k++) {
        for (m = 0; m < n; m++) {
          r = &rlist[list[m]];
          if (r->type != RECOMBINATION) continue;
          if (r->nproduct != 2 || r->products[1] < 0) continue;
          if (r->products[1] == k) {
            reactions[i][j].sp2recomb[k] = list[m];
            break;
          }
        }
        if (m < n) continue;

        for (m = 0; m < n; m++) {
          r = &rlist[list[m]];
          if (r->type != RECOMBINATION) continue;
          if (r->nproduct != 2 || r->products[1] >= 0) continue;
          if (r->products[1] == -1 && particle->species[k].vibdof == 0) {
            reactions[i][j].sp2recomb[k] = list[m];
            break;
          }
          if (r->products[1] == -2 && particle->species[k].vibdof > 0) {
            reactions[i][j].sp2recomb[k] = list[m];
            break;
          }
        }
        if (m < n) continue;

        for (m = 0; m < n; m++) {
          r = &rlist[list[m]];
          if (r->type != RECOMBINATION) continue;
          if (r->nproduct != 1) continue;
          reactions[i][j].sp2recomb[k] = list[m];
          break;
        }
        if (m < n) continue;

        for (m = 0; m < n; m++) {
          r = &rlist[list[m]];
          if (r->type != RECOMBINATION) continue;
          reactions[i][j].sp2recomb[k] = -1;
          break;
        }
      }
    }
}

/* ----------------------------------------------------------------------
   auto-generate a reverse (B-style) reaction for every eligible active
     forward reaction (react_modify reverse auto):
     EXCHANGE A + B -> C + D        gains  C + D -> A + B (exchange B)
     DISSOCIATION AB + M -> A+B+M   gains  A + B -> AB + M (recomb B)
       (explicit third body only; wildcard recombination cannot pair
       with a per-partner dissociation rate)
   a forward reaction is skipped if any active reaction already provides
     its reverse: an explicit B line, an independently fitted reverse, or
     a wildcard recombination covering the same product and third body;
     a recombination switched off only by react_modify recomb no counts,
     so switching recombination back on later cannot duplicate it
   recomb_off[m] = 1 for file reactions switched off by recomb no
   generated entries are marked and become inert if reverse_auto is
     turned off before a later run
------------------------------------------------------------------------- */

void ReactBird::generate_reverses(int *recomb_off)
{
  if (generated_flag) return;
  generated_flag = 1;

  Particle::Species *species = particle->species;

  int nforward = nlist;
  int ngen = 0;

  for (int m = 0; m < nforward; m++) {
    OneReaction *f = &rlist[m];
    if (!f->active || f->reverse) continue;
    if (f->style != ARRHENIUS) continue;      // QK reactions are not eligible

    int rtype;
    int reactants[2],products[2];
    if (f->type == EXCHANGE) {
      rtype = EXCHANGE;
      reactants[0] = f->products[0]; reactants[1] = f->products[1];
      products[0] = f->reactants[0]; products[1] = f->reactants[1];
    } else if (f->type == DISSOCIATION) {
      if (f->nproduct != 3) continue;
      if (f->reactants[1] < 0 || f->products[2] != f->reactants[1]) continue;
      rtype = RECOMBINATION;
      reactants[0] = f->products[0]; reactants[1] = f->products[1];
      products[0] = f->reactants[0]; products[1] = f->reactants[1];
    } else if (f->type == IONIZATION && f->nproduct == 2) {
      // associative ionization A + B -> AB+ + e: the reverse is
      // dissociative recombination AB+ + e -> A + B, an exchange
      rtype = EXCHANGE;
      reactants[0] = f->products[0]; reactants[1] = f->products[1];
      products[0] = f->reactants[0]; products[1] = f->reactants[1];
    } else continue;

    // skip if any active reaction already provides this reverse

    int exists = 0;
    for (int k = 0; k < nlist; k++) {
      OneReaction *r2 = &rlist[k];
      if (r2->type != rtype) continue;
      if (!r2->active && !(k < nforward && recomb_off[k])) continue;
      if (rtype == EXCHANGE) {
        if (set_match(r2->reactants,2,reactants,2) &&
            set_match(r2->products,2,products,2)) { exists = 1; break; }
      } else {
        if (!set_match(r2->reactants,2,reactants,2)) continue;
        if (r2->products[0] != products[0]) continue;
        // a wildcard third body covers this one only if sp2recomb would
        // select it for it: atom for species without vibration, mol for
        // species with vibration (see the end of init())
        int k = products[1];
        if (r2->nproduct < 2 || r2->products[1] == k ||
            (r2->products[1] == -1 && species[k].vibdof == 0) ||
            (r2->products[1] == -2 && species[k].vibdof > 0)) {
          exists = 1;
          break;
        }
      }
    }
    if (exists) continue;

    // append the generated reverse to rlist

    if (nlist == maxlist) {
      maxlist += DELTALIST;
      rlist = (OneReaction *)
        memory->srealloc(rlist,maxlist*sizeof(OneReaction),"react/bird:rlist");
      for (int i = nlist; i < maxlist; i++) {
        OneReaction *r = &rlist[i];
        r->nreactant = r->nproduct = 0;
        r->id_reactants = new char*[MAXREACTANT];
        r->id_products = new char*[MAXPRODUCT];
        r->reactants = new int[MAXREACTANT];
        r->products = new int[MAXPRODUCT];
        r->coeff = new double[MAXCOEFF];
        r->id = NULL;
        r->reverse = 0;
        r->reverse_partner = -1;
        r->reverse_bf = 0.0;
        r->a_raw = 0.0;
        r->reverse_A = 0.0;
        r->generated = 0;
        r->keq_flag = 0;
      }
      f = &rlist[m];   // rlist may have moved
    }

    OneReaction *b = &rlist[nlist];
    // a generated recombination obeys react_modify recomb no like any
    // other recombination (later runs apply it in the loop at the top of
    // init(), this covers the run that generates it)
    b->active = (rtype == RECOMBINATION && recombflag_user == 0) ? 0 : 1;
    b->initflag = 0;
    b->type = rtype;
    b->style = ARRHENIUS;
    b->ncoeff = 5;
    b->nreactant = 2;
    b->nproduct = 2;
    b->reverse = 1;
    b->reverse_partner = -1;
    b->reverse_bf = 0.0;
    b->a_raw = 0.0;
    b->reverse_A = 0.0;
    b->generated = 1;
    b->keq_flag = 0;
    for (int i = 0; i < 5; i++) b->coeff[i] = 0.0;

    char idbuf[MAXLINE];
    for (int i = 0; i < 2; i++) {
      char *name = species[reactants[i]].id;
      b->reactants[i] = reactants[i];
      b->id_reactants[i] = new char[strlen(name)+1];
      strcpy(b->id_reactants[i],name);
      name = species[products[i]].id;
      b->products[i] = products[i];
      b->id_products[i] = new char[strlen(name)+1];
      strcpy(b->id_products[i],name);
    }
    sprintf(idbuf,"%s + %s --> %s + %s",
            species[reactants[0]].id,species[reactants[1]].id,
            species[products[0]].id,species[products[1]].id);
    b->id = new char[strlen(idbuf)+1];
    strcpy(b->id,idbuf);

    nlist++;
    ngen++;
  }

  // tally arrays were sized for the file reactions: regrow
  // (virtual: the KOKKOS variant owns them through dual views)

  if (ngen) {
    grow_tallies();
    for (int i = 0; i < nlist; i++) tally_reactions[i] = 0;
  }

  if (comm->me == 0) {
    if (screen)
      fprintf(screen,"Generated %d reverse reaction(s)\n",ngen);
    if (logfile)
      fprintf(logfile,"Generated %d reverse reaction(s)\n",ngen);
  }
}

/* ----------------------------------------------------------------------
   resize the reaction tally arrays to the current nlist, after
   auto-generation appended reactions; the KOKKOS variant overrides this
   because its arrays are owned by dual views
------------------------------------------------------------------------- */

void ReactBird::grow_tallies()
{
  delete [] tally_reactions;
  delete [] tally_reactions_all;
  tally_reactions = new bigint[nlist];
  tally_reactions_all = new bigint[nlist];
}

/* ----------------------------------------------------------------------
   read equilibrium-constant curve fits from react_modify keq_file:
     2-line entries, comments and blank lines allowed:
       A + B --> C + D          (the FORWARD reaction the fit describes)
       park c0 c1 c2 c3 c4
     ln Keq = c0/Z + c1 + c2 ln(Z) + c3 Z + c4 Z^2,  Z = 10000 K / T,
     with Keq = k_f/k_b in SI units (dimensionless for an exchange pair,
     1/m^3 for a dissociation/recombination pair), also under units cgs:
     fit_keq_residual() converts it to the run's units
   entries whose species are not all defined are skipped, like reaction
     files; called every init so react_modify changes take effect
------------------------------------------------------------------------- */

void ReactBird::read_keq_file()
{
  memory->sfree(keqfits);
  keqfits = NULL;
  nkeqfits = 0;
  if (!keq_file) return;

  int maxfits = 0;
  char line1[MAXLINE],line2[MAXLINE];

  FILE *keqfp = NULL;
  if (comm->me == 0) {
    keqfp = fopen(keq_file,"r");
    if (!keqfp) {
      char str[MAXLINE+64];
      snprintf(str,sizeof(str),"Cannot open Keq file %s",keq_file);
      error->one(FLERR,str);
    }
  }

  while (1) {
    int eof = 0;
    if (comm->me == 0) {

      // read a 2-line entry, skipping blank and comment lines

      char *ptr;
      do {
        ptr = fgets(line1,MAXLINE,keqfp);
      } while (ptr && (strspn(line1," \t\n\r") == strlen(line1) ||
                       line1[strspn(line1," \t")] == '#'));
      if (!ptr) eof = 1;
      else if (!fgets(line2,MAXLINE,keqfp)) eof = 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(line1,MAXLINE,MPI_CHAR,0,world);
    MPI_Bcast(line2,MAXLINE,MPI_CHAR,0,world);

    // parse the formula: species names split by "+" and "-->"

    if (nkeqfits == maxfits) {
      maxfits += DELTALIST;
      keqfits = (KeqFit *)
        memory->srealloc(keqfits,maxfits*sizeof(KeqFit),"react/bird:keqfits");
    }
    KeqFit *fit = &keqfits[nkeqfits];
    fit->nreactant = fit->nproduct = 0;
    fit->used = 0;

    int side = 0;
    int skip = 0;
    char *word = strtok(line1," \t\n\r");
    while (word) {
      if (strcmp(word,"+") == 0) {
        word = strtok(NULL," \t\n\r");
        continue;
      }
      if (strcmp(word,"-->") == 0) {
        side = 1;
        word = strtok(NULL," \t\n\r");
        continue;
      }
      int isp = particle->find_species(word);
      if (isp < 0) skip = 1;
      if (side == 0) {
        if (fit->nreactant == 2)
          error->all(FLERR,"Too many reactants in Keq file entry");
        fit->reactants[fit->nreactant++] = isp;
      } else {
        if (fit->nproduct == 3)
          error->all(FLERR,"Too many products in Keq file entry");
        fit->products[fit->nproduct++] = isp;
      }
      word = strtok(NULL," \t\n\r");
    }
    if (fit->nreactant != 2 || fit->nproduct < 2)
      error->all(FLERR,"Invalid reaction formula in Keq file");

    word = strtok(line2," \t\n\r");
    if (!word || (strcmp(word,"park") != 0 && strcmp(word,"PARK") != 0))
      error->all(FLERR,"Invalid fit style in Keq file (expected park)");
    for (int i = 0; i < 5; i++) {
      word = strtok(NULL," \t\n\r");
      if (!word) error->all(FLERR,"Keq file park fit requires 5 coefficients");
      fit->coeff[i] = input->numeric(FLERR,word);
    }

    // entries naming absent species are skipped like reaction-file lines

    if (!skip) nkeqfits++;
  }

  if (comm->me == 0) fclose(keqfp);
}

/* ----------------------------------------------------------------------
   assign each Keq fit to the reverse reaction whose FORWARD partner it
   describes; a matched reverse evaluates k_b = k_f/Keq_fit at the cell
   temperature in place of the internal partition functions
------------------------------------------------------------------------- */

void ReactBird::assign_keq_fits()
{
  if (!nkeqfits) return;

  for (int m = 0; m < nlist; m++) {
    OneReaction *b = &rlist[m];
    if (!b->active || !b->reverse || b->reverse_partner < 0) continue;
    OneReaction *f = &rlist[b->reverse_partner];

    for (int k = 0; k < nkeqfits; k++) {
      KeqFit *fit = &keqfits[k];
      if (fit->nreactant != f->nreactant) continue;
      if (fit->nproduct != f->nproduct) continue;
      if (!set_match(fit->reactants,fit->nreactant,
                     f->reactants,f->nreactant)) continue;
      if (!set_match(fit->products,fit->nproduct,
                     f->products,f->nproduct)) continue;
      b->keq_flag = 1;
      for (int i = 0; i < 5; i++) b->keq_coeff[i] = fit->coeff[i];
      fit->used = 1;
      fit_keq_residual(m);
      break;
    }
  }

  int nunused = 0;
  for (int k = 0; k < nkeqfits; k++)
    if (!keqfits[k].used) nunused++;
  if (nunused && comm->me == 0) {
    char str[160];
    sprintf(str,"%d of %d Keq fit(s) matched no active reverse reaction",
            nunused,nkeqfits);
    error->warning(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   fit the residual thermal correction of a reverse reaction matched to an
   external Keq curve fit:
     R(T) = Keq_statmech(T) / Keq_ext(T)
   The detailed-balance table (build_db_table / build_db3_table) already
   reproduces the reverse rate for the statistical-mechanics equilibrium
   constant Keq_statmech, so multiplying its per-collision probability by
   R(T) at run time makes the thermal average reproduce the EXTERNAL Keq
   while keeping the energy-resolved, microscopically-reversible collision-
   energy shape (only the -- usually small and smooth -- Keq discrepancy is
   thermal).  ln R is fit to the same Park form as Keq so both the CPU and
   KOKKOS runtimes evaluate it identically from stored coefficients, with no
   partition functions on the device.  Keq_statmech of the FORWARD reaction
   is q_products/q_reactants times exp(dHrxn/kT); a spectator third body
   cancels between the two products/reactants sets.
------------------------------------------------------------------------- */

void ReactBird::fit_keq_residual(int i)
{
  OneReaction *b = &rlist[i];
  OneReaction *f = &rlist[b->reverse_partner];
  const double boltz = update->boltz;

  // least-squares fit of ln R(T) over a log-spaced temperature range wide
  // enough to cover any cell temperature the reaction will see

  const int NP = 16;
  const double Tlo = 1000.0, Thi = 60000.0;
  double ata[5][5], atb[5];

  // the external fit is in SI (Keq in 1/m^3 for a dissociation/recombination
  // pair), the partition functions are in the run's units (per cm^3 in cgs):
  // convert the fit, one factor 1e-6 per net number-density power dn of the
  // forward reaction (dn = 1 for dissociation, 0 for exchange/ionization)

  const int cgs = (strcmp(update->unit_style,"cgs") == 0);
  const int dn = f->nproduct - f->nreactant;
  const double keq_units = cgs ? pow(1.0e-6,dn) : 1.0;
  double lnr_min = HUGE_VAL;
  for (int a = 0; a < 5; a++) {
    atb[a] = 0.0;
    for (int c = 0; c < 5; c++) ata[a][c] = 0.0;
  }

  for (int p = 0; p < NP; p++) {
    double T = Tlo * pow(Thi/Tlo, p/(double)(NP-1));

    double keq_sm = exp(f->coeff[4]/(boltz*T));
    for (int j = 0; j < f->nproduct; j++)
      keq_sm *= partition_function(f->products[j],T);
    for (int j = 0; j < f->nreactant; j++)
      keq_sm /= partition_function(f->reactants[j],T);

    double lnR = log(keq_sm) - log(keq_units*keq_eval(b->keq_coeff,T));
    if (fabs(lnR) < lnr_min) lnr_min = fabs(lnR);

    double Z = 10000.0/T;
    double basis[5] = {1.0/Z, 1.0, log(Z), Z, Z*Z};
    for (int a = 0; a < 5; a++) {
      atb[a] += basis[a]*lnR;
      for (int c = 0; c < 5; c++) ata[a][c] += basis[a]*basis[c];
    }
  }

  // solve the 5x5 normal equations by Gaussian elimination with partial
  // pivoting; store the coefficients so R(T) = keq_eval(keq_resid_coeff,T)

  double m[5][6];
  for (int a = 0; a < 5; a++) {
    for (int c = 0; c < 5; c++) m[a][c] = ata[a][c];
    m[a][5] = atb[a];
  }
  for (int col = 0; col < 5; col++) {
    int piv = col;
    for (int r = col+1; r < 5; r++)
      if (fabs(m[r][col]) > fabs(m[piv][col])) piv = r;
    for (int c = 0; c < 6; c++) {
      double t = m[col][c]; m[col][c] = m[piv][c]; m[piv][c] = t;
    }
    // a zero (or non-finite) pivot means the normal equations are singular
    // -- degenerate/insufficient sample points, or a Keq curve that made a
    // basis column non-finite.  Dividing by it would poison every
    // coefficient with inf/NaN, and a NaN keq_resid_coeff silently turns
    // the reverse probability into NaN at run time, which removes the whole
    // species pair from the chemistry without any diagnostic.  Fail here
    // instead, where the cause is still identifiable.
    double d = m[col][col];
    if (!(fabs(d) > 0.0) || !std::isfinite(d)) {
      char str[MAXLINE+128];
      sprintf(str,"Reverse reaction %s: external-Keq residual fit is "
              "singular (zero pivot); check the keq_file fit coefficients "
              "and the species partition-function data",b->id);
      error->all(FLERR,str);
    }
    for (int c = col; c < 6; c++) m[col][c] /= d;
    for (int r = 0; r < 5; r++) {
      if (r == col) continue;
      double fac = m[r][col];
      for (int c = col; c < 6; c++) m[r][c] -= fac*m[col][c];
    }
  }
  for (int a = 0; a < 5; a++) {
    if (!std::isfinite(m[a][5])) {
      char str[MAXLINE+128];
      sprintf(str,"Reverse reaction %s: external-Keq residual fit produced "
              "a non-finite coefficient; check the keq_file fit and the "
              "species partition-function data",b->id);
      error->all(FLERR,str);
    }
    b->keq_resid_coeff[a] = m[a][5];
  }

  // goodness-of-fit self-check: the Park basis reproduces a statistical-
  // mechanics / Park-fit Keq ratio to well under 1% in practice, but a
  // pathological external fit (different functional shape, or wildly
  // inconsistent partition functions) could defeat the 5-term fit; warn if
  // the residual R(T) is off by more than 2% anywhere in the fit window so
  // the reproduced Keq is not silently wrong

  double maxerr = 0.0;
  for (int p = 0; p < NP; p++) {
    double T = Tlo * pow(Thi/Tlo, p/(double)(NP-1));
    double keq_sm = exp(f->coeff[4]/(boltz*T));
    for (int j = 0; j < f->nproduct; j++)
      keq_sm *= partition_function(f->products[j],T);
    for (int j = 0; j < f->nreactant; j++)
      keq_sm /= partition_function(f->reactants[j],T);
    double lnR = log(keq_sm) - log(keq_units*keq_eval(b->keq_coeff,T));
    double err = fabs(keq_eval(b->keq_resid_coeff,T)/exp(lnR) - 1.0);
    // NOTE: written so that a NaN err is caught (NaN fails every ordered
    // comparison, so "err > maxerr" alone would silently ignore it)
    if (!(err <= maxerr)) maxerr = err;
  }
  // the residual fit absorbs any constant factor between the two Keq (a
  // units mismatch, e.g. a fit not in SI, would vanish into coefficient 1),
  // so check its size separately: the external and partition-function Keq
  // of the same pair agree somewhere in the window to well within a factor
  // of ten (they part only where a fit is extrapolated), while a units
  // mismatch is a constant factor at every temperature

  if (!(lnr_min <= log(10.0)) && comm->me == 0) {
    char str[MAXLINE+192];
    snprintf(str,sizeof(str),"Reverse reaction %s: the keq_file fit and the "
             "partition-function Keq differ by more than a factor of %g at "
             "every temperature in %g-%g K; check that the fit is in SI units "
             "and the species partition-function data",b->id,exp(lnr_min),
             Tlo,Thi);
    error->warning(FLERR,str);
  }

  if (!(maxerr <= 0.02) && comm->me == 0) {
    char str[MAXLINE+128];
    sprintf(str,"Reverse reaction %s: external-Keq residual fit is off by "
            "%g%% within %g-%g K; the reproduced Keq will be biased by that "
            "much (check the Keq file fit and species partition-function data)",
            b->id,100.0*maxerr,Tlo,Thi);
    error->warning(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   check that the temperature exponent eta = coeff[3] of each active
     Arrhenius reaction is within the bound for which the TCE reaction
     probability
       P = C1 * Gamma(z+5/2-omega) / Gamma(z+eta+3/2) *
           (Ec-Ea)^(eta-1+omega) * (1-Ea/Ec)^(z+3/2-omega)
     reproduces the Arrhenius rate of the reaction
   the bound is eta > -(z + 3/2), see tce_eta_bound() for z
     below it the TCE inversion of the Arrhenius rate does not exist
     and the argument of Gamma(z+eta+3/2) is non-positive
     ReactTCE::attempt() clamps the gamma function to 1.0e-6,
     so the reaction still runs but its rate is incorrect, and
     ReactBird::init() switches off any reverse reaction derived from it
   other properties of the probability, e.g. whether it diverges as Ec
     approaches Ea or infinity, do not change its average, so the rate
     is still reproduced; where the probability exceeds 1 it is clipped,
     which ReactTCE::end_of_run() warns about when it happens
   B-style reverse reactions are seeded with eta = 0, so never flagged
   all violating reactions are reported in one warning, on the first
     call only, so that it is not repeated for each run
   called from ReactTCE::init() and ReactTCEKokkos::init(),
     after ReactBird::init()
------------------------------------------------------------------------- */

void ReactBird::check_tce_bounds()
{
  if (tce_bounds_checked) return;
  tce_bounds_checked = 1;

  char str[MAXLINE+64];

  int nbad = 0;
  std::string mesg;

  for (int m = 0; m < nlist; m++) {
    OneReaction *r = &rlist[m];
    if (!r->active) continue;

    // QK reactions do not use the TCE probability

    if (r->style != ARRHENIUS) continue;

    double eta = r->coeff[3];
    double bound = tce_eta_bound(r);

    if (eta <= bound) {
      snprintf(str,sizeof(str),"\n  %s: temperature exponent %g <= %g",
               r->id,eta,bound);
      mesg += str;
      nbad++;
    }
  }

  if (nbad && comm->me == 0) {
    std::string head = std::to_string(nbad) + " TCE reaction(s) have a "
      "temperature exponent below the bound -(z+3/2) for which the TCE "
      "model can reproduce their Arrhenius rate. Their gamma function is "
      "non-positive and clamped, so their reaction rates are incorrect";
    if (partialEnergy == 0) head += " and any reverse reaction derived "
      "from them is switched off";
    head += ":";
    error->warning(FLERR,(head + mesg).c_str());
  }
}

/* ----------------------------------------------------------------------
   warn once, over all procs, about invalid TCE reaction probabilities
     recorded on this proc during the run by the tce styles
   negflag = 1 if a negative probability occurred
   bigindex = rlist index of the first reaction at which a species pair's
     summed probability exceeded 1, -1 if none
   called on all procs by ReactTCE::end_of_run() and
     ReactTCEKokkos::end_of_run()
------------------------------------------------------------------------- */

void ReactBird::report_tce_prob_flags(int negflag, int bigindex)
{
  int negall;
  MPI_Allreduce(&negflag,&negall,1,MPI_INT,MPI_MAX,world);

  // report the lowest-index reaction flagged by any proc

  int mine = (bigindex >= 0) ? bigindex : nlist;
  int bigall;
  MPI_Allreduce(&mine,&bigall,1,MPI_INT,MPI_MIN,world);

  if (comm->me == 0) {
    if (negall)
      error->warning(FLERR,"Negative TCE reaction probability occurred "
                     "during this run, check reaction file coefficients");
    if (bigall < nlist) {
      std::string mesg = "Summed TCE reaction probability exceeded 1.0 "
        "during this run, e.g. at reaction " +
        std::string(rlist[bigall].id) + ": the Arrhenius rates exceed "
        "the collision rate at some collision energies, where the "
        "reactions of that species pair are limited to the collision rate "
        "(each selected in proportion to its probability), so their rates "
        "are under-predicted; this is independent of timestep and fnum";
      error->warning(FLERR,mesg.c_str());
    }
  }
}

/* ----------------------------------------------------------------------
   lower bound -(z+3/2) on the temperature exponent of Arrhenius reaction R
     for which the TCE probability reproduces its rate, see
     check_tce_bounds()
   z = continuum internal DOF contributing to the collision energy Ec,
     constant for a given reaction:
     partial energy (rDOF): z = coeff[0]
     total energy: z = average rotational DOF of the reactants, plus the
       average vibrational DOF when vibstyle = smooth; the discrete
       vibrational and electronic ladders enter through the microcanonical
       energy-factor tables, not through z
------------------------------------------------------------------------- */

double ReactBird::tce_eta_bound(OneReaction *r)
{
  Particle::Species *species = particle->species;
  int isp = r->reactants[0];
  int jsp = r->reactants[1];

  double z;
  if (partialEnergy) z = r->coeff[0];
  else {
    z = 0.5 * (species[isp].rotdof + species[jsp].rotdof);
    if (collide->vibstyle == SMOOTH)
      z += 0.5 * (species[isp].vibdof + species[jsp].vibdof);
  }
  return -(z+1.5);
}

/* ----------------------------------------------------------------------
   return 1 if any recombination reactions are defined for species pair ISP,JSP
   else return 0
   called from Collide::init(), after React::init() has been performed
------------------------------------------------------------------------- */

int ReactBird::recomb_exist(int isp, int jsp)
{
  if (reactions[isp][jsp].sp2recomb) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   check active reactions that include ambi ion or electron especies
   their format must be correct to work with ambi_reset()
   called after init() from collide::init()
------------------------------------------------------------------------- */

void ReactBird::ambi_check()
{
  int flag;
  OneReaction *r;

  // fix ambipolar must exist since collide caller extracted ambi vector/array

  int ifix;
  for (ifix = 0; ifix < modify->nfix; ifix++)
    if (strcmp(modify->fix[ifix]->style,"ambipolar") == 0) break;
  FixAmbipolar *afix = (FixAmbipolar *) modify->fix[ifix];
  int especies = afix->especies;
  int *ions = afix->ions;

  // loop over active reactions

  for (int i = 0; i < nlist; i++) {
    r = &rlist[i];
    if (!r->active) continue;

    // skip reaction if no ambipolar ions or electrons as reactant or product
    // r->products[j] can be < 0 for atom or mol

    flag = 0;
    for (int j = 0; j < r->nreactant; j++)
      if (r->reactants[j] == especies || ions[r->reactants[j]]) flag = 1;
    for (int j = 0; j < r->nproduct; j++) {
      if (r->products[j] < 0) continue;
      if (r->products[j] == especies || ions[r->products[j]]) flag = 1;
    }
    if (!flag) continue;

    // dissociation must match one of these orders
    // D: AB + e -> A + e + B
    // D: AB+ + e -> A+ + e + B

    flag = 1;

    if (r->type == DISSOCIATION) {
      if (r->nreactant == 2 && r->nproduct == 3) {
        if (ions[r->reactants[0]] == 0 && r->reactants[1] == especies &&
            ions[r->products[0]] == 0 && r->products[1] == especies &&
            ions[r->products[2]] == 0) flag = 0;
        else if (ions[r->reactants[0]] == 1 && r->reactants[1] == especies &&
                 ions[r->products[0]] == 1 && r->products[1] == especies &&
                 ions[r->products[2]] == 0) flag = 0;

        // ambipolar ion as inert third body, dissociating species is neutral
        // D: AB + C+ -> A + C+ + B
        //   third body C+ is conserved as middle product (products[1])
        //   so it maps to the surviving J particle in react attempt()

        else if (ions[r->reactants[0]] == 0 && ions[r->reactants[1]] == 1 &&
                 ions[r->products[0]] == 0 &&
                 r->products[1] == r->reactants[1] &&
                 ions[r->products[2]] == 0) flag = 0;

        // ambipolar ion as inert third body, dissociating species is an ion
        // D: AB+ + C+ -> A+ + C+ + B

        else if (ions[r->reactants[0]] == 1 && ions[r->reactants[1]] == 1 &&
                 ions[r->products[0]] == 1 &&
                 r->products[1] == r->reactants[1] &&
                 ions[r->products[2]] == 0) flag = 0;
      }
    }

    // ionization with 3 products must match this
    // I: A + e -> A+ + e + e

    else if (r->type == IONIZATION && r->nproduct == 3) {
      if (r->nreactant == 2 && r->nproduct == 3) {
        if (ions[r->reactants[0]] == 0 && r->reactants[1] == especies &&
            ions[r->products[0]] == 1 && r->products[1] == especies &&
            r->products[2] == especies) flag = 0;
      }
    }

    // ionization with 2 products must match this
    // I: A + B -> AB+ + e

    else if (r->type == IONIZATION && r->nproduct == 2) {
      if (r->nreactant == 2 && r->nproduct == 2) {
        if (ions[r->reactants[0]] == 0 && ions[r->reactants[1]] == 0 &&
            ions[r->products[0]] == 1 && r->products[1] == especies) flag = 0;
      }
    }

    // exchange must match one of these
    // E: AB+ + e -> A + B
    // E: AB+ + C -> A + BC+
    // E: C + AB+ -> A + BC+

    else if (r->type == EXCHANGE) {
      if (r->nreactant == 2 && r->nproduct == 2) {
        if (ions[r->reactants[0]] == 1 && r->reactants[1] == especies &&
            ions[r->products[0]] == 0 && ions[r->products[1]] == 0) flag = 0;
        else if (ions[r->reactants[0]] == 1 && ions[r->reactants[1]] == 0 &&
            ions[r->products[0]] == 0 && ions[r->products[1]] == 1) flag = 0;
        else if (ions[r->reactants[0]] == 0 && ions[r->reactants[1]] == 1 &&
            ions[r->products[0]] == 0 && ions[r->products[1]] == 1) flag = 0;
      }
    }

    // recombination must match one of these
    // R: A+ + e -> A
    // R: A + B+ -> AB+
    // R: A+ + B -> AB+

    else if (r->type == RECOMBINATION) {
      if (r->nreactant == 2 && r->nproduct == 1) {
        if (ions[r->reactants[0]] == 1 && r->reactants[1] == especies &&
            ions[r->products[0]] == 0) flag = 0;
        else if (ions[r->reactants[0]] == 0 && ions[r->reactants[1]] == 1 &&
            ions[r->products[0]] == 1) flag = 0;
        else if (ions[r->reactants[0]] == 1 && ions[r->reactants[1]] == 0 &&
            ions[r->products[0]] == 1) flag = 0;
      }

      // ambipolar ion as inert third body
      // R: A + B + C+ -> AB + C+
      //   third body C+ is the 2nd product (products[1]), used only to
      //   select the recombination rate; it is not a colliding particle
      //   reactants and recombined molecule (products[0]) are not electrons

      else if (r->nreactant == 2 && r->nproduct == 2) {
        if (r->reactants[0] != especies && r->reactants[1] != especies &&
            r->products[0] != especies && ions[r->products[1]] == 1) flag = 0;
      }
    }

    // flag = 1 means unrecognized reaction

    if (flag) {
      print_reaction_ambipolar(r);
      error->all(FLERR,"Invalid ambipolar reaction");
    }
  }
}

/* ---------------------------------------------------------------------- */

void ReactBird::readfile(char *fname)
{
  int n,n1,n2,eof;
  char line1[MAXLINE],line2[MAXLINE];
  char copy1[MAXLINE],copy2[MAXLINE];
  char *word;
  OneReaction *r;

  // proc 0 opens file

  if (comm->me == 0) {
    fp = fopen(fname,"r");
    if (fp == NULL) {
      char str[128];
      snprintf(str,sizeof(str),"Cannot open reaction file %s",fname);
      error->one(FLERR,str);
    }
  }

  // read reactions one at a time and store their info in rlist

  while (1) {
    if (comm->me == 0) eof = readone(line1,line2,n1,n2);
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;

    MPI_Bcast(&n1,1,MPI_INT,0,world);
    MPI_Bcast(&n2,1,MPI_INT,0,world);
    MPI_Bcast(line1,n1,MPI_CHAR,0,world);
    MPI_Bcast(line2,n2,MPI_CHAR,0,world);

    if (nlist == maxlist) {
      maxlist += DELTALIST;
      rlist = (OneReaction *)
        memory->srealloc(rlist,maxlist*sizeof(OneReaction),"react/bird:rlist");
      for (int i = nlist; i < maxlist; i++) {
        r = &rlist[i];
        r->nreactant = r->nproduct = 0;
        r->id_reactants = new char*[MAXREACTANT];
        r->id_products = new char*[MAXPRODUCT];
        r->reactants = new int[MAXREACTANT];
        r->products = new int[MAXPRODUCT];
        r->coeff = new double[MAXCOEFF];
        r->id = NULL;
        r->reverse = 0;
        r->reverse_partner = -1;
        r->reverse_bf = 0.0;
        r->a_raw = 0.0;
        r->reverse_A = 0.0;
        r->generated = 0;
        r->keq_flag = 0;
      }
    }

    strcpy(copy1,line1);
    strcpy(copy2,line2);

    r = &rlist[nlist];
    r->initflag = 0;

    int side = 0;
    int species = 1;

    n = strlen(line1) - 1;
    r->id = new char[n+1];
    strncpy(r->id,line1,n);
    r->id[n] = '\0';

    word = strtok(line1," \t\n\r");

    while (1) {
      if (!word) {
        if (side == 0) {
          print_reaction(copy1,copy2);
          error->all(FLERR,"Invalid reaction formula in file");
        }
        if (species) {
          print_reaction(copy1,copy2);
          error->all(FLERR,"Invalid reaction formula in file");
        }
        break;
      }
      if (species) {
        species = 0;
        if (side == 0) {
          if (r->nreactant == MAXREACTANT) {
            print_reaction(copy1,copy2);
            error->all(FLERR,"Too many reactants in a reaction formula");
          }
          n = strlen(word) + 1;
          r->id_reactants[r->nreactant] = new char[n];
          strcpy(r->id_reactants[r->nreactant],word);
          r->nreactant++;
        } else {
          if (r->nproduct == MAXPRODUCT) {
            print_reaction(copy1,copy2);
            error->all(FLERR,"Too many products in a reaction formula");
          }
          n = strlen(word) + 1;
          r->id_products[r->nproduct] = new char[n];
          strcpy(r->id_products[r->nproduct],word);
          r->nproduct++;
        }
      } else {
        species = 1;
        if (strcmp(word,"+") == 0) {
          word = strtok(NULL," \t\n\r");
          continue;
        }
        if (strcmp(word,"-->") != 0) {
          print_reaction(copy1,copy2);
          error->all(FLERR,"Invalid reaction formula in file");
        }
        side = 1;
      }
      word = strtok(NULL," \t\n\r");
    }

    word = strtok(line2," \t\n\r");
    if (!word) {
      print_reaction(copy1,copy2);
      error->all(FLERR,"Invalid reaction type in file");
    }
    if (word[0] == 'D' || word[0] == 'd') r->type = DISSOCIATION;
    else if (word[0] == 'E' || word[0] == 'e') r->type = EXCHANGE;
    else if (word[0] == 'I' || word[0] == 'i') r->type = IONIZATION;
    else if (word[0] == 'R' || word[0] == 'r') r->type = RECOMBINATION;
    else {
      print_reaction(copy1,copy2);
      error->all(FLERR,"Invalid reaction type in file");
    }

    // check that reactant/product counts are consistent with type

    if (r->type == DISSOCIATION) {
      if (r->nreactant != 2 || r->nproduct != 3) {
        print_reaction(copy1,copy2);
        error->all(FLERR,"Invalid dissociation reaction");
      }
    } else if (r->type == EXCHANGE) {
      if (r->nreactant != 2 || r->nproduct != 2) {
        print_reaction(copy1,copy2);
        error->all(FLERR,"Invalid exchange reaction");
      }
    } else if (r->type == IONIZATION) {
      if (r->nreactant != 2 || (r->nproduct != 2 && r->nproduct != 3)) {
        print_reaction(copy1,copy2);
        error->all(FLERR,"Invalid ionization reaction");
      }
    } else if (r->type == RECOMBINATION) {
      if (r->nreactant != 2 || (r->nproduct != 1 && r->nproduct != 2)) {
        print_reaction(copy1,copy2);
        error->all(FLERR,"Invalid recombination reaction");
      }
    }

    word = strtok(NULL," \t\n\r");
    if (!word) {
      print_reaction(copy1,copy2);
      error->all(FLERR,"Invalid reaction style in file");
    }
    if (word[0] == 'A' || word[0] == 'a') r->style = ARRHENIUS;
    else if (word[0] == 'Q' || word[0] == 'q') r->style = QUANTUM;
    else if (word[0] == 'B' || word[0] == 'b') {
      // 'B' = Arrhenius Backward: TCE reaction whose rate is derived from a
      //   forward reaction via detailed balance (PROTOTYPE, see ReactTCE).
      //   Treated as an Arrhenius reaction; its Arrhenius prefactor (C3) is
      //   overwritten in init() by the paired forward rate and the
      //   temperature-dependent partition-function ratio.
      r->style = ARRHENIUS;
      r->reverse = 1;
    } else {
      print_reaction(copy1,copy2);
      error->all(FLERR,"Invalid reaction style in file");
    }

    if (r->style == ARRHENIUS || r->style == QUANTUM) r->ncoeff = 5;

    for (int i = 0; i < r->ncoeff; i++) {
      word = strtok(NULL," \t\n\r");
      if (!word) {
        print_reaction(copy1,copy2);
        error->all(FLERR,"Invalid reaction coefficients in file");
      }
      r->coeff[i] = input->numeric(FLERR,word);
    }

    word = strtok(NULL," \t\n\r");
    if (word) {
      print_reaction(copy1,copy2);
      error->all(FLERR,"Too many coefficients in a reaction formula");
    }

    nlist++;
  }

  if (comm->me == 0) fclose(fp);
}

/* ----------------------------------------------------------------------
   read one reaction from file
   reaction = 2 lines
   return 1 if end-of-file, else return 0
------------------------------------------------------------------------- */

int ReactBird::readone(char *line1, char *line2, int &n1, int &n2)
{
  char *eof;
  while ((eof = fgets(line1,MAXLINE,fp))) {
    size_t pre = strspn(line1," \t\n\r");
    if (pre == strlen(line1) || line1[pre] == '#') continue;
    eof = fgets(line2,MAXLINE,fp);
    if (!eof) break;
    n1 = strlen(line1) + 1;
    n2 = strlen(line2) + 1;
    return 0;
  }

  return 1;
}

/* ----------------------------------------------------------------------
   check for duplicates in list of reactions read from file
   error if any exist
------------------------------------------------------------------------- */

void ReactBird::check_duplicate()
{
  OneReaction *r,*s;

  for (int i = 0; i < nlist; i++) {
    r = &rlist[i];

    for (int j = i+1; j < nlist; j++) {
      s = &rlist[j];

      if (r->type != s->type) continue;
      if (r->style != s->style) continue;
      if (r->nreactant != s->nreactant) continue;
      if (r->nproduct != s->nproduct) continue;

      int reactant_match = 0;
      if (strcmp(r->id_reactants[0],s->id_reactants[0]) == 0 &&
          strcmp(r->id_reactants[1],s->id_reactants[1]) == 0)
        reactant_match = 1;
      else if (strcmp(r->id_reactants[0],s->id_reactants[1]) == 0 &&
               strcmp(r->id_reactants[1],s->id_reactants[0]) == 0)
        reactant_match = 2;
      if (!reactant_match) continue;

      int product_match = 0;
      if (r->nproduct == 1) {
        if (strcmp(r->id_products[0],s->id_products[0]) == 0)
          product_match = 1;
      } else if (r->nproduct >= 2) {
        if (strcmp(r->id_products[0],s->id_products[0]) == 0 &&
            strcmp(r->id_products[1],s->id_products[1]) == 0)
          product_match = 1;
        else if (strcmp(r->id_products[0],s->id_products[1]) == 0 &&
                 strcmp(r->id_products[1],s->id_products[0]) == 0)
          product_match = 2;
      }
      if (!product_match) continue;

      if (comm->me == 0) {
        printf("MATCH %d %d %d: %d\n",i,j,nlist,product_match);
        printf("MATCH %d %d %d %d\n",
               r->products[0],r->products[1],s->products[0],s->products[1]);
      }
      print_reaction(r);
      print_reaction(s);
      error->all(FLERR,"Duplicate reactions in reaction file");
    }
  }
}

/* ----------------------------------------------------------------------
   print reaction as read from file
   only proc 0 performs output
------------------------------------------------------------------------- */

void ReactBird::print_reaction(char *line1, char *line2)
{
  if (comm->me) return;
  printf("Bad reaction format:\n");
  printf("%s\n%s\n",line1,line2);
};

/* ----------------------------------------------------------------------
   print reaction as stored in rlist
   only proc 0 performs output
------------------------------------------------------------------------- */

void ReactBird::print_reaction(OneReaction *r)
{
  if (comm->me) return;
  printf("Bad reaction:\n");

  char type;
  if (r->type == DISSOCIATION) type = 'D';
  else if (r->type == EXCHANGE) type = 'E';
  else if (r->type == IONIZATION) type = 'I';
  else if (r->type == RECOMBINATION) type = 'R';

  char style;
  if (r->style == ARRHENIUS) style = 'A';
  else if (r->style == QUANTUM) style = 'Q';

  if (r->nproduct == 1)
    printf("  %c %c: %s + %s --> %s\n",type,style,
           r->id_reactants[0],r->id_reactants[1],
           r->id_products[0]);
  else if (r->nproduct == 2)
    printf("  %c %c: %s + %s --> %s %s\n",type,style,
           r->id_reactants[0],r->id_reactants[1],
           r->id_products[0],r->id_products[1]);
  else if (r->nproduct == 3)
    printf("  %c %c: %s + %s --> %s %s %s\n",type,style,
           r->id_reactants[0],r->id_reactants[1],
           r->id_products[0],r->id_products[1],r->id_products[2]);
};

/* ----------------------------------------------------------------------
   print reaction as stored in rlist
   only proc 0 performs output
------------------------------------------------------------------------- */

void ReactBird::print_reaction_ambipolar(OneReaction *r)
{
  if (comm->me) return;
  printf("Bad ambipolar reaction format:\n");
  printf("  type %d style %d\n",r->type,r->style);
  printf("  nreactant %d:",r->nreactant);
  for (int i = 0; i < r->nreactant; i++)
    printf(" %s",r->id_reactants[i]);
  printf("\n");
  printf("  nproduct %d:",r->nproduct);
  for (int i = 0; i < r->nproduct; i++)
    printf(" %s",r->id_products[i]);
  printf("\n");
  printf("  ncoeff %d:",r->ncoeff);
  for (int i = 0; i < r->ncoeff; i++)
    printf(" %g",r->coeff[i]);
  printf("\n");
};

/* ----------------------------------------------------------------------
   return reaction ID = chemical formula
------------------------------------------------------------------------- */

char *ReactBird::reactionID(int m)
{
  return rlist[m].id;
};

/* ----------------------------------------------------------------------
   return tally associated with a reaction
------------------------------------------------------------------------- */

double ReactBird::extract_tally(int m)
{
  if (!tally_flag) {
    tally_flag = 1;
    MPI_Allreduce(tally_reactions,tally_reactions_all,nlist,
                  MPI_SPARTA_BIGINT,MPI_SUM,world);
  }

  return 1.0*tally_reactions_all[m];
};

/* ----------------------------------------------------------------------
   convolve arr (length n, grid spacing du) with one discrete ladder:
   arr_new(u) = sum_levels g * arr_old(u - eps), each level split linearly
   onto the two neighboring grid points; work = scratch of length n
------------------------------------------------------------------------- */

static void ladder_convolve(double *arr, double *work, int n, double du,
                            int nlevel, const double *eps, const double *g)
{
  memcpy(work,arr,n*sizeof(double));
  memset(arr,0,n*sizeof(double));
  for (int m = 0; m < nlevel; m++) {
    if (g[m] == 0.0) continue;
    double sh = eps[m]/du;
    int i0 = (int) sh;
    double frac = sh - i0;
    if (i0 < n) {
      double w0 = g[m]*(1.0-frac);
      for (int k = i0; k < n; k++) arr[k] += w0*work[k-i0];
    }
    if (i0+1 < n) {
      double w1 = g[m]*frac;
      for (int k = i0+1; k < n; k++) arr[k] += w1*work[k-i0-1];
    }
  }
}

/* ----------------------------------------------------------------------
   cap of a forward reaction's numerator, on the energy grid of its reverse
   the forward probability is capped at 1 at run time (ReactTCE::attempt);
     where that cap acts the forward rate falls below its Arrhenius fit,
     and a reverse table built from the uncapped forward numerator would
     run faster than detailed balance with the realized forward rate
   the forward probability at total energy E is
     P_F(E) = pref_F * Gamma_F * num_F(E)/den_F(E)
     so P_F <= 1 is num_F <= den_F/(pref_F*Gamma_F)
   on return cap[k] = den_F(E_k)/(pref_F*Gamma_F) with E_k = k*du - C5_F
     the forward total energy that corresponds to reverse-grid energy
     k*du (energy conservation: the reverse pair holds the forward
     energy plus the forward reaction energy C5_F)
   den_F(E) = E^(zcontF+3/2-omegaF) (x) discrete ladders of the forward
     pair is built on its own grid, which starts at E = 0, and
     interpolated at E_k
------------------------------------------------------------------------- */

void ReactBird::forward_cap(OneReaction *f, double du, int n, double *cap)
{
  Particle::Species *species = particle->species;
  double boltz = update->boltz;

  int sps[2];
  sps[0] = f->reactants[0];
  sps[1] = f->reactants[1];

  int vibdiscrete = (collide->vibstyle == DISCRETE);
  int elecdiscrete = (collide->elecstyle == DISCRETE);

  double zcont = 0.5*(species[sps[0]].rotdof + species[sps[1]].rotdof);
  if (collide->vibstyle == SMOOTH)
    zcont += 0.5*(species[sps[0]].vibdof + species[sps[1]].vibdof);
  double omega = collide->extract(sps[0],sps[1],"omega");

  double gden = tgamma(zcont + f->coeff[3] + 1.5);
  if (gden < 1.0e-6) gden = 1.0e-6;
  double pg = f->coeff[2] * tgamma(zcont + 2.5 - omega) / gden;

  if (!(pg > 0.0)) {
    for (int k = 0; k < n; k++) cap[k] = 1.0e300;
    return;
  }

  double shift = -f->coeff[4]/du;
  int ne = n + (int) fabs(shift) + 3;
  double emax = ne*du;

  double *den,*work,*leps,*lg;
  memory->create(den,ne,"react:cap_den");
  memory->create(work,ne,"react:cap_work");

  double expden = zcont + 1.5 - omega;
  den[0] = 0.0;
  for (int j = 1; j < ne; j++) den[j] = pow(j*du,expden);

  double theta_min = 0.0;
  if (vibdiscrete)
    for (int s = 0; s < 2; s++)
      for (int m = 0; m < species[sps[s]].nvibmode; m++) {
        double th = species[sps[s]].vibtemp[m];
        if (th > 0.0 && (theta_min == 0.0 || th < theta_min)) theta_min = th;
      }
  int maxlev = (theta_min > 0.0) ? (int) (emax/(boltz*theta_min)) + 2 : 0;
  if (max_nelecstate() + 2 > maxlev) maxlev = max_nelecstate() + 2;
  memory->create(leps,maxlev,"react:cap_leps");
  memory->create(lg,maxlev,"react:cap_lg");

  for (int s = 0; s < 2; s++) {
    int sp = sps[s];
    if (vibdiscrete)
      for (int m = 0; m < species[sp].nvibmode; m++) {
        double th = species[sp].vibtemp[m];
        if (th <= 0.0) continue;
        int d = species[sp].vibdegen[m] > 1 ? species[sp].vibdegen[m] : 1;
        int nlev = 0;
        double g = 1.0;
        for (int l = 0; l*th*boltz < emax && nlev < maxlev; l++) {
          leps[nlev] = l*th*boltz;
          if (l > 0) g = g*(l+d-1)/l;
          lg[nlev] = g;
          nlev++;
        }
        ladder_convolve(den,work,ne,du,nlev,leps,lg);
      }
    if (elecdiscrete && species[sp].elecdat) {
      int nlev = species[sp].elecdat->nelecstate;
      for (int l = 0; l < nlev; l++) {
        leps[l] = boltz*species[sp].elecdat->states[l].temp;
        lg[l] = species[sp].elecdat->states[l].degen;
      }
      ladder_convolve(den,work,ne,du,nlev,leps,lg);
    }
  }

  for (int k = 0; k < n; k++) {
    double x = k + shift;
    if (x <= 0.0) cap[k] = 0.0;
    else if (x >= ne-1) cap[k] = den[ne-1]/pg;
    else {
      int j = (int) x;
      double fr = x - j;
      cap[k] = ((1.0-fr)*den[j] + fr*den[j+1])/pg;
    }
  }

  memory->destroy(den);
  memory->destroy(work);
  memory->destroy(leps);
  memory->destroy(lg);
}

/* ----------------------------------------------------------------------
   build the per-reaction tables of the microcanonical TCE energy factor
   used with react_modify partial_energy no:
     factor(E_tot) = sum_p g_p (E_tot - eps_p - Ea)_+^(zrot+eta+1/2)
                   / sum_p g_p (E_tot - eps_p)_+^(zrot+3/2-omega)
   where p runs over the joint discrete states of the two reactants: the
   SHO levels of every vibrational mode (level degeneracy of a d-fold
   degenerate mode = C(l+d-1,d-1)) when vibration is discrete, and the
   electronic states when electronic energy is discrete. The continuum is translation + rotation only
   (Gamma-distributed at equilibrium), so this factor keeps the
   equilibrium TCE rate on the input Arrhenius rate with DISCRETE
   vibration, in place of the instantaneous-vibrational-DOF heuristic.
   The numerator/denominator sums are accumulated on an energy grid by
   convolution (shift-and-add per ladder level); runtime evaluation is a
   single linear interpolation (vib_micro_factor() in react_bird.h).
   Tables are clamped at umax = Ea + 60 eV; collisions beyond use the
   last table value (their probability weight is negligible below
   ~100000 K, and clamping keeps CPU/Kokkos evaluation identical).
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   largest electronic-state count over all species, used to size the
   per-ladder scratch buffers: the electronic ladder of any one species
   can be longer than the vibrational-grid depth that otherwise bounds
   the scratch arrays, so the buffers must cover it to avoid a heap
   overflow when a species carries a densely resolved elecfile
------------------------------------------------------------------------- */

int ReactBird::max_nelecstate()
{
  Particle::Species *species = particle->species;
  int maxn = 0;
  for (int sp = 0; sp < particle->nspecies; sp++)
    if (species[sp].elecdat && species[sp].elecdat->nelecstate > maxn)
      maxn = species[sp].elecdat->nelecstate;
  return maxn;
}

/* ----------------------------------------------------------------------
   warn, on the first run only, about inputs that the detailed-balance
     reverse tables silently mis-handle:
   (a) a rotating species in an active reverse whose rotfile line omits
       the symmetry number: sigma = 1 is assumed, a constant factor-sigma
       error in Keq per such species (2 for a homonuclear molecule) that
       the tables' self-consistency checks cannot detect
   (b) collide_modify rotate no: the tables count rotational degrees of
       freedom, but particles then carry no rotational energy, so reverse
       rates involving rotating species are biased, and a reverse
       recombination with a nonlinear third body never fires (its
       continuum weight erot3^(1/2) is zero)
   (c) a reaction whose products and reactants differ in mass: every
       event then changes the centre-of-mass kinetic energy, and a gas
       cycling through forward and reverse reactions gains or loses
       energy steadily
   (d) reaction energies that do not sum to zero around reaction cycles
       (e.g. N2 + O -> NO + N, NO + M -> N + O + M, N + N + M -> N2 + M):
       each loop then creates or destroys energy, and the equilibrium
       constants admit no common equilibrium.  Tested by fitting one
       formation energy per species to the forward reaction energies
       (least squares) and reporting the largest residual
------------------------------------------------------------------------- */

void ReactBird::check_reverse_modes()
{
  if (reverse_modes_checked) return;
  reverse_modes_checked = 1;
  if (comm->me) return;

  Particle::Species *species = particle->species;
  int nspecies = particle->nspecies;
  int *flag = new int[nspecies];
  for (int i = 0; i < nspecies; i++) flag[i] = 0;

  int nreverse = 0;
  int nonlinear3 = 0;

  for (int m = 0; m < nlist; m++) {
    OneReaction *r = &rlist[m];
    if (!r->active || !r->reverse || r->reverse_partner < 0) continue;
    nreverse++;
    for (int i = 0; i < r->nreactant; i++)
      if (r->reactants[i] >= 0) flag[r->reactants[i]] = 1;
    for (int i = 0; i < r->nproduct; i++)
      if (r->products[i] >= 0) flag[r->products[i]] = 1;
    if (r->type == RECOMBINATION && r->products[1] >= 0 &&
        species[r->products[1]].rotdof == 3) nonlinear3 = 1;
  }

  if (nreverse) {
    std::string missing;
    int rotating = 0;
    for (int i = 0; i < nspecies; i++) {
      if (!flag[i] || species[i].rotdof == 0) continue;
      rotating = 1;
      if (species[i].nrottemp > 0 && !species[i].rotsymm_read) {
        if (!missing.empty()) missing += " ";
        missing += species[i].id;
      }
    }

    if (!missing.empty()) {
      std::string mesg = "Rotation file gives no symmetry number for "
        "species used by detailed-balance reverse reactions, so sigma = 1 "
        "is assumed (use 2 for homonuclear molecules): " + missing;
      error->warning(FLERR,mesg.c_str());
    }

    if (rotating && collide && collide->rotstyle == NONE) {
      std::string mesg = "Detailed-balance reverse reactions with "
        "collide_modify rotate no: particles carry no rotational energy, "
        "so reverse rates involving rotating species are biased";
      if (nonlinear3) mesg += ", and recombination with a nonlinear "
                        "third body never occurs";
      error->warning(FLERR,mesg.c_str());
    }

    // (c) mass conservation of every active reaction
    // the third body of a recombination is not listed as a reactant

    int nmass = 0, mworst = -1;
    double dmworst = 0.0;
    for (int m = 0; m < nlist; m++) {
      OneReaction *r = &rlist[m];
      if (!r->active) continue;
      double mr = 0.0, mp = 0.0;
      for (int i = 0; i < r->nreactant; i++)
        if (r->reactants[i] >= 0) mr += species[r->reactants[i]].mass;
      int nprod = (r->type == RECOMBINATION) ? 1 : r->nproduct;
      for (int i = 0; i < nprod; i++)
        if (r->products[i] >= 0) mp += species[r->products[i]].mass;
      double dm = fabs(mp-mr) / mr;
      if (dm > 1.0e-6) {
        nmass++;
        if (dm > dmworst) { dmworst = dm; mworst = m; }
      }
    }
    if (nmass) {
      char str[512];
      snprintf(str,sizeof(str),"%d reaction(s) do not conserve mass "
               "(largest relative defect %.2g, reaction %s): with "
               "detailed-balance reverse reactions energy is then not "
               "conserved; make the species masses consistent",
               nmass,dmworst,rlist[mworst].id);
      error->warning(FLERR,str);
    }

    // (d) closure of the reaction energies around reaction cycles
    // unknowns: one formation energy per species; equations: for each
    //   active forward reaction, sum(E_reactants) - sum(E_products) = C5
    //   (a reverse reaction carries -C5 of its partner by construction)
    // normal equations solved by Gaussian elimination with pivoting; the
    //   null space (one free constant per conserved element and charge)
    //   is fixed by zeroing the variables with vanishing pivots

    int n = nspecies;
    double **a = new double*[n];
    double *rhs = new double[n];
    double *ef = new double[n];
    for (int i = 0; i < n; i++) {
      a[i] = new double[n];
      for (int j = 0; j < n; j++) a[i][j] = 0.0;
      rhs[i] = 0.0;
      ef[i] = 0.0;
    }
    double *nu = new double[n];
    double c5max = 0.0;

    for (int m = 0; m < nlist; m++) {
      OneReaction *r = &rlist[m];
      if (!r->active || r->reverse) continue;
      for (int i = 0; i < n; i++) nu[i] = 0.0;
      for (int i = 0; i < r->nreactant; i++)
        if (r->reactants[i] >= 0) nu[r->reactants[i]] += 1.0;
      int nprod = (r->type == RECOMBINATION) ? 1 : r->nproduct;
      for (int i = 0; i < nprod; i++)
        if (r->products[i] >= 0) nu[r->products[i]] -= 1.0;
      for (int i = 0; i < n; i++) {
        if (nu[i] == 0.0) continue;
        rhs[i] += nu[i]*r->coeff[4];
        for (int j = 0; j < n; j++) a[i][j] += nu[i]*nu[j];
      }
      if (fabs(r->coeff[4]) > c5max) c5max = fabs(r->coeff[4]);
    }

    double amax = 0.0;
    for (int i = 0; i < n; i++) if (a[i][i] > amax) amax = a[i][i];
    int *free = new int[n];
    for (int i = 0; i < n; i++) free[i] = 0;

    for (int k = 0; k < n && amax > 0.0; k++) {
      int piv = k;
      for (int i = k+1; i < n; i++)
        if (fabs(a[i][k]) > fabs(a[piv][k])) piv = i;
      if (fabs(a[piv][k]) < 1.0e-9*amax) { free[k] = 1; continue; }
      if (piv != k) {
        double *t = a[piv]; a[piv] = a[k]; a[k] = t;
        double tr = rhs[piv]; rhs[piv] = rhs[k]; rhs[k] = tr;
      }
      for (int i = k+1; i < n; i++) {
        double f = a[i][k]/a[k][k];
        if (f == 0.0) continue;
        for (int j = k; j < n; j++) a[i][j] -= f*a[k][j];
        rhs[i] -= f*rhs[k];
      }
    }
    for (int k = n-1; k >= 0 && amax > 0.0; k--) {
      if (free[k]) { ef[k] = 0.0; continue; }
      double sum = rhs[k];
      for (int j = k+1; j < n; j++) sum -= a[k][j]*ef[j];
      ef[k] = sum/a[k][k];
    }

    double resmax = 0.0;
    int rworst = -1;
    for (int m = 0; m < nlist; m++) {
      OneReaction *r = &rlist[m];
      if (!r->active || r->reverse) continue;
      double sum = 0.0;
      for (int i = 0; i < r->nreactant; i++)
        if (r->reactants[i] >= 0) sum += ef[r->reactants[i]];
      int nprod = (r->type == RECOMBINATION) ? 1 : r->nproduct;
      for (int i = 0; i < nprod; i++)
        if (r->products[i] >= 0) sum -= ef[r->products[i]];
      double res = fabs(sum - r->coeff[4]);
      if (res > resmax) { resmax = res; rworst = m; }
    }
    if (c5max > 0.0 && resmax > 1.0e-5*c5max) {
      char str[512];
      snprintf(str,sizeof(str),"Reaction energies are not consistent "
               "around reaction cycles (largest residual %.3g, %.2g of "
               "the largest reaction energy, at reaction %s): with "
               "detailed-balance reverse reactions energy is then not "
               "conserved; derive all reaction energies from one set "
               "of species formation energies",
               resmax,resmax/c5max,rlist[rworst].id);
      error->warning(FLERR,str);
    }

    for (int i = 0; i < n; i++) delete [] a[i];
    delete [] a;
    delete [] rhs;
    delete [] ef;
    delete [] nu;
    delete [] free;
  }

  delete [] flag;
}

/* ---------------------------------------------------------------------- */

void ReactBird::build_micro_tables()
{
  free_micro_tables();
  check_reverse_modes();

  Particle::Species *species = particle->species;
  double boltz = update->boltz;
  // 60 eV in the energy unit of the unit style (J for si, erg for cgs)
  double ev60 = 60.0 * (strcmp(update->unit_style,"cgs") == 0 ?
                        1.602176634e-12 : 1.602176634e-19);

  mtab_nlist = nlist;
  mtab = new double*[nlist];
  mtab_num = new double*[nlist];
  memory->create(mtab_du,nlist,"react:mtab_du");
  memory->create(mtab_n,nlist,"react:mtab_n");

  int sps[2];
  double *num,*den,*work,*leps,*lg;

  for (int i = 0; i < nlist; i++) {
    mtab[i] = NULL;
    mtab_num[i] = NULL;
    mtab_du[i] = 0.0;
    mtab_n[i] = 0;

    OneReaction *r = &rlist[i];
    if (!r->active || r->nreactant != 2) continue;

    // a reverse (B-style) exchange reaction gets a temperature-free
    // detailed-balance table (the energy-resolved microscopic-reversibility
    // shape) instead of the standard energy factor.  With an external Keq
    // fit it keeps the same table -- calibrated to the statistical-mechanics
    // Keq -- and applies only the residual thermal correction
    // R(T) = Keq_statmech(T)/Keq_ext(T) at the cell temperature (see
    // assign_keq_fits / ReactTCE::attempt), so the collision-energy
    // selectivity stays microscopically reversible and only the (usually
    // small) Keq discrepancy is thermal

    if (r->reverse && r->type == EXCHANGE && r->reverse_partner >= 0) {
      build_db_table(i);
      continue;
    }

    // a reverse (B-style) recombination gets 3-body detailed-balance
    // tables: the pair density of states in its mtab slot plus the
    // calibrated forward-numerator/flat-measure ratio in mtab_num; an
    // external Keq fit uses the same tables times the residual R(T)

    if (r->reverse && r->type == RECOMBINATION && r->reverse_partner >= 0) {
      build_db3_table(i);
      continue;
    }

    sps[0] = r->reactants[0];
    sps[1] = r->reactants[1];

    // grid resolution from the smallest vibrational quantum;
    // skip the table entirely if the reactants carry no discrete ladders
    // (runtime then falls back to the standard / elec-micro factor)
    //
    // rotation stays in the continuum even with collide_modify rotate
    // discrete, for two reasons: (a) the rigid-rotor ladder is so dense
    // (theta_rot of order 1 K) that its density of states differs from
    // the continuum limit only to O(theta_rot/T); and (b) removing
    // rotation from the continuum makes the TCE calibration singular for
    // temperature exponents eta <= zrot - 3/2 (the numerator seed
    // exponent zcont+eta+1/2 reaches -1 and the required energy factor
    // degenerates to a delta comb), which excludes common rate sets
    // (e.g. eta = -1.5 dissociation): Bird's TCE validity constraint
    // zbar + eta + 3/2 > 0 is satisfied by keeping rotation continuous

    int vibdiscrete = (collide->vibstyle == DISCRETE);
    int elecdiscrete = (collide->elecstyle == DISCRETE);

    double theta_min = 0.0;
    int nladder = 0;
    for (int s = 0; s < 2; s++) {
      int sp = sps[s];
      if (vibdiscrete)
        for (int m = 0; m < species[sp].nvibmode; m++) {
          nladder++;
          double th = species[sp].vibtemp[m];
          if (theta_min == 0.0 || th < theta_min) theta_min = th;
        }
      if (elecdiscrete && species[sp].elecdat) nladder++;
    }
    if (nladder == 0) continue;

    // threshold of the numerator: the signed activation energy, as in the
    // analytic factor tce_efactor() (a negative Ea shifts the numerator,
    // (u - Ea)^a = (u + |Ea|)^a, rather than removing the barrier; the
    // derived-reverse tables seed their numerators unclamped likewise);
    // the grid extent uses Ea clamped at 0

    double ea = r->coeff[1];
    double eaP = ea > 0.0 ? ea : 0.0;
    double umax = eaP + ev60;   // + 60 eV: Boltzmann tail to ~1e5 K
    double du;
    if (theta_min > 0.0) du = boltz*theta_min/16.0;
    else du = umax/20000.0;
    if (umax/du > 200000.0) du = umax/200000.0;
    int n = (int) (umax/du) + 2;

    mtab_du[i] = du;
    mtab_n[i] = n;

    memory->create(num,n,"react:mtab_num");
    memory->create(den,n,"react:mtab_den");
    memory->create(work,n,"react:mtab_work");

    // continuum internal DOF: rotation, plus vibration when it is a
    // continuous (smooth) mode; must match the z used at runtime in
    // ReactTCE::attempt so the Gamma prefactor stays consistent with
    // the table seeds

    double zcont = 0.5*(species[sps[0]].rotdof + species[sps[1]].rotdof);
    if (collide->vibstyle == SMOOTH)
      zcont += 0.5*(species[sps[0]].vibdof + species[sps[1]].vibdof);
    double omega = collide->extract(sps[0],sps[1],"omega");
    double exp_num = zcont + r->coeff[3] + 0.5;
    double exp_den = zcont + 1.5 - omega;

    for (int k = 0; k < n; k++) {
      double u = k*du;
      den[k] = (u > 0.0) ? pow(u,exp_den) : 0.0;
      num[k] = (u > ea) ? pow(u-ea,exp_num) : 0.0;
    }

    // convolve numerator and denominator with each ladder

    // with no discrete vibrational ladder (theta_min = 0) only the
    // electronic ladders below use leps/lg
    int maxlev = (theta_min > 0.0) ? (int) (umax/(boltz*theta_min)) + 2 : 0;
    if (max_nelecstate() + 2 > maxlev) maxlev = max_nelecstate() + 2;
    memory->create(leps,maxlev,"react:mtab_leps");
    memory->create(lg,maxlev,"react:mtab_lg");

    for (int s = 0; s < 2; s++) {
      int sp = sps[s];

      if (vibdiscrete)
        for (int m = 0; m < species[sp].nvibmode; m++) {
          double th = species[sp].vibtemp[m];
          int d = species[sp].vibdegen[m] > 1 ? species[sp].vibdegen[m] : 1;
          int nlev = 0;
          double g = 1.0;
          for (int l = 0; l*th*boltz < umax && nlev < maxlev; l++) {
            leps[nlev] = l*th*boltz;
            // level degeneracy of a d-fold degenerate SHO mode:
            // C(l+d-1,d-1), built iteratively
            if (l > 0) g = g*(l+d-1)/l;
            lg[nlev] = g;
            nlev++;
          }
          ladder_convolve(num,work,n,du,nlev,leps,lg);
          ladder_convolve(den,work,n,du,nlev,leps,lg);
        }

      if (elecdiscrete && species[sp].elecdat) {
        int nlev = species[sp].elecdat->nelecstate;
        for (int l = 0; l < nlev; l++) {
          leps[l] = boltz*species[sp].elecdat->states[l].temp;
          lg[l] = species[sp].elecdat->states[l].degen;
        }
        ladder_convolve(num,work,n,du,nlev,leps,lg);
        ladder_convolve(den,work,n,du,nlev,leps,lg);
      }
    }

    mtab[i] = new double[n];
    for (int k = 0; k < n; k++)
      mtab[i][k] = (den[k] > 0.0) ? num[k]/den[k] : 0.0;

    memory->destroy(leps);
    memory->destroy(lg);
    memory->destroy(num);
    memory->destroy(den);
    memory->destroy(work);
  }

  // second pass over the reverse tables: when the capped probabilities
  // of all channels of a forward pair sum to more than 1, the run-time
  // selection scales each of them by 1/sum (ReactTCE::attempt), so the
  // forward reaction realizes min(P,1)/sum.  Apply the same factor to
  // its reverse table, at the forward pair's total energy, so that the
  // pair stays reversible under saturation.  The sum uses the tables as
  // built above (a reverse channel in the sum enters with its own
  // single-channel cap, which is exact to first order).  All factors are
  // computed before any table is changed: a reverse exchange can itself
  // be a channel of another reverse's forward pair, and dividing its
  // table in the same loop would make the result depend on the order of
  // the reactions in the input file.

  double **satfac = new double*[nlist];      // pair sums S_F > 1, or 1
  for (int i = 0; i < nlist; i++) satfac[i] = NULL;

  for (int i = 0; i < nlist; i++) {
    OneReaction *r = &rlist[i];
    if (!r->active || !r->reverse || r->reverse_partner < 0) continue;
    if (!mtab[i]) continue;
    double *tab = (r->type == RECOMBINATION) ? mtab_num[i] : mtab[i];
    if (!tab) continue;
    OneReaction *f = &rlist[r->reverse_partner];
    int isp = f->reactants[0];
    int jsp = f->reactants[1];
    int nch = reactions[isp][jsp].n;
    int *list = reactions[isp][jsp].list;
    if (nch < 2) continue;
    // per-channel constants of pair_channel_prob(), hoisted out of the
    // grid loop (same values, same operation order)
    double *cz = new double[nch];
    double *cpre = new double[nch];
    int *clive = new int[nch];
    for (int c = 0; c < nch; c++) {
      OneReaction *rc = &rlist[list[c]];
      clive[c] = (rc->active && rc->type != RECOMBINATION);
      cz[c] = clive[c] ? pair_channel_z(rc) : 0.0;
      cpre[c] = clive[c] ? tce_prefactor(rc,cz[c],rc->coeff[2]) : 0.0;
    }
    satfac[i] = new double[mtab_n[i]];
    for (int k = 0; k < mtab_n[i]; k++) {
      satfac[i][k] = 1.0;
      if (tab[k] <= 0.0) continue;
      double e = k*mtab_du[i] - f->coeff[4];
      if (e <= 0.0) continue;
      double sum = 0.0;
      for (int c = 0; c < nch; c++) {
        if (!clive[c] || tce_excess(&rlist[list[c]],e) <= 0.0) continue;
        double p = cpre[c] * tce_efactor(list[c],e,cz[c],1);
        if (p > 1.0) p = 1.0;
        sum += p;
      }
      if (sum > 1.0) satfac[i][k] = sum;
    }
    delete [] cz;
    delete [] cpre;
    delete [] clive;
  }

  for (int i = 0; i < nlist; i++) {
    if (!satfac[i]) continue;
    double *tab = (rlist[i].type == RECOMBINATION) ? mtab_num[i] : mtab[i];
    for (int k = 0; k < mtab_n[i]; k++)
      if (satfac[i][k] > 1.0) tab[k] /= satfac[i][k];
    delete [] satfac[i];
  }
  delete [] satfac;
}

/* ----------------------------------------------------------------------
   probability of two-body reaction j as a function of the total collision
   energy of its reactant pair alone (total-energy model), uncapped;
   mirrors ReactTCE::channel_prob() for dissociation, ionization and
   exchange, and returns 0 for recombination, whose probability also
   depends on the third particle
------------------------------------------------------------------------- */

double ReactBird::pair_channel_prob(int j, double ecc)
{
  OneReaction *r = &rlist[j];
  if (!r->active || r->type == RECOMBINATION) return 0.0;
  double z = pair_channel_z(r);
  if (tce_excess(r,ecc) <= 0.0) return 0.0;
  return tce_prefactor(r,z,r->coeff[2]) * tce_efactor(j,ecc,z,1);
}

/* ----------------------------------------------------------------------
   continuum z of a two-body channel in the total-energy model, as
   ReactTCE::channel_prob() forms it at run time
------------------------------------------------------------------------- */

double ReactBird::pair_channel_z(OneReaction *r)
{
  Particle::Species *species = particle->species;
  int isp = r->reactants[0];
  int jsp = r->reactants[1];
  double z = (species[isp].rotdof + species[jsp].rotdof)/2.0;
  if (collide->vibstyle == SMOOTH)
    z += (species[isp].vibdof + species[jsp].vibdof)/2.0;
  return z;
}

/* ----------------------------------------------------------------------
   total partition function per unit volume for a species at temperature T,
     q = q_trans * q_rot * q_vib * q_elec   (issue #472 reverse reactions)
   - translational: (2 pi m kB T / h^2)^{3/2}  (per unit volume)
   - rotational:    rigid rotor, linear molecule, T/(sigma theta_r), with
                    the symmetry number sigma from the species rotfile
                    (nonlinear molecules use the classical
                    sqrt(pi/(sigma^2) * T^3/(tA tB tC)) form)
   - vibrational:   harmonic oscillator, ground-state referenced, with
                    d-fold degenerate modes contributing (1-x)^-d
   - electronic:    sum over the species elecfile ladder, g_i e^(-theta_i/T);
                    1 if the species defines no electronic data
   energies are referenced to each species' own ground state, consistent
   with the reaction energy (C5) used to seed the backward coefficients
------------------------------------------------------------------------- */

double ReactBird::partition_function(int isp, double T)
{
  Particle::Species *sp = &particle->species[isp];
  const double kb = update->boltz;
  const int cgs = (strcmp(update->unit_style,"cgs") == 0);
  const double h = cgs ? 6.62607015e-27 : 6.62607015e-34;  // erg s or J s

  // translational partition function per unit volume

  double qtrans = pow(2.0*MY_PI*sp->mass*kb*T/(h*h), 1.5);

  // rotational partition function (rigid rotor), high-temperature form:
  // exact to O(theta_r/T), i.e. to ~1e-4 at any temperature where the
  // reaction rates themselves are non-negligible

  double qrot = 1.0;
  if (sp->rotdof == 2 && sp->nrottemp >= 1 && sp->rottemp[0] > 0.0)
    qrot = T / (sp->rotsymm * sp->rottemp[0]);
  else if (sp->rotdof == 3 && sp->nrottemp == 3 &&
           sp->rottemp[0] > 0.0 && sp->rottemp[1] > 0.0 &&
           sp->rottemp[2] > 0.0)
    qrot = sqrt(MY_PI*T*T*T /
                (sp->rotsymm*sp->rotsymm *
                 sp->rottemp[0]*sp->rottemp[1]*sp->rottemp[2]));

  // vibrational partition function: quantum harmonic oscillator
  // (ground-state referenced) with DISCRETE vibration, matching the
  // energy-resolved SHO ladders convolved into the detailed-balance
  // tables; classical harmonic (q = T/theta per mode) with SMOOTH
  // vibration, matching the continuum vibrational DOF the tables and the
  // forward TCE rate use in that mode -- so the calibration target and
  // the table DOS carry the same vibrational temperature dependence

  // the partition functions must contain only the internal states the
  //   simulation models, else the calibration target and the table
  //   densities of states differ by a T-dependent factor (the table drift
  //   check then fires, or the reverse rate is biased); a mode the
  //   collision model does not carry (collide_modify vibrate no or
  //   electronic no) keeps only its ground-state contribution: 1 for a
  //   ground-referenced oscillator, the ground-state degeneracy for an
  //   electronic ladder -- T-independent constants that belong in Keq and
  //   cannot cause drift

  double qvib = 1.0;
  int vibsmooth = (collide && collide->vibstyle == SMOOTH);
  int vibnone = (collide && collide->vibstyle == NONE);
  for (int m = 0; m < sp->nvibmode && !vibnone; m++) {
    if (sp->vibtemp[m] <= 0.0) continue;
    int g = sp->vibdegen[m] > 0 ? sp->vibdegen[m] : 1;
    if (vibsmooth)
      qvib *= pow(T/sp->vibtemp[m], g);
    else {
      double x = exp(-sp->vibtemp[m]/T);
      qvib *= pow(1.0/(1.0-x), g);
    }
  }

  // electronic partition function from the species elecfile ladder

  double qelec = 1.0;
  int elecdiscrete = (collide && collide->elecstyle == DISCRETE);
  if (sp->elecdat && !elecdiscrete) qelec = sp->elecdat->states[0].degen;
  else if (sp->elecdat) {
    qelec = 0.0;
    for (int i = 0; i < sp->elecdat->nelecstate; i++)
      qelec += sp->elecdat->states[i].degen *
        exp(-sp->elecdat->states[i].temp/T);
  }

  // a free electron carries spin degeneracy 2: it is structureless (no
  // elecfile ladder), so without this factor the Saha-type equilibrium
  // constant of an electron-producing reaction (associative ionization
  // vs dissociative recombination) would be off by exactly 2x.  The
  // electron is identified by its physics (negative charge and a mass
  // far below any atom -- the electron is ~9.1e-31 kg, ~1800x lighter
  // than the lightest nucleus; the threshold follows the unit style)
  // rather than by fix ambipolar's especies: the reverse machinery also
  // serves reactions with real electron particles and no fix ambipolar,
  // so no especies is available here.

  else if (sp->charge < 0.0 && sp->mass < (cgs ? 1.0e-26 : 1.0e-29))
    qelec = 2.0;

  return qtrans*qrot*qvib*qelec;
}

/* ----------------------------------------------------------------------
   build the temperature-free detailed-balance table of a reverse
   (B-style) EXCHANGE reaction, stored in its mtab slot:
     P_b(u) = prefactor * Gamma-ratio * table(u)
     table(u) = scale * num_F(u) / den_B(u)
   where u is the total collision energy of the backward reactant pair,
     den_B(u) = u^(zcontB+3/2-omegaB) (x) discrete ladders of the pair
       = collision-weighted density of states of the colliding pair
       (identical to the standard table denominator), and
     num_F(u) = (u - ea_eff)_+^(zcontF+etaF+1/2) (x) discrete ladders of
       the FORWARD reactant pair, with ea_eff = Ea_F + dHf the backward
       threshold: num_F is the forward reaction's microcanonical
       numerator expressed in the backward channel's energy variable
       (energy conservation shifts the argument, which folds into the
       threshold), so the ratio enforces energy-resolved microscopic
       reversibility between the two channels.
   the constant `scale` collects all channel constants (VSS collision
   rates, reduced masses, rotational 1/(sigma k theta) factors, Gamma
   normalizations): it is fixed by calibrating the thermal average of
   P_b against the exact detailed-balance target
     k_b(T) = A_F T^bF e^(-ea_eff/kT) q_prod(T)/q_react(T)
   at one temperature; the ratio of the two sides is temperature
   independent by construction (both are Laplace transforms of the same
   energy profile), so one calibration temperature suffices, and doing
   it on the discrete grid absorbs quadrature bias as well.
   the backward rate then satisfies k_b(T) = k_f(T)/K_eq(T) at every
   temperature with no temperature evaluated at run time.
------------------------------------------------------------------------- */

void ReactBird::build_db_table(int i)
{
  Particle::Species *species = particle->species;
  double boltz = update->boltz;
  // 60 eV in the energy unit of the unit style (J for si, erg for cgs)
  double ev60 = 60.0 * (strcmp(update->unit_style,"cgs") == 0 ?
                        1.602176634e-12 : 1.602176634e-19);

  OneReaction *r = &rlist[i];
  OneReaction *f = &rlist[r->reverse_partner];

  int spsB[2],spsF[2];
  spsB[0] = r->reactants[0]; spsB[1] = r->reactants[1];
  spsF[0] = f->reactants[0]; spsF[1] = f->reactants[1];

  int vibdiscrete = (collide->vibstyle == DISCRETE);
  int elecdiscrete = (collide->elecstyle == DISCRETE);

  // backward threshold: forward barrier shifted into the backward
  // channel by the reaction energy (unclamped: a negative value means
  // the backward reaction is barrierless with excess energy)

  double ea_eff = f->coeff[1] + f->coeff[4];
  double eaP = ea_eff > 0.0 ? ea_eff : 0.0;

  // energy grid: resolve the smallest vibrational quantum of EITHER
  // pair; always build (a ladder-free pair still needs the continuum
  // density-of-states ratio)

  double theta_min = 0.0;
  for (int s = 0; s < 4; s++) {
    int sp = s < 2 ? spsB[s] : spsF[s-2];
    if (vibdiscrete)
      for (int m = 0; m < species[sp].nvibmode; m++) {
        double th = species[sp].vibtemp[m];
        if (theta_min == 0.0 || th < theta_min) theta_min = th;
      }
  }

  double umax = eaP + ev60;   // + 60 eV: Boltzmann tail to ~1e5 K
  double du;
  if (theta_min > 0.0) du = boltz*theta_min/16.0;
  else du = umax/20000.0;
  if (umax/du > 200000.0) du = umax/200000.0;
  int n = (int) (umax/du) + 2;

  mtab_du[i] = du;
  mtab_n[i] = n;

  double *num,*den,*work,*leps,*lg;
  memory->create(num,n,"react:mtab_num");
  memory->create(den,n,"react:mtab_den");
  memory->create(work,n,"react:mtab_work");

  // continuum internal DOF and omega of each channel, matching the
  // runtime z of the colliding (backward) pair and the forward table
  // construction respectively

  double zcontB = 0.5*(species[spsB[0]].rotdof + species[spsB[1]].rotdof);
  double zcontF = 0.5*(species[spsF[0]].rotdof + species[spsF[1]].rotdof);
  if (collide->vibstyle == SMOOTH) {
    zcontB += 0.5*(species[spsB[0]].vibdof + species[spsB[1]].vibdof);
    zcontF += 0.5*(species[spsF[0]].vibdof + species[spsF[1]].vibdof);
  }
  double omegaB = collide->extract(spsB[0],spsB[1],"omega");

  double exp_num = zcontF + f->coeff[3] + 0.5;
  double exp_den = zcontB + 1.5 - omegaB;

  for (int k = 0; k < n; k++) {
    double u = k*du;
    den[k] = (u > 0.0) ? pow(u,exp_den) : 0.0;
    num[k] = (u > ea_eff) ? pow(u-ea_eff,exp_num) : 0.0;
  }

  // convolve den with the backward pair's ladders and num with the
  // forward pair's ladders

  // with no discrete vibrational ladder (theta_min = 0) only the
  // electronic ladders below use leps/lg
  int maxlev = (theta_min > 0.0) ? (int) (umax/(boltz*theta_min)) + 2 : 0;
  if (max_nelecstate() + 2 > maxlev) maxlev = max_nelecstate() + 2;
  memory->create(leps,maxlev,"react:mtab_leps");
  memory->create(lg,maxlev,"react:mtab_lg");

  for (int which = 0; which < 2; which++) {
    double *arr = which ? num : den;
    int *sps = which ? spsF : spsB;

    for (int s = 0; s < 2; s++) {
      int sp = sps[s];

      if (vibdiscrete)
        for (int m = 0; m < species[sp].nvibmode; m++) {
          double th = species[sp].vibtemp[m];
          int d = species[sp].vibdegen[m] > 1 ? species[sp].vibdegen[m] : 1;
          int nlev = 0;
          double g = 1.0;
          for (int l = 0; l*th*boltz < umax && nlev < maxlev; l++) {
            leps[nlev] = l*th*boltz;
            if (l > 0) g = g*(l+d-1)/l;
            lg[nlev] = g;
            nlev++;
          }
          ladder_convolve(arr,work,n,du,nlev,leps,lg);
        }

      if (elecdiscrete && species[sp].elecdat) {
        int nlev = species[sp].elecdat->nelecstate;
        for (int l = 0; l < nlev; l++) {
          leps[l] = boltz*species[sp].elecdat->states[l].temp;
          lg[l] = species[sp].elecdat->states[l].degen;
        }
        ladder_convolve(arr,work,n,du,nlev,leps,lg);
      }
    }
  }

  mtab[i] = new double[n];
  for (int k = 0; k < n; k++)
    mtab[i][k] = (den[k] > 0.0) ? num[k]/den[k] : 0.0;

  // calibrate the channel constant: match the thermal average of the
  // runtime probability, times the VSS collision rate of the backward
  // pair (in SPARTA's own convention, inverted from the C1 transform:
  // kcoll = 2 sqrt(pi) d^2 sqrt(2 k Tref/mr) (T/Tref)^(1-omega) / eps),
  // to the exact detailed-balance rate at one temperature; the common
  // factor e^(-eaP/kT) is removed from both sides so large thresholds
  // cannot underflow

  double tcal = 5000.0;
  if (ea_eff > 0.0 && ea_eff/(5.0*boltz) > tcal) tcal = ea_eff/(5.0*boltz);
  if (tcal > umax/(20.0*boltz)) tcal = umax/(20.0*boltz);

  double diam = collide->extract(spsB[0],spsB[1],"diam");
  double tref = collide->extract(spsB[0],spsB[1],"tref");
  double mr = species[spsB[0]].mass*species[spsB[1]].mass /
    (species[spsB[0]].mass + species[spsB[1]].mass);
  double epsB = (spsB[0] == spsB[1]) ? 2.0 : 1.0;
  double kcoll = 2.0*MY_PIS*diam*diam*sqrt(2.0*boltz*tref/mr) *
    pow(tcal/tref,1.0-omegaB) / epsB;

  double cpre = r->coeff[2] *
    tgamma(zcontB+2.5-omegaB)/tgamma(zcontB+r->coeff[3]+1.5);

  double s1 = 0.0, s0 = 0.0;
  for (int k = 0; k < n; k++) {
    double u = k*du;
    s0 += den[k]*exp(-u/(boltz*tcal));
    if (mtab[i][k] > 0.0)
      s1 += den[k]*mtab[i][k]*exp(-(u-eaP)/(boltz*tcal));
  }
  double kb_model = kcoll*cpre*s1/s0;   // = model k_b * e^(+eaP/kTcal)

  double qratio =
    partition_function(r->products[0],tcal) *
    partition_function(r->products[1],tcal) /
    (partition_function(r->reactants[0],tcal) *
     partition_function(r->reactants[1],tcal));
  double target = r->reverse_A * pow(tcal,r->reverse_bf) * qratio *
    exp(-(ea_eff-eaP)/(boltz*tcal));   // = exact k_b * e^(+eaP/kTcal)

  char str[MAXLINE+128];
  if (!(kb_model > 0.0) || !(target > 0.0)) {
    sprintf(str,"Reverse reaction %s: detailed-balance table "
            "calibration failed",r->id);
    error->all(FLERR,str);
  }

  double scale = target/kb_model;
  for (int k = 0; k < n; k++) mtab[i][k] *= scale;

  // self-check: the calibrated table must reproduce the detailed-balance
  // target at EVERY temperature (the ratio of the two sides is
  // temperature independent in the continuum limit); verify at a second
  // temperature and warn if the drift exceeds 2%, which indicates
  // inconsistent level data or an under-resolved energy grid

  double tcal2 = 2.0*tcal;
  if (tcal2 > umax/(10.0*boltz)) tcal2 = umax/(10.0*boltz);
  if (tcal2 != tcal) {
    double s1b = 0.0, s0b = 0.0;
    for (int k = 0; k < n; k++) {
      double u = k*du;
      s0b += den[k]*exp(-u/(boltz*tcal2));
      if (mtab[i][k] > 0.0)
        s1b += den[k]*mtab[i][k]*exp(-(u-eaP)/(boltz*tcal2));
    }
    double kcoll2 = 2.0*MY_PIS*diam*diam*sqrt(2.0*boltz*tref/mr) *
      pow(tcal2/tref,1.0-omegaB) / epsB;
    double kb2 = kcoll2*cpre*s1b/s0b;
    double qratio2 =
      partition_function(r->products[0],tcal2) *
      partition_function(r->products[1],tcal2) /
      (partition_function(r->reactants[0],tcal2) *
       partition_function(r->reactants[1],tcal2));
    double target2 = r->reverse_A * pow(tcal2,r->reverse_bf) * qratio2 *
      exp(-(ea_eff-eaP)/(boltz*tcal2));
    double drift = kb2/target2 - 1.0;
    // the calibrated table must reproduce k_b = k_f/Keq at EVERY
    // temperature, so a large drift means the reverse rate actually
    // produced is wrong by that much -- the guarantee the B style exists
    // to provide.  Warn on a small drift (grid discretization), refuse to
    // run on a gross one, which in practice means missing or inconsistent
    // level data (e.g. no rotfile).  Written NaN-safely.
    if (!(fabs(drift) <= DB_DRIFT_FATAL)) {
      sprintf(str,"Reverse reaction %s: detailed-balance table drifts "
              "%g%% between %g and %g K, so its reverse rate cannot "
              "reproduce k_f/Keq; the species level data "
              "(elecfile/rotfile/vibfile) are missing or inconsistent",
              r->id,100.0*drift,tcal,tcal2);
      error->all(FLERR,str);
    } else if (fabs(drift) > DB_DRIFT_WARN && comm->me == 0) {
      sprintf(str,"Reverse reaction %s: detailed-balance table drifts "
              "%g%% between %g and %g K; check the species level data "
              "(elecfile/rotfile/vibfile) for consistency",
              r->id,100.0*drift,tcal,tcal2);
      error->warning(FLERR,str);
    }
  }

  // keep the pair reversible where the forward probability is capped at
  // 1: limit the forward numerator in the reverse table to the value the
  // capped forward realizes (see forward_cap); done after calibration
  // and self-check, which refer to the uncapped Arrhenius rate

  forward_cap(f,du,n,work);
  for (int k = 0; k < n; k++)
    if (den[k] > 0.0 && scale*work[k]/den[k] < mtab[i][k])
      mtab[i][k] = scale*work[k]/den[k];

  memory->destroy(leps);
  memory->destroy(lg);
  memory->destroy(num);
  memory->destroy(den);
  memory->destroy(work);
}

/* ----------------------------------------------------------------------
   build the temperature-free detailed-balance tables of a reverse
   (B-style) RECOMBINATION reaction A + B -> AB + M, the 3-body analog
   of build_db_table: the third particle M participates in the energy
   balance, so the backward probability is resolved in the TOTAL
   available energy
     w = u + e3,   e3 = eps_t + erot3 + evib3 + eelec3
   where u is the collision energy of the A,B pair and e3 collects the
   third particle's energies: eps_t is its translational energy relative
   to the pair's center of mass (with reduced mass m3*(mA+mB)/(m3+mA+mB),
   which is exactly the relative translational energy of the forward
   AB + M collision), plus its internal energies.  the probability is
     P_r(u,e3) = num(w) / (x_AB(u) * c3 * V3(w))
   with
     x_AB(u) = u^(zcontB+3/2-omegaB) (x) discrete ladders of the A,B
       pair = collision-weighted density of states of the colliding
       pair, stored in the reaction's mtab slot,
     c3 = eps_t^(1/2) * erot3^(rotdof3/2-1) * evib3^(vibdof3/2-1) the
       continuum density weights of the third particle's energies
       (evaluated at run time; the vib factor only for SMOOTH vibration),
     num(w) = (w - ea_eff)_+^(zcontF+etaF+1/2) (x) discrete ladders of
       the FORWARD pair (AB,M), with ea_eff = Ea_F + dHf: the forward
       dissociation's microcanonical numerator expressed in w (energy
       conservation folds the released energy into the threshold), and
     V3(w) = w^(nflat-1) (x) discrete ladders of M: the flat measure of
       the (u, eps_t, erot3, evib3) decomposition of w at fixed total -
       dividing by the continuum weights c3 leaves each continuous
       variable flat, so nflat counts them: the two always-present
       variables u and eps_t, plus M's rotational energy and (SMOOTH
       only) M's vibrational energy.
   dividing num by V3 makes the thermal average of P_r collapse to a
   single Laplace transform of num, so the average reproduces
     k_r(T) = A_F T^bF e^(-ea_eff/kT) q_AB(T)/(q_A(T) q_B(T))
       = k_f(T)/K_eq(T)
   at every temperature simultaneously; as in build_db_table one
   calibration temperature fixes the overall constant (VSS collision
   rates, Gamma normalizations, grid quadrature), verified at a second
   temperature.  num/V3 is stored in mtab_num on the same energy grid
   as mtab, and no temperature is evaluated at run time.
------------------------------------------------------------------------- */

void ReactBird::build_db3_table(int i)
{
  Particle::Species *species = particle->species;
  double boltz = update->boltz;
  // 60 eV in the energy unit of the unit style (J for si, erg for cgs)
  double ev60 = 60.0 * (strcmp(update->unit_style,"cgs") == 0 ?
                        1.602176634e-12 : 1.602176634e-19);

  OneReaction *r = &rlist[i];
  OneReaction *f = &rlist[r->reverse_partner];

  int spsB[2],spsF[2];
  spsB[0] = r->reactants[0]; spsB[1] = r->reactants[1];
  spsF[0] = f->reactants[0]; spsF[1] = f->reactants[1];
  int spM = r->products[1];

  int vibdiscrete = (collide->vibstyle == DISCRETE);
  int elecdiscrete = (collide->elecstyle == DISCRETE);

  // backward threshold in the total available energy w: the forward
  // barrier shifted by the reaction energy (unclamped; for a
  // dissociation whose barrier equals the well depth this is zero)

  double ea_eff = f->coeff[1] + f->coeff[4];
  double eaP = ea_eff > 0.0 ? ea_eff : 0.0;

  // energy grid: resolve the smallest vibrational quantum of either
  // pair; always build (the continuum density-of-states tables are
  // needed even without ladders)

  double theta_min = 0.0;
  for (int s = 0; s < 4; s++) {
    int sp = s < 2 ? spsB[s] : spsF[s-2];
    if (vibdiscrete)
      for (int m = 0; m < species[sp].nvibmode; m++) {
        double th = species[sp].vibtemp[m];
        if (theta_min == 0.0 || th < theta_min) theta_min = th;
      }
  }

  double umax = eaP + ev60;   // + 60 eV: Boltzmann tail to ~1e5 K
  double du;
  if (theta_min > 0.0) du = boltz*theta_min/16.0;
  else du = umax/20000.0;
  if (umax/du > 200000.0) du = umax/200000.0;
  int n = (int) (umax/du) + 2;

  mtab_du[i] = du;
  mtab_n[i] = n;

  double *num,*den,*v3,*work,*leps,*lg;
  memory->create(num,n,"react:mtab_num");
  memory->create(den,n,"react:mtab_den");
  memory->create(v3,n,"react:mtab_v3");
  memory->create(work,n,"react:mtab_work");

  // continuum internal DOF and omega of each channel; the flat-measure
  // dimension nflat and the continuum normalization exponent pcont of
  // the third particle must match the weights divided out at run time
  // in ReactTCE::attempt

  double zcontB = 0.5*(species[spsB[0]].rotdof + species[spsB[1]].rotdof);
  double zcontF = 0.5*(species[spsF[0]].rotdof + species[spsF[1]].rotdof);
  if (collide->vibstyle == SMOOTH) {
    zcontB += 0.5*(species[spsB[0]].vibdof + species[spsB[1]].vibdof);
    zcontF += 0.5*(species[spsF[0]].vibdof + species[spsF[1]].vibdof);
  }
  double omegaB = collide->extract(spsB[0],spsB[1],"omega");

  double exp_num = zcontF + f->coeff[3] + 0.5;
  double exp_den = zcontB + 1.5 - omegaB;

  int nflat = 2;
  double pcont = 1.5;
  if (species[spM].rotdof > 0) {
    nflat++;
    pcont += 0.5*species[spM].rotdof;
  }
  if (collide->vibstyle == SMOOTH && species[spM].vibdof > 0) {
    nflat++;
    pcont += 0.5*species[spM].vibdof;
  }

  for (int k = 0; k < n; k++) {
    double u = k*du;
    den[k] = (u > 0.0) ? pow(u,exp_den) : 0.0;
    num[k] = (u > ea_eff) ? pow(u-ea_eff,exp_num) : 0.0;
    v3[k] = (u > 0.0) ? pow(u,(double) (nflat-1)) : 0.0;
  }

  // convolve den with the backward pair's ladders, num with the
  // forward pair's ladders, and v3 with the third body's ladders

  // with no discrete vibrational ladder (theta_min = 0) only the
  // electronic ladders below use leps/lg
  int maxlev = (theta_min > 0.0) ? (int) (umax/(boltz*theta_min)) + 2 : 0;
  if (max_nelecstate() + 2 > maxlev) maxlev = max_nelecstate() + 2;
  memory->create(leps,maxlev,"react:mtab_leps");
  memory->create(lg,maxlev,"react:mtab_lg");

  for (int which = 0; which < 3; which++) {
    double *arr;
    int *sps,nsp;
    int spm[1] = {spM};
    if (which == 0) { arr = den; sps = spsB; nsp = 2; }
    else if (which == 1) { arr = num; sps = spsF; nsp = 2; }
    else { arr = v3; sps = spm; nsp = 1; }

    for (int s = 0; s < nsp; s++) {
      int sp = sps[s];

      if (vibdiscrete)
        for (int m = 0; m < species[sp].nvibmode; m++) {
          double th = species[sp].vibtemp[m];
          int d = species[sp].vibdegen[m] > 1 ? species[sp].vibdegen[m] : 1;
          int nlev = 0;
          double g = 1.0;
          for (int l = 0; l*th*boltz < umax && nlev < maxlev; l++) {
            leps[nlev] = l*th*boltz;
            if (l > 0) g = g*(l+d-1)/l;
            lg[nlev] = g;
            nlev++;
          }
          ladder_convolve(arr,work,n,du,nlev,leps,lg);
        }

      if (elecdiscrete && species[sp].elecdat) {
        int nlev = species[sp].elecdat->nelecstate;
        for (int l = 0; l < nlev; l++) {
          leps[l] = boltz*species[sp].elecdat->states[l].temp;
          lg[l] = species[sp].elecdat->states[l].degen;
        }
        ladder_convolve(arr,work,n,du,nlev,leps,lg);
      }
    }
  }

  mtab[i] = new double[n];
  mtab_num[i] = new double[n];
  for (int k = 0; k < n; k++) {
    mtab[i][k] = den[k];
    mtab_num[i][k] = (v3[k] > 0.0) ? num[k]/v3[k] : 0.0;
  }

  // calibrate the channel constant against the exact detailed-balance
  // rate at one temperature, as in build_db_table: the model rate is
  //   k_r(T) = kcoll_AB(T) * I1(T) / ((nflat-1)! * ZB(T) * M3(T))
  // where I1 is the Laplace transform of num, ZB the partition sum of
  // the colliding pair over its collision-weighted density of states,
  // and M3 the third particle's normalization: Gamma(d/2) (kT)^(d/2)
  // for each of its continuum weights times the partition sums of its
  // discrete ladders.  the constants must be EXACT here, not merely
  // consistent between the two calibration temperatures: the v3 table
  // divides the runtime probability by (nflat-1)! times the true flat
  // measure (its kernel omits the factorial), and the third particle's
  // energies are Gamma-distributed with the Gamma-function norms, so
  // any constant omitted from this model would NOT cancel out of the
  // thermal average of the runtime probability and would directly
  // scale the realized rate (it cancels between tcal and tcal2, so the
  // self-check below cannot catch it).  the common factor e^(-eaP/kT)
  // is removed from both sides so large thresholds cannot underflow

  double tcal = 5000.0;
  if (ea_eff > 0.0 && ea_eff/(5.0*boltz) > tcal) tcal = ea_eff/(5.0*boltz);
  if (tcal > umax/(20.0*boltz)) tcal = umax/(20.0*boltz);

  double diam = collide->extract(spsB[0],spsB[1],"diam");
  double tref = collide->extract(spsB[0],spsB[1],"tref");
  double mr = species[spsB[0]].mass*species[spsB[1]].mass /
    (species[spsB[0]].mass + species[spsB[1]].mass);
  double epsB = (spsB[0] == spsB[1]) ? 2.0 : 1.0;

  // partition sums of the third body's discrete ladders (must count
  // exactly the ladders convolved into v3)

  auto q3disc = [&](double T) {
    double q = 1.0;
    if (vibdiscrete)
      for (int m = 0; m < species[spM].nvibmode; m++) {
        if (species[spM].vibtemp[m] <= 0.0) continue;
        double x = exp(-species[spM].vibtemp[m]/T);
        int g = species[spM].vibdegen[m] > 0 ? species[spM].vibdegen[m] : 1;
        q *= pow(1.0/(1.0-x), g);
      }
    if (elecdiscrete && species[spM].elecdat) {
      double qe = 0.0;
      for (int l = 0; l < species[spM].elecdat->nelecstate; l++)
        qe += species[spM].elecdat->states[l].degen *
          exp(-species[spM].elecdat->states[l].temp/T);
      q *= qe;
    }
    return q;
  };

  // exact continuum normalization constant of the third particle's
  // Gamma-distributed energies, and the factorial by which the v3
  // kernel exceeds the true flat measure (see the comment above)

  double m3const = tgamma(1.5);
  if (species[spM].rotdof > 0) m3const *= tgamma(0.5*species[spM].rotdof);
  if (collide->vibstyle == SMOOTH && species[spM].vibdof > 0)
    m3const *= tgamma(0.5*species[spM].vibdof);
  for (int k = 2; k <= nflat-1; k++) m3const *= k;

  auto kr_model = [&](double T) {
    double s1 = 0.0, s0 = 0.0;
    for (int k = 0; k < n; k++) {
      double u = k*du;
      s0 += den[k]*exp(-u/(boltz*T));
      if (num[k] > 0.0) s1 += num[k]*exp(-(u-eaP)/(boltz*T));
    }
    double kcoll = 2.0*MY_PIS*diam*diam*sqrt(2.0*boltz*tref/mr) *
      pow(T/tref,1.0-omegaB) / epsB;
    return kcoll*s1/(s0*m3const*pow(boltz*T,pcont)*q3disc(T));
  };

  auto kr_target = [&](double T) {
    double qratio = partition_function(r->products[0],T) /
      (partition_function(r->reactants[0],T) *
       partition_function(r->reactants[1],T));
    return r->reverse_A * pow(T,r->reverse_bf) * qratio *
      exp(-(ea_eff-eaP)/(boltz*T));
  };

  double kb_model = kr_model(tcal);    // both sides carry e^(+eaP/kT)
  double target = kr_target(tcal);

  char str[MAXLINE+128];
  if (!(kb_model > 0.0) || !(target > 0.0)) {
    sprintf(str,"Reverse reaction %s: detailed-balance table "
            "calibration failed",r->id);
    error->all(FLERR,str);
  }

  double scale = target/kb_model;
  for (int k = 0; k < n; k++) mtab_num[i][k] *= scale;

  // self-check at a second temperature, as in build_db_table: the
  // model/target ratio is temperature independent by construction, so
  // drift indicates inconsistent level data or an under-resolved grid

  double tcal2 = 2.0*tcal;
  if (tcal2 > umax/(10.0*boltz)) tcal2 = umax/(10.0*boltz);
  if (tcal2 != tcal) {
    double drift = scale*kr_model(tcal2)/kr_target(tcal2) - 1.0;
    // the calibrated table must reproduce k_b = k_f/Keq at EVERY
    // temperature, so a large drift means the reverse rate actually
    // produced is wrong by that much -- the guarantee the B style exists
    // to provide.  Warn on a small drift (grid discretization), refuse to
    // run on a gross one, which in practice means missing or inconsistent
    // level data (e.g. no rotfile).  Written NaN-safely.
    if (!(fabs(drift) <= DB_DRIFT_FATAL)) {
      sprintf(str,"Reverse reaction %s: detailed-balance table drifts "
              "%g%% between %g and %g K, so its reverse rate cannot "
              "reproduce k_f/Keq; the species level data "
              "(elecfile/rotfile/vibfile) are missing or inconsistent",
              r->id,100.0*drift,tcal,tcal2);
      error->all(FLERR,str);
    } else if (fabs(drift) > DB_DRIFT_WARN && comm->me == 0) {
      sprintf(str,"Reverse reaction %s: detailed-balance table drifts "
              "%g%% between %g and %g K; check the species level data "
              "(elecfile/rotfile/vibfile) for consistency",
              r->id,100.0*drift,tcal,tcal2);
      error->warning(FLERR,str);
    }
  }

  // keep the pair reversible where the forward dissociation probability
  // is capped at 1 (see forward_cap and build_db_table)

  forward_cap(f,du,n,work);
  for (int k = 0; k < n; k++)
    if (v3[k] > 0.0 && scale*work[k]/v3[k] < mtab_num[i][k])
      mtab_num[i][k] = scale*work[k]/v3[k];

  memory->destroy(leps);
  memory->destroy(lg);
  memory->destroy(num);
  memory->destroy(den);
  memory->destroy(v3);
  memory->destroy(work);
}

/* ---------------------------------------------------------------------- */

void ReactBird::free_micro_tables()
{
  if (mtab) {
    for (int i = 0; i < mtab_nlist; i++) delete [] mtab[i];
    delete [] mtab;
  }
  if (mtab_num) {
    for (int i = 0; i < mtab_nlist; i++) delete [] mtab_num[i];
    delete [] mtab_num;
  }
  memory->destroy(mtab_du);
  memory->destroy(mtab_n);
  mtab = NULL;
  mtab_num = NULL;
  mtab_du = NULL;
  mtab_n = NULL;
  mtab_nlist = 0;
}

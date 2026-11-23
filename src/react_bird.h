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

#ifndef SPARTA_REACT_BIRD_H
#define SPARTA_REACT_BIRD_H

#include "stdio.h"
#include "math.h"
#include "react.h"
#include "particle.h"

namespace SPARTA_NS {

class ReactBird : public React {
 public:
  ReactBird(class SPARTA *, int, char **);
  ReactBird(class SPARTA *);
  virtual ~ReactBird();
  virtual void init();
  int recomb_exist(int, int);
  void ambi_check();
  virtual int attempt(Particle::OnePart *, Particle::OnePart *,
                      double, double, double, double, double &, int &) = 0;
  char *reactionID(int);
  virtual double extract_tally(int);

  // per-reaction tables of the microcanonical TCE energy factor over the
  // joint vibrational (x electronic) ladder, used with react_modify
  // partial_energy no; built by build_micro_tables() at init, shared by
  // the CPU and Kokkos TCE variants

  double **mtab;            // factor vs total collision energy, per rlist
                            // entry; NULL if the reactants have no ladders
  double *mtab_du;          // energy grid spacing (J) per table
  int *mtab_n;              // # of grid points per table
  int mtab_nlist;           // # of rlist entries tables were built for

  // 3-body detailed-balance table of a reverse recombination: the
  // calibrated forward numerator divided by the flat measure of the
  // (pair energy, third-body energy) decomposition, as a function of
  // the total available energy w = u + e3; the pair density of states
  // x_AB(u) is stored in the reaction's mtab slot; both share the mtab
  // grid; NULL for all other reactions (see build_db3_table)

  double **mtab_num;

  void build_micro_tables();
  void free_micro_tables();

  inline double vib_micro_factor(int ireact, double ecc) const {
    const double *tab = mtab[ireact];
    const int n = mtab_n[ireact];
    const double x = ecc / mtab_du[ireact];
    const int k = (int) x;
    if (k >= n-1) return tab[n-1];
    const double f = x - k;
    return (1.0-f)*tab[k] + f*tab[k+1];
  }

  inline double db3_num_factor(int ireact, double ew) const {
    const double *tab = mtab_num[ireact];
    const int n = mtab_n[ireact];
    const double x = ew / mtab_du[ireact];
    const int k = (int) x;
    if (k >= n-1) return tab[n-1];
    const double f = x - k;
    return (1.0-f)*tab[k] + f*tab[k+1];
  }

 protected:
  FILE *fp;

  // tallies for reactions

  bigint *tally_reactions,*tally_reactions_all;
  int tally_flag;

  struct OneReaction {
    int active;                    // 1 if reaction is active
    int initflag;                  // 1 if reaction params have been init
    int type;                      // reaction type = DISSOCIATION, etc
    int style;                     // reaction style = ARRHENIUS, etc
    int ncoeff;                    // # of numerical coeffs
    int nreactant,nproduct;        // # of reactants and products
    char **id_reactants,**id_products;  // species IDs of reactants/products
    int *reactants,*products;      // species indices of reactants/products
    double *coeff;                 // numerical coeffs for reaction
    char *id;                      // reaction ID (formula)
    int reverse;                   // 1 if backward rate is derived from a
                                   //   forward reaction by detailed balance
                                   //   (reaction style 'B'), else 0
    int reverse_partner;           // rlist index of the forward reaction whose
                                   //   Arrhenius rate this reverse rate derives
                                   //   from, or -1 if none found
    double reverse_bf;             // temperature exponent of the forward
                                   //   reaction, applied at the cell
                                   //   temperature in the backward
                                   //   recombination prefactor
    double a_raw;                  // Arrhenius prefactor as read from the
                                   //   file, saved before the TCE transform
                                   //   overwrites coeff[2], so a reverse
                                   //   generated or paired in a later run
                                   //   is still seeded with the raw value
    double reverse_A;              // raw Arrhenius prefactor of the forward
                                   //   reaction (stashed before the TCE
                                   //   transform), used to calibrate the
                                   //   detailed-balance table of a B-style
                                   //   exchange reaction
    int generated;                 // 1 if this reverse reaction was
                                   //   auto-generated (react_modify
                                   //   reverse auto), else 0
    int keq_flag;                  // 1 if this reverse reaction uses an
                                   //   external equilibrium-constant curve
                                   //   fit (react_modify keq_file), else 0
    double keq_coeff[5];           // Park-form fit coefficients of the
                                   //   external equilibrium constant:
                                   //   ln Keq = c0/Z + c1 + c2 ln(Z) +
                                   //   c3 Z + c4 Z^2 with Z = 10000 K / T
    double keq_resid_coeff[5];     // Park-form fit of the residual thermal
                                   //   correction ln R(T) = ln Keq_statmech
                                   //   - ln Keq_ext, applied on top of the
                                   //   detailed-balance table at run time so
                                   //   the external Keq keeps the energy-
                                   //   resolved shape (only R is thermal)
  };

  OneReaction *rlist;              // list of all reactions read from file

  // pieces of the TCE probability of a two-body channel shared by the
  // run-time ReactTCE::channel_prob() and the table-time
  // pair_channel_prob(), so the two cannot drift apart; the operation
  // order is that of the original expressions, so results are unchanged
  // (the KOKKOS twin, ReactTCEKokkos::channel_prob_kk, mirrors them)

  // energy above the threshold: the activation energy, or the
  // endothermicity where it is larger
  inline double tce_excess(const OneReaction *r, double ecc) const {
    if (r->coeff[1] > ((-1)*r->coeff[4])) return ecc - r->coeff[1];
    return ecc + r->coeff[4];
  }

  // energy factor: the microcanonical table when one is built and used,
  // else the standard (ecc-Ea)^(eta-1+omega) (1-Ea/ecc)^(z+3/2-omega)
  inline double tce_efactor(int ireact, double ecc, double z,
                            int use_table) const {
    const OneReaction *r = &rlist[ireact];
    if (use_table && mtab && mtab[ireact])
      return vib_micro_factor(ireact,ecc);
    return pow(ecc-r->coeff[1],r->coeff[3]-1+r->coeff[5]) *
      pow(1.0-r->coeff[1]/ecc,z+1.5-r->coeff[5]);
  }

  // prefactor times the Gamma-function ratio; the denominator is clamped
  // at 1e-6 below the TCE exponent bound (warned by check_tce_bounds)
  inline double tce_prefactor(const OneReaction *r, double z,
                              double prefactor) const {
    double g = tgamma(z+r->coeff[3]+1.5);
    if (g < 1.0e-6) g = 1.0e-6;
    return prefactor * tgamma(z+2.5-r->coeff[5]) / g;
  }

  int maxlist;                     // max # of reactions in rlist

  // all reactions a pair of reactant species is part of

  struct ReactionIJ {
    int *list;       // N-length list of rlist indices
                     //   for reactions defined for this IJ pair,
                     //   just a ptr into sub-section of long list_ij vector
                     //   for all pairs
    int *sp2recomb;  // Nspecies-length list of rlist indices
                     //   for recomb reactions defined for this IJ pair,
                     //   one index for all 3rd particle species,
                     //   just a ptr into sub-section of long sp2recomb_ij
                     //   vector for all pairs which have recomb reactions
    int n;           // # of reactions in list
  };

  ReactionIJ **reactions;     // reaction info for all IJ pairs of species
  int *list_ij;               // chunks of rlist indices,
                              //   one chunk per IJ pair,
                              //   stored in contiguous vector,
                              //   length of each chunk is # of IJ reactions
                              // pointed into by reactions[i][k].list
  int *sp2recomb_ij;          // chunks of rlist indices,
                              //   one chunk per IJ pair that has
                              //     recombination reactions,
                              //   stored in contiguous vector
                              //   length of each chunk is # of species
                              // pointed into by reactions[i][k].sp2recomb

  // equilibrium-constant curve fits read from react_modify keq_file

  struct KeqFit {
    int reactants[2];                // species indices of the FORWARD
    int products[3];                 //   reaction the fit belongs to
    int nreactant,nproduct;
    double coeff[5];                 // Park-form coefficients
    int used;                        // 1 once matched to a reaction
  };

  KeqFit *keqfits;
  int nkeqfits;

  int generated_flag;              // 1 once auto-reverses were generated
  int tce_style;                   // 1 for the tce styles (set by their
                                   //   constructors), whose Arrhenius
                                   //   probability carries the TCE Gamma
                                   //   normalization and which implement
                                   //   detailed-balance reverse reactions

  void readfile(char *);
  int readone(char *, char *, int &, int &);
  void check_duplicate();
  int tce_bounds_checked;          // 1 if check_tce_bounds() already ran
  int reverse_off_warned;          // # of B lines switched off at the last
                                   //   warning, to warn once per change
  void check_tce_bounds();
  int reverse_modes_checked;       // 1 if check_reverse_modes() already ran
  void check_reverse_modes();
  double tce_eta_bound(OneReaction *);
  void report_tce_prob_flags(int, int, int);
  virtual void grow_tallies();
  double partition_function(int, double);
  void forward_cap(OneReaction *, double, int, double *);
  double pair_channel_prob(int, double);
  double pair_channel_z(OneReaction *);
  int max_nelecstate();
  void build_db_table(int);
  void build_db3_table(int);
  void generate_reverses(int *);
  void read_keq_file();
  void assign_keq_fits();
  void fit_keq_residual(int);

  inline double keq_eval(const double *c, double T) const {
    double z = 10000.0/T;
    return exp(c[0]/z + c[1] + c[2]*log(z) + c[3]*z + c[4]*z*z);
  }
  void print_reaction(char *, char *);
  void print_reaction(OneReaction *);
  void print_reaction_ambipolar(OneReaction *);
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: React tce can only be used with collide vss

Self-explanatory.

E: Ionization and recombination reactions are not yet implemented

This error conditions will be removed after those reaction styles are
fully implemented.

E: Unknown outcome in reaction

The specified type of the reaction is not encoded in the reaction
style.

E: Cannot open reaction file %s

Self-explanatory.

E: Invalid reaction formula in file

Self-explanatory.

E: Invalid reaction type in file

Self-explanatory.

E: Invalid reaction style in file

Self-explanatory.

E: Invalid reaction coefficients in file

Self-explanatory.

*/

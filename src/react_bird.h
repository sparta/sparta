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

#ifndef SPARTA_REACT_BIRD_H
#define SPARTA_REACT_BIRD_H

#include "stdio.h"
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
                      double, double, double, double &, int &) = 0;
  char *reactionID(int);
  virtual double extract_tally(int);

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
  };

  OneReaction *rlist;              // list of all reactions read from file
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

  void readfile(char *);
  int readone(char *, char *, int &, int &);
  void check_duplicate();
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

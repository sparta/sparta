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

#ifdef SURF_REACT_CLASS

SurfReactStyle(prob,SurfReactProb)

#else

#ifndef SPARTA_SURF_REACT_PROB_H
#define SPARTA_SURF_REACT_PROB_H

#include "surf_react.h"

namespace SPARTA_NS {

class SurfReactProb : public SurfReact {
 public:
  SurfReactProb(class SPARTA *, int, char **);
  SurfReactProb(class SPARTA *sparta) : SurfReact(sparta) {}
  virtual ~SurfReactProb();
  virtual void init();
  int react(Particle::OnePart *&, int, double *, Particle::OnePart *&, int &);
  char *reactionID(int);
  int match_reactant(char *, int);
  int match_product(char *, int);

  // reaction info, as read from file

  struct OneReaction {
    int active;                    // 1 if reaction is active
    int type;                      // reaction type = DISSOCIATION, etc
    int style;                     // reaction style = ARRHENIUS, etc
    int ncoeff;                    // # of numerical coeffs
    int nreactant,nproduct;        // # of reactants and products
    char **id_reactants,**id_products;  // species IDs of reactants/products
    int *reactants,*products;      // species indices of reactants/products
    double *coeff;                 // numerical coeffs for reaction
    char *id;                      // reaction ID (formula)
  };

 protected:
  class RanKnuth *random;     // RNG for reaction probabilities

  OneReaction *rlist;              // list of all reactions read from file
  int nlist_prob;                  // # of reactions read from file
  int maxlist_prob;                // max # of reactions in rlist

  // possible reactions a reactant species is part of

  struct ReactionI {
    int *list;           // list of indices into rlist, ptr into indices
    int n;               // # of reactions in list
  };

  ReactionI *reactions;       // reactions for all species
  int *indices;               // master list of indices

  virtual void init_reactions();
  void readfile(char *);
  int readone(char *, char *, int &, int &);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

*/

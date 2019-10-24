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

#ifndef SPARTA_SURF_REACT_H
#define SPARTA_SURF_REACT_H

#include "pointers.h"
#include "particle.h"

namespace SPARTA_NS {

class SurfReact : protected Pointers {
 public:
  char *id;
  char *style;

  int vector_flag;          // 0/1 if compute_vector() function exists
  int size_vector;          // length of global vector
  int nlist;                // # of reactions defined or read from file

  SurfReact(class SPARTA *, int, char **);
  virtual ~SurfReact();
  virtual void init();
  virtual int react(Particle::OnePart *&, double *, Particle::OnePart *&) = 0;

  void tally_update();
  virtual char *reactionID(int);
  double compute_vector(int i);

 protected:
  FILE *fp;

  // tallies for reactions

  int nsingle,ntotal;
  double one[2],all[2];
  int *tally_single,*tally_total;
  int *tally_single_all,*tally_total_all;
  int tally_two_flag,tally_single_flag,tally_total_flag;

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

  OneReaction *rlist;              // list of all reactions read from file
  int maxlist;                     // max # of reactions in rlist

  // possible reactions a reactant species is part of

  struct ReactionI {
    int *list;           // list of indices into rlist, ptr into indices
    int n;               // # of reactions in list
  };

  ReactionI *reactions;       // reactions for all species
  int *indices;               // master list of indices

  void init_reactions();
  void readfile(char *);
  int readone(char *, char *, int &, int &);
};

}

#endif

/* ERROR/WARNING messages:

E: Surf_react ID must be alphanumeric or underscore characters

Self-explanatory.

*/

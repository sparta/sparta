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

#ifndef SPARTA_SURF_REACT_H
#define SPARTA_SURF_REACT_H

#include "pointers.h"
#include "particle.h"

namespace SPARTA_NS {

class SurfReact : protected Pointers {
 public:
  char *id;
  char *style;

  int nlist;                // # of reactions defined or read from file
  int vector_flag;          // 0/1 if compute_vector() function exists
  int size_vector;          // length of global vector
  int kokkosable;           // 1 if Kokkos version
  int copy,copymode;        // 1 if copy of class, used by Kokkos

  SurfReact(class SPARTA *, int, char **);
  SurfReact(class SPARTA *sparta) : Pointers(sparta) {}
  virtual ~SurfReact();
  virtual void init();
  virtual int react(Particle::OnePart *&, int, double *,
                    Particle::OnePart *&, int &) = 0;
  virtual char *reactionID(int) = 0;
  virtual int match_reactant(char *, int) = 0;
  virtual int match_product(char *, int) = 0;

  virtual void tally_reset();
  virtual void tally_update();
  double compute_vector(int i);

 protected:
  FILE *fp;

  // tallies for reactions
  // nsingle = all reactions in one step
  // ntotal = cumulative nsingle across all steps
  // tally_single = per-reaction counts in one step
  // tally_all = cumulative tally_single across all steps
  // 3 flags used in compute_vector() to minimize AllReduce calls

  int nsingle,ntotal;
  double one[2],all[2];
  int *tally_single,*tally_total;
  int *tally_single_all,*tally_total_all;
  int tally_two_flag,tally_single_flag,tally_total_flag;
};

}

#endif

/* ERROR/WARNING messages:

E: Surf_react ID must be alphanumeric or underscore characters

Self-explanatory.

*/

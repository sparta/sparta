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

SurfReactStyle(global,SurfReactGlobal)

#else

#ifndef SPARTA_SURF_REACT_GLOBAL_H
#define SPARTA_SURF_REACT_GLOBAL_H

#include "surf_react.h"

namespace SPARTA_NS {

class SurfReactGlobal : public SurfReact {
 public:
  SurfReactGlobal(class SPARTA *, int, char **);
  SurfReactGlobal(class SPARTA *sparta) : SurfReact(sparta) {}
  virtual ~SurfReactGlobal();
  int react(Particle::OnePart *&, int, double *, Particle::OnePart *&, int &);
  char *reactionID(int);
  int match_reactant(char *, int) {return 1;}
  int match_product(char *, int) {return 1;}

 protected:
  double prob_create,prob_destroy;
  class RanKnuth *random;     // RNG for reaction probabilities
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

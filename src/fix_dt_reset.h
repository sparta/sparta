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

#ifdef FIX_CLASS

FixStyle(dt/reset,FixDtReset)

#else

#ifndef SPARTA_FIX_DT_RESET_H
#define SPARTA_FIX_DT_RESET_H

#include "fix.h"

namespace SPARTA_NS {

class FixDtReset : public Fix {
 public:
  FixDtReset(class SPARTA *, int, char **);
  ~FixDtReset();
  int setmask();
  void init();
  void end_of_step();
  double compute_scalar();
  double compute_vector(int);
  void reallocate();

 protected:
  int step_which,step_index,resetflag;
  int maxgrid;
  double weight;
  double dtmin,dtmax,dtave,dtnew;

  char *id_step;
  class Compute *cstep;
  class Fix *fstep;

  double *gridstep;
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

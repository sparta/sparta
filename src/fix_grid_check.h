/* ----------------------------------------------------------------------
   DSMC - Sandia parallel DSMC code
   www.sandia.gov/~sjplimp/dsmc.html
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2011) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level DSMC directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(grid/check,FixGridCheck)

#else

#ifndef DSMC_FIX_GRID_CHECK_H
#define DSMC_FIX_GRID_CHECK_H

#include "fix.h"

namespace DSMC_NS {

class FixGridCheck : public Fix {
 public:
  FixGridCheck(class DSMC *, int, char **);
  ~FixGridCheck() {}
  int setmask();
  void init();
  void end_of_step();
  double compute_scalar();

 private:
  int nevery,nflag;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running DSMC to see the offending line.

*/

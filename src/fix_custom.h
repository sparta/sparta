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

#ifdef FIX_CLASS

FixStyle(custom,FixCustom)

#else

#ifndef SPARTA_FIX_CUSTOM_H
#define SPARTA_FIX_CUSTOM_H

#include "fix.h"
#include "surf.h"

namespace SPARTA_NS {

class FixCustom : public Fix {
 public:
  FixCustom(class SPARTA *, int, char **);
  virtual ~FixCustom();
  int setmask();
  virtual void init();
  virtual void end_of_step();

 private:
  int mode,action;
  int vindex,vstyle;
  int groupbit;
  char *fname,*aname,*vname;

  class Variable *variable;
  class Mixture *mixture;
  class Region *region;
  class Custom *custom;

  int ctype,csize,cindex,ccol;
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

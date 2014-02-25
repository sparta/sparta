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

#ifdef FIX_CLASS

FixStyle(print,FixPrint)

#else

#ifndef SPARTA_FIX_PRINT_H
#define SPARTA_FIX_PRINT_H

#include "stdio.h"
#include "fix.h"

namespace SPARTA_NS {

class FixPrint : public Fix {
 public:
  FixPrint(class SPARTA *, int, char **);
  ~FixPrint();
  int setmask();
  void end_of_step();

 private:
  int me,screenflag;
  FILE *fp;
  char *string,*copy,*work;
  int maxcopy,maxwork;
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/

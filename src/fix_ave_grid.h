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

FixStyle(ave/grid,FixAveGrid)

#else

#ifndef LMP_FIX_AVE_GRID_H
#define LMP_FIX_AVE_GRID_H

#include "stdio.h"
#include "fix.h"

namespace DSMC_NS {

class FixAveGrid : public Fix {
 public:
  FixAveGrid(class DSMC *, int, char **);
  ~FixAveGrid();
  int setmask();
  void init();
  void setup(int);
  void end_of_step();
  double memory_usage();

 private:
  int nvalues;
  int nevery,nrepeat,irepeat;
  bigint nvalid;
  int *which,*argindex,*value2index;
  char **ids;

  int *pcount;

  bigint nextvalid();
};

}

#endif
#endif

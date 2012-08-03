/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2012) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(ave/grid2,FixAveGrid2)

#else

#ifndef LMP_FIX_AVE_GRID2_H
#define LMP_FIX_AVE_GRID2_H

#include "fix.h"

namespace SPARTA_NS {

class FixAveGrid2 : public Fix {
 public:
  FixAveGrid2(class SPARTA *, int, char **);
  ~FixAveGrid2();
  int setmask();
  void init();
  void setup();
  void end_of_step();
  double memory_usage();

 private:
  int nvalues,maxvalues;
  int nevery,nrepeat,irepeat,nsample,ave;
  bigint nvalid;
  int *which,*argindex,*value2index;
  char **ids;

  int nglocal;
  double *accvec;
  double **accarray;

  int *normacc;        // 1 if Ith value triggers one-time norm accumulation
  int *normindex;      // index of norm vector for Ith value, -1 if none
  double **norms;      // pointers to accumulated norms
  double **cfv_norms;  // pointers to snapshot norms by compute,fix,variable
  int nnorm;           // # of norm pointers in norms and cfv_norms

  void options(int, char **);
  void grow();
  bigint nextvalid();
};

}

#endif
#endif

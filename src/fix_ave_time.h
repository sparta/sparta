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

FixStyle(ave/time,FixAveTime)

#else

#ifndef SPARTA_FIX_AVE_TIME_H
#define SPARTA_FIX_AVE_TIME_H

#include "stdio.h"
#include "fix.h"

namespace SPARTA_NS {

class FixAveTime : public Fix {
 public:
  FixAveTime(class SPARTA *, int, char **);
  ~FixAveTime();
  int setmask();
  void init();
  void setup();
  void end_of_step();
  double compute_scalar();
  double compute_vector(int);
  double compute_array(int,int);

 private:
  int me,nvalues,maxvalues;
  int nrepeat,nfreq,irepeat;
  bigint nvalid;
  int *which,*argindex,*value2index,*offcol;
  char **ids;
  FILE *fp;
  int nrows;

  int ave,nwindow,nsum,startstep,mode;
  int noff;
  int *offlist;
  char *title1,*title2,*title3;

  int norm,iwindow,window_limit;
  double *vector;
  double *vector_total;
  double **vector_list;
  double *column;
  double **array;
  double **array_total;
  double ***array_list;
  double **norms;

  void invoke_scalar(bigint);
  void invoke_vector(bigint);
  void options(int, char **);
  void grow();
  bigint nextvalid();
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/

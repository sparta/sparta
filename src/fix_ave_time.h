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
  void options(int, int, char **);
  void grow();
  bigint nextvalid();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Invalid fix ave/time off column

Self-explanatory.

E: Compute ID for fix ave/time does not exist

Self-explanatory.

E: Fix ave/time compute does not calculate a scalar

Self-explanatory.

E: Fix ave/time compute does not calculate a vector

Self-explanatory.

E: Fix ave/time compute vector is accessed out-of-range

The index for the vector is out of bounds.

E: Fix ave/time compute does not calculate an array

Self-explanatory.

E: Fix ave/time compute array is accessed out-of-range

An index for the array is out of bounds.

E: Fix ID for fix ave/time does not exist

Self-explanatory.

E: Fix ave/time fix does not calculate a scalar

Self-explanatory.

E: Fix ave/time fix does not calculate a vector

Self-explanatory.

E: Fix ave/time fix vector is accessed out-of-range

The index for the vector is out of bounds.

E: Fix for fix ave/time not computed at compatible time

Fixes generate their values on specific timesteps.  Fix ave/time
is requesting a value on a non-allowed timestep.

E: Fix ave/time fix does not calculate an array

Self-explanatory.

E: Fix ave/time fix array is accessed out-of-range

An index for the array is out of bounds.

E: Variable name for fix ave/time does not exist

Self-explanatory.

E: Fix ave/time variable is not equal-style variable

Self-explanatory.

E: Fix ave/time cannot use variable with vector mode

Variables produce scalar values.

E: Fix ave/time columns are inconsistent lengths

Self-explanatory.

E: Cannot open fix ave/time file %s

The specified file cannot be opened.  Check that the path and name are
correct.

*/

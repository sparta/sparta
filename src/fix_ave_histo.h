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

FixStyle(ave/histo,FixAveHisto)

#else

#ifndef SPARTA_FIX_AVE_HISTO_H
#define SPARTA_FIX_AVE_HISTO_H

#include <stdio.h>
#include "fix.h"

namespace SPARTA_NS {

class FixAveHisto : public Fix {
 public:
  FixAveHisto(class SPARTA *, int, char **);
  virtual ~FixAveHisto();
  int setmask();
  void init();
  void setup();
  virtual void end_of_step();
  double compute_vector(int);
  double compute_array(int,int);

 protected:
  int me,nvalues;
  int nrepeat,nfreq,irepeat;
  bigint nvalid;
  int *which,*argindex,*value2index;
  char **ids;
  FILE *fp;
  double lo,hi,binsize,bininv;
  int regionflag,mixflag,groupflag;
  int iregion,imix,groupbit;
  int kind,beyond,overwrite;
  long filepos;

  double stats[4],stats_total[4],stats_all[4];
  double **stats_list;

  int nbins;
  double *bin,*bin_total,*bin_all;
  double **bin_list;
  double *coord;

  double *vector;
  int maxvector;

  int ave,nwindow,startstep,mode;
  char *title1,*title2,*title3;
  int iwindow,window_limit;

  // data used by ave/histo/weight

  int weightflag;
  double weight;
  double *weights;
  double *vectorwt;
  int maxvectorwt;

  // methods

  virtual void bin_one(double);
  virtual void bin_vector(int, double *, int);
  virtual void bin_particles(int, int);
  virtual void bin_particles(double *, int);
  virtual void bin_grid_cells(double *, int);

  virtual void calculate_weights() {}

  void options(int, int, char **);
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

E: Compute ID for fix ave/histo does not exist

Self-explanatory.

E: Fix ID for fix ave/histo does not exist

Self-explanatory.

E: Fix ave/histo input is invalid compute

Self-explanatory.

E: Fix ave/histo input is invalid fix

Self-explanatory.

E: Fix ave/histo input is invalid variable

Self-explanatory.

E: Fix ave/histo inputs are not all global, peratom, or local

All inputs in a single fix ave/histo command must be of the
same style.

E: Fix ave/histo cannot input per-atom values in scalar mode

Self-explanatory.

E: Fix ave/histo cannot input local values in scalar mode

Self-explanatory.

E: Fix ave/histo compute does not calculate a global scalar

Self-explanatory.

E: Fix ave/histo compute does not calculate a global vector

Self-explanatory.

E: Fix ave/histo compute vector is accessed out-of-range

Self-explanatory.

E: Fix ave/histo compute does not calculate a global array

Self-explanatory.

E: Fix ave/histo compute array is accessed out-of-range

Self-explanatory.

E: Fix ave/histo compute does not calculate per-atom values

Self-explanatory.

E: Fix ave/histo compute does not calculate a per-atom vector

Self-explanatory.

E: Fix ave/histo compute does not calculate a per-atom array

Self-explanatory.

E: Fix ave/histo compute does not calculate local values

Self-explanatory.

E: Fix ave/histo compute does not calculate a local vector

Self-explanatory.

E: Fix ave/histo compute does not calculate a local array

Self-explanatory.

E: Fix ave/histo fix does not calculate a global scalar

Self-explanatory.

E: Fix ave/histo fix does not calculate a global vector

Self-explanatory.

E: Fix ave/histo fix vector is accessed out-of-range

Self-explanatory.

E: Fix for fix ave/histo not computed at compatible time

Fixes generate their values on specific timesteps.  Fix ave/histo is
requesting a value on a non-allowed timestep.

E: Fix ave/histo fix does not calculate a global array

Self-explanatory.

E: Fix ave/histo fix array is accessed out-of-range

Self-explanatory.

E: Fix ave/histo fix does not calculate per-atom values

Self-explanatory.

E: Fix ave/histo fix does not calculate a per-atom vector

Self-explanatory.

E: Fix ave/histo fix does not calculate a per-atom array

Self-explanatory.

E: Fix ave/histo fix does not calculate local values

Self-explanatory.

E: Fix ave/histo fix does not calculate a local vector

Self-explanatory.

E: Fix ave/histo fix does not calculate a local array

Self-explanatory.

E: Variable name for fix ave/histo does not exist

Self-explanatory.

E: Error writing file header

Something in the output to the file triggered an error.

E: Invalid timestep reset for fix ave/histo

Resetting the timestep has invalidated the sequence of timesteps this
fix needs to process.

E: Error writing out histogram data

Something in the output to the file triggered an error.

E: Cannot open fix ave/histo file %s

The specified file cannot be opened.  Check that the path and name are
correct.

*/

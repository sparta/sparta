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

FixStyle(ave/grid,FixAveGrid)

#else

#ifndef LMP_FIX_AVE_GRID_H
#define LMP_FIX_AVE_GRID_H

#include "fix.h"

namespace SPARTA_NS {

class FixAveGrid : public Fix {
 public:
  FixAveGrid(class SPARTA *, int, char **);
  ~FixAveGrid();
  int setmask();
  void init();
  void setup();
  void end_of_step();
  int pack_grid_one(int, char *, int);
  int unpack_grid_one(int, char *);
  void compress_grid();
  double memory_usage();

 private:
  int nvalues,maxvalues;
  int nrepeat,irepeat,nsample,ave;
  bigint nvalid;
  int *which,*argindex,*value2index;
  char **ids;

  int nglocal;               // # of owned grid cells
  int nglocalmax;            // max size of per-cell vectors/arrays
  double *vector;            // tally vector when ave = RUNNING
  double **array;            // tally array when ave = RUNNING
  double **array_extra;      // extra tally array for computes that need it
  double **norm_extra;       // extra norm array for computes that need it

  int *normacc;        // 1 if Ith value triggers one-time norm accumulation
  int *normindex;      // index of norm vector for Ith value, -1 if none
  double **norms;      // pointers to accumulated norms
  int nnorm;           // # of norm pointers in norms

  int pack_one(int, char *, int);
  int unpack_one(char *, int);
  void options(int, char **);
  void grow();
  bigint nextvalid();
  void allocate(int);
  void grow_percell(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Compute ID for fix ave/grid does not exist

Self-explanatory.

E: Fix ave/grid compute does not calculate per-grid values

Self-explanatory.

E: Fix ave/grid compute does not calculate a per-grid vector

Self-explanatory.

E: Fix ave/grid compute does not calculate a per-grid array

Self-explanatory.

E: Fix ave/grid compute array is accessed out-of-range

Self-explanatory.

E: Fix ID for fix ave/grid does not exist

Self-explanatory.

E: Fix ave/grid fix does not calculate per-grid values

Self-explanatory.

E: Fix ave/grid fix does not calculate a per-grid vector

Self-explanatory.

E: Fix ave/grid fix does not calculate a per-grid array

Self-explanatory.

E: Fix ave/grid fix array is accessed out-of-range

Self-explanatory.

E: Fix for fix ave/grid not computed at compatible time

Fixes generate values on specific timesteps.  Fix ave/grid is
requesting a value on a non-allowed timestep.

E: Variable name for fix ave/grid does not exist

Self-explanatory.

E: Fix ave/grid variable is not grid-style variable

Self-explanatory.

*/

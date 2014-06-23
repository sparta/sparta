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

FixStyle(ave/surf,FixAveSurf)

#else

#ifndef LMP_FIX_AVE_SURF_H
#define LMP_FIX_AVE_SURF_H

#include "fix.h"

namespace SPARTA_NS {

class FixAveSurf : public Fix {
 public:
  FixAveSurf(class SPARTA *, int, char **);
  ~FixAveSurf();
  int setmask();
  void init();
  void setup();
  void end_of_step();
  double memory_usage();

 private:
  int nvalues,maxvalues;
  int nrepeat,irepeat,nsample,ave;
  bigint nvalid;
  int *which,*argindex,*value2index;
  char **ids;

  int nslocal;               // # of surfs I own
  double *accvec;            // accumulation vector
  double **accarray;         // accumulation array

  int nsurf;               // # of global surfs, lines or triangles
  int nlocal;              // # of local surfs
  int maxlocal;            // # of local surfs currently allocated
  int *glob2loc;           // glob2loc[I] = local index of Ith global surf
  int *loc2glob;           // loc2glob[I] = global index of Ith local surf

  double *vec_local;
  double **array_local;

  void options(int, char **);
  void grow();
  void grow_local();
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

E: Compute ID for fix ave/surf does not exist

Self-explanatory.  

E: Fix ave/surf compute does not calculate per-surf values

Self-explanatory.  

E: Fix ave/surf compute does not calculate a per-surf vector

Self-explanatory.  

E: Fix ave/surf compute does not calculate a per-surf array

Self-explanatory.  

E: Fix ave/surf compute array is accessed out-of-range

Self-explanatory.  

E: Fix ID for fix ave/surf does not exist

Self-explanatory.  

E: Fix ave/surf fix does not calculate per-surf values

Self-explanatory.  

E: Fix ave/surf fix does not calculate a per-surf vector

Self-explanatory.  

E: Fix ave/surf fix does not calculate a per-surf array

Self-explanatory.  

E: Fix ave/surf fix array is accessed out-of-range

Self-explanatory.  

E: Fix for fix ave/surf not computed at compatible time

Fixes generate their values on specific timesteps.  Fix ave/surf is
requesting a value on a non-allowed timestep.

E: Variable name for fix ave/surf does not exist

Self-explanatory.  

E: Fix ave/surf variable is not surf-style variable

Self-explanatory.

*/

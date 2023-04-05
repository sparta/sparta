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

#ifdef COMPUTE_CLASS

ComputeStyle(reduce,ComputeReduce)

#else

#ifndef SPARTA_COMPUTE_REDUCE_H
#define SPARTA_COMPUTE_REDUCE_H

#include "compute.h"

namespace SPARTA_NS {

class ComputeReduce : public Compute {
 public:
  ComputeReduce(class SPARTA *, int, char **);
  ~ComputeReduce();
  void init();
  double compute_scalar();
  void compute_vector();
  bigint memory_usage();

 protected:
  int me,nprocs;
  int mode,nvalues,iregion;
  double area_total;
  int *which,*argindex,*flavor,*value2index;
  char **ids;
  double *onevec;
  int *replace,*indices,*owner;
  int index;

  char *subsetID;
  int isubset;
  int *s2g;
  int gridgroupbit,surfgroupbit;

  int maxparticle,maxgrid;
  double *varparticle,*vargrid;
  double *areasurf;

  struct Pair {
    double value;
    int proc;
  };
  Pair pairme,pairall;

  double compute_one(int, int);
  bigint count_included();
  double area_per_surf();
  void combine(double &, double, int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Compute reduce replace requires min or max mode

Self-explanatory.

E: Invalid replace values in compute reduce

Self-explanatory.

E: Compute ID for compute reduce does not exist

Self-explanatory.

E: Compute reduce compute does not calculate a per-particle vector

This is necessary if no column index is used to specify the compute.

E: Compute reduce compute does not calculate a per-particle array

This is necessary if a column index is used to specify the compute.

E: Compute reduce compute array is accessed out-of-range

An index for the array is out of bounds.

E: Compute reduce compute does not calculate a per-grid vector

This is necessary if no column index is used to specify the compute.

E: Compute reduce compute does not calculate a per-grid array

This is necessary if a column index is used to specify the compute.

E: Compute reduce compute calculates global or surf values

The compute reduce command does not operate on this kind of values.
The variable command has special functions that can reduce global
values.

E: Fix ID for compute reduce does not exist

Self-explanatory.

E: Compute reduce fix does not calculate a per-particle vector

This is necessary if no column index is used to specify the fix.

E: Compute reduce fix does not calculate a per-particle array

This is necessary if a column index is used to specify the fix.

E: Compute reduce fix array is accessed out-of-range

An index for the array is out of bounds.

E: Compute reduce fix does not calculate a per-grid vector

This is necessary if no column index is used to specify the fix.

E: Compute reduce fix does not calculate a per-grid array

This is necessary if a column index is used to specify the fix.

E: Compute reduce fix does not calculate a per-surf vector

This is necessary if no column index is used to specify the fix.

E: Compute reduce fix does not calculate a per-surf array

This is necessary if a column index is used to specify the fix.

E: Compute reduce fix calculates global values

A fix that calculates peratom or local values is required.

E: Variable name for compute reduce does not exist

Self-explanatory.

E: Compute reduce variable is not particle-style variable

This is the only style of variable that can be reduced.

E: Fix used in compute reduce not computed at compatible time

Fixes generate their values on specific timesteps.  Compute reduce is
requesting a value on a non-allowed timestep.

*/

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

ComputeStyle(count,ComputeCount)

#else

#ifndef SPARTA_COMPUTE_COUNT_H
#define SPARTA_COMPUTE_COUNT_H

#include "compute.h"

namespace SPARTA_NS {

class ComputeCount : public Compute {
 public:
  ComputeCount(class SPARTA *, int, char **);
  ~ComputeCount();
  void init();
  double compute_scalar();
  void compute_vector();

 protected:
  int nvalues,maxvalues;     // number of values to count
  int *spmix;                // value is SPECIES or MIXTURE group
  int *index;                // index of species or mixture
  int *indexgroup;           // index of mixture group
  int *mixgroups;            // # of groups in mixture at constructor time

  int maxspecies;            // size of count vector
  int *count;                // per-species count of particles on this proc
  bigint *onevec;            // this proc's contribution to each value
  bigint *sumvec;            // values summed across all procs

  bigint lasttally;          // last timestep particles were counted

  void per_species_tally();
  void allocate(int);
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

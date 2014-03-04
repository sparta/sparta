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
  int me;
  int mode,nvalues,iregion;
  int *which,*argindex,*flavor,*value2index;
  char **ids;
  double *onevec;
  int *replace,*indices,*owner;
  int index;
  char *idregion;

  int maxparticle;
  double *varparticle;

  struct Pair {
    double value;
    int proc;
  };
  Pair pairme,pairall;

  double compute_one(int, int);
  bigint count(int);
  void combine(double &, double, int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/

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

ComputeStyle(lambda/grid,ComputeLambdaGrid)

#else

#ifndef SPARTA_COMPUTE_LAMBDA_GRID_H
#define SPARTA_COMPUTE_LAMBDA_GRID_H

#include "compute.h"

namespace SPARTA_NS {

class ComputeLambdaGrid : public Compute {
 public:
  ComputeLambdaGrid(class SPARTA *, int, char **);
  ~ComputeLambdaGrid();
  void init();
  virtual void compute_per_grid();
  virtual void reallocate();
  bigint memory_usage();

 protected:
  int nglocal;
  int nrhowhich,tempwhich,kflag;

  char *id_nrho,*id_temp;
  int nrhoindex,tempindex;
  class Compute *cnrho,*ctemp;
  class Fix *fnrho,*ftemp;
  double *nrho,*temp;

  char *species;
  double dref,tref,omega,prefactor;
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

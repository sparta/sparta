/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
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
  int nglocal,nvalues,nparams,tmax,noutputs;
  int tempwhich;
  int lambdaflag,tauflag;
  int knallflag,knxflag,knyflag,knzflag,knanyflag;

  char **ids_nrho;           // ID/name of compute,fix,variable to access
  int *nrhowhich;            // COMPUTE or FIX or VARIABLE
  int *nrhoindex;            // which column from compute or fix to access
  int *value2index;          // index of compute,fix,variable
  int *post_process;         // 1 if need compute->post_process() on value

  int ntotal;                // total # of columns in tally array

                             // used when normalizing tallies
  int *nmap;                 // # of tally quantities for each value
                             //   these may not be unique
  int **map;                 // indices of non-unique tally quantities
                             //   in tally, for each value

                             // used when accumulating tallies
  int *numap;                // # of unique tally quantities for each value
  int **umap;                // indices of tally quants in tally for each value
  int **uomap;               // indices of corresponding quantities (0 to N-1)
                             //   in compute/fix tally array, for each value
  int *output_order;

  char *id_temp;
  int tempindex;
  class Compute *cnrho,*ctemp;
  class Fix *fnrho,*ftemp;
  double **nrho,*temp,**lambdainv,**tauinv;
  double **array_grid1;

  char *species;
  double dref,tref,omega,mj,mk,mr;
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

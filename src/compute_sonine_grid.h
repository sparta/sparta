/* ----------------------------------------------------------------------
   DSMC - Sandia parallel DSMC code
   www.sandia.gov/~sjplimp/dsmc.html
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2011) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level DSMC directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(sonine/grid,ComputeSonineGrid)

#else

#ifndef DSMC_COMPUTE_SONINE_GRID_H
#define DSMC_COMPUTE_SONINE_GRID_H

#include "compute.h"

namespace DSMC_NS {

class ComputeSonineGrid : public Compute {
 public:
  ComputeSonineGrid(class DSMC *, int, char **);
  ~ComputeSonineGrid();
  void init();
  void compute_per_grid();
  bigint memory_usage();

 private:
  int imix,nvalues,ngroups,npergroup,ntotal;
  int *which,*moment,*order;

  double **vave,***avave;
  int *count,**acount;

  double *sonine_vector;
  double **sonine_array;
};

}

#endif
#endif

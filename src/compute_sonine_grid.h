/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2012) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(sonine/grid,ComputeSonineGrid)

#else

#ifndef SPARTA_COMPUTE_SONINE_GRID_H
#define SPARTA_COMPUTE_SONINE_GRID_H

#include "compute.h"

namespace SPARTA_NS {

class ComputeSonineGrid : public Compute {
 public:
  ComputeSonineGrid(class SPARTA *, int, char **);
  ~ComputeSonineGrid();
  void init();
  void compute_per_grid();
  double *normptr(int);
  bigint memory_usage();

 private:
  int imix,nvalue,ngroup,npergroup,ntotal;
  int *which,*moment,*order;

  double ***vave;
  double **sonine;

  int *group_norm_style;      // norm style of all values in Ith group
  double **norms;             // per-group ptr to norm vector by count or mass
  double **count;
};

}

#endif
#endif

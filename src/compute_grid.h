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

ComputeStyle(grid,ComputeGrid)

#else

#ifndef SPARTA_COMPUTE_GRID_H
#define SPARTA_COMPUTE_GRID_H

#include "compute.h"

namespace SPARTA_NS {

class ComputeGrid : public Compute {
 public:
  ComputeGrid(class SPARTA *, int, char **);
  ~ComputeGrid();
  void init();
  void compute_per_grid();
  int normflag(int);
  double *normptr(int);
  void reallocate();
  bigint memory_usage();

 private:
  int imix,nvalue,ngroup,ntotal;
  int *which;

  int nglocal;               // # of owned grid cells
  int *norm_style;           // norm style for each value
  double **norm_count;       // per-group ptr to norm vector, by count
  double **norm_mass;        // per-group ptr to norm vector, by mass
  double **norm_dof;         // per-group ptr to norm vector, by dof
};

}

#endif
#endif

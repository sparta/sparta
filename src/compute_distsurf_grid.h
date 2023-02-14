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

ComputeStyle(distsurf/grid,ComputeDistSurfGrid)

#else

#ifndef SPARTA_COMPUTE_DISTSURF_GRID_H
#define SPARTA_COMPUTE_DISTSURF_GRID_H

#include "compute.h"

namespace SPARTA_NS {

class ComputeDistSurfGrid : public Compute {
 public:
  ComputeDistSurfGrid(class SPARTA *, int, char **);
  ~ComputeDistSurfGrid();
  void init();
  void compute_per_grid();
  void reallocate();
  bigint memory_usage();

 protected:
  int nglocal,groupbit,sgroupbit;
  double sdir[3];
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

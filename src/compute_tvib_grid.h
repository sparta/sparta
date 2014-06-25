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

ComputeStyle(tvib/grid,ComputeTvibGrid)

#else

#ifndef SPARTA_COMPUTE_TVIB_GRID_H
#define SPARTA_COMPUTE_TVIB_GRID_H

#include "compute.h"

namespace SPARTA_NS {

class ComputeTvibGrid : public Compute {
 public:
  ComputeTvibGrid(class SPARTA *, int, char **);
  ~ComputeTvibGrid();
  void init();
  void compute_per_grid();
  void post_process_grid(double **, double **);
  void normwhich(int, int &, int &);
  double *normptr(int);
  void reallocate();
  bigint memory_usage();

 private:
  int imix,ngroup,nspecies;

  int nglocal;               // # of owned grid cells
  double **norm_count;       // per-group ptr to norm vector
  double *tspecies;          // per-species vibrational temperature
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Compute grid mixture ID does not exist

Self-explanatory.

E: Number of groups in compute grid mixture has changed

This mixture property cannot be changed after this compute command is
issued.

*/

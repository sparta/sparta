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
  void post_process_grid(void *, void *, int, int, double *, int);
  void normwhich(int, int &, int &);
  double *normptr(int);
  void reallocate();
  bigint memory_usage();

 private:
  int imix,nvalue,ngroup,npergroup,ntotal;
  int *which,*moment,*order;

  int nglocal;                  // # of owned grid cells
  double ***vcom;               // COM velocity per group and per cell
  double **masstot;             // total mass per group and per cell
  double **sonine;              // accumulator array

  int **value_norm_style;       // I,J = norm style of Jth value in Ith group
  double **norm_count;          // per-group ptr to norm vector, by count
  double **norm_mass;           // per-group ptr to norm vector, by mass
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Compute sonine/grid mixture ID does not exist

Self-explanatory.

E: Number of groups in compute sonine/grid mixture has changed

This mixture property cannot be changed after this compute command is
issued.

E: Invalid call to ComputeSonineGrid::post_process_grid()

This indicates a coding error.  Please report the issue to the SPARTA
developers.

*/

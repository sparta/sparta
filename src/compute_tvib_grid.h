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
  int query_tally_grid(int, double **&, int *&);
  double post_process_grid(int, int, int, double **, int *, double *, int);
  void reallocate();
  bigint memory_usage();

 private:
  int imix,ngroup,mixspecies,nspecies;

  int ntotal;                // total # of columns in tally array
  int nglocal;               // # of owned grid cells

  int *nmap;                 // # of tally quantities each group value uses
  int **map;                 // which tally columns each group value uses
  double **tally;            // array of tally quantities, cells by ntotal

  double *tspecies;          // per-species vibrational temperature
  int *s2t;                  // convert particle species to tally column
  int *t2s;                  // convert tally column to particle species
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

E: Number of groups in compute tvib/grid mixture has changed

This mixture property cannot be changed after this compute command is
issued.

E: Number of species in compute tvib/grid mixture has changed

This mixture property cannot be changed after this compute command is
issued.

*/

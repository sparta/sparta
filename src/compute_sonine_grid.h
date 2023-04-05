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
  int query_tally_grid(int, double **&, int *&);
  void post_process_grid(int, int, double **, int *, double *, int);
  void reallocate();
  bigint memory_usage();

 protected:
  int groupbit,imix,ngroup;

  int nvalue;                // # of user-requested values
                             // each may generate multiple output values

  int *which;                // keyword for each user requested value
  int *moment;               // moment for each user requested value
  int *order;                // order (1-5) for each user requested value

  int noutpergroup;          // # of output values per group
  int npergroup;             // # of unique tally quantities per group
  int ntotal;                // total # of columns in tally array
  int nglocal;               // # of owned grid cells

  int nmap;                  // # of tally quantities every output value uses
  int **map;                 // which tally columns each output value uses
  double **tally;            // array of tally quantities, cells by ntotal

  double ***vcom;            // COM velocity and mass per group and per cell
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

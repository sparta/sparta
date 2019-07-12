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

ComputeStyle(isurf/grid,ComputeISurfGrid)

#else

#ifndef SPARTA_COMPUTE_ISURF_GRID_H
#define SPARTA_COMPUTE_ISURF_GRID_H

#include "compute.h"
#include "grid.h"
#include "surf.h"

namespace SPARTA_NS {

class ComputeISurfGrid : public Compute {
 public:
  ComputeISurfGrid(class SPARTA *, int, char **);
  ~ComputeISurfGrid();
  void init();
  void compute_per_grid();
  double post_process_grid(int, int, int, double **, int *, double *, int);
  virtual void clear();
  void surf_tally(int, int, Particle::OnePart *, 
                  Particle::OnePart *, Particle::OnePart *);
  void reallocate();
  bigint memory_usage();

 protected:
  int groupbit,imix,nvalue,ngroup;
  int collapsed;
  int maxsend;

  int *proclist;
  double *sbuf;

  int *which;              // keyword for each user requested value

  int nsurf;               // # of lines/tris I store
                           // surf->nlocal+nghost for explicit all or distributed

  int npergroup;           // # of unique tally quantities per group
  int ntotal;              // total # of columns in tally array
  int ngtotal;             // # of owned + ghost grid cells

  int dimension;           // local copies
  Grid::ChildCell *cells;
  Grid::ChildInfo *cinfo;
  Grid::SplitInfo *sinfo;
  Surf::Line *lines;
  Surf::Tri *tris;

  int weightflag;          // 1 if cell weighting is enabled
  double weight;           // particle weight, based on initial cell
  double nfactor;          // dt/fnum for normalization
  double nfactor_inverse;  // fnum/dt for normalization
  double *normflux;        // normalization factor for each surf element
  double nfactor_previous; // nfactor from previous run

  void grow_tally();
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

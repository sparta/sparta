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

ComputeStyle(property/grid,ComputePropertyGrid)

#else

#ifndef SPARTA_COMPUTE_PROPERTY_GRID_H
#define SPARTA_COMPUTE_PROPERTY_GRID_H

#include "compute.h"

namespace SPARTA_NS {

class ComputePropertyGrid : public Compute {
 public:
  ComputePropertyGrid(class SPARTA *, int, char **);
  ~ComputePropertyGrid();
  void init();
  void compute_per_grid();
  void reallocate();
  bigint memory_usage();

 protected:
  int groupbit,nvalues,nglocal;
  int *index;
  double *buf;

  typedef void (ComputePropertyGrid::*FnPtrPack)(int);
  FnPtrPack *pack_choice;              // ptrs to pack functions

  void pack_id(int);
  void pack_proc(int);
  void pack_xlo(int);
  void pack_ylo(int);
  void pack_zlo(int);
  void pack_xhi(int);
  void pack_yhi(int);
  void pack_zhi(int);
  void pack_xc(int);
  void pack_yc(int);
  void pack_zc(int);
  void pack_vol(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Invalid compute property/grid field for 2d simulation

Fields that reference z-dimension properties cannot be used
in a 2d simulation.

E: Invalid keyword in compute property/grid command

Self-explanatory.

*/

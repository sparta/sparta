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

#ifdef FIX_CLASS

FixStyle(grid/check,FixGridCheck)

#else

#ifndef SPARTA_FIX_GRID_CHECK_H
#define SPARTA_FIX_GRID_CHECK_H

#include "fix.h"

namespace SPARTA_NS {

class FixGridCheck : public Fix {
 public:
  FixGridCheck(class SPARTA *, int, char **);
  virtual ~FixGridCheck() {}
  int setmask();
  void init();
  void setup();
  virtual void end_of_step();
  double compute_scalar();

 protected:
  int ntotal,outflag;
  int outside_check;
  int dim;

  class Cut2d *cut2d;
  class Cut3d *cut3d;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Particle %d,%d on proc %d is in invalid cell  on timestep %ld

The particle is in a cell indexed by a value that is out-of-bounds for
the cells owned by this processor.

E: Particle %d,%d on proc %d is in split cell  on timestep %ld

This should not happend.  The particle should be in one
of the sub-cells of the split cell.

E: Particle %d,%d on proc %d is outside cell  on timestep %ld

The particle's coordinates are not within the grid cell
it is supposed to be in.

W: %d particles were in wrong cells on timestep %ld

This is the total number of particles that are incorrectly
matched to their grid cell.

*/

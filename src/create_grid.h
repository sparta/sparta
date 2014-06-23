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

#ifdef COMMAND_CLASS

CommandStyle(create_grid,CreateGrid)

#else

#ifndef SPARTA_CREATE_GRID_H
#define SPARTA_CREATE_GRID_H

#include "pointers.h"

namespace SPARTA_NS {

class CreateGrid : protected Pointers {
 public:
  CreateGrid(class SPARTA *);
  void command(int, char **);

 private:
  void bounds(char *, int, int &, int &);
  void procs2grid(int, int, int, int &, int &, int &);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot create grid before simulation box is defined

Self-explanatory.

E: Cannot create grid when grid is already defined

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Create_grid nz value must be 1 for a 2d simulation

Self-explanatory.

E: Cannot use specified create_grid options with more than one level

When defining a grid with more than one level, the other create_grid
keywords (stride, clump, block, etc) cannot be used.  The child grid
cells will be assigned to processors in round-robin order as explained
on the create_grid doc page.

E: Bad grid of processors for create_grid

For block style, product of Px,Py,Pz must equal total number of
processors.

E: Numeric index is out of bounds

A command with an argument that specifies an integer or range of
integers is using a value that is less than 1 or greater than the
maximum allowed limit.

*/

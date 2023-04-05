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

#ifdef COMMAND_CLASS

CommandStyle(balance_grid,BalanceGrid)

#else

#ifndef SPARTA_BALANCE_GRID_H
#define SPARTA_BALANCE_GRID_H

#include "pointers.h"

namespace SPARTA_NS {

class BalanceGrid : protected Pointers {
 public:
  BalanceGrid(class SPARTA *);
  void command(int, char **, int outflag=1);

 private:
  double last;

  void procs2grid(int, int, int, int &, int &, int &);
  void timer_cell_weights(double *&);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot balance grid before grid is defined

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Invalid balance_grid style for non-uniform grid

Some balance styles can only be used when the grid is uniform.  See
the command doc page for details.

E: Bad grid of processors for balance_grid block

Product of Px,Py,Pz must equal total number of processors.

*/

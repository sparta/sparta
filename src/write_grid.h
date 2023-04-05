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

CommandStyle(write_grid,WriteGrid)

#else

#ifndef SPARTA_WRITE_GRID_H
#define SPARTA_WRITE_GRID_H

#include "stdio.h"
#include "pointers.h"

namespace SPARTA_NS {

class WriteGrid : protected Pointers {
 public:
  int silent;

  WriteGrid(class SPARTA *);
  void command(int, char **);

 private:
  FILE *fp;

  void header();
  void write();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot write grid when grid is not defined

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Cannot open file %s

The specified file cannot be opened.  Check that the path and name are
correct. If the file is a compressed file, also check that the gzip
executable can be found and run.

*/

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

CommandStyle(read_grid,ReadGrid)

#else

#ifndef SPARTA_READ_GRID_H
#define SPARTA_READ_GRID_H

#include "stdio.h"
#include "pointers.h"

namespace SPARTA_NS {

class ReadGrid : protected Pointers {
 public:
  ReadGrid(class SPARTA *);
  ~ReadGrid();
  void command(int, char **);

 private:
  int me;
  char *line,*keyword,*buffer;
  FILE *fp;
  int compressed;

  int nparents;

  void create_parents(int, char *);
  void create_children();
  void open(char *);
  void header();
  void parse_keyword(int);
  int count_words(char *);
  void read_parents();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot read grid before simulation box is defined

UNDOCUMENTED

E: Cannot read grid when grid is already defined

UNDOCUMENTED

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Grid file cannot parse Parents section

UNDOCUMENTED

E: Incorrect parents format in grid file

UNDOCUMENTED

E: Invalid cell ID in read_grid file

UNDOCUMENTED

E: Duplicate cell ID in read_grid

UNDOCUMENTED

E: Parent cell's parent does not exist in read_grid

UNDOCUMENTED

E: 

UNDOCUMENTED

E: Nz value in read_grid file must be 1 for a 2d simulation

UNDOCUMENTED

E: Cannot open gzipped file

SPARTA was compiled without support for reading and writing gzipped
files through a pipeline to the gzip program with -DSPARTA_GZIP.

E: Cannot open file %s

The specified file cannot be opened.  Check that the path and name are
correct. If the file is a compressed file, also check that the gzip
executable can be found and run.

E: Unexpected end of grid file

UNDOCUMENTED

E: Grid file does not contain parents

UNDOCUMENTED

*/

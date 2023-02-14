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
  void read(char *, int);

 private:
  int me,nprocs;
  char *line,*keyword,*buffer;
  FILE *fp;
  int compressed;

  bigint ncell;
  int nlevels;
  int whichproc;

  // grid level info

  struct Level {
    int setflag;                         // setflag = 1 if specified
    int cx,cy,cz;                        // grid of child cells at this level
  };

  void read_cells();
  void create_cells(int, char *);
  void open(char *);
  void header();
  void parse_keyword(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot read grid before simulation box is defined

Self-explanatory.

E: Cannot read grid when grid is already defined

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Read_grid did not find parents section of grid file

Expected Parents section but did not find keyword.

E: Incorrect format of parent cell in grid file

Number of words in a parent cell line was not the expected number.

E: Invalid cell ID in grid file

A cell ID could not be converted into numeric format.

E: Duplicate cell ID in grid file

Parent cell IDs must be unique.

E: Parent cell's parent does not exist in grid file

Parent cells must be listed in order such that
each cell's parents have already appeared in the list.

E: Invalid Nx,Ny,Nz values in grid file

A Nx or Ny or Nz value for a parent cell is <= 0.

E: Nz value in read_grid file must be 1 for a 2d simulation

Self-explanatory.

E: Cannot open gzipped file

SPARTA was compiled without support for reading and writing gzipped
files through a pipeline to the gzip program with -DSPARTA_GZIP.

E: Cannot open file %s

The specified file cannot be opened.  Check that the path and name are
correct. If the file is a compressed file, also check that the gzip
executable can be found and run.

E: Unexpected end of grid file

Self-explanatory.

E: Grid file does not contain parents

No parent cells appeared in the grid file.

*/

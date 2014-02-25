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

*/

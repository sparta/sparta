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

CommandStyle(write_surf,WriteSurf)

#else

#ifndef SPARTA_WRITE_SURF_H
#define SPARTA_WRITE_SURF_H

#include "stdio.h"
#include "pointers.h"

namespace SPARTA_NS {

class WriteSurf : protected Pointers {
 public:
  WriteSurf(class SPARTA *);
  void command(int, char **);
  void write_file(FILE *, int);

 private:
  int me,nprocs;

  struct SurfIDType {
    surfint id;
    int type;
  };

  void write_file_all_points(FILE *);
  void write_file_all_nopoints(FILE *);
  void write_file_distributed_points(FILE *);
  void write_file_distributed_nopoints(FILE *);
  void write_file_implicit_points(FILE *);
  void write_file_implicit_nopoints(FILE *);
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

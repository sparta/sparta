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

CommandStyle(write_isurf,WriteISurf)

#else

#ifndef SPARTA_WRITE_ISURF_H
#define SPARTA_WRITE_ISURF_H

#include "stdio.h"
#include "pointers.h"

namespace SPARTA_NS {

class WriteISurf : protected Pointers {
 public:
  WriteISurf(class SPARTA *);
  void command(int, char **);

 private:
  int me,nprocs;
  int dim;
  int ggroup,groupbit;
  int nx,ny,nz,precision;
  int ncorner;

  class FixAblate *ablate;
  double *dbuf,*dbufall;

  void collect_values();
  void write_file(FILE *);
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

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

CommandStyle(write_surf,WriteSurf)

#else

#ifndef SPARTA_WRITE_SURF_H
#define SPARTA_WRITE_SURF_H

#include "stdio.h"
#include "pointers.h"

namespace SPARTA_NS {

class WriteSurf : protected Pointers {
 public:
  int statflag;

  WriteSurf(class SPARTA *);
  void command(int, char **);

 private:
  int me,nprocs;
  int dim;
  FILE *fp;

  int pointflag;             // 1/0 to include/exclude Points section in file
  int multiproc;             // 0 = proc 0 writes for all
                             // else # of procs writing files
  int filewriter;            // 1 if this proc writes to file, else 0
  int icluster;              // which cluster I am in
  int nclusterprocs;         // # of procs in my cluster that write to one file
  int fileproc;              // ID of proc in my cluster who writes to file

  struct SurfIDType {
    surfint id;
    int type;
  };

  void write_file(char *);
  void write_file_all_points(char *);
  void write_file_all_nopoints(char *);
  void write_file_distributed_points(char *);
  void write_file_distributed_nopoints(char *);

  void write_base(char *);
  void open(char *);
  void write_header(int);
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

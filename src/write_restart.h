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

CommandStyle(write_restart,WriteRestart)

#else

#ifndef SPARTA_WRITE_RESTART_H
#define SPARTA_WRITE_RESTART_H

#include "stdio.h"
#include "pointers.h"

namespace SPARTA_NS {

class WriteRestart : protected Pointers {
 public:
  WriteRestart(class SPARTA *);
  void command(int, char **);
  void multiproc_options(int, int, char **);
  void write(char *);
  void write_less_memory(char *);

 private:
  int me,nprocs;
  FILE *fp;

  int multiproc;             // 0 = proc 0 writes for all
                             // else # of procs writing files
  int nclusterprocs;         // # of procs in my cluster that write to one file
  int filewriter;            // 1 if this proc writes a file, else 0
  int fileproc;              // ID of proc in my cluster who writes to file
  int icluster;              // which cluster I am in

  void header();
  void box_params();
  void particle_params();
  void grid_params();
  void surf_params();
  void file_layout(int);

  void magic_string();
  void endian();
  void version_numeric();

  void write_int(int, int);
  void write_bigint(int, bigint);
  void write_double(int, double);
  void write_string(int, char *);
  void write_int_vec(int, int, int *);
  void write_double_vec(int, int, double *);
  void write_char_vec(int, bigint, char *);
  void write_char_vec(int, bigint, int, char *);
  void write_char_vec(int, char *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot write restart file before grid is defined

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Cannot (yet) use global mem/limit without % in restart file name

This feature is not yet implemented.

E: Cannot use write_restart fileper without % in restart file name

Self-explanatory.

E: Cannot use write_restart nfile without % in restart file name

Self-explanatory.

E: Cannot open restart file %s

The specified file cannot be opened.  Check that the path and name are
correct.  If the file is a compressed file, also check that the gzip
executable can be found and run.

*/

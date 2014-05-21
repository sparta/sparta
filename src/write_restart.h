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
};

}

#endif
#endif

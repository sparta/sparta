/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   www/sandia.gov/sparta.html
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2012) Sandia Corporation.  Under the terms of Contract
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
  void write(char *);

 private:
  int me,nprocs;
  FILE *fp;

  void write_int(int, int);
  void write_double(int, double);
  void write_char(int, char *);
  void write_bigint(int, bigint);
};

}

#endif
#endif

/* ----------------------------------------------------------------------
   DSMC - Sandia parallel DSMC code
   www.sandia.gov/~sjplimp/dsmc.html
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2011) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level DSMC directory.
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS

CommandStyle(write_restart,WriteRestart)

#else

#ifndef DSMC_WRITE_RESTART_H
#define DSMC_WRITE_RESTART_H

#include "stdio.h"
#include "pointers.h"

namespace DSMC_NS {

class WriteRestart : protected Pointers {
 public:
  WriteRestart(class DSMC *);
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

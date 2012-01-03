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

#ifndef DSMC_MODIFY_H
#define DSMC_MODIFY_H

#include "pointers.h"

namespace DSMC_NS {

class Modify : protected Pointers {
 public:
  int nfix,maxfix;
  int n_start_of_step,n_end_of_step;

  class Fix **fix;           // list of fixes
  int *fmask;                // bit mask for when each fix is applied

  int ncompute,maxcompute;   // list of computes
  class Compute **compute;

  Modify(class DSMC *);
  ~Modify();
  void init();
  void start_of_step();
  void end_of_step();

  void add_fix(int, char **, char *suffix = NULL);
  void delete_fix(const char *);
  int find_fix(const char *);

  void add_compute(int, char **, char *suffix = NULL);
  void delete_compute(const char *);
  int find_compute(const char *);
  void clearstep_compute();

  bigint memory_usage();

 protected:

  // lists of fixes to apply at different stages of timestep

  int *list_start_of_step,*list_end_of_step;

  void list_init(int, int &, int *&);
};

}

#endif

/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2012) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifndef SPARTA_MODIFY_H
#define SPARTA_MODIFY_H

#include "pointers.h"

namespace SPARTA_NS {

class Modify : protected Pointers {
 public:
  int nfix,maxfix;
  int n_start_of_step,n_end_of_step;

  class Fix **fix;           // list of fixes
  int *fmask;                // bit mask for when each fix is applied

  int ncompute,maxcompute;   // list of computes
  class Compute **compute;

  Modify(class SPARTA *);
  ~Modify();
  void init();
  void start_of_step();
  void end_of_step();

  void add_fix(int, char **);
  void delete_fix(const char *);
  int find_fix(const char *);

  void add_compute(int, char **);
  void delete_compute(const char *);
  int find_compute(const char *);

  void clearstep_compute();
  void addstep_compute(bigint);
  void addstep_compute_all(bigint);

  bigint memory_usage();

 protected:

  // lists of fixes to apply at different stages of timestep

  int *list_start_of_step,*list_end_of_step;

  int n_timeflag;            // list of computes that store time invocation
  int *list_timeflag;

  void list_init(int, int &, int *&);
  void list_init_compute();
};

}

#endif

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

#ifndef SPARTA_MODIFY_H
#define SPARTA_MODIFY_H

#include "pointers.h"
#include "particle.h"

namespace SPARTA_NS {

class Modify : protected Pointers {
 public:
  int nfix,maxfix;
  int n_start_of_step,n_end_of_step;
  int n_pergrid,n_update_custom,n_gas_react,n_surf_react;

  class Fix **fix;           // list of fixes
  int *fmask;                // bit mask for when each fix is applied

  int ncompute,maxcompute;   // list of computes
  class Compute **compute;

  Modify(class SPARTA *);
  ~Modify();
  void init();
  void setup();
  virtual void start_of_step();
  virtual void end_of_step();

  virtual int pack_grid_one(int, char *, int);
  virtual int unpack_grid_one(int, char *);
  virtual void copy_grid_one(int, int);
  virtual void add_grid_one();
  virtual void reset_grid_count(int);
  virtual void grid_changed();

  void add_fix(int, char **);
  void delete_fix(const char *);
  int find_fix(const char *);

  void add_compute(int, char **);
  void delete_compute(const char *);
  int find_compute(const char *);

  void clearstep_compute();
  void addstep_compute(bigint);
  void addstep_compute_all(bigint);

  void list_init_fixes();
  void list_init_computes();

  virtual void update_custom(int, double, double, double, double *);
  virtual void gas_react(int);
  virtual void surf_react(Particle::OnePart *, int &, int &);

  bigint memory_usage();

 protected:

  // lists of fixes to apply at different stages of timestep

  int *list_start_of_step,*list_end_of_step;

  int *end_of_step_every;

  int *list_pergrid;         // list of fixes that store per grid cell info
  int *list_update_custom;    // list of fixes with update_custom() method
  int *list_gas_react;       // list of fixes with gas_react() method
  int *list_surf_react;      // list of fixes with surf_react() method

  int n_timeflag;            // list of computes that store time invocation
  int *list_timeflag;

  void list_init(int, int &, int *&);
  void list_init_end_of_step(int, int &, int *&);
};

}

#endif

/* ERROR/WARNING messages:

E: Fix command before simulation box is defined

The fix command cannot be used before a read_data, read_restart, or
create_box command.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Replacing a fix, but new style != old style

A fix ID can be used a 2nd time, but only if the style matches the
previous fix.  In this case it is assumed you with to reset a fix's
parameters.  This error may mean you are mistakenly re-using a fix ID
when you do not intend to.

E: Invalid fix style

The choice of fix style is unknown.

E: Could not find fix ID to delete

Self-explanatory.

E: Reuse of compute ID

A compute ID cannot be used twice.

E: Invalid compute style

Self-explanatory.

E: Could not find compute ID to delete

Self-explanatory.

*/

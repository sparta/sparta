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

#ifndef SPARTA_STATS_H
#define SPARTA_STATS_H

#include "pointers.h"

namespace SPARTA_NS {

class Stats : protected Pointers {
 public:
  Stats(class SPARTA *);
  ~Stats();
  void init();
  void modify_params(int, char **);
  void set_fields(int, char **);
  void header();
  void compute(int);
  int evaluate_keyword(char *, double *);

 private:
  char *line;
  char **keyword;
  int *vtype;

  int nfield;
  int me;

  char **format;
  char *format_line_user;
  char *format_float_user,*format_int_user,*format_bigint_user;
  char **format_column_user;

  char *format_float_def,*format_int_def;
  char format_bigint_def[8];

  int firststep;
  int flushflag,lineflag;

  double wall0;
  double last_tpcpu,last_spcpu;
  double last_time;
  bigint last_step;
                         // data used by routines that compute single values
  int ivalue;            // integer value to print
  double dvalue;         // double value to print
  bigint bivalue;        // big integer value to print
  int ifield;            // which field in thermo output is being computed
  int *field2index;      // which compute,fix,variable calcs this field
  int *argindex1;        // indices into compute,fix scalar,vector
  int *argindex2;

  int ncompute;                // # of Compute objects called by stats
  char **id_compute;           // their IDs
  int *compute_which;          // 0/1/2 if should call scalar,vector,array
  class Compute **computes;    // list of ptrs to the Compute objects

  int nfix;                    // # of Fix objects called by stats
  char **id_fix;               // their IDs
  class Fix **fixes;           // list of ptrs to the Fix objects

  int nsurfcollide;            // # of SurfCollide objs called by stats
  char **id_surf_collide;      // their IDs
  class SurfCollide **sc;      // list of ptrs to SurfCollide objects

  int nsurfreact;             // # of SurfReact objects called by stats
  char **id_surf_react;       // their IDs
  class SurfReact **sr;       // list of ptrs to SurfReact objects

  int nvariable;               // # of variables evaulated by stats
  char **id_variable;          // list of variable names
  int *variables;              // list of Variable indices

  // private methods

  void allocate();
  void deallocate();

  int add_compute(const char *, int);
  int add_fix(const char *);
  int add_surf_collide(const char *);
  int add_surf_react(const char *);
  int add_variable(const char *);

  typedef void (Stats::*FnPtr)();
  void addfield(const char *, FnPtr, int);
  FnPtr *vfunc;                // list of ptrs to functions

  void compute_compute();        // functions that compute a single value
  void compute_fix();            // via calls to Compute,Fix,
  void compute_surf_collide();   //   SurfCollide,SurfReact,Variable classes
  void compute_surf_react();
  void compute_variable();

  // functions that compute a single value
  // customize a new keyword by adding a method prototype

  void compute_step();
  void compute_elapsed();
  void compute_elaplong();
  void compute_dt();
  void compute_time();
  void compute_cpu();
  void compute_tpcpu();
  void compute_spcpu();
  void compute_wall();

  void compute_np();
  void compute_ntouch();
  void compute_ncomm();
  void compute_nbound();
  void compute_nexit();
  void compute_nscoll();
  void compute_nscheck();
  void compute_ncoll();
  void compute_nattempt();
  void compute_nreact();
  void compute_nsreact();

  void compute_npave();
  void compute_ntouchave();
  void compute_ncommave();
  void compute_nboundave();
  void compute_nexitave();
  void compute_nscollave();
  void compute_nscheckave();
  void compute_ncollave();
  void compute_nattemptave();
  void compute_nreactave();
  void compute_nsreactave();

  void compute_ngrid();
  void compute_nsplit();
  void compute_maxlevel();

  void compute_vol();
  void compute_lx();
  void compute_ly();
  void compute_lz();

  void compute_xlo();
  void compute_xhi();
  void compute_ylo();
  void compute_yhi();
  void compute_zlo();
  void compute_zhi();
};

}

#endif

/* ERROR/WARNING messages:

E: Could not find stats compute ID

Compute ID specified in stats_style command does not exist.

E: Could not find stats fix ID

Fix ID specified in stats_style command does not exist.

E: Stats and fix not computed at compatible times

Fixes generate values on specific timesteps.  The stats output
does not match these timesteps.

E: Could not find stats variable name

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Stats_modify int format does not contain d character

Self-explanatory.

E: Stats compute does not compute scalar

Self-explanatory.

E: Stats compute does not compute vector

Self-explanatory.

E: Stats compute vector is accessed out-of-range

Self-explanatory.

E: Stats compute does not compute array

Self-explanatory.

E: Stats compute array is accessed out-of-range

Self-explanatory.

E: Stats fix does not compute scalar

Self-explanatory.

E: Stats fix does not compute vector

Self-explanatory.

E: Stats fix vector is accessed out-of-range

Self-explanatory.

E: Stats fix does not compute array

Self-explanatory.

E: Stats fix array is accessed out-of-range

Self-explanatory.

E: Stats variable is not equal-style variable

Only equal-style variables can be output with stats output, not
particle-style or grid-style or surf-style variables.

E: Stats variable cannot be indexed

A variable used as a stats keyword cannot be indexed.
E.g. v_foo must be used, not v_foo[100].

E: Invalid keyword in stats_style command

One or more specified keywords are not recognized.

E: Variable stats keyword cannot be used between runs

Stats keywords that refer to time (such as cpu, elapsed) do not make
sense in between runs.

*/

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

#ifndef SPARTA_OUTPUT_H
#define SPARTA_OUTPUT_H

#include "pointers.h"

namespace SPARTA_NS {

class Output : protected Pointers {
 public:
  bigint next;                 // next timestep for any kind of output

  bigint next_stats;           // next timestep for stats output
  int stats_every;             // stats output every this many steps
  bigint last_stats;           // last timestep stats was output
  char *var_stats;             // variable name for stats frequency
  int ivar_stats;              // variable index for stats frequency
  class Stats *stats;          // statistical output

  int ndump;                   // # of Dumps defined
  int max_dump;                // max size of Dump list
  bigint next_dump_any;        // next timestep for any Dump
  int *every_dump;             // output of each Dump every this many steps
  bigint *next_dump;           // next timestep to do each Dump
  bigint *last_dump;           // last timestep each snapshot was output
  char **var_dump;             // variable name for dump frequency
  int *ivar_dump;              // variable index for dump frequency
  class Dump **dump;           // list of defined Dumps

  int restart_flag;            // 1 if any restart files are written
  int restart_flag_single;     // 1 if single restart files are written
  int restart_flag_double;     // 1 if double restart files are written
  bigint next_restart;         // next timestep to write any restart file
  bigint next_restart_single;  // next timestep to write a single restart file
  bigint next_restart_double;  // next timestep to write a double restart file
  int restart_every_single;    // single restart file write freq, 0 if var
  int restart_every_double;    // double restart file write freq, 0 if var
  bigint last_restart;         // last timestep any restart file was output
  int restart_toggle;          // 0 if use restart2a as prefix, 1 if restart2b
  char *var_restart_single;    // variable name for single restart freq
  char *var_restart_double;    // variable name for double restart freq
  int ivar_restart_single;     // index of var_restart_single
  int ivar_restart_double;     // index of var_restart_double
  char *restart1;              // name single restart file
  char *restart2a,*restart2b;  // names of double restart files
  class WriteRestart *restart; // class for writing restart files

  Output(class SPARTA *);
  ~Output();
  void init();
  void setup(int);                   // initial output before run/min
  void write(bigint);                // output for current timestep
  void write_dump(bigint);           // force output of dump snapshots
  void write_restart(bigint);        // force output of a restart file
  void reset_timestep(bigint);       // reset next timestep for all output

  void add_dump(int, char **);       // add a Dump to Dump list
  void modify_dump(int, char **);    // modify a Dump
  void delete_dump(char *);          // delete a Dump from Dump list

  void set_stats(int, char **);      // set stats output frequency
  void create_stats(int, char **);   // create a Stats style
  void create_restart(int, char **); // create Restart and restart files

  void memory_usage();               // print out memory usage
};

}

#endif

/* ERROR/WARNING messages:

E: Variable name for stats every does not exist

Self-explanatory.

E: Variable for stats every is invalid style

It must be an equal-style variable.

E: Variable name for dump every does not exist

Self-explanatory.

E: Variable for dump every is invalid style

Only equal-style variables can be used.

E: Variable name for restart does not exist

Self-explanatory.

E: Variable for restart is invalid style

It must be an equal-style variable.

E: Dump every variable returned a bad timestep

The variable must return a timestep greater than the current timestep.

E: Restart variable returned a bad timestep

The variable must return a timestep greater than the current timestep.

E: Stats every variable returned a bad timestep

The variable must return a timestep greater than the current timestep.

E: Stats_modify every variable returned a bad timestep

The variable must return a timestep greater than the current timestep.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Reuse of dump ID

A dump ID cannot be used twice.

E: Invalid dump frequency

Dump frequency must be 1 or greater.

E: Invalid dump style

The choice of dump style is unknown.

E: Cound not find dump_modify ID

Self-explanatory.

E: Could not find undump ID

A dump ID used in the undump command does not exist.

E: Both restart files must use % or neither

Self-explanatory.

*/

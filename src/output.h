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

  bigint next_restart;         // next timestep to write a restart file
  int restart_every;           // write a restart file every this many steps
  bigint last_restart;         // last timestep a restart file was output
  int restart_toggle;          // 0 if use restart1 as prefix
                               // 1 if use restart1 as file, 2 for restart2
  char *restart1,*restart2;    // names for restart files
  class WriteRestart *restart; // Restart output

  Output(class SPARTA *);
  ~Output();
  void init();
  void setup(int);                   // initial output before run/min
  void write(bigint);                // output for current timestep
  void write_dump(bigint);           // force output of dump snapshots
  void write_restart(bigint);        // force output of a restart file

  void add_dump(int, char **);       // add a Dump to Dump list
  void modify_dump(int, char **);    // modify a Dump
  void delete_dump(char *);          // delete a Dump from Dump list

  void create_stats(int, char **);   // create a Stats style
  void create_restart(int, char **); // create Restart and restart files

  void memory_usage();               // print out memory usage
};

}

#endif

/* ERROR/WARNING messages:

E: Variable name for stats every does not exist

UNDOCUMENTED

E: Variable for stats every is invalid style

UNDOCUMENTED

E: Variable name for dump every does not exist

Self-explanatory.

E: Variable for dump every is invalid style

Only equal-style variables can be used.

E: Dump every variable returned a bad timestep

The variable must return a timestep greater than the current timestep.

E: stats every variable returned a bad timestep

UNDOCUMENTED

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

*/

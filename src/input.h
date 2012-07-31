/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   www.sandia.gov/sparta.html
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2012) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifndef SPARTA_INPUT_H
#define SPARTA_INPUT_H

#include "stdio.h"
#include "pointers.h"

namespace SPARTA_NS {

class Input : protected Pointers {
 public:
  int narg;                    // # of command args
  char **arg;                  // parsed args for command
  class Variable *variable;    // defined variables

  Input(class SPARTA *, int, char **);
  ~Input();
  void file();                   // process all input
  void file(const char *);       // process an input script
  char *one(const char *);       // process a single command
  void substitute(char *, int);  // substitute for variables in a string

 private:
  int me;                      // proc ID
  char *command;               // ptr to current command
  int maxarg;                  // max # of args in arg
  char *line,*copy,*work;      // input line & copy of it
  int echo_screen;             // 0 = no, 1 = yes
  int echo_log;                // 0 = no, 1 = yes
  int nfile,maxfile;           // current # and max # of open input files
  int label_active;            // 0 = no label, 1 = looking for label
  char *labelstr;              // label string being looked for
  int jump_skip;               // 1 if skipping next jump, 0 otherwise

  FILE **infiles;              // list of open input files

  void parse();                // parse an input text line
  int execute_command();       // execute a single command

  void clear();                // input script commands
  void echo();
  void ifthenelse();
  void include();
  void jump();
  void label();
  void log();
  void next_command();
  void partition();
  void print();
  void shell();
  void variable_command();

  // SPARTA commands

  void boundary();
  void bound_modify();
  void collide_command();
  void compute();
  void dimension();
  void dump();
  void dump_modify();
  void fix();
  void global();
  void mixture();
  void seed();
  void species();
  void stats();
  void stats_modify();
  void stats_style();
  void surf_collide();
  void timestep();
  void uncompute();
  void undump();
  void unfix();
  void units();
};

}

#endif

/* ERROR/WARNING messages:

E: Label wasn't found in input script

Self-explanatory.

E: Input line too long: %s

This is a hard (very large) limit defined in the input.cpp file.

E: Unknown command: %s

The command is not known to SPARTA.  Check the input script.

E: Another input script is already being processed

Cannot attempt to open a 2nd input script, when the original file is
still being processed.

E: Cannot open input script %s

Self-explanatory.

E: Unbalanced quotes in input line

No matching end double quote was found following a leading double
quote.

E: Invalid variable name

Variable name used in an input script line is invalid.

E: Substitution for illegal variable

Input script line contained a variable that could not be substituted
for.

E: Input line too long after variable substitution

This is a hard (very large) limit defined in the input.cpp file.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Cannot open logfile %s

The SPARTA log file specified in the input script cannot be opened.
Check that the path and name are correct.

E: Partition numeric index is out of bounds

UNDOCUMENTED

E: Invalid collide style

UNDOCUMENTED

E: Dimension command after simulation box is defined

The dimension command cannot be used after a read_data,
read_restart, or create_box command.

*/

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
  void substitute(char *&, char *&, int &, int &, int);
                                 // substitute for variables in a string
  int expand_args(int, char **, int, char **&);  // expand args due to wildcard

  double numeric(const char *, int, const char *);    // arg checking
  int inumeric(const char *, int, char *);
  bigint bnumeric(const char *, int, char *);
  void bounds(char *, int, int &, int &, int nmin=1);
  int count_words(char *);

 private:
  int me;                      // proc ID
  char *command;               // ptr to current command
  int maxarg;                  // max # of args in arg
  char *line,*copy,*work;      // input line & copy and work string
  int maxline,maxcopy,maxwork; // max lengths of char strings
  int echo_screen;             // 0 = no, 1 = yes
  int echo_log;                // 0 = no, 1 = yes
  int nfile,maxfile;           // current # and max # of open input files
  int label_active;            // 0 = no label, 1 = looking for label
  char *labelstr;              // label string being looked for
  int jump_skip;               // 1 if skipping next jump, 0 otherwise
  int ifthenelse_flag;         // 1 if executing commands inside an if-then-else

  FILE **infiles;              // list of open input files

  void parse();                          // parse an input text line
  char *nextword(char *, char **);       // find next word in string with quotes
  void reallocate(char *&, int &, int);  // reallocate a char string
  int execute_command();                 // execute a single command

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
  void quit();
  void shell();
  void variable_command();

  // SPARTA commands

  void boundary();
  void bound_modify();
  void collide_command();
  void collide_modify();
  void compute();
  void dimension();
  void dump();
  void dump_modify();
  void fix();
  void global();
  void group();
  void mixture();
  void package();
  void react_command();
  void react_modify();
  void region();
  void reset_timestep();
  void restart();
  void seed();
  void species();
  void species_modify();
  void stats();
  void stats_modify();
  void stats_style();
  void surf_collide();
  void surf_modify();
  void surf_react();
  void timestep();
  void uncompute();
  void undump();
  void unfix();
  void units();
  void weight();
};

}

#endif

/* ERROR/WARNING messages:

E: Label wasn't found in input script

Self-explanatory.

E: Unknown command: %s

The command is not known to SPARTA.  Check the input script.

E: Invalid use of library file() function

This function is called thru the library interface.  This
error should not occur.  Contact the developers if it does.

E: Cannot open input script %s

Self-explanatory.

E: Unbalanced quotes in input line

No matching end double quote was found following a leading double
quote.

E: Input line quote not followed by whitespace

An end quote must be followed by whitespace.

E: Invalid variable name

Variable name used in an input script line is invalid.

E: Invalid immediate variable

Syntax of immediate value is incorrect.

E: Substitution for illegal variable

Input script line contained a variable that could not be substituted
for.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Cannot use include command within an if command

Self-explanatory.

E: Cannot open logfile %s

The SPARTA log file specified in the input script cannot be opened.
Check that the path and name are correct.

E: Partition numeric index is out of bounds

It must be an integer from 1 to the number of partitions.

E: Cannot open print file %s

Self-explanatory.

E: Invalid collide style

The choice of collision style is unknown.

E: Cannot use collide_modify with no collisions defined

A collision style must be specified first.

E: Dimension command after simulation box is defined

The dimension command cannot be used after a read_data,
read_restart, or create_box command.

E: Invalid react style

The choice of reaction style is unknown.

E: Units command after simulation box is defined

The units command cannot be used after a read_data, read_restart, or
create_box command.

*/

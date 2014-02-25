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
  void mixture();
  void react_command();
  void region();
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
  void weight();
};

}

#endif

/* ERROR/WARNING messages:

*/

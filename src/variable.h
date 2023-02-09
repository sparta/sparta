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

#ifndef SPARTA_VARIABLE_H
#define SPARTA_VARIABLE_H

#include "stdio.h"
#include "pointers.h"

namespace SPARTA_NS {

class Variable : protected Pointers {
 public:
  Variable(class SPARTA *);
  ~Variable();
  void set(int, char **);
  void set(char *, int, char **);
  int next(int, char **);
  int find(char *);

  int equal_style(int);
  int particle_style(int);
  int grid_style(int);
  int surf_style(int);
  int internal_style(int);

  char *retrieve(char *);
  double compute_equal(int);
  double compute_equal(char *);
  void compute_particle(int, double *, int, int);
  void compute_grid(int, double *, int, int);
  void compute_surf(int, double *, int, int) {}  // not yet supported
  void internal_set(int, double);

  int int_between_brackets(char *&, int);
  double evaluate_boolean(char *);

 private:
  int me;
  int nvar;                // # of defined variables
  int maxvar;              // max # of variables following lists can hold
  char **names;            // name of each variable
  int *style;              // style of each variable
  int *num;                // # of values for each variable
  int *which;              // next available value for each variable
  int *pad;                // 1 = pad loop/uloop variables with 0s, 0 = no pad
  class VarReader **reader;   // variable that reads from file
  char ***data;            // str value of each variable's values
  double *dvalue;          // single numeric value for internal variables

  int *eval_in_progress;   // flag if evaluation of variable is in progress

  class RanKnuth *randomequal;     // RNG for equal-style vars
  class RanKnuth *randomparticle;  // RNG for particle-style vars

  int precedence[17];      // precedence level of math operators
                           // set length to include up to OR in enum

  int treestyle;           // tree being used for particle or grid-style var

                           // local copies of compute vector_grid vectors
  int nvec_storage;        // # of vectors currently stored locally
  int maxvec_storage;      // max # of vectors in vec_storage
  double **vec_storage;    // list of vector copies
  int *maxlen_storage;     // allocated length of each vector

  struct Tree {            // parse tree for particle-style variables
    double value;          // single scalar
    double *array;         // per-atom or per-type list of doubles
    char *carray;          // ptr into data struct with nstride = sizeof(struct)
    int type;              // operation, see enum{} in variable.cpp
    int nstride;           // stride between atoms if array is a 2d array
    int selfalloc;         // 1 if array is allocated here, else 0
    int ivalue1,ivalue2;   // extra values for needed for gmask,rmask,grmask
    Tree *left,*middle,*right;    // ptrs further down tree
  };

  void remove(int);
  void grow();
  void copy(int, char **, char **);
  double evaluate(char *, Tree **);
  double collapse_tree(Tree *);
  double eval_tree(Tree *, int);
  void free_tree(Tree *);
  int find_matching_paren(char *, int, char *&);
  int math_function(char *, char *, Tree **, Tree **, int &, double *, int &);
  int special_function(char *, char *, Tree **, Tree **,
                       int &, double *, int &);
  int is_particle_vector(char *);
  void particle_vector(char *, Tree **, Tree **, int &);
  int is_grid_vector(char *);
  void grid_vector(char *, Tree **, Tree **, int &);
  int is_constant(char *);
  double constant(char *);
  char *find_next_comma(char *);
  void print_tree(Tree *, int);
  double *add_storage(double *);
};

class VarReader : protected Pointers {
 public:
  VarReader(class SPARTA *, char *, char *, int);
  ~VarReader();
  int read_scalar(char *);

 private:
  int me,style;
  FILE *fp;
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: World variable count doesn't match # of partitions

A world-style variable must specify a number of values equal to the
number of processor partitions.

E: Universe/uloop variable count < # of partitions

A universe or uloop style variable must specify a number of values >= to the
number of processor partitions.

E: All universe/uloop variables must have same # of values

Self-explanatory.

E: Cannot redefine variable as a different style

An equal-style variable can be re-defined but only if it was
originally an equal-style variable.

E: File variable could not read value

Check the file assigned to the variable.

E: Grid-style variables are not yet implemented

Self-explanatory.

E: Surf-style variables are not yet implemented

Self-explanatory.

E: Variable name must be alphanumeric or underscore characters

Self-explanatory.

E: Invalid variable in next command

Self-explanatory.

E: All variables in next command must be same style

Self-explanatory.

E: Invalid variable style with next command

Variable styles {equal} and {world} cannot be used in a next
command.

E: Next command must list all universe and uloop variables

This is to insure they stay in sync.

E: Variable has circular dependency

A circular dependency is when variable "a" in used by variable "b" and
variable "b" is also used by varaible "a".  Circular dependencies with
longer chains of dependence are also not allowed.

E: Invalid syntax in variable formula

Self-explanatory.

E: Variable evaluation before simulation box is defined

Cannot evaluate a compute or fix or atom-based value in a variable
before the simulation has been setup.

E: Invalid compute ID in variable formula

The compute is not recognized.

E: Compute used in variable between runs is not current

Computes cannot be invoked by a variable in between runs.  Thus they
must have been evaluated on the last timestep of the previous run in
order for their value(s) to be accessed.  See the doc page for the
variable command for more info.

E: Variable formula compute vector is accessed out-of-range

Self-explanatory.

E: Variable formula compute array is accessed out-of-range

Self-explanatory.

E: Per-particle compute in equal-style variable formula

Equal-style variables cannot use per-particle quantities.

E: Mismatched compute in variable formula

A compute is referenced incorrectly or a compute that produces per-atom
values is used in an equal-style variable formula.

E: Invalid fix ID in variable formula

The fix is not recognized.

E: Fix in variable not computed at compatible time

Fixes generate their values on specific timesteps.  The variable is
requesting the values on a non-allowed timestep.

E: Variable formula fix vector is accessed out-of-range

Self-explanatory.

E: Variable formula fix array is accessed out-of-range

Self-explanatory.

E: Per-particle fix in equal-style variable formula

Equal-style variables cannot use per-particle quantities.

E: Mismatched fix in variable formula

A fix is referenced incorrectly or a fix that produces per-atom
values is used in an equal-style variable formula.

E: Invalid variable name in variable formula

Variable name is not recognized.

E: Invalid variable evaluation in variable formula

A variable used in a formula could not be evaluated.

E: Particle-style variable in equal-style variable formula

Equal-style variables cannot use per-particle quantities.

E: Mismatched variable in variable formula

A variable is referenced incorrectly or an atom-style variable that
produces per-atom values is used in an equal-style variable
formula.

E: Invalid math/special function in variable formula

Self-explanatory.

E: Too many particles per processor for particle-style variable

The number of particles per MPI rank is too large, increase the
number of MPI ranks.

E: Too many grid cells per processor for grid-style variable

The number of grid cells per MPI rank is too large, increase the
number of MPI ranks.

E: Invalid stats keyword in variable formula

The keyword is not recognized.

E: Divide by 0 in variable formula

Self-explanatory.

E: Modulo 0 in variable formula

Self-explanatory.

E: Power by 0 in variable formula

Self-explanatory.

E: Sqrt of negative value in variable formula

Self-explanatory.

E: Log of zero/negative value in variable formula

Self-explanatory.

E: Arcsin of invalid value in variable formula

Argument of arcsin() must be between -1 and 1.

E: Arccos of invalid value in variable formula

Argument of arccos() must be between -1 and 1.

E: Invalid math function in variable formula

Self-explanatory.

E: Variable name between brackets must be alphanumeric or underscore characters

Self-explanatory.

E: Non digit character between brackets in variable

Self-explanatory.

E: Mismatched brackets in variable

Self-explanatory.

E: Empty brackets in variable

There is no variable syntax that uses empty brackets.  Check
the variable doc page.

E: Index between variable brackets must be positive

Self-explanatory.

E: Cannot use ramp in variable formula between runs

This is because the ramp() function is time dependent.

E: Cannot use vdisplace in variable formula between runs

This is a function of elapsed time.

E: Cannot use swiggle in variable formula between runs

This is a function of elapsed time.

E: Cannot use cwiggle in variable formula between runs

This is a function of elapsed time.

E: Invalid special function in variable formula

Self-explanatory.

E: Variable ID in variable formula does not exist

Self-explanatory.

E: Invalid variable style in special function next

Only file-style or atomfile-style variables can be used with next().

E: Particle vector in equal-style variable formula

Equal-style variables cannot use per-particle quantities.

E: Invalid Boolean syntax in if command

Self-explanatory.

E: Cannot open file variable file %s

The specified file cannot be opened.  Check that the path and name are
correct.

*/

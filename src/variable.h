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

#ifndef DSMC_VARIABLE_H
#define DSMC_VARIABLE_H

#include "pointers.h"

namespace DSMC_NS {

class Variable : protected Pointers {
 public:
  Variable(class DSMC *);
  ~Variable();
  void set(int, char **);
  void set(char *, int, char **);
  int next(int, char **);
  int find(char *);
  int equal_style(int);
  int particle_style(int);
  int cell_style(int);
  char *retrieve(char *);
  double compute_equal(int);
  void compute_particle(int, double *, int, int);
  void compute_grid(int, double *, int, int);
  int int_between_brackets(char *&);
  double evaluate_boolean(char *);

 private:
  int me;
  int nvar;                // # of defined variables
  int maxvar;              // max # of variables arrays can hold
  char **names;            // name of each variable
  int *style;              // style of each variable
  int *num;                // # of values for each variable
  int *which;              // next available value for each variable
  int *pad;                // 1 = pad loop/uloop variables with 0s, 0 = no pad
  char ***data;            // str value of each variable's values
  double PI;

  class RanPark *randomequal;   // RNG for equal-style vars
  class RanPark *randompart;    // RNG for particle-style vars

  int precedence[16];      // precedence level of math operators
                           // set length to include OR in enum

  struct Tree {            // parse tree for particle-style variables
    double value;
    double *array;
    char *carray;
    int type;
    int nstride;
    int ivalue1,ivalue2;
    Tree *left,*middle,*right;
  };

  void remove(int);
  void extend();
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
  int is_constant(char *);
  double constant(char *);
  double numeric(char *);
  int inumeric(char *);
  char *find_next_comma(char *);
  void print_tree(Tree *, int);
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running DSMC to see the offending line.

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

E: Variable name must be alphanumeric or underscore characters

Self-explanatory.

E: Invalid variable in next command

Self-explanatory.

E: All variables in next command must be same style

Self-explanatory.

E: Invalid variable style with next command

Variable styles {equal} and {world} cannot be used in a next
command.

E: Invalid syntax in variable formula

Self-explanatory.

E: Invalid variable name in variable formula

Variable name is not recognized.

E: Invalid variable evaluation in variable formula

A variable used in a formula could not be evaluated.

E: Atom-style variable in equal-style variable formula

Atom-style variables generate one value per atom which is not allowed
in an equal-style variable.

E: Mismatched variable in variable formula

A variable is referenced incorrectly or an atom-style variable that
produces per-atom values is used in an equal-style variable
formula.

E: Invalid math/group/special function in variable formula

Self-explanatory.

E: Divide by 0 in variable formula

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

E: Non digit character between brackets in variable

Self-explantory.

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

E: Group ID in variable formula does not exist

Self-explanatory.

E: Invalid group function in variable formula

Group function is not recognized.

E: Region ID in variable formula does not exist

Self-explanatory.

E: Invalid special function in variable formula

Self-explanatory.

E: Invalid compute ID in variable formula

The compute is not recognized.

E: Compute used in variable between runs is not current

Computes cannot be invoked by a variable in between runs.  Thus they
must have been evaluated on the last timestep of the previous run in
order for their value(s) to be accessed.  See the doc page for the
variable command for more info.

E: Variable formula compute array is accessed out-of-range

Self-explanatory.

E: Mismatched compute in variable formula

A compute is referenced incorrectly or a compute that produces per-atom
values is used in an equal-style variable formula.

E: Invalid fix ID in variable formula

The fix is not recognized.

E: Fix in variable not computed at compatible time

Fixes generate their values on specific timesteps.  The variable is
requesting the values on a non-allowed timestep.

E: Variable formula fix array is accessed out-of-range

Self-explanatory.

E: Mismatched fix in variable formula

A fix is referenced incorrectly or a fix that produces per-atom
values is used in an equal-style variable formula.

E: Gmask function in equal-style variable formula

Gmask is per-atom operation.

E: Rmask function in equal-style variable formula

Rmask is per-atom operation.

E: Grmask function in equal-style variable formula

Grmask is per-atom operation.

E: Indexed per-atom vector in variable formula without atom map

Accessing a value from an atom vector requires the ability to lookup
an atom index, which is provided by an atom map.  An atom map does not
exist (by default) for non-molecular problems.  Using the atom_modify
map command will force an atom map to be created.

E: Invalid atom vector in variable formula

The atom vector is not recognized.

E: Atom vector in equal-style variable formula

Atom vectors generate one value per atom which is not allowed
in an equal-style variable.

E: Expected floating point parameter in variable definition

The quantity being read is a non-numeric value.

E: Expected integer parameter in variable definition

The quantity being read is a floating point or non-numeric value.

E: Invalid Boolean syntax in if command

Self-explanatory.

*/

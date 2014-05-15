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

#ifndef SPARTA_VARIABLE_H
#define SPARTA_VARIABLE_H

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
  char *retrieve(char *);
  double compute_equal(int);
  double compute_equal(char *);
  void compute_particle(int, double *, int, int);
  void compute_grid(int, double *, int, int) {}
  void compute_surf(int, double *, int, int) {}
  int int_between_brackets(char *&);
  double evaluate_boolean(char *);

 private:
  int nvar;                // # of defined variables
  int maxvar;              // max # of variables following lists can hold
  char **names;            // name of each variable
  int *style;              // style of each variable
  int *num;                // # of values for each variable
  int *which;              // next available value for each variable
  int *pad;                // 1 = pad loop/uloop variables with 0s, 0 = no pad
  class VarReader **reader;   // variable that reads from file
  char ***data;            // str value of each variable's values

  int *eval_in_progress;   // flag if evaluation of variable is in progress

  class RanPark *randomequal;     // RNG for equal-style vars
  class RanPark *randomparticle;  // RNG for particle-style vars

  int precedence[17];      // precedence level of math operators
                           // set length to include up to OR in enum
  int me;

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
  int is_constant(char *);
  double constant(char *);
  char *find_next_comma(char *);
  void print_tree(Tree *, int);
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

*/

/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "ctype.h"
#include "unistd.h"
#include "variable_kokkos.h"
#include "universe.h"
#include "update.h"
#include "domain.h"
#include "comm.h"
#include "particle_kokkos.h"
#include "grid_kokkos.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "surf_kokkos.h"
#include "surf_collide.h"
#include "surf_react.h"
#include "input.h"
#include "output.h"
#include "spapython.h"
#include "stats.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "utils.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "sparta_masks.h"

using namespace SPARTA_NS;
using namespace MathConst;

#define VARDELTA 4
#define MAXLEVEL 4
#define MAXLINE 256
#define CHUNK 1024
#define MAXFUNCARG 6

#define MYROUND(a) (( a-floor(a) ) >= .5) ? ceil(a) : floor(a)

enum{INDEX,LOOP,WORLD,UNIVERSE,ULOOP,STRING,GETENV,
  SCALARFILE,FORMAT,EQUAL,PARTICLE,GRID,SURF,INTERNAL,PYTHON};
enum{ARG,OP};

enum{INT,DOUBLE};                       // several files

enum{PARTICLE_CUSTOM,GRID_CUSTOM,SURF_CUSTOM};

// customize by adding a function
// if add before OR,
//   also set precedence level in constructor and precedence length in *.h

enum{DONE,ADD,SUBTRACT,MULTIPLY,DIVIDE,CARAT,MODULO,UNARY,
     NOT,EQ,NE,LT,LE,GT,GE,AND,OR,
     SQRT,EXP,LN,LOG,ABS,SIN,COS,TAN,ASIN,ACOS,ATAN,ATAN2,ERF,
     RANDOM,NORMAL,CEIL,FLOOR,ROUND,RAMP,STAGGER,LOGFREQ,STRIDE,
     VDISPLACE,SWIGGLE,CWIGGLE,
     PYWRAPPER,
     VALUE,ARRAY,ARRAYINT,PARTARRAYDOUBLE,PARTARRAYINT,SPECARRAY,PARTGRIDARRAY};

// customize by adding a special function

enum{SUM,XMIN,XMAX,AVE,TRAP,SLOPE,GRID2PART};

#define INVOKED_SCALAR 1
#define INVOKED_VECTOR 2
#define INVOKED_ARRAY 4
#define INVOKED_PER_PARTICLE 8
#define INVOKED_PER_GRID 16
#define INVOKED_PER_SURF 32

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

VariableKokkos::VariableKokkos(SPARTA *sparta) : Variable(sparta)
{

}

/* ----------------------------------------------------------------------
   recursive evaluation of a string str
   str is an equal-style or particle-style or grid-style formula
     containing one or more items:
     number = 0.0, -5.45, 2.8e-4, ...
     constant = PI
     stats keyword = step, np, vol, ...
     math operation = (),-x,x+y,x-y,x*y,x/y,x^y,
                      x==y,x!=y,x<y,x<=y,x>y,x>=y,x&&y,x||y,
                      sqrt(x),exp(x),ln(x),log(x),abs(x),
                      sin(x),cos(x),tan(x),asin(x),atan2(y,x),...
     special function = sum(x),min(x), ...
     python function wrapper = py_varname(x,y,z,...)
     particle vector = x, y, vx, ...
     grid vector = cxlo, cxhi, cylo, cyhi, czlo, czhi
     compute = c_ID, c_ID[i], c_ID[i][j]
     fix = f_ID, f_ID[i], f_ID[i][j]
     variable = v_name
   equal-style variable passes in tree = NULL:
     evaluate the formula, return result as a double
   particle-style or grid-style or surf-style variable passes in tree = non-NULL:
     parse the formula but do not evaluate it
     create a parse tree and return it
------------------------------------------------------------------------- */

double VariableKokkos::evaluate(char *str, Tree **tree)
{
  ParticleKokkos* particle_kk = ((ParticleKokkos*)particle);
  GridKokkos* grid_kk = ((GridKokkos*)grid);
  SurfKokkos* surf_kk = ((SurfKokkos*)surf);

  int op,opprevious;
  double value1,value2;
  char onechar;
  char *ptr;

  double argstack[MAXLEVEL];
  Tree *treestack[MAXLEVEL];
  int opstack[MAXLEVEL];
  int nargstack = 0;
  int ntreestack = 0;
  int nopstack = 0;

  int i = 0;
  int expect = ARG;

  while (1) {
    onechar = str[i];

    // whitespace: just skip

    if (isspace(onechar)) i++;

    // ----------------
    // parentheses: recursively evaluate contents of parens
    // ----------------

    else if (onechar == '(') {
      if (expect == OP) error->all(FLERR,"Invalid syntax in variable formula");
      expect = OP;

      char *contents;
      i = find_matching_paren(str,i,contents);
      i++;

      // evaluate contents and push on stack

      if (tree) {
        Tree *newtree;
        evaluate(contents,&newtree);
        treestack[ntreestack++] = newtree;
      } else argstack[nargstack++] = evaluate(contents,NULL);

      delete [] contents;

    // ----------------
    // number: push value onto stack
    // ----------------

    } else if (isdigit(onechar) || onechar == '.') {
      if (expect == OP) error->all(FLERR,"Invalid syntax in variable formula");
      expect = OP;

      // istop = end of number, including scientific notation

      int istart = i;
      while (isdigit(str[i]) || str[i] == '.') i++;
      if (str[i] == 'e' || str[i] == 'E') {
        i++;
        if (str[i] == '+' || str[i] == '-') i++;
        while (isdigit(str[i])) i++;
      }
      int istop = i - 1;

      int n = istop - istart + 1;
      char *number = new char[n+1];
      strncpy(number,&str[istart],n);
      number[n] = '\0';

      if (tree) {
        Tree *newtree = new Tree();
        newtree->type = VALUE;
        newtree->value = atof(number);
        treestack[ntreestack++] = newtree;
      } else argstack[nargstack++] = atof(number);

      delete [] number;

    // ----------------
    // letter: c_ID, c_ID[], c_ID[][], f_ID, f_ID[], f_ID[][],
    //         p_ID, p_ID[], g_ID, g_ID[], s_ID, s_ID[],
    //         sc_ID[], sr_ID[], v_name, exp(), x, PI, vol
    // ----------------

    } else if (isalpha(onechar)) {
      if (expect == OP) error->all(FLERR,"Invalid syntax in variable formula");
      expect = OP;

      // istop = end of word
      // word = all alphanumeric or underscore

      int istart = i;
      while (isalnum(str[i]) || str[i] == '_') i++;
      int istop = i-1;

      int n = istop - istart + 1;
      char *word = new char[n+1];
      strncpy(word,&str[istart],n);
      word[n] = '\0';

      // ----------------
      // compute
      // ----------------

      if (strncmp(word,"c_",2) == 0) {
        if (domain->box_exist == 0)
          error->all(FLERR,
                     "Variable evaluation before simulation box is defined");

        n = strlen(word) - 2 + 1;
        char *id = new char[n];
        strcpy(id,&word[2]);

        int icompute = modify->find_compute(id);
        if (icompute < 0)
          error->all(FLERR,"Invalid compute ID in variable formula");
        Compute *compute = modify->compute[icompute];
        delete [] id;

        // parse zero or one or two trailing brackets
        // point i beyond last bracket
        // nbracket = # of bracket pairs
        // index1,index2 = int inside each bracket pair

        int nbracket,index1,index2;
        if (str[i] != '[') nbracket = 0;
        else {
          nbracket = 1;
          ptr = &str[i];
          index1 = int_between_brackets(ptr,1);
          i = ptr-str+1;
          if (str[i] == '[') {
            nbracket = 2;
            ptr = &str[i];
            index2 = int_between_brackets(ptr,1);
            i = ptr-str+1;
          }
        }

        // c_ID = scalar from global scalar

        if (nbracket == 0 && compute->scalar_flag) {

          if (!compute->first_init)
            error->all(FLERR,"Variable formula compute cannot be invoked"
                       " before initialized");
          if (!(compute->invoked_flag & INVOKED_SCALAR)) {
            compute->compute_scalar();
            compute->invoked_flag |= INVOKED_SCALAR;
          }

          value1 = compute->scalar;
          if (tree) {
            Tree *newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;

        // c_ID[i] = scalar from global vector

        } else if (nbracket == 1 && compute->vector_flag) {

          if (index1 > compute->size_vector)
            error->all(FLERR,"Variable formula compute vector "
                       "is accessed out-of-range");
          if (!compute->first_init)
            error->all(FLERR,"Variable formula compute cannot be invoked"
                       " before initialized");
          if (!(compute->invoked_flag & INVOKED_VECTOR)) {
            compute->compute_vector();
            compute->invoked_flag |= INVOKED_VECTOR;
          }

          value1 = compute->vector[index1-1];
          if (tree) {
            Tree *newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;

        // c_ID[i][j] = scalar from global array

        } else if (nbracket == 2 && compute->array_flag) {

          if (index1 > compute->size_array_rows)
            error->all(FLERR,"Variable formula compute array "
                       "is accessed out-of-range");
          if (index2 > compute->size_array_cols)
            error->all(FLERR,"Variable formula compute array "
                       "is accessed out-of-range");
          if (!compute->first_init)
            error->all(FLERR,"Variable formula compute cannot be invoked"
                       " before initialized");
          if (!(compute->invoked_flag & INVOKED_ARRAY)) {
            compute->compute_array();
            compute->invoked_flag |= INVOKED_ARRAY;
          }

          value1 = compute->array[index1-1][index2-1];
          if (tree) {
            Tree *newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;

        // c_ID = vector from per-particle vector

        } else if (nbracket == 0 && compute->per_particle_flag &&
                   compute->size_per_particle_cols == 0) {

          if (tree == NULL || treestyle != PARTICLE)
            error->all(FLERR,"Per-particle compute in "
                       "non particle-style variable formula");
          if (!compute->first_init)
            error->all(FLERR,"Variable formula compute cannot be invoked"
                       " before initialized");
          if (!(compute->invoked_flag & INVOKED_PER_PARTICLE)) {
            compute->compute_per_particle();
            compute->invoked_flag |= INVOKED_PER_PARTICLE;
          }

          Tree *newtree = new Tree();
          newtree->type = ARRAY;
          newtree->array = compute->vector_particle;
          newtree->nstride = 1;
          treestack[ntreestack++] = newtree;

        // c_ID[i] = vector from per-particle array

        } else if (nbracket == 1 && compute->per_particle_flag &&
                   compute->size_per_particle_cols > 0) {

          if (tree == NULL || treestyle != PARTICLE)
            error->all(FLERR,"Per-particle compute in "
                       "non particle-style variable formula");
          if (index1 > compute->size_per_particle_cols)
            error->all(FLERR,"Variable formula compute array "
                       "is accessed out-of-range");
          if (!compute->first_init)
            error->all(FLERR,"Variable formula compute cannot be invoked"
                       " before initialized");
          if (!(compute->invoked_flag & INVOKED_PER_PARTICLE)) {
            compute->compute_per_particle();
            compute->invoked_flag |= INVOKED_PER_PARTICLE;
          }

          Tree *newtree = new Tree();
          newtree->type = ARRAY;
          newtree->array = &compute->array_particle[0][index1-1];
          newtree->nstride = compute->size_per_particle_cols;
          treestack[ntreestack++] = newtree;

        // c_ID = vector from per-grid vector

        } else if (nbracket == 0 && compute->per_grid_flag &&
                   compute->size_per_grid_cols == 0) {

          if (tree == NULL || treestyle != GRID)
            error->all(FLERR,"Per-grid compute in "
                       "non grid-style variable formula");
          if (!compute->first_init)
            error->all(FLERR,"Variable formula compute cannot be invoked"
                       " before initialized");
          if (!(compute->invoked_flag & INVOKED_PER_GRID)) {
            compute->compute_per_grid();
            compute->invoked_flag |= INVOKED_PER_GRID;
          }

          if (compute->post_process_grid_flag)
            compute->post_process_grid(0,1,NULL,NULL,NULL,1);
          else if (compute->post_process_isurf_grid_flag)
            compute->post_process_isurf_grid();

          Tree *newtree = new Tree();
          newtree->type = ARRAY;
          newtree->array = compute->vector_grid;
          newtree->nstride = 1;
          treestack[ntreestack++] = newtree;

        // c_ID[i] = vector from per-grid array
        // if compute sets post_process_grid_flag:
        //   then values are in computes's vector_grid,
        //   must store them locally via add_storage()
        //   since compute's vector_grid will be overwritten
        //   if this variable accesses multiple columns from compute's array

        } else if (nbracket == 1 && compute->per_grid_flag &&
                   compute->size_per_grid_cols > 0) {

          if (tree == NULL || treestyle != GRID)
            error->all(FLERR,"Per-grid compute in "
                       "non grid-style variable formula");
          if (index1 > compute->size_per_grid_cols)
            error->all(FLERR,"Variable formula compute array "
                       "is accessed out-of-range");
          if (!compute->first_init)
            error->all(FLERR,"Variable formula compute cannot be invoked"
                       " before initialized");
          if (!(compute->invoked_flag & INVOKED_PER_GRID)) {
            compute->compute_per_grid();
            compute->invoked_flag |= INVOKED_PER_GRID;
          }

          if (compute->post_process_grid_flag)
            compute->post_process_grid(index1,1,NULL,NULL,NULL,1);
          else if (compute->post_process_isurf_grid_flag)
            compute->post_process_isurf_grid();

          Tree *newtree = new Tree();
          newtree->type = ARRAY;
          if (compute->post_process_grid_flag) {
            newtree->array = add_storage(compute->vector_grid);
            newtree->nstride = 1;
          } else {
            newtree->array = &compute->array_grid[0][index1-1];
            newtree->nstride = compute->size_per_grid_cols;
          }
          treestack[ntreestack++] = newtree;

	// c_ID = vector from per-surf vector

        } else if (nbracket == 0 && compute->per_surf_flag &&
                   compute->size_per_surf_cols == 0) {

          if (tree == NULL || treestyle != SURF)
            error->all(FLERR,"Per-surf compute in "
                       "non surf-style variable formula");
          if (!compute->first_init)
            error->all(FLERR,"Variable formula compute cannot be invoked"
                       " before initialized");
          if (!(compute->invoked_flag & INVOKED_PER_SURF)) {
            compute->compute_per_surf();
            compute->invoked_flag |= INVOKED_PER_SURF;
          }

	  compute->post_process_surf();

          Tree *newtree = new Tree();
          newtree->type = ARRAY;
          newtree->array = compute->vector_surf;
          newtree->nstride = 1;
          treestack[ntreestack++] = newtree;
	
	// c_ID[i] = vector from per-surf array

        } else if (nbracket == 1 && compute->per_surf_flag &&
                   compute->size_per_surf_cols > 0) {

          if (tree == NULL || treestyle != SURF)
            error->all(FLERR,"Per-surf compute in "
                       "non surf-style variable formula");
          if (index1 > compute->size_per_surf_cols)
            error->all(FLERR,"Variable formula compute array "
                       "is accessed out-of-range");
          if (!compute->first_init)
            error->all(FLERR,"Variable formula compute cannot be invoked"
                       " before initialized");
          if (!(compute->invoked_flag & INVOKED_PER_SURF)) {
            compute->compute_per_surf();
            compute->invoked_flag |= INVOKED_PER_SURF;
          }

	  compute->post_process_surf();

          Tree *newtree = new Tree();
          newtree->type = ARRAY;
	  newtree->array = &compute->array_surf[0][index1-1];
	  newtree->nstride = compute->size_per_surf_cols;
          treestack[ntreestack++] = newtree;

	// unrecognized compute
	
        } else error->all(FLERR,"Mismatched compute in variable formula");

      // ----------------
      // fix
      // ----------------

      } else if (strncmp(word,"f_",2) == 0) {
        if (domain->box_exist == 0)
          error->all(FLERR,
                     "Variable evaluation before simulation box is defined");

        n = strlen(word) - 2 + 1;
        char *id = new char[n];
        strcpy(id,&word[2]);

        int ifix = modify->find_fix(id);
        if (ifix < 0) error->all(FLERR,"Invalid fix ID in variable formula");
        Fix *fix = modify->fix[ifix];
        delete [] id;

        // parse zero or one or two trailing brackets
        // point i beyond last bracket
        // nbracket = # of bracket pairs
        // index1,index2 = int inside each bracket pair

        int nbracket,index1,index2;
        if (str[i] != '[') nbracket = 0;
        else {
          nbracket = 1;
          ptr = &str[i];
          index1 = int_between_brackets(ptr,1);
          i = ptr-str+1;
          if (str[i] == '[') {
            nbracket = 2;
            ptr = &str[i];
            index2 = int_between_brackets(ptr,1);
            i = ptr-str+1;
          }
        }

        if (nbracket == 0 && fix->scalar_flag) {

          if (update->runflag > 0 && update->ntimestep % fix->global_freq)
            error->all(FLERR,"Fix in variable not computed at compatible time");

          value1 = fix->compute_scalar();
          if (tree) {
            Tree *newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;

        // f_ID[i] = scalar from global vector

        } else if (nbracket == 1 && fix->vector_flag) {

          if (index1 > fix->size_vector)
            error->all(FLERR,
                       "Variable formula fix vector is accessed out-of-range");
          if (update->runflag > 0 && update->ntimestep % fix->global_freq)
            error->all(FLERR,"Fix in variable not computed at compatible time");

          value1 = fix->compute_vector(index1-1);
          if (tree) {
            Tree *newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;

        // f_ID[i][j] = scalar from global array

        } else if (nbracket == 2 && fix->array_flag) {

          if (index1 > fix->size_array_rows)
            error->all(FLERR,
                       "Variable formula fix array is accessed out-of-range");
          if (index2 > fix->size_array_cols)
            error->all(FLERR,
                       "Variable formula fix array is accessed out-of-range");
          if (update->runflag > 0 && update->ntimestep % fix->global_freq)
            error->all(FLERR,"Fix in variable not computed at compatible time");

          value1 = fix->compute_array(index1-1,index2-1);
          if (tree) {
            Tree *newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;

        // f_ID = vector from per-particle vector

        } else if (nbracket == 0 && fix->per_particle_flag &&
                   fix->size_per_particle_cols == 0) {

          if (tree == NULL)
            error->all(FLERR,
                       "Per-particle fix in equal-style variable formula");
          if (update->runflag > 0 &&
              update->ntimestep % fix->per_particle_freq)
            error->all(FLERR,"Fix in variable not computed at compatible time");

          Tree *newtree = new Tree();
          newtree->type = ARRAY;
          newtree->array = fix->vector_particle;
          newtree->nstride = 1;
          treestack[ntreestack++] = newtree;

        // f_ID[i] = vector from per-particle array

        } else if (nbracket == 1 && fix->per_particle_flag &&
                   fix->size_per_particle_cols > 0) {

          if (tree == NULL)
            error->all(FLERR,
                       "Per-particle fix in equal-style variable formula");
          if (index1 > fix->size_per_particle_cols)
            error->all(FLERR,
                       "Variable formula fix array is accessed out-of-range");
          if (update->runflag > 0 &&
              update->ntimestep % fix->per_particle_freq)
            error->all(FLERR,"Fix in variable not computed at compatible time");

          Tree *newtree = new Tree();
          newtree->type = ARRAY;
          newtree->array = &fix->array_particle[0][index1-1];
          newtree->nstride = fix->size_per_particle_cols;
          treestack[ntreestack++] = newtree;

        // f_ID = vector from per-grid vector

        } else if (nbracket == 0 && fix->per_grid_flag &&
                   fix->size_per_grid_cols == 0) {

          if (tree == NULL || treestyle != GRID)
            error->all(FLERR,"Per-grid fix in "
                       "non grid-style variable formula");
          if (update->runflag > 0 &&
              update->ntimestep % fix->per_grid_freq)
            error->all(FLERR,"Fix in variable not computed at compatible time");

          Tree *newtree = new Tree();
          newtree->type = ARRAY;
          newtree->array = fix->vector_grid;
          newtree->nstride = 1;
          treestack[ntreestack++] = newtree;

        // f_ID[i] = vector from per-grid array

        } else if (nbracket == 1 && fix->per_grid_flag &&
                   fix->size_per_grid_cols > 0) {

          if (tree == NULL || treestyle != GRID)
            error->all(FLERR,"Per-grid fix in "
                       "non grid-style variable formula");
          if (index1 > fix->size_per_grid_cols)
            error->all(FLERR,
                       "Variable formula fix array is accessed out-of-range");
          if (update->runflag > 0 &&
              update->ntimestep % fix->per_grid_freq)
            error->all(FLERR,"Fix in variable not computed at compatible time");

          Tree *newtree = new Tree();
          newtree->type = ARRAY;
          newtree->array = &fix->array_grid[0][index1-1];
          newtree->nstride = fix->size_per_grid_cols;
          treestack[ntreestack++] = newtree;

        // f_ID = vector from per-surf vector

        } else if (nbracket == 0 && fix->per_surf_flag &&
                   fix->size_per_surf_cols == 0) {

          if (tree == NULL || treestyle != SURF)
            error->all(FLERR,"Per-surf fix in "
                       "non surf-style variable formula");
          if (update->runflag > 0 &&
              update->ntimestep % fix->per_surf_freq)
            error->all(FLERR,"Fix in variable not computed at compatible time");

          Tree *newtree = new Tree();
          newtree->type = ARRAY;
          newtree->array = fix->vector_surf;
          newtree->nstride = 1;
          treestack[ntreestack++] = newtree;

        // f_ID[i] = vector from per-surf array

        } else if (nbracket == 1 && fix->per_surf_flag &&
                   fix->size_per_surf_cols > 0) {

          if (tree == NULL || treestyle != SURF)
            error->all(FLERR,"Per-surf fix in "
                       "non surf-style variable formula");
          if (index1 > fix->size_per_surf_cols)
            error->all(FLERR,
                       "Variable formula fix array is accessed out-of-range");
          if (update->runflag > 0 &&
              update->ntimestep % fix->per_surf_freq)
            error->all(FLERR,"Fix in variable not computed at compatible time");

          Tree *newtree = new Tree();
          newtree->type = ARRAY;
          newtree->array = &fix->array_surf[0][index1-1];
          newtree->nstride = fix->size_per_surf_cols;
          treestack[ntreestack++] = newtree;

	// unrecognized fix
	
	} else error->all(FLERR,"Mismatched fix in variable formula");

      // ----------------
      // custom per-particle, per-grid, per-surf data
      // ----------------

      } else if ((strncmp(word,"p_",2) == 0) ||
                 (strncmp(word,"g_",2) == 0) ||
                 (strncmp(word,"s_",2) == 0)) {

        if (domain->box_exist == 0)
          error->all(FLERR,
                     "Custom attribute evaluation before simulation box is defined");

	int cwhich;
	if (strncmp(word,"p_",2) == 0) cwhich = PARTICLE_CUSTOM;
	else if (strncmp(word,"g_",2) == 0) cwhich = GRID_CUSTOM;
	else if (strncmp(word,"s_",2) == 0) cwhich = SURF_CUSTOM;

        if (cwhich == PARTICLE_CUSTOM)
          particle_kk->sync(Host,CUSTOM_MASK);	
        else if (cwhich == GRID_CUSTOM)
          grid_kk->sync(Host,CUSTOM_MASK);	
        else if (cwhich == SURF_CUSTOM)
          surf_kk->sync(Host,CUSTOM_MASK);	
	
        n = strlen(word) - 2 + 1;
        char *id = new char[n];
        strcpy(id,&word[2]);

	int icustom,size,type;
	if (cwhich == PARTICLE_CUSTOM) {
	  if (tree == NULL || treestyle != PARTICLE)
	    error->all(FLERR,"Per-particle custom attribute in "
		       "non particle-style variable formula");
	  icustom = particle->find_custom(id);
	  if (icustom < 0)
	    error->all(FLERR,"Invalid custom attribute ID in variable formula");
	  size = particle->esize[icustom];
	  type = particle->etype[icustom];
	} else if (cwhich == GRID_CUSTOM) {
	  if (tree == NULL || treestyle != GRID)
	    error->all(FLERR,"Per-grid custom attribute in "
		       "non grid-style variable formula");
	  icustom = grid->find_custom(id);
	  if (icustom < 0)
	    error->all(FLERR,"Invalid custom attribute ID in variable formula");
	  size = grid->esize[icustom];
	  type = grid->etype[icustom];
	} else if (cwhich == SURF_CUSTOM) {
	  if (tree == NULL || treestyle != SURF)
	    error->all(FLERR,"Per-surf custom attribute in "
		       "non surf-style variable formula");
	  icustom = surf->find_custom(id);
	  if (icustom < 0)
	    error->all(FLERR,"Invalid custom attribute ID in variable formula");
	  size = surf->esize[icustom];
	  type = surf->etype[icustom];
	}
	
        delete [] id;

        // parse zero or one or two trailing brackets
        // point i beyond last bracket
        // nbracket = # of bracket pairs
        // index1,index2 = int inside each bracket pair

        int nbracket,index1;
        if (str[i] != '[') nbracket = 0;
        else {
          nbracket = 1;
          ptr = &str[i];
          index1 = int_between_brackets(ptr,1);
          i = ptr-str+1;
          if (str[i] == '[') {
            nbracket = 2;
            ptr = &str[i];
            i = ptr-str+1;
          }
        }
	
	if (nbracket == 0 && size == 0) {

	  Tree *newtree = new Tree();
	  if (type == INT) {
	    newtree->type = ARRAYINT;
	    if (cwhich == PARTICLE_CUSTOM)
	      newtree->iarray = particle->eivec[particle->ewhich[icustom]];
	    else if (cwhich == GRID_CUSTOM)
	      newtree->iarray = grid->eivec[grid->ewhich[icustom]];
	    else if (cwhich == SURF_CUSTOM)
	      newtree->iarray = surf->eivec[surf->ewhich[icustom]];
	  } else if (type == DOUBLE) {
	    newtree->type = ARRAY;
	    if (cwhich == PARTICLE_CUSTOM)
	      newtree->array = particle->edvec[particle->ewhich[icustom]];
	    else if (cwhich == GRID_CUSTOM)
	      newtree->array = grid->edvec[grid->ewhich[icustom]];
	    else if (cwhich == SURF_CUSTOM)
	      newtree->array = surf->edvec[surf->ewhich[icustom]];
	  }
	  newtree->nstride = 1;
	  treestack[ntreestack++] = newtree;
	
	} else if (nbracket == 1 && size > 0) {
	
	  Tree *newtree = new Tree();
	  if (type == INT) {
	    newtree->type = ARRAYINT;
	    if (cwhich == PARTICLE_CUSTOM)
	      newtree->iarray = particle->eiarray[particle->ewhich[icustom]][index1-1];
	    else if (cwhich == GRID_CUSTOM)
	      newtree->iarray = grid->eiarray[grid->ewhich[icustom]][index1-1];
	    else if (cwhich == SURF_CUSTOM)
	      newtree->iarray = surf->eiarray[surf->ewhich[icustom]][index1-1];
	  } else if (type == DOUBLE) {
	    newtree->type = ARRAY;
	    if (cwhich == PARTICLE_CUSTOM)
	      newtree->array = particle->edvec[particle->ewhich[icustom]];
	    else if (cwhich == GRID_CUSTOM)
	      newtree->array = grid->edvec[grid->ewhich[icustom]];
	    else if (cwhich == SURF_CUSTOM)
	      newtree->array = surf->edvec[surf->ewhich[icustom]];
	  }
	  newtree->nstride = size;
	  treestack[ntreestack++] = newtree;

	// unrecognized custom attribute
	
	} else error->all(FLERR,"Mismatched custom attribute in variable formula");

      // ----------------
      // surface collide model
      // ----------------

      } else if (strncmp(word,"sc_",3) == 0) {
        if (domain->box_exist == 0)
          error->all(FLERR,
                     "Variable evaluation before simulation box is defined");

        n = strlen(word) - 3 + 1;
        char *id = new char[n];
        strcpy(id,&word[3]);

        int isc = surf->find_collide(id);
        if (isc < 0)
          error->all(FLERR,"Invalid surf collide ID in variable formula");
        SurfCollide *sc = surf->sc[isc];
        delete [] id;

        // parse zero or one or two trailing brackets
        // point i beyond last bracket
        // nbracket = # of bracket pairs
        // index1 = int inside each bracket pair

        int nbracket,index1;
        if (str[i] != '[') nbracket = 0;
        else {
          nbracket = 1;
          ptr = &str[i];
          index1 = int_between_brackets(ptr,1);
          i = ptr-str+1;
          if (str[i] == '[') {
            nbracket = 2;
            ptr = &str[i];
            i = ptr-str+1;
          }
        }

        // sc_ID[i] = scalar from global vector

        if (nbracket == 1 && sc->vector_flag) {
          if (index1 > sc->size_vector)
            error->all(FLERR,"Variable formula surf collide vector "
                       "is accessed out-of-range");

          value1 = sc->compute_vector(index1-1);
          if (tree) {
            Tree *newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;

        } else error->all(FLERR,"Mismatched surf collide in variable formula");

      // ----------------
      // surface reaction model
      // ----------------

      } else if (strncmp(word,"sr_",3) == 0) {
        if (domain->box_exist == 0)
          error->all(FLERR,
                     "Variable evaluation before simulation box is defined");

        n = strlen(word) - 3 + 1;
        char *id = new char[n];
        strcpy(id,&word[3]);

        int isr = surf->find_react(id);
        if (isr < 0)
          error->all(FLERR,"Invalid surf reaction ID in variable formula");
        SurfReact *sr = surf->sr[isr];
        delete [] id;

        // parse zero or one or two trailing brackets
        // point i beyond last bracket
        // nbracket = # of bracket pairs
        // index1 = int inside each bracket pair

        int nbracket,index1;
        if (str[i] != '[') nbracket = 0;
        else {
          nbracket = 1;
          ptr = &str[i];
          index1 = int_between_brackets(ptr,1);
          i = ptr-str+1;
          if (str[i] == '[') {
            nbracket = 2;
            ptr = &str[i];
            int_between_brackets(ptr,1);
            i = ptr-str+1;
          }
        }

        // sr_ID[i] = scalar from global vector

        if (nbracket == 1 && sr->vector_flag) {
          if (index1 > sr->size_vector)
            error->all(FLERR,"Variable formula surf reaction vector "
                       "is accessed out-of-range");

          value1 = sr->compute_vector(index1-1);
          if (tree) {
            Tree *newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;

        } else error->all(FLERR,"Mismatched surf reaction in variable formula");

      // ----------------
      // variable
      // ----------------

      } else if (strncmp(word,"v_",2) == 0) {
        n = strlen(word) - 2 + 1;
        char *id = new char[n];
        strcpy(id,&word[2]);

        int ivar = find(id);
        if (ivar < 0)
          error->all(FLERR,"Invalid variable name in variable formula");
        if (eval_in_progress[ivar])
          error->all(FLERR,"Variable has circular dependency");

        // parse zero or one trailing brackets
        // point i beyond last bracket
        // nbracket = # of bracket pairs

        int nbracket;
        if (str[i] != '[') nbracket = 0;
        else {
          nbracket = 1;
          ptr = &str[i];
          int_between_brackets(ptr,1);
          i = ptr-str+1;
        }

        // v_name = scalar from internal-style variable
        // access value directly

        if (nbracket == 0 && style[ivar] == INTERNAL) {

          value1 = dvalue[ivar];
          if (tree) {
            Tree *newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;

        // v_name = scalar from non particle/grid/surf variable
        // access value via retrieve()

        } else if (nbracket == 0 && style[ivar] != PARTICLE &&
                   style[ivar] != GRID && style[ivar] != SURF) {

          char *var = retrieve(id);
          if (var == NULL)
            error->all(FLERR,"Invalid variable evaluation in variable formula");
          if (tree) {
            Tree *newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = atof(var);
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = atof(var);

        // v_name = per-particle vector from particle-style variable
        // evaluate the particle-style variable as newtree

        } else if (nbracket == 0 && style[ivar] == PARTICLE) {

          if (tree == NULL || treestyle != PARTICLE)
            error->all(FLERR,"Per-particle variable in "
                       "non particle-style variable formula");
          Tree *newtree;
          evaluate(data[ivar][0],&newtree);
          treestack[ntreestack++] = newtree;

        // v_name = per-grid vector from grid-style variable
        // evaluate the grid-style variable as newtree

        } else if (nbracket == 0 && style[ivar] == GRID) {

          if (tree == NULL || treestyle != GRID)
            error->all(FLERR,"Per-grid variable in "
                       "non grid-style variable formula");
          Tree *newtree;
          evaluate(data[ivar][0],&newtree);
          treestack[ntreestack++] = newtree;

        // v_name = per-surf vector from surf-style variable
        // evaluate the surf-style variable as newtree

        } else if (nbracket == 0 && style[ivar] == SURF) {

          if (tree == NULL || treestyle != SURF)
            error->all(FLERR,"Per-surf variable in "
                       "non surf-style variable formula");
          Tree *newtree;
          evaluate(data[ivar][0],&newtree);
          treestack[ntreestack++] = newtree;

	// unrecognized variable
	
        } else error->all(FLERR,"Mismatched variable in variable formula");

        delete [] id;

      // ----------------
      // math/special function or particle vector or grid vector or
      //   constant or stats keyword
      // ----------------

      } else {

        // ----------------
        // math or special function
        // math_function() includes Python function wrapper
        // ----------------

        if (str[i] == '(') {
          char *contents;
          i = find_matching_paren(str,i,contents);
          i++;

          if (math_function(word,contents,tree,
                            treestack,ntreestack,argstack,nargstack));
          else if (special_function(word,contents,tree,
                                    treestack,ntreestack,argstack,nargstack));
          else error->all(FLERR,"Invalid math/special function "
                          "in variable formula");
          delete [] contents;

        // ----------------
        // particle vector
        // ----------------

        } else if (is_particle_vector(word)) {
          if (domain->box_exist == 0)
            error->all(FLERR,
                       "Variable evaluation before simulation box is defined");
          particle_vector(word,tree,treestack,ntreestack);

        // ----------------
        // grid vector
        // ----------------

        } else if (is_grid_vector(word)) {
          if (domain->box_exist == 0)
            error->all(FLERR,
                       "Variable evaluation before simulation box is defined");
          grid_vector(word,tree,treestack,ntreestack);

        // ----------------
        // constant
        // ----------------

        } else if (is_constant(word)) {
          value1 = constant(word);
          if (tree) {
            Tree *newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;

        // ----------------
        // stats keyword
        // ----------------

        } else {
          if (domain->box_exist == 0)
            error->all(FLERR,
                       "Variable evaluation before simulation box is defined");

          int flag = output->stats->evaluate_keyword(word,&value1);
          if (flag)
            error->all(FLERR,"Invalid stats keyword in variable formula");
          if (tree) {
            Tree *newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;
        }
      }

      delete [] word;

    // ----------------
    // math operator, including end-of-string
    // ----------------

    } else if (strchr("+-*/^<>=!&|%\0",onechar)) {
      if (onechar == '+') op = ADD;
      else if (onechar == '-') op = SUBTRACT;
      else if (onechar == '*') op = MULTIPLY;
      else if (onechar == '/') op = DIVIDE;
      else if (onechar == '%') op = MODULO;
      else if (onechar == '^') op = CARAT;
      else if (onechar == '=') {
        if (str[i+1] != '=')
          error->all(FLERR,"Invalid syntax in variable formula");
        op = EQ;
        i++;
      } else if (onechar == '!') {
        if (str[i+1] == '=') {
          op = NE;
          i++;
        } else op = NOT;
      } else if (onechar == '<') {
        if (str[i+1] != '=') op = LT;
        else {
          op = LE;
          i++;
        }
      } else if (onechar == '>') {
        if (str[i+1] != '=') op = GT;
        else {
          op = GE;
          i++;
        }
      } else if (onechar == '&') {
        if (str[i+1] != '&')
          error->all(FLERR,"Invalid syntax in variable formula");
        op = AND;
        i++;
      } else if (onechar == '|') {
        if (str[i+1] != '|')
          error->all(FLERR,"Invalid syntax in variable formula");
        op = OR;
        i++;
      } else op = DONE;

      i++;

      if (op == SUBTRACT && expect == ARG) {
        opstack[nopstack++] = UNARY;
        continue;
      }
      if (op == NOT && expect == ARG) {
        opstack[nopstack++] = op;
        continue;
      }

      if (expect == ARG) error->all(FLERR,"Invalid syntax in variable formula");
      expect = ARG;

      // evaluate stack as deep as possible while respecting precedence
      // before pushing current op onto stack

      while (nopstack && precedence[opstack[nopstack-1]] >= precedence[op]) {
        opprevious = opstack[--nopstack];

        if (tree) {
          Tree *newtree = new Tree();
          newtree->type = opprevious;
          if (opprevious == UNARY) {
            newtree->first = treestack[--ntreestack];
          } else {
            newtree->second = treestack[--ntreestack];
            newtree->first = treestack[--ntreestack];
          }
          treestack[ntreestack++] = newtree;

        } else {
          value2 = argstack[--nargstack];
          if (opprevious != UNARY && opprevious != NOT)
            value1 = argstack[--nargstack];

          if (opprevious == ADD)
            argstack[nargstack++] = value1 + value2;
          else if (opprevious == SUBTRACT)
            argstack[nargstack++] = value1 - value2;
          else if (opprevious == MULTIPLY)
            argstack[nargstack++] = value1 * value2;
          else if (opprevious == DIVIDE) {
            if (value2 == 0.0)
              error->one(FLERR,"Divide by 0 in variable formula");
            argstack[nargstack++] = value1 / value2;
          } else if (opprevious == MODULO) {
            if (value2 == 0.0)
              error->one(FLERR,"Modulo 0 in variable formula");
            argstack[nargstack++] = fmod(value1,value2);
          } else if (opprevious == CARAT) {
            if (value2 == 0.0)
              error->one(FLERR,"Power by 0 in variable formula");
            argstack[nargstack++] = pow(value1,value2);
          } else if (opprevious == UNARY) {
            argstack[nargstack++] = -value2;
          } else if (opprevious == NOT) {
            if (value2 == 0.0) argstack[nargstack++] = 1.0;
            else argstack[nargstack++] = 0.0;
          } else if (opprevious == EQ) {
            if (value1 == value2) argstack[nargstack++] = 1.0;
            else argstack[nargstack++] = 0.0;
          } else if (opprevious == NE) {
            if (value1 != value2) argstack[nargstack++] = 1.0;
            else argstack[nargstack++] = 0.0;
          } else if (opprevious == LT) {
            if (value1 < value2) argstack[nargstack++] = 1.0;
            else argstack[nargstack++] = 0.0;
          } else if (opprevious == LE) {
            if (value1 <= value2) argstack[nargstack++] = 1.0;
            else argstack[nargstack++] = 0.0;
          } else if (opprevious == GT) {
            if (value1 > value2) argstack[nargstack++] = 1.0;
            else argstack[nargstack++] = 0.0;
          } else if (opprevious == GE) {
            if (value1 >= value2) argstack[nargstack++] = 1.0;
            else argstack[nargstack++] = 0.0;
          } else if (opprevious == AND) {
            if (value1 != 0.0 && value2 != 0.0) argstack[nargstack++] = 1.0;
            else argstack[nargstack++] = 0.0;
          } else if (opprevious == OR) {
            if (value1 != 0.0 || value2 != 0.0) argstack[nargstack++] = 1.0;
            else argstack[nargstack++] = 0.0;
          }
        }
      }

      // if end-of-string, break out of entire formula evaluation loop

      if (op == DONE) break;

      // push current operation onto stack

      opstack[nopstack++] = op;

    } else error->all(FLERR,"Invalid syntax in variable formula");
  }

  if (nopstack) error->all(FLERR,"Invalid syntax in variable formula");

  // for particle-style, grid-style, or surf-style variable, return remaining tree
  // for equal-style variable, return remaining arg

  if (tree) {
    if (ntreestack != 1) error->all(FLERR,"Invalid syntax in variable formula");
    *tree = treestack[0];
    return 0.0;
  } else {
    if (nargstack != 1) error->all(FLERR,"Invalid syntax in variable formula");
    return argstack[0];
  }
}

/* ----------------------------------------------------------------------
   evaluate a particle-style variable parse tree for particle I
     or a grid-style variable parse tree for grid cell I
     or a surf-style variable parse tree for surf element I
   tree was created by one-time parsing of formula string via evaulate()
   customize by adding a function:
     sqrt(),exp(),ln(),log(),sin(),cos(),tan(),asin(),acos(),atan(),
     atan2(y,x),random(x,y),normal(x,y),ceil(),floor(),round(),
     ramp(x,y),stagger(x,y),logfreq(x,y,z),stride(x,y,z),
     vdisplace(x,y),swiggle(x,y,z),cwiggle(x,y,z)
---------------------------------------------------------------------- */

double VariableKokkos::eval_tree(Tree *tree, int i)
{
  ParticleKokkos* particle_kk = ((ParticleKokkos*)particle);

  if (tree->type == SPECARRAY)
    particle_kk->sync(Host,PARTICLE_MASK);

  if (tree->type == PARTGRIDARRAY)
    particle_kk->sync(Host,PARTICLE_MASK);

  return Variable::eval_tree(tree,i);
}

/* ----------------------------------------------------------------------
   process a particle vector in formula
   push result onto tree
   word = particle vector
   customize by adding a particle vector:
     id,x,y,z,vx,vy,vz,type,mass,q,mu
------------------------------------------------------------------------- */

void VariableKokkos::particle_vector(char *word, Tree **tree,
                               Tree **treestack, int &ntreestack)
{
  ParticleKokkos* particle_kk = ((ParticleKokkos*)particle);

  if (strcmp(word,"x") == 0)
    particle_kk->sync(Host,PARTICLE_MASK);
  else if (strcmp(word,"y") == 0)
    particle_kk->sync(Host,PARTICLE_MASK);
  else if (strcmp(word,"z") == 0)
    particle_kk->sync(Host,PARTICLE_MASK);
  else if (strcmp(word,"vx") == 0)
    particle_kk->sync(Host,PARTICLE_MASK);
  else if (strcmp(word,"vy") == 0)
    particle_kk->sync(Host,PARTICLE_MASK);
  else if (strcmp(word,"vz") == 0)
    particle_kk->sync(Host,PARTICLE_MASK);

  else if (strcmp(word,"id") == 0)
    particle_kk->sync(Host,PARTICLE_MASK);
  else if (strcmp(word,"type") == 0)
    particle_kk->sync(Host,PARTICLE_MASK);
  else if (strcmp(word,"mass") == 0)
    particle_kk->sync(Host,SPECIES_MASK);
  else if (strcmp(word,"q") == 0)
    particle_kk->sync(Host,SPECIES_MASK);
  else if (strcmp(word,"mu") == 0)
    particle_kk->sync(Host,SPECIES_MASK);

  Variable::particle_vector(word,tree,treestack,ntreestack);
}

/* ----------------------------------------------------------------------
   process a grid vector in formula
   push result onto tree
   word = grid vector
   customize by adding a grid vector:
     cxlo,cxhi,cylo,cyhi,czlo,czhi
------------------------------------------------------------------------- */

void VariableKokkos::grid_vector(char *word, Tree **tree,
                           Tree **treestack, int &ntreestack)
{
  GridKokkos* grid_kk = ((GridKokkos*)grid);
  grid_kk->sync(Host,CELL_MASK);

  Variable::grid_vector(word,tree,treestack,ntreestack);
}

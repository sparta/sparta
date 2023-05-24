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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "ctype.h"
#include "unistd.h"
#include "variable.h"
#include "universe.h"
#include "update.h"
#include "domain.h"
#include "comm.h"
#include "particle.h"
#include "grid.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "surf.h"
#include "surf_collide.h"
#include "surf_react.h"
#include "input.h"
#include "output.h"
#include "stats.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

#define VARDELTA 4
#define MAXLEVEL 4
#define MAXLINE 256
#define CHUNK 1024
#define VALUELENGTH 64

#define MYROUND(a) (( a-floor(a) ) >= .5) ? ceil(a) : floor(a)

enum{INDEX,LOOP,WORLD,UNIVERSE,ULOOP,STRING,GETENV,
     SCALARFILE,FORMAT,EQUAL,PARTICLE,GRID,SURF,INTERNAL};
enum{ARG,OP};

// customize by adding a function
// if add before OR,
// also set precedence level in constructor and precedence length in *.h

enum{DONE,ADD,SUBTRACT,MULTIPLY,DIVIDE,CARAT,MODULO,UNARY,
     NOT,EQ,NE,LT,LE,GT,GE,AND,OR,
     SQRT,EXP,LN,LOG,ABS,SIN,COS,TAN,ASIN,ACOS,ATAN,ATAN2,ERF,
     RANDOM,NORMAL,CEIL,FLOOR,ROUND,RAMP,STAGGER,LOGFREQ,STRIDE,
     VDISPLACE,SWIGGLE,CWIGGLE,
     VALUE,ARRAY,PARTARRAYDOUBLE,PARTARRAYINT,SPECARRAY};

// customize by adding a special function

enum{SUM,XMIN,XMAX,AVE,TRAP,SLOPE};

#define INVOKED_SCALAR 1
#define INVOKED_VECTOR 2
#define INVOKED_ARRAY 4
#define INVOKED_PER_PARTICLE 8
#define INVOKED_PER_GRID 16
#define INVOKED_PER_SURF 32

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

Variable::Variable(SPARTA *sparta) : Pointers(sparta)
{
  MPI_Comm_rank(world,&me);

  nvar = maxvar = 0;
  names = NULL;
  style = NULL;
  num = NULL;
  which = NULL;
  pad = NULL;
  reader = NULL;
  data = NULL;
  dvalue = NULL;

  eval_in_progress = NULL;

  randomequal = NULL;
  randomparticle = NULL;

  precedence[DONE] = 0;
  precedence[OR] = 1;
  precedence[AND] = 2;
  precedence[EQ] = precedence[NE] = 3;
  precedence[LT] = precedence[LE] = precedence[GT] = precedence[GE] = 4;
  precedence[ADD] = precedence[SUBTRACT] = 5;
  precedence[MULTIPLY] = precedence[DIVIDE] = precedence[MODULO] = 6;
  precedence[CARAT] = 7;
  precedence[UNARY] = precedence[NOT] = 8;

  // local storage of compute vector_grid values
  // stored when a grid-style variable compute uses post_process() method

  maxvec_storage = 0;
  vec_storage = NULL;
  maxlen_storage = NULL;
}

/* ---------------------------------------------------------------------- */

Variable::~Variable()
{
  for (int i = 0; i < nvar; i++) {
    delete [] names[i];
    delete reader[i];
    if (style[i] == LOOP || style[i] == ULOOP) delete [] data[i][0];
    else for (int j = 0; j < num[i]; j++) delete [] data[i][j];
    delete [] data[i];
  }
  memory->sfree(names);
  memory->destroy(style);
  memory->destroy(num);
  memory->destroy(which);
  memory->destroy(pad);
  memory->sfree(reader);
  memory->sfree(data);
  memory->sfree(dvalue);

  memory->destroy(eval_in_progress);

  delete randomequal;
  delete randomparticle;

  for (int i = 0; i < maxvec_storage; i++)
    memory->destroy(vec_storage[i]);
  memory->sfree(vec_storage);
  memory->sfree(maxlen_storage);
}

/* ----------------------------------------------------------------------
   called by variable command in input script
------------------------------------------------------------------------- */

void Variable::set(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal variable command");

  int replaceflag = 0;

  // DELETE
  // doesn't matter if variable no longer exists

  if (strcmp(arg[1],"delete") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal variable command");
    if (find(arg[0]) >= 0) remove(find(arg[0]));
    return;

  // INDEX
  // num = listed args, which = 1st value, data = copied args

  } else if (strcmp(arg[1],"index") == 0) {
    if (narg < 3) error->all(FLERR,"Illegal variable command");
    if (find(arg[0]) >= 0) return;
    if (nvar == maxvar) grow();
    style[nvar] = INDEX;
    num[nvar] = narg - 2;
    which[nvar] = 0;
    pad[nvar] = 0;
    data[nvar] = new char*[num[nvar]];
    copy(num[nvar],&arg[2],data[nvar]);

  // LOOP
  // 1 arg + pad: num = N, which = 1st value, data = single string
  // 2 args + pad: num = N2, which = N1, data = single string

  } else if (strcmp(arg[1],"loop") == 0) {
    if (find(arg[0]) >= 0) return;
    if (nvar == maxvar) grow();
    style[nvar] = LOOP;
    int nfirst,nlast;
    if (narg == 3 || (narg == 4 && strcmp(arg[3],"pad") == 0)) {
      nfirst = 1;
      nlast = atoi(arg[2]);
      if (nlast <= 0) error->all(FLERR,"Illegal variable command");
      if (narg == 4 && strcmp(arg[3],"pad") == 0) {
        char digits[12];
        sprintf(digits,"%d",nlast);
        pad[nvar] = strlen(digits);
      } else pad[nvar] = 0;
    } else if (narg == 4 || (narg == 5 && strcmp(arg[4],"pad") == 0)) {
      nfirst = atoi(arg[2]);
      nlast = atoi(arg[3]);
      if (nfirst > nlast || nlast < 0)
        error->all(FLERR,"Illegal variable command");
      if (narg == 5 && strcmp(arg[4],"pad") == 0) {
        char digits[12];
        sprintf(digits,"%d",nlast);
        pad[nvar] = strlen(digits);
      } else pad[nvar] = 0;
    } else error->all(FLERR,"Illegal variable command");
    num[nvar] = nlast;
    which[nvar] = nfirst-1;
    data[nvar] = new char*[1];
    data[nvar][0] = NULL;

  // WORLD
  // num = listed args, which = partition this proc is in, data = copied args
  // error check that num = # of worlds in universe

  } else if (strcmp(arg[1],"world") == 0) {
    if (narg < 3) error->all(FLERR,"Illegal variable command");
    if (find(arg[0]) >= 0) return;
    if (nvar == maxvar) grow();
    style[nvar] = WORLD;
    num[nvar] = narg - 2;
    if (num[nvar] != universe->nworlds)
      error->all(FLERR,"World variable count doesn't match # of partitions");
    which[nvar] = universe->iworld;
    pad[nvar] = 0;
    data[nvar] = new char*[num[nvar]];
    copy(num[nvar],&arg[2],data[nvar]);

  // UNIVERSE and ULOOP
  // for UNIVERSE: num = listed args, data = copied args
  // for ULOOP: num = N, data = single string
  // which = partition this proc is in
  // universe proc 0 creates lock file
  // error check that all other universe/uloop variables are same length

  } else if (strcmp(arg[1],"universe") == 0 || strcmp(arg[1],"uloop") == 0) {
    if (strcmp(arg[1],"universe") == 0) {
      if (narg < 3) error->all(FLERR,"Illegal variable command");
      if (find(arg[0]) >= 0) return;
      if (nvar == maxvar) grow();
      style[nvar] = UNIVERSE;
      num[nvar] = narg - 2;
      pad[nvar] = 0;
      data[nvar] = new char*[num[nvar]];
      copy(num[nvar],&arg[2],data[nvar]);
    } else if (strcmp(arg[1],"uloop") == 0) {
      if (narg < 3 || narg > 4 || (narg == 4 && strcmp(arg[3],"pad") != 0))
        error->all(FLERR,"Illegal variable command");
      if (find(arg[0]) >= 0) return;
      if (nvar == maxvar) grow();
      style[nvar] = ULOOP;
      num[nvar] = atoi(arg[2]);
      data[nvar] = new char*[1];
      data[nvar][0] = NULL;
      if (narg == 4) {
        char digits[12];
        sprintf(digits,"%d",num[nvar]);
        pad[nvar] = strlen(digits);
      } else pad[nvar] = 0;
    }

    if (num[nvar] < universe->nworlds)
      error->all(FLERR,"Universe/uloop variable count < # of partitions");
    which[nvar] = universe->iworld;

    if (universe->me == 0) {
      FILE *fp = fopen("tmp.sparta.variable","w");
      fprintf(fp,"%d\n",universe->nworlds);
      fclose(fp);
    }

    for (int jvar = 0; jvar < nvar; jvar++)
      if (num[jvar] && (style[jvar] == UNIVERSE || style[jvar] == ULOOP) &&
          num[nvar] != num[jvar])
        error->all(FLERR,
                   "All universe/uloop variables must have same # of values");

  // STRING
  // replace pre-existing var if also style STRING (allows it to be reset)
  // num = 1, which = 1st value
  // data = 1 value, string to eval

  } else if (strcmp(arg[1],"string") == 0) {
    if (narg != 3) error->all(FLERR,"Illegal variable command");
    int ivar = find(arg[0]);
    if (ivar >= 0) {
      if (style[find(arg[0])] != STRING)
        error->all(FLERR,"Cannot redefine variable as a different style");
      delete [] data[ivar][0];
      copy(1,&arg[2],data[ivar]);
      replaceflag = 1;
    } else {
      if (nvar == maxvar) grow();
      style[nvar] = STRING;
      num[nvar] = 1;
      which[nvar] = 0;
      pad[nvar] = 0;
      data[nvar] = new char*[num[nvar]];
      copy(1,&arg[2],data[nvar]);
    }

  // GETENV
  // remove pre-existing var if also style GETENV (allows it to be reset)
  // num = 1, which = 1st value
  // data = 1 value, string to eval

  } else if (strcmp(arg[1],"getenv") == 0) {
    if (narg != 3) error->all(FLERR,"Illegal variable command");
    if (find(arg[0]) >= 0) {
      if (style[find(arg[0])] != GETENV)
        error->all(FLERR,"Cannot redefine variable as a different style");
      remove(find(arg[0]));
    }
    if (nvar == maxvar) grow();
    style[nvar] = GETENV;
    num[nvar] = 1;
    which[nvar] = 0;
    pad[nvar] = 0;
    data[nvar] = new char*[num[nvar]];
    copy(1,&arg[2],data[nvar]);
    data[nvar][1] = new char[VALUELENGTH];
    strcpy(data[nvar][1],"(undefined)");

  // SCALARFILE for strings or numbers
  // which = 1st value
  // data = 1 value, string to eval

  } else if (strcmp(arg[1],"file") == 0) {
    if (narg != 3) error->all(FLERR,"Illegal variable command");
    if (find(arg[0]) >= 0) return;
    if (nvar == maxvar) grow();
    style[nvar] = SCALARFILE;
    num[nvar] = 1;
    which[nvar] = 0;
    pad[nvar] = 0;
    data[nvar] = new char*[num[nvar]];
    data[nvar][0] = new char[MAXLINE];
    reader[nvar] = new VarReader(sparta,arg[0],arg[2],SCALARFILE);
    int flag = reader[nvar]->read_scalar(data[nvar][0]);
    if (flag) error->all(FLERR,"File variable could not read value");

  // FORMAT
  // num = 3, which = 1st value
  // data = 3 values
  //   1st is name of variable to eval, 2nd is format string,
  //   3rd is filled on retrieval

  } else if (strcmp(arg[1],"format") == 0) {
    if (narg != 4) error->all(FLERR,"Illegal variable command");
    if (find(arg[0]) >= 0) return;
    if (nvar == maxvar) grow();
    style[nvar] = FORMAT;
    num[nvar] = 3;
    which[nvar] = 0;
    pad[nvar] = 0;
    data[nvar] = new char*[num[nvar]];
    copy(2,&arg[2],data[nvar]);
    data[nvar][2] = new char[VALUELENGTH];
    strcpy(data[nvar][2],"(undefined)");

  // EQUAL
  // replace pre-existing var if also style EQUAL (allows it to be reset)
  // num = 2, which = 1st value
  // data = 2 values, 1st is string to eval, 2nd is filled on retrieval

  } else if (strcmp(arg[1],"equal") == 0) {
    if (narg != 3) error->all(FLERR,"Illegal variable command");
    int ivar = find(arg[0]);
    if (ivar >= 0) {
      if (style[find(arg[0])] != EQUAL)
        error->all(FLERR,"Cannot redefine variable as a different style");
      delete [] data[ivar][0];
      copy(1,&arg[2],data[ivar]);
      replaceflag = 1;
    } else {
      if (nvar == maxvar) grow();
      style[nvar] = EQUAL;
      num[nvar] = 2;
      which[nvar] = 0;
      pad[nvar] = 0;
      data[nvar] = new char*[num[nvar]];
      copy(1,&arg[2],data[nvar]);
      data[nvar][1] = new char[VALUELENGTH];
      strcpy(data[nvar][1],"(undefined)");
    }

  // PARTICLE
  // replace pre-existing var if also style PARTICLE (allows it to be reset)
  // num = 1, which = 1st value
  // data = 1 value, string to eval

  } else if (strcmp(arg[1],"particle") == 0) {
    if (narg != 3) error->all(FLERR,"Illegal variable command");
    int ivar = find(arg[0]);
    if (ivar >= 0) {
      if (style[find(arg[0])] != PARTICLE)
        error->all(FLERR,"Cannot redefine variable as a different style");
      delete [] data[ivar][0];
      copy(1,&arg[2],data[ivar]);
      replaceflag = 1;
    } else {
      if (nvar == maxvar) grow();
      style[nvar] = PARTICLE;
      num[nvar] = 1;
      which[nvar] = 0;
      pad[nvar] = 0;
      data[nvar] = new char*[num[nvar]];
      copy(1,&arg[2],data[nvar]);
    }

  // GRID
  // replace pre-existing var if also style GRID (allows it to be reset)
  // num = 1, which = 1st value
  // data = 1 value, string to eval

  } else if (strcmp(arg[1],"grid") == 0) {
    if (narg != 3) error->all(FLERR,"Illegal variable command");
    int ivar = find(arg[0]);
    if (ivar >= 0) {
      if (style[find(arg[0])] != GRID)
        error->all(FLERR,"Cannot redefine variable as a different style");
      delete [] data[ivar][0];
      copy(1,&arg[2],data[ivar]);
      replaceflag = 1;
    } else {
      if (nvar == maxvar) grow();
      style[nvar] = GRID;
      num[nvar] = 1;
      which[nvar] = 0;
      pad[nvar] = 0;
      data[nvar] = new char*[num[nvar]];
      copy(1,&arg[2],data[nvar]);
    }

  // SURF (not implemented yet)
  // replace pre-existing var if also style SURF (allows it to be reset)
  // num = 1, which = 1st value
  // data = 1 value, string to eval

  } else if (strcmp(arg[1],"surf") == 0) {
    error->all(FLERR,"Surf-style variables are not yet implemented");

  // INTERNAL
  // replace pre-existing var if also style INTERNAL (allows it to be reset)
  // num = 1, for string representation of dvalue, set by retrieve()
  // dvalue = numeric initialization from 2nd arg, reset by internal_set()

  } else if (strcmp(arg[1],"internal") == 0) {
    if (narg != 3) error->all(FLERR,"Illegal variable command");
    int ivar = find(arg[0]);
    if (ivar >= 0) {
      if (style[ivar] != INTERNAL)
        error->all(FLERR,"Cannot redefine variable as a different style");
      dvalue[nvar] = input->numeric(FLERR,arg[2]);
      replaceflag = 1;
    } else {
      if (nvar == maxvar) grow();
      style[nvar] = INTERNAL;
      num[nvar] = 1;
      which[nvar] = 0;
      pad[nvar] = 0;
      data[nvar] = new char*[num[nvar]];
      data[nvar][0] = new char[VALUELENGTH];
      dvalue[nvar] = input->numeric(FLERR,arg[2]);
    }

  } else error->all(FLERR,"Illegal variable command");

  // set name of variable, if not replacing one flagged with replaceflag
  // name must be all alphanumeric chars or underscores

  if (replaceflag) return;

  int n = strlen(arg[0]) + 1;
  names[nvar] = new char[n];
  strcpy(names[nvar],arg[0]);

  for (int i = 0; i < n-1; i++)
    if (!isalnum(names[nvar][i]) && names[nvar][i] != '_')
      error->all(FLERR,"Variable name must be alphanumeric or "
                 "underscore characters");
  nvar++;
}

/* ----------------------------------------------------------------------
   INDEX variable created by command-line argument
   make it INDEX rather than STRING so cannot be re-defined in input script
------------------------------------------------------------------------- */

void Variable::set(char *name, int narg, char **arg)
{
  char **newarg = new char*[2+narg];
  newarg[0] = name;
  newarg[1] = (char *) "index";
  for (int i = 0; i < narg; i++) newarg[2+i] = arg[i];
  set(2+narg,newarg);
  delete [] newarg;
}

/* ----------------------------------------------------------------------
   increment variable(s)
   return 0 if OK if successfully incremented
   return 1 if any variable is exhausted, free the variable to allow re-use
------------------------------------------------------------------------- */

int Variable::next(int narg, char **arg)
{
  int ivar;

  if (narg == 0) error->all(FLERR,"Illegal next command");

  // check that variables exist and are all the same style
  // exception: UNIVERSE and ULOOP variables can be mixed in same next command

  for (int iarg = 0; iarg < narg; iarg++) {
    ivar = find(arg[iarg]);
    if (ivar == -1) error->all(FLERR,"Invalid variable in next command");
    if (style[ivar] == ULOOP && style[find(arg[0])] == UNIVERSE) continue;
    else if (style[ivar] == UNIVERSE && style[find(arg[0])] == ULOOP) continue;
    else if (style[ivar] != style[find(arg[0])])
      error->all(FLERR,"All variables in next command must be same style");
  }

  // invalid styles: STRING, EQUAL, WORLD, PARTICLE, GRID, GETENV,
  //                 FORMAT, INTERNAL

  int istyle = style[find(arg[0])];
  if (istyle == STRING || istyle == EQUAL || istyle == WORLD ||
      istyle == GETENV || istyle == PARTICLE || istyle == GRID ||
      istyle == FORMAT || istyle == INTERNAL)
    error->all(FLERR,"Invalid variable style with next command");

  // if istyle = UNIVERSE or ULOOP, insure all such variables are incremented

  if (istyle == UNIVERSE || istyle == ULOOP)
    for (int i = 0; i < nvar; i++) {
      if (style[i] != UNIVERSE && style[i] != ULOOP) continue;
      int iarg = 0;
      for (iarg = 0; iarg < narg; iarg++)
        if (strcmp(arg[iarg],names[i]) == 0) break;
      if (iarg == narg)
        error->universe_one(FLERR,"Next command must list all "
                            "universe and uloop variables");
    }

  // increment all variables in list
  // if any variable is exhausted, set flag = 1 and remove var to allow re-use

  int flag = 0;

  if (istyle == INDEX || istyle == LOOP) {
    for (int iarg = 0; iarg < narg; iarg++) {
      ivar = find(arg[iarg]);
      which[ivar]++;
      if (which[ivar] >= num[ivar]) {
        flag = 1;
        remove(ivar);
      }
    }

  } else if (istyle == SCALARFILE) {

    for (int iarg = 0; iarg < narg; iarg++) {
      ivar = find(arg[iarg]);
      int done = reader[ivar]->read_scalar(data[ivar][0]);
      if (done) {
        flag = 1;
        remove(ivar);
      }
    }

  } else if (istyle == UNIVERSE || istyle == ULOOP) {

    // wait until lock file can be created and owned by proc 0 of this world
    // rename() is not atomic in practice, but no known simple fix
    //   means multiple procs can read/write file at the same time (bad!)
    // random delays help
    // delay for random fraction of 1 second before first rename() call
    // delay for random fraction of 1 second before subsequent tries
    // when successful, read next available index and Bcast it within my world

    int nextindex;
    if (me == 0) {
      int seed = 12345 + universe->me + which[find(arg[0])];
      RanMars *random = new RanMars(sparta);
      random->init(seed);
      int delay = (int) (1000000*random->uniform());
      usleep(delay);
      while (1) {
        if (!rename("tmp.sparta.variable","tmp.sparta.variable.lock")) break;
        delay = (int) (1000000*random->uniform());
        usleep(delay);
      }
      delete random;

      FILE *fp = fopen("tmp.sparta.variable.lock","r");
      int tmp = fscanf(fp,"%d",&nextindex);
      //printf("READ %d %d\n",universe->me,nextindex);
      fclose(fp);
      fp = fopen("tmp.sparta.variable.lock","w");
      fprintf(fp,"%d\n",nextindex+1);
      //printf("WRITE %d %d\n",universe->me,nextindex+1);
      fclose(fp);
      rename("tmp.sparta.variable.lock","tmp.sparta.variable");
      if (universe->uscreen)
        fprintf(universe->uscreen,
                "Increment via next: value %d on partition %d\n",
                nextindex+1,universe->iworld);
      if (universe->ulogfile)
        fprintf(universe->ulogfile,
                "Increment via next: value %d on partition %d\n",
                nextindex+1,universe->iworld);
    }
    MPI_Bcast(&nextindex,1,MPI_INT,0,world);

    // set all variables in list to nextindex
    // must increment all UNIVERSE and ULOOP variables here
    // error check above tested for this

    for (int iarg = 0; iarg < narg; iarg++) {
      ivar = find(arg[iarg]);
      which[ivar] = nextindex;
      if (which[ivar] >= num[ivar]) {
        flag = 1;
        remove(ivar);
      }
    }
  }

  return flag;
}

/* ----------------------------------------------------------------------
   return ptr to the data text associated with a variable
   if INDEX or WORLD or UNIVERSE or STRING or SCALARFILE var,
     return ptr to stored string
   if LOOP or ULOOP var, write int to data[0] and return ptr to string
   if EQUAL var, evaluate variable and put result in str
   if FORMAT var, evaluate its variable and put formatted result in str
   if GETENV var, query environment and put result in str
   if PARTICLE or GRID var, return NULL
   if INTERNAL, convert dvalue and put result in str
   return NULL if no variable with name or which value is bad,
     caller must respond
------------------------------------------------------------------------- */

char *Variable::retrieve(char *name)
{
  int ivar = find(name);
  if (ivar == -1) return NULL;
  if (which[ivar] >= num[ivar]) return NULL;

  char *str;
  if (style[ivar] == INDEX || style[ivar] == WORLD ||
      style[ivar] == UNIVERSE || style[ivar] == STRING ||
      style[ivar] == SCALARFILE) {
    str = data[ivar][which[ivar]];
  } else if (style[ivar] == LOOP || style[ivar] == ULOOP) {
    char result[16];
    if (pad[ivar] == 0) sprintf(result,"%d",which[ivar]+1);
    else {
      char padstr[16];
      sprintf(padstr,"%%0%dd",pad[ivar]);
      sprintf(result,padstr,which[ivar]+1);
    }
    int n = strlen(result) + 1;
    delete [] data[ivar][0];
    data[ivar][0] = new char[n];
    strcpy(data[ivar][0],result);
    str = data[ivar][0];
  } else if (style[ivar] == EQUAL) {
    double answer = evaluate(data[ivar][0],NULL);
    sprintf(data[ivar][1],"%.15g",answer);
    str = data[ivar][1];
  } else if (style[ivar] == FORMAT) {
    int jvar = find(data[ivar][0]);
    if (jvar == -1) return NULL;
    if (!equal_style(jvar)) return NULL;
    double answer = compute_equal(jvar);
    sprintf(data[ivar][2],data[ivar][1],answer);
    str = data[ivar][2];
  } else if (style[ivar] == GETENV) {
    const char *result = getenv(data[ivar][0]);
    if (result == NULL) result = (const char *)"";
    int n = strlen(result) + 1;
    if (n > VALUELENGTH) {
      delete [] data[ivar][1];
      data[ivar][1] = new char[n];
    }
    strcpy(data[ivar][1],result);
    str = data[ivar][1];
  } else if (style[ivar] == INTERNAL) {
    sprintf(data[ivar][0],"%.15g",dvalue[ivar]);
    str = data[ivar][0];
  } else if (style[ivar] == PARTICLE || style[ivar] == GRID) return NULL;

  return str;
}

/* ----------------------------------------------------------------------
   return result of equal-style variable evaluation
   can be EQUAL or INTERNAL style
------------------------------------------------------------------------- */

double Variable::compute_equal(int ivar)
{
  if (eval_in_progress[ivar])
    error->all(FLERR,"Variable has circular dependency");

  eval_in_progress[ivar] = 1;

  double value;
  if (style[ivar] == EQUAL) value = evaluate(data[ivar][0],NULL);
  else if (style[ivar] == INTERNAL) value = dvalue[ivar];

  eval_in_progress[ivar] = 0;
  return value;
}

/* ----------------------------------------------------------------------
   return result of immediate equal-style variable evaluation
   called from Input::substitute()
------------------------------------------------------------------------- */

double Variable::compute_equal(char *str)
{
  return evaluate(str,NULL);
}

/* ----------------------------------------------------------------------
   compute result of particle-style variable evaluation
   answers are placed every stride locations into result
   if sumflag, add variable values to existing result
------------------------------------------------------------------------- */

void Variable::compute_particle(int ivar, double *result,
                                int stride, int sumflag)
{
  Tree *tree;

  if (eval_in_progress[ivar])
    error->all(FLERR,"Variable has circular dependency");
  eval_in_progress[ivar] = 1;

  treestyle = PARTICLE;
  evaluate(data[ivar][0],&tree);
  collapse_tree(tree);

  int nlocal = particle->nlocal;

  if (sumflag == 0) {
    int m = 0;
    for (int i = 0; i < nlocal; i++) {
      result[m] = eval_tree(tree,i);
      m += stride;
    }

  } else {
    int m = 0;
    for (int i = 0; i < nlocal; i++) {
      result[m] += eval_tree(tree,i);
      m += stride;
    }
  }

  free_tree(tree);

  eval_in_progress[ivar] = 0;
}

/* ----------------------------------------------------------------------
   compute result of grid-style variable evaluation
   answers are placed every stride locations into result
   if sumflag, add variable values to existing result
------------------------------------------------------------------------- */

void Variable::compute_grid(int ivar, double *result,
                            int stride, int sumflag)
{
  Tree *tree;

  if (eval_in_progress[ivar])
    error->all(FLERR,"Variable has circular dependency");
  eval_in_progress[ivar] = 1;

  nvec_storage = 0;
  treestyle = GRID;
  evaluate(data[ivar][0],&tree);
  collapse_tree(tree);

  int nglocal = grid->nlocal;

  if (sumflag == 0) {
    int m = 0;
    for (int i = 0; i < nglocal; i++) {
      result[m] = eval_tree(tree,i);
      m += stride;
    }

  } else {
    int m = 0;
    for (int i = 0; i < nglocal; i++) {
      result[m] += eval_tree(tree,i);
      m += stride;
    }
  }

  free_tree(tree);

  eval_in_progress[ivar] = 0;
}

/* ----------------------------------------------------------------------
   set value stored by INTERNAL style ivar
------------------------------------------------------------------------- */

void Variable::internal_set(int ivar, double value)
{
  dvalue[ivar] = value;
}

/* ----------------------------------------------------------------------
   search for name in list of variables names
   return index or -1 if not found
------------------------------------------------------------------------- */

int Variable::find(char *name)
{
  for (int i = 0; i < nvar; i++)
    if (strcmp(name,names[i]) == 0) return i;
  return -1;
}

/* ----------------------------------------------------------------------
   return 1 if variable is EQUAL or INTERNAL style, 0 if not
------------------------------------------------------------------------- */

int Variable::equal_style(int ivar)
{
  if (style[ivar] == EQUAL || style[ivar] == INTERNAL) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   return 1 if variable is PARTICLE style, 0 if not
------------------------------------------------------------------------- */

int Variable::particle_style(int ivar)
{
  if (style[ivar] == PARTICLE) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   return 1 if variable is GRID style, 0 if not
------------------------------------------------------------------------- */

int Variable::grid_style(int ivar)
{
  if (style[ivar] == GRID) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   return 1 if variable is SURF style, 0 if not
------------------------------------------------------------------------- */

int Variable::surf_style(int ivar)
{
  if (style[ivar] == SURF) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   return 1 if variable is INTERNAL style, 0 if not
   this is checked before call to set_internal() to assure it can be set
------------------------------------------------------------------------- */

int Variable::internal_style(int ivar)
{
  if (style[ivar] == INTERNAL) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   remove Nth variable from list and compact list
   delete reader explicitly if it exists
------------------------------------------------------------------------- */

void Variable::remove(int n)
{
  delete [] names[n];
  if (style[n] == LOOP || style[n] == ULOOP) delete [] data[n][0];
  else for (int i = 0; i < num[n]; i++) delete [] data[n][i];
  delete [] data[n];
  delete reader[n];

  for (int i = n+1; i < nvar; i++) {
    names[i-1] = names[i];
    style[i-1] = style[i];
    num[i-1] = num[i];
    which[i-1] = which[i];
    pad[i-1] = pad[i];
    reader[i-1] = reader[i];
    data[i-1] = data[i];
  }
  nvar--;
}

/* ----------------------------------------------------------------------
  make space in arrays for new variable
------------------------------------------------------------------------- */

void Variable::grow()
{
  int old = maxvar;
  maxvar += VARDELTA;
  names = (char **) memory->srealloc(names,maxvar*sizeof(char *),"var:names");
  memory->grow(style,maxvar,"var:style");
  memory->grow(num,maxvar,"var:num");
  memory->grow(which,maxvar,"var:which");
  memory->grow(pad,maxvar,"var:pad");

  reader = (VarReader **)
    memory->srealloc(reader,maxvar*sizeof(VarReader *),"var:reader");
  for (int i = old; i < maxvar; i++) reader[i] = NULL;

  data = (char ***) memory->srealloc(data,maxvar*sizeof(char **),"var:data");
  memory->grow(dvalue,maxvar,"var:dvalue");

  memory->grow(eval_in_progress,maxvar,"var:eval_in_progress");
  for (int i = 0; i < maxvar; i++) eval_in_progress[i] = 0;
}

/* ----------------------------------------------------------------------
   copy narg strings from **from to **to, and allocate space for them
------------------------------------------------------------------------- */

void Variable::copy(int narg, char **from, char **to)
{
  int n;
  for (int i = 0; i < narg; i++) {
    n = strlen(from[i]) + 1;
    to[i] = new char[n];
    strcpy(to[i],from[i]);
  }
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
     particle vector = x, y, vx, ...
     grid vector = cxlo, cxhi, cylo, cyhi, czlo, czhi
     compute = c_ID, c_ID[i], c_ID[i][j]
     fix = f_ID, f_ID[i], f_ID[i][j]
     variable = v_name
   equal-style variable passes in tree = NULL:
     evaluate the formula, return result as a double
   particle-style or grid-style variable passes in tree = non-NULL:
     parse the formula but do not evaluate it
     create a parse tree and return it
------------------------------------------------------------------------- */

double Variable::evaluate(char *str, Tree **tree)
{
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
        newtree->left = newtree->middle = newtree->right = NULL;
        treestack[ntreestack++] = newtree;
      } else argstack[nargstack++] = atof(number);

      delete [] number;

    // ----------------
    // letter: c_ID, c_ID[], c_ID[][], f_ID, f_ID[], f_ID[][], s_ID[], r_ID[],
    //         v_name, exp(), x, PI, vol
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

          if (update->runflag == 0) {
            if (compute->invoked_scalar != update->ntimestep)
              error->all(FLERR,"Compute used in variable between runs "
                         "is not current");
          } else if (!(compute->invoked_flag & INVOKED_SCALAR)) {
            compute->compute_scalar();
            compute->invoked_flag |= INVOKED_SCALAR;
          }

          value1 = compute->scalar;
          if (tree) {
            Tree *newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            newtree->left = newtree->middle = newtree->right = NULL;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;

        // c_ID[i] = scalar from global vector

        } else if (nbracket == 1 && compute->vector_flag) {

          if (index1 > compute->size_vector)
            error->all(FLERR,"Variable formula compute vector "
                       "is accessed out-of-range");
          if (update->runflag == 0) {
            if (compute->invoked_vector != update->ntimestep)
              error->all(FLERR,"Compute used in variable between runs "
                         "is not current");
          } else if (!(compute->invoked_flag & INVOKED_VECTOR)) {
            compute->compute_vector();
            compute->invoked_flag |= INVOKED_VECTOR;
          }

          value1 = compute->vector[index1-1];
          if (tree) {
            Tree *newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            newtree->left = newtree->middle = newtree->right = NULL;
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
          if (update->runflag == 0) {
            if (compute->invoked_array != update->ntimestep)
              error->all(FLERR,"Compute used in variable between runs "
                         "is not current");
          } else if (!(compute->invoked_flag & INVOKED_ARRAY)) {
            compute->compute_array();
            compute->invoked_flag |= INVOKED_ARRAY;
          }

          value1 = compute->array[index1-1][index2-1];
          if (tree) {
            Tree *newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            newtree->left = newtree->middle = newtree->right = NULL;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;

        // c_ID = vector from per-particle vector

        } else if (nbracket == 0 && compute->per_particle_flag &&
                   compute->size_per_particle_cols == 0) {

          if (tree == NULL || treestyle != PARTICLE)
            error->all(FLERR,"Per-particle compute in "
                       "non particle-style variable formula");
          if (update->runflag == 0) {
            if (compute->invoked_per_particle != update->ntimestep)
              error->all(FLERR,"Compute used in variable between runs "
                         "is not current");
          } else if (!(compute->invoked_flag & INVOKED_PER_PARTICLE)) {
            compute->compute_per_particle();
            compute->invoked_flag |= INVOKED_PER_PARTICLE;
          }

          Tree *newtree = new Tree();
          newtree->type = ARRAY;
          newtree->array = compute->vector_particle;
          newtree->nstride = 1;
          newtree->selfalloc = 0;
          newtree->left = newtree->middle = newtree->right = NULL;
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
          if (update->runflag == 0) {
            if (compute->invoked_per_particle != update->ntimestep)
              error->all(FLERR,"Compute used in variable between runs "
                         "is not current");
          } else if (!(compute->invoked_flag & INVOKED_PER_PARTICLE)) {
            compute->compute_per_particle();
            compute->invoked_flag |= INVOKED_PER_PARTICLE;
          }

          Tree *newtree = new Tree();
          newtree->type = ARRAY;
          newtree->array = &compute->array_particle[0][index1-1];
          newtree->nstride = compute->size_per_particle_cols;
          newtree->selfalloc = 0;
          newtree->left = newtree->middle = newtree->right = NULL;
          treestack[ntreestack++] = newtree;

        // c_ID = vector from per-grid vector

        } else if (nbracket == 0 && compute->per_grid_flag &&
                   compute->size_per_grid_cols == 0) {

          if (tree == NULL || treestyle != GRID)
            error->all(FLERR,"Per-grid compute in "
                       "non grid-style variable formula");
          if (update->runflag == 0) {
            if (compute->invoked_per_grid != update->ntimestep)
              error->all(FLERR,"Compute used in variable between runs "
                         "is not current");
          } else if (!(compute->invoked_flag & INVOKED_PER_GRID)) {
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
          newtree->selfalloc = 0;
          newtree->left = newtree->middle = newtree->right = NULL;
          treestack[ntreestack++] = newtree;

        // c_ID[i] = vector from per-grid array
        // if compute sets post_process_grid_flag:
        //   then values are in computes's vector_grid,
        //   must store them locally via add_vector()
        //   since compute's vector_grid be overwritten
        //   if this variable accesses multiple columns from compute's array

        } else if (nbracket == 1 && compute->per_grid_flag &&
                   compute->size_per_grid_cols > 0) {

          if (tree == NULL || treestyle != GRID)
            error->all(FLERR,"Per-grid compute in "
                       "non grid-style variable formula");
          if (index1 > compute->size_per_grid_cols)
            error->all(FLERR,"Variable formula compute array "
                       "is accessed out-of-range");
          if (update->runflag == 0) {
            if (compute->invoked_per_grid != update->ntimestep)
              error->all(FLERR,"Compute used in variable between runs "
                         "is not current");
          } else if (!(compute->invoked_flag & INVOKED_PER_GRID)) {
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
          newtree->selfalloc = 0;
          newtree->left = newtree->middle = newtree->right = NULL;
          treestack[ntreestack++] = newtree;

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
            newtree->left = newtree->middle = newtree->right = NULL;
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
            newtree->left = newtree->middle = newtree->right = NULL;
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
            newtree->left = newtree->middle = newtree->right = NULL;
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
          newtree->selfalloc = 0;
          newtree->left = newtree->middle = newtree->right = NULL;
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
          newtree->selfalloc = 0;
          newtree->left = newtree->middle = newtree->right = NULL;
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
          newtree->selfalloc = 0;
          newtree->left = newtree->middle = newtree->right = NULL;
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
          newtree->selfalloc = 0;
          newtree->left = newtree->middle = newtree->right = NULL;
          treestack[ntreestack++] = newtree;

        } else error->all(FLERR,"Mismatched fix in variable formula");

      // ----------------
      // surface collide model
      // ----------------

      } else if (strncmp(word,"s_",2) == 0) {
        if (domain->box_exist == 0)
          error->all(FLERR,
                     "Variable evaluation before simulation box is defined");

        n = strlen(word) - 2 + 1;
        char *id = new char[n];
        strcpy(id,&word[2]);

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

        // s_ID[i] = scalar from global vector

        if (nbracket == 1 && sc->vector_flag) {
          if (index1 > sc->size_vector)
            error->all(FLERR,"Variable formula surf collide vector "
                       "is accessed out-of-range");

          value1 = sc->compute_vector(index1-1);
          if (tree) {
            Tree *newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            newtree->left = newtree->middle = newtree->right = NULL;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;

        } else error->all(FLERR,"Mismatched surf collide in variable formula");

      // ----------------
      // surface reaction model
      // ----------------

      } else if (strncmp(word,"r_",2) == 0) {
        if (domain->box_exist == 0)
          error->all(FLERR,
                     "Variable evaluation before simulation box is defined");

        n = strlen(word) - 2 + 1;
        char *id = new char[n];
        strcpy(id,&word[2]);

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

        // r_ID[i] = scalar from global vector

        if (nbracket == 1 && sr->vector_flag) {
          if (index1 > sr->size_vector)
            error->all(FLERR,"Variable formula surf reaction vector "
                       "is accessed out-of-range");

          value1 = sr->compute_vector(index1-1);
          if (tree) {
            Tree *newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = value1;
            newtree->left = newtree->middle = newtree->right = NULL;
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
            newtree->left = newtree->middle = newtree->right = NULL;
            treestack[ntreestack++] = newtree;
          } else argstack[nargstack++] = value1;

        // v_name = scalar from non particle/grid variable
        // access value via retrieve()

        } else if (nbracket == 0 && style[ivar] != PARTICLE &&
                   style[ivar] != GRID) {

          char *var = retrieve(id);
          if (var == NULL)
            error->all(FLERR,"Invalid variable evaluation in variable formula");
          if (tree) {
            Tree *newtree = new Tree();
            newtree->type = VALUE;
            newtree->value = atof(var);
            newtree->left = newtree->middle = newtree->right = NULL;
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

        } else error->all(FLERR,"Mismatched variable in variable formula");

        delete [] id;

      // ----------------
      // math/special function or particle vector or constant or stats keyword
      // ----------------

      } else {

        // ----------------
        // math or special function
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
            newtree->left = newtree->middle = newtree->right = NULL;
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
            newtree->left = newtree->middle = newtree->right = NULL;
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
            newtree->left = treestack[--ntreestack];
            newtree->middle = newtree->right = NULL;
          } else {
            newtree->right = treestack[--ntreestack];
            newtree->middle = NULL;
            newtree->left = treestack[--ntreestack];
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

  // for particle-style variable, return remaining tree
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
   one-time collapse of a particle-style or grid-style variable parse tree
   tree was created by one-time parsing of formula string via evaulate()
   only keep tree nodes that depend on ARRAYs
   remainder is converted to single VALUE
   this enables optimal eval_tree loop over particles or grid cells
   customize by adding a function:
     sqrt(),exp(),ln(),log(),abs(),sin(),cos(),tan(),asin(),acos(),atan(),
     atan2(y,x),random(x,y),normal(x,y),ceil(),floor(),round(),
     ramp(x,y),stagger(x,y),logfreq(x,y,z),stride(x,y,z),
     vdisplace(x,y),swiggle(x,y,z),cwiggle(x,y,z)
---------------------------------------------------------------------- */

double Variable::collapse_tree(Tree *tree)
{
  double arg1,arg2;

  if (tree->type == VALUE) return tree->value;
  if (tree->type == ARRAY) return 0.0;
  if (tree->type == PARTARRAYDOUBLE) return 0.0;
  if (tree->type == PARTARRAYINT) return 0.0;
  if (tree->type == SPECARRAY) return 0.0;

  if (tree->type == ADD) {
    arg1 = collapse_tree(tree->left);
    arg2 = collapse_tree(tree->right);
    if (tree->left->type != VALUE || tree->right->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = arg1 + arg2;
    return tree->value;
  }

  if (tree->type == SUBTRACT) {
    arg1 = collapse_tree(tree->left);
    arg2 = collapse_tree(tree->right);
    if (tree->left->type != VALUE || tree->right->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = arg1 - arg2;
    return tree->value;
  }

  if (tree->type == MULTIPLY) {
    arg1 = collapse_tree(tree->left);
    arg2 = collapse_tree(tree->right);
    if (tree->left->type != VALUE || tree->right->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = arg1 * arg2;
    return tree->value;
  }

  if (tree->type == DIVIDE) {
    arg1 = collapse_tree(tree->left);
    arg2 = collapse_tree(tree->right);
    if (tree->left->type != VALUE || tree->right->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg2 == 0.0) error->one(FLERR,"Divide by 0 in variable formula");
    tree->value = arg1 / arg2;
    return tree->value;
  }

  if (tree->type == MODULO) {
    arg1 = collapse_tree(tree->left);
    arg2 = collapse_tree(tree->right);
    if (tree->left->type != VALUE || tree->right->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg2 == 0.0) error->one(FLERR,"Modulo 0 in variable formula");
    tree->value = fmod(arg1,arg2);
    return tree->value;
  }

  if (tree->type == CARAT) {
    arg1 = collapse_tree(tree->left);
    arg2 = collapse_tree(tree->right);
    if (tree->left->type != VALUE || tree->right->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg2 == 0.0) error->one(FLERR,"Power by 0 in variable formula");
    tree->value = pow(arg1,arg2);
    return tree->value;
  }

  if (tree->type == UNARY) {
    arg1 = collapse_tree(tree->left);
    if (tree->left->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = -arg1;
    return tree->value;
  }

  if (tree->type == NOT) {
    arg1 = collapse_tree(tree->left);
    if (tree->left->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 == 0.0) tree->value = 1.0;
    else tree->value = 0.0;
    return tree->value;
  }

  if (tree->type == EQ) {
    arg1 = collapse_tree(tree->left);
    arg2 = collapse_tree(tree->right);
    if (tree->left->type != VALUE || tree->right->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 == arg2) tree->value = 1.0;
    else tree->value = 0.0;
    return tree->value;
  }

  if (tree->type == NE) {
    arg1 = collapse_tree(tree->left);
    arg2 = collapse_tree(tree->right);
    if (tree->left->type != VALUE || tree->right->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 != arg2) tree->value = 1.0;
    else tree->value = 0.0;
    return tree->value;
  }

  if (tree->type == LT) {
    arg1 = collapse_tree(tree->left);
    arg2 = collapse_tree(tree->right);
    if (tree->left->type != VALUE || tree->right->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 < arg2) tree->value = 1.0;
    else tree->value = 0.0;
    return tree->value;
  }

  if (tree->type == LE) {
    arg1 = collapse_tree(tree->left);
    arg2 = collapse_tree(tree->right);
    if (tree->left->type != VALUE || tree->right->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 <= arg2) tree->value = 1.0;
    else tree->value = 0.0;
    return tree->value;
  }

  if (tree->type == GT) {
    arg1 = collapse_tree(tree->left);
    arg2 = collapse_tree(tree->right);
    if (tree->left->type != VALUE || tree->right->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 > arg2) tree->value = 1.0;
    else tree->value = 0.0;
    return tree->value;
  }

  if (tree->type == GE) {
    arg1 = collapse_tree(tree->left);
    arg2 = collapse_tree(tree->right);
    if (tree->left->type != VALUE || tree->right->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 >= arg2) tree->value = 1.0;
    else tree->value = 0.0;
    return tree->value;
  }

  if (tree->type == AND) {
    arg1 = collapse_tree(tree->left);
    arg2 = collapse_tree(tree->right);
    if (tree->left->type != VALUE || tree->right->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 != 0.0 && arg2 != 0.0) tree->value = 1.0;
    else tree->value = 0.0;
    return tree->value;
  }

  if (tree->type == OR) {
    arg1 = collapse_tree(tree->left);
    arg2 = collapse_tree(tree->right);
    if (tree->left->type != VALUE || tree->right->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 != 0.0 || arg2 != 0.0) tree->value = 1.0;
    else tree->value = 0.0;
    return tree->value;
  }

  if (tree->type == SQRT) {
    arg1 = collapse_tree(tree->left);
    if (tree->left->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 < 0.0)
      error->one(FLERR,"Sqrt of negative value in variable formula");
    tree->value = sqrt(arg1);
    return tree->value;
  }

  if (tree->type == EXP) {
    arg1 = collapse_tree(tree->left);
    if (tree->left->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = exp(arg1);
    return tree->value;
  }

  if (tree->type == LN) {
    arg1 = collapse_tree(tree->left);
    if (tree->left->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 <= 0.0)
      error->one(FLERR,"Log of zero/negative value in variable formula");
    tree->value = log(arg1);
    return tree->value;
  }

  if (tree->type == LOG) {
    arg1 = collapse_tree(tree->left);
    if (tree->left->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 <= 0.0)
      error->one(FLERR,"Log of zero/negative value in variable formula");
    tree->value = log10(arg1);
    return tree->value;
  }

  if (tree->type == ABS) {
    arg1 = collapse_tree(tree->left);
    if (tree->left->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = fabs(arg1);
    return tree->value;
  }

  if (tree->type == SIN) {
    arg1 = collapse_tree(tree->left);
    if (tree->left->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = sin(arg1);
    return tree->value;
  }

  if (tree->type == COS) {
    arg1 = collapse_tree(tree->left);
    if (tree->left->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = cos(arg1);
    return tree->value;
  }

  if (tree->type == TAN) {
    arg1 = collapse_tree(tree->left);
    if (tree->left->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = tan(arg1);
    return tree->value;
  }

  if (tree->type == ASIN) {
    arg1 = collapse_tree(tree->left);
    if (tree->left->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 < -1.0 || arg1 > 1.0)
      error->one(FLERR,"Arcsin of invalid value in variable formula");
    tree->value = asin(arg1);
    return tree->value;
  }

  if (tree->type == ACOS) {
    arg1 = collapse_tree(tree->left);
    if (tree->left->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg1 < -1.0 || arg1 > 1.0)
      error->one(FLERR,"Arccos of invalid value in variable formula");
    tree->value = acos(arg1);
    return tree->value;
  }

  if (tree->type == ATAN) {
    arg1 = collapse_tree(tree->left);
    if (tree->left->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = atan(arg1);
    return tree->value;
  }

  if (tree->type == ATAN2) {
    arg1 = collapse_tree(tree->left);
    arg2 = collapse_tree(tree->right);
    if (tree->left->type != VALUE || tree->right->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = atan2(arg1,arg2);
    return tree->value;
  }

  if (tree->type == ERF) {
    arg1 = collapse_tree(tree->left);
    if (tree->left->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = erf(arg1);
    return tree->value;
  }

  // random() or normal() do not become a single collapsed value

  if (tree->type == RANDOM) {
    collapse_tree(tree->left);
    collapse_tree(tree->right);
    if (randomparticle == NULL) {
      randomparticle = new RanKnuth(update->ranmaster->uniform());
      double seed = update->ranmaster->uniform();
      randomparticle->reset(seed,me,100);
    }
    return 0.0;
  }

  if (tree->type == NORMAL) {
    collapse_tree(tree->left);
    double sigma = collapse_tree(tree->right);
    if (sigma < 0.0)
      error->one(FLERR,"Invalid math function in variable formula");
    if (randomparticle == NULL) {
      randomparticle = new RanKnuth(update->ranmaster->uniform());
      double seed = update->ranmaster->uniform();
      randomparticle->reset(seed,me,100);
    }
    return 0.0;
  }

  if (tree->type == CEIL) {
    arg1 = collapse_tree(tree->left);
    if (tree->left->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = ceil(arg1);
    return tree->value;
  }

  if (tree->type == FLOOR) {
    arg1 = collapse_tree(tree->left);
    if (tree->left->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = floor(arg1);
    return tree->value;
  }

  if (tree->type == ROUND) {
    arg1 = collapse_tree(tree->left);
    if (tree->left->type != VALUE) return 0.0;
    tree->type = VALUE;
    tree->value = MYROUND(arg1);
    return tree->value;
  }

  if (tree->type == RAMP) {
    arg1 = collapse_tree(tree->left);
    arg2 = collapse_tree(tree->right);
    if (tree->left->type != VALUE || tree->right->type != VALUE) return 0.0;
    tree->type = VALUE;
    double delta = update->ntimestep - update->beginstep;
    if (delta != 0.0) delta /= update->endstep - update->beginstep;
    tree->value = arg1 + delta*(arg2-arg1);
    return tree->value;
  }

  if (tree->type == STAGGER) {
    int ivalue1 = static_cast<int> (collapse_tree(tree->left));
    int ivalue2 = static_cast<int> (collapse_tree(tree->right));
    if (tree->left->type != VALUE || tree->right->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (ivalue1 <= 0 || ivalue2 <= 0 || ivalue1 <= ivalue2)
      error->one(FLERR,"Invalid math function in variable formula");
    int lower = update->ntimestep/ivalue1 * ivalue1;
    int delta = update->ntimestep - lower;
    if (delta < ivalue2) tree->value = lower+ivalue2;
    else tree->value = lower+ivalue1;
    return tree->value;
  }

  if (tree->type == LOGFREQ) {
    int ivalue1 = static_cast<int> (collapse_tree(tree->left));
    int ivalue2 = static_cast<int> (collapse_tree(tree->middle));
    int ivalue3 = static_cast<int> (collapse_tree(tree->right));
    if (tree->left->type != VALUE || tree->middle->type != VALUE ||
        tree->right->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (ivalue1 <= 0 || ivalue2 <= 0 || ivalue3 <= 0 || ivalue2 >= ivalue3)
      error->one(FLERR,"Invalid math function in variable formula");
    if (update->ntimestep < ivalue1) tree->value = ivalue1;
    else {
      int lower = ivalue1;
      while (update->ntimestep >= ivalue3*lower) lower *= ivalue3;
      int multiple = update->ntimestep/lower;
      if (multiple < ivalue2) tree->value = (multiple+1)*lower;
      else tree->value = lower*ivalue3;
    }
    return tree->value;
  }

  if (tree->type == STRIDE) {
    int ivalue1 = static_cast<int> (collapse_tree(tree->left));
    int ivalue2 = static_cast<int> (collapse_tree(tree->middle));
    int ivalue3 = static_cast<int> (collapse_tree(tree->right));
    if (tree->left->type != VALUE || tree->middle->type != VALUE ||
        tree->right->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (ivalue1 < 0 || ivalue2 < 0 || ivalue3 <= 0 || ivalue1 > ivalue2)
      error->one(FLERR,"Invalid math function in variable formula");
    if (update->ntimestep < ivalue1) tree->value = ivalue1;
    else if (update->ntimestep < ivalue2) {
      int offset = update->ntimestep - ivalue1;
      tree->value = ivalue1 + (offset/ivalue3)*ivalue3 + ivalue3;
      if (tree->value > ivalue2) tree->value = 9.0e18;
    } else tree->value = 9.0e18;
    return tree->value;
  }

  if (tree->type == VDISPLACE) {
    double arg1 = collapse_tree(tree->left);
    double arg2 = collapse_tree(tree->right);
    if (tree->left->type != VALUE || tree->right->type != VALUE) return 0.0;
    tree->type = VALUE;
    double delta = update->ntimestep - update->beginstep;
    tree->value = arg1 + arg2*delta*update->dt;
    return tree->value;
  }

  if (tree->type == SWIGGLE) {
    double arg1 = collapse_tree(tree->left);
    double arg2 = collapse_tree(tree->middle);
    double arg3 = collapse_tree(tree->right);
    if (tree->left->type != VALUE || tree->middle->type != VALUE ||
        tree->right->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg3 == 0.0)
      error->one(FLERR,"Invalid math function in variable formula");
    double delta = update->ntimestep - update->beginstep;
    double omega = 2.0*MY_PI/arg3;
    tree->value = arg1 + arg2*sin(omega*delta*update->dt);
    return tree->value;
  }

  if (tree->type == CWIGGLE) {
    double arg1 = collapse_tree(tree->left);
    double arg2 = collapse_tree(tree->middle);
    double arg3 = collapse_tree(tree->right);
    if (tree->left->type != VALUE || tree->middle->type != VALUE ||
        tree->right->type != VALUE) return 0.0;
    tree->type = VALUE;
    if (arg3 == 0.0)
      error->one(FLERR,"Invalid math function in variable formula");
    double delta = update->ntimestep - update->beginstep;
    double omega = 2.0*MY_PI/arg3;
    tree->value = arg1 + arg2*(1.0-cos(omega*delta*update->dt));
    return tree->value;
  }

  return 0.0;
}

/* ----------------------------------------------------------------------
   evaluate a particle-style variable parse tree for particle I
     or a grid-style variable parse tree for grid cell I
   tree was created by one-time parsing of formula string via evaulate()
   customize by adding a function:
     sqrt(),exp(),ln(),log(),sin(),cos(),tan(),asin(),acos(),atan(),
     atan2(y,x),random(x,y),normal(x,y),ceil(),floor(),round(),
     ramp(x,y),stagger(x,y),logfreq(x,y,z),stride(x,y,z),
     vdisplace(x,y),swiggle(x,y,z),cwiggle(x,y,z)
---------------------------------------------------------------------- */

double Variable::eval_tree(Tree *tree, int i)
{
  double arg,arg1,arg2,arg3;

  if (tree->type == VALUE) return tree->value;
  if (tree->type == ARRAY) return tree->array[i*tree->nstride];
  if (tree->type == PARTARRAYDOUBLE)
    return *((double *) &tree->carray[i*tree->nstride]);
  if (tree->type == PARTARRAYINT)
    return *((int *) &tree->carray[i*tree->nstride]);
  if (tree->type == SPECARRAY)
    return *((double *)
             &tree->carray[particle->particles[i].ispecies*tree->nstride]);

  if (tree->type == ADD)
    return eval_tree(tree->left,i) + eval_tree(tree->right,i);
  if (tree->type == SUBTRACT)
    return eval_tree(tree->left,i) - eval_tree(tree->right,i);
  if (tree->type == MULTIPLY)
    return eval_tree(tree->left,i) * eval_tree(tree->right,i);
  if (tree->type == DIVIDE) {
    double denom = eval_tree(tree->right,i);
    if (denom == 0.0) error->one(FLERR,"Divide by 0 in variable formula");
    return eval_tree(tree->left,i) / denom;
  }
  if (tree->type == MODULO) {
    double denom = eval_tree(tree->right,i);
    if (denom == 0.0) error->one(FLERR,"Modulo 0 in variable formula");
    return fmod(eval_tree(tree->left,i),denom);
  }
  if (tree->type == CARAT) {
    double exponent = eval_tree(tree->right,i);
    if (exponent == 0.0) error->one(FLERR,"Power by 0 in variable formula");
    return pow(eval_tree(tree->left,i),exponent);
  }
  if (tree->type == UNARY) return -eval_tree(tree->left,i);

  if (tree->type == NOT) {
    if (eval_tree(tree->left,i) == 0.0) return 1.0;
    else return 0.0;
  }
  if (tree->type == EQ) {
    if (eval_tree(tree->left,i) == eval_tree(tree->right,i)) return 1.0;
    else return 0.0;
  }
  if (tree->type == NE) {
    if (eval_tree(tree->left,i) != eval_tree(tree->right,i)) return 1.0;
    else return 0.0;
  }
  if (tree->type == LT) {
    if (eval_tree(tree->left,i) < eval_tree(tree->right,i)) return 1.0;
    else return 0.0;
  }
  if (tree->type == LE) {
    if (eval_tree(tree->left,i) <= eval_tree(tree->right,i)) return 1.0;
    else return 0.0;
  }
  if (tree->type == GT) {
    if (eval_tree(tree->left,i) > eval_tree(tree->right,i)) return 1.0;
    else return 0.0;
  }
  if (tree->type == GE) {
    if (eval_tree(tree->left,i) >= eval_tree(tree->right,i)) return 1.0;
    else return 0.0;
  }
  if (tree->type == AND) {
    if (eval_tree(tree->left,i) != 0.0 && eval_tree(tree->right,i) != 0.0)
      return 1.0;
    else return 0.0;
  }
  if (tree->type == OR) {
    if (eval_tree(tree->left,i) != 0.0 || eval_tree(tree->right,i) != 0.0)
      return 1.0;
    else return 0.0;
  }

  if (tree->type == SQRT) {
    arg1 = eval_tree(tree->left,i);
    if (arg1 < 0.0)
      error->one(FLERR,"Sqrt of negative value in variable formula");
    return sqrt(arg1);
  }
  if (tree->type == EXP)
    return exp(eval_tree(tree->left,i));
  if (tree->type == LN) {
    arg1 = eval_tree(tree->left,i);
    if (arg1 <= 0.0)
      error->one(FLERR,"Log of zero/negative value in variable formula");
    return log(arg1);
  }
  if (tree->type == LOG) {
    arg1 = eval_tree(tree->left,i);
    if (arg1 <= 0.0)
      error->one(FLERR,"Log of zero/negative value in variable formula");
    return log10(arg1);
  }
  if (tree->type == ABS)
    return fabs(eval_tree(tree->left,i));

  if (tree->type == SIN)
    return sin(eval_tree(tree->left,i));
  if (tree->type == COS)
    return cos(eval_tree(tree->left,i));
  if (tree->type == TAN)
    return tan(eval_tree(tree->left,i));

  if (tree->type == ASIN) {
    arg1 = eval_tree(tree->left,i);
    if (arg1 < -1.0 || arg1 > 1.0)
      error->one(FLERR,"Arcsin of invalid value in variable formula");
    return asin(arg1);
  }
  if (tree->type == ACOS) {
    arg1 = eval_tree(tree->left,i);
    if (arg1 < -1.0 || arg1 > 1.0)
      error->one(FLERR,"Arccos of invalid value in variable formula");
    return acos(arg1);
  }
  if (tree->type == ATAN)
    return atan(eval_tree(tree->left,i));
  if (tree->type == ATAN2)
    return atan2(eval_tree(tree->left,i),eval_tree(tree->right,i));
  if (tree->type == ERF)
    return erf(eval_tree(tree->left,i));

  if (tree->type == RANDOM) {
    double lower = eval_tree(tree->left,i);
    double upper = eval_tree(tree->right,i);
    if (randomparticle == NULL) {
      randomparticle = new RanKnuth(update->ranmaster->uniform());
      double seed = update->ranmaster->uniform();
      randomparticle->reset(seed,me,100);
    }
    return randomparticle->uniform()*(upper-lower)+lower;
  }
  if (tree->type == NORMAL) {
    double mu = eval_tree(tree->left,i);
    double sigma = eval_tree(tree->right,i);
    if (sigma < 0.0)
      error->one(FLERR,"Invalid math function in variable formula");
    if (randomparticle == NULL) {
      randomparticle = new RanKnuth(update->ranmaster->uniform());
      double seed = update->ranmaster->uniform();
      randomparticle->reset(seed,me,100);
    }
    return mu + sigma*randomparticle->gaussian();
  }

  if (tree->type == CEIL)
    return ceil(eval_tree(tree->left,i));
  if (tree->type == FLOOR)
    return floor(eval_tree(tree->left,i));
  if (tree->type == ROUND)
    return MYROUND(eval_tree(tree->left,i));

  if (tree->type == RAMP) {
    arg1 = eval_tree(tree->left,i);
    arg2 = eval_tree(tree->right,i);
    double delta = update->ntimestep - update->beginstep;
    if (delta != 0.0) delta /= update->endstep - update->beginstep;
    arg = arg1 + delta*(arg2-arg1);
    return arg;
  }

  if (tree->type == STAGGER) {
    int ivalue1 = static_cast<int> (eval_tree(tree->left,i));
    int ivalue2 = static_cast<int> (eval_tree(tree->right,i));
    if (ivalue1 <= 0 || ivalue2 <= 0 || ivalue1 <= ivalue2)
      error->one(FLERR,"Invalid math function in variable formula");
    int lower = update->ntimestep/ivalue1 * ivalue1;
    int delta = update->ntimestep - lower;
    if (delta < ivalue2) arg = lower+ivalue2;
    else arg = lower+ivalue1;
    return arg;
  }

  if (tree->type == LOGFREQ) {
    int ivalue1 = static_cast<int> (eval_tree(tree->left,i));
    int ivalue2 = static_cast<int> (eval_tree(tree->middle,i));
    int ivalue3 = static_cast<int> (eval_tree(tree->right,i));
    if (ivalue1 <= 0 || ivalue2 <= 0 || ivalue3 <= 0 || ivalue2 >= ivalue3)
      error->one(FLERR,"Invalid math function in variable formula");
    if (update->ntimestep < ivalue1) arg = ivalue1;
    else {
      int lower = ivalue1;
      while (update->ntimestep >= ivalue3*lower) lower *= ivalue3;
      int multiple = update->ntimestep/lower;
      if (multiple < ivalue2) arg = (multiple+1)*lower;
      else arg = lower*ivalue3;
    }
    return arg;
  }

  if (tree->type == STRIDE) {
    int ivalue1 = static_cast<int> (eval_tree(tree->left,i));
    int ivalue2 = static_cast<int> (eval_tree(tree->middle,i));
    int ivalue3 = static_cast<int> (eval_tree(tree->right,i));
    if (ivalue1 < 0 || ivalue2 < 0 || ivalue3 <= 0 || ivalue1 > ivalue2)
      error->one(FLERR,"Invalid math function in variable formula");
    if (update->ntimestep < ivalue1) arg = ivalue1;
    else if (update->ntimestep < ivalue2) {
      int offset = update->ntimestep - ivalue1;
      arg = ivalue1 + (offset/ivalue3)*ivalue3 + ivalue3;
      if (arg > ivalue2) arg = 9.0e18;
    } else arg = 9.0e18;
    return arg;
  }

  if (tree->type == VDISPLACE) {
    arg1 = eval_tree(tree->left,i);
    arg2 = eval_tree(tree->right,i);
    double delta = update->ntimestep - update->beginstep;
    arg = arg1 + arg2*delta*update->dt;
    return arg;
  }

  if (tree->type == SWIGGLE) {
    arg1 = eval_tree(tree->left,i);
    arg2 = eval_tree(tree->middle,i);
    arg3 = eval_tree(tree->right,i);
    if (arg3 == 0.0)
      error->one(FLERR,"Invalid math function in variable formula");
    double delta = update->ntimestep - update->beginstep;
    double omega = 2.0*MY_PI/arg3;
    arg = arg1 + arg2*sin(omega*delta*update->dt);
    return arg;
  }

  if (tree->type == CWIGGLE) {
    arg1 = eval_tree(tree->left,i);
    arg2 = eval_tree(tree->middle,i);
    arg3 = eval_tree(tree->right,i);
    if (arg3 == 0.0)
      error->one(FLERR,"Invalid math function in variable formula");
    double delta = update->ntimestep - update->beginstep;
    double omega = 2.0*MY_PI/arg3;
    arg = arg1 + arg2*(1.0-cos(omega*delta*update->dt));
    return arg;
  }

  return 0.0;
}

/* ---------------------------------------------------------------------- */

void Variable::free_tree(Tree *tree)
{
  if (tree->left) free_tree(tree->left);
  if (tree->middle) free_tree(tree->middle);
  if (tree->right) free_tree(tree->right);
  delete tree;
}

/* ----------------------------------------------------------------------
   find matching parenthesis in str, allocate contents = str between parens
   i = left paren
   return loc or right paren
------------------------------------------------------------------------- */

int Variable::find_matching_paren(char *str, int i,char *&contents)
{
  // istop = matching ')' at same level, allowing for nested parens

  int istart = i;
  int ilevel = 0;
  while (1) {
    i++;
    if (!str[i]) break;
    if (str[i] == '(') ilevel++;
    else if (str[i] == ')' && ilevel) ilevel--;
    else if (str[i] == ')') break;
  }
  if (!str[i]) error->all(FLERR,"Invalid syntax in variable formula");
  int istop = i;

  int n = istop - istart - 1;
  contents = new char[n+1];
  strncpy(contents,&str[istart+1],n);
  contents[n] = '\0';

  return istop;
}

/* ----------------------------------------------------------------------
   find int between brackets and return it
   ptr initially points to left bracket
   return it pointing to right bracket
   error if no right bracket or brackets are empty or index = 0
   if varallow = 0: error if any between-bracket chars are non-digits
   if varallow = 1: also allow for v_name, where name is variable name
------------------------------------------------------------------------- */

int Variable::int_between_brackets(char *&ptr, int varallow)
{
  int varflag,index;

  char *start = ++ptr;

  if (varallow && strstr(ptr,"v_") == ptr) {
    varflag = 1;
    while (*ptr && *ptr != ']') {
      if (!isalnum(*ptr) && *ptr != '_')
        error->all(FLERR,"Variable name between brackets must be "
                   "alphanumeric or underscore characters");
      ptr++;
    }

  } else {
    varflag = 0;
    while (*ptr && *ptr != ']') {
      if (!isdigit(*ptr))
        error->all(FLERR,"Non digit character between brackets in variable");
      ptr++;
    }
  }

  if (*ptr != ']') error->all(FLERR,"Mismatched brackets in variable");
  if (ptr == start) error->all(FLERR,"Empty brackets in variable");

  *ptr = '\0';

  // evaluate index as variable or as simple integer via atoi()

  if (varflag) {
    char *id = start+2;
    int ivar = find(id);
    if (ivar < 0)
      error->all(FLERR,"Invalid variable name in variable formula");
    if (eval_in_progress[ivar])
      error->all(FLERR,"Variable has circular dependency");

    char *var = retrieve(id);
    if (var == NULL)
      error->all(FLERR,"Invalid variable evaluation in variable formula");
    index = static_cast<int> (atof(var));

  } else {
    index = atoi(start);
  }

  *ptr = ']';

  if (index == 0)
    error->all(FLERR,"Index between variable brackets must be positive");
  return index;
}

/* ----------------------------------------------------------------------
   process a math function in formula
   push result onto tree or arg stack
   word = math function
   contents = str between parentheses with one,two,three args
   return 0 if not a match, 1 if successfully processed
   customize by adding a math function:
     sqrt(),exp(),ln(),log(),abs(),sin(),cos(),tan(),asin(),acos(),atan(),
     atan2(y,x),random(x,y),normal(x,y),ceil(),floor(),round(),
     ramp(x,y),stagger(x,y),logfreq(x,y,z),stride(x,y,z),
     vdisplace(x,y),swiggle(x,y,z),cwiggle(x,y,z)
------------------------------------------------------------------------- */

int Variable::math_function(char *word, char *contents, Tree **tree,
                            Tree **treestack, int &ntreestack,
                            double *argstack, int &nargstack)
{
  // word not a match to any math function

  if (strcmp(word,"sqrt") && strcmp(word,"exp") &&
      strcmp(word,"ln") && strcmp(word,"log") &&
      strcmp(word,"abs") &&
      strcmp(word,"sin") && strcmp(word,"cos") &&
      strcmp(word,"tan") && strcmp(word,"asin") &&
      strcmp(word,"acos") && strcmp(word,"atan") &&
      strcmp(word,"atan2") && strcmp(word,"erf") &&
      strcmp(word,"random") &&
      strcmp(word,"normal") && strcmp(word,"ceil") &&
      strcmp(word,"floor") && strcmp(word,"round") &&
      strcmp(word,"ramp") && strcmp(word,"stagger") &&
      strcmp(word,"logfreq") && strcmp(word,"stride") &&
      strcmp(word,"vdisplace") &&
      strcmp(word,"swiggle") && strcmp(word,"cwiggle"))
    return 0;

  // parse contents for arg1,arg2,arg3 separated by commas
  // ptr1,ptr2 = location of 1st and 2nd comma, NULL if none

  char *arg1,*arg2,*arg3;
  char *ptr1,*ptr2;

  ptr1 = find_next_comma(contents);
  if (ptr1) {
    *ptr1 = '\0';
    ptr2 = find_next_comma(ptr1+1);
    if (ptr2) *ptr2 = '\0';
  } else ptr2 = NULL;

  int n = strlen(contents) + 1;
  arg1 = new char[n];
  strcpy(arg1,contents);
  int narg = 1;
  if (ptr1) {
    n = strlen(ptr1+1) + 1;
    arg2 = new char[n];
    strcpy(arg2,ptr1+1);
    narg = 2;
  } else arg2 = NULL;
  if (ptr2) {
    n = strlen(ptr2+1) + 1;
    arg3 = new char[n];
    strcpy(arg3,ptr2+1);
    narg = 3;
  } else arg3 = NULL;

  // evaluate args

  Tree *newtree;
  double value1,value2,value3;

  if (tree) {
    newtree = new Tree();
    Tree *argtree;
    if (narg == 1) {
      evaluate(arg1,&argtree);
      newtree->left = argtree;
      newtree->middle = newtree->right = NULL;
    } else if (narg == 2) {
      evaluate(arg1,&argtree);
      newtree->left = argtree;
      newtree->middle = NULL;
      evaluate(arg2,&argtree);
      newtree->right = argtree;
    } else if (narg == 3) {
      evaluate(arg1,&argtree);
      newtree->left = argtree;
      evaluate(arg2,&argtree);
      newtree->middle = argtree;
      evaluate(arg3,&argtree);
      newtree->right = argtree;
    }
    treestack[ntreestack++] = newtree;
  } else {
    if (narg == 1) {
      value1 = evaluate(arg1,NULL);
    } else if (narg == 2) {
      value1 = evaluate(arg1,NULL);
      value2 = evaluate(arg2,NULL);
    } else if (narg == 3) {
      value1 = evaluate(arg1,NULL);
      value2 = evaluate(arg2,NULL);
      value3 = evaluate(arg3,NULL);
    }
  }

  if (strcmp(word,"sqrt") == 0) {
    if (narg != 1)
      error->all(FLERR,"Invalid math function in variable formula");
    if (tree) newtree->type = SQRT;
    else {
      if (value1 < 0.0)
        error->all(FLERR,"Sqrt of negative value in variable formula");
      argstack[nargstack++] = sqrt(value1);
    }

  } else if (strcmp(word,"exp") == 0) {
    if (narg != 1)
      error->all(FLERR,"Invalid math function in variable formula");
    if (tree) newtree->type = EXP;
    else argstack[nargstack++] = exp(value1);
  } else if (strcmp(word,"ln") == 0) {
    if (narg != 1)
      error->all(FLERR,"Invalid math function in variable formula");
    if (tree) newtree->type = LN;
    else {
      if (value1 <= 0.0)
        error->all(FLERR,"Log of zero/negative value in variable formula");
      argstack[nargstack++] = log(value1);
    }
  } else if (strcmp(word,"log") == 0) {
    if (narg != 1)
      error->all(FLERR,"Invalid math function in variable formula");
    if (tree) newtree->type = LOG;
    else {
      if (value1 <= 0.0)
        error->all(FLERR,"Log of zero/negative value in variable formula");
      argstack[nargstack++] = log10(value1);
    }
  } else if (strcmp(word,"abs") == 0) {
    if (narg != 1)
      error->all(FLERR,"Invalid math function in variable formula");
    if (tree) newtree->type = ABS;
    else argstack[nargstack++] = fabs(value1);

  } else if (strcmp(word,"sin") == 0) {
    if (narg != 1)
      error->all(FLERR,"Invalid math function in variable formula");
    if (tree) newtree->type = SIN;
    else argstack[nargstack++] = sin(value1);
  } else if (strcmp(word,"cos") == 0) {
    if (narg != 1)
      error->all(FLERR,"Invalid math function in variable formula");
    if (tree) newtree->type = COS;
    else argstack[nargstack++] = cos(value1);
  } else if (strcmp(word,"tan") == 0) {
    if (narg != 1)
      error->all(FLERR,"Invalid math function in variable formula");
    if (tree) newtree->type = TAN;
    else argstack[nargstack++] = tan(value1);

  } else if (strcmp(word,"asin") == 0) {
    if (narg != 1)
      error->all(FLERR,"Invalid math function in variable formula");
    if (tree) newtree->type = ASIN;
    else {
      if (value1 < -1.0 || value1 > 1.0)
        error->all(FLERR,"Arcsin of invalid value in variable formula");
      argstack[nargstack++] = asin(value1);
    }
  } else if (strcmp(word,"acos") == 0) {
    if (narg != 1)
      error->all(FLERR,"Invalid math function in variable formula");
    if (tree) newtree->type = ACOS;
    else {
      if (value1 < -1.0 || value1 > 1.0)
        error->all(FLERR,"Arccos of invalid value in variable formula");
      argstack[nargstack++] = acos(value1);
    }
  } else if (strcmp(word,"atan") == 0) {
    if (narg != 1)
      error->all(FLERR,"Invalid math function in variable formula");
    if (tree) newtree->type = ATAN;
    else argstack[nargstack++] = atan(value1);
  } else if (strcmp(word,"atan2") == 0) {
    if (narg != 2)
      error->all(FLERR,"Invalid math function in variable formula");
    if (tree) newtree->type = ATAN2;
    else argstack[nargstack++] = atan2(value1,value2);
  } else if (strcmp(word,"erf") == 0) {
    if (narg != 1)
      error->all(FLERR,"Invalid math function in variable formula");
    if (tree) newtree->type = ERF;
    else argstack[nargstack++] = erf(value1);

  } else if (strcmp(word,"random") == 0) {
    if (narg != 2)
      error->all(FLERR,"Invalid math function in variable formula");
    if (tree) newtree->type = RANDOM;
    else {
      if (randomequal == NULL) {
        randomequal = new RanKnuth(update->ranmaster->uniform());
        double seed = update->ranmaster->uniform();
        randomequal->reset(seed,me,100);
      }
      argstack[nargstack++] = randomequal->uniform()*(value2-value1) + value1;
    }
  } else if (strcmp(word,"normal") == 0) {
    if (narg != 2)
      error->all(FLERR,"Invalid math function in variable formula");
    if (tree) newtree->type = NORMAL;
    else {
      if (value2 < 0.0)
        error->all(FLERR,"Invalid math function in variable formula");
      if (randomequal == NULL) {
        randomequal = new RanKnuth(update->ranmaster->uniform());
        double seed = update->ranmaster->uniform();
        randomequal->reset(seed,me,100);
      }
      argstack[nargstack++] = value1 + value2*randomequal->gaussian();
    }

  } else if (strcmp(word,"ceil") == 0) {
    if (narg != 1)
      error->all(FLERR,"Invalid math function in variable formula");
    if (tree) newtree->type = CEIL;
    else argstack[nargstack++] = ceil(value1);

  } else if (strcmp(word,"floor") == 0) {
    if (narg != 1)
      error->all(FLERR,"Invalid math function in variable formula");
    if (tree) newtree->type = FLOOR;
    else argstack[nargstack++] = floor(value1);

  } else if (strcmp(word,"round") == 0) {
    if (narg != 1)
      error->all(FLERR,"Invalid math function in variable formula");
    if (tree) newtree->type = ROUND;
    else argstack[nargstack++] = MYROUND(value1);

  } else if (strcmp(word,"ramp") == 0) {
    if (narg != 2)
      error->all(FLERR,"Invalid math function in variable formula");
    if (update->runflag == 0)
      error->all(FLERR,"Cannot use ramp in variable formula between runs");
    if (tree) newtree->type = RAMP;
    else {
      double delta = update->ntimestep - update->beginstep;
      if (delta != 0.0) delta /= update->endstep - update->beginstep;
      double value = value1 + delta*(value2-value1);
      argstack[nargstack++] = value;
    }

  } else if (strcmp(word,"stagger") == 0) {
    if (narg != 2)
      error->all(FLERR,"Invalid math function in variable formula");
    if (tree) newtree->type = STAGGER;
    else {
      int ivalue1 = static_cast<int> (value1);
      int ivalue2 = static_cast<int> (value2);
      if (ivalue1 <= 0 || ivalue2 <= 0 || ivalue1 <= ivalue2)
        error->all(FLERR,"Invalid math function in variable formula");
      int lower = update->ntimestep/ivalue1 * ivalue1;
      int delta = update->ntimestep - lower;
      double value;
      if (delta < ivalue2) value = lower+ivalue2;
      else value = lower+ivalue1;
      argstack[nargstack++] = value;
    }

  } else if (strcmp(word,"logfreq") == 0) {
    if (narg != 3)
      error->all(FLERR,"Invalid math function in variable formula");
    if (tree) newtree->type = LOGFREQ;
    else {
      int ivalue1 = static_cast<int> (value1);
      int ivalue2 = static_cast<int> (value2);
      int ivalue3 = static_cast<int> (value3);
      if (ivalue1 <= 0 || ivalue2 <= 0 || ivalue3 <= 0 || ivalue2 >= ivalue3)
        error->all(FLERR,"Invalid math function in variable formula");
      double value;
      if (update->ntimestep < ivalue1) value = ivalue1;
      else {
        int lower = ivalue1;
        while (update->ntimestep >= ivalue3*lower) lower *= ivalue3;
        int multiple = update->ntimestep/lower;
        if (multiple < ivalue2) value = (multiple+1)*lower;
        else value = lower*ivalue3;
      }
      argstack[nargstack++] = value;
    }

  } else if (strcmp(word,"stride") == 0) {
    if (narg != 3)
      error->all(FLERR,"Invalid math function in variable formula");
    if (tree) newtree->type = STRIDE;
    else {
      int ivalue1 = static_cast<int> (value1);
      int ivalue2 = static_cast<int> (value2);
      int ivalue3 = static_cast<int> (value3);
      if (ivalue1 < 0 || ivalue2 < 0 || ivalue3 <= 0 || ivalue1 > ivalue2)
        error->one(FLERR,"Invalid math function in variable formula");
      double value;
      if (update->ntimestep < ivalue1) value = ivalue1;
      else if (update->ntimestep < ivalue2) {
        int offset = update->ntimestep - ivalue1;
        value = ivalue1 + (offset/ivalue3)*ivalue3 + ivalue3;
        if (value > ivalue2) value = 9.0e18;
      } else value = 9.0e18;
      argstack[nargstack++] = value;
    }

  } else if (strcmp(word,"vdisplace") == 0) {
    if (narg != 2)
      error->all(FLERR,"Invalid math function in variable formula");
    if (update->runflag == 0)
      error->all(FLERR,"Cannot use vdisplace in variable formula between runs");
    if (tree) newtree->type = VDISPLACE;
    else {
      double delta = update->ntimestep - update->beginstep;
      double value = value1 + value2*delta*update->dt;
      argstack[nargstack++] = value;
    }

  } else if (strcmp(word,"swiggle") == 0) {
    if (narg != 3)
      error->all(FLERR,"Invalid math function in variable formula");
    if (update->runflag == 0)
      error->all(FLERR,"Cannot use swiggle in variable formula between runs");
    if (tree) newtree->type = CWIGGLE;
    else {
      if (value3 == 0.0)
        error->all(FLERR,"Invalid math function in variable formula");
      double delta = update->ntimestep - update->beginstep;
      double omega = 2.0*MY_PI/value3;
      double value = value1 + value2*sin(omega*delta*update->dt);
      argstack[nargstack++] = value;
    }

  } else if (strcmp(word,"cwiggle") == 0) {
    if (narg != 3)
      error->all(FLERR,"Invalid math function in variable formula");
    if (update->runflag == 0)
      error->all(FLERR,"Cannot use cwiggle in variable formula between runs");
    if (tree) newtree->type = CWIGGLE;
    else {
      if (value3 == 0.0)
        error->all(FLERR,"Invalid math function in variable formula");
      double delta = update->ntimestep - update->beginstep;
      double omega = 2.0*MY_PI/value3;
      double value = value1 + value2*(1.0-cos(omega*delta*update->dt));
      argstack[nargstack++] = value;
    }
  }

  delete [] arg1;
  delete [] arg2;
  delete [] arg3;

  return 1;
}

/* ----------------------------------------------------------------------
   process a special function in formula
   push result onto tree or arg stack
   word = special function
   contents = str between parentheses with one,two,three args
   return 0 if not a match, 1 if successfully processed
   customize by adding a special function:
     sum(x),min(x),max(x),ave(x),trap(x),slope(x),next(x)
------------------------------------------------------------------------- */

int Variable::special_function(char *word, char *contents, Tree **tree,
                               Tree **treestack, int &ntreestack,
                               double *argstack, int &nargstack)
{
  double value,xvalue,sx,sy,sxx,sxy;

  // word not a match to any special function

  if (strcmp(word,"sum") && strcmp(word,"min") && strcmp(word,"max") &&
      strcmp(word,"ave") && strcmp(word,"trap") && strcmp(word,"slope") &&
      strcmp(word,"next"))
    return 0;

  // parse contents for arg1,arg2,arg3 separated by commas
  // ptr1,ptr2 = location of 1st and 2nd comma, NULL if none

  char *arg1,*arg2,*arg3;
  char *ptr1,*ptr2;

  ptr1 = find_next_comma(contents);
  if (ptr1) {
    *ptr1 = '\0';
    ptr2 = find_next_comma(ptr1+1);
    if (ptr2) *ptr2 = '\0';
  } else ptr2 = NULL;

  int n = strlen(contents) + 1;
  arg1 = new char[n];
  strcpy(arg1,contents);
  int narg = 1;
  if (ptr1) {
    n = strlen(ptr1+1) + 1;
    arg2 = new char[n];
    strcpy(arg2,ptr1+1);
    narg = 2;
  } else arg2 = NULL;
  if (ptr2) {
    n = strlen(ptr2+1) + 1;
    arg3 = new char[n];
    strcpy(arg3,ptr2+1);
    narg = 3;
  } else arg3 = NULL;

  // special functions that operate on global vectors

  if (strcmp(word,"sum") == 0 || strcmp(word,"min") == 0 ||
      strcmp(word,"max") == 0 || strcmp(word,"ave") == 0 ||
      strcmp(word,"trap") == 0 || strcmp(word,"slope") == 0) {

    int method;
    if (strcmp(word,"sum") == 0) method = SUM;
    else if (strcmp(word,"min") == 0) method = XMIN;
    else if (strcmp(word,"max") == 0) method = XMAX;
    else if (strcmp(word,"ave") == 0) method = AVE;
    else if (strcmp(word,"trap") == 0) method = TRAP;
    else if (strcmp(word,"slope") == 0) method = SLOPE;

    if (narg != 1)
      error->all(FLERR,"Invalid special function in variable formula");

    Compute *compute = NULL;
    Fix *fix = NULL;
    int index,nvec,nstride;

    if (strstr(arg1,"c_") == arg1) {
      ptr1 = strchr(arg1,'[');
      if (ptr1) {
        ptr2 = ptr1;
        index = int_between_brackets(ptr2,0);
        *ptr1 = '\0';
      } else index = 0;

      int icompute = modify->find_compute(&arg1[2]);
      if (icompute < 0)
        error->all(FLERR,"Invalid compute ID in variable formula");
      compute = modify->compute[icompute];
      if (index == 0 && compute->vector_flag) {
        if (update->runflag == 0) {
          if (compute->invoked_vector != update->ntimestep)
            error->all(FLERR,
                       "Compute used in variable between runs is not current");
        } else if (!(compute->invoked_flag & INVOKED_VECTOR)) {
          compute->compute_vector();
          compute->invoked_flag |= INVOKED_VECTOR;
        }
        nvec = compute->size_vector;
        nstride = 1;
      } else if (index && compute->array_flag) {
        if (index > compute->size_array_cols)
          error->all(FLERR,"Variable formula compute array "
                     "is accessed out-of-range");
        if (update->runflag == 0) {
          if (compute->invoked_array != update->ntimestep)
            error->all(FLERR,
                       "Compute used in variable between runs is not current");
        } else if (!(compute->invoked_flag & INVOKED_ARRAY)) {
          compute->compute_array();
          compute->invoked_flag |= INVOKED_ARRAY;
        }
        nvec = compute->size_array_rows;
        nstride = compute->size_array_cols;
      } else error->all(FLERR,"Mismatched compute in variable formula");

    } else if (strstr(arg1,"f_") == arg1) {
      ptr1 = strchr(arg1,'[');
      if (ptr1) {
        ptr2 = ptr1;
        index = int_between_brackets(ptr2,0);
        *ptr1 = '\0';
      } else index = 0;

      int ifix = modify->find_fix(&arg1[2]);
      if (ifix < 0) error->all(FLERR,"Invalid fix ID in variable formula");
      fix = modify->fix[ifix];
      if (index == 0 && fix->vector_flag) {
        if (update->runflag > 0 && update->ntimestep % fix->global_freq)
          error->all(FLERR,"Fix in variable not computed at compatible time");
        nvec = fix->size_vector;
        nstride = 1;
      } else if (index && fix->array_flag) {
        if (index > fix->size_array_cols)
          error->all(FLERR,
                     "Variable formula fix array is accessed out-of-range");
        if (update->runflag > 0 && update->ntimestep % fix->global_freq)
          error->all(FLERR,"Fix in variable not computed at compatible time");
        nvec = fix->size_array_rows;
        nstride = fix->size_array_cols;
      } else error->all(FLERR,"Mismatched fix in variable formula");

    } else error->all(FLERR,"Invalid special function in variable formula");

    value = 0.0;
    if (method == SLOPE) sx = sy = sxx = sxy = 0.0;
    if (method == XMIN) value = BIG;
    if (method == XMAX) value = -BIG;

    if (compute) {
      double *vec;
      if (index) {
        if (compute->array) vec = &compute->array[0][index-1];
        else vec = NULL;
      } else vec = compute->vector;

      int j = 0;
      for (int i = 0; i < nvec; i++) {
        if (method == SUM) value += vec[j];
        else if (method == XMIN) value = MIN(value,vec[j]);
        else if (method == XMAX) value = MAX(value,vec[j]);
        else if (method == AVE) value += vec[j];
        else if (method == TRAP) value += vec[j];
        else if (method == SLOPE) {
          if (nvec > 1) xvalue = (double) i / (nvec-1);
          else xvalue = 0.0;
          sx += xvalue;
          sy += vec[j];
          sxx += xvalue*xvalue;
          sxy += xvalue*vec[j];
        }
        j += nstride;
      }
      if (method == TRAP) value -= 0.5*vec[0] + 0.5*vec[nvec-1];
    }

    if (fix) {
      double one;
      for (int i = 0; i < nvec; i++) {
        if (index) one = fix->compute_array(i,index-1);
        else one = fix->compute_vector(i);
        if (method == SUM) value += one;
        else if (method == XMIN) value = MIN(value,one);
        else if (method == XMAX) value = MAX(value,one);
        else if (method == AVE) value += one;
        else if (method == TRAP) value += one;
        else if (method == SLOPE) {
          if (nvec > 1) xvalue = (double) i / (nvec-1);
          else xvalue = 0.0;
          sx += xvalue;
          sy += one;
          sxx += xvalue*xvalue;
          sxy += xvalue*one;
        }
      }
      if (method == TRAP) {
        if (index) value -= 0.5*fix->compute_array(0,index-1) +
                     0.5*fix->compute_array(nvec-1,index-1);
        else value -= 0.5*fix->compute_vector(0) +
               0.5*fix->compute_vector(nvec-1);
      }
    }

    if (method == AVE) value /= nvec;

    if (method == SLOPE) {
      double numerator = sxy - sx*sy;
      double denominator = sxx - sx*sx;
      if (denominator != 0.0) value = numerator/denominator / nvec;
      else value = BIG;
    }

    // save value in tree or on argstack

    if (tree) {
      Tree *newtree = new Tree();
      newtree->type = VALUE;
      newtree->value = value;
      newtree->left = newtree->middle = newtree->right = NULL;
      treestack[ntreestack++] = newtree;
    } else argstack[nargstack++] = value;

  // special function for file-style variable

  } else if (strcmp(word,"next") == 0) {
    if (narg != 1)
      error->all(FLERR,"Invalid special function in variable formula");

    int ivar = find(arg1);
    if (ivar == -1)
      error->all(FLERR,"Variable ID in variable formula does not exist");

    // SCALARFILE has single current value, read next one
    // save value in tree or on argstack

    if (style[ivar] == SCALARFILE) {
      double value = atof(data[ivar][0]);
      int done = reader[ivar]->read_scalar(data[ivar][0]);
      if (done) remove(ivar);

      if (tree) {
        Tree *newtree = new Tree();
        newtree->type = VALUE;
        newtree->value = value;
        newtree->left = newtree->middle = newtree->right = NULL;
        treestack[ntreestack++] = newtree;
      } else argstack[nargstack++] = value;

    } else error->all(FLERR,"Invalid variable style in special function next");
  }

  delete [] arg1;
  delete [] arg2;
  delete [] arg3;

  return 1;
}

/* ----------------------------------------------------------------------
   check if word matches a particle vector
   return 1 if yes, else 0
   customize by adding a particle vector:
     x,y,z,vx,vy,vz,type,mass,q,mu
------------------------------------------------------------------------- */

int Variable::is_particle_vector(char *word)
{
  if (strcmp(word,"x") == 0) return 1;
  if (strcmp(word,"y") == 0) return 1;
  if (strcmp(word,"z") == 0) return 1;
  if (strcmp(word,"vx") == 0) return 1;
  if (strcmp(word,"vy") == 0) return 1;
  if (strcmp(word,"vz") == 0) return 1;
  if (strcmp(word,"type") == 0) return 1;
  if (strcmp(word,"mass") == 0) return 1;
  if (strcmp(word,"q") == 0) return 1;
  if (strcmp(word,"mu") == 0) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   process a particle vector in formula
   push result onto tree
   word = particle vector
   customize by adding a particle vector:
     x,y,z,vx,vy,vz,type,mass,q,mu
------------------------------------------------------------------------- */

void Variable::particle_vector(char *word, Tree **tree,
                               Tree **treestack, int &ntreestack)
{
  if (tree == NULL || treestyle != PARTICLE)
    error->all(FLERR,"Particle vector in non particle-style variable formula");

  Particle::OnePart *particles = particle->particles;
  Particle::Species *species = particle->species;

  Tree *newtree = new Tree();
  newtree->type = PARTARRAYDOUBLE;
  newtree->nstride = sizeof(Particle::OnePart);
  newtree->left = newtree->middle = newtree->right = NULL;
  treestack[ntreestack++] = newtree;

  if (strcmp(word,"x") == 0)
    newtree->carray = (char *) &particles[0].x[0];
  else if (strcmp(word,"y") == 0)
    newtree->carray = (char *) &particles[0].x[1];
  else if (strcmp(word,"z") == 0)
    newtree->carray = (char *) &particles[0].x[2];
  else if (strcmp(word,"vx") == 0)
    newtree->carray = (char *) &particles[0].v[0];
  else if (strcmp(word,"vy") == 0)
    newtree->carray = (char *) &particles[0].v[1];
  else if (strcmp(word,"vz") == 0)
    newtree->carray = (char *) &particles[0].v[2];

  else if (strcmp(word,"type") == 0) {
    newtree->type = PARTARRAYINT;
    newtree->carray = (char *) &particles[0].ispecies;
  } else if (strcmp(word,"mass") == 0) {
    newtree->type = SPECARRAY;
    newtree->nstride = sizeof(Particle::Species);
    newtree->carray = (char *) &species[0].mass;
  } else if (strcmp(word,"q") == 0) {
    newtree->type = SPECARRAY;
    newtree->nstride = sizeof(Particle::Species);
    newtree->carray = (char *) &species[0].charge;
  } else if (strcmp(word,"mu") == 0) {
    newtree->type = SPECARRAY;
    newtree->nstride = sizeof(Particle::Species);
    newtree->carray = (char *) &species[0].magmoment;
  }

  if ((bigint)particle->nlocal*newtree->nstride > MAXSMALLINT)
    error->all(FLERR,"Too many particles per processor for particle-style variable");
}

/* ----------------------------------------------------------------------
   check if word matches a grid vector
   return 1 if yes, else 0
   customize by adding a grid vector:
     cxlo,cxhi,cylo,cyhi,czlo,czhi
------------------------------------------------------------------------- */

int Variable::is_grid_vector(char *word)
{
  if (strcmp(word,"cxlo") == 0) return 1;
  if (strcmp(word,"cxhi") == 0) return 1;
  if (strcmp(word,"cylo") == 0) return 1;
  if (strcmp(word,"cyhi") == 0) return 1;
  if (strcmp(word,"czlo") == 0) return 1;
  if (strcmp(word,"czhi") == 0) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   process a grid vector in formula
   push result onto tree
   word = grid vector
   customize by adding a grid vector:
     cxlo,cxhi,cylo,cyhi,czlo,czhi
------------------------------------------------------------------------- */

void Variable::grid_vector(char *word, Tree **tree,
                           Tree **treestack, int &ntreestack)
{
  if (tree == NULL || treestyle != GRID)
    error->all(FLERR,"Grid vector in non grid-style variable formula");

  Grid::ChildCell *cells = grid->cells;

  Tree *newtree = new Tree();
  newtree->type = PARTARRAYDOUBLE;
  newtree->nstride = sizeof(Grid::ChildCell);
  newtree->left = newtree->middle = newtree->right = NULL;
  treestack[ntreestack++] = newtree;

  if (strcmp(word,"cxlo") == 0)
    newtree->carray = (char *) &cells[0].lo[0];
  else if (strcmp(word,"cxhi") == 0)
    newtree->carray = (char *) &cells[0].hi[0];
  else if (strcmp(word,"cylo") == 0)
    newtree->carray = (char *) &cells[0].lo[1];
  else if (strcmp(word,"cyhi") == 0)
    newtree->carray = (char *) &cells[0].hi[1];
  else if (strcmp(word,"czlo") == 0)
    newtree->carray = (char *) &cells[0].lo[2];
  else if (strcmp(word,"czhi") == 0)
    newtree->carray = (char *) &cells[0].hi[2];

  if ((bigint)grid->nlocal*newtree->nstride > MAXSMALLINT)
    error->all(FLERR,"Too many grid cells per processor for grid-style variable");
}

/* ----------------------------------------------------------------------
   check if word matches a constant
   return 1 if yes, else 0
   customize by adding a constant: PI
------------------------------------------------------------------------- */

int Variable::is_constant(char *word)
{
  if (strcmp(word,"PI") == 0) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   process a constant in formula
   customize by adding a constant: PI
------------------------------------------------------------------------- */

double Variable::constant(char *word)
{
  if (strcmp(word,"PI") == 0) return MY_PI;
  return 0.0;
}

/* ----------------------------------------------------------------------
   find next comma in str
   skip commas inside one or more nested parenthesis
   only return ptr to comma at level 0, else NULL if not found
------------------------------------------------------------------------- */

char *Variable::find_next_comma(char *str)
{
  int level = 0;
  for (char *p = str; *p; ++p) {
    if ('(' == *p) level++;
    else if (')' == *p) level--;
    else if (',' == *p && !level) return p;
  }
  return NULL;
}

/* ----------------------------------------------------------------------
   copy of cvec in local storage
   return local ptr to copied vec
------------------------------------------------------------------------- */

double *Variable::add_storage(double *cvec)
{
  if (nvec_storage == maxvec_storage) {
    maxvec_storage++;
    vec_storage = (double **)
      memory->srealloc(vec_storage,maxvec_storage*sizeof(double *),
                       "variable:vec_storage");
    maxlen_storage = (int *)
      memory->srealloc(maxlen_storage,maxvec_storage*sizeof(int),
                       "variable:maxlen_storage");
    vec_storage[nvec_storage] = NULL;
    maxlen_storage[nvec_storage] = 0;
  }

  int n = grid->nlocal;
  if (maxlen_storage[nvec_storage] < n) {
    maxlen_storage[nvec_storage] = n;
    memory->destroy(vec_storage[nvec_storage]);
    memory->create(vec_storage[nvec_storage],n,"variable:vec_storage");
  }

  memcpy(vec_storage[nvec_storage],cvec,n*sizeof(double));
  nvec_storage++;

  return vec_storage[nvec_storage-1];
}

/* ----------------------------------------------------------------------
   debug routine for printing formula tree recursively
------------------------------------------------------------------------- */

void Variable::print_tree(Tree *tree, int level)
{
  printf("TREE %d: %d %g\n",level,tree->type,tree->value);
  if (tree->left) print_tree(tree->left,level+1);
  if (tree->middle) print_tree(tree->middle,level+1);
  if (tree->right) print_tree(tree->right,level+1);
  return;
}

/* ----------------------------------------------------------------------
   recursive evaluation of string str
   called from "if" command in input script
   str is a boolean expression containing one or more items:
     number = 0.0, -5.45, 2.8e-4, ...
     math operation = (),x==y,x!=y,x<y,x<=y,x>y,x>=y,x&&y,x||y
------------------------------------------------------------------------- */

double Variable::evaluate_boolean(char *str)
{
  int op,opprevious,flag1,flag2;
  double value1,value2;
  char onechar;
  char *str1,*str2;

  struct Arg {
    int flag;          // 0 for numeric value, 1 for string
    double value;      // stored numeric value
    char *str;         // stored string
  };

  Arg argstack[MAXLEVEL];
  int opstack[MAXLEVEL];
  int nargstack = 0;
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
      if (expect == OP)
        error->all(FLERR,"Invalid Boolean syntax in if command");
      expect = OP;

      char *contents;
      i = find_matching_paren(str,i,contents);
      i++;

      // evaluate contents and push on stack

      argstack[nargstack].value = evaluate_boolean(contents);
      argstack[nargstack].flag = 0;
      nargstack++;

      delete [] contents;

    // ----------------
    // number: push value onto stack
    // ----------------

    } else if (isdigit(onechar) || onechar == '.' || onechar == '-') {
      if (expect == OP)
        error->all(FLERR,"Invalid Boolean syntax in if command");
      expect = OP;

      // set I to end of number, including scientific notation

      int istart = i++;
      while (isdigit(str[i]) || str[i] == '.') i++;
      if (str[i] == 'e' || str[i] == 'E') {
        i++;
        if (str[i] == '+' || str[i] == '-') i++;
        while (isdigit(str[i])) i++;
      }

      onechar = str[i];
      str[i] = '\0';
      argstack[nargstack].value = atof(&str[istart]);
      str[i] = onechar;

      argstack[nargstack++].flag = 0;

    // ----------------
    // string: push string onto stack
    // ----------------

    } else if (isalpha(onechar)) {
      if (expect == OP)
        error->all(FLERR,"Invalid Boolean syntax in if command");
      expect = OP;

      // set I to end of string

      int istart = i++;
      while (isalnum(str[i]) || str[i] == '_') i++;

      int n = i - istart + 1;
      argstack[nargstack].str = new char[n];
      onechar = str[i];
      str[i] = '\0';
      strcpy(argstack[nargstack].str,&str[istart]);
      str[i] = onechar;

      argstack[nargstack++].flag = 1;

    // ----------------
    // Boolean operator, including end-of-string
    // ----------------

    } else if (strchr("<>=!&|\0",onechar)) {
      if (onechar == '=') {
        if (str[i+1] != '=')
          error->all(FLERR,"Invalid Boolean syntax in if command");
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
          error->all(FLERR,"Invalid Boolean syntax in if command");
        op = AND;
        i++;
      } else if (onechar == '|') {
        if (str[i+1] != '|')
          error->all(FLERR,"Invalid Boolean syntax in if command");
        op = OR;
        i++;
      } else op = DONE;

      i++;

      if (op == NOT && expect == ARG) {
        opstack[nopstack++] = op;
        continue;
      }

      if (expect == ARG)
        error->all(FLERR,"Invalid Boolean syntax in if command");
      expect = ARG;

      // evaluate stack as deep as possible while respecting precedence
      // before pushing current op onto stack

      while (nopstack && precedence[opstack[nopstack-1]] >= precedence[op]) {
        opprevious = opstack[--nopstack];

        nargstack--;
        flag2 = argstack[nargstack].flag;
        value2 = argstack[nargstack].value;
        str2 = argstack[nargstack].str;
        if (opprevious != NOT) {
          nargstack--;
          flag1 = argstack[nargstack].flag;
          value1 = argstack[nargstack].value;
          str1 = argstack[nargstack].str;
        }

        if (opprevious == NOT) {
          if (flag2) error->all(FLERR,"Invalid Boolean syntax in if command");
          if (value2 == 0.0) argstack[nargstack].value = 1.0;
          else argstack[nargstack].value = 0.0;
        } else if (opprevious == EQ) {
          if (flag1 != flag2)
            error->all(FLERR,"Invalid Boolean syntax in if command");
          if (flag2 == 0) {
            if (value1 == value2) argstack[nargstack].value = 1.0;
            else argstack[nargstack].value = 0.0;
          } else {
            if (strcmp(str1,str2) == 0) argstack[nargstack].value = 1.0;
            else argstack[nargstack].value = 0.0;
            delete [] str1;
            delete [] str2;
          }
        } else if (opprevious == NE) {
          if (flag1 != flag2)
            error->all(FLERR,"Invalid Boolean syntax in if command");
          if (flag2 == 0) {
            if (value1 != value2) argstack[nargstack].value = 1.0;
            else argstack[nargstack].value = 0.0;
          } else {
            if (strcmp(str1,str2) != 0) argstack[nargstack].value = 1.0;
            else argstack[nargstack].value = 0.0;
            delete [] str1;
            delete [] str2;
          }
        } else if (opprevious == LT) {
          if (flag2) error->all(FLERR,"Invalid Boolean syntax in if command");
          if (value1 < value2) argstack[nargstack].value = 1.0;
          else argstack[nargstack].value = 0.0;
        } else if (opprevious == LE) {
          if (flag2) error->all(FLERR,"Invalid Boolean syntax in if command");
          if (value1 <= value2) argstack[nargstack].value = 1.0;
          else argstack[nargstack].value = 0.0;
        } else if (opprevious == GT) {
          if (flag2) error->all(FLERR,"Invalid Boolean syntax in if command");
          if (value1 > value2) argstack[nargstack].value = 1.0;
          else argstack[nargstack].value = 0.0;
        } else if (opprevious == GE) {
          if (flag2) error->all(FLERR,"Invalid Boolean syntax in if command");
          if (value1 >= value2) argstack[nargstack].value = 1.0;
          else argstack[nargstack].value = 0.0;
        } else if (opprevious == AND) {
          if (flag2) error->all(FLERR,"Invalid Boolean syntax in if command");
          if (value1 != 0.0 && value2 != 0.0) argstack[nargstack].value = 1.0;
          else argstack[nargstack].value = 0.0;
        } else if (opprevious == OR) {
          if (flag2) error->all(FLERR,"Invalid Boolean syntax in if command");
          if (value1 != 0.0 || value2 != 0.0) argstack[nargstack].value = 1.0;
          else argstack[nargstack].value = 0.0;
        }

        argstack[nargstack++].flag = 0;
      }

      // if end-of-string, break out of entire formula evaluation loop

      if (op == DONE) break;

      // push current operation onto stack

      opstack[nopstack++] = op;

    } else error->all(FLERR,"Invalid Boolean syntax in if command");
  }

  if (nopstack) error->all(FLERR,"Invalid Boolean syntax in if command");
  if (nargstack != 1) error->all(FLERR,"Invalid Boolean syntax in if command");
  return argstack[0].value;
}

/* ----------------------------------------------------------------------
   class to read variable values from a file
   for flag = SCALARFILE, reads one value per line
------------------------------------------------------------------------- */

VarReader::VarReader(SPARTA *sparta, char *, char *file, int flag) :
  Pointers(sparta)
{
  me = comm->me;
  style = flag;

  if (me == 0) {
    fp = fopen(file,"r");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open file variable file %s",file);
      error->one(FLERR,str);
    }
  } else fp = NULL;
}

/* ---------------------------------------------------------------------- */

VarReader::~VarReader()
{
  if (me == 0) fclose(fp);
}

/* ----------------------------------------------------------------------
   read for SCALARFILE style
   read next value from file into str for file-style variable
   strip comments, skip blank lines
   return 0 if successful, 1 if end-of-file
------------------------------------------------------------------------- */

int VarReader::read_scalar(char *str)
{
  int n;
  char *ptr;

  // read one string from file

  if (me == 0) {
    while (1) {
      if (fgets(str,MAXLINE,fp) == NULL) n = 0;
      else n = strlen(str);
      if (n == 0) break;                                 // end of file
      str[n-1] = '\0';                                   // strip newline
      if ((ptr = strchr(str,'#'))) *ptr = '\0';          // strip comment
      if (strtok(str," \t\n\r\f") == NULL) continue;     // skip if blank
      n = strlen(str) + 1;
      break;
    }
  }

  MPI_Bcast(&n,1,MPI_INT,0,world);
  if (n == 0) return 1;
  MPI_Bcast(str,n,MPI_CHAR,0,world);
  return 0;
}

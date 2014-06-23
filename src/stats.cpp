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

#include "spatype.h"
#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "stats.h"
#include "update.h"
#include "particle.h"
#include "collide.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "input.h"
#include "variable.h"
#include "output.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

// customize a new keyword by adding to this list:

// step,elapsed,elaplong,dt,cpu,tpcpu,spcpu
// np,ntouch,ncomm,nbound,nexit,nscoll,nscheck,ncoll,nattempt,nreact,
// npave,ntouchave,ncommave,nboundave,nexitave,nscollave,nscheckave,
// ncollave,nattemptave,nreactave,
// vol,lx,ly,lz,xlo,xhi,ylo,yhi,zlo,zhi

enum{INT,FLOAT,BIGINT};
enum{SCALAR,VECTOR,ARRAY};

#define INVOKED_SCALAR 1
#define INVOKED_VECTOR 2
#define INVOKED_ARRAY 4

#define MAXLINE 8192               // make this 4x longer than Input::MAXLINE
#define DELTA 8

/* ---------------------------------------------------------------------- */

Stats::Stats(SPARTA *sparta) : Pointers(sparta)
{
  MPI_Comm_rank(world,&me);

  line = new char[MAXLINE];

  keyword = NULL;
  vfunc = NULL;
  vtype = NULL;

  format = NULL;
  format_user = NULL;

  field2index = NULL;
  argindex1 = NULL;
  argindex2 = NULL;

  // default args

  char **arg = new char*[3];
  arg[0] = (char *) "step";
  arg[1] = (char *) "cpu";
  arg[2] = (char *) "np";

  nfield = 3;
  allocate();
  set_fields(3,arg);

  delete [] arg;

  // stats_modify defaults

  flushflag = 0;

  // format strings

  char *bigint_format = (char *) BIGINT_FORMAT;

  format_float_one_def = (char *) "%12.8g";
  format_int_one_def = (char *) "%8d";
  sprintf(format_bigint_one_def,"%%8%s",&bigint_format[1]);

  format_float_user = NULL;
  format_int_user = NULL;
  format_bigint_user = NULL;
}

/* ---------------------------------------------------------------------- */

Stats::~Stats()
{
  delete [] line;
  deallocate();

  // format strings

  delete [] format_float_user;
  delete [] format_int_user;
  delete [] format_bigint_user;
}

/* ---------------------------------------------------------------------- */

void Stats::init()
{
  // set format string for each field
  // add trailing '/n' to last value

  char *ptr;
  for (int i = 0; i < nfield; i++) {
    format[i][0] = '\0';

    if (format_user[i]) ptr = format_user[i];
    else if (vtype[i] == FLOAT) {
      if (format_float_user) ptr = format_float_user;
      else ptr = format_float_one_def;
    } else if (vtype[i] == INT) {
      if (format_int_user) ptr = format_int_user;
      else ptr = format_int_one_def;
    } else if (vtype[i] == BIGINT) {
      if (format_bigint_user) ptr = format_bigint_user;
      else ptr = format_bigint_one_def;
    }

    int n = strlen(format[i]);
    sprintf(&format[i][n],"%s ",ptr);

    if (i == nfield-1) strcat(format[i],"\n");
  }

  // find current ptr for each Compute ID

  int icompute;
  for (int i = 0; i < ncompute; i++) {
    icompute = modify->find_compute(id_compute[i]);
    if (icompute < 0) error->all(FLERR,"Could not find stats compute ID");
    computes[i] = modify->compute[icompute];
  }

  // find current ptr for each Fix ID
  // check that fix frequency is acceptable with stats output frequency

  int ifix;
  for (int i = 0; i < nfix; i++) {
    ifix = modify->find_fix(id_fix[i]);
    if (ifix < 0) error->all(FLERR,"Could not find stats fix ID");
    fixes[i] = modify->fix[ifix];
    if (output->stats_every % fixes[i]->global_freq)
      error->all(FLERR,"Stats and fix not computed at compatible times");
  }

  // find current ptr for each Variable ID

  int ivariable;
  for (int i = 0; i < nvariable; i++) {
    ivariable = input->variable->find(id_variable[i]);
    if (ivariable < 0) 
      error->all(FLERR,"Could not find stats variable name");
    variables[i] = ivariable;
  }
}

/* ---------------------------------------------------------------------- */

void Stats::header()
{
  int loc = 0;
  for (int i = 0; i < nfield; i++)
    loc += sprintf(&line[loc],"%s ",keyword[i]);
  sprintf(&line[loc],"\n");
  
  if (me == 0) {
    if (screen) fprintf(screen,line);
    if (logfile) fprintf(logfile,line);
  }
}

/* ---------------------------------------------------------------------- */

void Stats::compute(int flag)
{
  int i;

  firststep = flag;
  bigint ntimestep = update->ntimestep;

  // invoke Compute methods needed for stats keywords

  for (i = 0; i < ncompute; i++)
    if (compute_which[i] == SCALAR) {
      if (!(computes[i]->invoked_flag & INVOKED_SCALAR)) {
	computes[i]->compute_scalar();
	computes[i]->invoked_flag |= INVOKED_SCALAR;
      }
    } else if (compute_which[i] == VECTOR) {
      if (!(computes[i]->invoked_flag & INVOKED_VECTOR)) {
	computes[i]->compute_vector();
	computes[i]->invoked_flag |= INVOKED_VECTOR;
      }
    } else if (compute_which[i] == ARRAY) {
      if (!(computes[i]->invoked_flag & INVOKED_ARRAY)) {
	computes[i]->compute_array();
	computes[i]->invoked_flag |= INVOKED_ARRAY;
      }
    }

  // add each stat value to line with its specific format

  int loc = 0;
  for (ifield = 0; ifield < nfield; ifield++) {
    (this->*vfunc[ifield])();
    if (vtype[ifield] == FLOAT)
      loc += sprintf(&line[loc],format[ifield],dvalue);
    else if (vtype[ifield] == INT) 
      loc += sprintf(&line[loc],format[ifield],ivalue);
    else if (vtype[ifield] == BIGINT) {
      loc += sprintf(&line[loc],format[ifield],bivalue);
    }
  }

  // print line to screen and logfile

  if (me == 0) {
    if (screen) fprintf(screen,line);
    if (logfile) {
      fprintf(logfile,line);
      if (flushflag) fflush(logfile);
    }
  }
}

/* ----------------------------------------------------------------------
   modify stats parameters
------------------------------------------------------------------------- */

void Stats::modify_params(int narg, char **arg)
{
  if (narg == 0) error->all(FLERR,"Illegal stats_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal stats_modify command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
	delete [] output->var_stats;
	int n = strlen(&arg[iarg+1][2]) + 1;
	output->var_stats = new char[n];
	strcpy(output->var_stats,&arg[iarg+1][2]);
      } else error->all(FLERR,"Illegal stats_modify command");
      output->stats_every = 0;
      iarg += 2;

    } else if (strcmp(arg[iarg],"flush") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal stats_modify command");
      if (strcmp(arg[iarg+1],"no") == 0) flushflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) flushflag = 1;
      else error->all(FLERR,"Illegal stats_modify command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"format") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal stats_modify command");
      if (strcmp(arg[iarg+1],"int") == 0) {
	if (format_int_user) delete [] format_int_user;
	int n = strlen(arg[iarg+2]) + 1;
	format_int_user = new char[n];
	strcpy(format_int_user,arg[iarg+2]);
	if (format_bigint_user) delete [] format_bigint_user;
	n = strlen(format_int_user) + 3;
	format_bigint_user = new char[n];
	char *ptr = strchr(format_int_user,'d');
	if (ptr == NULL) 
	  error->all(FLERR,
		     "Stats_modify int format does not contain d character");
	*ptr = '\0';
	sprintf(format_bigint_user,"%s%s%s",format_int_user,
		BIGINT_FORMAT,ptr+1);
	*ptr = 'd';
      } else if (strcmp(arg[iarg+1],"float") == 0) {
	if (format_float_user) delete [] format_float_user;
	int n = strlen(arg[iarg+2]) + 1;
	format_float_user = new char[n];
	strcpy(format_float_user,arg[iarg+2]);
      } else {
	int i = atoi(arg[iarg+1]) - 1;
	if (i < 0 || i >= nfield)
	  error->all(FLERR,"Illegal stats_modify command");
	if (format_user[i]) delete [] format_user[i];
	int n = strlen(arg[iarg+2]) + 1;
	format_user[i] = new char[n];
	strcpy(format_user[i],arg[iarg+2]);
      }
      iarg += 3;

    } else error->all(FLERR,"Illegal stats_modify command");
  }
}

/* ----------------------------------------------------------------------
   allocate all per-field memory
------------------------------------------------------------------------- */

void Stats::allocate()
{
  int n = nfield;

  keyword = new char*[n];
  for (int i = 0; i < n; i++) keyword[i] = new char[32];
  vfunc = new FnPtr[n];
  vtype = new int[n];

  format = new char*[n];
  for (int i = 0; i < n; i++) format[i] = new char[32];
  format_user = new char*[n];
  for (int i = 0; i < n; i++) format_user[i] = NULL;

  field2index = new int[n];
  argindex1 = new int[n];
  argindex2 = new int[n];

  // memory for computes, fixes, variables

  ncompute = 0;
  id_compute = new char*[n];
  compute_which = new int[n];
  computes = new Compute*[n];

  nfix = 0;
  id_fix = new char*[n];
  fixes = new Fix*[n];

  nvariable = 0;
  id_variable = new char*[n];
  variables = new int[n];
}

/* ----------------------------------------------------------------------
   deallocate all per-field memory
------------------------------------------------------------------------- */

void Stats::deallocate()
{
  int n = nfield;

  for (int i = 0; i < n; i++) delete [] keyword[i];
  delete [] keyword;
  delete [] vfunc;
  delete [] vtype;

  for (int i = 0; i < n; i++) delete [] format[i];
  delete [] format;
  for (int i = 0; i < n; i++) delete [] format_user[i];
  delete [] format_user;

  delete [] field2index;
  delete [] argindex1;
  delete [] argindex2;

  for (int i = 0; i < ncompute; i++) delete [] id_compute[i];
  delete [] id_compute;
  delete [] compute_which;
  delete [] computes;

  for (int i = 0; i < nfix; i++) delete [] id_fix[i];
  delete [] id_fix;
  delete [] fixes;

  for (int i = 0; i < nvariable; i++) delete [] id_variable[i];
  delete [] id_variable;
  delete [] variables;
}

/* ----------------------------------------------------------------------
   set fields of stats output from args
------------------------------------------------------------------------- */

void Stats::set_fields(int narg, char **arg)
{
  deallocate();
  nfield = narg;
  allocate();

  nfield = 0;

  // customize a new keyword by adding to if statement

  for (int i = 0; i < narg; i++) {
    if (strcmp(arg[i],"step") == 0) {
      addfield("Step",&Stats::compute_step,BIGINT);
    } else if (strcmp(arg[i],"elapsed") == 0) {
      addfield("Elapsed",&Stats::compute_elapsed,BIGINT);
    } else if (strcmp(arg[i],"elaplong") == 0) {
      addfield("Elapsed",&Stats::compute_elaplong,BIGINT);
    } else if (strcmp(arg[i],"dt") == 0) {
      addfield("Dt",&Stats::compute_dt,FLOAT);
    } else if (strcmp(arg[i],"cpu") == 0) {
      addfield("CPU",&Stats::compute_cpu,FLOAT);
    } else if (strcmp(arg[i],"tpcpu") == 0) {
      addfield("T/CPU",&Stats::compute_tpcpu,FLOAT);
    } else if (strcmp(arg[i],"spcpu") == 0) {
      addfield("S/CPU",&Stats::compute_spcpu,FLOAT);

    } else if (strcmp(arg[i],"np") == 0) {
      addfield("Np",&Stats::compute_np,BIGINT);
    } else if (strcmp(arg[i],"ntouch") == 0) {
      addfield("Ntouch",&Stats::compute_ntouch,BIGINT);
    } else if (strcmp(arg[i],"ncomm") == 0) {
      addfield("Ncomm",&Stats::compute_ncomm,BIGINT);
    } else if (strcmp(arg[i],"nbound") == 0) {
      addfield("Nbound",&Stats::compute_nbound,BIGINT);
    } else if (strcmp(arg[i],"nexit") == 0) {
      addfield("Nexit",&Stats::compute_nexit,BIGINT);
    } else if (strcmp(arg[i],"nscoll") == 0) {
      addfield("Nscoll",&Stats::compute_nscoll,BIGINT);
    } else if (strcmp(arg[i],"nscheck") == 0) {
      addfield("Nscheck",&Stats::compute_nscheck,BIGINT);
    } else if (strcmp(arg[i],"ncoll") == 0) {
      addfield("Ncoll",&Stats::compute_ncoll,BIGINT);
    } else if (strcmp(arg[i],"nattempt") == 0) {
      addfield("Natt",&Stats::compute_nattempt,BIGINT);
    } else if (strcmp(arg[i],"nreact") == 0) {
      addfield("Nreact",&Stats::compute_nreact,BIGINT);

    } else if (strcmp(arg[i],"npave") == 0) {
      addfield("Npave",&Stats::compute_npave,FLOAT);
    } else if (strcmp(arg[i],"ntouchave") == 0) {
      addfield("Ntouchave",&Stats::compute_ntouchave,FLOAT);
    } else if (strcmp(arg[i],"ncommave") == 0) {
      addfield("Ncommave",&Stats::compute_ncommave,FLOAT);
    } else if (strcmp(arg[i],"nboundave") == 0) {
      addfield("Nboundave",&Stats::compute_nboundave,FLOAT);
    } else if (strcmp(arg[i],"nexitave") == 0) {
      addfield("Nexitave",&Stats::compute_nexitave,FLOAT);
    } else if (strcmp(arg[i],"nscollave") == 0) {
      addfield("Nscollave",&Stats::compute_nscollave,FLOAT);
    } else if (strcmp(arg[i],"nscheckave") == 0) {
      addfield("Nschckave",&Stats::compute_nscheckave,FLOAT);
    } else if (strcmp(arg[i],"ncollave") == 0) {
      addfield("Ncollave",&Stats::compute_ncollave,FLOAT);
    } else if (strcmp(arg[i],"nattemptave") == 0) {
      addfield("Nattave",&Stats::compute_nattemptave,FLOAT);
    } else if (strcmp(arg[i],"nreactave") == 0) {
      addfield("Nattave",&Stats::compute_nreactave,FLOAT);

    } else if (strcmp(arg[i],"vol") == 0) {
      addfield("Volume",&Stats::compute_vol,FLOAT);
    } else if (strcmp(arg[i],"lx") == 0) {
      addfield("Lx",&Stats::compute_lx,FLOAT);
    } else if (strcmp(arg[i],"ly") == 0) {
      addfield("Ly",&Stats::compute_ly,FLOAT);
    } else if (strcmp(arg[i],"lz") == 0) {
      addfield("Lz",&Stats::compute_lz,FLOAT);

    } else if (strcmp(arg[i],"xlo") == 0) {
      addfield("Xlo",&Stats::compute_xlo,FLOAT);
    } else if (strcmp(arg[i],"xhi") == 0) {
      addfield("Xhi",&Stats::compute_xhi,FLOAT);
    } else if (strcmp(arg[i],"ylo") == 0) {
      addfield("Ylo",&Stats::compute_ylo,FLOAT);
    } else if (strcmp(arg[i],"yhi") == 0) {
      addfield("Yhi",&Stats::compute_yhi,FLOAT);
    } else if (strcmp(arg[i],"zlo") == 0) {
      addfield("Zlo",&Stats::compute_zlo,FLOAT);
    } else if (strcmp(arg[i],"zhi") == 0) {
      addfield("Zhi",&Stats::compute_zhi,FLOAT);

    // compute value = c_ID, fix value = f_ID, variable value = v_ID
    // count trailing [] and store int arguments
    // copy = at most 8 chars of ID to pass to addfield

    } else if ((strncmp(arg[i],"c_",2) == 0) || 
	       (strncmp(arg[i],"f_",2) == 0) ||
	       (strncmp(arg[i],"v_",2) == 0)) {

      int n = strlen(arg[i]);
      char *id = new char[n];
      strcpy(id,&arg[i][2]);
      char copy[9];
      strncpy(copy,id,8);
      copy[8] = '\0';

      // parse zero or one or two trailing brackets from ID
      // argindex1,argindex2 = int inside each bracket pair, 0 if no bracket

      char *ptr = strchr(id,'[');
      if (ptr == NULL) argindex1[nfield] = 0;
      else {
	*ptr = '\0';
	argindex1[nfield] = input->variable->int_between_brackets(ptr,0);
	ptr++;
	if (*ptr == '[') {
	  argindex2[nfield] = input->variable->int_between_brackets(ptr,0);
	  ptr++;
	} else argindex2[nfield] = 0;
      }

      if (arg[i][0] == 'c') {
	n = modify->find_compute(id);
	if (n < 0) error->all(FLERR,"Could not find stats compute ID");
	if (argindex1[nfield] == 0 && modify->compute[n]->scalar_flag == 0)
	  error->all(FLERR,"Stats compute does not compute scalar");
	if (argindex1[nfield] > 0 && argindex2[nfield] == 0) {
	  if (modify->compute[n]->vector_flag == 0)
	    error->all(FLERR,"Stats compute does not compute vector");
	  if (argindex1[nfield] > modify->compute[n]->size_vector)
	    error->all(FLERR,"Stats compute vector is accessed out-of-range");
	}
	if (argindex1[nfield] > 0 && argindex2[nfield] > 0) {
	  if (modify->compute[n]->array_flag == 0)
	    error->all(FLERR,"Stats compute does not compute array");
	  if (argindex1[nfield] > modify->compute[n]->size_array_rows ||
	      argindex2[nfield] > modify->compute[n]->size_array_cols)
	    error->all(FLERR,"Stats compute array is accessed out-of-range");
	}

	if (argindex1[nfield] == 0)
	  field2index[nfield] = add_compute(id,SCALAR);
	else if (argindex2[nfield] == 0)
	  field2index[nfield] = add_compute(id,VECTOR);
	else 
	  field2index[nfield] = add_compute(id,ARRAY);
	addfield(copy,&Stats::compute_compute,FLOAT);

      } else if (arg[i][0] == 'f') {
	n = modify->find_fix(id);
	if (n < 0) error->all(FLERR,"Could not find stats fix ID");
	if (argindex1[nfield] == 0 && modify->fix[n]->scalar_flag == 0)
	  error->all(FLERR,"Stats fix does not compute scalar");
	if (argindex1[nfield] > 0 && argindex2[nfield] == 0) {
	  if (modify->fix[n]->vector_flag == 0)
	    error->all(FLERR,"Stats fix does not compute vector");
	  if (argindex1[nfield] > modify->fix[n]->size_vector)
	    error->all(FLERR,"Stats fix vector is accessed out-of-range");
	}
	if (argindex1[nfield] > 0 && argindex2[nfield] > 0) {
	  if (modify->fix[n]->array_flag == 0)
	    error->all(FLERR,"Stats fix does not compute array");
	  if (argindex1[nfield] > modify->fix[n]->size_array_rows ||
	      argindex2[nfield] > modify->fix[n]->size_array_cols)
	    error->all(FLERR,"Stats fix array is accessed out-of-range");
	}

	field2index[nfield] = add_fix(id);
	addfield(copy,&Stats::compute_fix,FLOAT);

      } else if (arg[i][0] == 'v') {
	n = input->variable->find(id);
	if (n < 0) error->all(FLERR,"Could not find stats variable name");
	if (input->variable->equal_style(n) == 0)
	  error->all(FLERR,"Stats variable is not equal-style variable");
	if (argindex1[nfield]) 
	  error->all(FLERR,"Stats variable cannot be indexed");

	field2index[nfield] = add_variable(id);
	addfield(copy,&Stats::compute_variable,FLOAT);
      }

      delete [] id;

    } else error->all(FLERR,"Invalid keyword in stats_style command");
  }
}

/* ----------------------------------------------------------------------
   add field to list of quantities to print
------------------------------------------------------------------------- */

void Stats::addfield(const char *key, FnPtr func, int typeflag)
{
  strcpy(keyword[nfield],key);
  vfunc[nfield] = func;
  vtype[nfield] = typeflag;
  nfield++;
}

/* ----------------------------------------------------------------------
   add compute ID to list of Compute objects to call
   return location of where this Compute is in list
   if already in list with same which, do not add, just return index
------------------------------------------------------------------------- */

int Stats::add_compute(const char *id, int which)
{
  int icompute;
  for (icompute = 0; icompute < ncompute; icompute++)
    if ((strcmp(id,id_compute[icompute]) == 0) && 
	which == compute_which[icompute]) break;
  if (icompute < ncompute) return icompute;

  int n = strlen(id) + 1;
  id_compute[ncompute] = new char[n];
  strcpy(id_compute[ncompute],id);
  compute_which[ncompute] = which;
  ncompute++;
  return ncompute-1;
}

/* ----------------------------------------------------------------------
   add fix ID to list of Fix objects to call
------------------------------------------------------------------------- */

int Stats::add_fix(const char *id)
{
  int n = strlen(id) + 1;
  id_fix[nfix] = new char[n];
  strcpy(id_fix[nfix],id);
  nfix++;
  return nfix-1;
}

/* ----------------------------------------------------------------------
   add variable ID to list of Variables to evaluate
------------------------------------------------------------------------- */

int Stats::add_variable(const char *id)
{
  int n = strlen(id) + 1;
  id_variable[nvariable] = new char[n];
  strcpy(id_variable[nvariable],id);
  nvariable++;
  return nvariable-1;
}

/* ----------------------------------------------------------------------
   compute a single stats value, word is any stats keyword
   called when a variable is evaluated by Variable class
   return value as double in answer
   return 0 if str is recoginzed keyword, 1 if unrecognized
   customize a new keyword by adding to if statement
------------------------------------------------------------------------- */

int Stats::evaluate_keyword(char *word, double *answer)
{
  // invoke a lo-level stats routine to compute the variable value

  if (strcmp(word,"step") == 0) {
    compute_step();
    dvalue = bivalue;

  } else if (strcmp(word,"elapsed") == 0) {
    if (update->runflag == 0) 
      error->all(FLERR,
		 "Variable stats keyword cannot be used between runs");
    compute_elapsed();
    dvalue = bivalue;

  } else if (strcmp(word,"elaplong") == 0) {
    if (update->runflag == 0) 
      error->all(FLERR,
		 "Variable stats keyword cannot be used between runs");
    compute_elaplong();
    dvalue = bivalue;

  } else if (strcmp(word,"dt") == 0) {
    compute_dt();

  } else if (strcmp(word,"cpu") == 0) {
    if (update->runflag == 0) 
      error->all(FLERR,
		 "Variable stats keyword cannot be used between runs");
    compute_cpu();

  } else if (strcmp(word,"tpcpu") == 0) {
    if (update->runflag == 0) 
      error->all(FLERR,
		 "Variable stats keyword cannot be used between runs");
    compute_tpcpu();

  } else if (strcmp(word,"spcpu") == 0) {
    if (update->runflag == 0) 
      error->all(FLERR,
		 "Variable stats keyword cannot be used between runs");
    compute_spcpu();

  } else if (strcmp(word,"np") == 0) {
    compute_np();
    dvalue = bivalue;
  } else if (strcmp(word,"ntouch") == 0) {
    compute_ntouch();
    dvalue = bivalue;
  } else if (strcmp(word,"ncomm") == 0) {
    compute_ncomm();
    dvalue = bivalue;
  } else if (strcmp(word,"nbound") == 0) {
    compute_nbound();
    dvalue = bivalue;
  } else if (strcmp(word,"nexit") == 0) {
    compute_nexit();
    dvalue = bivalue;
  } else if (strcmp(word,"nscoll") == 0) {
    compute_nscoll();
    dvalue = bivalue;
  } else if (strcmp(word,"nscheck") == 0) {
    compute_nscheck();
    dvalue = bivalue;
  } else if (strcmp(word,"ncoll") == 0) {
    compute_ncoll();
    dvalue = bivalue;
  } else if (strcmp(word,"nattempt") == 0) {
    compute_nattempt();
    dvalue = bivalue;
  } else if (strcmp(word,"nreact") == 0) {
    compute_nreact();
    dvalue = bivalue;
  }

  else if (strcmp(word,"npave") == 0) compute_npave();
  else if (strcmp(word,"ntouchave") == 0) compute_ntouchave();
  else if (strcmp(word,"ncommave") == 0) compute_ncommave();
  else if (strcmp(word,"nboundave") == 0) compute_nboundave();
  else if (strcmp(word,"nexitave") == 0) compute_nexitave();
  else if (strcmp(word,"nscollave") == 0) compute_nscollave();
  else if (strcmp(word,"nscheckave") == 0) compute_nscheckave();
  else if (strcmp(word,"ncollave") == 0) compute_ncollave();
  else if (strcmp(word,"nattemptave") == 0) compute_nattemptave();
  else if (strcmp(word,"nreactave") == 0) compute_nreactave();

  else if (strcmp(word,"vol") == 0) compute_vol();
  else if (strcmp(word,"lx") == 0) compute_lx();
  else if (strcmp(word,"ly") == 0) compute_ly();
  else if (strcmp(word,"lz") == 0) compute_lz();

  else if (strcmp(word,"xlo") == 0) compute_xlo();
  else if (strcmp(word,"xhi") == 0) compute_xhi();
  else if (strcmp(word,"ylo") == 0) compute_ylo();
  else if (strcmp(word,"yhi") == 0) compute_yhi();
  else if (strcmp(word,"zlo") == 0) compute_zlo();
  else if (strcmp(word,"zhi") == 0) compute_zhi();

  else return 1;

  *answer = dvalue;
  return 0;
}

/* ----------------------------------------------------------------------
   extraction of Compute, Fix, Variable results
------------------------------------------------------------------------- */

void Stats::compute_compute()
{
  int m = field2index[ifield];
  Compute *compute = computes[m];

  if (compute_which[m] == SCALAR)
    dvalue = compute->scalar;
  else if (compute_which[m] == VECTOR)
    dvalue = compute->vector[argindex1[ifield]-1];
  else
    dvalue = compute->array[argindex1[ifield]-1][argindex2[ifield]-1];
}

/* ---------------------------------------------------------------------- */

void Stats::compute_fix()
{
  int m = field2index[ifield];
  Fix *fix = fixes[m];

  if (argindex1[ifield] == 0)
    dvalue = fix->compute_scalar();
  else if (argindex2[ifield] == 0)
    dvalue = fix->compute_vector(argindex1[ifield]-1);
  else
    dvalue = fix->compute_array(argindex1[ifield]-1,argindex2[ifield]-1);
}

/* ---------------------------------------------------------------------- */

void Stats::compute_variable()
{
  dvalue = input->variable->compute_equal(variables[field2index[ifield]]);
}

/* ----------------------------------------------------------------------
   one method for every keyword stats can output
   called by compute() or evaluate_keyword()
   compute will have already been called
   set ivalue/dvalue/bivalue if value is int/double/bigint
   customize a new keyword by adding a method
------------------------------------------------------------------------- */

void Stats::compute_step()
{
  bivalue = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void Stats::compute_elapsed()
{
  bivalue = update->ntimestep - update->firststep;
}

/* ---------------------------------------------------------------------- */

void Stats::compute_elaplong()
{
  bivalue = update->ntimestep - update->beginstep;
}

/* ---------------------------------------------------------------------- */

void Stats::compute_dt()
{
  dvalue = update->dt;
}

/* ---------------------------------------------------------------------- */

void Stats::compute_cpu()
{
  if (firststep == 0) dvalue = 0.0;
  else dvalue = timer->elapsed(TIME_LOOP);
}

/* ---------------------------------------------------------------------- */

void Stats::compute_tpcpu()
{
  double new_cpu;
  double new_time = update->ntimestep * update->dt;

  if (firststep == 0) {
    new_cpu = 0.0;
    dvalue = 0.0;
  } else {
    new_cpu = timer->elapsed(TIME_LOOP);
    double cpu_diff = new_cpu - last_tpcpu;
    double time_diff = new_time - last_time;
    if (time_diff > 0.0 && cpu_diff > 0.0) dvalue = time_diff/cpu_diff;
    else dvalue = 0.0;
  }

  last_time = new_time;
  last_tpcpu = new_cpu;
}

/* ---------------------------------------------------------------------- */

void Stats::compute_spcpu()
{
  double new_cpu;
  int new_step = update->ntimestep;

  if (firststep == 0) {
    new_cpu = 0.0;
    dvalue = 0.0;
  } else {
    new_cpu = timer->elapsed(TIME_LOOP);
    double cpu_diff = new_cpu - last_spcpu;
    int step_diff = new_step - last_step;
    if (cpu_diff > 0.0) dvalue = step_diff/cpu_diff;
    else dvalue = 0.0;
  }

  last_step = new_step;
  last_spcpu = new_cpu;
}

/* ---------------------------------------------------------------------- */

void Stats::compute_np()
{
  bigint n = particle->nlocal;
  MPI_Allreduce(&n,&particle->nglobal,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  bivalue = particle->nglobal;
}

/* ---------------------------------------------------------------------- */

void Stats::compute_ntouch()
{
  bigint n = update->ntouch_one;
  MPI_Allreduce(&n,&bivalue,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void Stats::compute_ncomm()
{
  bigint n = update->ncomm_one;
  MPI_Allreduce(&n,&bivalue,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void Stats::compute_nbound()
{
  bigint n = update->nboundary_one;
  MPI_Allreduce(&n,&bivalue,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void Stats::compute_nexit()
{
  bigint n = update->nexit_one;
  MPI_Allreduce(&n,&bivalue,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void Stats::compute_nscoll()
{
  bigint n = update->nscollide_one;
  MPI_Allreduce(&n,&bivalue,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void Stats::compute_nscheck()
{
  bigint n = update->nscheck_one;
  MPI_Allreduce(&n,&bivalue,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void Stats::compute_ncoll()
{
  if (!collide) bivalue = 0;
  else {
    bigint n = collide->ncollide_one;
    MPI_Allreduce(&n,&bivalue,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  }
}

/* ---------------------------------------------------------------------- */

void Stats::compute_nattempt()
{
  if (!collide) bivalue = 0;
  else {
    bigint n = collide->nattempt_one;
    MPI_Allreduce(&n,&bivalue,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  }
}

/* ---------------------------------------------------------------------- */

void Stats::compute_nreact()
{
  if (!collide) bivalue = 0;
  else {
    bigint n = collide->nreact_one;
    MPI_Allreduce(&n,&bivalue,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  }
}

/* ---------------------------------------------------------------------- */

void Stats::compute_npave()
{
  MPI_Allreduce(&update->nmove_running,&bivalue,1,MPI_SPARTA_BIGINT,
		MPI_SUM,world);
  if (update->ntimestep == update->firststep) dvalue = 0.0;
  else dvalue = 1.0*bivalue / (update->ntimestep - update->firststep);
}

/* ---------------------------------------------------------------------- */

void Stats::compute_ntouchave()
{
  MPI_Allreduce(&update->ntouch_running,&bivalue,1,MPI_SPARTA_BIGINT,
		MPI_SUM,world);
  if (update->ntimestep == update->firststep) dvalue = 0.0;
  else dvalue = 1.0*bivalue / (update->ntimestep - update->firststep);
}

/* ---------------------------------------------------------------------- */

void Stats::compute_ncommave()
{
  MPI_Allreduce(&update->ncomm_running,&bivalue,1,MPI_SPARTA_BIGINT,
		MPI_SUM,world);
  if (update->ntimestep == update->firststep) dvalue = 0.0;
  else dvalue = 1.0*bivalue / (update->ntimestep - update->firststep);
}

/* ---------------------------------------------------------------------- */

void Stats::compute_nboundave()
{
  MPI_Allreduce(&update->nboundary_running,&bivalue,1,MPI_SPARTA_BIGINT,
		MPI_SUM,world);
  if (update->ntimestep == update->firststep) dvalue = 0.0;
  else dvalue = 1.0*bivalue / (update->ntimestep - update->firststep);
}

/* ---------------------------------------------------------------------- */

void Stats::compute_nexitave()
{
  MPI_Allreduce(&update->nexit_running,&bivalue,1,MPI_SPARTA_BIGINT,
		MPI_SUM,world);
  if (update->ntimestep == update->firststep) dvalue = 0.0;
  else dvalue = 1.0*bivalue / (update->ntimestep - update->firststep);
}

/* ---------------------------------------------------------------------- */

void Stats::compute_nscollave()
{
  MPI_Allreduce(&update->nscollide_running,&bivalue,1,MPI_SPARTA_BIGINT,
		MPI_SUM,world);
  if (update->ntimestep == update->firststep) dvalue = 0.0;
  else dvalue = 1.0*bivalue / (update->ntimestep - update->firststep);
}

/* ---------------------------------------------------------------------- */

void Stats::compute_nscheckave()
{
  MPI_Allreduce(&update->nscheck_running,&bivalue,1,MPI_SPARTA_BIGINT,
		MPI_SUM,world);
  if (update->ntimestep == update->firststep) dvalue = 0.0;
  else dvalue = 1.0*bivalue / (update->ntimestep - update->firststep);
}

/* ---------------------------------------------------------------------- */

void Stats::compute_ncollave()
{
  if (!collide) dvalue = 0.0;
  else {
    MPI_Allreduce(&collide->ncollide_running,&bivalue,1,MPI_SPARTA_BIGINT,
		  MPI_SUM,world);
    if (update->ntimestep == update->firststep) dvalue = 0.0;
    else dvalue = 1.0*bivalue / (update->ntimestep - update->firststep);
  }
}

/* ---------------------------------------------------------------------- */

void Stats::compute_nattemptave()
{
  if (!collide) dvalue = 0.0;
  else {
    MPI_Allreduce(&collide->nattempt_running,&bivalue,1,MPI_SPARTA_BIGINT,
		  MPI_SUM,world);
    if (update->ntimestep == update->firststep) dvalue = 0.0;
    else dvalue = 1.0*bivalue / (update->ntimestep - update->firststep);
  }
}

/* ---------------------------------------------------------------------- */

void Stats::compute_nreactave()
{
  if (!collide) dvalue = 0.0;
  else {
    MPI_Allreduce(&collide->nreact_running,&bivalue,1,MPI_SPARTA_BIGINT,
		  MPI_SUM,world);
    if (update->ntimestep == update->firststep) dvalue = 0.0;
    else dvalue = 1.0*bivalue / (update->ntimestep - update->firststep);
  }
}

/* ---------------------------------------------------------------------- */

void Stats::compute_vol()
{
  if (domain->dimension == 3)
    dvalue = domain->xprd * domain->yprd * domain->zprd;
  else
    dvalue = domain->xprd * domain->yprd;
}

/* ---------------------------------------------------------------------- */

void Stats::compute_lx()
{
  dvalue = domain->xprd;
}

/* ---------------------------------------------------------------------- */

void Stats::compute_ly()
{
  dvalue = domain->yprd;
}

/* ---------------------------------------------------------------------- */

void Stats::compute_lz()
{
  dvalue = domain->zprd;
}

/* ---------------------------------------------------------------------- */

void Stats::compute_xlo()
{
  dvalue = domain->boxlo[0];
}

/* ---------------------------------------------------------------------- */

void Stats::compute_xhi()
{
  dvalue = domain->boxhi[0];
}

/* ---------------------------------------------------------------------- */

void Stats::compute_ylo()
{
  dvalue = domain->boxlo[1];
}

/* ---------------------------------------------------------------------- */

void Stats::compute_yhi()
{
  dvalue = domain->boxhi[1];
}

/* ---------------------------------------------------------------------- */

void Stats::compute_zlo()
{
  dvalue = domain->boxlo[2];
}

/* ---------------------------------------------------------------------- */

void Stats::compute_zhi()
{
  dvalue = domain->boxhi[2];
}

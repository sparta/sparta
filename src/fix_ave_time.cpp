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
#include "stdlib.h"
#include "string.h"
#include "fix_ave_time.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{COMPUTE,FIX,VARIABLE};
enum{ONE,RUNNING,WINDOW};
enum{SCALAR,VECTOR};

#define INVOKED_SCALAR 1
#define INVOKED_VECTOR 2
#define INVOKED_ARRAY 4
#define DELTA 8;

/* ---------------------------------------------------------------------- */

FixAveTime::FixAveTime(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg < 6) error->all(FLERR,"Illegal fix ave/time command");

  MPI_Comm_rank(world,&me);

  nevery = atoi(arg[2]);
  nrepeat = atoi(arg[3]);
  nfreq = atoi(arg[4]);

  time_depend = 1;
  global_freq = nfreq;

  // scan values to count them
  // then read options so know mode = SCALAR/VECTOR before re-reading values

  nvalues = 0;

  int iarg = 5;
  while (iarg < narg) {
    if ((strncmp(arg[iarg],"c_",2) == 0) || 
	(strncmp(arg[iarg],"f_",2) == 0) || 
	(strncmp(arg[iarg],"v_",2) == 0)) {
      nvalues++;
      iarg++;
    } else break;
  }

  if (nvalues == 0) error->all(FLERR,"No values in fix ave/time command");

  options(iarg,narg,arg);

  // expand args if any have wildcard character "*"
  // this can reset nvalues

  int expand = 0;
  char **earg;
  nvalues = input->expand_args(nvalues,&arg[5],mode,earg);

  if (earg != &arg[5]) expand = 1;
  arg = earg;

  // parse values

  which = new int[nvalues];
  argindex = new int[nvalues];
  value2index = new int[nvalues];
  offcol = new int[nvalues];
  ids = new char*[nvalues];

  for (int i = 0; i < nvalues; i++) {
    if (arg[i][0] == 'c') which[i] = COMPUTE;
    else if (arg[i][0] == 'f') which[i] = FIX;
    else if (arg[i][0] == 'v') which[i] = VARIABLE;

    int n = strlen(arg[i]);
    char *suffix = new char[n];
    strcpy(suffix,&arg[i][2]);

    char *ptr = strchr(suffix,'[');
    if (ptr) {
      if (suffix[strlen(suffix)-1] != ']')
        error->all(FLERR,"Illegal fix ave/time command");
      argindex[i] = atoi(ptr+1);
      *ptr = '\0';
    } else argindex[i] = 0;

    n = strlen(suffix) + 1;
    ids[i] = new char[n];
    strcpy(ids[i],suffix);
    delete [] suffix;
  }

  // set off columns now that nvalues is finalized

  for (int i = 0; i < nvalues; i++) offcol[i] = 0;
  for (int i = 0; i < noff; i++) {
    if (offlist[i] < 1 || offlist[i] > nvalues)
      error->all(FLERR,"Invalid fix ave/time off column");
    offcol[offlist[i]-1] = 1;
  }

  // setup and error check
  // for fix inputs, check that fix frequency is acceptable

  if (nevery <= 0 || nrepeat <= 0 || nfreq <= 0)
    error->all(FLERR,"Illegal fix ave/time command");
  if (nfreq % nevery || (nrepeat-1)*nevery >= nfreq)
    error->all(FLERR,"Illegal fix ave/time command");

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE && mode == SCALAR) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
	error->all(FLERR,"Compute ID for fix ave/time does not exist");
      if (argindex[i] == 0 && modify->compute[icompute]->scalar_flag == 0)
	error->all(FLERR,"Fix ave/time compute does not calculate a scalar");
      if (argindex[i] && modify->compute[icompute]->vector_flag == 0)
	error->all(FLERR,"Fix ave/time compute does not calculate a vector");
      if (argindex[i] && argindex[i] > modify->compute[icompute]->size_vector)
	error->all(FLERR,
		   "Fix ave/time compute vector is accessed out-of-range");

    } else if (which[i] == COMPUTE && mode == VECTOR) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
	error->all(FLERR,"Compute ID for fix ave/time does not exist");
      if (argindex[i] == 0 && modify->compute[icompute]->vector_flag == 0)
	error->all(FLERR,"Fix ave/time compute does not calculate a vector");
      if (argindex[i] && modify->compute[icompute]->array_flag == 0)
	error->all(FLERR,"Fix ave/time compute does not calculate an array");
      if (argindex[i] && 
	  argindex[i] > modify->compute[icompute]->size_array_cols)
	error->all(FLERR,"Fix ave/time compute array is accessed out-of-range");

    } else if (which[i] == FIX && mode == SCALAR) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
	error->all(FLERR,"Fix ID for fix ave/time does not exist");
      if (argindex[i] == 0 && modify->fix[ifix]->scalar_flag == 0)
	error->all(FLERR,"Fix ave/time fix does not calculate a scalar");
      if (argindex[i] && modify->fix[ifix]->vector_flag == 0)
	error->all(FLERR,"Fix ave/time fix does not calculate a vector");
      if (argindex[i] && argindex[i] > modify->fix[ifix]->size_vector)
	error->all(FLERR,"Fix ave/time fix vector is accessed out-of-range");
      if (nevery % modify->fix[ifix]->global_freq)
	error->all(FLERR,
		   "Fix for fix ave/time not computed at compatible time");

    } else if (which[i] == FIX && mode == VECTOR) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
	error->all(FLERR,"Fix ID for fix ave/time does not exist");
      if (argindex[i] == 0 && modify->fix[ifix]->vector_flag == 0)
	error->all(FLERR,"Fix ave/time fix does not calculate a vector");
      if (argindex[i] && modify->fix[ifix]->array_flag == 0)
	error->all(FLERR,"Fix ave/time fix does not calculate an array");
      if (argindex[i] && argindex[i] > modify->fix[ifix]->size_array_cols)
	error->all(FLERR,"Fix ave/time fix array is accessed out-of-range");
      if (nevery % modify->fix[ifix]->global_freq)
	error->all(FLERR,
		   "Fix for fix ave/time not computed at compatible time");

    } else if (which[i] == VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
	error->all(FLERR,"Variable name for fix ave/time does not exist");
      if (input->variable->equal_style(ivariable) == 0)
	error->all(FLERR,"Fix ave/time variable is not equal-style variable");
      if (mode == VECTOR)
	error->all(FLERR,"Fix ave/time cannot use variable with vector mode");
    }
  }

  // if VECTOR mode, check that all columns are same length
  // nrows = # of rows in output array

  if (mode == VECTOR) {
    int length;

    for (int i = 0; i < nvalues; i++) {
      if (which[i] == COMPUTE) {
	int icompute = modify->find_compute(ids[i]);
	if (argindex[i] == 0) length = modify->compute[icompute]->size_vector;
	else length = modify->compute[icompute]->size_array_rows;
      } else if (which[i] == FIX) {
	int ifix = modify->find_fix(ids[i]);
	if (argindex[i] == 0) length = modify->fix[ifix]->size_vector;
	else length = modify->fix[ifix]->size_array_rows;
      }
      if (i == 0) nrows = length;
      else if (length != nrows)
	error->all(FLERR,"Fix ave/time columns are inconsistent lengths");
    }

    column = new double[nrows];
  } else column = NULL;

  // print file comment lines
  // for mode = VECTOR, cannot use arg to print
  // since array args may have been expanded to multiple vectors

  if (fp && me == 0) {
    if (title1) fprintf(fp,"%s\n",title1);
    else fprintf(fp,"# Time-averaged data for fix %s\n",id);
    if (title2) fprintf(fp,"%s\n",title2);
    else if (mode == SCALAR) {
      fprintf(fp,"# TimeStep");
      for (int i = 0; i < nvalues; i++) fprintf(fp," %s",earg[i]);
      fprintf(fp,"\n");
    } else fprintf(fp,"# TimeStep Number-of-rows\n");
    if (title3 && mode == VECTOR) fprintf(fp,"%s\n",title3);
    else if (mode == VECTOR) {
      fprintf(fp,"# Row");
      for (int i = 0; i < nvalues; i++) fprintf(fp," %s",earg[i]);
      fprintf(fp,"\n");
    }
  }

  delete [] title1;
  delete [] title2;
  delete [] title3;

  // if wildcard expansion occurred, free earg memory from expand_args()
  // wait to do this until after file comment lines are printed

  if (expand) {
    for (int i = 0; i < nvalues; i++) delete [] earg[i];
    memory->sfree(earg);
  }

  // allocate accumulators

  vector = vector_total = NULL;
  vector_list = NULL;
  array = array_total = NULL;
  array_list = NULL;

  if (mode == SCALAR) {
    vector = new double[nvalues];
    vector_total = new double[nvalues];
    if (ave == WINDOW)
      memory->create(vector_list,nwindow,nvalues,"ave/time:vector_list");
  } else {
    memory->create(array,nrows,nvalues,"ave/time:array");
    memory->create(array_total,nrows,nvalues,"ave/time:array_total");
    if (ave == WINDOW)
      memory->create(array_list,nwindow,nrows,nvalues,"ave/time:array_list");
  }

  // this fix produces either a global scalar or vector or array
  // SCALAR mode produces either a scalar or vector
  // VECTOR mode produces either a vector or array
  
  if (mode == SCALAR) {
    if (nvalues == 1) scalar_flag = 1;
    else {
      vector_flag = 1;
      size_vector = nvalues;
    }
  } else {
    if (nvalues == 1) {
      vector_flag = 1;
      size_vector = nrows;
    } else {
      array_flag = 1;
      size_array_rows = nrows;
      size_array_cols = nvalues;
    }
  }

  // initializations
  // set vector_total/array_total to zero since it accumulates

  irepeat = 0;
  iwindow = window_limit = 0;
  norm = 0;

  if (mode == SCALAR)
    for (int i = 0; i < nvalues; i++) vector_total[i] = 0.0;
  else
    for (int i = 0; i < nrows; i++)
      for (int j = 0; j < nvalues; j++)
        array_total[i][j] = 0.0;

  // nvalid = next step on which end_of_step does something
  // add nvalid to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  nvalid = nextvalid();
  modify->addstep_compute_all(nvalid);
}

/* ---------------------------------------------------------------------- */

FixAveTime::~FixAveTime()
{
  memory->destroy(which);
  memory->destroy(argindex);
  memory->destroy(value2index);
  memory->destroy(offcol);
  for (int i = 0; i < nvalues; i++) delete [] ids[i];
  memory->sfree(ids);

  if (fp && me == 0) fclose(fp);

  delete [] vector;
  delete [] vector_total;
  delete [] column;
  memory->destroy(array);
  memory->destroy(array_total);
  memory->destroy(array_list);
}

/* ---------------------------------------------------------------------- */

int FixAveTime::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAveTime::init()
{
  // set current indices for all computes,fixes,variables

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
	error->all(FLERR,"Compute ID for fix ave/time does not exist");
      value2index[i] = icompute;

    } else if (which[i] == FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0) 
	error->all(FLERR,"Fix ID for fix ave/time does not exist");
      value2index[i] = ifix;

    } else if (which[i] == VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0) 
	error->all(FLERR,"Variable name for fix ave/time does not exist");
      value2index[i] = ivariable;
    }
  }
}

/* ----------------------------------------------------------------------
   only does something if nvalid = current timestep
------------------------------------------------------------------------- */

void FixAveTime::setup()
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveTime::end_of_step()
{
  // skip if not step which requires doing something

  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;

  if (mode == SCALAR) invoke_scalar(ntimestep);
  else invoke_vector(ntimestep);
}

/* ---------------------------------------------------------------------- */

void FixAveTime::invoke_scalar(bigint ntimestep)
{
  int i,m;
  double scalar;

  // zero if first step

  if (irepeat == 0)
    for (i = 0; i < nvalues; i++) vector[i] = 0.0;

  // accumulate results of computes,fixes,variables to local copy
  // compute/fix/variable may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  for (i = 0; i < nvalues; i++) {
    m = value2index[i];
    
    // invoke compute if not previously invoked
    
    if (which[i] == COMPUTE) {
      Compute *compute = modify->compute[m];

      if (argindex[i] == 0) {
	if (!(compute->invoked_flag & INVOKED_SCALAR)) {
	  compute->compute_scalar();
	  compute->invoked_flag |= INVOKED_SCALAR;
	}
	scalar = compute->scalar;
      } else {
	if (!(compute->invoked_flag & INVOKED_VECTOR)) {
	  compute->compute_vector();
	  compute->invoked_flag |= INVOKED_VECTOR;
	}
	scalar = compute->vector[argindex[i]-1];
      }
      
    // access fix fields, guaranteed to be ready
      
    } else if (which[i] == FIX) {
      if (argindex[i] == 0) 
	scalar = modify->fix[m]->compute_scalar();
      else
	scalar = modify->fix[m]->compute_vector(argindex[i]-1);
      
    // evaluate equal-style variable
      
    } else if (which[i] == VARIABLE)
      scalar = input->variable->compute_equal(m);
    
    // add value to vector or just set directly if offcol is set
    
    if (offcol[i]) vector[i] = scalar;
    else vector[i] += scalar;
  }

  // done if irepeat < nrepeat
  // else reset irepeat and nvalid

  irepeat++;
  if (irepeat < nrepeat) {
    nvalid += nevery;
    modify->addstep_compute(nvalid);
    return;
  }

  irepeat = 0;
  nvalid = ntimestep + nfreq - (nrepeat-1)*nevery;
  modify->addstep_compute(nvalid);

  // average the final result for the Nfreq timestep
  // no other normalization factor used

  double repeat = nrepeat;
  for (i = 0; i < nvalues; i++)
    if (offcol[i] == 0) vector[i] /= repeat;

  // if ave = ONE, only single Nfreq timestep value is needed
  // if ave = RUNNING, combine with all previous Nfreq timestep values
  // if ave = WINDOW, combine with nwindow most recent Nfreq timestep values

  if (ave == ONE) {
    for (i = 0; i < nvalues; i++) vector_total[i] = vector[i];
    norm = 1;
    
  } else if (ave == RUNNING) {
    for (i = 0; i < nvalues; i++) vector_total[i] += vector[i];
    norm++;
    
  } else if (ave == WINDOW) {
    for (i = 0; i < nvalues; i++) {
      vector_total[i] += vector[i];
      if (window_limit) vector_total[i] -= vector_list[iwindow][i];
      vector_list[iwindow][i] = vector[i];
    }
    
    iwindow++;
    if (iwindow == nwindow) {
      iwindow = 0;
      window_limit = 1;
    }
    if (window_limit) norm = nwindow;
    else norm = iwindow;
  }
  
  // insure any columns with offcol set are effectively set to last value

  for (i = 0; i < nvalues; i++)
    if (offcol[i]) vector_total[i] = norm*vector[i];

  // output result to file

  if (fp && me == 0) {
    fprintf(fp,BIGINT_FORMAT,ntimestep);
    for (i = 0; i < nvalues; i++) fprintf(fp," %g",vector_total[i]/norm);
    fprintf(fp,"\n");
    fflush(fp);
  }
}

/* ---------------------------------------------------------------------- */

void FixAveTime::invoke_vector(bigint ntimestep)
{
  int i,j,m;

  // zero if first step

  if (irepeat == 0)
    for (i = 0; i < nrows; i++)
      for (j = 0; j < nvalues; j++) array[i][j] = 0.0;
  
  // accumulate results of computes,fixes,variables to local copy
  // compute/fix/variable may invoke computes so wrap with clear/add
  
  modify->clearstep_compute();
  
  for (j = 0; j < nvalues; j++) {
    m = value2index[j];
    
    // invoke compute if not previously invoked
    
    if (which[j] == COMPUTE) {
      Compute *compute = modify->compute[m];
      
      if (argindex[j] == 0) {
	if (!(compute->invoked_flag & INVOKED_VECTOR)) {
	  compute->compute_vector();
	  compute->invoked_flag |= INVOKED_VECTOR;
	}
	double *cvector = compute->vector;
	for (i = 0; i < nrows; i++)
	  column[i] = cvector[i];
	
      } else {
	if (!(compute->invoked_flag & INVOKED_ARRAY)) {
	  compute->compute_array();
	  compute->invoked_flag |= INVOKED_ARRAY;
	}
	double **carray = compute->array;
	int icol = argindex[j]-1;
	for (i = 0; i < nrows; i++)
	  column[i] = carray[i][icol];
      }
      
    // access fix fields, guaranteed to be ready
      
    } else if (which[j] == FIX) {
      Fix *fix = modify->fix[m];
      if (argindex[j] == 0)
	for (i = 0; i < nrows; i++)
	  column[i] = fix->compute_vector(i);
      else {
	int icol = argindex[j]-1;
	for (i = 0; i < nrows; i++)
	  column[i] = fix->compute_array(i,icol);
      }
    }
    
    // add columns of values to array or just set directly if offcol is set
    
    if (offcol[j])
      for (i = 0; i < nrows; i++)
	array[i][j] = column[i];
    else
      for (i = 0; i < nrows; i++)
	array[i][j] += column[i];
  }
  
  // done if irepeat < nrepeat
  // else reset irepeat and nvalid

  irepeat++;
  if (irepeat < nrepeat) {
    nvalid += nevery;
    modify->addstep_compute(nvalid);
    return;
  }

  irepeat = 0;
  nvalid = ntimestep + nfreq - (nrepeat-1)*nevery;
  modify->addstep_compute(nvalid);

  // average the final result for the Nfreq timestep
  
  double repeat = nrepeat;

  for (m = 0; m < nvalues; m++) 
    if (offcol[m] == 0)
      for (i = 0; i < nrows; i++) array[i][m] /= repeat;
  
  // if ave = ONE, only single Nfreq timestep value is needed
  // if ave = RUNNING, combine with all previous Nfreq timestep values
  // if ave = WINDOW, combine with nwindow most recent Nfreq timestep values
  
  if (ave == ONE) {
    for (i = 0; i < nrows; i++)
      for (j = 0; j < nvalues; j++) array_total[i][j] = array[i][j];
    norm = 1;
    
  } else if (ave == RUNNING) {
    for (i = 0; i < nrows; i++)
      for (j = 0; j < nvalues; j++) array_total[i][j] += array[i][j];
    norm++;
    
  } else if (ave == WINDOW) {
    for (i = 0; i < nrows; i++)
      for (j = 0; j < nvalues; j++) {
	array_total[i][j] += array[i][j];
	if (window_limit) array_total[i][j] -= array_list[iwindow][i][j];
	array_list[iwindow][i][j] = array[i][j];
      }
    
    iwindow++;
    if (iwindow == nwindow) {
      iwindow = 0;
      window_limit = 1;
    }
    if (window_limit) norm = nwindow;
    else norm = iwindow;
  }
  
  // insure any columns with offcol set are effectively set to last value

  for (i = 0; i < nrows; i++)
    for (j = 0; j < nvalues; j++)
      if (offcol[j]) array_total[i][j] = norm*array[i][j];
  
  // output result to file

  if (fp && me == 0) {
    fprintf(fp,BIGINT_FORMAT " %d\n",ntimestep,nrows);
    for (i = 0; i < nrows; i++) {
      fprintf(fp,"%d",i+1);
      for (j = 0; j < nvalues; j++) fprintf(fp," %g",array_total[i][j]/norm);
      fprintf(fp,"\n");
    }
    fflush(fp);
  }
}

/* ----------------------------------------------------------------------
   return scalar value
------------------------------------------------------------------------- */

double FixAveTime::compute_scalar()
{
  if (norm) return vector_total[0]/norm;
  return 0.0;
}

/* ----------------------------------------------------------------------
   return Ith vector value
------------------------------------------------------------------------- */

double FixAveTime::compute_vector(int i)
{
  if (norm) {
    if (mode == SCALAR) return vector_total[i]/norm;
    if (mode == VECTOR) return array_total[i][0];
  }
  return 0.0;
}

/* ----------------------------------------------------------------------
   return I,J array value
------------------------------------------------------------------------- */

double FixAveTime::compute_array(int i, int j)
{
  if (norm) return array_total[i][j]/norm;
  return 0.0;
}

/* ----------------------------------------------------------------------
   parse optional args
------------------------------------------------------------------------- */

void FixAveTime::options(int iarg, int narg, char **arg)
{
  // option defaults

  fp = NULL;
  ave = ONE;
  startstep = 0;
  mode = SCALAR;
  noff = 0;
  offlist = NULL;
  title1 = NULL;
  title2 = NULL;
  title3 = NULL;

  // optional args

  while (iarg < narg) {
    if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/time command");
      if (me == 0) {
	fp = fopen(arg[iarg+1],"w");
	if (fp == NULL) {
	  char str[128];
	  sprintf(str,"Cannot open fix ave/time file %s",arg[iarg+1]);
	  error->one(FLERR,str);
	}
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"ave") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/time command");
      if (strcmp(arg[iarg+1],"one") == 0) ave = ONE;
      else if (strcmp(arg[iarg+1],"running") == 0) ave = RUNNING;
      else if (strcmp(arg[iarg+1],"window") == 0) ave = WINDOW;
      else error->all(FLERR,"Illegal fix ave/time command");
      if (ave == WINDOW) {
	if (iarg+3 > narg) error->all(FLERR,"Illegal fix ave/time command");
	nwindow = atoi(arg[iarg+2]);
	if (nwindow <= 0) error->all(FLERR,"Illegal fix ave/time command");
      }
      iarg += 2;
      if (ave == WINDOW) iarg++;
    } else if (strcmp(arg[iarg],"start") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/time command");
      startstep = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"mode") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/time command");
      if (strcmp(arg[iarg+1],"scalar") == 0) mode = SCALAR;
      else if (strcmp(arg[iarg+1],"vector") == 0) mode = VECTOR;
      else error->all(FLERR,"Illegal fix ave/time command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"off") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/time command");
      memory->grow(offlist,noff+1,"ave/time:offlist");
      offlist[noff++] = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title1") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/spatial command");
      delete [] title1;
      int n = strlen(arg[iarg+1]) + 1;
      title1 = new char[n];
      strcpy(title1,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title2") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/spatial command");
      delete [] title2;
      int n = strlen(arg[iarg+1]) + 1;
      title2 = new char[n];
      strcpy(title2,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title3") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/spatial command");
      delete [] title3;
      int n = strlen(arg[iarg+1]) + 1;
      title3 = new char[n];
      strcpy(title3,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix ave/time command");
  }
}

/* ----------------------------------------------------------------------
   reallocate vectors for each input value, of length N
------------------------------------------------------------------------- */

void FixAveTime::grow()
{
  maxvalues += DELTA;
  memory->grow(which,maxvalues,"ave/time:which");
  memory->grow(argindex,maxvalues,"ave/time:argindex");
  memory->grow(value2index,maxvalues,"ave/time:value2index");
  memory->grow(offcol,maxvalues,"ave/time:offcol");
  ids = (char **) memory->srealloc(ids,maxvalues*sizeof(char *),"ave/time:ids");
}

/* ----------------------------------------------------------------------
   calculate nvalid = next step on which end_of_step does something
   can be this timestep if multiple of nfreq and nrepeat = 1
   else backup from next multiple of nfreq
   startstep is lower bound on nfreq multiple
------------------------------------------------------------------------- */

bigint FixAveTime::nextvalid()
{
  bigint nvalid = (update->ntimestep/nfreq)*nfreq + nfreq;
  while (nvalid < startstep) nvalid += nfreq;
  if (nvalid-nfreq == update->ntimestep && nrepeat == 1)
    nvalid = update->ntimestep;
  else
    nvalid -= (nrepeat-1)*nevery;
  if (nvalid < update->ntimestep) nvalid += nfreq;
  return nvalid;
}

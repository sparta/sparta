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

#include "sptype.h"
#include "stdlib.h"
#include "string.h"
#include "fix_ave_grid.h"
#include "grid.h"
#include "particle.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{STANDARD,COMPUTE,FIX,VARIABLE};
enum{ONE,RUNNING};

#define INVOKED_PER_GRID 16

/* ---------------------------------------------------------------------- */

FixAveGrid::FixAveGrid(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal fix ave/grid command");

  nevery = atoi(arg[2]);
  nrepeat = atoi(arg[3]);
  per_grid_freq = atoi(arg[4]);

  // scan values, then read options

  nvalues = 0;

  int iarg = 5;
  while (iarg < narg) {
    if ((strcmp(arg[iarg],"standard") == 0) ||
	(strncmp(arg[iarg],"c_",2) == 0) || 
	(strncmp(arg[iarg],"f_",2) == 0) || 
	(strncmp(arg[iarg],"v_",2) == 0)) {
      nvalues++;
      iarg++;
    } else break;
  }
  if (nvalues == 0) error->all(FLERR,"Illegal fix ave/grid command");

  options(narg,arg);

  // parse values until one isn't recognized
  // treat value = "standard" as a special case

  which = argindex = value2index = NULL;
  ids = NULL;
  allocate_values(nvalues);
  nvalues = 0;
  standard = 0;

  iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"standard") == 0) {
      if (nvalues != 0) error->all(FLERR,"Illegal fix ave/grid command");
      standard = 1;
      nvalues = 7 * particle->nspecies;
      allocate_values(nvalues);
      for (int i = 0; i < nvalues; i++) {
	which[i] = STANDARD;
	ids[i] = NULL;
      }

    } else if (strncmp(arg[iarg],"c_",2) == 0 || 
	       strncmp(arg[iarg],"f_",2) == 0 || 
	       strncmp(arg[iarg],"v_",2) == 0) {

      if (standard) error->all(FLERR,"Illegal fix ave/grid command");
      if (arg[iarg][0] == 'c') which[nvalues] = COMPUTE;
      else if (arg[iarg][0] == 'f') which[nvalues] = FIX;
      else if (arg[iarg][0] == 'v') which[nvalues] = VARIABLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
	if (suffix[strlen(suffix)-1] != ']')
	  error->all(FLERR,"Illegal fix ave/grid command");
	argindex[nvalues] = atoi(ptr+1);
	*ptr = '\0';
      } else argindex[nvalues] = 0;

      n = strlen(suffix) + 1;
      ids[nvalues] = new char[n];
      strcpy(ids[nvalues],suffix);
      nvalues++;
      delete [] suffix;

    } else break;

    iarg++;
  }

  // setup and error check
  // for fix inputs, check that fix frequency is acceptable

  if (nevery <= 0 || nrepeat <= 0 || per_grid_freq <= 0)
    error->all(FLERR,"Illegal fix ave/grid command");
  if (per_grid_freq % nevery || (nrepeat-1)*nevery >= per_grid_freq)
    error->all(FLERR,"Illegal fix ave/grid command");

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
	error->all(FLERR,"Compute ID for fix ave/grid does not exist");
      if (modify->compute[icompute]->per_grid_flag == 0)
	error->all(FLERR,
		   "Fix ave/grid compute does not calculate per-grid values");
      if (argindex[i] == 0 && 
	  modify->compute[icompute]->size_per_grid_cols != 0)
	error->all(FLERR,"Fix ave/grid compute does not "
		   "calculate a per-grid vector");
      if (argindex[i] && modify->compute[icompute]->size_per_grid_cols == 0)
	error->all(FLERR,"Fix ave/grid compute does not "
		   "calculate a per-grid array");
      if (argindex[i] && 
	  argindex[i] > modify->compute[icompute]->size_per_grid_cols)
	error->all(FLERR,"Fix ave/grid compute array is accessed out-of-range");

    } else if (which[i] == FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
	error->all(FLERR,"Fix ID for fix ave/grid does not exist");
      if (modify->fix[ifix]->per_grid_flag == 0)
	error->all(FLERR,"Fix ave/grid fix does not calculate per-grid values");
      if (argindex[i] == 0 && modify->fix[ifix]->size_per_grid_cols != 0)
	error->all(FLERR,
		   "Fix ave/grid fix does not calculate a per-grid vector");
      if (argindex[i] && modify->fix[ifix]->size_per_grid_cols == 0)
	error->all(FLERR,
		   "Fix ave/grid fix does not calculate a per-grid array");
      if (argindex[i] && argindex[i] > modify->fix[ifix]->size_per_grid_cols)
	error->all(FLERR,"Fix ave/grid fix array is accessed out-of-range");
      if (nevery % modify->fix[ifix]->per_grid_freq)
	error->all(FLERR,
		   "Fix for fix ave/grid not computed at compatible time");

    } else if (which[i] == VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
	error->all(FLERR,"Variable name for fix ave/grid does not exist");
      if (input->variable->grid_style(ivariable) == 0)
	error->all(FLERR,"Fix ave/grid variable is not grid-style variable");
    }
  }

  // this fix produces either a per-grid vector or array

  per_grid_flag = 1;
  if (nvalues == 1) size_per_grid_cols = 0;
  else size_per_grid_cols = nvalues;

  // allocate accumulators
  // if ave = RUNNING, allocate extra set of accvec/accarray
  // mcount = molecule counts = only needed for standard

  nspecies = particle->nspecies;
  nglocal = grid->nlocal;

  if (standard) memory->create(mcount,nspecies,nglocal,"ave/grid:mcount");
  else mcount = NULL;

  if (nvalues == 1) memory->create(vector_grid,nglocal,"ave/grid:vector_grid");
  else memory->create(array_grid,nglocal,nvalues,"ave/grid:array_grid");
  
  if (ave == RUNNING) {
    if (nvalues == 1) memory->create(accvec,nglocal,"ave/grid:accvec");
    else memory->create(accarray,nglocal,nvalues,"ave/grid:accarray");
  } else {
    if (nvalues == 1) accvec = vector_grid;
    else accarray = array_grid;
  }

  // zero mcount,accvec,accarray one time if ave = RUNNING

  int i,isp,m;

  if (ave == RUNNING) {
    if (standard) {
      for (isp = 0; isp < nspecies; isp++)
	for (i = 0; i < nglocal; i++)
	  mcount[isp][i] = 0;
    }
    if (nvalues == 1) {
      for (i = 0; i < nglocal; i++)
	accvec[i] = 0.0;
    } else {
      for (i = 0; i < nglocal; i++)
	for (m = 0; m < nvalues; m++)
	  accarray[i][m] = 0.0;
    }
  }

  // zero vector/array since dump may access it on timestep 0
  // zero vector/array since a variable may access it before first run

  if (nvalues == 0) {
    for (i = 0; i < nglocal; i++)
      vector_grid[i] = 0.0;
  } else {
    for (i = 0; i < nglocal; i++)
      for (m = 0; m < nvalues; m++)
	array_grid[i][m] = 0.0;
  }

  // nvalid = next step on which end_of_step does something
  // add nvalid to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  nsample = 0;
  irepeat = 0;
  nvalid = nextvalid();
  modify->addstep_compute_all(nvalid);
}

/* ---------------------------------------------------------------------- */

FixAveGrid::~FixAveGrid()
{
  memory->destroy(which);
  memory->destroy(argindex);
  memory->destroy(value2index);
  for (int i = 0; i < nvalues; i++) delete [] ids[i];
  memory->sfree(ids);

  memory->destroy(mcount);
  if (nvalues == 1) memory->destroy(vector_grid);
  else memory->destroy(array_grid);
  if (ave == RUNNING) {
    if (nvalues == 1) memory->destroy(accvec);
    else memory->destroy(accarray);
  }
}

/* ---------------------------------------------------------------------- */

int FixAveGrid::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAveGrid::init()
{
  // nspecies cannot have changed since constructor

  if (standard && particle->nspecies != nspecies)
	error->all(FLERR,"Species count changed for fix ave/grid");

  // set indices and check validity of all computes,fixes,variables

  for (int m = 0; m < nvalues; m++) {
    if (which[m] == COMPUTE) {
      int icompute = modify->find_compute(ids[m]);
      if (icompute < 0)
	error->all(FLERR,"Compute ID for fix ave/grid does not exist");
      value2index[m] = icompute;
      
    } else if (which[m] == FIX) {
      int ifix = modify->find_fix(ids[m]);
      if (ifix < 0) 
	error->all(FLERR,"Fix ID for fix ave/grid does not exist");
      value2index[m] = ifix;

    } else if (which[m] == VARIABLE) {
      int ivariable = input->variable->find(ids[m]);
      if (ivariable < 0) 
	error->all(FLERR,"Variable name for fix ave/grid does not exist");
      value2index[m] = ivariable;

    } else value2index[m] = -1;
  }
}

/* ----------------------------------------------------------------------
   only does something if nvalid = current timestep
------------------------------------------------------------------------- */

void FixAveGrid::setup(int vflag)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveGrid::end_of_step()
{
  int i,j,m,n,isp,icol;

  // skip if not step which requires doing something

  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;

  // zero accumulators if ave = ONE and first sample

  if (ave == ONE && irepeat == 0) {
    if (standard) {
      for (isp = 0; isp < nspecies; isp++)
	for (i = 0; i < nglocal; i++)
	  mcount[isp][i] = 0;
    }
    if (nvalues == 1) {
      for (i = 0; i < nglocal; i++)
	accvec[i] = 0.0;
    } else {
      for (i = 0; i < nglocal; i++)
	for (m = 0; m < nvalues; m++)
	  accarray[i][m] = 0.0;
    }
  }

  // accumulate results of attributes,computes,fixes,variables
  // compute/fix/variable may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  if (standard) {
    Grid::OneCell *cells = grid->cells;
    Particle::Species *species = particle->species;
    Particle::OnePart *particles = particle->particles;
    int nlocal = particle->nlocal;

    int icell,ilocal;
    double *v;

    for (int i = 0; i < nlocal; i++) {
      isp = particles[i].ispecies;
      icell = particles[i].icell;
      ilocal = cells[icell].local;
      mcount[isp][ilocal]++;
      icol = isp*7;

      v = particles[i].v;
      accarray[ilocal][icol+1] += v[0];
      accarray[ilocal][icol+2] += v[1];
      accarray[ilocal][icol+3] += v[2];
      accarray[ilocal][icol+4] += v[0]*v[0];
      accarray[ilocal][icol+5] += v[1]*v[1];
      accarray[ilocal][icol+6] += v[2]*v[2];
    }

  } else {
    for (m = 0; m < nvalues; m++) {
      n = value2index[m];
      j = argindex[m];

      // invoke compute if not previously invoked

      if (which[m] == COMPUTE) {
	Compute *compute = modify->compute[n];
	if (!(compute->invoked_flag & INVOKED_PER_GRID)) {
	  compute->compute_per_grid();
	  compute->invoked_flag |= INVOKED_PER_GRID;
	}
	
	if (j == 0) {
	  double *compute_vector = compute->vector_grid;
	  if (nvalues == 1) {
	    for (i = 0; i < nglocal; i++)
	      accvec[i] += compute_vector[i];
	  } else {
	    for (i = 0; i < nglocal; i++)
	      accarray[i][m] += compute_vector[i];
	  }
	} else {
	  int jm1 = j - 1;
	  double **compute_array = compute->array_grid;
	  if (nvalues == 1) {
	    for (i = 0; i < nglocal; i++)
	      accvec[i] += compute_array[i][jm1];
	  } else {
	    for (i = 0; i < nglocal; i++)
	      accarray[i][m] += compute_array[i][jm1];
	  }
	}
	
      // access fix fields, guaranteed to be ready
	
      } else if (which[m] == FIX) {
	if (j == 0) {
	  double *fix_vector = modify->fix[n]->vector_grid;
	  if (nvalues == 1) {
	    for (i = 0; i < nglocal; i++)
	      accvec[i] += fix_vector[i];
	  } else {
	    for (i = 0; i < nglocal; i++)
	      accarray[i][m] += fix_vector[i];
	  }
	} else {
	  int jm1 = j - 1;
	  double **fix_array = modify->fix[n]->array_grid;
	  if (nvalues == 1) {
	    for (i = 0; i < nglocal; i++)
	      accvec[i] += fix_array[i][jm1];
	  } else {
	    for (i = 0; i < nglocal; i++)
	      accarray[i][m] += fix_array[i][jm1];
	  }
	}

      // evaluete grid-style variable

      } else if (which[m] == VARIABLE) {
      }
    }
  }

  // done if irepeat < nrepeat
  // else reset irepeat and nvalid

  nsample++;
  irepeat++;
  if (irepeat < nrepeat) {
    nvalid += nevery;
    modify->addstep_compute(nvalid);
    return;
  }

  irepeat = 0;
  nvalid = ntimestep+per_grid_freq - (nrepeat-1)*nevery;
  modify->addstep_compute(nvalid);

  // average the accumulators for output on Nfreq timestep
  // molecule count is normalized by nsample
  // molecule properties are normalized by molecule count
  // reset nsample if ave = ONE
  // NOTE: have to decide how to do averaging for non-standard case

  if (ave == ONE) {
    if (standard) {
      for (i = 0; i < nglocal; i++)
	for (isp = 0; isp < nspecies; isp++) {
	  icol = isp*7;
	  if (mcount[isp][i]) {
	    array_grid[i][icol+1] /= mcount[isp][i];
	    array_grid[i][icol+2] /= mcount[isp][i];
	    array_grid[i][icol+3] /= mcount[isp][i];
	    array_grid[i][icol+4] /= mcount[isp][i];
	    array_grid[i][icol+5] /= mcount[isp][i];
	    array_grid[i][icol+6] /= mcount[isp][i];
	  }
	}
      for (i = 0; i < nglocal; i++)
	for (isp = 0; isp < nspecies; isp++)
	  array_grid[i][isp*7] = 1.0*mcount[isp][i] / nsample;
    }
  } else {
    if (standard) {
      for (i = 0; i < nglocal; i++)
	for (isp = 0; isp < nspecies; isp++) {
	  icol = isp*7;
	  if (mcount[isp][i]) {
	    array_grid[i][icol+1] = accarray[i][icol+1] / mcount[isp][i];
	    array_grid[i][icol+2] = accarray[i][icol+2] / mcount[isp][i];
	    array_grid[i][icol+3] = accarray[i][icol+3] / mcount[isp][i];
	    array_grid[i][icol+4] = accarray[i][icol+4] / mcount[isp][i];
	    array_grid[i][icol+5] = accarray[i][icol+5] / mcount[isp][i];
	    array_grid[i][icol+6] = accarray[i][icol+6] / mcount[isp][i];
	  } else {
	    array_grid[i][icol+1] = 0.0;
	    array_grid[i][icol+2] = 0.0;
	    array_grid[i][icol+3] = 0.0;
	    array_grid[i][icol+4] = 0.0;
	    array_grid[i][icol+5] = 0.0;
	    array_grid[i][icol+6] = 0.0;
	  }
	}
      for (i = 0; i < nglocal; i++)
	for (isp = 0; isp < nspecies; isp++)
	  array_grid[i][isp*7] = 1.0*mcount[isp][i] / nsample;
    }
  }

  if (ave == ONE) nsample = 0;
}

/* ----------------------------------------------------------------------
   parse optional args
------------------------------------------------------------------------- */

void FixAveGrid::options(int narg, char **arg)
{
  // option defaults

  ave = ONE;

  // optional args

  int iarg = 5 + nvalues;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"ave") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/grid command");
      if (strcmp(arg[iarg+1],"one") == 0) ave = ONE;
      else if (strcmp(arg[iarg+1],"running") == 0) ave = RUNNING;
      else error->all(FLERR,"Illegal fix ave/grid command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix ave/grid command");
  }
}

/* ----------------------------------------------------------------------
   reallocate vectors for each input value, of length N
------------------------------------------------------------------------- */

void FixAveGrid::allocate_values(int n)
{
  memory->grow(which,n,"ave/grid:which");
  memory->grow(argindex,n,"ave/grid:argindex");
  memory->grow(value2index,n,"ave/grid:value2index");
  ids = (char **) memory->srealloc(ids,n*sizeof(char *),"ave/grid:ids");
}

/* ----------------------------------------------------------------------
   memory usage of accumulators
------------------------------------------------------------------------- */

double FixAveGrid::memory_usage()
{
  double bytes = 0.0;
  if (standard) bytes += nglocal*nspecies * sizeof(int);
  bytes += nglocal*nvalues * sizeof(double);
  if (ave == RUNNING) bytes += nglocal*nvalues * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   calculate nvalid = next step on which end_of_step does something
   can be this timestep if multiple of nfreq and nrepeat = 1
   else backup from next multiple of nfreq
------------------------------------------------------------------------- */

bigint FixAveGrid::nextvalid()
{
  bigint nvalid = (update->ntimestep/per_grid_freq)*per_grid_freq + 
    per_grid_freq;
  if (nvalid-per_grid_freq == update->ntimestep && nrepeat == 1)
    nvalid = update->ntimestep;
  else
    nvalid -= (nrepeat-1)*nevery;
  if (nvalid < update->ntimestep) nvalid += per_grid_freq;
  return nvalid;
}

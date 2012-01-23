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

#include "dsmctype.h"
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

using namespace DSMC_NS;

enum{COMPUTE,FIX,VARIABLE};

#define STANDARD 6

#define INVOKED_PER_GRID 16

/* ---------------------------------------------------------------------- */

FixAveGrid::FixAveGrid(DSMC *dsmc, int narg, char **arg) :
  Fix(dsmc, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal fix ave/grid command");

  nevery = atoi(arg[2]);
  nrepeat = atoi(arg[3]);
  per_grid_freq = atoi(arg[4]);

  // parse remaining values

  which = new int[narg-5];
  argindex = new int[narg-5];
  ids = new char*[narg-5];
  value2index = new int[narg-5];
  nvalues = 0;

  int iarg = 5;
  while (iarg < narg) {
    ids[nvalues] = NULL;

    if (strncmp(arg[iarg],"c_",2) == 0 || 
	       strncmp(arg[iarg],"f_",2) == 0 || 
	       strncmp(arg[iarg],"v_",2) == 0) {
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

    } else error->all(FLERR,"Illegal fix ave/grid command");

    iarg++;
  }

  // for now, allow no explicit values

  if (nvalues) error->all(FLERR,"Illegal fix ave/grid command");

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
      if (input->variable->cell_style(ivariable) == 0)
	error->all(FLERR,"Fix ave/grid variable is not cell-style variable");
    }
  }

  // this fix produces either a per-grid vector or array

  per_grid_flag = 1;
  if (nvalues == 0) size_per_grid_cols = STANDARD;
  else if (nvalues == 1) size_per_grid_cols = 0;
  else size_per_grid_cols = nvalues;

  // perform initial allocation of cell-based count and vector/array

  int nglocal = grid->nlocal;
  memory->create(pcount,nglocal,"ave/time:pcount");

  if (nvalues == 0)
    memory->create(array_grid,nglocal,STANDARD,"ave/time:array_grid");
  else if (nvalues == 1) 
    memory->create(vector_grid,nglocal,"ave/time:vector_grid");
  else
    memory->create(array_grid,nglocal,nvalues,"ave/time:array_grid");
  
  // zero vector/array since dump may access it on timestep 0
  // zero vector/array since a variable may access it before first run

  if (nvalues == 0) {
    for (int i = 0; i < nglocal; i++)
      for (int m = 0; m < STANDARD; m++)
	array_grid[i][m] = 0.0;
  } else if (nvalues == 0) {
    for (int i = 0; i < nglocal; i++)
      vector_grid[i] = 0.0;
  } else {
    for (int i = 0; i < nglocal; i++)
      for (int m = 0; m < nvalues; m++)
	array_grid[i][m] = 0.0;
  }

  // nvalid = next step on which end_of_step does something

  irepeat = 0;
  nvalid = nextvalid();
}

/* ---------------------------------------------------------------------- */

FixAveGrid::~FixAveGrid()
{
  delete [] which;
  delete [] argindex;
  for (int m = 0; m < nvalues; m++) delete [] ids[m];
  delete [] ids;
  delete [] value2index;

  memory->destroy(pcount);
  if (nvalues == 0) memory->destroy(array_grid);
  else if (nvalues == 1) memory->destroy(vector_grid);
  else memory->destroy(array_grid);
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
  int i,j,m,n;

  // skip if not step which requires doing something

  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;

  // zero vector/array if first step

  int nglocal = grid->nlocal;

  if (irepeat == 0) {
    for (i = 0; i < nglocal; i++) pcount[i] = 0;
    if (nvalues == 0) {
      for (i = 0; i < nglocal; i++)
	for (m = 0; m < STANDARD; m++)
	  array_grid[i][m] = 0.0;
    } else if (nvalues == 1) {
      for (i = 0; i < nglocal; i++)
	vector_grid[i] = 0.0;
    } else {
      for (i = 0; i < nglocal; i++)
	for (m = 0; m < nvalues; m++)
	  array_grid[i][m] = 0.0;
    }
  }

  // accumulate results of attributes,computes,fixes,variables to local copy

  if (nvalues == 0) {
    Grid::OneCell *cells = grid->cells;
    Particle::Species *species = particle->species;
    Particle::OnePart *particles = particle->particles;
    int nlocal = particle->nlocal;

    int icell,ilocal;
    double *v,*row;

    for (int i = 0; i < nlocal; i++) {
      icell = particles[i].icell;
      ilocal = cells[icell].local;
      pcount[ilocal]++;

      v = particles[i].v;
      row = array_grid[ilocal];
      row[0] += v[0];
      row[1] += v[1];
      row[2] += v[2];
      row[3] += v[0]*v[0];
      row[4] += v[1]*v[1];
      row[5] += v[2]*v[2];
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
	      vector_grid[i] += compute_vector[i];
	  } else {
	    for (i = 0; i < nglocal; i++)
	      array_grid[i][m] += compute_vector[i];
	  }
	} else {
	  int jm1 = j - 1;
	  double **compute_array = compute->array_grid;
	  if (nvalues == 1) {
	    for (i = 0; i < nglocal; i++)
	      vector_grid[i] += compute_array[i][jm1];
	  } else {
	    for (i = 0; i < nglocal; i++)
	      array_grid[i][m] += compute_array[i][jm1];
	  }
	}
	
      // access fix fields, guaranteed to be ready
	
      } else if (which[m] == FIX) {
	if (j == 0) {
	  double *fix_vector = modify->fix[n]->vector_grid;
	  if (nvalues == 1) {
	    for (i = 0; i < nglocal; i++)
	      vector_grid[i] += fix_vector[i];
	  } else {
	    for (i = 0; i < nglocal; i++)
	      array_grid[i][m] += fix_vector[i];
	  }
	} else {
	  int jm1 = j - 1;
	  double **fix_array = modify->fix[n]->array_grid;
	  if (nvalues == 1) {
	    for (i = 0; i < nglocal; i++)
	      vector_grid[i] += fix_array[i][jm1];
	  } else {
	    for (i = 0; i < nglocal; i++)
	      array_grid[i][m] += fix_array[i][jm1];
	  }
	}

      // evaluete cell-style variable

      } else if (which[m] == VARIABLE) {
      }
    }
  }

  // done if irepeat < nrepeat
  // else reset irepeat and nvalid

  irepeat++;
  if (irepeat < nrepeat) {
    nvalid += nevery;
    return;
  }

  irepeat = 0;
  nvalid = ntimestep+per_grid_freq - (nrepeat-1)*nevery;

  // average the final result for the Nfreq timestep

  if (nvalues == 0) {
    for (i = 0; i < nglocal; i++)
      for (m = 0; m < STANDARD; m++)
	array_grid[i][m] /= pcount[i];
  } else if (nvalues == 1) {
    for (i = 0; i < nglocal; i++)
      vector_grid[i] /= pcount[i];
  } else {
    for (i = 0; i < nglocal; i++)
      for (m = 0; m < nvalues; m++)
	array_grid[i][m] /= pcount[i];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local cell-based array
------------------------------------------------------------------------- */

double FixAveGrid::memory_usage()
{
  double bytes = grid->nlocal * sizeof(int);
  if (nvalues == 0) bytes += grid->nlocal*STANDARD * sizeof(double);
  else if (nvalues == 1) bytes += grid->nlocal * sizeof(double);
  else bytes += grid->nlocal*nvalues * sizeof(double);
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

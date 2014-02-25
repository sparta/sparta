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
#include "fix_ave_grid.h"
#include "grid.h"
#include "particle.h"
#include "comm.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{COMPUTE,FIX,VARIABLE};
enum{ONE,RUNNING};

// NOTE: set DELTAGRID to big value when done debugging

#define INVOKED_PER_GRID 16
#define DELTAINPUT 8
#define DELTAGRID 10              // must be bigger than split cells per cell

/* ---------------------------------------------------------------------- */

FixAveGrid::FixAveGrid(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal fix ave/grid command");

  nevery = atoi(arg[2]);
  nrepeat = atoi(arg[3]);
  per_grid_freq = atoi(arg[4]);
  gridmigrate = 1;

  // scan values, then read options

  int iarg = 5;
  while (iarg < narg) {
    if ((strncmp(arg[iarg],"c_",2) == 0) || 
	(strncmp(arg[iarg],"f_",2) == 0) || 
	(strncmp(arg[iarg],"v_",2) == 0)) iarg++;
    else break;
  }

  options(narg-iarg,&arg[iarg]);

  // parse values until one isn't recognized
  // expand compute or fix array into full set of columns

  which = argindex = value2index = NULL;
  ids = NULL;
  nvalues = maxvalues = 0;

  iarg = 5;
  while (iarg < narg) {
    if (strncmp(arg[iarg],"c_",2) == 0 || 
        strncmp(arg[iarg],"f_",2) == 0 || 
        strncmp(arg[iarg],"v_",2) == 0) {

      if (nvalues == maxvalues) grow();
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

      if ((which[nvalues] == COMPUTE || which[nvalues] == FIX) && 
          argindex[nvalues] == 0) {
        int ndup = 1;
        if (which[nvalues] == COMPUTE) {
          int icompute = modify->find_compute(ids[nvalues]);
          if (icompute >= 0)
            if (modify->compute[icompute]->per_grid_flag)
              ndup = modify->compute[icompute]->size_per_grid_cols;
        }
        if (which[nvalues] == FIX) {
          int ifix = modify->find_fix(ids[nvalues]);
          if (ifix >= 0)
            if (modify->fix[ifix]->per_grid_flag)
              ndup = modify->fix[ifix]->size_per_grid_cols;
        }
        if (ndup > 1) {
          argindex[nvalues] = 1;
          nvalues++;
          for (int icol = 2; icol <= ndup; icol++) {
            if (nvalues == maxvalues) grow();
            which[nvalues] = which[nvalues-1];
            argindex[nvalues] = icol;
            n = strlen(suffix) + 1;
            ids[nvalues] = new char[n];
            strcpy(ids[nvalues],suffix);
            nvalues++;
          }
          nvalues--;
        }
      }

      nvalues++;
      delete [] suffix;

      iarg++;
    } else break;
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

  vector_grid = vector = NULL;
  array_grid = array = NULL;

  // setup norm vectors and norm pointers
  // only store unique norms by checking if 
  //   current compute & flag matches previous computes & flags
  // NOTE: add more logic for fixes and variables if enable them

  normacc = new int[nvalues];
  normindex = new int[nvalues];
  norms = new double*[nvalues];
  int *list_compute = new int[nvalues];
  int *list_flag = new int[nvalues];

  nnorm = 0;
  for (int m = 0; m < nvalues; m++) {
    int j = argindex[m];
    
    normacc[m] = 0;
    normindex[m] = -1;
    
    if (which[m] == COMPUTE) {
      int icompute = modify->find_compute(ids[m]);
      Compute *compute = modify->compute[icompute];
      int flag = compute->normflag(j);
      if (!flag) continue;
      
      int i;
      for (i = 0; i < nnorm; i++)
        if (icompute == list_compute[i] && flag == list_flag[i]) break;
      if (i < nnorm) normindex[m] = i;
      else {
        normacc[m] = 1;
        normindex[m] = nnorm;
        norms[nnorm] = NULL;
        list_compute[nnorm] = icompute;
        list_flag[nnorm] = flag;
        nnorm++;
      }
    }
  }

  delete [] list_compute;
  delete [] list_flag;

  // allocate per-grid cell memory for vectors/arrays and norms

  nglocal = nglocalmax = grid->nlocal;
  allocate(nglocal);

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

  if (nvalues == 1) memory->destroy(vector_grid);
  else memory->destroy(array_grid);
  if (ave == RUNNING) {
    if (nvalues == 1) memory->destroy(vector);
    else memory->destroy(array);
  }

  delete [] normacc;
  delete [] normindex;
  for (int i = 0; i < nnorm; i++) memory->destroy(norms[i]);
  delete [] norms;
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

void FixAveGrid::setup()
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveGrid::end_of_step()
{
  int i,j,m,n,isp,icol;
  double *norm,*cfv_norm;

  // skip if not step which requires doing something

  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;

  // zero vector/array and norms if ave = ONE and first sample

  if (ave == ONE && irepeat == 0) {
    if (nvalues == 1)
      for (i = 0; i < nglocal; i++)
	vector[i] = 0.0;
    else
      for (i = 0; i < nglocal; i++)
	for (m = 0; m < nvalues; m++)
	  array[i][m] = 0.0;
    for (int m = 0; m < nnorm; m++) {
      norm = norms[m];
      for (i = 0; i < nglocal; i++) norm[i] = 0.0;
    }
  }

  // accumulate results of computes,fixes,variables
  // compute/fix/variable may invoke computes so wrap with clear/add
  // accumulate norm values if necessary, only done once per unique norm
  // NOTE: add more logic for fixes and variables if enable them

  modify->clearstep_compute();

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
        if (nvalues == 1)
          for (i = 0; i < nglocal; i++)
            vector[i] += compute_vector[i];
        else
          for (i = 0; i < nglocal; i++)
            array[i][m] += compute_vector[i];
      } else {
        int jm1 = j - 1;
        double **compute_array = compute->array_grid;
        if (nvalues == 1)
          for (i = 0; i < nglocal; i++)
            vector[i] += compute_array[i][jm1];
        else
          for (i = 0; i < nglocal; i++)
            array[i][m] += compute_array[i][jm1];
      }

      if (normacc[m]) {
        cfv_norm = compute->normptr(j);
        norm = norms[normindex[m]];
        for (i = 0; i < nglocal; i++) norm[i] += cfv_norm[i];
      }
      
    // access fix fields, guaranteed to be ready
      
    } else if (which[m] == FIX) {
      if (j == 0) {
        double *fix_vector = modify->fix[n]->vector_grid;
        if (nvalues == 1)
          for (i = 0; i < nglocal; i++)
            vector[i] += fix_vector[i];
        else
          for (i = 0; i < nglocal; i++)
            array[i][m] += fix_vector[i];
      } else {
        int jm1 = j - 1;
        double **fix_array = modify->fix[n]->array_grid;
        if (nvalues == 1)
          for (i = 0; i < nglocal; i++)
            vector[i] += fix_array[i][jm1];
        else
          for (i = 0; i < nglocal; i++)
            array[i][m] += fix_array[i][jm1];
      }
      
      // evaluate grid-style variable
      
    } else if (which[m] == VARIABLE) {
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

  // normalize for output on Nfreq timestep
  // normindex < 0, just normalize by # of samples
  // normindex >= 0, normalize by accumulated norm vector

  if (ave == ONE) {
    if (nvalues == 1) {
      if (normindex[0] < 0)
        for (i = 0; i < nglocal; i++) vector_grid[i] /= nsample;
      else {
        norm = norms[normindex[0]];
        for (i = 0; i < nglocal; i++)
          if (norm[i] > 0.0) vector_grid[i] /= norm[i];
      }
    } else {
      for (m = 0; m < nvalues; m++) {
        if (normindex[m] < 0)
          for (i = 0; i < nglocal; i++) array_grid[i][m] /= nsample;
        else {
          norm = norms[normindex[m]];
          for (i = 0; i < nglocal; i++)
            if (norm[i] > 0.0) array_grid[i][m] /= norm[i];
        }
      }
    }

  } else {
    if (nvalues == 1) {
      if (normindex[0] < 0)
        for (i = 0; i < nglocal; i++) vector_grid[i] = vector[i]/nsample;
      else {
        norm = norms[normindex[0]];
        for (i = 0; i < nglocal; i++) 
          if (norm[i] > 0.0) vector_grid[i] = vector[i]/norm[i];
      }
    } else {
      for (m = 0; m < nvalues; m++) {
        if (normindex[m] < 0)
          for (i = 0; i < nglocal; i++)
            array_grid[i][m] = array[i][m]/nsample;
        else {
          norm = norms[normindex[m]];
          for (i = 0; i < nglocal; i++) 
            if (norm[i] > 0.0) array_grid[i][m] = array[i][m]/norm[i];
        }
      }
    }
  }

  // reset nsample if ave = ONE

  if (ave == ONE) nsample = 0;
}

/* ----------------------------------------------------------------------
   pack icell values for per-cell arrays into buf
   if icell is a split cell, also pack all sub cell values 
   return byte count of amount packed
   if memflag, only return count, do not fill buf
------------------------------------------------------------------------- */

int FixAveGrid::pack_grid_one(int icell, char *buf, int memflag)
{
  char *ptr = buf;

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  Grid::SplitInfo *sinfo = grid->sinfo;

  ptr += pack_one(icell,ptr,memflag);

  if (cells[icell].nsplit > 1) {
    int isplit = cells[icell].isplit;
    int nsplit = cells[icell].nsplit;
    for (int i = 0; i < nsplit; i++)
      ptr += pack_one(sinfo[isplit].csubs[i],ptr,memflag);
  }

  return ptr-buf;
}
 
/* ----------------------------------------------------------------------
   pack one set of values into buf from icell
------------------------------------------------------------------------- */

int FixAveGrid::pack_one(int icell, char *buf, int memflag)
{
  char *ptr = buf;

  if (memflag) {
    if (nvalues == 1) *((double *) ptr) = vector_grid[icell];
    else memcpy(ptr,array_grid[icell],nvalues*sizeof(double));
  }
  ptr += nvalues*sizeof(double);

  if (ave == RUNNING) {
    if (memflag) {
      if (nvalues == 1) *((double *) ptr) = vector[icell];
      else memcpy(ptr,array[icell],nvalues*sizeof(double));
    }
    ptr += nvalues*sizeof(double);
  }

  if (memflag) {
    double *dptr = (double *) ptr;
    for (int i = 0; i < nnorm; i++) {
      *dptr = norms[i][icell];
      dptr++;
    }
  }
  ptr += nnorm*sizeof(double);

  return ptr-buf;
}

/* ----------------------------------------------------------------------
   unpack icell values for per-cell arrays from buf
   if icell is a split cell, also unpack all sub cell values 
   return byte count of amount unpacked
------------------------------------------------------------------------- */

int FixAveGrid::unpack_grid_one(int icell, char *buf)
{
  char *ptr = buf;

  Grid::ChildCell *cells = grid->cells;
  Grid::SplitInfo *sinfo = grid->sinfo;

  grow_percell(1);
  ptr += unpack_one(ptr,icell);
  nglocal++;

  if (cells[icell].nsplit > 1) {
    int isplit = cells[icell].isplit;
    int nsplit = cells[icell].nsplit;
    grow_percell(nsplit);
    for (int i = 0; i < nsplit; i++)
      ptr += unpack_one(ptr,sinfo[isplit].csubs[i]);
    nglocal += nsplit;
  }

  return ptr-buf;
}

/* ----------------------------------------------------------------------
   unpack one set of values from buf into icell
------------------------------------------------------------------------- */

int FixAveGrid::unpack_one(char *buf, int icell)
{
  char *ptr = buf;

  if (nvalues == 1) vector_grid[icell] = *((double *) ptr);
  else memcpy(array_grid[icell],ptr,nvalues*sizeof(double));
  ptr += nvalues*sizeof(double);

  if (ave == RUNNING) {
    if (nvalues == 1) vector[icell] = *((double *) ptr);
    else memcpy(array[icell],ptr,nvalues*sizeof(double));
    ptr += nvalues*sizeof(double);
  }

  double *dptr = (double *) ptr;
  for (int i = 0; i < nnorm; i++) {
    norms[i][icell] = *dptr;
    dptr++;
  }
  ptr += nnorm*sizeof(double);

  return ptr-buf;
}

/* ----------------------------------------------------------------------
   compress per-cell arrays due to cells migrating to new procs
   criteria for keeping/discarding a cell is same as in Grid::compress()
   this keeps final ordering of per-cell arrays consistent with Grid class
------------------------------------------------------------------------- */

void FixAveGrid::compress_grid()
{
  int me = comm->me;
  Grid::ChildCell *cells = grid->cells;

  // keep an unsplit or split cell if staying on this proc
  // keep a sub cell if its split cell is staying on this proc

  int ncurrent = nglocal;
  nglocal = 0;
  for (int icell = 0; icell < ncurrent; icell++) {
    if (cells[icell].nsplit >= 1) {
      if (cells[icell].proc != me) continue;
    } else {
      int isplit = cells[icell].isplit;
      if (cells[grid->sinfo[isplit].icell].proc != me) continue;
    }

    if (nglocal != icell)  {
      if (nvalues == 1) vector_grid[nglocal] = vector_grid[icell];
      else memcpy(array_grid[nglocal],array_grid[icell],nvalues*sizeof(double));
      if (ave == RUNNING) {
        if (nvalues == 1) vector[nglocal] = vector[icell];
        else memcpy(array[nglocal],array[icell],nvalues*sizeof(double));
      }
      for (int m = 0; m < nnorm; m++) norms[m][nglocal] = norms[m][icell];
    }
    nglocal++;
  }
}

/* ----------------------------------------------------------------------
   memory usage of accumulators
------------------------------------------------------------------------- */

double FixAveGrid::memory_usage()
{
  double bytes = 0.0;
  bytes += nglocalmax*nvalues * sizeof(double);
  if (ave == RUNNING) bytes += nglocalmax*nvalues * sizeof(double);
  bytes += nnorm*nglocalmax * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   parse optional args
------------------------------------------------------------------------- */

void FixAveGrid::options(int narg, char **arg)
{
  // option defaults

  ave = ONE;

  // optional args

  int iarg = 0;
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
   grow vectors for each input value
------------------------------------------------------------------------- */

void FixAveGrid::grow()
{
  maxvalues += DELTAINPUT;
  memory->grow(which,maxvalues,"ave/grid:which");
  memory->grow(argindex,maxvalues,"ave/grid:argindex");
  memory->grow(value2index,maxvalues,"ave/grid:value2index");
  ids = (char **) memory->srealloc(ids,maxvalues*sizeof(char *),"ave/grid:ids");
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

/* ----------------------------------------------------------------------
   allocate per-grid vectors/arrays/norms
   zero newly allocated values in case used by dump or load balancer
------------------------------------------------------------------------- */

void FixAveGrid::allocate(int n)
{
  // if ave = ONE, allocate only vector_grid and array_grid
  // if ave = RUNNING, also allocate vector and array

  int m;

  if (nvalues == 1) {
    memory->grow(vector_grid,n,"ave/grid:vector_grid");
    for (int i = 0; i < n; i++) vector_grid[i] = 0.0;
  } else {
    memory->grow(array_grid,n,nvalues,"ave/grid:array_grid");
    for (int i = 0; i < n; i++)
      for (m = 0; m < nvalues; m++) array_grid[i][m] = 0.0;
  }
  
  if (ave == RUNNING) {
    if (nvalues == 1) {
      memory->grow(vector,n,"ave/grid:vector");
      for (int i = 0; i < n; i++) vector[i] = 0.0;
    } else {
      memory->grow(array,n,nvalues,"ave/grid:array");
      for (int i = 0; i < n; i++)
        for (m = 0; m < nvalues; m++) array[i][m] = 0.0;
    }
  } else {
    if (nvalues == 1) vector = vector_grid;
    else array = array_grid;
  }
  
  for (m = 0; m < nnorm; m++) {
    memory->grow(norms[m],n,"ave/grid:norms");
    double *norm = norms[m];
    for (int i = 0; i < n; i++) norm[i] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   insure per-cell arrays are allocated long enough for N new cells
------------------------------------------------------------------------- */

void FixAveGrid::grow_percell(int nnew)
{
  if (nglocal+nnew < nglocalmax) return;
  nglocalmax += DELTAGRID;
  int n = nglocalmax;

  if (nvalues == 1) memory->grow(vector_grid,n,"ave/grid:vector_grid");
  else memory->grow(array_grid,n,nvalues,"ave/grid:array_grid");
  
  if (ave == RUNNING) {
    if (nvalues == 1) memory->grow(vector,n,"ave/grid:vector");
    else memory->grow(array,n,nvalues,"ave/grid:array");
  } else {
    if (nvalues == 1) vector = vector_grid;
    else array = array_grid;
  }
  
  for (int m = 0; m < nnorm; m++)
    memory->grow(norms[m],n,"ave/grid:norms");
}

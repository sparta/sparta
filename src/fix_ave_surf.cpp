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
#include "stdlib.h"
#include "string.h"
#include "fix_ave_surf.h"
#include "surf.h"
#include "particle.h"
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "compute.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{COMPUTE,FIX,VARIABLE};
enum{ONE,RUNNING};

#define INVOKED_PER_SURF 32
#define DELTA 8;

/* ---------------------------------------------------------------------- */

FixAveSurf::FixAveSurf(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal fix ave/surf command");

  int igroup = surf->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Could not find fix ave/surf group ID");
  groupbit = surf->bitmask[igroup];

  nevery = atoi(arg[3]);
  nrepeat = atoi(arg[4]);
  per_surf_freq = atoi(arg[5]);

  time_depend = 1;

  // scan values, then read options

  nvalues = 0;

  int iarg = 6;
  while (iarg < narg) {
    if ((strncmp(arg[iarg],"c_",2) == 0) || 
	(strncmp(arg[iarg],"f_",2) == 0) || 
	(strncmp(arg[iarg],"v_",2) == 0)) {
      nvalues++;
      iarg++;
    } else break;
  }

  if (nvalues == 0) error->all(FLERR,"No values in fix ave/surf command");

  options(iarg,narg,arg);

  // expand args if any have wildcard character "*"
  // this can reset nvalues

  int expand = 0;
  char **earg;
  nvalues = input->expand_args(nvalues,&arg[6],1,earg);

  if (earg != &arg[6]) expand = 1;
  arg = earg;

  // parse values

  which = new int[nvalues];
  argindex = new int[nvalues];
  value2index = new int[nvalues];
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
        error->all(FLERR,"Illegal fix ave/surf command");
      argindex[i] = atoi(ptr+1);
      *ptr = '\0';
    } else argindex[i] = 0;

    n = strlen(suffix) + 1;
    ids[i] = new char[n];
    strcpy(ids[i],suffix);
    delete [] suffix;
  }

  // if wildcard expansion occurred, free earg memory from expand_args()
  // wait to do this until after file comment lines are printed

  if (expand) {
    for (int i = 0; i < nvalues; i++) delete [] earg[i];
    memory->sfree(earg);
  }

  // setup and error check
  // for fix inputs, check that fix frequency is acceptable

  if (nevery <= 0 || nrepeat <= 0 || per_surf_freq <= 0)
    error->all(FLERR,"Illegal fix ave/surf command");
  if (per_surf_freq % nevery || (nrepeat-1)*nevery >= per_surf_freq)
    error->all(FLERR,"Illegal fix ave/surf command");

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
	error->all(FLERR,"Compute ID for fix ave/surf does not exist");
      if (modify->compute[icompute]->per_surf_flag == 0)
	error->all(FLERR,
		   "Fix ave/surf compute does not calculate per-surf values");
      if (argindex[i] == 0 && 
	  modify->compute[icompute]->size_per_surf_cols != 0)
	error->all(FLERR,"Fix ave/surf compute does not "
		   "calculate a per-surf vector");
      if (argindex[i] && modify->compute[icompute]->size_per_surf_cols == 0)
	error->all(FLERR,"Fix ave/surf compute does not "
		   "calculate a per-surf array");
      if (argindex[i] && 
	  argindex[i] > modify->compute[icompute]->size_per_surf_cols)
	error->all(FLERR,"Fix ave/surf compute array is accessed out-of-range");

    } else if (which[i] == FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
	error->all(FLERR,"Fix ID for fix ave/surf does not exist");
      if (modify->fix[ifix]->per_surf_flag == 0)
	error->all(FLERR,"Fix ave/surf fix does not calculate per-surf values");
      if (argindex[i] == 0 && modify->fix[ifix]->size_per_surf_cols != 0)
	error->all(FLERR,
		   "Fix ave/surf fix does not calculate a per-surf vector");
      if (argindex[i] && modify->fix[ifix]->size_per_surf_cols == 0)
	error->all(FLERR,
		   "Fix ave/surf fix does not calculate a per-surf array");
      if (argindex[i] && argindex[i] > modify->fix[ifix]->size_per_surf_cols)
	error->all(FLERR,"Fix ave/surf fix array is accessed out-of-range");
      if (nevery % modify->fix[ifix]->per_surf_freq)
	error->all(FLERR,
		   "Fix for fix ave/surf not computed at compatible time");

    } else if (which[i] == VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
	error->all(FLERR,"Variable name for fix ave/surf does not exist");
      if (input->variable->surf_style(ivariable) == 0)
	error->all(FLERR,"Fix ave/surf variable is not surf-style variable");
    }
  }

  // this fix produces either a per-surf vector or array

  per_surf_flag = 1;
  if (nvalues == 1) size_per_surf_cols = 0;
  else size_per_surf_cols = nvalues;

  // allocate accumulators for owned surfaces
  // if ave = RUNNING, allocate extra set of accvec/accarray

  nslocal = surf->nlocal;

  memory->create(buflocal,nslocal,"ave/surf:buflocal");
  memory->create(masks,nslocal,"ave/surf:masks");

  if (nvalues == 1) memory->create(vector_surf,nslocal,"ave/surf:vector_surf");
  else memory->create(array_surf,nslocal,nvalues,"ave/surf:array_surf");

  if (ave == RUNNING) {
    if (nvalues == 1) memory->create(accvec,nslocal,"ave/surf:accvec");
    else memory->create(accarray,nslocal,nvalues,"ave/surf:accarray");
  } else {
    if (nvalues == 1) accvec = vector_surf;
    else accarray = array_surf;
  }

  // allocate accumulators for local surfaces

  if (domain->dimension == 2) nsurf = surf->nline;
  else nsurf = surf->ntri;

  memory->create(glob2loc,nsurf,"surf:glob2loc");
  for (int i = 0; i < nsurf; i++) glob2loc[i] = -1;

  nlocal = maxlocal = 0;
  loc2glob = NULL;
  vec_local = NULL;
  array_local = NULL;

  // zero accumulators one time if ave = RUNNING

  if (ave == RUNNING) {
    if (nvalues == 1)
      for (int i = 0; i < nslocal; i++)
	accvec[i] = 0.0;
    else {
      int m;
      for (int i = 0; i < nslocal; i++)
	for (m = 0; m < nvalues; m++)
	  accarray[i][m] = 0.0;
    }
  }

  // zero vector/array since dump may access it on timestep 0
  // zero vector/array since a variable may access it before first run

  if (nvalues == 1) {
    for (int i = 0; i < nslocal; i++)
      vector_surf[i] = 0.0;
  } else {
    int m;
    for (int i = 0; i < nslocal; i++)
      for (m = 0; m < nvalues; m++)
	array_surf[i][m] = 0.0;
  }

  // local storage of surf element masks

  int *mysurfs = surf->mysurfs;

  if (domain->dimension == 2) {
    Surf::Line *lines = surf->lines;
    for (int i = 0; i < nslocal; i++)
      masks[i] = lines[mysurfs[i]].mask;
  } else {
    Surf::Tri *tris = surf->tris;
    for (int i = 0; i < nslocal; i++)
      masks[i] = tris[mysurfs[i]].mask;
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

FixAveSurf::~FixAveSurf()
{
  delete [] which;
  delete [] argindex;
  delete [] value2index;
  for (int i = 0; i < nvalues; i++) delete [] ids[i];
  delete [] ids;

  memory->destroy(buflocal);
  memory->destroy(masks);

  if (nvalues == 1) memory->destroy(vector_surf);
  else memory->destroy(array_surf);
  if (ave == RUNNING) {
    if (nvalues == 1) memory->destroy(accvec);
    else memory->destroy(accarray);
  }

  memory->destroy(glob2loc);
  memory->destroy(loc2glob);
  memory->destroy(vec_local);
  memory->destroy(array_local);
}

/* ---------------------------------------------------------------------- */

int FixAveSurf::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAveSurf::init()
{
  if ((domain->dimension == 2 && nsurf != surf->nline) || 
      (domain->dimension == 3 && nsurf != surf->ntri)) 
    error->all(FLERR,"Number of surface elements changed in dump surf");

  // set indices and check validity of all computes,fixes,variables

  for (int m = 0; m < nvalues; m++) {
    if (which[m] == COMPUTE) {
      int icompute = modify->find_compute(ids[m]);
      if (icompute < 0)
	error->all(FLERR,"Compute ID for fix ave/surf does not exist");
      value2index[m] = icompute;
      
    } else if (which[m] == FIX) {
      int ifix = modify->find_fix(ids[m]);
      if (ifix < 0) 
	error->all(FLERR,"Fix ID for fix ave/surf does not exist");
      value2index[m] = ifix;

    } else if (which[m] == VARIABLE) {
      int ivariable = input->variable->find(ids[m]);
      if (ivariable < 0) 
	error->all(FLERR,"Variable name for fix ave/surf does not exist");
      value2index[m] = ivariable;

    } else value2index[m] = -1;
  }
}

/* ----------------------------------------------------------------------
   only does something if nvalid = current timestep
------------------------------------------------------------------------- */

void FixAveSurf::setup()
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveSurf::end_of_step()
{
  int i,j,k,m,n,isurf,ilocal;
  double *vec;

  // skip if not step which requires doing something

  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;

  // zero accumulators and norms if ave = ONE and first sample

  if (ave == ONE && irepeat == 0) {
    if (nvalues == 1)
      for (i = 0; i < nslocal; i++)
	accvec[i] = 0.0;
    else
      for (i = 0; i < nslocal; i++)
	for (m = 0; m < nvalues; m++)
	  accarray[i][m] = 0.0;
  }

  // reset all set glob2loc values to -1 and nlocal to 0 if first sample

  if (irepeat == 0) {
    for (i = 0; i < nlocal; i++) glob2loc[loc2glob[i]] = -1;
    nlocal = 0;
  }

  // accumulate results of computes,fixes,variables
  // compute/fix/variable may invoke computes so wrap with clear/add
  // NOTE: need to add logic for fixes and variables if enable them

  modify->clearstep_compute();

  for (m = 0; m < nvalues; m++) {
    n = value2index[m];
    j = argindex[m];

    // invoke compute if not previously invoked
    
    if (which[m] == COMPUTE) {
      Compute *compute = modify->compute[n];
      if (!(compute->invoked_flag & INVOKED_PER_SURF)) {
        compute->compute_per_surf();
        compute->invoked_flag |= INVOKED_PER_SURF;
      }
      int *loc2glob_compute;
      int nlocal_compute = compute->surfinfo(loc2glob_compute);
      
      if (j == 0) {
        double *vector = compute->vector_surf_tally;
        if (nvalues == 1) {
          for (i = 0; i < nlocal_compute; i++) {
            isurf = loc2glob_compute[i];
            ilocal = glob2loc[isurf];
            if (ilocal < 0) {
              if (nlocal == maxlocal) grow_local();
              ilocal = nlocal++;
              loc2glob[ilocal] = isurf;
              glob2loc[isurf] = ilocal;
              vec_local[ilocal] = 0.0;
            }
            vec_local[ilocal] += vector[i];
          }
        } else {
          for (i = 0; i < nlocal_compute; i++) {
            isurf = loc2glob_compute[i];
            ilocal = glob2loc[isurf];
            if (ilocal < 0) {
              if (nlocal == maxlocal) grow_local();
              ilocal = nlocal++;
              loc2glob[ilocal] = isurf;
              glob2loc[isurf] = ilocal;
              vec = array_local[ilocal];
              for (k = 0; k < nvalues; k++) vec[k] = 0.0;
            }
            array_local[ilocal][m] += vector[i];
          }
        }
      } else {
        int jm1 = j - 1;
        double **array = compute->array_surf_tally;
        if (nvalues == 1) {
          for (i = 0; i < nlocal_compute; i++) {
            isurf = loc2glob_compute[i];
            ilocal = glob2loc[isurf];
            if (ilocal < 0) {
              if (nlocal == maxlocal) grow_local();
              ilocal = nlocal++;
              loc2glob[ilocal] = isurf;
              glob2loc[isurf] = ilocal;
              vec_local[ilocal] = 0.0;
            }
            vec_local[ilocal] += array[i][jm1];
          }
        } else {
          for (i = 0; i < nlocal_compute; i++) {
            isurf = loc2glob_compute[i];
            ilocal = glob2loc[isurf];
            if (ilocal < 0) {
              if (nlocal == maxlocal) grow_local();
              ilocal = nlocal++;
              loc2glob[ilocal] = isurf;
              glob2loc[isurf] = ilocal;
              vec = array_local[ilocal];
              for (k = 0; k < nvalues; k++) vec[k] = 0.0;
            }
            array_local[ilocal][m] += array[i][jm1];
          }
        }
      }
      
    // access fix fields, guaranteed to be ready
    // check group mask in case other fix uses a different surf group

    } else if (which[m] == FIX) {
      if (j == 0) {
        double *fix_vector = modify->fix[n]->vector_surf;
        if (nvalues == 1) {
          for (i = 0; i < nslocal; i++)
            if (masks[i] & groupbit) accvec[i] += fix_vector[i];
        } else {
          for (i = 0; i < nslocal; i++)
            if (masks[i] & groupbit) accarray[i][m] += fix_vector[i];
        }
      } else {
        int jm1 = j - 1;
        double **fix_array = modify->fix[n]->array_surf;
        if (nvalues == 1) {
          for (i = 0; i < nslocal; i++)
            if (masks[i] & groupbit) accvec[i] += fix_array[i][jm1];
        } else {
          for (i = 0; i < nslocal; i++)
            if (masks[i] & groupbit) accarray[i][m] += fix_array[i][jm1];
        }
      }

    // evaluete surf-style variable
      
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
  nvalid = ntimestep+per_surf_freq - (nrepeat-1)*nevery;
  modify->addstep_compute(nvalid);

  // for values from computes, merge local values across all procs
  // returns surf element values I own in buflocal
  // then sum buflocal into accvec and accarray
  //   do not sum for surf elements not in surf group
  //   so those values stay 0.0 even if compute's surf group is different
  // NOTE: could test if all values are from COMPUTE and do collate_array()

  for (m = 0; m < nvalues; m++) {
    if (which[m] != COMPUTE) continue;
    if (nvalues == 1) {
      surf->collate_vector(nlocal,loc2glob,vec_local,1,buflocal);
      for (i = 0; i < nslocal; i++)
        if (masks[i] & groupbit) accvec[i] += buflocal[i];
    } else {
      if (nlocal)
        surf->collate_vector(nlocal,loc2glob,&array_local[0][m],
                             nvalues,buflocal);
      else
        surf->collate_vector(nlocal,loc2glob,NULL,nvalues,buflocal);
      for (i = 0; i < nslocal; i++)
        if (masks[i] & groupbit) accarray[i][m] += buflocal[i];
    }
  }

  // normalize the accumulators for output on Nfreq timestep
  // normindex < 0, just normalize by # of samples
  // normindex >= 0, normalize by accumulated norm vector

  if (ave == ONE) {
    if (nvalues == 1)
      for (i = 0; i < nslocal; i++) vector_surf[i] /= nsample;
    else
      for (m = 0; m < nvalues; m++)
        for (i = 0; i < nslocal; i++)
          array_surf[i][m] /= nsample;

  } else {
    if (nvalues == 1)
      for (i = 0; i < nslocal; i++) vector_surf[i] = accvec[i]/nsample;
    else
      for (m = 0; m < nvalues; m++)
        for (i = 0; i < nslocal; i++)
          array_surf[i][m] = accarray[i][m]/nsample;
  }

  // reset nsample if ave = ONE

  if (ave == ONE) nsample = 0;
}

/* ----------------------------------------------------------------------
   parse optional args
------------------------------------------------------------------------- */

void FixAveSurf::options(int iarg, int narg, char **arg)
{
  // option defaults

  ave = ONE;

  // optional args

  while (iarg < narg) {
    if (strcmp(arg[iarg],"ave") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/surf command");
      if (strcmp(arg[iarg+1],"one") == 0) ave = ONE;
      else if (strcmp(arg[iarg+1],"running") == 0) ave = RUNNING;
      else error->all(FLERR,"Illegal fix ave/surf command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix ave/surf command");
  }
}

/* ---------------------------------------------------------------------- */

void FixAveSurf::grow_local()
{
  maxlocal += DELTA;
  memory->grow(loc2glob,maxlocal,"ave/surf:loc2glob");
  if (nvalues == 1)
    memory->grow(vec_local,maxlocal,"ave/surf:vec_local");
  else
    memory->grow(array_local,maxlocal,nvalues,"ave/surf:array_local");
}

/* ----------------------------------------------------------------------
   memory usage of accumulators
------------------------------------------------------------------------- */

double FixAveSurf::memory_usage()
{
  double bytes = 0.0;
  bytes += nslocal*nvalues * sizeof(double);
  if (ave == RUNNING) bytes += nslocal*nvalues * sizeof(double);
  //bytes += nslocal*nnorm * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   calculate nvalid = next step on which end_of_step does something
   can be this timestep if multiple of nfreq and nrepeat = 1
   else backup from next multiple of nfreq
------------------------------------------------------------------------- */

bigint FixAveSurf::nextvalid()
{
  bigint nvalid = (update->ntimestep/per_surf_freq)*per_surf_freq + 
    per_surf_freq;
  if (nvalid-per_surf_freq == update->ntimestep && nrepeat == 1)
    nvalid = update->ntimestep;
  else
    nvalid -= (nrepeat-1)*nevery;
  if (nvalid < update->ntimestep) nvalid += per_surf_freq;
  return nvalid;
}

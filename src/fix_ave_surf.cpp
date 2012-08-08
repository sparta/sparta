/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2012) Sandia Corporation.  Under the terms of Contract
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
  if (narg < 5) error->all(FLERR,"Illegal fix ave/surf command");

  nevery = atoi(arg[2]);
  nrepeat = atoi(arg[3]);
  per_surf_freq = atoi(arg[4]);

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
	  error->all(FLERR,"Illegal fix ave/surf command");
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
            if (modify->compute[icompute]->per_surf_flag)
              ndup = modify->compute[icompute]->size_per_surf_cols;
        }
        if (which[nvalues] == FIX) {
          int ifix = modify->find_fix(ids[nvalues]);
          if (ifix >= 0)
            if (modify->fix[ifix]->per_surf_flag)
              ndup = modify->fix[ifix]->size_per_surf_cols;
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

  // allocate accumulators and norm vectors
  // if ave = RUNNING, allocate extra set of accvec/accarray

  if (domain->dimension == 2) nsurf = surf->nline;
  else nsurf = surf->ntri;
  nslocal = surf->nlocal;

  if (nvalues == 1) memory->create(vector_surf,nslocal,"ave/surf:vector_surf");
  else memory->create(array_surf,nslocal,nvalues,"ave/surf:array_surf");

  if (ave == RUNNING) {
    if (nvalues == 1) memory->create(accvec,nslocal,"ave/surf:accvec");
    else memory->create(accarray,nslocal,nvalues,"ave/surf:accarray");
  } else {
    if (nvalues == 1) accvec = vector_surf;
    else accarray = array_surf;
  }

  nnorm = 0;
  normacc = new int[nvalues];
  normindex = new int[nvalues];
  norms = new double*[nvalues];
  cfv_norms = new double*[nvalues];

  // zero mcount,accvec,accarray one time if ave = RUNNING

  int i,m;

  if (ave == RUNNING) {
    if (nvalues == 1)
      for (i = 0; i < nslocal; i++)
	accvec[i] = 0.0;
    else
      for (i = 0; i < nslocal; i++)
	for (m = 0; m < nvalues; m++)
	  accarray[i][m] = 0.0;
  }

  // zero vector/array since dump may access it on timestep 0
  // zero vector/array since a variable may access it before first run

  if (nvalues == 0) {
    for (i = 0; i < nslocal; i++)
      vector_surf[i] = 0.0;
  } else {
    for (i = 0; i < nslocal; i++)
      for (m = 0; m < nvalues; m++)
	array_surf[i][m] = 0.0;
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
  memory->destroy(which);
  memory->destroy(argindex);
  memory->destroy(value2index);
  for (int i = 0; i < nvalues; i++) delete [] ids[i];
  memory->sfree(ids);

  if (nvalues == 1) memory->destroy(vector_surf);
  else memory->destroy(array_surf);
  if (ave == RUNNING) {
    if (nvalues == 1) memory->destroy(accvec);
    else memory->destroy(accarray);
  }

  delete [] normacc;
  delete [] normindex;
  for (int i = 0; i < nnorm; i++) memory->sfree(norms[i]);
  delete [] norms;
  delete [] cfv_norms;
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
  // ont-time setup of norm pointers and nnorm
  // must do here after computes have been initialized
  // only store unique norms by checking if returned ptr matches previous ptr
  // zero norm vectors one time if ave = RUNNING
  // NOTE: need to add logic for fixes and variables if enable them

  if (nnorm == 0) {
    double *ptr;
    int j,n,iptr;

    for (int m = 0; m < nvalues; m++) {
      n = value2index[m];
      j = argindex[m];
      
      if (which[m] == COMPUTE) {
        Compute *compute = modify->compute[n];
        ptr = compute->normptr(j-1);
        if (ptr == NULL) {
          normacc[m] = 0;
          normindex[m] = -1;
        } else {
          for (iptr = 0; iptr < nnorm; iptr++)
            if (ptr == cfv_norms[iptr]) break;
          if (iptr < nnorm) {
            normacc[m] = 0;
            normindex[m] = iptr;
          } else {
            normacc[m] = 1;
            normindex[m] = nnorm;
            cfv_norms[nnorm] = ptr;
            memory->create(norms[nnorm],nslocal,"ave/surf:norms");
            if (ave == RUNNING) {
              double *norm = norms[nnorm];
              for (int i = 0; i < nslocal; i++) norm[i] = 0.0;
            }
            nnorm++;
          }
        }
      }
    }
  }

  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveSurf::end_of_step()
{
  int i,j,m,n,isp,icol;
  double *norm,*cfv_norm;

  // skip if not step which requires doing something

  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;

  // zero accumulators and norms if ave = ONE and first sample
  // mpivecone, mpiearrayone are global nsurf in length
  // NOTE: change this when remove Allreduce

  if (ave == ONE && irepeat == 0) {
    if (nvalues == 1)
      for (i = 0; i < nsurf; i++)
	mpivecone[i] = 0.0;
    else
      for (i = 0; i < nsurf; i++)
	for (m = 0; m < nvalues; m++)
	  mpiarrayone[i][m] = 0.0;
    for (int m = 0; m < nnorm; m++) {
      norm = norms[m];
      for (i = 0; i < nslocal; i++) norm[i] = 0.0;
    }
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
      
      if (j == 0) {
        //int nlocal = compute->nlocal;
        //int *loc2glob = compute->loc2glob;
        int nlocal = 0;
        int *loc2glob = NULL;
        double *compute_vector = compute->vector_surf;
        if (nvalues == 1) {
          for (i = 0; i < nlocal; i++)
            mpivecone[loc2glob[i]] += compute_vector[i];
        } else {
          for (i = 0; i < nlocal; i++)
            mpiarrayone[loc2glob[i]][m] += compute_vector[i];
        }
      } else {
        int jm1 = j - 1;
        //int nlocal = compute->nlocal;
        //int *loc2glob = compute->loc2glob;
        int nlocal = 0;
        int *loc2glob = NULL;
        double **compute_array = compute->array_surf;
        if (nvalues == 1) {
          for (i = 0; i < nlocal; i++)
            mpivecone[loc2glob[i]] += compute_array[i][jm1];
        } else {
          for (i = 0; i < nlocal; i++)
            mpiarrayone[loc2glob[i]][m] += compute_array[i][jm1];
        }
      }
      
      // access fix fields, guaranteed to be ready
      
    } else if (which[m] == FIX) {
      if (j == 0) {
        double *fix_vector = modify->fix[n]->vector_surf;
        if (nvalues == 1) {
          for (i = 0; i < nslocal; i++)
            accvec[i] += fix_vector[i];
        } else {
          for (i = 0; i < nslocal; i++)
            accarray[i][m] += fix_vector[i];
        }
      } else {
        int jm1 = j - 1;
        double **fix_array = modify->fix[n]->array_surf;
        if (nvalues == 1) {
          for (i = 0; i < nslocal; i++)
            accvec[i] += fix_array[i][jm1];
        } else {
          for (i = 0; i < nslocal; i++)
            accarray[i][m] += fix_array[i][jm1];
        }
      }
      
      // evaluete surf-style variable
      
    } else if (which[m] == VARIABLE) {
    }

    // accumulate norm values if necessary
    // only done once per unique norm

    if (normacc[m]) {
      norm = norms[normindex[m]];
      cfv_norm = cfv_norms[normindex[m]];
      for (i = 0; i < nslocal; i++) norm[i] += cfv_norm[i];
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

  // sum mpi arrays across all procs via Allreduce
  // copy summed value into owned-surf vec, array
  // NOTE: remove this when use remove Allreduce

  if (nvalues == 1)
    MPI_Allreduce(mpivecone,mpivec,nsurf,MPI_DOUBLE,MPI_SUM,world);
  else
    MPI_Allreduce(&mpiarrayone[0][0],&mpiarray[0][0],nsurf*nvalues,
                  MPI_DOUBLE,MPI_SUM,world);

  int *mysurfs = surf->mysurfs;

  if (nvalues == 1)
    for (i = 0; i < nslocal; i++) vector_surf[i] = mpivec[mysurfs[i]];
  else
    for (i = 0; i < nslocal; i++) {
      m = mysurfs[i];
      for (j = 0; j < nvalues; j++)
        array_surf[i][j] = mpiarray[m][j];
    }

  // normalize the accumulators for output on Nfreq timestep
  // normindex = 0, just normalize by # of samples
  // normindex = 1, normalize by accumulated norm vector
  // reset nsample if ave = ONE

  if (ave == ONE) {
    if (nvalues == 1) {
      n = value2index[0];
      j = argindex[0];
      if (normindex[0] < 0) {
        for (i = 0; i < nslocal; i++) vector_surf[i] /= nsample;
      } else {
        norm = norms[normindex[0]];
        for (i = 0; i < nslocal; i++)
          if (norm[i] > 0.0) vector_surf[i] /= norm[i];
      }
    } else {
      for (m = 0; m < nvalues; m++) {
        n = value2index[m];
        j = argindex[m];

        if (normindex[m] < 0) {
          for (i = 0; i < nslocal; i++) array_surf[i][m] /= nsample;
        } else {
          norm = norms[normindex[m]];
          for (i = 0; i < nslocal; i++)
            if (norm[i] > 0.0) array_surf[i][m] /= norm[i];
        }
      }
    }

  } else {
    if (nvalues == 1) {
      n = value2index[0];
      j = argindex[0];
      if (normindex[0] < 0) {
        for (i = 0; i < nslocal; i++) vector_surf[i] = accvec[i]/nsample;
      } else {
        norm = norms[normindex[0]];
        for (i = 0; i < nslocal; i++) vector_surf[i] = accvec[i]/norm[i];
      }
    } else {
      for (m = 0; m < nvalues; m++) {
        n = value2index[m];
        j = argindex[m];
        if (normindex[m] < 0) {
          for (i = 0; i < nslocal; i++)
            array_surf[i][m] = accarray[i][m]/nsample;
        } else {
          norm = norms[normindex[m]];
          for (i = 0; i < nslocal; i++) 
            array_surf[i][m] = accarray[i][m]/norm[i];
        }
      }
    }
  }

  if (ave == ONE) nsample = 0;
}

/* ----------------------------------------------------------------------
   parse optional args
------------------------------------------------------------------------- */

void FixAveSurf::options(int narg, char **arg)
{
  // option defaults

  ave = ONE;

  // optional args

  int iarg = 0;
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

/* ----------------------------------------------------------------------
   grow vectors for each input value
------------------------------------------------------------------------- */

void FixAveSurf::grow()
{
  maxvalues += DELTA;
  memory->grow(which,maxvalues,"ave/surf:which");
  memory->grow(argindex,maxvalues,"ave/surf:argindex");
  memory->grow(value2index,maxvalues,"ave/surf:value2index");
  ids = (char **) memory->srealloc(ids,maxvalues*sizeof(char *),"ave/surf:ids");
}

/* ----------------------------------------------------------------------
   memory usage of accumulators
------------------------------------------------------------------------- */

double FixAveSurf::memory_usage()
{
  double bytes = 0.0;
  bytes += nslocal*nvalues * sizeof(double);
  if (ave == RUNNING) bytes += nslocal*nvalues * sizeof(double);
  bytes += nslocal*nnorm * sizeof(double);
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

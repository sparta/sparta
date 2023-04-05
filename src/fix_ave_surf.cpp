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

#include "spatype.h"
#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "fix_ave_surf.h"
#include "surf.h"
#include "particle.h"
#include "update.h"
#include "domain.h"
#include "comm.h"
#include "modify.h"
#include "compute.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "sparta_masks.h"

using namespace SPARTA_NS;

enum{COMPUTE,FIX,VARIABLE};
enum{ONE,RUNNING};

#define INVOKED_PER_SURF 32
#define DELTA 1024;

/* ---------------------------------------------------------------------- */

FixAveSurf::FixAveSurf(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal fix ave/surf command");

  if (surf->implicit)
    error->all(FLERR,"Cannot use fix ave/surf with implicit surfs");

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

  // if any input is a compute, all must be

  int cflag = 0;
  for (int i = 0; i < nvalues; i++)
    if (which[i] == COMPUTE) cflag++;

  if (cflag && cflag != nvalues)
    error->all(FLERR,"Fix ave/surf inputs must be all computes or no computes");

  // this fix produces either a per-surf vector or array

  per_surf_flag = 1;
  if (nvalues == 1) size_per_surf_cols = 0;
  else size_per_surf_cols = nvalues;

  // allocate accumulators for owned surfaces
  // if ave = RUNNING, allocate extra set of accvec/accarray

  nown = surf->nown;
  memory->create(masks,nown,"ave/surf:masks");

  if (cflag) {
    bufvec = NULL;
    bufarray = NULL;
    if (nvalues == 1) memory->create(bufvec,nown,"ave/surf:bufvec");
    else memory->create(bufarray,nown,nvalues,"ave/surf:bufarray");
  }

  if (nvalues == 1) memory->create(vector_surf,nown,"ave/surf:vector_surf");
  else memory->create(array_surf,nown,nvalues,"ave/surf:array_surf");

  if (ave == RUNNING) {
    if (nvalues == 1) memory->create(accvec,nown,"ave/surf:accvec");
    else memory->create(accarray,nown,nvalues,"ave/surf:accarray");
  } else {
    if (nvalues == 1) accvec = vector_surf;
    else accarray = array_surf;
  }

  // tally accumulators

  ntally = maxtally = 0;
  tally2surf = NULL;
  vec_tally = NULL;
  array_tally = NULL;

  // zero accumulators one time if ave = RUNNING

  if (ave == RUNNING) {
    if (nvalues == 1)
      for (int i = 0; i < nown; i++)
        accvec[i] = 0.0;
    else {
      int m;
      for (int i = 0; i < nown; i++)
        for (m = 0; m < nvalues; m++)
          accarray[i][m] = 0.0;
    }
  }

  // zero vector/array since dump may access it on timestep 0
  // zero vector/array since a variable may access it before first run

  if (nvalues == 1) {
    for (int i = 0; i < nown; i++)
      vector_surf[i] = 0.0;
  } else {
    int m;
    for (int i = 0; i < nown; i++)
      for (m = 0; m < nvalues; m++)
        array_surf[i][m] = 0.0;
  }

  // set surf element masks for owned surfs

  if (surf->distributed) {
    if (domain->dimension == 2) {
      Surf::Line *lines = surf->mylines;
      for (int i = 0; i < nown; i++)
        masks[i] = lines[i].mask;
    } else {
      Surf::Tri *tris = surf->mytris;
      for (int i = 0; i < nown; i++)
      masks[i] = tris[i].mask;
    }

  } else {
    int me = comm->me;
    int nprocs = comm->nprocs;
    int nsurf = surf->nsurf;
    int m = 0;
    if (domain->dimension == 2) {
      Surf::Line *lines = surf->lines;
      for (int i = me; i < nsurf; i += nprocs)
        masks[m++] = lines[i].mask;
    } else {
      Surf::Tri *tris = surf->tris;
      for (int i = me; i < nsurf; i += nprocs)
        masks[m++] = tris[i].mask;
    }
  }

  // nvalid = next step on which end_of_step does something
  // add nvalid to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  nsample = 0;
  irepeat = 0;
  nvalid = nextvalid();
  modify->addstep_compute_all(nvalid);

  // hash for mapping surfIDs to tally indices

  hash = new MyHash;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

FixAveSurf::~FixAveSurf()
{
  delete [] which;
  delete [] argindex;
  delete [] value2index;
  for (int i = 0; i < nvalues; i++) delete [] ids[i];
  delete [] ids;

  memory->destroy(bufvec);
  memory->destroy(bufarray);
  memory->destroy(masks);

  if (nvalues == 1) memory->destroy(vector_surf);
  else memory->destroy(array_surf);
  if (ave == RUNNING) {
    if (nvalues == 1) memory->destroy(accvec);
    else memory->destroy(accarray);
  }

  memory->destroy(tally2surf);
  memory->destroy(vec_tally);
  memory->destroy(array_tally);

  delete hash;
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
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveSurf::end_of_step()
{
  int i,j,k,m,n,isurf,itally;
  surfint surfID;
  double *vec;

  // skip if not step which requires doing something

  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;

  // zero accumulators if ave = ONE and first sample

  if (ave == ONE && irepeat == 0) {
    if (nvalues == 1)
      for (i = 0; i < nown; i++)
        accvec[i] = 0.0;
    else
      for (i = 0; i < nown; i++)
        for (m = 0; m < nvalues; m++)
          accarray[i][m] = 0.0;
  }

  // clear hash of tallied surf IDs if first sample

  if (irepeat == 0) {
    hash->clear();
    ntally = 0;
  }

  // accumulate results of computes,fixes,variables
  // compute/fix/variable may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  for (m = 0; m < nvalues; m++) {
    n = value2index[m];
    j = argindex[m];

    // access list of tallies from compute, add to my list

    if (which[m] == COMPUTE) {
      Compute *compute = modify->compute[n];
      if (!(compute->invoked_flag & INVOKED_PER_SURF)) {
        compute->compute_per_surf();
        compute->invoked_flag |= INVOKED_PER_SURF;
      }
      surfint *tally2surf_compute;
      int ntally_compute = compute->tallyinfo(tally2surf_compute);

      if (j == 0) {
        double *vector = compute->vector_surf_tally;
        if (nvalues == 1) {
          for (i = 0; i < ntally_compute; i++) {
            surfID = tally2surf_compute[i];
            if (hash->find(surfID) != hash->end()) itally = (*hash)[surfID];
            else {
              if (ntally == maxtally) grow_tally();
              itally = ntally;
              (*hash)[surfID] = itally;
              tally2surf[itally] = surfID;
              vec_tally[itally] = 0.0;
              ntally++;
            }
            vec_tally[itally] += vector[i];
          }
        } else {
          for (i = 0; i < ntally_compute; i++) {
            surfID = tally2surf_compute[i];
            if (hash->find(surfID) != hash->end()) itally = (*hash)[surfID];
            else {
              if (ntally == maxtally) grow_tally();
              itally = ntally;
              (*hash)[surfID] = itally;
              tally2surf[itally] = surfID;
              vec = array_tally[itally];
              for (k = 0; k < nvalues; k++) vec[k] = 0.0;
              ntally++;
            }
            array_tally[itally][m] += vector[i];
          }
        }
      } else {
        int jm1 = j - 1;
        double **array = compute->array_surf_tally;
        if (nvalues == 1) {
          for (i = 0; i < ntally_compute; i++) {
            surfID = tally2surf_compute[i];
            if (hash->find(surfID) != hash->end()) itally = (*hash)[surfID];
            else {
              if (ntally == maxtally) grow_tally();
              itally = ntally;
              (*hash)[surfID] = itally;
              tally2surf[itally] = surfID;
              vec_tally[itally] = 0.0;
              ntally++;
            }
            vec_tally[itally] += array[i][jm1];
          }
        } else {
          for (i = 0; i < ntally_compute; i++) {
            surfID = tally2surf_compute[i];
            if (hash->find(surfID) != hash->end()) itally = (*hash)[surfID];
            else {
              if (ntally == maxtally) grow_tally();
              itally = ntally;
              (*hash)[surfID] = itally;
              tally2surf[itally] = surfID;
              vec = array_tally[itally];
              for (k = 0; k < nvalues; k++) vec[k] = 0.0;
              ntally++;
            }
            array_tally[itally][m] += array[i][jm1];
          }
        }
      }

    // access fix fields, guaranteed to be ready

    } else if (which[m] == FIX) {
      if (j == 0) {
        double *fix_vector = modify->fix[n]->vector_surf;
        if (nvalues == 1)
          for (i = 0; i < nown; i++) accvec[i] += fix_vector[i];
        else
          for (i = 0; i < nown; i++) accarray[i][m] += fix_vector[i];
      } else {
        int jm1 = j - 1;
        double **fix_array = modify->fix[n]->array_surf;
        if (nvalues == 1)
          for (i = 0; i < nown; i++) accvec[i] += fix_array[i][jm1];
        else
          for (i = 0; i < nown; i++) accarray[i][m] += fix_array[i][jm1];
      }

    // evaluate surf-style variable

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

  // invoke surf->collate() on tallies this fix stores for multiple steps
  // this merges tallies to owned surfs
  // NOTE: this should only be done if source is a COMPUTE ?

  if (nvalues == 1) {
    surf->collate_vector(ntally,tally2surf,vec_tally,1,bufvec);
    for (i = 0; i < nown; i++) accvec[i] += bufvec[i];
  } else {
    surf->collate_array(ntally,nvalues,tally2surf,array_tally,bufarray);
    for (i = 0; i < nown; i++)
      for (m = 0; m < nvalues; m++)
        accarray[i][m] += bufarray[i][m];
  }

  // normalize the accumulators for output on Nfreq timestep
  // everything is just normalized by # of samples

  if (ave == ONE) {
    if (nvalues == 1)
      for (i = 0; i < nown; i++) vector_surf[i] /= nsample;
    else
      for (i = 0; i < nown; i++)
        for (m = 0; m < nvalues; m++)
          array_surf[i][m] /= nsample;
  } else {
    if (nvalues == 1)
      for (i = 0; i < nown; i++)
        vector_surf[i] = accvec[i]/nsample;
    else
      for (i = 0; i < nown; i++)
        for (m = 0; m < nvalues; m++)
          array_surf[i][m] = accarray[i][m]/nsample;
  }

  // set values for surfs not in group to zero

  if (groupbit != 1) {
    if (nvalues == 1) {
      for (i = 0; i < nown; i++)
        if (!(masks[i] & groupbit)) vector_surf[i] = 0.0;
    } else {
      for (i = 0; i < nown; i++)
        if (!(masks[i] & groupbit))
          for (m = 0; m < nvalues; m++) array_surf[i][m] = 0.0;
    }
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

void FixAveSurf::grow_tally()
{
  maxtally += DELTA;
  memory->grow(tally2surf,maxtally,"ave/surf:tally2surf");
  if (nvalues == 1)
    memory->grow(vec_tally,maxtally,"ave/surf:vec_tally");
  else
    memory->grow(array_tally,maxtally,nvalues,"ave/surf:array_tally");
}

/* ----------------------------------------------------------------------
   memory usage of accumulators
------------------------------------------------------------------------- */

double FixAveSurf::memory_usage()
{
  double bytes = 0.0;
  bytes += nown*nvalues * sizeof(double);
  if (ave == RUNNING) bytes += nown*nvalues * sizeof(double);
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

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

enum{PERGRID,PERGRIDSURF};
enum{COMPUTE,FIX,VARIABLE};
enum{ONE,RUNNING};                // multiple files

#define INVOKED_PER_GRID 16
#define DELTAGRID 1024            // must be bigger than split cells per cell
#define DELTASURF 1024;

/* ---------------------------------------------------------------------- */

FixAveGrid::FixAveGrid(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal fix ave/grid command");

  int igroup = grid->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Could not find fix ave/grid group ID");
  groupbit = grid->bitmask[igroup];

  nevery = atoi(arg[3]);
  nrepeat = atoi(arg[4]);
  per_grid_freq = atoi(arg[5]);

  time_depend = 1;
  gridmigrate = 1;

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

  if (nvalues == 0) error->all(FLERR,"No values in fix ave/grid command");

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
  post_process = new int[nvalues];
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
        error->all(FLERR,"Illegal fix ave/grid command");
      argindex[i] = atoi(ptr+1);
      *ptr = '\0';
    } else argindex[i] = 0;

    n = strlen(suffix) + 1;
    ids[i] = new char[n];
    strcpy(ids[i],suffix);
    delete [] suffix;

    post_process[i] = 0;
    if (which[i] == COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix ave/grid does not exist");
      post_process[i] =
        modify->compute[icompute]->post_process_grid_flag;
    }
  }

  // if wildcard expansion occurred, free earg memory from expand_args()
  // wait to do this until after file comment lines are printed

  if (expand) {
    for (int i = 0; i < nvalues; i++) delete [] earg[i];
    memory->sfree(earg);
  }

  // setup and error check
  // for fix inputs, check that fix frequency is acceptable

  int flavor_pergrid = 0;
  int flavor_pergridsurf = 0;

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
                   "calculate per-grid vector");
      if (argindex[i] && modify->compute[icompute]->size_per_grid_cols == 0)
        error->all(FLERR,"Fix ave/grid compute does not "
                   "calculate per-grid array");
      if (argindex[i] &&
          argindex[i] > modify->compute[icompute]->size_per_grid_cols)
        error->all(FLERR,"Fix ave/grid compute array is accessed out-of-range");
      if (modify->compute[icompute]->post_process_isurf_grid_flag)
        flavor_pergridsurf = 1;
      else flavor_pergrid = 1;

    } else if (which[i] == FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for fix ave/grid does not exist");
      if (modify->fix[ifix]->per_grid_flag == 0)
        error->all(FLERR,"Fix ave/grid fix does not calculate per-grid values");
      if (argindex[i] == 0 && modify->fix[ifix]->size_per_grid_cols != 0)
        error->all(FLERR,
                   "Fix ave/grid fix does not calculate per-grid vector");
      if (argindex[i] && modify->fix[ifix]->size_per_grid_cols == 0)
        error->all(FLERR,
                   "Fix ave/grid fix does not calculate per-grid array");
      if (argindex[i] && argindex[i] > modify->fix[ifix]->size_per_grid_cols)
        error->all(FLERR,"Fix ave/grid fix array is accessed out-of-range");
      if (nevery % modify->fix[ifix]->per_grid_freq)
        error->all(FLERR,
                   "Fix for fix ave/grid not computed at compatible time");
      flavor_pergrid = 1;

    } else if (which[i] == VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for fix ave/grid does not exist");
      if (input->variable->grid_style(ivariable) == 0)
        error->all(FLERR,"Fix ave/grid variable is not grid-style variable");
      flavor_pergrid = 1;
    }
  }

  // require all inputs be flavor PERGRID or PERGRIDSURF

  if (flavor_pergrid && flavor_pergridsurf)
    error->all(FLERR,"Fix ave/grid cannot mix grid and grid/surf inputs");

  if (flavor_pergrid) flavor = PERGRID;
  if (flavor_pergridsurf) flavor = PERGRIDSURF;

  // this could be allowed but is dangerous if load-balancing
  // per-cellID tallies stored by this proc might grow enormous

  if (flavor == PERGRIDSURF && ave == RUNNING)
    error->all(FLERR,"Fix ave/grid for grid/surf inputs cannot use ave running");

  // this fix produces either a per-grid vector or array

  per_grid_flag = 1;
  if (nvalues == 1) size_per_grid_cols = 0;
  else size_per_grid_cols = nvalues;

  nglocal = maxgrid = grid->nlocal;

  // allocate per-grid cell data storage
  // zero vector/array grid in case used by dump or load balancer

  vector_grid = NULL;
  array_grid = NULL;

  if (nvalues == 1) {
    memory->create(vector_grid,nglocal,"ave/grid:vector_grid");
    for (int i = 0; i < nglocal; i++) vector_grid[i] = 0.0;
  } else {
    memory->create(array_grid,nglocal,nvalues,"ave/grid:array_grid");
    for (int i = 0; i < nglocal; i++)
      for (int m = 0; m < nvalues; m++) array_grid[i][m] = 0.0;
  }

  // nvalid = next step on which end_of_step does something
  // add nvalid to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  nsample = 0;
  irepeat = 0;
  nvalid = nextvalid();
  modify->addstep_compute_all(nvalid);

  // determine size of map/umap/uomap data structs and allocate them
  // tmax = max # of tally quantities for any value

  tmax = 1;
  for (int m = 0; m < nvalues; m++) {
    int n = -1;
    int j = argindex[m];
    if (which[m] != COMPUTE) continue;
    n = modify->find_compute(ids[m]);
    if (!modify->compute[n]->post_process_grid_flag) continue;
    double **array;
    int *cmap;
    int ncount = modify->compute[n]->query_tally_grid(j,array,cmap);
    tmax = MAX(tmax,ncount);
  }

  nmap = new int[nvalues];
  memory->create(map,nvalues,tmax,"ave/grid:map");
  numap = new int[nvalues];
  memory->create(umap,nvalues,tmax,"ave/grid:umap");
  memory->create(uomap,nvalues,tmax,"ave/grid:uomap");

  // setup nmap/map and numap/umap/uomap data structs for all values
  // ntotal = total # of unique tally quantities = columns in tally array

  ntotal = 0;
  for (int m = 0; m < nvalues; m++) {
    int n = -1;
    int j = argindex[m];

    int pflag = 0;
    if (which[m] == COMPUTE) {
      n = modify->find_compute(ids[m]);
      if (modify->compute[n]->post_process_grid_flag) pflag = 1;
    }

    // if not a compute that post-processes,
    // add single new tally to nmap/map and numap/umap

    if (!pflag) {
      nmap[m] = 1;
      map[m][0] = ntotal;
      numap[m] = 1;
      umap[m][0] = ntotal;
      ntotal++;

    // else add all compute tallies to nmap/map
    // and only unique compute tallies to numap/umap/uomap

    } else {
      double **array;
      int *cmap;
      int ncount = modify->compute[n]->query_tally_grid(j,array,cmap);
      nmap[m] = numap[m] = 0;
      for (int i = 0; i < ncount; i++) {

        // set ucol = -1 if first time this compute quantity is tallied
        // else set to tally column that already tallies it

        int col = -1;
        for (int mm = 0; mm <= m; mm++) {
          if (which[mm] != COMPUTE) continue;
          int nn = modify->find_compute(ids[mm]);
          if (!modify->compute[nn]->post_process_grid_flag) continue;
          if (n != nn) continue;  // not same compute
          for (int kk = 0; kk < numap[mm]; kk++)
            if (cmap[i] == uomap[mm][kk]) col = umap[mm][kk];
        }

        // if this quantity already tallied, just point nmap/map to it
        // else add to nmap/map and numap/umap/uomap

        if (col >= 0) {
          map[m][nmap[m]] = col;
          nmap[m]++;
        } else {
          map[m][nmap[m]] = ntotal;
          nmap[m]++;
          umap[m][numap[m]] = ntotal;
          uomap[m][numap[m]] = cmap[i];
          numap[m]++;
          ntotal++;
        }
      }
    }
  }

  // allocate tally array for flavor = PERGRID
  // zero in case used by ave = RUNNING or accessed for immediate output

  if (flavor == PERGRID) {
    memory->create(tally,nglocal,ntotal,"ave/grid:tally");
    for (int i = 0; i < nglocal; i++)
      for (int j = 0; j < ntotal; j++)
        tally[i][j] = 0.0;
  } else tally = NULL;

  // tally accumulators for flavor = PERGRIDSURF

  ntallyID = maxtallyID = 0;
  tally2cell = NULL;
  vec_tally = NULL;
  array_tally = NULL;

  // hash for mapping cellIDs to tally indices

  hash = new MyHash;
}

/* ---------------------------------------------------------------------- */

FixAveGrid::~FixAveGrid()
{
  if (copymode) return;

  delete [] which;
  delete [] argindex;
  delete [] value2index;
  delete [] post_process;
  for (int i = 0; i < nvalues; i++) delete [] ids[i];
  delete [] ids;

  delete [] nmap;
  memory->destroy(map);
  delete [] numap;
  memory->destroy(umap);
  memory->destroy(uomap);

  if (nvalues == 1) memory->destroy(vector_grid);
  else memory->destroy(array_grid);

  memory->destroy(tally);

  memory->destroy(tally2cell);
  memory->destroy(vec_tally);
  memory->destroy(array_tally);

  delete hash;
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

  // reallocate per-grid data if necessary

  nglocal = grid->nlocal;
  grow_percell(0);
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
  int i,j,k,m,n,itally;
  int ntally_col,kk;
  cellint cellID;
  int *itmp;
  double *vec;
  double **ctally;

  // skip if not step which requires doing something

  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;

  // zero grid tallies if ave = ONE and first sample
  // could do this with memset()

  if (flavor == PERGRID && ave == ONE && irepeat == 0) {
    for (i = 0; i < nglocal; i++)
      for (j = 0; j < ntotal; j++)
        tally[i][j] = 0.0;
  }

  // clear hash of cellID tallies if ave = ONE and first sample

  if (flavor == PERGRIDSURF && ave == ONE && irepeat == 0) {
    hash->clear();
    ntallyID = 0;
  }

  // accumulate results of computes,fixes,variables
  // compute/fix/variable may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  // PERGRID = all values are computes,fixes,variable which
  //           calculate per-grid values in standard manner

  if (flavor == PERGRID) {

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

        // accumulate one or more compute values to umap columns of tally array
        // if compute does not post-process, access its vec/array grid directly
        // else access uomap columns in its ctally array

        if (post_process[m]) {
          ntally_col = numap[m];
          compute->query_tally_grid(j,ctally,itmp);
          for (i = 0; i < nglocal; i++)
            for (itally = 0; itally < ntally_col; itally++) {
              k = umap[m][itally];
              kk = uomap[m][itally];
              tally[i][k] += ctally[i][kk];
            }
        } else {
          k = umap[m][0];
          if (j == 0) {
            double *compute_vector = compute->vector_grid;
            for (i = 0; i < nglocal; i++)
              tally[i][k] += compute_vector[i];
          } else {
            int jm1 = j - 1;
            double **compute_array = compute->array_grid;
            for (i = 0; i < nglocal; i++)
              tally[i][k] += compute_array[i][jm1];
          }
        }

      // access fix fields, guaranteed to be ready

      } else if (which[m] == FIX) {
        k = umap[m][0];
        if (j == 0) {
          double *fix_vector = modify->fix[n]->vector_grid;
          for (i = 0; i < nglocal; i++)
            tally[i][k] += fix_vector[i];
        } else {
          int jm1 = j - 1;
        double **fix_array = modify->fix[n]->array_grid;
        for (i = 0; i < nglocal; i++)
          tally[i][k] += fix_array[i][jm1];
        }

      // evaluate grid-style variable, sum values to Kth column of tally array

      } else if (which[m] == VARIABLE) {
        k = umap[m][0];
        input->variable->compute_grid(n,&tally[0][k],ntotal,1);
      }
    }

  // PERGRIDSURF = all values are computes which tally info on collisions
  //               with implicit surfs and store them as per-grid-cell tallies

  } else if (flavor == PERGRIDSURF) {

    for (m = 0; m < nvalues; m++) {
      n = value2index[m];
      j = argindex[m];

      // invoke compute if not previously invoked

      Compute *compute = modify->compute[n];
      if (!(compute->invoked_flag & INVOKED_PER_GRID)) {
        compute->compute_per_grid();
        compute->invoked_flag |= INVOKED_PER_GRID;
      }

      // similar logic to fix ave/surf and its per-surf tally accumulation

      surfint *tally2surf_compute;
      int ntallyID_compute = compute->tallyinfo(tally2surf_compute);
      cellint *tally2cell_compute = (cellint *) tally2surf_compute;

      if (j == 0) {
        double *vector = compute->vector_surf_tally;
        if (nvalues == 1) {
          for (i = 0; i < ntallyID_compute; i++) {
            cellID = tally2cell_compute[i];
            if (hash->find(cellID) != hash->end()) itally = (*hash)[cellID];
            else {
              if (ntallyID == maxtallyID) grow_tally();
              itally = ntallyID;
              (*hash)[cellID] = itally;
              tally2cell[itally] = cellID;
              vec_tally[itally] = 0.0;
              ntallyID++;
            }
            vec_tally[itally] += vector[i];
          }
        } else {
          for (i = 0; i < ntallyID_compute; i++) {
            cellID = tally2cell_compute[i];
            if (hash->find(cellID) != hash->end()) itally = (*hash)[cellID];
            else {
              if (ntallyID == maxtallyID) grow_tally();
              itally = ntallyID;
              (*hash)[cellID] = itally;
              tally2cell[itally] = cellID;
              vec = array_tally[itally];
              for (k = 0; k < nvalues; k++) vec[k] = 0.0;
              ntallyID++;
            }
            array_tally[itally][m] += vector[i];
          }
        }
      } else {
        int jm1 = j - 1;
        double **array = compute->array_surf_tally;
        if (nvalues == 1) {
          for (i = 0; i < ntallyID_compute; i++) {
            cellID = tally2cell_compute[i];
            if (hash->find(cellID) != hash->end()) itally = (*hash)[cellID];
            else {
              if (ntallyID == maxtallyID) grow_tally();
              itally = ntallyID;
              (*hash)[cellID] = itally;
              tally2cell[itally] = cellID;
              vec_tally[itally] = 0.0;
              ntallyID++;
            }
            vec_tally[itally] += array[i][jm1];
          }
        } else {
          for (i = 0; i < ntallyID_compute; i++) {
            cellID = tally2cell_compute[i];
            if (hash->find(cellID) != hash->end()) itally = (*hash)[cellID];
            else {
              if (ntallyID == maxtallyID) grow_tally();
              itally = ntallyID;
              (*hash)[cellID] = itally;
              tally2cell[itally] = cellID;
              vec = array_tally[itally];
              for (k = 0; k < nvalues; k++) vec[k] = 0.0;
              ntallyID++;
            }
            array_tally[itally][m] += array[i][jm1];
          }
        }
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

  // final PERGRID values for output
  // normalize the accumulators for output on Nfreq timestep
  // works for ave = ONE or RUNNING
  // if post_process flag set, compute performs normalization via pp_grid()
  // else just divide by nsample

  if (flavor == PERGRID) {
    if (nvalues == 1) {
      if (post_process[0]) {
        n = value2index[0];
        j = argindex[0];
        Compute *c = modify->compute[n];
        c->post_process_grid(j,nsample,tally,map[0],vector_grid,1);
      } else {
        k = map[0][0];
        for (i = 0; i < nglocal; i++) vector_grid[i] = tally[i][k] / nsample;
      }
    } else {
      for (m = 0; m < nvalues; m++) {
        if (post_process[m]) {
          n = value2index[m];
          j = argindex[m];
          Compute *c = modify->compute[n];
          if (array_grid) c->post_process_grid(j,nsample,tally,map[m],
                                               &array_grid[0][m],nvalues);
        } else {
          k = map[m][0];
          for (i = 0; i < nglocal; i++) array_grid[i][m] = tally[i][k] / nsample;
        }
      }
    }

  // final PERGRIDSURF values for output
  // invoke surf->collate() on cellID tallies this fix stores for multiple steps
  //   this merges tallies to owned grid cells
  // divide results by nsample
  // copy split cell values to their sub cells, used by dump grid

  } else if (flavor == PERGRIDSURF) {
    if (nvalues == 1)
      surf->collate_vector_implicit(ntallyID,tally2cell,vec_tally,vector_grid);
    else
      surf->collate_array_implicit(ntallyID,nvalues,tally2cell,
                                   array_tally,array_grid);

    Grid::ChildCell *cells = grid->cells;

    if (nvalues == 1) {
      for (int icell = 0; icell < nglocal; icell++)
        vector_grid[icell] /= nsample;
    } else {
      for (int icell = 0; icell < nglocal; icell++)
        for (m = 0; m < nvalues; m++)
          array_grid[icell][m] /= nsample;
    }

    Grid::SplitInfo *sinfo = grid->sinfo;

    int jcell,nsplit;
    int *csubs;

    for (int icell = 0; icell < nglocal; icell++) {
      if (cells[icell].nsplit <= 1) continue;
      nsplit = cells[icell].nsplit;
      csubs = sinfo[cells[icell].isplit].csubs;
      if (nvalues == 1) {
        for (int j = 0; j < nsplit; j++) {
          jcell = csubs[j];
          vector_grid[jcell] = vector_grid[icell];
        }
      } else {
        for (int j = 0; j < nsplit; j++) {
          jcell = csubs[j];
          memcpy(array_grid[jcell],array_grid[icell],nvalues*sizeof(double));
        }
      }
    }
  }

  // set values for grid cells not in group to zero

  if (groupbit != 1) {
    Grid::ChildInfo *cinfo = grid->cinfo;
    if (nvalues == 1) {
      for (i = 0; i < nglocal; i++)
        if (!(cinfo[i].mask & groupbit)) vector_grid[i] = 0.0;
    } else {
      for (i = 0; i < nglocal; i++)
        if (!(cinfo[i].mask & groupbit))
          for (m = 0; m < nvalues; m++) array_grid[i][m] = 0.0;
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

  if (flavor == PERGRID) {
    if (memflag) memcpy(ptr,tally[icell],ntotal*sizeof(double));
    ptr += ntotal*sizeof(double);
  }

  return ptr-buf;
}

/* ----------------------------------------------------------------------
   unpack icell values for per-cell arrays from buf
   if icell is a split cell, also unpack all sub cell values
   return byte count of amount unpacked
   NOTE: why packing/unpacking parent cell if a split cell?
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

  if (flavor == PERGRID) {
    memcpy(tally[icell],ptr,ntotal*sizeof(double));
    ptr += ntotal*sizeof(double);
  }

  return ptr-buf;
}

/* ----------------------------------------------------------------------
   copy per-cell info from Icell to Jcell
   called when a grid cell is removed from this processor's list
   caller checks that Icell != Jcell
------------------------------------------------------------------------- */

void FixAveGrid::copy_grid_one(int icell, int jcell)
{
  if (nvalues == 1) vector_grid[jcell] = vector_grid[icell];
  else memcpy(array_grid[jcell],array_grid[icell],nvalues*sizeof(double));
  if (flavor == PERGRID)
    memcpy(tally[jcell],tally[icell],ntotal*sizeof(double));
}

/* ----------------------------------------------------------------------
   add a grid cell
   called when a grid cell is added to this processor's list
   initialize values to 0.0
------------------------------------------------------------------------- */

void FixAveGrid::add_grid_one()
{
  grow_percell(1);

  if (nvalues == 1) vector_grid[nglocal] = 0.0;
  else
    for (int i = 0; i < nvalues; i++) array_grid[nglocal][i] = 0.0;

  if (flavor == PERGRID)
    for (int i = 0; i < ntotal; i++) tally[nglocal][i] = 0.0;

  nglocal++;
}

/* ----------------------------------------------------------------------
   reset final grid cell count after grid cell removals
------------------------------------------------------------------------- */

void FixAveGrid::reset_grid_count(int nlocal)
{
  nglocal = nlocal;
}

/* ----------------------------------------------------------------------
   parse optional args
------------------------------------------------------------------------- */

void FixAveGrid::options(int iarg, int narg, char **arg)
{
  // option defaults

  ave = ONE;

  // optional args

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
   insure per-cell arrays are allocated long enough for N new cells
------------------------------------------------------------------------- */

void FixAveGrid::grow_percell(int nnew)
{
  if (nglocal+nnew < maxgrid) return;

  int maxgridold = maxgrid;
  while (maxgrid < nglocal+nnew) maxgrid += DELTAGRID;

  if (nvalues == 1) memory->grow(vector_grid,maxgrid,"ave/grid:vector_grid");
  else memory->grow(array_grid,maxgrid,nvalues,"ave/grid:array_grid");
  if (flavor == PERGRID) memory->grow(tally,maxgrid,ntotal,"ave/grid:tally");

  if (nvalues == 1)
    for (int i = maxgridold; i < maxgrid; i++)
      vector_grid[i] = 0.0;
  else
    for (int i = maxgridold; i < maxgrid; i++)
      for (int j = 0; j < nvalues; j++)
        array_grid[i][j] = 0.0;
}

/* ---------------------------------------------------------------------- */

void FixAveGrid::grow_tally()
{
  maxtallyID += DELTASURF;
  memory->grow(tally2cell,maxtallyID,"ave/grid:tally2cell");
  if (nvalues == 1)
    memory->grow(vec_tally,maxtallyID,"ave/grid:vec_tally");
  else
    memory->grow(array_tally,maxtallyID,nvalues,"ave/grid:array_tally");
}

/* ----------------------------------------------------------------------
   memory usage of accumulators
------------------------------------------------------------------------- */

double FixAveGrid::memory_usage()
{
  double bytes = 0.0;
  bytes += maxgrid*nvalues * sizeof(double);    // vector or array grid
  if (flavor == PERGRID) bytes += ntotal*maxgrid * sizeof(double);
  if (flavor == PERGRIDSURF) {
    bytes += maxtallyID * sizeof(cellint);
    bytes += nvalues*maxtallyID * sizeof(double);
  }
  return bytes;
}

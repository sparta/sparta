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
#include "fix_ablate.h"
#include "grid.h"
#include "domain.h"
#include "modify.h"
#include "compute.h"
#include "memory.h"
#include "error.h"

// DEBUG
#include "comm.h"

using namespace SPARTA_NS;

enum{COMPUTE,FIX};

#define DELTAGRID 1024

// NOTES
// use grid group to only ablate some cells?
// should I store one value or 8 per cell
// should I store as bytes or doubles and only convert to double for output
// allow Nevery = 0 and no end_of_step mask
// which pack/unpack methods needed for grid migration, not adapt
// enforce existence of implicit surfs?
// how to restart and create correct array_grid?

/* ---------------------------------------------------------------------- */

FixAblate::FixAblate(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg != 4) error->all(FLERR,"Illegal fix ablate command");

  igroup = grid->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Could not find fix ablate group ID");
  groupbit = grid->bitmask[igroup];

  nevery = atoi(arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix ablate command");

  idsource = NULL;

  /*
  if ((strncmp(arg[4],"c_",2) == 0) || (strncmp(arg[4],"f_",2) == 0)) {
    if (arg[4][0] == 'c') which = COMPUTE;
    else if (arg[4][0] == 'f') which = FIX;

    int n = strlen(arg[4]);
    char *suffix = new char[n];
    strcpy(suffix,&arg[4][2]);

    char *ptr = strchr(suffix,'[');
    if (ptr) {
      if (suffix[strlen(suffix)-1] != ']')
        error->all(FLERR,"Illegal fix ablate command");
      argindex = atoi(ptr+1);
      *ptr = '\0';
    } else argindex = 0;

    n = strlen(suffix) + 1;
    idsource = new char[n];
    strcpy(idsource,suffix);
    delete [] suffix;

  } else error->all(FLERR,"Illegal fix ablate command");
  */

  // error check

  /*
  if (which == COMPUTE) {
    icompute = modify->find_compute(idsource);
    if (icompute < 0)
      error->all(FLERR,"Compute ID for fix ablate does not exist");
    if (modify->compute[icompute]->per_grid_flag == 0)
      error->all(FLERR,
                 "Fix ablate compute does not calculate per-grid values");
    if (argindex == 0 && 
        modify->compute[icompute]->size_per_grid_cols != 0)
      error->all(FLERR,"Fix ablate compute does not "
                 "calculate per-grid vector");
    if (argindex && modify->compute[icompute]->size_per_grid_cols == 0)
      error->all(FLERR,"Fix ablate compute does not "
                 "calculate per-grid array");
    if (argindex && argindex > modify->compute[icompute]->size_per_grid_cols)
      error->all(FLERR,"Fix ablate compute array is accessed out-of-range");

  } else if (which == FIX) {
    ifix = modify->find_fix(idsource);
    if (ifix < 0)
      error->all(FLERR,"Fix ID for fix ablate does not exist");
    if (modify->fix[ifix]->per_grid_flag == 0)
      error->all(FLERR,"Fix ablate fix does not calculate per-grid values");
    if (argindex == 0 && modify->fix[ifix]->size_per_grid_cols != 0)
      error->all(FLERR,
                 "Fix ablate fix does not calculate per-grid vector");
    if (argindex && modify->fix[ifix]->size_per_grid_cols == 0)
      error->all(FLERR,
                 "Fix ablate fix does not calculate per-grid array");
    if (argindex && argindex > modify->fix[ifix]->size_per_grid_cols)
      error->all(FLERR,"Fix ablate fix array is accessed out-of-range");
    if (nevery % modify->fix[ifix]->per_grid_freq)
      error->all(FLERR,
                 "Fix for fix ablate not computed at compatible time");
  }
  */

  // this fix produces a per-grid array

  per_grid_flag = 1;
  if (domain->dimension == 2) size_per_grid_cols = 4;
  else size_per_grid_cols = 8;
  per_grid_freq = 1;
  gridmigrate = 1;

  storeflag = 0;
  array_grid = NULL;
}

/* ---------------------------------------------------------------------- */

FixAblate::~FixAblate()
{
  delete [] idsource;
  memory->destroy(array_grid);
}

/* ---------------------------------------------------------------------- */

int FixAblate::setmask()
{
  int mask = 0;
  if (nevery) mask |= END_OF_STEP;
  return mask;
}

/* ----------------------------------------------------------------------
   store grid corner point values in array_grid
   called by ReadIsurf when corner point grid is read in
------------------------------------------------------------------------- */

void FixAblate::store_corners(int **corners)
{
  storeflag = 1;

  // allocate per-grid cell data storage
  // zero array in case used by dump or load balancer

  Grid::ChildCell *cells = grid->cells;
  Grid::SplitInfo *sinfo = grid->sinfo;
  nglocal = nglocalmax = grid->nlocal;
  ncols = size_per_grid_cols;

  memory->grow(array_grid,nglocal,ncols,"ablate:array_grid");
  for (int i = 0; i < nglocal; i++)
    for (int m = 0; m < ncols; m++) array_grid[i][m] = 0.0;

  // NOTE: assuming that corners is in same order as original cells
  //       new split cells were added to end

  for (int icell = 0; icell < nglocal; icell++) {
    if (cells[icell].nsplit > 0) {
      for (int m = 0; m < ncols; m++) 
        array_grid[icell][m] = corners[icell][m];
    } else {
      int splitcell = sinfo[cells[icell].isplit].icell;
      for (int m = 0; m < ncols; m++) 
        array_grid[icell][m] = array_grid[splitcell][m];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixAblate::init()
{
  if (!storeflag) 
    error->all(FLERR,"Fix ablate corner point values not yet stored");

  MPI_Comm_rank(world,&me);

  /*
  if (which == COMPUTE) {
    icompute = modify->find_compute(idsource);
    if (icompute < 0) 
      error->all(FLERR,"Compute ID for fix ablate does not exist");
  } else if (which == FIX) {
    ifix = modify->find_fix(idsource);
    if (ifix < 0) 
      error->all(FLERR,"Fix ID for fix ablate does not exist");
  }
  */
}

/* ---------------------------------------------------------------------- */

void FixAblate::end_of_step()
{
  /*

  // perform ablation
  // read from source
  // decrement corner point values
  // sync shared corner point values
  // invoke marching cubes
  // recompute all tris

  // count # of face datums to send
  // skip if face has no neigh cell within ISurf grid
  // skip if adjacent cell is owned by this proc

  int nsend = 0;
  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue
    if (cells[icell].nsurf == 0) continue;
    if (cells[icell].nsplit <= 0) continue;
    // NOTE: need to determine i,j,k of this cell in Isurf grid
    //       so no which faces do not matter
    if (dim == 2) {
      for (iface = 0; iface < 4; iface++) {
        // check 2 skip conditions
        nsend++;
      }
    } else {
      for (iface = 0; iface < 6; iface++) {
        // check 2 skip conditions
        nsend++;
      }
    }
  }


  // sync corner point values
  // ncomm = ilocal + iface + 2/4 corner points

  int ncomm;
  if (dim == 2) ncomm = 4;
  else ncomm = 6;

  if (nsend > maxsend) {
    maxsend = nsend;
    memory->destroy(proclist);
    memory->create(proclist,maxsend,"ablate:proclist");
    memory->destroy(sbuf);
    memory->create(sbuf,maxsend*ncomm,"ablate:sbuf");
  }

  // pack face datums to send
  // skip cells not in group or with no surfs or sub-cells
  // data = ilocal (of icell on other proc) + Ntotal values

  int m = 0;
  nsend = 0;
  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue
    if (cells[icell].nsurf == 0) continue;
    if (cells[icell].nsplit <= 0) continue;
    if (dim == 2) {
      for (iface = 0; iface < 4; iface++) {
        // check 2 skip conditions
        proclist[nsend] = cells[icell].proc;
        sbuf[m++] = cells[icell].ilocal;    // NOTE: need ubuf
        sbuf[m++] = iface;    // NOTE: needs to be inverse for receiving cell
        sbuf[m++] = delta[icell][f2c[iface][0]];   // NOTE: map face to 2 pts
        sbuf[m++] = delta[icell][f2c[iface][1]];
        nsend++;
      }
    } else {
      for (iface = 0; iface < 6; iface++) {
        // check 2 skip conditions
        proclist[nsend] = cells[icell].proc;
        sbuf[m++] = cells[icell].ilocal;    // NOTE: need ubuf
        sbuf[m++] = iface;    // NOTE: needs to be inverse for receiving cell
        sbuf[m++] = delta[icell][f2c[iface][0]];   // NOTE: map face to 4 pts
        sbuf[m++] = delta[icell][f2c[iface][1]];
        sbuf[m++] = delta[icell][f2c[iface][2]];
        sbuf[m++] = delta[icell][f2c[iface][3]];
        nsend++;
      }
    }
  }

  double *rbuf;
  int nrecv = comm->irregular_uniform_neighs(nsend,proclist,(char *) sbuf,
                                             ncomm*sizeof(double),
                                             (char **) &rbuf);
  
  // unpack received data and sum Ntotal values into array_grid

  m = 0;
  for (int i = 0; i < nrecv; i++) {
    ilocal = static_cast<int> (rbuf[m++]);   // NOTE: need ubuf
    iface = static_cast<int> (rbuf[m++]);
    if (dim == 2) {
      for (iface = 0; iface < 4; iface++) {
        delta[icell][f2c[iface][0]] += rbuf[m++];  // NOTE: map face to 2 pts
        delta[icell][f2c[iface][1]] += rbuf[m++];  // NOTE: is unorder the same
      }
    } else {
      for (iface = 0; iface < 4; iface++) {
        delta[icell][f2c[iface][0]] += rbuf[m++];  // NOTE: map face to 2 pts
        delta[icell][f2c[iface][1]] += rbuf[m++];  // NOTE: is unorder the same
        delta[icell][f2c[iface][2]] += rbuf[m++];
        delta[icell][f2c[iface][3]] += rbuf[m++];
      }
    }
  }

  // NOTE: also have to do internal face sharing
  // NOTE: how to insure no round-off errors
  //       so final shared delta values are exactly the same

  */
}

/* ----------------------------------------------------------------------
   pack icell values for per-cell arrays into buf
   if icell is a split cell, also pack all sub cell values 
   return byte count of amount packed
   if memflag, only return count, do not fill buf
------------------------------------------------------------------------- */

int FixAblate::pack_grid_one(int icell, char *buf, int memflag)
{
  char *ptr = buf;
  Grid::ChildCell *cells = grid->cells;
  Grid::SplitInfo *sinfo = grid->sinfo;

  if (memflag) memcpy(ptr,array_grid[icell],ncols*sizeof(double));
  ptr += ncols*sizeof(double);

  if (cells[icell].nsplit > 1) {
    int isplit = cells[icell].isplit;
    int nsplit = cells[icell].nsplit;
    for (int i = 0; i < nsplit; i++) {
      int jcell = sinfo[isplit].csubs[i];
      if (memflag) memcpy(ptr,array_grid[jcell],ncols*sizeof(double));
      ptr += ncols*sizeof(double);
    }
  }

  return ptr-buf;
}

/* ----------------------------------------------------------------------
   unpack icell values for per-cell array from buf
   return byte count of amount unpacked
------------------------------------------------------------------------- */

int FixAblate::unpack_grid_one(int icell, char *buf)
{
  char *ptr = buf;
  Grid::ChildCell *cells = grid->cells;
  Grid::SplitInfo *sinfo = grid->sinfo;

  grow_percell(1);
  memcpy(array_grid[icell],ptr,ncols*sizeof(double));
  ptr += ncols*sizeof(double);
  nglocal++;

 if (cells[icell].nsplit > 1) {
    int isplit = cells[icell].isplit;
    int nsplit = cells[icell].nsplit;
    grow_percell(nsplit);
    for (int i = 0; i < nsplit; i++) {
      int jcell = sinfo[isplit].csubs[i];
      memcpy(array_grid[jcell],ptr,ncols*sizeof(double));
      ptr += ncols*sizeof(double);
    }
    nglocal += nsplit;
  }

  return ptr-buf;
}

/* ----------------------------------------------------------------------
   compress per-cell arrays due to cells migrating to new procs
   criteria for keeping/discarding a cell is same as in Grid::compress()
   this keeps final ordering of per-cell arrays consistent with Grid class
------------------------------------------------------------------------- */

void FixAblate::compress_grid()
{
  Grid::ChildCell *cells = grid->cells;

  int ncurrent = nglocal;
  nglocal = 0;
  for (int icell = 0; icell < ncurrent; icell++) {
    if (cells[icell].nsplit >= 1) {
      if (cells[icell].proc != me) continue;
    } else {
      int isplit = cells[icell].isplit;
      if (cells[grid->sinfo[isplit].icell].proc != me) continue;
    }

    if (nglocal != icell)
      memcpy(array_grid[nglocal],array_grid[icell],ncols*sizeof(double));

    nglocal++;
  }
}

/* ----------------------------------------------------------------------
   memory usage of accumulators
------------------------------------------------------------------------- */

double FixAblate::memory_usage()
{
  double bytes = 0.0;
  bytes += nglocalmax*ncols * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   insure per-cell arrays are allocated long enough for N new cells
------------------------------------------------------------------------- */

void FixAblate::grow_percell(int nnew)
{
  if (nglocal+nnew < nglocalmax) return;
  nglocalmax += DELTAGRID;
  memory->grow(array_grid,nglocalmax,ncols,"ablate:array_grid");
}

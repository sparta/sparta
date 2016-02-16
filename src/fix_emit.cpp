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

#include "stdlib.h"
#include "string.h"
#include "fix_emit.h"
#include "update.h"
#include "domain.h"
#include "region.h"
#include "grid.h"
#include "comm.h"
#include "random_mars.h"
#include "random_park.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};           // same as Grid

#define DELTAGRID 1024
#define DELTACELL 1024

/* ---------------------------------------------------------------------- */

FixEmit::FixEmit(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  gridmigrate = 1;

  // RNG

  int me = comm->me;
  random = new RanPark(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,me,100);

  // local storage of emit data structures

  c2list = NULL;
  nglocal = nglocalmax = 0;
  clist = clistnum = clistfirst = NULL;
  nlist = nlistmax = 0;

  // counters common to all emit styles for output from fix

  nsingle = ntotal = 0;
}

/* ---------------------------------------------------------------------- */

FixEmit::~FixEmit()
{
  delete random;

  memory->destroy(c2list);
  memory->destroy(clist);
  memory->destroy(clistnum);
  memory->destroy(clistfirst);
}

/* ---------------------------------------------------------------------- */

int FixEmit::setmask()
{
  int mask = 0;
  mask |= START_OF_STEP;
  return mask;
}

/* ----------------------------------------------------------------------
   initialize per grid cell task data for grid cells I own
   do this from scratch every run in case grid or emit properties change
   callback to child emit style onecell() for each cell
   called from init() of child emit styles
------------------------------------------------------------------------- */

void FixEmit::init()
{
  particle->exist = 1;
  ntotal = 0;

  int dimension = domain->dimension;
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  nglocal = grid->nlocal;

  if (nglocal > nglocalmax) {
    memory->destroy(c2list);
    nglocalmax = nglocal;
    memory->create(c2list,nglocalmax,"emit:c2list");
  }

  // upsplit, split, sub cells store c2list flag
  // upsplit, split cells can store clist data, but only if have tasks
  // no tasks for a cell inside surface
  // no tasks if cell is entirely outside region bounding box

  int flag,ntaskcell,ntaskfirst;

  nlist = 0;
  ntaskfirst = 0;

  for (int icell = 0; icell < nglocal; icell++) {
    c2list[icell] = -1;
    if (cells[icell].nsplit <= 0) continue;
    if (cinfo[icell].type == INSIDE) continue;
    if (region && region->bboxflag) {
      flag = 1;
      if (cells[icell].hi[0] > region->extent_xlo &&
          cells[icell].lo[0] < region->extent_xhi) flag = 0;
      if (cells[icell].hi[1] > region->extent_ylo &&
          cells[icell].lo[1] < region->extent_yhi) flag = 0;
      if (dimension == 3) {
        if (cells[icell].hi[2] > region->extent_zlo &&
            cells[icell].lo[2] < region->extent_zhi) flag = 0;
      }
      if (flag) continue;
    }

    ntaskcell = create_task(icell);
    if (ntaskcell) {
      if (nlist == nlistmax) {
	nlistmax += DELTAGRID;
	memory->grow(clist,nlistmax,"emit:clist");
	memory->grow(clistnum,nlistmax,"emit:clistnum");
	memory->grow(clistfirst,nlistmax,"emit:clistfirst");
      }
      c2list[icell] = nlist;
      clist[nlist] = icell;
      clistnum[nlist] = ntaskcell;
      clistfirst[nlist] = ntaskfirst;
      ntaskfirst += ntaskcell;
      nlist++;
    }
  }

  active_current = 0;
}

/* ----------------------------------------------------------------------
   perform all particle emit tasks for cells I own
------------------------------------------------------------------------- */

void FixEmit::start_of_step()
{
  if (update->ntimestep % nevery) return;

  nsingle = 0;
  perform_task();
  ntotal += nsingle;
}

/* ----------------------------------------------------------------------
   add tasks for a new child cell added by adapt_grid or fix adapt
   similar logic to init()
------------------------------------------------------------------------- */

void FixEmit::add_grid_one(int icell, int flag)
{
  active_current = 0;

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  if (flag == 0) {
    if (nglocal == nglocalmax) {
      nglocalmax += DELTACELL;
      memory->grow(c2list,nglocalmax,"emit:c2list");
    }
    nglocal++;

    c2list[icell] = -1;
  }

  if (flag == 1) {
    if (cells[icell].nsplit <= 0) return;
    if (cinfo[icell].type == INSIDE) return;
    if (region && region->bboxflag) {
      int rflag = 1;
      if (cells[icell].hi[0] > region->extent_xlo &&
          cells[icell].lo[0] < region->extent_xhi) rflag = 0;
      if (cells[icell].hi[1] > region->extent_ylo &&
          cells[icell].lo[1] < region->extent_yhi) rflag = 0;
      if (domain->dimension == 3) {
        if (cells[icell].hi[2] > region->extent_zlo &&
            cells[icell].lo[2] < region->extent_zhi) rflag = 0;
      }
      if (rflag) return;
    }
    
    int ntaskcell = create_task(icell);
    
    if (ntaskcell) {
      if (nlist == nlistmax) {
        nlistmax += DELTAGRID;
        memory->grow(clist,nlistmax,"emit:clist");
        memory->grow(clistnum,nlistmax,"emit:clistnum");
        memory->grow(clistfirst,nlistmax,"emit:clistfirst");
      }
      c2list[icell] = nlist;
      clist[nlist] = icell;
      clistnum[nlist] = ntaskcell;
      clistfirst[nlist] = ntask - ntaskcell;
      nlist++;
    }
  }
}

/* ----------------------------------------------------------------------
   pack task count and tasks for grid cell icell into buf
   return byte count of amount packed
   if not memflag, only return count, do not fill buf
------------------------------------------------------------------------- */

int FixEmit::pack_grid_one(int icell, char *buf, int memflag)
{
  active_current = 0;

  char *ptr = buf;

  // ntaskcell = task count for icell

  int ilist = c2list[icell];

  int ntaskcell = 0;
  if (ilist >= 0) ntaskcell = clistnum[ilist];

  if (memflag) memcpy(ptr,&ntaskcell,sizeof(int));
  ptr += sizeof(int);
  ptr = ROUNDUP(ptr);

  if (!ntaskcell) return ptr-buf;

  // pack each insert task via child class pack_task()

  int itask = clistfirst[ilist];
  for (int i = 0; i < ntaskcell; i++) {
    ptr += pack_task(itask,ptr,memflag);
    ptr = ROUNDUP(ptr);
    itask++;
  }

  return ptr-buf;
}

/* ----------------------------------------------------------------------
   unpack task flag and tasks for grid cell icell from buf
   also unpack cellface data for flagged faces
   return byte count of amount unpacked
------------------------------------------------------------------------- */

int FixEmit::unpack_grid_one(int icell, char *buf)
{
  active_current = 0;

  char *ptr = buf;

  if (nglocal == nglocalmax) grow_percell(1);

  int ntaskcell;
  memcpy(&ntaskcell,ptr,sizeof(int));
  ptr += sizeof(int);
  ptr = ROUNDUP(ptr);

  if (ntaskcell) c2list[icell] = nlist;
  else c2list[icell] = -1;
  nglocal++;

  // if new cell is split cell, set c2list for subcells to -1

  int nsplit = grid->cells[icell].nsplit;
  if (nsplit > 1) {
    if (nglocal+nsplit > nglocalmax) grow_percell(nsplit);
    for (int i = 0; i < nsplit; i++) c2list[nglocal++] = -1;
  }

  if (!ntaskcell) return ptr-buf;

  // add list entry to all clist vectors

  if (nlist == nlistmax) grow_list();
  clist[nlist] = icell;
  clistnum[nlist] = ntaskcell;
  if (nlist) clistfirst[nlist] = clistfirst[nlist-1] + clistnum[nlist-1];
  else clistfirst[nlist] = 0;
  nlist++;

  // unpack each insert task via unpack_task() provided by child class

  for (int i = 0; i < ntaskcell; i++) {
    ptr += unpack_task(ptr,icell);
    ptr = ROUNDUP(ptr);
  }

  return ptr-buf;
}

/* ----------------------------------------------------------------------
   compress per-cell data due to cells migrating to new procs
   criteria for keeping/discarding a cell is same as in Grid::compress()
   this keeps final ordering of per-cell data consistent with Grid class
------------------------------------------------------------------------- */

void FixEmit::compress_grid()
{
  active_current = 0;

  int me = comm->me;
  Grid::ChildCell *cells = grid->cells;

  // keep an unsplit or split cell if staying on this proc
  // keep a sub cell if its split cell is staying on this proc
  // kept upsplit, split, sub cells store c2list flag
  // kept upsplit, split cells store clist data only if have tasks

  int ncurrent = nglocal;
  nglocal = 0;

  int oldlist,ntaskcell,oldntaskfirst;
  int ntaskfirst = 0;
  nlist = ntask = 0;

  for (int icell = 0; icell < ncurrent; icell++) {
    if (cells[icell].nsplit >= 1) {
      if (cells[icell].proc != me) continue;
      if (c2list[icell] < 0) {
	c2list[nglocal++] = -1;
	continue;
      }
    } else {
      int isplit = cells[icell].isplit;
      if (cells[grid->sinfo[isplit].icell].proc != me) continue;
      c2list[nglocal++] = -1;
      continue;
    }

    oldlist = c2list[icell];
    ntaskcell = clistnum[oldlist];
    oldntaskfirst = clistfirst[oldlist];

    c2list[nglocal] = nlist;
    clist[nlist] = nglocal;
    clistnum[nlist] = ntaskcell;
    clistfirst[nlist] = ntaskfirst;
    copy_task(nglocal,ntaskcell,ntaskfirst,oldntaskfirst);
    nglocal++;
    ntaskfirst += ntaskcell;
    nlist++;
  }
}

/* ----------------------------------------------------------------------
   insure c2list allocated long enough for N new cells
------------------------------------------------------------------------- */

void FixEmit::grow_percell(int n)
{
  while (nglocal+n > nglocalmax) nglocalmax += DELTAGRID;
  memory->grow(c2list,nglocalmax,"emit:c2list");
}

/* ----------------------------------------------------------------------
   grow all Nlist vectors
------------------------------------------------------------------------- */

void FixEmit::grow_list()
{
  nlistmax += DELTAGRID;
  memory->grow(clist,nlistmax,"emit:clist");
  memory->grow(clistnum,nlistmax,"emit:clistnum");
  memory->grow(clistfirst,nlistmax,"emit:clistfirst");
}

/* ----------------------------------------------------------------------
   calculate flux of particles of a species with vscale/fraction
     entering a grid cell
   indot = vstream dotted into face normal, assumed to be >= 0.0
   scosine = s cos(theta) in Bird notation where vscale = 1/beta
   see Bird 1994, eq 4.22
   bounding by -3.0 allows backflow influx of particles opposite to
     streaming velocity up to a reasonable limit
   if did not bound, would rarely emit a particle, but when do,
     could take too many iterations of double do while loop,
     e.g. in FixEmitFace::perform_task(), 
     to generate an inward velocity for the particle
------------------------------------------------------------------------- */

double FixEmit::mol_inflow(double indot, double vscale, double fraction)
{
  double scosine = indot / vscale;
  if (scosine < -3.0) return 0.0;
  double inward_number_flux = vscale*fraction *
    (exp(-scosine*scosine) + sqrt(MY_PI)*scosine*(1.0 + erf(scosine))) / 
    (2*sqrt(MY_PI));
  return inward_number_flux;
}

/* ----------------------------------------------------------------------
   process optional keywords common to all emit styles
   pass unrecognized keyword back to child option() method
   called from constructor of child emit styles
------------------------------------------------------------------------- */

int FixEmit::subsonic_temperature_check(int flag, double tempmax)
{
  int allflag;
  MPI_Allreduce(&flag,&allflag,1,MPI_INT,MPI_SUM,world);
  if (allflag) {
    double allmax;
    MPI_Allreduce(&tempmax,&allmax,1,MPI_DOUBLE,MPI_MAX,world);
    if (comm->me == 0) {
      char str[128];
      sprintf(str,"Excessive subsonic thermal temp = %g",allmax);
      error->warning(FLERR,str);
    }
    return 1;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   process optional keywords common to all emit styles
   pass unrecognized keyword back to child option() method
   called from constructor of child emit styles
------------------------------------------------------------------------- */

void FixEmit::options(int narg, char **arg)
{
  nevery = 1;
  perspecies = 1;
  region = NULL;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"nevery") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix emit command");
      nevery = atoi(arg[iarg+1]);
      if (nevery <= 0) error->all(FLERR,"Illegal fix emit command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"perspecies") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix emit command");
      if (strcmp(arg[iarg+1],"yes") == 0) perspecies = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) perspecies = 0;
      else error->all(FLERR,"Illegal fix emit command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix emit command");
      int iregion = domain->find_region(arg[iarg+1]);
      if (iregion < 0) 
        error->all(FLERR,"Fix emit region does not exist");
      region = domain->regions[iregion];
      iarg += 2;

    } else iarg += option(narg-iarg,&arg[iarg]);
  }
}

/* ----------------------------------------------------------------------
   process unknown keyword
------------------------------------------------------------------------- */

int FixEmit::option(int narg, char **arg)
{
  error->all(FLERR,"Illegal fix emit command");
  return 0;
}

/* ----------------------------------------------------------------------
   return one-step or total count of particle insertions
------------------------------------------------------------------------- */

double FixEmit::compute_vector(int i)
{
  double one,all;
  
  if (i == 0) one = nsingle;
  else one = ntotal;
  MPI_Allreduce(&one,&all,1,MPI_DOUBLE,MPI_SUM,world);
  return all;
}

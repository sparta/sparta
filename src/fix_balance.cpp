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

#include "string.h"
#include "stdlib.h"
#include "fix_balance.h"
#include "balance_grid.h"
#include "update.h"
#include "grid.h"
#include "particle.h"
#include "comm.h"
#include "rcb.h"
#include "modify.h"
#include "compute.h"
#include "random_mars.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{RANDOM,PROC,BISECTION};
enum{CELL,PARTICLE};

#define ZEROPARTICLE 0.1

/* ---------------------------------------------------------------------- */

FixBalance::FixBalance(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal fix balance command");

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;

  // parse arguments

  nevery = atoi(arg[2]);
  thresh = atof(arg[3]);

  if (strcmp(arg[4],"random") == 0) {
    if (narg != 5) error->all(FLERR,"Illegal fix balance command");
    bstyle = RANDOM;
  } else if (strcmp(arg[4],"proc") == 0) {
    if (narg != 5) error->all(FLERR,"Illegal fix balance command");
    bstyle = PROC;
  } else if (strcmp(arg[4],"rcb") == 0) {
    if (narg != 6) error->all(FLERR,"Illegal fix balance command");
    bstyle = BISECTION;
    if (strcmp(arg[5],"cell") == 0) rcbwt = CELL;
    else if (strcmp(arg[5],"mol") == 0) rcbwt = PARTICLE;
  } else error->all(FLERR,"Illegal fix balance command");

  // error check

  if (nevery < 0 || thresh < 1.0)
    error->all(FLERR,"Illegal fix balance command");

  me = comm->me;
  nprocs = comm->nprocs;

  // create instance of RNG or RCB

  random = NULL;
  rcb = NULL;

  if (bstyle == RANDOM || bstyle == PROC) 
    random = new RanPark(update->ranmaster->uniform()); 
  if (bstyle == BISECTION) rcb = new RCB(sparta);

  // compute initial outputs

  imbfinal = imbprev = imbalance_nlocal(maxperproc);
}

/* ---------------------------------------------------------------------- */

FixBalance::~FixBalance()
{
  delete random;
  delete rcb;
}

/* ---------------------------------------------------------------------- */

int FixBalance::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBalance::init()
{
  // error b/c acquire_ghosts() is a no-op in this case

  if (bstyle != BISECTION && grid->cutoff >= 0.0)
    error->all(FLERR,"Cannot use non-rcb fix balance with a grid cutoff");
}

/* ----------------------------------------------------------------------
   perform dynamic load balancing
------------------------------------------------------------------------- */

void FixBalance::end_of_step()
{
  // return if imbalance < threshhold

  imbnow = imbalance_nlocal(maxperproc);
  if (imbnow <= thresh) return;
  imbprev = imbnow;

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  int nglocal = grid->nlocal;

  // re-assign each of my local child cells to a proc
  // only assign unsplit and split cells
  // do not assign sub-cells since they migrate with their split cell
  // set nmigrate = # of cells that will migrate to a new proc
  // reset proc field in cells for migrating cells

  int nmigrate = 0;

  if (bstyle == RANDOM) {
    int newproc;
    for (int icell = 0; icell < nglocal; icell++) {
      if (cells[icell].nsplit <= 0) continue;
      newproc = nprocs * random->uniform();
      if (newproc != cells[icell].proc) nmigrate++;
      cells[icell].proc = newproc;
    }

  } else if (bstyle == PROC) {
    int newproc = nprocs * random->uniform();
    for (int icell = 0; icell < nglocal; icell++) {
      if (cells[icell].nsplit <= 0) continue;
      if (newproc != cells[icell].proc) nmigrate++;
      cells[icell].proc = newproc;
      newproc++;
      if (newproc == nprocs) newproc = 0;
    }

  } else if (bstyle == BISECTION) {
    double **x;
    memory->create(x,nglocal,3,"balance:x");

    double *lo,*hi;

    int nbalance = 0;
    for (int icell = 0; icell < nglocal; icell++) {
      if (cells[icell].nsplit <= 0) continue;
      lo = cells[icell].lo;
      hi = cells[icell].hi;
      x[nbalance][0] = 0.5*(lo[0]+hi[0]);
      x[nbalance][1] = 0.5*(lo[1]+hi[1]);
      x[nbalance][2] = 0.5*(lo[2]+hi[2]);
      nbalance++;
    }

    double *wt = NULL;
    if (rcbwt == PARTICLE) {
      int n;
      memory->create(wt,nglocal,"balance:wt");
      nbalance = 0;
      for (int icell = 0; icell < nglocal; icell++) {
        if (cells[icell].nsplit <= 0) continue;
        n = cinfo[icell].count;
        if (n) wt[nbalance++] = n;
        else wt[nbalance++] = ZEROPARTICLE;
      }
    }

    rcb->compute(nbalance,x,wt);
    rcb->invert();

    nbalance = 0;
    int *sendproc = rcb->sendproc;
    for (int icell = 0; icell < nglocal; icell++) {
      if (cells[icell].nsplit <= 0) continue;
      cells[icell].proc = sendproc[nbalance++];
    }
    nmigrate = nbalance - rcb->nkeep;

    memory->destroy(x);
    memory->destroy(wt);
  }

  if (bstyle == BISECTION) grid->clumped = 1;
  else grid->clumped = 0;

  // sort particles

  particle->sort();

  // migrate grid cells and their particles to new owners
  // invoke grid methods to complete grid setup

  grid->unset_neighbors();
  grid->remove_ghosts();
  comm->migrate_cells(nmigrate);

  grid->setup_owned();
  grid->acquire_ghosts();
  grid->reset_neighbors();
  comm->reset_neighbors();

  // reallocate per grid cell arrays in any relevant computes

  Compute **compute = modify->compute;
  for (int i = 0; i < modify->ncompute; i++)
    if (compute[i]->per_grid_flag)
      compute[i]->reallocate();

  // final imbalance factor

  imbfinal = imbalance_nlocal(maxperproc);
}

/* ----------------------------------------------------------------------
   calculate imbalance based on current particle count
   return max = max particles per proc
   return imbalance factor = max per proc / ave per proc
------------------------------------------------------------------------- */

double FixBalance::imbalance_nlocal(int &max)
{
  bigint n,nglobal;
  n = particle->nlocal;
  MPI_Allreduce(&n,&nglobal,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  MPI_Allreduce(&particle->nlocal,&max,1,MPI_INT,MPI_MAX,world);
  double imbalance = 1.0;
  if (max) imbalance = max / (1.0 * nglobal / nprocs);
  return imbalance;
}

/* ----------------------------------------------------------------------
   return imbalance factor after last rebalance
------------------------------------------------------------------------- */

double FixBalance::compute_scalar()
{
  return imbfinal;
}

/* ----------------------------------------------------------------------
   return stats for last rebalance
------------------------------------------------------------------------- */

double FixBalance::compute_vector(int i)
{
  if (i == 0) return (double) maxperproc;
  return imbprev;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

double FixBalance::memory_usage()
{
  double bytes = 0.0;
  //double bytes = irregular->memory_usage();
  return bytes;
}

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
#include "output.h"
#include "dump.h"
#include "random_mars.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"
#include "timer.h"

// DEBUG
#include "surf.h"

using namespace SPARTA_NS;

enum{RANDOM,PROC,BISECTION};
enum{CELL,PARTICLE,TIME};

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
    else if (strcmp(arg[5],"part") == 0) rcbwt = PARTICLE;
    else if (strcmp(arg[5],"time") == 0) rcbwt = TIME;
    else error->all(FLERR,"Illegal fix balance command");
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

  last = 0.0;
  imbfinal = imbprev = imbalance_factor(maxperproc);
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

  last = 0.0;
  timer->init();
}

/* ----------------------------------------------------------------------
   perform dynamic load balancing
------------------------------------------------------------------------- */

void FixBalance::end_of_step()
{
  // return if imbalance < threshhold

  imbnow = imbalance_factor(maxperproc);
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
      if (!particle->sorted) particle->sort();
      memory->create(wt,nglocal,"balance:wt");
      int n;
      nbalance = 0;
      for (int icell = 0; icell < nglocal; icell++) {
        if (cells[icell].nsplit <= 0) continue;
        n = cinfo[icell].count;
        if (n) wt[nbalance++] = n;
        else wt[nbalance++] = ZEROPARTICLE;
      }
    } else if (rcbwt == TIME) {
      memory->create(wt,nglocal,"balance:wt");
      timer_cell_weights(wt);
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

  if (nprocs == 1 || bstyle == BISECTION) grid->clumped = 1;
  else grid->clumped = 0;

  // sort particles

  if (!particle->sorted) particle->sort();

  // migrate grid cells and their particles to new owners
  // invoke grid methods to complete grid setup

  grid->unset_neighbors();
  grid->remove_ghosts();

  comm->migrate_cells(nmigrate);
  grid->hashfilled = 0;

  grid->setup_owned();
  grid->acquire_ghosts();

  grid->reset_neighbors();
  comm->reset_neighbors();

  // notify all classes that store per-grid data that grid may have changed
  
  grid->notify_changed();

  // final imbalance factor

  if (bstyle == BISECTION && rcbwt == TIME)
    imbfinal = 0.0; // can't compute imbalance from timers since grid cells moved
  else
    imbfinal = imbalance_factor(maxperproc);
}

/* ----------------------------------------------------------------------
   calculate imbalance based on current particle count
   return maxcost = max particles per proc or CPU time per proc
   return imbalance factor = max per proc / ave per proc
------------------------------------------------------------------------- */

double FixBalance::imbalance_factor(double &maxcost)
{
  double mycost,totalcost;
  double mycost_proc_weighted,maxcost_proc_weighted,nprocs_weighted;

  if (bstyle == BISECTION && rcbwt == TIME) {
    timer_cost();
    mycost = my_timer_cost;
  } else mycost = particle->nlocal;

  MPI_Allreduce(&mycost,&totalcost,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&mycost,&maxcost,1,MPI_DOUBLE,MPI_MAX,world);

  double imbalance = 1.0;
  if (maxcost) imbalance = maxcost / (totalcost / nprocs);
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
  if (i == 0) return maxperproc;
  return imbprev;
}

/* -------------------------------------------------------------------- */

void FixBalance::timer_cost()
{
  // my_timer_cost = CPU time for relevant timers since last invocation

  my_timer_cost = -last;
  my_timer_cost += timer->array[TIME_MOVE];
  my_timer_cost += timer->array[TIME_SORT];
  my_timer_cost += timer->array[TIME_COLLIDE];
  my_timer_cost += timer->array[TIME_MODIFY];

  // last = time up to this point

  last += my_timer_cost;
}

/* -------------------------------------------------------------------- */

void FixBalance::timer_cell_weights(double *weight)
{
  // localwt = weight assigned to each owned grid cell
  // just return if no time yet tallied

  double maxcost;
  MPI_Allreduce(&my_timer_cost,&maxcost,1,MPI_DOUBLE,MPI_MAX,world);
  if (maxcost <= 0.0) {
    memory->destroy(weight);
    weight = NULL;
    return;
  }

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  int nglocal = grid->nlocal;

  double localwt_total = 0.0;
  if (nglocal) localwt_total = my_timer_cost/nglocal;
  if (nglocal && localwt_total <= 0.0) error->one(FLERR,"Balance weight <= 0.0");

  if (!particle->sorted) particle->sort();
  double wttotal = 0;
  int nbalance = 0;
  double* localwt;
  memory->create(localwt,nglocal,"imbalance_time:localwt");
  for (int icell = 0; icell < nglocal; icell++) {
    localwt[icell] = 0.0;
    if (cells[icell].nsplit <= 0) continue;
    int n = cinfo[icell].count;
    if (n) localwt[nbalance++] = n;
    else localwt[nbalance++] = ZEROPARTICLE;
    wttotal += localwt[nbalance-1];
  }

  for (int icell = 0; icell < nglocal; icell++) 
    weight[icell] = my_timer_cost*localwt[icell]/wttotal;

  memory->destroy(localwt);
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

double FixBalance::memory_usage()
{
  double bytes = 0.0;
  // tally wt vector?
  return bytes;
}

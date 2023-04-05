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

#include "string.h"
#include "stdlib.h"
#include "fix_balance.h"
#include "balance_grid.h"
#include "update.h"
#include "grid.h"
#include "particle.h"
#include "domain.h"
#include "comm.h"
#include "rcb.h"
#include "modify.h"
#include "compute.h"
#include "output.h"
#include "dump.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "memory.h"
#include "error.h"
#include "timer.h"
#include "fix_adapt.h"

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

  int iarg;
  if (strcmp(arg[4],"random") == 0) {
    bstyle = RANDOM;
    iarg = 5;
  } else if (strcmp(arg[4],"proc") == 0) {
    bstyle = PROC;
    iarg = 5;
  } else if (strcmp(arg[4],"rcb") == 0) {
    if (narg < 6) error->all(FLERR,"Illegal fix balance command");
    bstyle = BISECTION;
    if (strcmp(arg[5],"cell") == 0) rcbwt = CELL;
    else if (strcmp(arg[5],"part") == 0) rcbwt = PARTICLE;
    else if (strcmp(arg[5],"time") == 0) rcbwt = TIME;
    else error->all(FLERR,"Illegal fix balance command");
    iarg = 6;
  } else error->all(FLERR,"Illegal fix balance command");

  // optional args

  strcpy(eligible,"xyz");
  rcbflip = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"axes") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix balance command");
      if (strlen(arg[iarg+1]) > 3)
        error->all(FLERR,"Illegal fix balance command");
      strcpy(eligible,arg[iarg+1]);
      int xdim = 0;
      int ydim = 0;
      int zdim = 0;
      if (strchr(eligible,'x')) xdim = 1;
      if (strchr(eligible,'y')) ydim = 1;
      if (strchr(eligible,'z')) zdim = 1;
      if (zdim && domain->dimension == 2)
        error->all(FLERR,"Illegal balance_grid command");
      if (xdim+ydim+zdim != strlen(eligible))
        error->all(FLERR,"Illegal fix balance command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"flip") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix balance command");
      if (strcmp(arg[iarg+1],"yes") == 0) rcbflip = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) rcbflip = 0;
      else error->all(FLERR,"Illegal fix balance command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix balance command");
  }

  // error check

  if (nevery < 0 || thresh < 1.0)
    error->all(FLERR,"Illegal fix balance command");

  me = comm->me;
  nprocs = comm->nprocs;

  // create instance of RNG or RCB

  random = NULL;
  rcb = NULL;

  if (bstyle == RANDOM || bstyle == PROC)
    random = new RanKnuth(update->ranmaster->uniform());
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

  // check if fix balance rcb time is after fix adapt with coarsening

  if (rcbwt == TIME) {
    int coarsen_flag = 0;
    int after_flag = 0;
    for (int ifix = 0; ifix < modify->nfix; ifix++) {
      if (strstr(modify->fix[ifix]->style,"adapt") != NULL)
        if (((FixAdapt*)modify->fix[ifix])->coarsen_flag)
          coarsen_flag = 1;

      if (strstr(modify->fix[ifix]->style,"balance") != NULL)
        if (coarsen_flag) {
          after_flag = 1;
          break;
        }
    }

    if (after_flag && comm->me == 0) {
      error->warning(FLERR,"Using fix_adapt coarsen before fix balance "
        "rcb time may make the accumulated timing data less accurate");
    }
  }

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

    rcb->compute(nbalance,x,wt,eligible,rcbflip);
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
  // some fixes have post migration operations to perform

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

void FixBalance::timer_cell_weights(double* &weight)
{
  // localwt = weight assigned to each owned grid cell
  // just return if no time yet tallied

  double maxcost;
  MPI_Allreduce(&my_timer_cost,&maxcost,1,MPI_DOUBLE,MPI_MAX,world);
  if (maxcost <= 0.0) {
    memory->destroy(weight);
    weight = NULL;
    if (comm->me == 0) {
      error->warning(FLERR,"No time history accumulated for fix balance "
        "rcb time, using rcb cell option instead");
    }
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

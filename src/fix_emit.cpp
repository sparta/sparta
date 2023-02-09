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

#include "stdlib.h"
#include "string.h"
#include "fix_emit.h"
#include "update.h"
#include "domain.h"
#include "region.h"
#include "grid.h"
#include "comm.h"
#include "random_mars.h"
#include "random_knuth.h"
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
  random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,me,100);

  // counters common to all emit styles for output from fix

  nsingle = ntotal = 0;
}

/* ---------------------------------------------------------------------- */

FixEmit::~FixEmit()
{
  if (copymode) return;

  delete random;
}

/* ---------------------------------------------------------------------- */

int FixEmit::setmask()
{
  int mask = 0;
  mask |= START_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixEmit::init()
{
  particle->exist = 1;
  ntotal = 0;
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
   grid changed operation
   invoke create_tasks() to rebuild entire task list
   invoked after per-processor list of grid cells has changed
------------------------------------------------------------------------- */

void FixEmit::grid_changed()
{
  create_tasks();
}

/* ----------------------------------------------------------------------
   initialize per grid cell task data for grid cells I own
   calls back to child create_task() for each cell
   called from init() of child emit styles, in case grid or emit parama change
   also called after on-the-fly grid adaptation or load balancing
------------------------------------------------------------------------- */

void FixEmit::create_tasks()
{
  int dimension = domain->dimension;
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  int nglocal = grid->nlocal;

  // no tasks for a cell inside surface
  // no tasks if cell is entirely outside region bounding box

  ntask = 0;

  int rflag;

  for (int icell = 0; icell < nglocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    if (cinfo[icell].type == INSIDE) continue;
    if (region && region->bboxflag) {
      rflag = 1;
      if (cells[icell].hi[0] > region->extent_xlo &&
          cells[icell].lo[0] < region->extent_xhi) rflag = 0;
      if (cells[icell].hi[1] > region->extent_ylo &&
          cells[icell].lo[1] < region->extent_yhi) rflag = 0;
      if (dimension == 3) {
        if (cells[icell].hi[2] > region->extent_zlo &&
            cells[icell].lo[2] < region->extent_zhi) rflag = 0;
      }
      if (rflag) continue;
    }

    create_task(icell);
  }

  active_current = 0;
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
    (exp(-scosine*scosine) + MY_PIS*scosine*(1.0 + erf(scosine))) /
    (2*MY_PIS);
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

int FixEmit::option(int, char **)
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

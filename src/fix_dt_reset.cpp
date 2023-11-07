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

/* ----------------------------------------------------------------------
   Contributing author: Alan Stagg (Sandia)
------------------------------------------------------------------------- */

#include "fix_dt_reset.h"
#include "string.h"
#include "update.h"
#include "grid.h"
#include "modify.h"
#include "compute.h"
#include "memory.h"
#include "error.h"
#include <iostream>
using namespace SPARTA_NS;

enum{COMPUTE,FIX};
enum{NORESET,RESET,RESETWARN};

#define INVOKED_PER_GRID 16
#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

FixDtReset::FixDtReset(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg != 6) error->all(FLERR,"Illegal fix dt/reset command");

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;

  nevery = atoi(arg[2]);
  weight = atof(arg[4]);
  resetflag = atoi(arg[5]);

  // error checks

  if (nevery <= 0) error->all(FLERR,"Illegal fix dt/reset command");
  if (weight < 0.0 || weight > 1.0)
    error->all(FLERR,"Fix dt/reset weight is not within [0,1]");
  if (resetflag < 0 || resetflag > 2) error->all(FLERR,"Illegal fix dt/reset command");

  // process per grid timestep arg as compute or fix

  if (strncmp(arg[3],"c_",2) == 0) step_which = COMPUTE;
  else if (strncmp(arg[3],"f_",2) == 0) step_which = FIX;
  else error->all(FLERR,"Invalid step in fix dt/reset command");

  int n = strlen(arg[3]);
  id_step = new char[n];
  strcpy(id_step,&arg[3][2]);
  char *ptr = strchr(id_step,'[');
  if (ptr) {
    if (id_step[strlen(id_step)-1] != ']')
      error->all(FLERR,"Invalid step in fix dt/reset command");
    step_index = atoi(ptr+1);
    *ptr = '\0';
  } else step_index = 0;

  if (step_which == COMPUTE) {
    int m = modify->find_compute(id_step);
    if (m < 0)
      error->all(FLERR,"Could not find fix dt/reset compute ID");
    if (modify->compute[m]->per_grid_flag == 0)
      error->all(FLERR,"Fix dt/reset compute does not "
                 "compute per-grid info");
    if (step_index == 0 && modify->compute[m]->size_per_grid_cols > 0)
      error->all(FLERR,"Fix dt/reset compute does not"
                 "compute per-grid vector");
    if (step_index > 0 && modify->compute[m]->size_per_grid_cols == 0)
      error->all(FLERR,"Fix dt/reset compute does not "
                 "compute per-grid array");
    if (step_index > 0 && step_index > modify->compute[m]->size_per_grid_cols)
      error->all(FLERR,"Fix dt/reset compute array is "
                 "accessed out-of-range");
  } else if (step_which == FIX) {
    int m = modify->find_fix(id_step);
    if (m < 0) error->all(FLERR,"Could not find fix dt/reset fix ID");
    if (modify->fix[m]->per_grid_flag == 0)
      error->all(FLERR,"Fix dt/reset fix does not "
                 "compute per-grid info");
    if (step_index == 0 && modify->fix[m]->size_per_grid_cols > 0)
      error->all(FLERR,"Fix dt/reset fix does not "
                 "compute per-grid vector");
    if (step_index > 0 && modify->fix[m]->size_per_grid_cols == 0)
      error->all(FLERR,"Fix dt/reset fix does not "
                 "compute per-grid array");
    if (step_index > 0 && step_index > modify->fix[m]->size_per_grid_cols)
      error->all(FLERR,"Fix dt/reset fix array is "
                 "accessed out-of-range");
  }

  // initializations

  dtmin = 0.0;
  dtmax = 0.0;
  dtave = 0.0;
  dtnew = 0.0;

  maxgrid = 0;
  gridstep = NULL;
}

/* ---------------------------------------------------------------------- */

FixDtReset::~FixDtReset()
{
  if (copymode) return;

  delete [] id_step;
  memory->destroy(gridstep);
}

/* ---------------------------------------------------------------------- */

void FixDtReset::init()
{
  // cstep/fstep = compute or fix that calculates per grid cell timestep

  if (step_which == COMPUTE) {
    int icompute = modify->find_compute(id_step);
    if (icompute < 0)
      error->all(FLERR,"Could not find fix dt/reset compute ID");
    cstep = modify->compute[icompute];
  } else if (step_which == FIX) {
    int ifix = modify->find_fix(id_step);
    if (ifix < 0)
      error->all(FLERR,"Could not find fix dt/reset fix ID");
    fstep = modify->fix[ifix];
  }

  for (int i = 0; i < modify->nfix; i++)
    if (strstr(modify->fix[i]->style,"adapt") != NULL)
      if (modify->fix[i]->nevery < nevery)
        error->warning(FLERR,"Adapting grid more often than resetting global timestep");
}

/* ---------------------------------------------------------------------- */

int FixDtReset::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDtReset::end_of_step()
{
  // reallocate per-grid storage if necessary

  int nglocal = grid->nlocal;

  if (nglocal > maxgrid) {
    maxgrid = grid->nlocal;
    memory->destroy(gridstep);
    memory->create(gridstep,maxgrid,"dt/reset:gridstep");
  }

  // check that fix is computed at compatible time

  if (step_which == FIX && update->ntimestep % fstep->per_grid_freq)
    error->all(FLERR,"Fix dt/reset fix not computed at compatible time");

  // grab per grid cell timestep from compute or fix, invoke compute if needed

  if (step_which == COMPUTE) {
    if (!(cstep->invoked_flag & INVOKED_PER_GRID)) {
      cstep->compute_per_grid();
      cstep->invoked_flag |= INVOKED_PER_GRID;
    }

    if (step_index == 0)
      memcpy(gridstep,cstep->vector_grid,nglocal*sizeof(double));
    else {
      int index = step_index-1;
      double **array = cstep->array_grid;
      for (int i = 0; i < nglocal; i++)
        gridstep[i] = array[i][index];
    }
  } else if (step_which == FIX) {
    if (step_index == 0) {
      memcpy(gridstep,fstep->vector_grid,nglocal*sizeof(double));
    } else {
      double **array = fstep->array_grid;
      int index = step_index-1;
      for (int i = 0; i < nglocal; i++)
        gridstep[i] = array[i][index];
    }
  }

  // set dtmin,dtmax,dtave
  // skip cells whose timestep is zero (e.g. no particles)

  double dtmin_me = BIG;
  double dtmax_me = 0.0;
  double dtsum_me = 0.0;
  int count = 0;

  for (int i = 0; i < nglocal; i++) {
    if (gridstep[i] == 0.0) continue;
    dtmin_me = MIN(dtmin_me,gridstep[i]);
    dtmax_me = MAX(dtmax_me,gridstep[i]);
    dtsum_me += gridstep[i];
    count++;
  }

  bigint bcount_me = count;
  bigint bcount;
  MPI_Allreduce(&bcount_me,&bcount,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  if (bcount == 0) return;  // done if no cell computed a timestep

  double min_max_data_me[2];
  double min_max_data[2];
  min_max_data_me[0] = dtmin_me;
  min_max_data_me[1] = 1./dtmax_me;
  MPI_Allreduce(&min_max_data_me,&min_max_data,2,MPI_DOUBLE,MPI_MIN,world);
  dtmin = min_max_data[0];
  dtmax = 1./min_max_data[1];
  MPI_Allreduce(&dtsum_me,&dtave,1,MPI_DOUBLE,MPI_SUM,world);
  dtave /= bcount;

  // calculate new global timestep

  dtnew = (1.0-weight) * dtmin + weight * dtave;

  // reset global timestep if requested
  // also reset the global time

  if (resetflag) {
    update->time += (update->ntimestep - update->time_last_update) * update->dt;
    update->time_last_update = update->ntimestep;
    update->dt = dtnew;
  }
}

/* ----------------------------------------------------------------------
   return optimal calculated timestep
------------------------------------------------------------------------- */

double FixDtReset::compute_scalar()
{
  return dtnew;
}

/* ----------------------------------------------------------------------
   return min/max/avg per grid cell timesteps calaculated by compute
------------------------------------------------------------------------- */

double FixDtReset::compute_vector(int index)
{
  if (index == 0) return dtmin;
  else if (index == 1) return dtmax;
  else if (index == 2) return dtave;
  return 0.0;
}

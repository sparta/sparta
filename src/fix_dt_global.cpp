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
#include "fix_dt_global.h"
#include "update.h"
#include "input.h"
#include "grid.h"
#include "modify.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include <iostream>

using namespace SPARTA_NS;

#define INVOKED_PER_GRID 16
#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

FixDtGlobal::FixDtGlobal(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg != 7) error->all(FLERR,"Illegal fix dt global command");
  nevery = atoi(arg[2]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix dt global command");

  id_lambda = NULL;
  id_usq = NULL;
  id_vsq = NULL;
  id_wsq = NULL;

  // lambda compute
  if (strncmp(arg[3],"c_",2) == 0) {
    int n = strlen(arg[3]);
    id_lambda = new char[n];
    strcpy(id_lambda,&arg[3][2]);
    char *ptr = strchr(id_lambda,'[');
    if (ptr) {
      if (id_lambda[strlen(id_lambda)-1] != ']')
        error->all(FLERR,"Invalid lambda in fix dt_global command");
      lambdaindex = atoi(ptr+1);
      *ptr = '\0';
    } else lambdaindex = 0;

    n = modify->find_compute(id_lambda);
    if (n < 0)
      error->all(FLERR,"Could not find fix dt_global compute ID");
    if (modify->compute[n]->per_grid_flag == 0)
      error->all(FLERR,"Fix dt_global compute does not "
                 "compute per-grid info");
    if (lambdaindex == 0 && modify->compute[n]->size_per_grid_cols > 0)
      error->all(FLERR,
                 "Fix dt_global compute does not "
                 "compute per-grid vector");
    if (lambdaindex > 0 && modify->compute[n]->size_per_grid_cols == 0)
      error->all(FLERR,
                 "Fix dt_global compute does not "
                 "compute per-grid array");
    if (lambdaindex > 0 && lambdaindex > modify->compute[n]->size_per_grid_cols)
      error->all(FLERR,"Fix dt_global compute vector is "
                 "accessed out-of-range");
  }
  else {
    error->all(FLERR,"Illegal fix dt global command");
  }

  // squared velocity fixes
  if ((strncmp(arg[4],"f_",2) == 0) &&
      (strncmp(arg[5],"f_",2) == 0) &&
      (strncmp(arg[6],"f_",2) == 0)) {

    // usq----------------------------
    int n = strlen(arg[4]);
    id_usq = new char[n];
    strcpy(id_usq,&arg[4][2]);

    char *ptr = strchr(id_usq,'[');
    if (ptr) {
      if (id_usq[strlen(id_usq)-1] != ']')
        error->all(FLERR,"Invalid usq in fix dt_global command");
      usqindex = atoi(ptr+1);
      *ptr = '\0';
    } else usqindex = 0;

    n = modify->find_fix(id_usq);
    if (n < 0) error->all(FLERR,"Could not find fix dt_global fix ID");
    if (modify->fix[n]->per_grid_flag == 0)
      error->all(FLERR,"Fix dt_global fix does not "
                 "compute per-grid info");
    if (usqindex == 0 && modify->fix[n]->size_per_grid_cols > 0)
      error->all(FLERR,"Fix dt_global fix does not "
                 "compute per-grid vector");
    if (usqindex > 0 && modify->fix[n]->size_per_grid_cols == 0)
      error->all(FLERR,"Fix dt_global fix does not "
                 "compute per-grid array");
    if (usqindex > 0 && usqindex > modify->fix[n]->size_per_grid_cols)
      error->all(FLERR,"Fix dt_global fix array is "
                 "accessed out-of-range");

    // vsq----------------------------
    n = strlen(arg[5]);
    id_vsq = new char[n];
    strcpy(id_vsq,&arg[5][2]);

    ptr = strchr(id_vsq,'[');
    if (ptr) {
      if (id_vsq[strlen(id_vsq)-1] != ']')
        error->all(FLERR,"Invalid vsq in fix dt_global command");
      vsqindex = atoi(ptr+1);
      *ptr = '\0';
    } else vsqindex = 0;

    n = modify->find_fix(id_vsq);
    if (n < 0) error->all(FLERR,"Could not find fix dt_global fix ID");
    if (modify->fix[n]->per_grid_flag == 0)
      error->all(FLERR,"Fix dt_global fix does not "
                 "compute per-grid info");
    if (vsqindex == 0 && modify->fix[n]->size_per_grid_cols > 0)
      error->all(FLERR,"Fix dt_global fix does not "
                 "compute per-grid vector");
    if (vsqindex > 0 && modify->fix[n]->size_per_grid_cols == 0)
      error->all(FLERR,"Fix dt_global fix does not "
                 "compute per-grid array");
    if (vsqindex > 0 && vsqindex > modify->fix[n]->size_per_grid_cols)
      error->all(FLERR,"Fix dt_global fix array is "
                 "accessed out-of-range");

    // wsq----------------------------
    n = strlen(arg[6]);
    id_wsq = new char[n];
    strcpy(id_wsq,&arg[6][2]);

    ptr = strchr(id_wsq,'[');
    if (ptr) {
      if (id_wsq[strlen(id_wsq)-1] != ']')
        error->all(FLERR,"Invalid wsq in fix dt_global command");
      wsqindex = atoi(ptr+1);
      *ptr = '\0';
    } else wsqindex = 0;

    n = modify->find_fix(id_wsq);
    if (n < 0) error->all(FLERR,"Could not find fix dt_global fix ID");
    if (modify->fix[n]->per_grid_flag == 0)
      error->all(FLERR,"Fix dt_global fix does not "
                 "compute per-grid info");
    if (wsqindex == 0 && modify->fix[n]->size_per_grid_cols > 0)
      error->all(FLERR,"Fix dt_global fix does not "
                 "compute per-grid vector");
    if (wsqindex > 0 && modify->fix[n]->size_per_grid_cols == 0)
      error->all(FLERR,"Fix dt_global fix does not "
                 "compute per-grid array");
    if (wsqindex > 0 && wsqindex > modify->fix[n]->size_per_grid_cols)
      error->all(FLERR,"Fix dt_global fix array is "
                 "accessed out-of-range");
  }
  else {
    error->all(FLERR,"Illegal fix dt global command");
  }

  MPI_Comm_rank(world,&me);

  nglocal = 0;
  lambda = usq = vsq = wsq = NULL;
}

/* ---------------------------------------------------------------------- */

FixDtGlobal::~FixDtGlobal()
{
  delete [] id_lambda;
  delete [] id_usq;
  delete [] id_vsq;
  delete [] id_wsq;
  memory->destroy(lambda);
  memory->destroy(usq);
  memory->destroy(vsq);
  memory->destroy(wsq);
}

/* ---------------------------------------------------------------------- */

void FixDtGlobal::init()
{
  if (me == 0)
    std::cout << "NEED TO ADD INITIAL ESTIMATION OF GLOBAL TIMESTEP HERE (TO BE "
              << "USED BEFORE FIXES ARE AVAILABLE)!!!!!!!!!!!!!!!!!\n";

  reallocate();

  // lambda compute
  int icompute = modify->find_compute(id_lambda);
  if (icompute < 0)
    error->all(FLERR,"Could not find fix dt_global compute ID");
  clambda = modify->compute[icompute];

  // usq fix
  int ifix = modify->find_fix(id_usq);
  if (ifix < 0)
    error->all(FLERR,"Could not find usq fix ID");
  fusq = modify->fix[ifix];

  // vsq fix
  ifix = modify->find_fix(id_vsq);
  if (ifix < 0)
    error->all(FLERR,"Could not find vsq fix ID");
  fvsq = modify->fix[ifix];

  // wsq fix
  ifix = modify->find_fix(id_wsq);
  if (ifix < 0)
    error->all(FLERR,"Could not find wsq fix ID");
  fwsq = modify->fix[ifix];
}

/* ---------------------------------------------------------------------- */

int FixDtGlobal::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDtGlobal::end_of_step()
{
  // check that required fixes are computed at compatible times
  if (update->ntimestep % fusq->per_grid_freq)
    error->all(FLERR,"usq fix not computed at compatible time");
  if (update->ntimestep % fvsq->per_grid_freq)
    error->all(FLERR,"vsq fix not computed at compatible time");
  if (update->ntimestep % fwsq->per_grid_freq)
    error->all(FLERR,"wsq fix not computed at compatible time");

  // get lambda from compute
  if (!(clambda->invoked_flag & INVOKED_PER_GRID)) {
    clambda->compute_per_grid();
    clambda->invoked_flag |= INVOKED_PER_GRID;
  }

  if (clambda->post_process_grid_flag)
    clambda->post_process_grid(lambdaindex,1,NULL,NULL,NULL,1);

  if (lambdaindex == 0 || clambda->post_process_grid_flag)
    memcpy(lambda,clambda->vector_grid,nglocal*sizeof(double));
  else {
    int index = lambdaindex-1;
    double **array = clambda->array_grid;
    for (int i = 0; i < nglocal; i++)
      lambda[i] = array[i][index];
  }

  // get usq from fix
  if (usqindex == 0) {
    memcpy(usq,fusq->vector_grid,nglocal*sizeof(double));
  } else {
    double **array = fusq->array_grid;
    int index = usqindex-1;
    for (int i = 0; i < nglocal; i++)
      usq[i] = array[i][index];
  }

  // get vsq from fix
  if (vsqindex == 0) {
    memcpy(vsq,fvsq->vector_grid,nglocal*sizeof(double));
  } else {
    double **array = fvsq->array_grid;
    int index = vsqindex-1;
    for (int i = 0; i < nglocal; i++)
      vsq[i] = array[i][index];
  }

  // get wsq from fix
  if (wsqindex == 0) {
    memcpy(wsq,fwsq->vector_grid,nglocal*sizeof(double));
  } else {
    double **array = fwsq->array_grid;
    int index = wsqindex-1;
    for (int i = 0; i < nglocal; i++)
      wsq[i] = array[i][index];
  }

  // compute cell desired timestep
  Grid::ChildCell *cells = grid->cells;
  double dtmin_this_pe = BIG;
  for (int i = 0; i < nglocal; ++i) {
    double mean_collision_time = lambda[i]/sqrt(usq[i] + vsq[i] + wsq[i]);
    cells[i].dt_desired = mean_collision_time;
    if (cells[i].dt_desired < dtmin_this_pe)
      dtmin_this_pe = cells[i].dt_desired;
  }

  double dtmin;
  MPI_Allreduce(&dtmin_this_pe,&dtmin,1,MPI_DOUBLE,MPI_MIN,world);
  grid->dt_global = dtmin;
  if (me == 0) {
    std::cout << "WARNING!!!!! values stored in dt_global and cells[i].dt_desired are placeholders!!!!!\n";
    std::cout << "...also figure out what it means for index values to be zero, e.g. lambdaindex\n";
  }
}

/* ----------------------------------------------------------------------
   reallocate arrays if nglocal has changed
   called by init() and load balancer
------------------------------------------------------------------------- */

void FixDtGlobal::reallocate()
{
  if (grid->nlocal == nglocal) return;

  nglocal = grid->nlocal;

  memory->destroy(lambda);
  memory->destroy(usq);
  memory->destroy(vsq);
  memory->destroy(wsq);
  memory->create(lambda,nglocal,"dt_global:lambda");
  memory->create(usq,nglocal,"dt_global:usq");
  memory->create(vsq,nglocal,"dt_global:vsq");
  memory->create(wsq,nglocal,"dt_global:wsq");

}

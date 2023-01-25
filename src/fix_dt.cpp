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

#include "fix_dt.h"
#include "string.h"
#include "update.h"
#include "grid.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "comm.h"
#include <iostream>

using namespace SPARTA_NS;

enum{COMPUTE,FIX};
enum{WARN, USE_CALCULATED_DT};
enum{INT,DOUBLE};

#define INVOKED_PER_GRID 16
#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

FixDt::FixDt(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  // arguments
  int narg_required = 12;
  if (narg < narg_required) error->all(FLERR,"Illegal fix dt command");
  scalar_flag = 1;
  global_freq = 1;

  nevery = atoi(arg[2]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix dt command");

  dt_global_weight = atof(arg[3]);
  transit_fraction = atof(arg[4]);
  collision_fraction = atof(arg[5]);

  if (strncmp(arg[6],"c_",2) != 0 && strncmp(arg[6],"f_",2) != 0)
    error->all(FLERR,"Illegal fix dt command");

  if ( !((strncmp(arg[7],"f_",2) == 0) &&
         (strncmp(arg[8],"f_",2) == 0) &&
         (strncmp(arg[9],"f_",2) == 0) &&
         (strncmp(arg[10],"f_",2) == 0)) )
    error->all(FLERR,"Illegal fix dt command");

  id_lambda = NULL;
  id_temp = NULL;
  id_usq = NULL;
  id_vsq = NULL;
  id_wsq = NULL;

  // lambda compute or fix
  int n = strlen(arg[6]);
  id_lambda = new char[n];
  strcpy(id_lambda,&arg[6][2]);
  char *ptr = strchr(id_lambda,'[');
  if (ptr) {
    if (id_lambda[strlen(id_lambda)-1] != ']')
      error->all(FLERR,"Invalid lambda in fix dt command");
    lambdaindex = atoi(ptr+1);
    *ptr = '\0';
  } else lambdaindex = 0;

  if (strncmp(arg[6],"c_",2) == 0) lambdawhich = COMPUTE;
  else lambdawhich = FIX;

  if (lambdawhich == COMPUTE) { // lambda compute
    n = modify->find_compute(id_lambda);
    if (n < 0)
      error->all(FLERR,"Could not find fix dt compute ID");
    if (modify->compute[n]->per_grid_flag == 0)
      error->all(FLERR,"Fix dt compute does not "
                 "compute per-grid info");
    if (lambdaindex == 0 && modify->compute[n]->size_per_grid_cols > 0)
      error->all(FLERR,
                 "Fix dt compute does not "
                 "compute per-grid vector");
    if (lambdaindex > 0 && modify->compute[n]->size_per_grid_cols == 0)
      error->all(FLERR,
                 "Fix dt compute does not "
                 "compute per-grid array");
    if (lambdaindex > 0 && lambdaindex > modify->compute[n]->size_per_grid_cols)
      error->all(FLERR,"Fix dt compute vector is "
                 "accessed out-of-range");
  } else { // lambda fix
    n = modify->find_fix(id_lambda);
    if (n < 0) error->all(FLERR,"Could not find fix dt fix ID");
    if (modify->fix[n]->per_grid_flag == 0)
      error->all(FLERR,"Fix dt fix does not "
                 "compute per-grid info");
    if (lambdaindex == 0 && modify->fix[n]->size_per_grid_cols > 0)
      error->all(FLERR,"Fix dt fix does not "
                 "compute per-grid vector");
    if (lambdaindex > 0 && modify->fix[n]->size_per_grid_cols == 0)
      error->all(FLERR,"Fix dt fix does not "
                 "compute per-grid array");
    if (lambdaindex > 0 && lambdaindex > modify->fix[n]->size_per_grid_cols)
      error->all(FLERR,"Fix dt fix array is "
                 "accessed out-of-range");
  }

  // temperature fix
  n = strlen(arg[7]);
  id_temp = new char[n];
  strcpy(id_temp,&arg[7][2]);

  ptr = strchr(id_temp,'[');
  if (ptr) {
    if (id_temp[strlen(id_temp)-1] != ']')
      error->all(FLERR,"Invalid temperature in fix dt command");
    tempindex = atoi(ptr+1);
    *ptr = '\0';
  } else tempindex = 0;

  n = modify->find_fix(id_temp);
  if (n < 0) error->all(FLERR,"Could not find fix dt tempeature fix ID");
  if (modify->fix[n]->per_grid_flag == 0)
    error->all(FLERR,"Fix dt fix does not "
               "compute per-grid info");
  if (tempindex == 0 && modify->fix[n]->size_per_grid_cols > 0)
    error->all(FLERR,"Fix dt fix does not "
               "compute per-grid vector");
  if (tempindex > 0 && modify->fix[n]->size_per_grid_cols == 0)
    error->all(FLERR,"Fix dt fix does not "
               "compute per-grid array");
  if (tempindex > 0 && tempindex > modify->fix[n]->size_per_grid_cols)
    error->all(FLERR,"Fix dt fix array is "
               "accessed out-of-range");

  // squared-velocity fixes

  // usq----------------------------
  n = strlen(arg[8]);
  id_usq = new char[n];
  strcpy(id_usq,&arg[8][2]);

  ptr = strchr(id_usq,'[');
  if (ptr) {
    if (id_usq[strlen(id_usq)-1] != ']')
      error->all(FLERR,"Invalid usq in fix dt command");
    usqindex = atoi(ptr+1);
    *ptr = '\0';
  } else usqindex = 0;

  n = modify->find_fix(id_usq);
  if (n < 0) error->all(FLERR,"Could not find fix dt usq fix ID");
  if (modify->fix[n]->per_grid_flag == 0)
    error->all(FLERR,"Fix dt fix does not "
               "compute per-grid info");
  if (usqindex == 0 && modify->fix[n]->size_per_grid_cols > 0)
    error->all(FLERR,"Fix dt fix does not "
               "compute per-grid vector");
  if (usqindex > 0 && modify->fix[n]->size_per_grid_cols == 0)
    error->all(FLERR,"Fix dt fix does not "
               "compute per-grid array");
  if (usqindex > 0 && usqindex > modify->fix[n]->size_per_grid_cols)
    error->all(FLERR,"Fix dt fix array is "
               "accessed out-of-range");

  // vsq----------------------------
  n = strlen(arg[9]);
  id_vsq = new char[n];
  strcpy(id_vsq,&arg[9][2]);

  ptr = strchr(id_vsq,'[');
  if (ptr) {
    if (id_vsq[strlen(id_vsq)-1] != ']')
      error->all(FLERR,"Invalid vsq in fix dt command");
    vsqindex = atoi(ptr+1);
    *ptr = '\0';
  } else vsqindex = 0;

  n = modify->find_fix(id_vsq);
  if (n < 0) error->all(FLERR,"Could not find fix dt vsq fix ID");
  if (modify->fix[n]->per_grid_flag == 0)
    error->all(FLERR,"Fix dt fix does not "
               "compute per-grid info");
  if (vsqindex == 0 && modify->fix[n]->size_per_grid_cols > 0)
    error->all(FLERR,"Fix dt fix does not "
               "compute per-grid vector");
  if (vsqindex > 0 && modify->fix[n]->size_per_grid_cols == 0)
    error->all(FLERR,"Fix dt fix does not "
               "compute per-grid array");
  if (vsqindex > 0 && vsqindex > modify->fix[n]->size_per_grid_cols)
    error->all(FLERR,"Fix dt fix array is "
               "accessed out-of-range");

  // wsq----------------------------
  n = strlen(arg[10]);
  id_wsq = new char[n];
  strcpy(id_wsq,&arg[10][2]);

  ptr = strchr(id_wsq,'[');
  if (ptr) {
    if (id_wsq[strlen(id_wsq)-1] != ']')
      error->all(FLERR,"Invalid wsq in fix dt command");
    wsqindex = atoi(ptr+1);
    *ptr = '\0';
  } else wsqindex = 0;

  n = modify->find_fix(id_wsq);
  if (n < 0) error->all(FLERR,"Could not find fix dt wsq fix ID");
  if (modify->fix[n]->per_grid_flag == 0)
    error->all(FLERR,"Fix dt fix does not "
               "compute per-grid info");
  if (wsqindex == 0 && modify->fix[n]->size_per_grid_cols > 0)
    error->all(FLERR,"Fix dt fix does not "
               "compute per-grid vector");
  if (wsqindex > 0 && modify->fix[n]->size_per_grid_cols == 0)
    error->all(FLERR,"Fix dt fix does not "
               "compute per-grid array");
  if (wsqindex > 0 && wsqindex > modify->fix[n]->size_per_grid_cols)
    error->all(FLERR,"Fix dt fix array is "
               "accessed out-of-range");

  // mode specification
  if (strcmp(arg[11],"warn") == 0)
    mode = WARN;
  else if (strcmp(arg[11],"use_calculated_dt") == 0)
    mode = USE_CALCULATED_DT;
  else
    error->all(FLERR,"Illegal fix dt command: mode argument not recognized");

  // optional argument to store cell dt values
  per_grid_flag = 0;
  int iarg =  narg_required;
  if (narg == iarg+1) {
    if (strcmp(arg[iarg],"store_cell_dt") == 0)
      per_grid_flag = 1;
    else
      error->all(FLERR,"Illegal fix dt command: optional argument not recognized");
  }

  // initialize data structures
  dt_global_calculated = 0.;
  int me = comm->me;
  MPI_Comm_rank(world,&me);
  nglocal = 0;
  lambda = temp = usq = vsq = wsq = NULL;
  vector_grid = NULL;
  per_grid_freq = nevery;
  size_per_grid_cols = 0;


  // find minimum species mass
  Particle::Species *species = particle->species;
  min_species_mass = BIG;
  for (int s = 0; s < particle->nspecies; ++s)
    min_species_mass = MIN(species[s].mass,min_species_mass);

  reallocate();

  // initialize cell timestep vector
  // note:  a value of zero indicates that this fix has not calculated a
  //        timestep for the particular cell.  Such cells are not used in
  //        the global timestep calculation.
  if (per_grid_flag)
    for (int i = 0; i < nglocal; ++i)
      vector_grid[i] = 0.;
}

/* ---------------------------------------------------------------------- */

FixDt::~FixDt()
{
  if (copymode) return;

  delete [] id_lambda;
  delete [] id_temp;
  delete [] id_usq;
  delete [] id_vsq;
  delete [] id_wsq;
  memory->destroy(lambda);
  memory->destroy(temp);
  memory->destroy(usq);
  memory->destroy(vsq);
  memory->destroy(wsq);
  if (per_grid_flag)
    memory->destroy(vector_grid);
}

/* ---------------------------------------------------------------------- */

void FixDt::init()
{
  // setup computes and fixes

  // lambda
  if (lambdawhich == COMPUTE) {
    int icompute = modify->find_compute(id_lambda);
    if (icompute < 0)
      error->all(FLERR,"Could not find fix dt compute ID");
    clambda = modify->compute[icompute];
  } else if (lambdawhich == FIX) {
    int ifix = modify->find_fix(id_lambda);
    if (ifix < 0)
      error->all(FLERR,"Could not find fix dt fix ID");
    flambda = modify->fix[ifix];
  }

  // temperature fix
  int ifix = modify->find_fix(id_temp);
  if (ifix < 0)
    error->all(FLERR,"Could not find temperature fix ID");
  ftemp = modify->fix[ifix];

  // usq fix
  ifix = modify->find_fix(id_usq);
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

int FixDt::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDt::end_of_step()
{
  // check that fixes are computed at compatible times
  if (lambdawhich == FIX && update->ntimestep % flambda->per_grid_freq)
    error->all(FLERR,"Fix dt lambda fix not computed at compatible time");
  if (update->ntimestep % fusq->per_grid_freq)
    error->all(FLERR,"Fix dt usq fix not computed at compatible time");
  if (update->ntimestep % fvsq->per_grid_freq)
    error->all(FLERR,"Fix dt vsq fix not computed at compatible time");
  if (update->ntimestep % fwsq->per_grid_freq)
    error->all(FLERR,"Fix dt wsq fix not computed at compatible time");

  // get lambda from fix or by invoking compute
  if (lambdawhich == COMPUTE) {
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
  } else if (lambdawhich == FIX) {
    if (lambdaindex == 0) {
      memcpy(lambda,flambda->vector_grid,nglocal*sizeof(double));
    } else {
      double **array = flambda->array_grid;
      int index = lambdaindex-1;
      for (int i = 0; i < nglocal; i++)
        lambda[i] = array[i][index];
    }
  }

  // get temperature from fix
  if (tempindex == 0) {
    memcpy(temp,ftemp->vector_grid,nglocal*sizeof(double));
  } else {
    double **array = ftemp->array_grid;
    int index = tempindex-1;
    for (int i = 0; i < nglocal; i++)
      temp[i] = array[i][index];
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

  // compute cell desired timestep ---------------------------------
  Grid::ChildCell *cells = grid->cells;
  Particle::OnePart *particles = particle->particles;

  double dt_sum = 0.;
  double dtmin = BIG;

  bigint ncells_computing_dt = 0;
  for (int i = 0; i < nglocal; ++i) {

    if (per_grid_flag)
      vector_grid[i] = 0.;

    double cell_dt_desired = 0.;

    // cell dt based on mean collision time
    double velocitymag = sqrt(usq[i] + vsq[i] + wsq[i]);
    if (velocitymag > 0.) {
      double mean_collision_time = lambda[i]/velocitymag;
      if (mean_collision_time > 0.)
        cell_dt_desired = collision_fraction*mean_collision_time;
    }

    // cell size
    double dx = cells[i].hi[0] - cells[i].lo[0];
    double dy = cells[i].hi[1] - cells[i].lo[1];
    double dz = 0.;

    // cell dt based on transit time using average velocities
    double dt_candidate;
    double umag = sqrt(usq[i]);
    if (umag > 0.) {
      dt_candidate = transit_fraction*dx/umag;
      cell_dt_desired = MIN(dt_candidate,cell_dt_desired);
    }

    double vmag = sqrt(vsq[i]);
    if (vmag > 0.) {
      dt_candidate = transit_fraction*dy/vmag;
      cell_dt_desired = MIN(dt_candidate,cell_dt_desired);
    }

    if (domain->dimension == 3) {
      dz = cells[i].hi[2] - cells[i].lo[2];
      double wmag = sqrt(wsq[i]);
      if (wmag > 0.) {
        dt_candidate = transit_fraction*dz/wmag;
        cell_dt_desired = MIN(dt_candidate,cell_dt_desired);
      }
    }

    // cell dt based on transit time using maximum most probable speed
    double vrm_max = sqrt(2.0*update->boltz * temp[i] / min_species_mass);
    if (vrm_max > 0.) {
      dt_candidate = transit_fraction*dx/vrm_max;
      cell_dt_desired = MIN(dt_candidate,cell_dt_desired);

      dt_candidate = transit_fraction*dy/vrm_max;
      cell_dt_desired = MIN(dt_candidate,cell_dt_desired);

      if (domain->dimension == 3) {
        dt_candidate = transit_fraction*dz/vrm_max;
        cell_dt_desired = MIN(dt_candidate,cell_dt_desired);
      }
    }

    if (cell_dt_desired > 0.) {
      dtmin = MIN(dtmin, cell_dt_desired);
      dt_sum += cell_dt_desired;
      if (per_grid_flag)
        vector_grid[i] = cell_dt_desired;
      ncells_computing_dt++;
    }
  }

  // set global dtmin
  double dtmin_global;
  MPI_Allreduce(&dtmin,&dtmin_global,1,MPI_DOUBLE,MPI_MIN,world);
  dtmin = dtmin_global;

  if (dtmin < BIG) { // at least one cell computed a timestep

    // set global dtavg
    double cell_sums[2];
    double cell_sums_global[2];
    cell_sums[0] = dt_sum;
    cell_sums[1] = ncells_computing_dt; // implicit conversion of bigint to double to avoid 2 MPI_Allreduce calls
    MPI_Allreduce(cell_sums,cell_sums_global,2,MPI_DOUBLE,MPI_SUM,world);
    double dtavg = cell_sums_global[0]/cell_sums_global[1];

    // set calculated timestep based on user-specified weighting
    dt_global_calculated = (1.-dt_global_weight)*dtmin + dt_global_weight*dtavg;

    if (mode > WARN)
      grid->dt_global = dt_global_calculated;
    else if (mode == WARN && grid->dt_global > dt_global_calculated) {
      if (me == 0) {
        std::cout << std::endl;
        std::cout << "    WARNING: user-set global timestep(=" << grid->dt_global
                  << ") is greater than the calculated global timestep(=" << dt_global_calculated
                  << ")\n\n";
      }
    }
  }
}

/* ----------------------------------------------------------------------
   reallocate arrays if nglocal has changed
   called by constructor and load balancer
------------------------------------------------------------------------- */

void FixDt::reallocate()
{
  if (grid->nlocal == nglocal) return;

  nglocal = grid->nlocal;

  memory->destroy(lambda);
  memory->destroy(temp);
  memory->destroy(usq);
  memory->destroy(vsq);
  memory->destroy(wsq);
  if (per_grid_flag)
    memory->destroy(vector_grid);
  memory->create(lambda,nglocal,"dt:lambda");
  memory->create(temp,nglocal,"dt:temp");
  memory->create(usq,nglocal,"dt:usq");
  memory->create(vsq,nglocal,"dt:vsq");
  memory->create(wsq,nglocal,"dt:wsq");
  if (per_grid_flag)
    memory->create(vector_grid,nglocal,"dt:vector_grid");
}

/* ----------------------------------------------------------------------
   return optimal calculated timestep
------------------------------------------------------------------------- */
double FixDt::compute_scalar()
{
  return dt_global_calculated;
}


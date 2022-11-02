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
#include "fix_dt.h"
#include "update.h"
#include "input.h"
#include "grid.h"
#include "modify.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "particle.h"
#include "mixture.h"
#include "comm.h"
#include "random_mars.h"
#include "random_knuth.h"
#include <iostream>

using namespace SPARTA_NS;

enum{NONE,COMPUTE,FIX};
enum{INT,DOUBLE};

#define INVOKED_PER_GRID 16
#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

FixDt::FixDt(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{

  // RNG
  int me = comm->me;
  random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,me,100);

  // arguments
  if (narg < 12) error->all(FLERR,"Illegal fix dt command");
  scalar_flag = 1;
  global_freq = 1;

  nevery = atoi(arg[2]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix dt command");

  dt_global_weight = atof(arg[3]);
  transit_fraction = atof(arg[4]);
  collision_fraction = atof(arg[5]);

  imix = particle->find_mixture(arg[6]);
  if (imix < 0) error->all(FLERR,"fix dt mixture ID does not exist");

  if (strncmp(arg[7],"c_",2) != 0 && strncmp(arg[8],"f_",2) != 0)
    error->all(FLERR,"Illegal fix dt command");

  if ( !((strncmp(arg[8],"f_",2) == 0) &&
         (strncmp(arg[9],"f_",2) == 0) &&
         (strncmp(arg[10],"f_",2) == 0) &&
         (strncmp(arg[11],"f_",2) == 0)) )
    error->all(FLERR,"Illegal fix dt command");

  id_lambda = NULL;
  id_temp = NULL;
  id_usq = NULL;
  id_vsq = NULL;
  id_wsq = NULL;

  // lambda compute or fix
  int n = strlen(arg[7]);
  id_lambda = new char[n];
  strcpy(id_lambda,&arg[7][2]);
  char *ptr = strchr(id_lambda,'[');
  if (ptr) {
    if (id_lambda[strlen(id_lambda)-1] != ']')
      error->all(FLERR,"Invalid lambda in fix dt command");
    lambdaindex = atoi(ptr+1);
    *ptr = '\0';
  } else lambdaindex = 0;

  if (strncmp(arg[7],"c_",2) == 0) lambdawhich = COMPUTE;
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
  n = strlen(arg[8]);
  id_temp = new char[n];
  strcpy(id_temp,&arg[8][2]);

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
  n = strlen(arg[9]);
  id_usq = new char[n];
  strcpy(id_usq,&arg[9][2]);

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
  n = strlen(arg[10]);
  id_vsq = new char[n];
  strcpy(id_vsq,&arg[10][2]);

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
  n = strlen(arg[11]);
  id_wsq = new char[n];
  strcpy(id_wsq,&arg[11][2]);

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

  // optional args
  //  - prescribed temperature/velocity profiles for initial timestep calculation
  //  - mode specification for timestep utilization options
  tempflag = velflag = 0;

  int iarg = 12;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"temperature") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal fix dt command");
      tempflag = 1;
      tstr = arg[iarg+1];
      if (strcmp(arg[iarg+2],"NULL") == 0) txstr = NULL;
      else txstr = arg[iarg+2];
      if (strcmp(arg[iarg+3],"NULL") == 0) tystr = NULL;
      else tystr = arg[iarg+3];
      if (strcmp(arg[iarg+4],"NULL") == 0) tzstr = NULL;
      else tzstr = arg[iarg+4];
      iarg += 5;
    } else if (strcmp(arg[iarg],"velocity") == 0) {
      if (iarg+7 > narg) error->all(FLERR,"Illegal fix dt command");
      velflag = 1;
      if (strcmp(arg[iarg+1],"NULL") == 0) vxstr = NULL;
      else vxstr = arg[iarg+1];
      if (strcmp(arg[iarg+2],"NULL") == 0) vystr = NULL;
      else vystr = arg[iarg+2];
      if (strcmp(arg[iarg+3],"NULL") == 0) vzstr = NULL;
      else vzstr = arg[iarg+3];
      if (strcmp(arg[iarg+4],"NULL") == 0) vstrx = NULL;
      else vstrx = arg[iarg+4];
      if (strcmp(arg[iarg+5],"NULL") == 0) vstry = NULL;
      else vstry = arg[iarg+5];
      if (strcmp(arg[iarg+6],"NULL") == 0) vstrz = NULL;
      else vstrz = arg[iarg+6];
      iarg += 7;
    } else if (strcmp(arg[iarg],"none") == 0) {
      mode = FIXMODE::NONE;
      iarg += 1;
    } else if (strcmp(arg[iarg],"warn") == 0) {
      mode = FIXMODE::WARN;
      iarg += 1;
    } else if (strcmp(arg[iarg],"use_calculated_global_dt") == 0) {
      mode = FIXMODE::USE_CALCULATED_GLOBAL_DT;
      iarg += 1;
    }
    else if (strcmp(arg[iarg],"use_calculated_cell_dt") == 0) {
      mode = FIXMODE::USE_CALCULATED_CELL_DT;
      grid->use_cell_dt = true;
      iarg += 1;
    }
    else
      error->all(FLERR,"Illegal fix dt command: optional argument not recognized");
  }

  if (tempflag) {
    tvar = input->variable->find(tstr);
    if (tvar < 0)
      error->all(FLERR,"Variable name for fix dt does not exist");
    if (!input->variable->equal_style(tvar))
      error->all(FLERR,"Variable for fix dt is invalid style");
    if (txstr) {
      txvar = input->variable->find(txstr);
      if (txvar < 0)
        error->all(FLERR,"Variable name for fix dt does not exist");
      if (!input->variable->internal_style(txvar))
        error->all(FLERR,"Variable for fix dt is invalid style");
    }
    if (tystr) {
      tyvar = input->variable->find(tystr);
      if (tyvar < 0)
        error->all(FLERR,"Variable name for fix dt does not exist");
      if (!input->variable->internal_style(tyvar))
        error->all(FLERR,"Variable for fix dt is invalid style");
    }
    if (tzstr) {
      tzvar = input->variable->find(tzstr);
      if (tzvar < 0)
        error->all(FLERR,"Variable name for fix dt does not exist");
      if (!input->variable->internal_style(tzvar))
        error->all(FLERR,"Variable for fix dt is invalid style");
    }
  }

  if (velflag) {
    if (vxstr) {
      vxvar = input->variable->find(vxstr);
      if (vxvar < 0)
        error->all(FLERR,"Variable name for fix dt does not exist");
      if (!input->variable->equal_style(vxvar))
        error->all(FLERR,"Variable for fix dt is invalid style");
    }
    if (vystr) {
      vyvar = input->variable->find(vystr);
      if (vyvar < 0)
        error->all(FLERR,"Variable name for fix dt does not exist");
      if (!input->variable->equal_style(vyvar))
        error->all(FLERR,"Variable for fix dt is invalid style");
    }
    if (vzstr) {
      vzvar = input->variable->find(vzstr);
      if (vzvar < 0)
        error->all(FLERR,"Variable name for fix dt does not exist");
      if (!input->variable->equal_style(vzvar))
        error->all(FLERR,"Variable for fix dt is invalid style");
    }
    if (vstrx) {
      vvarx = input->variable->find(vstrx);
      if (vvarx < 0)
        error->all(FLERR,"Variable name for fix dt does not exist");
      if (!input->variable->internal_style(vvarx))
        error->all(FLERR,"Variable for fix dt is invalid style");
    }
    if (vstry) {
      vvary = input->variable->find(vstry);
      if (vvary < 0)
        error->all(FLERR,"Variable name for fix dt does not exist");
      if (!input->variable->internal_style(vvary))
        error->all(FLERR,"Variable for fix dt is invalid style");
    }
    if (vstrz) {
      vvarz = input->variable->find(vstrz);
      if (vvarz < 0)
        error->all(FLERR,"Variable name for fix dt does not exist");
      if (!input->variable->internal_style(vvarz))
        error->all(FLERR,"Variable for fix dt is invalid style");
    }
  }

  // initialize data structures
  MPI_Comm_rank(world,&me);
  nglocal = 0;
  lambda = temp = usq = vsq = wsq = NULL;
  vector_grid = NULL;
  per_grid_flag = 1;
  per_grid_freq = 1;
  size_per_grid_cols = 0;

  reallocate();

  // create per-particle vector for particle time
  particle_time_index = particle->add_custom((char *) "particle_time", DOUBLE, 0);

  // create per-cell array for cell time/desired timestep
  int ghostflag = 1;
  cell_time_index = grid->add_custom((char *) "cell_time", DOUBLE, 2, ghostflag);

  // set up initial timestep
  double x[3];
  Grid::ChildCell *cells = grid->cells;
  double mixture_temperature_thermal = particle->mixture[imix]->temp_thermal;
  double *vstream = particle->mixture[imix]->vstream;
  double vstream_variable[3];
  double v[3];
  double temperature;

  // find minimum species mass
  Particle::Species *species = particle->species;
  min_species_mass = BIG;
  for (int s = 0; s < particle->nspecies; ++s)
    min_species_mass = MIN(species[s].mass,min_species_mass);

  // custom cell time array
  double **cell_time = grid->edarray[grid->ewhich[cell_time_index]];

  for (int i = 0; i < nglocal; ++i) {

#if 0
    x[0] = 0.5*(cells[i].lo[0] + cells[i].hi[0]);
    x[1] = 0.5*(cells[i].lo[1] + cells[i].hi[1]);
    x[2] = 0.5*(cells[i].lo[2] + cells[i].hi[2]);

    if (tempflag)
      temperature = temperature_variable(x);
    else
      temperature = mixture_temperature_thermal;

    if (velflag) {
      velocity_variable(x,vstream,vstream_variable);
      v[0] = vstream_variable[0];
      v[1] = vstream_variable[1];
      v[2] = vstream_variable[2];
    }
    else {
      v[0] = vstream[0];
      v[1] = vstream[1];
      v[2] = vstream[2];
    }
    double vrm_max = sqrt(2.0*update->boltz * temperature / min_species_mass);
#endif

    // for now, set cell timestep to be global value, which was read in

    cell_time[i][0] = grid->time_global;
    double cell_dt_desired = grid->dt_global;
    cell_time[i][1] = cell_dt_desired;
    vector_grid[i] = cell_dt_desired;
  }
}

/* ---------------------------------------------------------------------- */

FixDt::~FixDt()
{
  if (copymode) return;

  delete random;
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
  memory->destroy(vector_grid);
  particle->remove_custom(particle_time_index);
  grid->remove_custom(cell_time_index);
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
  Grid::ChildInfo *cinfo = grid->cinfo;
  Particle::OnePart *particles = particle->particles;

  double dt_sum = 0.;
  double dtmin = BIG;

  // custom cell time array
  double **cell_time = grid->edarray[grid->ewhich[cell_time_index]];

  bigint ncells_with_a_particle = 0.;
  double cell_dt_desired;
  for (int i = 0; i < nglocal; ++i) {

    int np = cinfo[i].count;
    if (np < 1) continue;

    ncells_with_a_particle++;

    // cell dt based on mean collision time
    double mean_collision_time = lambda[i]/sqrt(usq[i] + vsq[i] + wsq[i]);
    if (mean_collision_time > 0.) {
      cell_dt_desired = collision_fraction*mean_collision_time;
    }

    // cell size
    double dx = cells[i].hi[0] - cells[i].lo[0];
    double dy = cells[i].hi[1] - cells[i].lo[1];
    double dz = 0.;

    // cell dt based on transit time using average velocities
    double dt_candidate = transit_fraction*dx/sqrt(usq[i]);
    cell_dt_desired = MIN(dt_candidate,cell_dt_desired);

    dt_candidate = transit_fraction*dy/sqrt(vsq[i]);
    cell_dt_desired = MIN(dt_candidate,cell_dt_desired);

    if (domain->dimension == 3) {
      dz = cells[i].hi[2] - cells[i].lo[2];
      dt_candidate = transit_fraction*dz/sqrt(wsq[i]);
      cell_dt_desired = MIN(dt_candidate,cell_dt_desired);
    }

    // cell dt based on transit time using maximum most probable speed
    double vrm_max = sqrt(2.0*update->boltz * temp[i] / min_species_mass);

    dt_candidate = transit_fraction*dx/vrm_max;
    cell_dt_desired = MIN(dt_candidate,cell_dt_desired);

    dt_candidate = transit_fraction*dy/vrm_max;
    cell_dt_desired = MIN(dt_candidate,cell_dt_desired);

    if (domain->dimension == 3) {
      dt_candidate = transit_fraction*dz/vrm_max;
      cell_dt_desired = MIN(dt_candidate,cell_dt_desired);
    }

    dtmin = MIN(dtmin, cell_dt_desired);
    dt_sum += cell_dt_desired;
    cell_time[i][1] = cell_dt_desired;
    vector_grid[i] = cell_dt_desired;
  }

  // set global dtmin
  double dtmin_global;
  MPI_Allreduce(&dtmin,&dtmin_global,1,MPI_DOUBLE,MPI_MIN,world);
  dtmin = dtmin_global;

  // set global dtavg
  double cell_sums[2];
  double cell_sums_global[2];
  cell_sums[0] = dt_sum;
  cell_sums[1] = ncells_with_a_particle; // implicit conversion of bigint to double to avoid 2 MPI_Allreduce calls
  MPI_Allreduce(cell_sums,cell_sums_global,2,MPI_DOUBLE,MPI_SUM,world);
  double dtavg = cell_sums_global[0]/cell_sums_global[1];

  // set calculated timestep based on user-specified weighting
  dt_global_calculated = (1.-dt_global_weight)*dtmin + dt_global_weight*dtavg;

  if (mode > FIXMODE::WARN)
    grid->dt_global = dt_global_calculated;
  else if (mode == FIXMODE::WARN && grid->dt_global > dt_global_calculated) {
    if (me == 0) {
      std::cout << std::endl;
      std::cout << "    WARNING: user-set global timestep(=" << grid->dt_global
                << ") is greater than the calculated global timestep(=" << dt_global_calculated
                << ")\n\n";
    }
  }
}

/* ----------------------------------------------------------------------
   called when a particle with index is created
   creation used cell desired timestep to set particle time
------------------------------------------------------------------------- */

void FixDt::update_custom(int index, double, double, double, double *)
{
  if (grid->use_cell_dt) {
    double *particle_time = particle->edvec[particle->ewhich[particle_time_index]];
    particle_time[index] = grid->time_global;
  }
  else
    return;
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
  memory->destroy(vector_grid);
  memory->create(lambda,nglocal,"dt:lambda");
  memory->create(temp,nglocal,"dt:temp");
  memory->create(usq,nglocal,"dt:usq");
  memory->create(vsq,nglocal,"dt:vsq");
  memory->create(wsq,nglocal,"dt:wsq");
  memory->create(vector_grid,nglocal,"dt:vector_grid");
}

/* ----------------------------------------------------------------------
   return optimal calculated timestep
------------------------------------------------------------------------- */
double FixDt::compute_scalar()
{
  return dt_global_calculated;
}

/* ----------------------------------------------------------------------
   use grid cell center in tvar variable to generate temperature scale factor
   first plug in cell x,y,z centroid values into txvar,tyvar,tzvar
------------------------------------------------------------------------- */

double FixDt::temperature_variable(double *x)
{

  if (txstr) input->variable->internal_set(txvar,x[0]);
  if (tystr) input->variable->internal_set(tyvar,x[1]);
  if (tzstr) input->variable->internal_set(tzvar,x[2]);

  double temperature = input->variable->compute_equal(tvar);
  return temperature;
}

/* ----------------------------------------------------------------------
   use cell position in vxvar,vyvar,vzvar variables to generate vel stream
   first plug in cell x,y,z centroid values into vvarx,vvary,vvarz
------------------------------------------------------------------------- */

void FixDt::velocity_variable(double *x, double *vstream,
                                    double *vstream_variable)
{
  if (vstrx) input->variable->internal_set(vvarx,x[0]);
  if (vstry) input->variable->internal_set(vvary,x[1]);
  if (vstrz) input->variable->internal_set(vvarz,x[2]);

  if (vxstr) vstream_variable[0] = input->variable->compute_equal(vxvar);
  else vstream_variable[0] = vstream[0];
  if (vystr) vstream_variable[1] = input->variable->compute_equal(vyvar);
  else vstream_variable[1] = vstream[1];
  if (vzstr) vstream_variable[2] = input->variable->compute_equal(vzvar);
  else vstream_variable[2] = vstream[2];
}


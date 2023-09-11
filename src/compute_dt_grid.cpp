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

#include "compute_dt_grid.h"
#include "update.h"
#include "grid.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "memory.h"
#include "error.h"
#include "string.h"

using namespace SPARTA_NS;

enum{NONE,COMPUTE,FIX};

#define INVOKED_PER_GRID 16
#define BIG 1.0e20
#define EPSILON 1.0e-6

/* ---------------------------------------------------------------------- */

ComputeDtGrid::ComputeDtGrid(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg != 10) error->all(FLERR,"Illegal compute dt/grid command");

  int igroup = grid->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Compute dt/grid group ID does not exist");
  groupbit = grid->bitmask[igroup];

  transit_fraction = atof(arg[3]);
  collision_fraction = atof(arg[4]);

  // error checks

  if (transit_fraction < 0.0)
    error->all(FLERR,"Compute dt/grid tfraction must be positive");
  if (collision_fraction < 0.0)
    error->all(FLERR,"Compute dt/grid cfraction must be positive");

  id_lambda = NULL;
  id_temp = NULL;
  id_usq = NULL;
  id_vsq = NULL;
  id_wsq = NULL;

  int n;
  char *ptr;

  // parse lambda as compute or fix

  if (strncmp(arg[5],"c_",2) == 0) lambda_which = COMPUTE;
  else if (strncmp(arg[5],"f_",2) == 0) lambda_which = FIX;
  else error->all(FLERR,"Invalid lambda in compute dt/grid command");

  n = strlen(arg[5]);
  id_lambda = new char[n];
  strcpy(id_lambda,&arg[5][2]);
  ptr = strchr(id_lambda,'[');
  if (ptr) {
    if (id_lambda[strlen(id_lambda)-1] != ']')
      error->all(FLERR,"Invalid lambda in compute dt/grid command");
    lambda_index = atoi(ptr+1);
    *ptr = '\0';
  } else lambda_index = 0;

  if (lambda_which == COMPUTE) {
    int m = modify->find_compute(id_lambda);
    if (m < 0)
      error->all(FLERR,"Could not find compute dt/grid compute ID");
    if (modify->compute[m]->per_grid_flag == 0)
      error->all(FLERR,"Compute dt/grid compute does not "
                 "compute per-grid info");
    if (lambda_index == 0 && modify->compute[m]->size_per_grid_cols > 0)
      error->all(FLERR,"Compute dt/grid compute does not"
                 "compute per-grid vector");
    if (lambda_index > 0 && modify->compute[m]->size_per_grid_cols == 0)
      error->all(FLERR,"Compute dt/grid compute does not "
                 "compute per-grid array");
    if (lambda_index > 0 && lambda_index > modify->compute[m]->size_per_grid_cols)
      error->all(FLERR,"Compute dt/grid compute array is "
                 "accessed out-of-range");
  } else if (lambda_which == FIX) {
    int m = modify->find_fix(id_lambda);
    if (m < 0) error->all(FLERR,"Could not find compute dt/grid fix ID");
    if (modify->fix[m]->per_grid_flag == 0)
      error->all(FLERR,"Compute dt/grid fix does not "
                 "compute per-grid info");
    if (lambda_index == 0 && modify->fix[m]->size_per_grid_cols > 0)
      error->all(FLERR,"Compute dt/grid fix does not "
                 "compute per-grid vector");
    if (lambda_index > 0 && modify->fix[m]->size_per_grid_cols == 0)
      error->all(FLERR,"Compute dt/grid fix does not "
                 "compute per-grid array");
    if (lambda_index > 0 && lambda_index > modify->fix[m]->size_per_grid_cols)
      error->all(FLERR,"Compute dt/grid fix array is "
                 "accessed out-of-range");
  }

  // parse temp as compute or fix

  if (strncmp(arg[6],"c_",2) == 0) temp_which = COMPUTE;
  else if (strncmp(arg[6],"f_",2) == 0) temp_which = FIX;
  else error->all(FLERR,"Invalid temp in compute dt/grid command");

  n = strlen(arg[6]);
  id_temp = new char[n];
  strcpy(id_temp,&arg[6][2]);
  ptr = strchr(id_temp,'[');
  if (ptr) {
    if (id_temp[strlen(id_temp)-1] != ']')
      error->all(FLERR,"Invalid temp in compute dt/grid command");
    temp_index = atoi(ptr+1);
    *ptr = '\0';
  } else temp_index = 0;

  if (temp_which == COMPUTE) {
    int m = modify->find_compute(id_temp);
    if (m < 0)
      error->all(FLERR,"Could not find compute dt/grid compute ID");
    if (modify->compute[m]->per_grid_flag == 0)
      error->all(FLERR,"Compute dt/grid compute does not "
                 "compute per-grid info");
    if (temp_index == 0 && modify->compute[m]->size_per_grid_cols > 0)
      error->all(FLERR,"Compute dt/grid compute does not"
                 "compute per-grid vector");
    if (temp_index > 0 && modify->compute[m]->size_per_grid_cols == 0)
      error->all(FLERR,"Compute dt/grid compute does not "
                 "compute per-grid array");
    if (temp_index > 0 && temp_index > modify->compute[m]->size_per_grid_cols)
      error->all(FLERR,"Compute dt/grid compute array is "
                 "accessed out-of-range");
  } else if (temp_which == FIX) {
    int m = modify->find_fix(id_temp);
    if (m < 0) error->all(FLERR,"Could not find compute dt/grid fix ID");
    if (modify->fix[m]->per_grid_flag == 0)
      error->all(FLERR,"Compute dt/grid fix does not "
                 "compute per-grid info");
    if (temp_index == 0 && modify->fix[m]->size_per_grid_cols > 0)
      error->all(FLERR,"Compute dt/grid fix does not "
                 "compute per-grid vector");
    if (temp_index > 0 && modify->fix[m]->size_per_grid_cols == 0)
      error->all(FLERR,"Compute dt/grid fix does not "
                 "compute per-grid array");
    if (temp_index > 0 && temp_index > modify->fix[m]->size_per_grid_cols)
      error->all(FLERR,"Compute dt/grid fix array is "
                 "accessed out-of-range");
  }

  // parse usq as compute or fix

  if (strncmp(arg[7],"c_",2) == 0) usq_which = COMPUTE;
  else if (strncmp(arg[7],"f_",2) == 0) usq_which = FIX;
  else error->all(FLERR,"Invalid usq in compute dt/grid command");

  n = strlen(arg[7]);
  id_usq = new char[n];
  strcpy(id_usq,&arg[7][2]);
  ptr = strchr(id_usq,'[');
  if (ptr) {
    if (id_usq[strlen(id_usq)-1] != ']')
      error->all(FLERR,"Invalid usq in compute dt/grid command");
    usq_index = atoi(ptr+1);
    *ptr = '\0';
  } else usq_index = 0;

  if (usq_which == COMPUTE) {
    int m = modify->find_compute(id_usq);
    if (m < 0)
      error->all(FLERR,"Could not find compute dt/grid compute ID");
    if (modify->compute[m]->per_grid_flag == 0)
      error->all(FLERR,"Compute dt/grid compute does not "
                 "compute per-grid info");
    if (usq_index == 0 && modify->compute[m]->size_per_grid_cols > 0)
      error->all(FLERR,"Compute dt/grid compute does not"
                 "compute per-grid vector");
    if (usq_index > 0 && modify->compute[m]->size_per_grid_cols == 0)
      error->all(FLERR,"Compute dt/grid compute does not "
                 "compute per-grid array");
    if (usq_index > 0 && usq_index > modify->compute[m]->size_per_grid_cols)
      error->all(FLERR,"Compute dt/grid compute array is "
                 "accessed out-of-range");
  } else if (usq_which == FIX) {
    int m = modify->find_fix(id_usq);
    if (m < 0) error->all(FLERR,"Could not find compute dt/grid fix ID");
    if (modify->fix[m]->per_grid_flag == 0)
      error->all(FLERR,"Compute dt/grid fix does not "
                 "compute per-grid info");
    if (usq_index == 0 && modify->fix[m]->size_per_grid_cols > 0)
      error->all(FLERR,"Compute dt/grid fix does not "
                 "compute per-grid vector");
    if (usq_index > 0 && modify->fix[m]->size_per_grid_cols == 0)
      error->all(FLERR,"Compute dt/grid fix does not "
                 "compute per-grid array");
    if (usq_index > 0 && usq_index > modify->fix[m]->size_per_grid_cols)
      error->all(FLERR,"Compute dt/grid fix array is "
                 "accessed out-of-range");
  }

  // parse vsq as compute or fix

  if (strncmp(arg[8],"c_",2) == 0) vsq_which = COMPUTE;
  else if (strncmp(arg[8],"f_",2) == 0) vsq_which = FIX;
  else error->all(FLERR,"Invalid vsq in compute dt/grid command");

  n = strlen(arg[8]);
  id_vsq = new char[n];
  strcpy(id_vsq,&arg[8][2]);
  ptr = strchr(id_vsq,'[');
  if (ptr) {
    if (id_vsq[strlen(id_vsq)-1] != ']')
      error->all(FLERR,"Invalid vsq in compute dt/grid command");
    vsq_index = atoi(ptr+1);
    *ptr = '\0';
  } else vsq_index = 0;

  if (vsq_which == COMPUTE) {
    int m = modify->find_compute(id_vsq);
    if (m < 0)
      error->all(FLERR,"Could not find compute dt/grid compute ID");
    if (modify->compute[m]->per_grid_flag == 0)
      error->all(FLERR,"Compute dt/grid compute does not "
                 "compute per-grid info");
    if (vsq_index == 0 && modify->compute[m]->size_per_grid_cols > 0)
      error->all(FLERR,"Compute dt/grid compute does not"
                 "compute per-grid vector");
    if (vsq_index > 0 && modify->compute[m]->size_per_grid_cols == 0)
      error->all(FLERR,"Compute dt/grid compute does not "
                 "compute per-grid array");
    if (vsq_index > 0 && vsq_index > modify->compute[m]->size_per_grid_cols)
      error->all(FLERR,"Compute dt/grid compute array is "
                 "accessed out-of-range");
  } else if (vsq_which == FIX) {
    int m = modify->find_fix(id_vsq);
    if (m < 0) error->all(FLERR,"Could not find compute dt/grid fix ID");
    if (modify->fix[m]->per_grid_flag == 0)
      error->all(FLERR,"Compute dt/grid fix does not "
                 "compute per-grid info");
    if (vsq_index == 0 && modify->fix[m]->size_per_grid_cols > 0)
      error->all(FLERR,"Compute dt/grid fix does not "
                 "compute per-grid vector");
    if (vsq_index > 0 && modify->fix[m]->size_per_grid_cols == 0)
      error->all(FLERR,"Compute dt/grid fix does not "
                 "compute per-grid array");
    if (vsq_index > 0 && vsq_index > modify->fix[m]->size_per_grid_cols)
      error->all(FLERR,"Compute dt/grid fix array is "
                 "accessed out-of-range");
  }

  // parse wsq as compute or fix

  if (strncmp(arg[9],"c_",2) == 0) wsq_which = COMPUTE;
  else if (strncmp(arg[9],"f_",2) == 0) wsq_which = FIX;
  else error->all(FLERR,"Invalid wsq in compute dt/grid command");

  n = strlen(arg[9]);
  id_wsq = new char[n];
  strcpy(id_wsq,&arg[9][2]);
  ptr = strchr(id_wsq,'[');
  if (ptr) {
    if (id_wsq[strlen(id_wsq)-1] != ']')
      error->all(FLERR,"Invalid wsq in compute dt/grid command");
    wsq_index = atoi(ptr+1);
    *ptr = '\0';
  } else wsq_index = 0;

  if (wsq_which == COMPUTE) {
    int m = modify->find_compute(id_wsq);
    if (m < 0)
      error->all(FLERR,"Could not find compute dt/grid compute ID");
    if (modify->compute[m]->per_grid_flag == 0)
      error->all(FLERR,"Compute dt/grid compute does not "
                 "compute per-grid info");
    if (wsq_index == 0 && modify->compute[m]->size_per_grid_cols > 0)
      error->all(FLERR,"Compute dt/grid compute does not"
                 "compute per-grid vector");
    if (wsq_index > 0 && modify->compute[m]->size_per_grid_cols == 0)
      error->all(FLERR,"Compute dt/grid compute does not "
                 "compute per-grid array");
    if (wsq_index > 0 && wsq_index > modify->compute[m]->size_per_grid_cols)
      error->all(FLERR,"Compute dt/grid compute array is "
                 "accessed out-of-range");
  } else if (wsq_which == FIX) {
    int m = modify->find_fix(id_wsq);
    if (m < 0) error->all(FLERR,"Could not find compute dt/grid fix ID");
    if (modify->fix[m]->per_grid_flag == 0)
      error->all(FLERR,"Compute dt/grid fix does not "
                 "compute per-grid info");
    if (wsq_index == 0 && modify->fix[m]->size_per_grid_cols > 0)
      error->all(FLERR,"Compute dt/grid fix does not "
                 "compute per-grid vector");
    if (wsq_index > 0 && modify->fix[m]->size_per_grid_cols == 0)
      error->all(FLERR,"Compute dt/grid fix does not "
                 "compute per-grid array");
    if (wsq_index > 0 && wsq_index > modify->fix[m]->size_per_grid_cols)
      error->all(FLERR,"Compute dt/grid fix array is "
                 "accessed out-of-range");
  }

  // find minimum species mass

  Particle::Species *species = particle->species;
  min_species_mass = BIG;
  for (int s = 0; s < particle->nspecies; ++s)
    min_species_mass = MIN(species[s].mass,min_species_mass);

  // initialize data structures

  per_grid_flag = 1;
  size_per_grid_cols = 0;

  nglocal = 0;
  vector_grid = NULL;
  lambda = temp = usq = vsq = wsq = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeDtGrid::~ComputeDtGrid()
{
  if (copymode) return;

  delete [] id_lambda;
  delete [] id_temp;
  delete [] id_usq;
  delete [] id_vsq;
  delete [] id_wsq;

  memory->destroy(vector_grid);
  memory->destroy(lambda);
  memory->destroy(temp);
  memory->destroy(usq);
  memory->destroy(vsq);
  memory->destroy(wsq);
}

/* ---------------------------------------------------------------------- */

void ComputeDtGrid::init()
{
  reallocate();

  // clambda/flambda = compute or fix that calculates per grid cell lambda
  // ditto for temp,usq,vsq,wsq

  if (lambda_which == COMPUTE) {
    int icompute = modify->find_compute(id_lambda);
    if (icompute < 0)
      error->all(FLERR,"Could not find compute dt/grid compute ID");
    clambda = modify->compute[icompute];
  } else if (lambda_which == FIX) {
    int ifix = modify->find_fix(id_lambda);
    if (ifix < 0)
      error->all(FLERR,"Could not find compute dt/grid fix ID");
    flambda = modify->fix[ifix];
  }

  if (temp_which == COMPUTE) {
    int icompute = modify->find_compute(id_temp);
    if (icompute < 0)
      error->all(FLERR,"Could not find compute dt/grid compute ID");
    ctemp = modify->compute[icompute];
  } else if (temp_which == FIX) {
    int ifix = modify->find_fix(id_temp);
    if (ifix < 0)
      error->all(FLERR,"Could not find compute dt/grid fix ID");
    ftemp = modify->fix[ifix];
  }

  if (usq_which == COMPUTE) {
    int icompute = modify->find_compute(id_usq);
    if (icompute < 0)
      error->all(FLERR,"Could not find compute dt/grid compute ID");
    cusq = modify->compute[icompute];
  } else if (usq_which == FIX) {
    int ifix = modify->find_fix(id_usq);
    if (ifix < 0)
      error->all(FLERR,"Could not find compute dt/grid fix ID");
    fusq = modify->fix[ifix];
  }

  if (vsq_which == COMPUTE) {
    int icompute = modify->find_compute(id_vsq);
    if (icompute < 0)
      error->all(FLERR,"Could not find compute dt/grid compute ID");
    cvsq = modify->compute[icompute];
  } else if (vsq_which == FIX) {
    int ifix = modify->find_fix(id_vsq);
    if (ifix < 0)
      error->all(FLERR,"Could not find compute dt/grid fix ID");
    fvsq = modify->fix[ifix];
  }

  if (wsq_which == COMPUTE) {
    int icompute = modify->find_compute(id_wsq);
    if (icompute < 0)
      error->all(FLERR,"Could not find compute dt/grid compute ID");
    cwsq = modify->compute[icompute];
  } else if (wsq_which == FIX) {
    int ifix = modify->find_fix(id_wsq);
    if (ifix < 0)
      error->all(FLERR,"Could not find compute dt/grid fix ID");
    fwsq = modify->fix[ifix];
  }
}

/* ---------------------------------------------------------------------- */

void ComputeDtGrid::compute_per_grid()
{
  invoked_per_grid = update->ntimestep;

  // check that fixes are computed at compatible times

  if (lambda_which == FIX && update->ntimestep % flambda->per_grid_freq)
    error->all(FLERR,"Fix dt lambda fix not computed at compatible time");
  if (temp_which == FIX && update->ntimestep % ftemp->per_grid_freq)
    error->all(FLERR,"Fix dt temp fix not computed at compatible time");
  if (usq_which == FIX && update->ntimestep % fusq->per_grid_freq)
    error->all(FLERR,"Fix dt usq fix not computed at compatible time");
  if (vsq_which == FIX && update->ntimestep % fvsq->per_grid_freq)
    error->all(FLERR,"Fix dt vsq fix not computed at compatible time");
  if (wsq_which == FIX && update->ntimestep % fwsq->per_grid_freq)
    error->all(FLERR,"Fix dt wsq fix not computed at compatible time");

  // grab per grid cell lambda from compute or fix, invoke compute if needed
  // ditto for temp,usq,vsq,wsq

  if (lambda_which == COMPUTE) {
    if (!(clambda->invoked_flag & INVOKED_PER_GRID)) {
      clambda->compute_per_grid();
      clambda->invoked_flag |= INVOKED_PER_GRID;
    }
    if (clambda->post_process_grid_flag)
      clambda->post_process_grid(lambda_index,1,NULL,NULL,NULL,1);

    if (lambda_index == 0)
      memcpy(lambda,clambda->vector_grid,nglocal*sizeof(double));
    else {
      int index = lambda_index-1;
      double **array = clambda->array_grid;
      for (int i = 0; i < nglocal; i++)
        lambda[i] = array[i][index];
    }
  } else if (lambda_which == FIX) {
    if (lambda_index == 0) {
      memcpy(lambda,flambda->vector_grid,nglocal*sizeof(double));
    } else {
      double **array = flambda->array_grid;
      int index = lambda_index-1;
      for (int i = 0; i < nglocal; i++)
        lambda[i] = array[i][index];
    }
  }

  if (temp_which == COMPUTE) {
    if (!(ctemp->invoked_flag & INVOKED_PER_GRID)) {
      ctemp->compute_per_grid();
      ctemp->invoked_flag |= INVOKED_PER_GRID;
    }
    if (ctemp->post_process_grid_flag)
      ctemp->post_process_grid(temp_index,1,NULL,NULL,NULL,1);

    if (temp_index == 0)
      memcpy(temp,ctemp->vector_grid,nglocal*sizeof(double));
    else {
      int index = temp_index-1;
      double **array = ctemp->array_grid;
      for (int i = 0; i < nglocal; i++)
        temp[i] = array[i][index];
    }
  } else if (temp_which == FIX) {
    if (temp_index == 0) {
      memcpy(temp,ftemp->vector_grid,nglocal*sizeof(double));
    } else {
      double **array = ftemp->array_grid;
      int index = temp_index-1;
      for (int i = 0; i < nglocal; i++)
        temp[i] = array[i][index];
    }
  }

  if (usq_which == COMPUTE) {
    if (!(cusq->invoked_flag & INVOKED_PER_GRID)) {
      cusq->compute_per_grid();
      cusq->invoked_flag |= INVOKED_PER_GRID;
    }
    if (cusq->post_process_grid_flag)
      cusq->post_process_grid(usq_index,1,NULL,NULL,NULL,1);

    if (usq_index == 0)
      memcpy(usq,cusq->vector_grid,nglocal*sizeof(double));
    else {
      int index = usq_index-1;
      double **array = cusq->array_grid;
      for (int i = 0; i < nglocal; i++)
        usq[i] = array[i][index];
    }
  } else if (usq_which == FIX) {
    if (usq_index == 0) {
      memcpy(usq,fusq->vector_grid,nglocal*sizeof(double));
    } else {
      double **array = fusq->array_grid;
      int index = usq_index-1;
      for (int i = 0; i < nglocal; i++)
        usq[i] = array[i][index];
    }
  }

  if (vsq_which == COMPUTE) {
    if (!(cvsq->invoked_flag & INVOKED_PER_GRID)) {
      cvsq->compute_per_grid();
      cvsq->invoked_flag |= INVOKED_PER_GRID;
    }
    if (cvsq->post_process_grid_flag)
      cvsq->post_process_grid(vsq_index,1,NULL,NULL,NULL,1);

    if (vsq_index == 0)
      memcpy(vsq,cvsq->vector_grid,nglocal*sizeof(double));
    else {
      int index = vsq_index-1;
      double **array = cvsq->array_grid;
      for (int i = 0; i < nglocal; i++)
        vsq[i] = array[i][index];
    }
  } else if (vsq_which == FIX) {
    if (vsq_index == 0) {
      memcpy(vsq,fvsq->vector_grid,nglocal*sizeof(double));
    } else {
      double **array = fvsq->array_grid;
      int index = vsq_index-1;
      for (int i = 0; i < nglocal; i++)
        vsq[i] = array[i][index];
    }
  }

  if (wsq_which == COMPUTE) {
    if (!(cwsq->invoked_flag & INVOKED_PER_GRID)) {
      cwsq->compute_per_grid();
      cwsq->invoked_flag |= INVOKED_PER_GRID;
    }
    if (cwsq->post_process_grid_flag)
      cwsq->post_process_grid(wsq_index,1,NULL,NULL,NULL,1);

    if (wsq_index == 0)
      memcpy(wsq,cwsq->vector_grid,nglocal*sizeof(double));
    else {
      int index = wsq_index-1;
      double **array = cwsq->array_grid;
      for (int i = 0; i < nglocal; i++)
        wsq[i] = array[i][index];
    }
  } else if (wsq_which == FIX) {
    if (wsq_index == 0) {
      memcpy(wsq,fwsq->vector_grid,nglocal*sizeof(double));
    } else {
      double **array = fwsq->array_grid;
      int index = wsq_index-1;
      for (int i = 0; i < nglocal; i++)
        wsq[i] = array[i][index];
    }
  }

  // calculate per grid cell timestep for cells in group
  // set timestep = 0.0 for cells not in group

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  Particle::OnePart *particles = particle->particles;
  int dimension = domain->dimension;

  double dx,dy,dz;
  double umag,vmag,wmag,velmag2;
  double vrm_max;
  double cell_dt_desired,mean_collision_time;
  double dt_candidate;

  for (int i = 0; i < nglocal; ++i) {
    vector_grid[i] = 0.0;

    // exclude cells not in the specified group
    if (!(cinfo[i].mask & groupbit)) continue;

    // exclude cells that have no particles
    //  (includes split cells and unsplit cells interior to surface objects)
    if (cinfo[i].count == 0) continue;

    // check sufficiency of cell data to calculate cell dt
    vrm_max = sqrt(2.0*update->boltz * temp[i] / min_species_mass);
    velmag2 = usq[i] + vsq[i] + wsq[i];
    if ( !((vrm_max > EPSILON) && (lambda[i] > EPSILON) && (velmag2 > EPSILON)) ) continue;

    // cell dt based on mean collision time
    mean_collision_time = lambda[i]/vrm_max;
    cell_dt_desired = collision_fraction*mean_collision_time;

    // cell size = dx,dy,dz
    dx = cells[i].hi[0] - cells[i].lo[0];
    dy = cells[i].hi[1] - cells[i].lo[1];
    dz = cells[i].hi[2] - cells[i].lo[2];

    // cell dt based on transit time using average velocities
    umag = sqrt(usq[i]);
    if (umag > 0.0) {
      dt_candidate = transit_fraction*dx/umag;
      cell_dt_desired = MIN(dt_candidate,cell_dt_desired);
    }
    vmag = sqrt(vsq[i]);
    if (vmag > 0.0) {
      dt_candidate = transit_fraction*dy/vmag;
      cell_dt_desired = MIN(dt_candidate,cell_dt_desired);
    }
    if (domain->dimension == 3) {
      wmag = sqrt(wsq[i]);
      if (wmag > 0.0) {
        dt_candidate = transit_fraction*dz/wmag;
        cell_dt_desired = MIN(dt_candidate,cell_dt_desired);
      }
    }

    // cell dt based on transit time using maximum most probable speed
    dt_candidate = transit_fraction*dx/vrm_max;
    cell_dt_desired = MIN(dt_candidate,cell_dt_desired);
    dt_candidate = transit_fraction*dy/vrm_max;
    cell_dt_desired = MIN(dt_candidate,cell_dt_desired);
    if (domain->dimension == 3) {
      dt_candidate = transit_fraction*dz/vrm_max;
      cell_dt_desired = MIN(dt_candidate,cell_dt_desired);
    }

    // per grid cell timestep = final cell_dt_desired for all criteria
    vector_grid[i] = cell_dt_desired;
  }
}

/* ----------------------------------------------------------------------
   reallocate arrays if nglocal has changed
   called by init() and load balancer
------------------------------------------------------------------------- */

void ComputeDtGrid::reallocate()
{
  if (grid->nlocal == nglocal) return;

  nglocal = grid->nlocal;

  memory->destroy(vector_grid);
  memory->destroy(lambda);
  memory->destroy(temp);
  memory->destroy(usq);
  memory->destroy(vsq);
  memory->destroy(wsq);

  memory->create(vector_grid,nglocal,"dt/grid:vector_grid");
  memory->create(lambda,nglocal,"dt/grid:lambda");
  memory->create(temp,nglocal,"dt/grid:temp");
  memory->create(usq,nglocal,"dt/grid:usq");
  memory->create(vsq,nglocal,"dt/grid:vsq");
  memory->create(wsq,nglocal,"dt/grid:wsq");
}

/* ----------------------------------------------------------------------
   memory usage of local grid-based array
------------------------------------------------------------------------- */

bigint ComputeDtGrid::memory_usage()
{
  bigint bytes;
  bytes = nglocal * sizeof(double);      // vector_grid
  bytes += 5*nglocal * sizeof(double);   // lambda,temp,usq,vsq,wsq
  return bytes;
}


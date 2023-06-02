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
#include "compute_lambda_grid.h"
#include "update.h"
#include "grid.h"
#include "domain.h"
#include "collide.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{NONE,COMPUTE,FIX};
enum{KNONE,KALL,KX,KY,KZ};

#define INVOKED_PER_GRID 16
#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

ComputeLambdaGrid::ComputeLambdaGrid(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg < 5 || narg > 6)
    error->all(FLERR,"Illegal compute lambda/grid command");

  // parse three required input fields
  // customize a new keyword by adding to if statement

  id_nrho = id_temp = NULL;

  if (strncmp(arg[2],"c_",2) == 0 || strncmp(arg[2],"f_",2) == 0) {
    int n = strlen(arg[2]);
    id_nrho = new char[n];
    strcpy(id_nrho,&arg[2][2]);

    char *ptr = strchr(id_nrho,'[');
    if (ptr) {
      if (id_nrho[strlen(id_nrho)-1] != ']')
        error->all(FLERR,"Invalid nrho in compute lambda/grid command");
      nrhoindex = atoi(ptr+1);
      *ptr = '\0';
    } else nrhoindex = 0;

    if (strncmp(arg[2],"c_",2) == 0) nrhowhich = COMPUTE;
    else nrhowhich = FIX;

    if (nrhowhich == COMPUTE) {
      int n = modify->find_compute(id_nrho);
      if (n < 0)
        error->all(FLERR,"Could not find compute lambda/grid compute ID");
      if (modify->compute[n]->per_grid_flag == 0)
        error->all(FLERR,"Compute lambda/grid compute does not "
                   "compute per-grid info");
      if (nrhoindex == 0 && modify->compute[n]->size_per_grid_cols > 0)
        error->all(FLERR,
                   "Compute lambda/grid compute does not "
                   "compute per-grid vector");
      if (nrhoindex > 0 && modify->compute[n]->size_per_grid_cols == 0)
        error->all(FLERR,
                   "Compute lambda/grid compute does not "
                   "compute per-grid array");
      if (nrhoindex > 0 && nrhoindex > modify->compute[n]->size_per_grid_cols)
        error->all(FLERR,"Compute lambda compute vector is "
                   "accessed out-of-range");
    } else {
      int n = modify->find_fix(id_nrho);
      if (n < 0) error->all(FLERR,"Could not find compute lambda/grid fix ID");
      if (modify->fix[n]->per_grid_flag == 0)
        error->all(FLERR,"Compute lambda/grid fix does not "
                   "compute per-grid info");
      if (nrhoindex == 0 && modify->fix[n]->size_per_grid_cols > 0)
        error->all(FLERR,"Compute lambda/grid fix does not "
                   "compute per-grid vector");
      if (nrhoindex > 0 && modify->fix[n]->size_per_grid_cols == 0)
        error->all(FLERR,"Compute lambda/grid fix does not "
                   "compute per-grid array");
      if (nrhoindex > 0 && nrhoindex > modify->fix[n]->size_per_grid_cols)
        error->all(FLERR,"Compute lambda/grid fix array is "
                   "accessed out-of-range");
    }
  } else error->all(FLERR,"Illegal compute lambda/grid command");

  if (strncmp(arg[3],"c_",2) == 0 || strncmp(arg[3],"f_",2) == 0) {
    int n = strlen(arg[3]);
    id_temp = new char[n];
    strcpy(id_temp,&arg[3][2]);

    char *ptr = strchr(id_temp,'[');
    if (ptr) {
      if (id_temp[strlen(id_temp)-1] != ']')
        error->all(FLERR,"Invalid temp in compute lambda/grid command");
      tempindex = atoi(ptr+1);
      *ptr = '\0';
    } else tempindex = 0;

    if (strncmp(arg[3],"c_",2) == 0) tempwhich = COMPUTE;
    else tempwhich = FIX;

    if (tempwhich == COMPUTE) {
      int n = modify->find_compute(id_temp);
      if (n < 0)
        error->all(FLERR,"Could not find compute lambda/grid compute ID");
      if (modify->compute[n]->per_grid_flag == 0)
        error->all(FLERR,"Compute lambda/grid compute does not "
                   "compute per-grid info");
      if (tempindex == 0 && modify->compute[n]->size_per_grid_cols > 0)
        error->all(FLERR,
                   "Compute lambda/grid compute does not "
                   "compute per-grid vector");
      if (tempindex > 0 && modify->compute[n]->size_per_grid_cols == 0)
        error->all(FLERR,
                   "Compute lambda/grid compute does not "
                   "compute per-grid array");
      if (tempindex > 0 && tempindex > modify->compute[n]->size_per_grid_cols)
        error->all(FLERR,"Compute lambda compute vector is "
                   "accessed out-of-range");
    } else {
      int n = modify->find_fix(id_temp);
      if (n < 0) error->all(FLERR,"Could not find compute lambda/grid fix ID");
      if (modify->fix[n]->per_grid_flag == 0)
        error->all(FLERR,"Compute lambda/grid fix does not "
                   "compute per-grid info");
      if (tempindex == 0 && modify->fix[n]->size_per_grid_cols > 0)
        error->all(FLERR,"Compute lambda/grid fix does not "
                   "compute per-grid vector");
      if (tempindex > 0 && modify->fix[n]->size_per_grid_cols == 0)
        error->all(FLERR,"Compute lambda/grid fix does not "
                   "compute per-grid array");
      if (tempindex > 0 && tempindex > modify->fix[n]->size_per_grid_cols)
        error->all(FLERR,"Compute lambda/grid fix array is "
                   "accessed out-of-range");
    }
  } else if (strcmp(arg[3],"NULL") == 0) tempwhich = NONE;
  else error->all(FLERR,"Illegal compute lambda/grid command");

  int n = strlen(arg[4]) + 1;
  species = new char[n];
  strcpy(species,arg[4]);

  // optional arg

  kflag = KNONE;
  if (narg == 6) {
    if (strcmp(arg[5],"kall") == 0) kflag = KALL;
    else if (strcmp(arg[5],"kx") == 0) kflag = KX;
    else if (strcmp(arg[5],"ky") == 0) kflag = KY;
    else if (strcmp(arg[5],"kz") == 0) kflag = KZ;
    else error->all(FLERR,"Illegal compute lambda/grid command");
  }

  if (kflag == KZ && domain->dimension == 2)
    error->all(FLERR,"Cannot use compute lambda/grid kz for 2d simulation");

  // initialize data structures

  per_grid_flag = 1;
  if (kflag == KNONE) size_per_grid_cols = 0;
  else size_per_grid_cols = 2;

  nglocal = 0;
  vector_grid = NULL;
  array_grid = NULL;
  nrho = temp = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeLambdaGrid::~ComputeLambdaGrid()
{
  if (copymode) return;

  delete [] id_nrho;
  delete [] id_temp;
  delete [] species;
  memory->destroy(vector_grid);
  memory->destroy(array_grid);
  memory->destroy(nrho);
  memory->destroy(temp);
}

/* ---------------------------------------------------------------------- */

void ComputeLambdaGrid::init()
{
  reallocate();

  // setup computes and fixes

  if (nrhowhich == COMPUTE) {
    int icompute = modify->find_compute(id_nrho);
    if (icompute < 0)
      error->all(FLERR,"Could not find compute lambda/grid compute ID");
    cnrho = modify->compute[icompute];
  } else if (nrhowhich == FIX) {
    int ifix = modify->find_fix(id_nrho);
    if (ifix < 0)
      error->all(FLERR,"Could not find compute lambda/grid fix ID");
    fnrho = modify->fix[ifix];
  }

  if (tempwhich == COMPUTE) {
    int icompute = modify->find_compute(id_temp);
    if (icompute < 0)
      error->all(FLERR,"Could not find compute lambda/grid compute ID");
    ctemp = modify->compute[icompute];
  } else if (tempwhich == FIX) {
    int ifix = modify->find_fix(id_temp);
    if (ifix < 0)
      error->all(FLERR,"Could not find compute lambda/grid fix ID");
    ftemp = modify->fix[ifix];
  }

  // setup collision parameter

  if (!collide)
    error->all(FLERR,"Compute lambda/grid requires a "
               "collision style be defined");

  int ispecies = particle->find_species(species);
  if (ispecies < 0)
    error->all(FLERR,"Compute lambda/grid species is not defined");

  dref = collide->extract(ispecies,ispecies,"diam");
  tref = collide->extract(ispecies,ispecies,"tref");
  omega = collide->extract(ispecies,ispecies,"omega");
  prefactor = sqrt(2.0) * MY_PI * dref*dref;
}

/* ---------------------------------------------------------------------- */

void ComputeLambdaGrid::compute_per_grid()
{
  invoked_per_grid = update->ntimestep;

  if (nrhowhich == FIX && update->ntimestep % fnrho->per_grid_freq)
    error->all(FLERR,"Compute lambda/grid fix not computed at compatible time");
  if (tempwhich == FIX && update->ntimestep % ftemp->per_grid_freq)
    error->all(FLERR,"Compute lambda/grid fix not computed at compatible time");

  // grab nrho and temp values from compute or fix
  // invoke nrho and temp computes as needed

  if (nrhowhich == COMPUTE) {
    if (!(cnrho->invoked_flag & INVOKED_PER_GRID)) {
      cnrho->compute_per_grid();
      cnrho->invoked_flag |= INVOKED_PER_GRID;
    }

    if (cnrho->post_process_grid_flag)
      cnrho->post_process_grid(nrhoindex,1,NULL,NULL,NULL,1);

    if (nrhoindex == 0 || cnrho->post_process_grid_flag)
      memcpy(nrho,cnrho->vector_grid,nglocal*sizeof(double));
    else {
      int index = nrhoindex-1;
      double **array = cnrho->array_grid;
      for (int i = 0; i < nglocal; i++)
        nrho[i] = array[i][index];
    }

  } else if (nrhowhich == FIX) {
    if (nrhoindex == 0) {
      memcpy(nrho,fnrho->vector_grid,nglocal*sizeof(double));
    } else {
      double **array = fnrho->array_grid;
      int index = nrhoindex-1;
      for (int i = 0; i < nglocal; i++)
        nrho[i] = array[i][index];
    }
  }

  if (tempwhich == COMPUTE) {
    if (!(ctemp->invoked_flag & INVOKED_PER_GRID)) {
      ctemp->compute_per_grid();
      ctemp->invoked_flag |= INVOKED_PER_GRID;
    }

    if (ctemp->post_process_grid_flag)
      ctemp->post_process_grid(tempindex,1,NULL,NULL,NULL,1);

    if (tempindex == 0 || ctemp->post_process_grid_flag)
      memcpy(temp,ctemp->vector_grid,nglocal*sizeof(double));
    else {
      int index = tempindex-1;
      double **array = ctemp->array_grid;
      for (int i = 0; i < nglocal; i++)
        temp[i] = array[i][index];
    }

  } else if (tempwhich == FIX) {
    if (tempindex == 0)
      memcpy(temp,ftemp->vector_grid,nglocal*sizeof(double));
    else {
      double **array = ftemp->array_grid;
      int index = tempindex-1;
      for (int i = 0; i < nglocal; i++)
        temp[i] = array[i][index];
    }
  }

  // compute mean free path for each grid cell
  // formula from Bird, eq 4.65

  double lambda;

  for (int i = 0; i < nglocal; i++) {
    if (nrho[i] == 0.0) lambda = BIG;
    else if (tempwhich == NONE || temp[i] == 0.0)
      lambda = 1.0 / (prefactor * nrho[i]);
    else
      lambda = 1.0 / (prefactor * nrho[i] * pow(tref/temp[i],omega-0.5));

    if (kflag == KNONE) vector_grid[i] = lambda;
    else array_grid[i][0] = lambda;
  }

  // calculate per-cell Knudsen number

  if (kflag == KNONE) return;

  Grid::ChildCell *cells = grid->cells;
  int dimension = domain->dimension;

  if (kflag == KALL) {
    double size;
    for (int i = 0; i < nglocal; i++) {
      size = (cells[i].hi[0] - cells[i].lo[0]);
      size += (cells[i].hi[1] - cells[i].lo[1]);
      if (dimension == 2) size *= 0.5;
      else {
        size += (cells[i].hi[2] - cells[i].lo[2]);
        size /= 3.0;
      }
      array_grid[i][1] = array_grid[i][0] / size;
    }
  } else if (kflag == KX) {
    for (int i = 0; i < nglocal; i++)
      array_grid[i][1] = array_grid[i][0] / (cells[i].hi[0] - cells[i].lo[0]);
  } else if (kflag == KY) {
    for (int i = 0; i < nglocal; i++)
      array_grid[i][1] = array_grid[i][0] / (cells[i].hi[1] - cells[i].lo[1]);
  } else if (kflag == KZ) {
    for (int i = 0; i < nglocal; i++)
      array_grid[i][1] = array_grid[i][0] / (cells[i].hi[2] - cells[i].lo[2]);
  }
}

/* ----------------------------------------------------------------------
   reallocate arrays if nglocal has changed
   called by init() and load balancer
------------------------------------------------------------------------- */

void ComputeLambdaGrid::reallocate()
{
  if (grid->nlocal == nglocal) return;

  nglocal = grid->nlocal;
  if (kflag == KNONE) {
    memory->destroy(vector_grid);
    memory->create(vector_grid,nglocal,"lambda/grid:vector_grid");
  } else {
    memory->destroy(array_grid);
    memory->create(array_grid,nglocal,2,"lambda/grid:array_grid");
  }

  memory->destroy(nrho);
  memory->create(nrho,nglocal,"lambda/grid:nrho");
  if (tempwhich != NONE) {
    memory->destroy(temp);
    memory->create(temp,nglocal,"lambda/grid:temp");
  }
}

/* ----------------------------------------------------------------------
   memory usage of local grid-based array
------------------------------------------------------------------------- */

bigint ComputeLambdaGrid::memory_usage()
{
  bigint bytes;
  bytes = nglocal * sizeof(double);
  if (kflag != KNONE) bytes += nglocal * sizeof(double);
  bytes += nglocal * sizeof(double);                            // nrho
  if (tempwhich != KNONE) bytes += nglocal * sizeof(double);    // temp
  return bytes;
}

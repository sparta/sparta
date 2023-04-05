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
#include "fix_field_grid.h"
#include "grid.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

FixFieldGrid::FixFieldGrid(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Illegal fix field/grid command");

  int ncols = 0;

  if (strcmp(arg[2],"NULL") == 0) axstr = NULL;
  else {
    int n = strlen(arg[2]) + 1;
    axstr = new char[n];
    strcpy(axstr,arg[2]);
    ncols++;
  }
  if (strcmp(arg[3],"NULL") == 0) aystr = NULL;
  else {
    int n = strlen(arg[3]) + 1;
    aystr = new char[n];
    strcpy(aystr,arg[3]);
    ncols++;
  }
  if (strcmp(arg[4],"NULL") == 0) azstr = NULL;
  else {
    int n = strlen(arg[4]) + 1;
    azstr = new char[n];
    strcpy(azstr,arg[4]);
    ncols++;
  }

  // fix settings

  per_grid_flag = 1;
  size_per_grid_cols = ncols;
  per_grid_freq = 1;
  per_grid_field = 1;

  field_active[0] = field_active[1] = field_active[2] = 0;
  if (axstr) field_active[0] = 1;
  if (aystr) field_active[1] = 1;
  if (azstr) field_active[2] = 1;

  // per-grid memory initialization

  maxgrid = 0;
  array_grid = NULL;
}

/* ---------------------------------------------------------------------- */

FixFieldGrid::~FixFieldGrid()
{
  delete [] axstr;
  delete [] aystr;
  delete [] azstr;

  memory->destroy(array_grid);
}

/* ---------------------------------------------------------------------- */

int FixFieldGrid::setmask()
{
  int mask = 0;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixFieldGrid::init()
{
  // check if all variables exist and are grid-style vars

  if (axstr) {
    axvar = input->variable->find(axstr);
    if (axvar < 0)
      error->all(FLERR,"Variable name for fix field/grid does not exist");
    if (!input->variable->grid_style(axvar))
      error->all(FLERR,"Variable for fix field/grid is invalid style");
  }
  if (aystr) {
    ayvar = input->variable->find(aystr);
    if (ayvar < 0)
      error->all(FLERR,"Variable name for fix field/grid does not exist");
    if (!input->variable->grid_style(ayvar))
      error->all(FLERR,"Variable for fix field/grid is invalid style");
  }
  if (azstr) {
    azvar = input->variable->find(azstr);
    if (azvar < 0)
      error->all(FLERR,"Variable name for fix field/grid does not exist");
    if (!input->variable->grid_style(azvar))
      error->all(FLERR,"Variable for fix field/grid is invalid style");
  }

  // set initial grid values to zero in case dump is performed at step 0

  if (grid->nlocal > maxgrid) {
    maxgrid = grid->maxlocal;
    memory->destroy(array_grid);
    memory->create(array_grid,maxgrid,size_per_grid_cols,"array_grid");
  }

  bigint nbytes = (bigint) grid->nlocal * size_per_grid_cols;
  memset(&array_grid[0][0],0,nbytes*sizeof(double));
}

/* ---------------------------------------------------------------------- */

void FixFieldGrid::compute_field()
{
  // reallocate array_grid if necessary

  if (grid->nlocal > maxgrid) {
    maxgrid = grid->maxlocal;
    memory->destroy(array_grid);
    memory->create(array_grid,maxgrid,size_per_grid_cols,"array_grid");
  }

  // evaluate each grid-style variable
  // results are put into strided array_grid

  int stride = size_per_grid_cols;
  int icol = 0;

  if (axstr) {
    input->variable->compute_grid(axvar,&array_grid[0][icol],stride,0);
    icol++;
  }

  if (aystr) {
    input->variable->compute_grid(ayvar,&array_grid[0][icol],stride,0);
    icol++;
  }

  if (azstr) {
    input->variable->compute_grid(azvar,&array_grid[0][icol],stride,0);
    icol++;
  }
}

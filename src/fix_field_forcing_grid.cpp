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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_field_forcing_grid.h"
#include "update.h"
#include "grid.h"
#include "input.h"
#include "variable.h"
#include "random_knuth.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

FixFieldForcingGrid::FixFieldForcingGrid(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix field/forcing/grid command");

  epsilon = atof(arg[2]);
  nkvec = atoi(arg[3]);

  if (nkvec == 0) error->all(FLERR,"Fix field/forcing/grid requires Kspace vecs");
  if (narg-4 != 3*nkvec) 
    error->all(FLERR,"Illegal fix field/forcing/grid command");

  memory->create(kvecs,nkvec,3,"");

  int iarg = 4;
  for (int i = 0; i < nkvec; i++) {
    kvecs[i][0] = atof(arg[iarg]);
    kvecs[i][1] = atof(arg[iarg+1]);
    kvecs[i][2] = atof(arg[iarg+2]);
    iarg += 3;
  }

  // fix settings

  per_grid_flag = 1;
  size_per_grid_cols = 3;
  per_grid_freq = 1;
  per_grid_field = 1;
  field_active[0] = field_active[1] = field_active[2] = 1;

  // random number generator, same on every proc

  random = new RanKnuth(update->ranmaster->uniform());

  // per-grid memory initialization

  maxgrid = 0;
  array_grid = NULL;
}

/* ---------------------------------------------------------------------- */

FixFieldForcingGrid::~FixFieldForcingGrid()
{  
  memory->destroy(kvecs);
  memory->destroy(array_grid);
  delete random;
}

/* ---------------------------------------------------------------------- */

int FixFieldForcingGrid::setmask()
{
  int mask = 0;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixFieldForcingGrid::init()
{
  // set initial grid values to zero in case dump is performed at step 0

  if (grid->nlocal > maxgrid) {
    maxgrid = grid->maxlocal;
    memory->destroy(array_grid);
    memory->create(array_grid,maxgrid,3,"field/forcing:array_grid");
  }

  bigint nbytes = (bigint) grid->nlocal * 3;
  memset(&array_grid[0][0],0,nbytes*sizeof(double));
}

/* ---------------------------------------------------------------------- */

void FixFieldForcingGrid::compute_field()
{
  // reallocate array_grid if necessary

  if (grid->nlocal > maxgrid) {
    maxgrid = grid->maxlocal;
    memory->destroy(array_grid);
    memory->create(array_grid,maxgrid,3,"field/forcing:array_grid");
  }

  // RYAN: this is code you need to write
  // set current 48 global prefactors for spatial/time varying field
  // can use epsilon and any other arguments to this fix
  // can use random->uniform() to get a uniform RN between 0 and 1

  double rn = random->uniform();
  double prefacN = rn;

  // RYAN: this is more code you need to write
  // array_grid = current field at each grid cell center point

  Grid::ChildCell *cells = grid->cells;
  int nglocal = grid->nlocal;

  double xc,yc,zc;

  for (int icell = 0; icell < nglocal; icell++) {

    // center point of grid cell

    xc = 0.5 * (cells[icell].lo[0] + cells[icell].hi[0]);
    yc = 0.5 * (cells[icell].lo[1] + cells[icell].hi[1]);
    zc = 0.5 * (cells[icell].lo[2] + cells[icell].hi[2]);

    // loop over Kspace vectors
    // RYAN: domain->boxlo and boxhi are global simulation box bounds
    //       domain->prd = simulation box size
    //       something like (xc-boxlo[0])/prd[0] 
    //         may be what you want for "x" in dot product ?

    array_grid[icell][0] = 0.0;
    array_grid[icell][1] = 0.0;
    array_grid[icell][2] = 0.0;

    for (int ik = 0; ik < nkvec; ik++) {
      double kdotx = kvecs[ik][0]*xc + kvecs[ik][1]*yc + kvecs[ik][2]*zc;
      double c = cos(kdotx);
      double s = sin(kdotx);

      array_grid[icell][0] += prefacN * c + prefacN * s;
      array_grid[icell][1] += prefacN * c + prefacN * s;
      array_grid[icell][2] += prefacN * c + prefacN * s;
    }
  }
}

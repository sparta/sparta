/* ----------------------------------------------------------------------
   DSMC - Sandia parallel DSMC code
   www.sandia.gov/~sjplimp/dsmc.html
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2011) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level DSMC directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"
#include "grid.h"
#include "domain.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace DSMC_NS;

enum{BLOCK,RANDOM};

/* ---------------------------------------------------------------------- */

Grid::Grid(DSMC *dsmc) : Pointers(dsmc)
{
  grid_exist = 0;

  ncell = 0;
  cells = NULL;
}

/* ---------------------------------------------------------------------- */

Grid::~Grid()
{
  memory->sfree(cells);
}

/* ---------------------------------------------------------------------- */

void Grid::create(int narg, char **arg)
{
  if (!domain->box_exist)
    error->all(FLERR,"Cannot create grid before simulation box is defined");
  if (grid_exist)
    error->all(FLERR,"Cannot create grid when grid is already defined");

  grid_exist = 1;

  if (narg != 4) error->all(FLERR,"Illegal create_grid command");

  nx = atoi(arg[0]);
  ny = atoi(arg[1]);
  nz = atoi(arg[2]);

  if (nx < 1 || ny < 1 || nz < 1)
    error->all(FLERR,"Illegal create_grid command");

  if (strcmp(arg[3],"block") == 0) bstyle = BLOCK;
  else if (strcmp(arg[3],"random") == 0) bstyle = RANDOM;
  else error->all(FLERR,"Illegal create_grid command");

  // box and grid cell geometry

  int dimension = domain->dimension;
  double xlo = domain->boxlo[0];
  double ylo = domain->boxlo[1];
  double zlo = domain->boxlo[2];
  double xhi = domain->boxhi[0];
  double yhi = domain->boxhi[1];
  double zhi = domain->boxhi[2];
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  if (dimension == 2 && nz != 1)
    error->all(FLERR,"Create_grid defines multiple z cells for 2d domain");

  xdelta = xprd / nx;
  ydelta = yprd / ny;
  zdelta = zprd / nz;

  // build a regular Nx x Ny x Nz global grid
  
  ncell = nx*ny*nz;
  // check on exceed smallint
  cells = (OneCell *) memory->smalloc(ncell*sizeof(OneCell),"grid:cells");

  int i,j,k,m;

  m = 0;
  for (k = 0; k < nz; k++) {
    for (j = 0; j < ny; j++) {
      for (i = 0; i < nx; i++) {
	cells[m].lo[0] = xlo + i*xdelta;
	cells[m].lo[1] = ylo + j*ydelta;
	cells[m].lo[2] = zlo + k*zdelta;
	cells[m].hi[0] = xlo + (i+1)*xdelta;
	cells[m].hi[1] = ylo + (j+1)*ydelta;
	cells[m].hi[2] = zlo + (k+1)*zdelta;
	if (i == nx-1) cells[m].hi[0] = xhi;
	if (j == ny-1) cells[m].hi[1] = yhi;
	if (k == nz-1) cells[m].hi[2] = zhi;

	cells[m].neigh[0] = m - 1;
	cells[m].neigh[1] = m + 1;
	cells[m].neigh[2] = m - nx;
	cells[m].neigh[3] = m + nx;
	cells[m].neigh[4] = m - nx*ny;
	cells[m].neigh[5] = m + nx*ny;

	if (i == 0) cells[m].neigh[0] = -1;
	if (i == nx-1) cells[m].neigh[1] = -1;
	if (j == 0) cells[m].neigh[2] = -1;
	if (j == ny-1) cells[m].neigh[3] = -1;
	if (k == 0) cells[m].neigh[4] = -1;
	if (k == nz-1) cells[m].neigh[5] = -1;

	cells[m].id = m+1;
	m++;
      }
    }
  }

  // assign owner to each grid cell, based on specified bstyle

  for (m = 0; m < ncell; m++) cells[m].proc = 0;

  // stats

  if (comm->me == 0) {
    if (screen) fprintf(screen,"Created %d grid cells\n",ncell);
    if (logfile) fprintf(logfile,"Created %d grid cells\n",ncell);
  }
}

/* ---------------------------------------------------------------------- */

int Grid::which_cell(double x, double y, double z)
{
  int ix = x / xdelta;
  int iy = y / ydelta;
  int iz = z / zdelta;
  int icell = iz*nx*ny + iy*nx + ix;
  return icell;
}

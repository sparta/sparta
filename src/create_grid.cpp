/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2012) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"
#include "create_grid.h"
#include "grid.h"
#include "domain.h"
#include "comm.h"
#include "error.h"

using namespace SPARTA_NS;

enum{STRIDE,BLOCK,RANDOM};
enum{XYZ,XZY,YXZ,YZX,ZXY,ZYX};

/* ---------------------------------------------------------------------- */

CreateGrid::CreateGrid(SPARTA *sparta) : Pointers(sparta) {}

/* ---------------------------------------------------------------------- */

void CreateGrid::command(int narg, char **arg)
{
  if (!domain->box_exist) 
    error->all(FLERR,"Cannot create_grid before simulation box is defined");
  if (grid->grid_exist)
    error->all(FLERR,"Cannot create grid when grid is already defined");

  grid->grid_exist = 1;

  if (narg < 3) error->all(FLERR,"Illegal create_grid command");

  int nx = atoi(arg[0]);
  int ny = atoi(arg[1]);
  int nz = atoi(arg[2]);

  if (nx < 1 || ny < 1 || nz < 1)
    error->all(FLERR,"Illegal create_grid command");
  if (domain->dimension == 2 && nz != 1)
    error->all(FLERR,"Create_grid nz value must be 1 for a 2d simulation");

  // optional args

  int bstyle = BLOCK;
  int px = 0;
  int py = 0;
  int pz = 0;
  int order;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"stride") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal create_grid command");
      bstyle = STRIDE;
      if (strlen(arg[iarg+1]) + strlen(arg[iarg+2]) + strlen(arg[iarg+3]) != 3)
	error->all(FLERR,"Illegal create_grid command");
      char str[4];
      str[0] = arg[iarg+1][0];
      str[1] = arg[iarg+2][0];
      str[2] = arg[iarg+3][0];
      str[3] = 0;
      if (strcmp(str,"xyz") == 0) order = XYZ;
      else if (strcmp(str,"xzy") == 0) order = XZY;
      else if (strcmp(str,"yxz") == 0) order = YXZ;
      else if (strcmp(str,"yzx") == 0) order = YZX;
      else if (strcmp(str,"zxy") == 0) order = ZXY;
      else if (strcmp(str,"zyx") == 0) order = ZYX;
      else error->all(FLERR,"Illegal create_grid command");
      iarg += 4;

    } else if (strcmp(arg[3],"block") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal create_grid command");
      bstyle = BLOCK;
      if (strcmp(arg[iarg+1],"*") == 0) px = 0;
      else px = atoi(arg[iarg+1]);
      if (strcmp(arg[iarg+1],"*") == 0) py = 0;
      else py = atoi(arg[iarg+2]);
      if (strcmp(arg[iarg+1],"*") == 0) pz = 0;
      else pz = atoi(arg[iarg+3]);
      iarg += 4;

    } else if (strcmp(arg[3],"random") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal create_grid command");
      bstyle = RANDOM;
      iarg += 1;

    } else error->all(FLERR,"Illegal create_grid command");
  }

  // box and grid cell geometry

  double xlo = domain->boxlo[0];
  double ylo = domain->boxlo[1];
  double zlo = domain->boxlo[2];
  double xhi = domain->boxhi[0];
  double yhi = domain->boxhi[1];
  double zhi = domain->boxhi[2];
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  double xdelta = xprd / nx;
  double ydelta = yprd / ny;
  double zdelta = zprd / nz;
  double xdeltainv = nx / xprd;
  double ydeltainv = ny / yprd;
  double zdeltainv = nz / zprd;

  // build a regular Nx x Ny x Nz global grid
  // neigh[face] = -1 if cell is adjacent to global boundary
  
  bigint ntotal = (bigint) nx * ny * nz;
  if (ntotal > MAXSMALLINT) 
    error->one(FLERR,"Per-processor grid count is too big");

  int i,j,k,m;
  int neigh[6];
  double lo[3],hi[3];

  m = 0;
  for (k = 0; k < nz; k++) {
    for (j = 0; j < ny; j++) {
      for (i = 0; i < nx; i++) {
	lo[0] = xlo + i*xdelta;
	lo[1] = ylo + j*ydelta;
	lo[2] = zlo + k*zdelta;
	hi[0] = xlo + (i+1)*xdelta;
	hi[1] = ylo + (j+1)*ydelta;
	hi[2] = zlo + (k+1)*zdelta;
	if (i == nx-1) hi[0] = xhi;
	if (j == ny-1) hi[1] = yhi;
	if (k == nz-1) hi[2] = zhi;

	neigh[0] = m - 1;
	neigh[1] = m + 1;
	neigh[2] = m - nx;
	neigh[3] = m + nx;
	neigh[4] = m - nx*ny;
	neigh[5] = m + nx*ny;

	if (i == 0) neigh[0] = -1;
	if (i == nx-1) neigh[1] = -1;
	if (j == 0) neigh[2] = -1;
	if (j == ny-1) neigh[3] = -1;
	if (k == 0) neigh[4] = -1;
	if (k == nz-1) neigh[5] = -1;

	grid->add_cell(m+1,lo,hi,neigh);
	m++;
      }
    }
  }

  // set grid geometry values within Grid class

  grid->nx = nx;
  grid->ny = ny;
  grid->nz = nz;
  grid->xdelta = xdelta;
  grid->ydelta = ydelta;
  grid->zdelta = zdelta;
  grid->xdeltainv = xdeltainv;
  grid->ydeltainv = ydeltainv;
  grid->zdeltainv = zdeltainv;

  // assign cells to processors based on bstyle

  if (bstyle == STRIDE) grid->assign_stride(order);
  else if (bstyle == BLOCK) grid->assign_block(px,py,pz);
  else if (bstyle == RANDOM) grid->assign_random();
  
  // make list of cells I own

  grid->setup_grid();

  // stats

  if (comm->me == 0) {
    if (screen) fprintf(screen,"Created %d grid cells\n",grid->ncell);
    if (logfile) fprintf(logfile,"Created %d grid cells\n",grid->ncell);
  }
}

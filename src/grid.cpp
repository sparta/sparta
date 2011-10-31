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
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace DSMC_NS;

enum{STRIDE,BLOCK,RANDOM};
enum{XYZ,XZY,YXZ,YZX,ZXY,ZYX};

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

  if (narg < 3) error->all(FLERR,"Illegal create_grid command");

  nx = atoi(arg[0]);
  ny = atoi(arg[1]);
  nz = atoi(arg[2]);

  if (nx < 1 || ny < 1 || nz < 1)
    error->all(FLERR,"Illegal create_grid command");
  if (domain->dimension == 2 && nz != 1)
    error->all(FLERR,"Create_grid nz value must be 1 for a 2d simulation");

  // optional args

  bstyle = BLOCK;
  user_procgrid[0] = user_procgrid[1] = user_procgrid[2] = 0;

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
      error->all(FLERR,"Illegal create_grid command");
      iarg += 4;

    } else if (strcmp(arg[3],"block") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal create_grid command");
      bstyle = BLOCK;
      if (strcmp(arg[iarg+1],"*") == 0) user_procgrid[0] = 0;
      else user_procgrid[0] = atoi(arg[iarg+1]);
      if (strcmp(arg[iarg+1],"*") == 0) user_procgrid[1] = 0;
      else user_procgrid[1] = atoi(arg[iarg+2]);
      if (strcmp(arg[iarg+1],"*") == 0) user_procgrid[2] = 0;
      else user_procgrid[2] = atoi(arg[iarg+3]);
      iarg += 4;

    } else if (strcmp(arg[3],"random") == 0) {
      if (narg != 5) error->all(FLERR,"Illegal create_grid command");
      bstyle = RANDOM;
      seed = atoi(arg[iarg+1]);
      if (seed <= 0) error->all(FLERR,"Illegal create_grid command");
      iarg += 2;

    } else error->all(FLERR,"Illegal create_grid command");
  }

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

  xdelta = xprd / nx;
  ydelta = yprd / ny;
  zdelta = zprd / nz;
  xdeltainv = nx / xprd;
  ydeltainv = ny / yprd;
  zdeltainv = nz / zprd;

  // build a regular Nx x Ny x Nz global grid
  
  bigint ntotal = (bigint) nx * ny * nz;
  if (ntotal > MAXSMALLINT) 
    error->one(FLERR,"Per-processor grid count is too big");
  ncell = ntotal;
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

  // assign cells to processors based on bstyle
  // calculates nlocal = # of cells I own

  if (bstyle == STRIDE) assign_stride();
  else if (bstyle == BLOCK) assign_block();
  else if (bstyle == RANDOM) assign_random();
  
  // stats

  if (comm->me == 0) {
    if (screen) fprintf(screen,"Created %d grid cells\n",ncell);
    if (logfile) fprintf(logfile,"Created %d grid cells\n",ncell);
  }
}

/* ---------------------------------------------------------------------- */

void Grid::assign_stride()
{
  int ix,iy,iz,nth;

  int me = comm->me;
  int nprocs = comm->nprocs;
  
  nlocal = 0;
  for (int m = 0; m < ncell; m++) {
    ix = m % nx;
    iy = (m / nx) % ny;
    iz = m / (nx*ny);
    
    if (order == XYZ) nth = iz*nx*ny + iy*nx + ix;
    else if (order == XZY) nth = iy*nx*nz + iz*nx + ix;
    else if (order == YXZ) nth = iz*ny*nx + ix*ny + iy;
    else if (order == YZX) nth = ix*ny*nz + iz*ny + iy;
    else if (order == ZXY) nth = iy*nz*nx + ix*nz + iz;
    else if (order == ZYX) nth = ix*nz*ny + iy*nz + iz;

    cells[m].proc = nth % nprocs;
    if (cells[m].proc == me) nlocal++;
  }
}

/* ---------------------------------------------------------------------- */

void Grid::assign_block()
{
  procs2grid();
  if (procgrid[0]*procgrid[1]*procgrid[2] != comm->nprocs)
    error->all(FLERR,"Bad grid of processors for create_grid");

  int ix,iy,iz,ipx,ipy,ipz,iproc;
  int me = comm->me;

  nlocal = 0;
  for (int m = 0; m < ncell; m++) {
    ix = m % nx;
    iy = (m / nx) % ny;
    iz = m / (nx*ny);
    ipx = ix*procgrid[0] / nx;
    ipy = iy*procgrid[1] / ny;
    ipz = iz*procgrid[2] / nz;
    iproc = ipz*procgrid[0]*procgrid[1] + ipy*procgrid[0] + ipx;
    cells[m].proc = iproc;
    if (cells[m].proc == me) nlocal++;
  }
}

/* ---------------------------------------------------------------------- */

void Grid::assign_random()
{
  int me = comm->me;
  int nprocs = comm->nprocs;
  RanPark *random = new RanPark(dsmc,seed);

  nlocal = 0;
  for (int m = 0; m < ncell; m++) {
    cells[m].proc = nprocs * random->uniform();
    if (cells[m].proc == me) nlocal++;
  }

  delete random;
}

/* ---------------------------------------------------------------------- */

int Grid::which_cell(double x, double y, double z)
{
  int ix = x * xdeltainv;
  int iy = y * ydeltainv;
  int iz = z * zdeltainv;
  int icell = iz*nx*ny + iy*nx + ix;
  return icell;
}

/* ----------------------------------------------------------------------
   assign nprocs to 3d grid so as to minimize surface area 
   area = surface area of each of 3 faces of simulation box
------------------------------------------------------------------------- */

void Grid::procs2grid()
{
  int nprocs = comm->nprocs;
  procgrid[0] = user_procgrid[0];
  procgrid[1] = user_procgrid[1];
  procgrid[2] = user_procgrid[2];

  // all 3 proc counts are specified

  if (procgrid[0] && procgrid[1] && procgrid[2]) return;

  // 2 out of 3 proc counts are specified

  if (procgrid[0] > 0 && procgrid[1] > 0) {
    procgrid[2] = nprocs/(procgrid[0]*procgrid[1]);
    return;
  } else if (procgrid[0] > 0 && procgrid[2] > 0) {
    procgrid[1] = nprocs/(procgrid[0]*procgrid[2]);
    return;
  } else if (procgrid[1] > 0 && procgrid[2] > 0) {
    procgrid[0] = nprocs/(procgrid[1]*procgrid[2]);
    return;
  } 

  // determine cross-sectional areas
  // area[0] = xy, area[1] = xz, area[2] = yz

  double area[3];
  area[0] = nx*ny;
  area[1] = nx*nz;
  area[2] = ny*nz;

  double bestsurf = 2.0 * (area[0]+area[1]+area[2]);

  // loop thru all possible factorizations of nprocs
  // only consider valid cases that match procgrid settings
  // surf = surface area of a proc sub-domain

  int ipx,ipy,ipz,valid;
  double surf;

  ipx = 1;
  while (ipx <= nprocs) {
    valid = 1;
    if (user_procgrid[0] && ipx != user_procgrid[0]) valid = 0;
    if (nprocs % ipx) valid = 0;
    if (!valid) {
      ipx++;
      continue;
    }

    ipy = 1;
    while (ipy <= nprocs/ipx) {
      valid = 1;
      if (user_procgrid[1] && ipy != user_procgrid[1]) valid = 0;
      if ((nprocs/ipx) % ipy) valid = 0;
      if (!valid) {
	ipy++;
	continue;
      }
      
      ipz = nprocs/ipx/ipy;
      valid = 1;
      if (user_procgrid[2] && ipz != user_procgrid[2]) valid = 0;
      if (domain->dimension == 2 && ipz != 1) valid = 0;
      if (!valid) {
	ipy++;
	continue;
      }
      
      surf = area[0]/ipx/ipy + area[1]/ipx/ipz + area[2]/ipy/ipz;
      if (surf < bestsurf) {
	bestsurf = surf;
	procgrid[0] = ipx;
	procgrid[1] = ipy;
	procgrid[2] = ipz;
      }
      ipy++;
    }

    ipx++;
  }
}

/* ---------------------------------------------------------------------- */

bigint Grid::memory_usage()
{
  bigint bytes = (bigint) ncell * sizeof(OneCell);
  return bytes;
}



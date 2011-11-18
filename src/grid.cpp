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

#include "grid.h"
#include "domain.h"
#include "comm.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace DSMC_NS;

#define DELTA 10000

enum{XYZ,XZY,YXZ,YZX,ZXY,ZYX};

/* ---------------------------------------------------------------------- */

Grid::Grid(DSMC *dsmc) : Pointers(dsmc)
{
  grid_exist = 0;

  ncell = maxcell = 0;
  cells = NULL;

  nlocal = 0;
  mycells = NULL;
}

/* ---------------------------------------------------------------------- */

Grid::~Grid()
{
  memory->sfree(cells);
  memory->sfree(mycells);
}

/* ---------------------------------------------------------------------- */

void Grid::add_cell(int id, double *lo, double *hi, int *neigh)
{
  if (ncell == maxcell) grow(1);

  OneCell *c = &cells[ncell];

  c->id = id;
  c->lo[0] = lo[0];
  c->lo[1] = lo[1];
  c->lo[2] = lo[2];
  c->hi[0] = hi[0];
  c->hi[1] = hi[1];
  c->hi[2] = hi[2];
  c->neigh[0] = neigh[0];
  c->neigh[1] = neigh[1];
  c->neigh[2] = neigh[2];
  c->neigh[3] = neigh[3];
  c->neigh[4] = neigh[4];
  c->neigh[5] = neigh[5];

  ncell++;
}

/* ----------------------------------------------------------------------
   setup owned grid cells
   create mycells list of owned cells
   compute volume of owned cells
------------------------------------------------------------------------- */

void Grid::setup_grid()
{
  // nlocal = # of cells I own
  // mycells = indices of cells I own

  int me = comm->me;

  nlocal = 0;
  for (int m = 0; m < ncell; m++)
    if (cells[m].proc == me) nlocal++;

  memory->destroy(mycells);
  memory->create(mycells,nlocal,"grid:mycells");

  nlocal = 0;
  for (int m = 0; m < ncell; m++)
    if (cells[m].proc == me) mycells[nlocal++] = m;

  // calculate volume of cells I own

  int icell;
  double dx,dy,dz;

  for (int m = 0; m < nlocal; m++) {
    icell = mycells[m];
    dx = cells[icell].hi[0] - cells[icell].lo[0];
    dy = cells[icell].hi[1] - cells[icell].lo[1];
    dz = cells[icell].hi[2] - cells[icell].lo[2];
    cells[icell].volume = dx*dy*dz;
  }
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

/* ---------------------------------------------------------------------- */

void Grid::assign_stride(int order)
{
  int ix,iy,iz,nth;

  int me = comm->me;
  int nprocs = comm->nprocs;
  
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
  }
}

/* ---------------------------------------------------------------------- */

void Grid::assign_block(int px, int py, int pz)
{
  procs2grid(px,py,pz);
  if (px*py*pz != comm->nprocs)
    error->all(FLERR,"Bad grid of processors for create_grid");

  int ix,iy,iz,ipx,ipy,ipz,iproc;
  int me = comm->me;

  for (int m = 0; m < ncell; m++) {
    ix = m % nx;
    iy = (m / nx) % ny;
    iz = m / (nx*ny);
    ipx = ix*px / nx;
    ipy = iy*py / ny;
    ipz = iz*pz / nz;
    iproc = ipz*px*py + ipy*px + ipx;
    cells[m].proc = iproc;
  }
}

/* ---------------------------------------------------------------------- */

void Grid::assign_random(int seed)
{
  int me = comm->me;
  int nprocs = comm->nprocs;
  RanPark *random = new RanPark(dsmc,seed);

  for (int m = 0; m < ncell; m++)
    cells[m].proc = nprocs * random->uniform();

  delete random;
}

/* ----------------------------------------------------------------------
   assign nprocs to 3d grid so as to minimize surface area 
   area = surface area of each of 3 faces of simulation box
------------------------------------------------------------------------- */

void Grid::procs2grid(int &px, int &py, int &pz)
{
  int upx = px;
  int upy = py;
  int upz = pz;

  int nprocs = comm->nprocs;

  // all 3 proc counts are specified

  if (px && py && pz) return;

  // 2 out of 3 proc counts are specified

  if (py > 0 && pz > 0) {
    px = nprocs/(py*pz);
    return;
  } else if (px > 0 && pz > 0) {
    py = nprocs/(px*pz);
    return;
  } else if (px > 0 && py > 0) {
    pz = nprocs/(px*py);
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
    if (upx && ipx != upx) valid = 0;
    if (nprocs % ipx) valid = 0;
    if (!valid) {
      ipx++;
      continue;
    }

    ipy = 1;
    while (ipy <= nprocs/ipx) {
      valid = 1;
      if (upy && ipy != upy) valid = 0;
      if ((nprocs/ipx) % ipy) valid = 0;
      if (!valid) {
	ipy++;
	continue;
      }
      
      ipz = nprocs/ipx/ipy;
      valid = 1;
      if (upz && ipz != upz) valid = 0;
      if (domain->dimension == 2 && ipz != 1) valid = 0;
      if (!valid) {
	ipy++;
	continue;
      }
      
      surf = area[0]/ipx/ipy + area[1]/ipx/ipz + area[2]/ipy/ipz;
      if (surf < bestsurf) {
	bestsurf = surf;
	px = ipx;
	py = ipy;
	pz = ipz;
      }
      ipy++;
    }

    ipx++;
  }
}

/* ----------------------------------------------------------------------
   insure cell list can hold nextra new grid cells
------------------------------------------------------------------------- */

void Grid::grow(int nextra)
{
  bigint target = (bigint) ncell + nextra;
  if (target <= maxcell) return;
  
  bigint newmax = maxcell;
  while (newmax < target) newmax += DELTA;
  
  if (newmax > MAXSMALLINT) 
    error->one(FLERR,"Per-processor grid count is too big");

  maxcell = newmax;
  cells = (OneCell *)
    memory->srealloc(cells,maxcell*sizeof(OneCell),"grid:cells");
}

/* ---------------------------------------------------------------------- */

bigint Grid::memory_usage()
{
  bigint bytes = (bigint) ncell * sizeof(OneCell);
  return bytes;
}

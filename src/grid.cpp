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

#include "math.h"
#include "grid.h"
#include "geometry.h"
#include "update.h"
#include "domain.h"
#include "comm.h"
#include "surf.h"
#include "random_mars.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace DSMC_NS;

#define DELTA 10000
#define EPSILON 1.0e-6
#define BIG 1.0e20

enum{XYZ,XZY,YXZ,YZX,ZXY,ZYX};
enum{SURFEXTERIOR,SURFINTERIOR,SURFCONTAIN};    // same as CreateMolecules
                                                // same as FixInflow

/* ---------------------------------------------------------------------- */

Grid::Grid(DSMC *dsmc) : Pointers(dsmc)
{
  grid_exist = 0;

  ncell = maxcell = 0;
  cells = NULL;

  nlocal = 0;
  mycells = NULL;
  csurfs = NULL;
}

/* ---------------------------------------------------------------------- */

Grid::~Grid()
{
  memory->sfree(cells);
  memory->destroy(mycells);
  memory->destroy(csurfs);
}

/* ---------------------------------------------------------------------- */

void Grid::init()
{
  // assign surf list to each global cell

  if (surf->surf_exist) {
    int i,m;
    double cmax,len,area;
    int dimension = domain->dimension;

    surf2grid();

    int icell;
    int stotal = 0;
    int smax = 0;
    double sratio = BIG;
    for (int m = 0; m < nlocal; m++) {
      icell = mycells[m];
      stotal += cells[icell].nsurf;
      smax = MAX(smax,cells[icell].nsurf);

      cmax = MAX(cells[icell].hi[0] - cells[icell].lo[0],
		 cells[icell].hi[1] - cells[icell].lo[1]);
      if (dimension == 3) 
	cmax = MAX(cmax,cells[icell].hi[2] - cells[icell].lo[2]);

      for (int i = 0; i < cells[icell].nsurf; i++) {
	surf->tri_size(csurfs[icell][i],len,area);
	sratio = MIN(sratio,len/cmax);
      }
    }

    int stotalall,smaxall;
    double sratioall;
    MPI_Allreduce(&stotal,&stotalall,1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&smax,&smaxall,1,MPI_INT,MPI_MAX,world);
    MPI_Allreduce(&sratio,&sratioall,1,MPI_DOUBLE,MPI_MIN,world);

    if (comm->me == 0) {
      if (screen) {
	fprintf(screen,"Grid/surf-element stats:\n");
	fprintf(screen,"  total surfs in all grid cells = %d\n",stotalall);
	fprintf(screen,"  max surfs in one grid cell = %d\n",smaxall);
	fprintf(screen,"  min surf-size/cell-size ratio = %g\n",sratioall);
      }
      if (logfile) {
	fprintf(logfile,"Grid/surf-element stats:\n");
	fprintf(logfile,"  total surfs in all grid cells = %d\n",stotalall);
	fprintf(logfile,"  max surfs in one grid cell = %d\n",smaxall);
	fprintf(logfile,"  min surf-size/cell-size ratio = %g\n",sratioall);
      }
    }
  }

  // set inflag for each owned cell

  if (surf->surf_exist) grid_inout();
  else {
    int icell;
    for (int m = 0; m < nlocal; m++) {
      icell = mycells[m];
      cells[icell].inflag = SURFEXTERIOR; 
    }
  }
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
  c->nsurf = 0;

  ncell++;
}

/* ----------------------------------------------------------------------
   setup owned grid cells
   create mycells list of owned cells
   compute local index for owned cells
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
    if (cells[m].proc == me) {
      cells[m].local = nlocal;
      mycells[nlocal++] = m;
    }

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

/* ----------------------------------------------------------------------
   map a particle coordinate into a grid cell
   NOTE: not currently used, see loop option in CreateMolecules
   NOTE: what if particle is at upper boundary of domain
   NOTE: assumes Nx by Ny by Nz grid
------------------------------------------------------------------------- */

int Grid::which_cell(double x, double y, double z)
{
  int ix = static_cast<int> (x * xdeltainv);
  int iy = static_cast<int> (y * ydeltainv);
  int iz = static_cast<int> (z * zdeltainv);
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

void Grid::assign_random()
{
  int me = comm->me;
  int nprocs = comm->nprocs;
  RanPark *random = new RanPark(update->ranmaster->uniform());

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
   assign surface elements (lines or triangles) to grid cells
   NOTE: no parallelism yet, since each proc owns entire grid & all surfs
------------------------------------------------------------------------- */

void Grid::surf2grid()
{
  int i,j,k,m,icell;
  int ilo,ihi,jlo,jhi,klo,khi;
  double lo[3],hi[3];
  double *x1,*x2,*x3;

  // epsilon = EPSILON fraction of largest box length

  int dimension = domain->dimension;
  double bmax = MAX(domain->xprd,domain->yprd);
  if (dimension == 3) bmax = MAX(bmax,domain->zprd);
  double epsilon = EPSILON * bmax;

  // count[M] = # of triangles overlapping grid cell M

  int *count;
  memory->create(count,ncell,"grid:count");
  for (m = 0; m < ncell; m++) count[m] = 0;

  // NOTE: this logic is specific to regular Nx by Ny by Nz grid
  // tally count by double loop over surfs and grid cells within surf bbox
  // lo/hi = bounding box around surf
  // ijk lo/hi = grid index bounding box around surf
  // add epsilon to insure surf is counted in any cell it touches
  // icell = index of a grid cell within bounding box

  Surf::Point *pts = surf->pts;
  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;

  int nsurf;
  if (dimension == 2) nsurf = surf->nline;
  else nsurf = surf->ntri;

  for (m = 0; m < nsurf; m++) {
    if (dimension == 2) {
      x1 = pts[lines[m].p1].x;
      x2 = pts[lines[m].p2].x;
      x3 = x2;
    } else {
      x1 = pts[tris[m].p1].x;
      x2 = pts[tris[m].p2].x;
      x3 = pts[tris[m].p3].x;
    }

    lo[0] = MIN(x1[0],x2[0]);
    lo[0] = MIN(lo[0],x3[0]);
    hi[0] = MAX(x1[0],x2[0]);
    hi[0] = MAX(hi[0],x3[0]);

    lo[1] = MIN(x1[1],x2[1]);
    lo[1] = MIN(lo[1],x3[1]);
    hi[1] = MAX(x1[1],x2[1]);
    hi[1] = MAX(hi[1],x3[1]);

    lo[2] = MIN(x1[2],x2[2]);
    lo[2] = MIN(lo[2],x3[2]);
    hi[2] = MAX(x1[2],x2[2]);
    hi[2] = MAX(hi[2],x3[2]);

    ilo = MAX(0,static_cast<int> ((lo[0]-epsilon)*xdeltainv));
    ihi = MIN(nx-1,static_cast<int> ((hi[0]+epsilon)*xdeltainv));
    jlo = MAX(0,static_cast<int> ((lo[1]-epsilon)*ydeltainv));
    jhi = MIN(ny-1,static_cast<int> ((hi[1]+epsilon)*ydeltainv));
    klo = MAX(0,static_cast<int> ((lo[2]-epsilon)*zdeltainv));
    khi = MIN(nz-1,static_cast<int> ((hi[2]+epsilon)*zdeltainv));

    if (dimension == 2) {
      for (i = ilo; i <= ihi; i++)
	for (j = jlo; j <= jhi; j++) {
	  icell = j*nx + i;
	  if (Geometry::line_quad_intersect(x1,x2,lines[m].norm,
					    cells[icell].lo,cells[icell].hi))
	    count[icell]++;
	}
    } else {
      for (i = ilo; i <= ihi; i++)
	for (j = jlo; j <= jhi; j++)
	  for (k = klo; k <= khi; k++) {
	    icell = k*nx*ny + j*nx + i;
	    if (Geometry::tri_hex_intersect(x1,x2,x3,tris[m].norm,
					    cells[icell].lo,cells[icell].hi))
	      count[icell]++;
	  }
    }
  }

  // (re)allocate ragged csurfs array
  // csurfs[I][J] = index of Jth tri in global cell I

  memory->destroy(csurfs);
  memory->create_ragged(csurfs,ncell,count,"grid:csurfs");

  // NOTE: this logic is specific to regular Nx by Ny by Nz grid
  // populate csurfs with same double loop

  for (m = 0; m < ncell; m++) count[m] = 0;
  
  for (m = 0; m < nsurf; m++) {
    if (dimension == 2) {
      x1 = pts[lines[m].p1].x;
      x2 = pts[lines[m].p2].x;
      x3 = x2;
    } else {
      x1 = pts[tris[m].p1].x;
      x2 = pts[tris[m].p2].x;
      x3 = pts[tris[m].p3].x;
    }

    lo[0] = MIN(x1[0],x2[0]);
    lo[0] = MIN(lo[0],x3[0]);
    hi[0] = MAX(x1[0],x2[0]);
    hi[0] = MAX(hi[0],x3[0]);

    lo[1] = MIN(x1[1],x2[1]);
    lo[1] = MIN(lo[1],x3[1]);
    hi[1] = MAX(x1[1],x2[1]);
    hi[1] = MAX(hi[1],x3[1]);

    lo[2] = MIN(x1[2],x2[2]);
    lo[2] = MIN(lo[2],x3[2]);
    hi[2] = MAX(x1[2],x2[2]);
    hi[2] = MAX(hi[2],x3[2]);

    ilo = MAX(0,static_cast<int> ((lo[0]-epsilon)*xdeltainv));
    ihi = MIN(nx-1,static_cast<int> ((hi[0]+epsilon)*xdeltainv));
    jlo = MAX(0,static_cast<int> ((lo[1]-epsilon)*ydeltainv));
    jhi = MIN(ny-1,static_cast<int> ((hi[1]+epsilon)*ydeltainv));
    klo = MAX(0,static_cast<int> ((lo[2]-epsilon)*zdeltainv));
    khi = MIN(nz-1,static_cast<int> ((hi[2]+epsilon)*zdeltainv));

    if (dimension == 2) {
      for (i = ilo; i <= ihi; i++)
	for (j = jlo; j <= jhi; j++) {
	  icell = j*nx + i;
	  if (Geometry::line_quad_intersect(x1,x2,lines[m].norm,
					    cells[icell].lo,cells[icell].hi))
	    csurfs[icell][count[icell]++] = m;
	}
    } else {
      for (i = ilo; i <= ihi; i++)
	for (j = jlo; j <= jhi; j++)
	  for (k = klo; k <= khi; k++) {
	    icell = k*nx*ny + j*nx + i;
	    if (Geometry::tri_hex_intersect(x1,x2,x3,tris[m].norm,
					    cells[icell].lo,cells[icell].hi))
	      csurfs[icell][count[icell]++] = m;
	  }
    }
  }

  // set per-cell surf count

  for (icell = 0; icell < ncell; icell++)
    cells[icell].nsurf = count[icell];

  // clean up
    
  memory->destroy(count);
}

/* ----------------------------------------------------------------------
   assign surface elements (lines or triangles) to grid cells
   NOTE: no parallelism yet, since each proc owns entire grid & all surfs
------------------------------------------------------------------------- */

void Grid::grid_inout()
{
  int icell;
  for (int m = 0; m < nlocal; m++) {
    icell = mycells[m];
    if (cells[icell].nsurf) cells[icell].inflag = SURFCONTAIN; 
    else cells[icell].inflag = SURFEXTERIOR; 
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
  bigint bytes = (bigint) ncell * sizeof(OneCell);   // cells
  bytes += nlocal*sizeof(int);                       // mycells

  if (surf->surf_exist) {                            // csurfs
    int n = 0;
    for (int i = 0; i < ncell; i++)
      n += cells[i].nsurf;
    bytes += n*sizeof(int);
  }
  
  return bytes;
}

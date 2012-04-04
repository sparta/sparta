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
#include "math_extra.h"
#include "update.h"
#include "domain.h"
#include "comm.h"
#include "surf.h"
#include "random_mars.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace DSMC_NS;
using namespace MathExtra;

#define DELTA 10000
#define EPSILON 1.0e-6

enum{XYZ,XZY,YXZ,YZX,ZXY,ZYX};
enum{NONE,OUTSIDE,INSIDE,BOTH};

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

  if (surf->npoint) {
    surf2grid();

    int icell;
    int stotal = 0;
    int smax = 0;
    for (int m = 0; m < nlocal; m++) {
      icell = mycells[m];
      stotal += cells[icell].nsurf;
      smax = MAX(smax,cells[icell].nsurf);
    }

    int stotalall,smaxall;
    MPI_Allreduce(&stotal,&stotalall,1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&smax,&smaxall,1,MPI_INT,MPI_MAX,world);

    if (comm->me == 0) {
      if (screen) {
	fprintf(screen,"Total surf elements in all cells = %d\n",stotalall);
	fprintf(screen,"Max surf elements in one cell = %d\n",smaxall);
      }
      if (logfile) {
	fprintf(logfile,"Total surf elements in all cells = %d\n",stotalall);
	fprintf(logfile,"Max surf elements in one cell = %d\n",smaxall);
      }
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

/* ---------------------------------------------------------------------- */

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
  // NOTE: no parallelism here yet, since each proc owns entire grid & surfs
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
	  if (line_quad_intersect(x1,x2,lines[m].norm,
				  cells[icell].lo,cells[icell].hi))
	    count[icell]++;
	}
    } else {
      for (i = ilo; i <= ihi; i++)
	for (j = jlo; j <= jhi; j++)
	  for (k = klo; k <= khi; k++) {
	    icell = k*nx*ny + j*nx + i;
	    if (tri_hex_intersect(x1,x2,x3,tris[m].norm,
				  cells[icell].lo,cells[icell].hi))
	      count[icell]++;
	  }
    }
  }

  // allocate ragged csurfs array
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
	  if (line_quad_intersect(x1,x2,lines[m].norm,
				  cells[icell].lo,cells[icell].hi))
	    csurfs[icell][count[icell]++] = m;
	}
    } else {
      for (i = ilo; i <= ihi; i++)
	for (j = jlo; j <= jhi; j++)
	  for (k = klo; k <= khi; k++) {
	    icell = k*nx*ny + j*nx + i;
	    if (tri_hex_intersect(x1,x2,x3,tris[m].norm,
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
   compute intersection of a line with a quad cell
   intersection is defined as
     any line pt (interior, vertex) in common with
     any rectangle pt (interior, edge, vertex)
   v0,v1 and norm = 2 vertices of line and unit normal vec
   lo,hi = opposite corner pts of rectangle
   return 1 if intersection, else 0
------------------------------------------------------------------------- */

int Grid::line_quad_intersect(double *v0, double *v1, double *norm,
			      double *lo, double *hi)
{
  double xlo,ylo,xhi,yhi,sum;
  double b[3],e[3],h0[3],h1[3],h2[3],h3[3],n[3],point[3];
  double param;
  int side;

  xlo = lo[0];
  ylo = lo[1];
  xhi = hi[0];
  yhi = hi[1];

  // if all 4 rectangle pts are on same side of line, no intersection

  sum = lineside(v0,norm,xlo,ylo);
  sum += lineside(v0,norm,xhi,ylo);
  sum += lineside(v0,norm,xlo,yhi);
  sum += lineside(v0,norm,xhi,yhi);
  
  if (sum == 4 || sum == -4) return 0;
	
  // if either of line vertices are inside quad, intersection
  // use <= and >= so touching hex surface is same as inside it

  if (v0[0] >= xlo && v0[0] <= xhi && v0[1] >= ylo && v0[1] <= yhi) return 1;
  if (v1[0] >= xlo && v1[0] <= xhi && v1[1] >= ylo && v1[1] <= yhi) return 1;

  // test 4 quad edges for intersection with line
  // b,e = begin/end of quad edge line segment
  // NOTE: need to write this method

  b[0] = xlo;   b[1] = ylo;   b[2] = 0.0;
  e[0] = xhi;   e[1] = ylo;   e[2] = 0.0;
  if (line_line_intersect(v0,v1,norm,b,e,point,param,side)) return 1;

  b[0] = xhi;   b[1] = ylo;   b[2] = 0.0;
  e[0] = xhi;   e[1] = yhi;   e[2] = 0.0;
  if (line_line_intersect(v0,v1,norm,b,e,point,param,side)) return 1;

  b[0] = xhi;   b[1] = yhi;   b[2] = 0.0;
  e[0] = xlo;   e[1] = yhi;   e[2] = 0.0;
  if (line_line_intersect(v0,v1,norm,b,e,point,param,side)) return 1;

  b[0] = xlo;   b[1] = yhi;   b[2] = 0.0;
  e[0] = xlo;   e[1] = ylo;   e[2] = 0.0;
  if (line_line_intersect(v0,v1,norm,b,e,point,param,side)) return 1;

  return 0;
}

/* ----------------------------------------------------------------------
   compute intersection of a triangle with a hex cell
   intersection is defined as
     any triangle pt (interior, edge, vertex) in common with
     any hex pt (interior, face, edge, vertex)
   v0,v1,v2 and norm = 3 vertices of triangle and unit normal vec
   lo,hi = opposite corner pts of hex
   return 1 if intersection, else 0
------------------------------------------------------------------------- */

int Grid::tri_hex_intersect(double *v0, double *v1, double *v2, double *norm,
			    double *lo, double *hi)
{
  double xlo,ylo,zlo,xhi,yhi,zhi,sum;
  double b[3],e[3],h0[3],h1[3],h2[3],h3[3],n[3],point[3];
  double param;
  int side;

  xlo = lo[0];
  ylo = lo[1];
  zlo = lo[2];
  xhi = hi[0];
  yhi = hi[1];
  zhi = hi[2];

  // if all 8 hex pts are on same side of tri plane, no intersection

  sum = triside(v0,norm,xlo,ylo,zlo);
  sum += triside(v0,norm,xhi,ylo,zlo);
  sum += triside(v0,norm,xlo,yhi,zlo);
  sum += triside(v0,norm,xhi,yhi,zlo);
  sum += triside(v0,norm,xlo,ylo,zhi);
  sum += triside(v0,norm,xhi,ylo,zhi);
  sum += triside(v0,norm,xlo,yhi,zhi);
  sum += triside(v0,norm,xhi,yhi,zhi);
  
  if (sum == 8 || sum == -8) return 0;
	
  // if any of 3 tri vertices are inside hex, intersection
  // use <= and >= so touching hex surface is same as inside it

  if (v0[0] >= xlo && v0[0] <= xhi && v0[1] >= ylo && v0[1] <= yhi &&
      v0[2] >= zlo && v0[2] <= zhi) return 1;

  if (v1[0] >= xlo && v1[0] <= xhi && v1[1] >= ylo && v1[1] <= yhi &&
      v1[2] >= zlo && v1[2] <= zhi) return 1;

  if (v2[0] >= xlo && v2[0] <= xhi && v2[1] >= ylo && v2[1] <= yhi &&
      v2[2] >= zlo && v2[2] <= zhi) return 1;

  // test 12 hex edges for intersection with tri
  // b,e = begin/end of hex edge line segment

  b[0] = xlo;   b[1] = ylo;   b[2] = zlo;
  e[0] = xhi;   e[1] = ylo;   e[2] = zlo;
  if (tri_line_intersect(v0,v1,v2,norm,b,e,point,param,side)) return 1;

  b[0] = xlo;   b[1] = yhi;   b[2] = zlo;
  e[0] = xhi;   e[1] = yhi;   e[2] = zlo;
  if (tri_line_intersect(v0,v1,v2,norm,b,e,point,param,side)) return 1;

  b[0] = xlo;   b[1] = ylo;   b[2] = zhi;
  e[0] = xhi;   e[1] = ylo;   e[2] = zhi;
  if (tri_line_intersect(v0,v1,v2,norm,b,e,point,param,side)) return 1;

  b[0] = xlo;   b[1] = yhi;   b[2] = zhi;
  e[0] = xhi;   e[1] = yhi;   e[2] = zhi;
  if (tri_line_intersect(v0,v1,v2,norm,b,e,point,param,side)) return 1;

  b[0] = xlo;   b[1] = ylo;   b[2] = zlo;
  e[0] = xlo;   e[1] = yhi;   e[2] = zlo;
  if (tri_line_intersect(v0,v1,v2,norm,b,e,point,param,side)) return 1;

  b[0] = xhi;   b[1] = ylo;   b[2] = zlo;
  e[0] = xhi;   e[1] = yhi;   e[2] = zlo;
  if (tri_line_intersect(v0,v1,v2,norm,b,e,point,param,side)) return 1;

  b[0] = xlo;   b[1] = ylo;   b[2] = zhi;
  e[0] = xlo;   e[1] = yhi;   e[2] = zhi;
  if (tri_line_intersect(v0,v1,v2,norm,b,e,point,param,side)) return 1;

  b[0] = xhi;   b[1] = ylo;   b[2] = zhi;
  e[0] = xhi;   e[1] = yhi;   e[2] = zhi;
  if (tri_line_intersect(v0,v1,v2,norm,b,e,point,param,side)) return 1;

  b[0] = xlo;   b[1] = ylo;   b[2] = zlo;
  e[0] = xlo;   e[1] = ylo;   e[2] = zhi;
  if (tri_line_intersect(v0,v1,v2,norm,b,e,point,param,side)) return 1;

  b[0] = xhi;   b[1] = ylo;   b[2] = zlo;
  e[0] = xhi;   e[1] = ylo;   e[2] = zhi;
  if (tri_line_intersect(v0,v1,v2,norm,b,e,point,param,side)) return 1;

  b[0] = xlo;   b[1] = yhi;   b[2] = zlo;
  e[0] = xlo;   e[1] = yhi;   e[2] = zhi;
  if (tri_line_intersect(v0,v1,v2,norm,b,e,point,param,side)) return 1;

  b[0] = xhi;   b[1] = yhi;   b[2] = zlo;
  e[0] = xhi;   e[1] = yhi;   e[2] = zhi;
  if (tri_line_intersect(v0,v1,v2,norm,b,e,point,param,side)) return 1;

  // test 3 tri edges for intersection with 6 faces of hex
  // h0,h1,h2,h3 = 4 corner pts of hex face
  // n = normal to xyz faces, depends on vertex ordering
  // each face is treated as 2 triangles -> 6 tests per face
  
  h0[0] = xlo;  h0[1] = ylo;  h0[2] = zlo;
  h1[0] = xlo;  h1[1] = yhi;  h1[2] = zlo;
  h2[0] = xlo;  h2[1] = yhi;  h2[2] = zhi;
  h3[0] = xlo;  h3[1] = ylo;  h3[2] = zhi;
  n[0]  = 1.0;  n[1]  = 0.0;  n[2]  = 0.0;
  
  if (tri_line_intersect(h0,h1,h2,n,v0,v1,point,param,side) ||
      tri_line_intersect(h0,h1,h2,n,v1,v2,point,param,side) ||
      tri_line_intersect(h0,h1,h2,n,v2,v0,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v0,v1,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v1,v2,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v2,v0,point,param,side)) return 1;

  h0[0] = h1[0] = h2[0] = h3[0] = xhi;

  if (tri_line_intersect(h0,h1,h2,n,v0,v1,point,param,side) ||
      tri_line_intersect(h0,h1,h2,n,v1,v2,point,param,side) ||
      tri_line_intersect(h0,h1,h2,n,v2,v0,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v0,v1,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v1,v2,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v2,v0,point,param,side)) return 1;

  h0[0] = xlo;  h0[1] = ylo;  h0[2] = zlo;
  h1[0] = xhi;  h1[1] = ylo;  h1[2] = zlo;
  h2[0] = xhi;  h2[1] = ylo;  h2[2] = zhi;
  h3[0] = xlo;  h3[1] = ylo;  h3[2] = zhi;
  n[0]  = 0.0;  n[1]  = -1.0;  n[2]  = 0.0;
  
  if (tri_line_intersect(h0,h1,h2,n,v0,v1,point,param,side) ||
      tri_line_intersect(h0,h1,h2,n,v1,v2,point,param,side) ||
      tri_line_intersect(h0,h1,h2,n,v2,v0,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v0,v1,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v1,v2,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v2,v0,point,param,side)) return 1;

  h0[1] = h1[1] = h2[1] = h3[1] = yhi;

  if (tri_line_intersect(h0,h1,h2,n,v0,v1,point,param,side) ||
      tri_line_intersect(h0,h1,h2,n,v1,v2,point,param,side) ||
      tri_line_intersect(h0,h1,h2,n,v2,v0,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v0,v1,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v1,v2,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v2,v0,point,param,side)) return 1;

  h0[0] = xlo;  h0[1] = ylo;  h0[2] = zlo;
  h1[0] = xhi;  h1[1] = ylo;  h1[2] = zlo;
  h2[0] = xhi;  h2[1] = yhi;  h2[2] = zlo;
  h3[0] = xlo;  h3[1] = yhi;  h3[2] = zlo;
  n[0]  = 0.0;  n[1]  = 0.0;  n[2]  = 1.0;
  
  if (tri_line_intersect(h0,h1,h2,n,v0,v1,point,param,side) ||
      tri_line_intersect(h0,h1,h2,n,v1,v2,point,param,side) ||
      tri_line_intersect(h0,h1,h2,n,v2,v0,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v0,v1,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v1,v2,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v2,v0,point,param,side)) return 1;

  h0[2] = h1[2] = h2[2] = h3[2] = zhi;
  
  if (tri_line_intersect(h0,h1,h2,n,v0,v1,point,param,side) ||
      tri_line_intersect(h0,h1,h2,n,v1,v2,point,param,side) ||
      tri_line_intersect(h0,h1,h2,n,v2,v0,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v0,v1,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v1,v2,point,param,side) ||
      tri_line_intersect(h0,h2,h3,n,v2,v0,point,param,side)) return 1;

  return 0;
}

/* ----------------------------------------------------------------------
   detect intersection between a line segment and a triangle
   intersection is defined as any line segment pt (including end pts)
     in common with any triangle pt (interior, edge, vertex)
   one exception is if both line end pts are in plane of triangle
     then is not an intersection
   v0,v1,v2 = 3 vertices of triangle
   norm = unit vector normal to triangle plane
     pointing OUTSIDE via right-hand rule
   start,stop = end points of directed line segment
   return TRUE if there is an intersection, else FALSE
   if TRUE also return:
     point = pt of intersection
     param = intersection pt is this fraction along segment (0-1 inclusive)
     side = OUTSIDE or INSIDE (enum value)
------------------------------------------------------------------------- */

bool Grid::line_line_intersect(double *v0, double *v1,
			      double *norm, double *start, double *stop,
			      double *point, double &param, int &side)
{
  return false;
}

/* ----------------------------------------------------------------------
   detect intersection between a line segment and a triangle
   intersection is defined as any line segment pt (including end pts)
     in common with any triangle pt (interior, edge, vertex)
   one exception is if both line end pts are in plane of triangle
     then is not an intersection
   v0,v1,v2 = 3 vertices of triangle
   norm = unit vector normal to triangle plane
     pointing OUTSIDE via right-hand rule
   start,stop = end points of directed line segment
   return TRUE if there is an intersection, else FALSE
   if TRUE also return:
     point = pt of intersection
     param = intersection pt is this fraction along segment (0-1 inclusive)
     side = OUTSIDE or INSIDE (enum value)
------------------------------------------------------------------------- */

bool Grid::tri_line_intersect(double *v0, double *v1, double *v2,
			      double *norm, double *start, double *stop,
			      double *point, double &param, int &side)
{
  double vec[3],start2stop[3],edge[3],pvec[3],xproduct[3];

  // if start,stop are on same side of triangle, no intersection
  // if start,stop are both in plane of triangle, no intersection

  MathExtra::sub3(start,v0,vec);
  double dotstart = MathExtra::dot3(norm,vec);
  MathExtra::sub3(stop,v0,vec);
  double dotstop = MathExtra::dot3(norm,vec);

  if (dotstart < 0.0 && dotstop < 0.0) return false;
  if (dotstart > 0.0 && dotstop > 0.0) return false;
  if (dotstart == 0.0 && dotstop == 0.0) return false;

  // param = parametric distance from start to stop
  //   at which tri plane is intersected
  // param = must be 0.0 to 1.0 inclusive, else no intersection

  MathExtra::sub3(v0,start,vec);
  MathExtra::sub3(stop,start,start2stop);
  param = MathExtra::dot3(norm,vec) / MathExtra::dot3(norm,start2stop);
  if (param < 0.0 || param > 1.0) return false;

  // point = intersection pt with plane of triangle

  point[0] = start[0] + param * start2stop[0];
  point[1] = start[1] + param * start2stop[1];
  point[2] = start[2] + param * start2stop[2];

  // test if intersection pt is inside triangle
  // edge = edge vector of triangle
  // pvec = vector from vertex to intersection point
  // xproduct = cross product of edge with pvec
  // if dot product of xproduct with norm < 0.0 for any of 3 edges,
  //   intersection point is outside tri

  MathExtra::sub3(v1,v0,edge);
  MathExtra::sub3(point,v0,pvec);
  MathExtra::cross3(edge,pvec,xproduct);
  if (MathExtra::dot3(xproduct,norm) < 0.0) return false;

  MathExtra::sub3(v2,v1,edge);
  MathExtra::sub3(point,v1,pvec);
  MathExtra::cross3(edge,pvec,xproduct);
  if (MathExtra::dot3(xproduct,norm) < 0.0) return false;

  MathExtra::sub3(v0,v2,edge);
  MathExtra::sub3(point,v2,pvec);
  MathExtra::cross3(edge,pvec,xproduct);
  if (MathExtra::dot3(xproduct,norm) < 0.0) return false;

  // there is a valid intersection inside the triangle
  // set side to INSIDE or OUTSIDE
  // if start point is inside or outside then is INSIDE or OUTSIDE
  // if particle started on surface, then is opposite of stop point

  if (dotstart < 0.0) side = INSIDE;
  else if (dotstart > 0.0) side = OUTSIDE;
  else if (dotstop > 0.0) side = INSIDE;
  else if (dotstop < 0.0) side = OUTSIDE;
  return true;
}

/* ----------------------------------------------------------------------
   determine which side of line the point x,y is on
   plane is defined by vertex pt v and unit normal vec
   return -1,0,1 for below,on,above line
------------------------------------------------------------------------- */

int Grid::lineside(double *v, double *norm, double x, double y)
{
  double vec[3];
  vec[0] = x - v[0];
  vec[1] = y - v[1];
  vec[2] = 0.0;

  double dotproduct = MathExtra::dot3(norm,vec);
  if (dotproduct < 0.0) return -1;
  else if (dotproduct > 0.0) return 1;
  else return 0;
}

/* ----------------------------------------------------------------------
   determine which side of plane the point x,y,z is on
   plane is defined by vertex pt v and unit normal vec
   return -1,0,1 for below,on,above plane
------------------------------------------------------------------------- */

int Grid::triside(double *v, double *norm, double x, double y, double z)
{
  double vec[3];
  vec[0] = x - v[0];
  vec[1] = y - v[1];
  vec[2] = z - v[2];

  double dotproduct = MathExtra::dot3(norm,vec);
  if (dotproduct < 0.0) return -1;
  else if (dotproduct > 0.0) return 1;
  else return 0;
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

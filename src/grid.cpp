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

#include "math.h"
#include "string.h"
#include "grid.h"
#include "geometry.h"
#include "update.h"
#include "domain.h"
#include "comm.h"
#include "irregular.h"
#include "surf.h"
#include "cut2d.h"
#include "random_mars.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

#define DELTA 10000
#define EPSILON 1.0e-6
#define BIG 1.0e20

enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};         // same as Domain
enum{XYZ,XZY,YXZ,YZX,ZXY,ZYX};

// cell is entirely outside/inside surfs or has any overlap with surfs

enum{CELLUNKNOWN,CELLOUTSIDE,CELLINSIDE,CELLOVERLAP};   // several files

// corner pt is outside/inside surfs or is on a surf

enum{CORNERUNKNOWN,CORNEROUTSIDE,CORNERINSIDE,CORNEROVERLAP};  // several files

/* ---------------------------------------------------------------------- */

Grid::Grid(SPARTA *sparta) : Pointers(sparta)
{
  grid_exist = 0;

  ncell = nsplit = maxcell = 0;
  cells = NULL;
  csurfs = NULL;
  csplits = NULL;

  nparent = 0;
  myparent = NULL;
  cflags = NULL;

  nchild = 0;
  mychild = NULL;
}

/* ---------------------------------------------------------------------- */

Grid::~Grid()
{
  memory->sfree(cells);
  memory->destroy(csurfs);
  memory->destroy(csplits);
  memory->destroy(myparent);
  memory->destroy(cflags);
  memory->destroy(mychild);
}

/* ---------------------------------------------------------------------- */

void Grid::init()
{
  int i,j,m,icell;

  // grid cell and surf geometry calculations

  if (comm->me == 0) {
    if (screen) fprintf(screen,"Grid/surf-element stats:\n");
    if (logfile) fprintf(logfile,"Grid/surf-element stats:\n");
  }
 
  // allocate cflags if necessary

  int ncorner = 4;
  if (domain->dimension == 3) ncorner = 8;

  if (cflags == NULL) memory->create(cflags,nparent,ncorner,"grid:cflags");

  // debug switch
  int old = 0;

  if (surf->changed) {
    surf->changed = 0;

    if (old || domain->dimension == 3) {
      surf2grid();
      surf2grid_stats();

      // unmark all local parent cells and corner flags

      for (m = 0; m < nparent; m++) {
        icell = myparent[m];
        cells[icell].type = CELLUNKNOWN;
        for (i = 0; i < ncorner; i++) cflags[m][i] = CORNERUNKNOWN;
        cells[icell].nsplit = 1;
        if (domain->dimension == 2)
          cells[icell].volume = (cells[icell].hi[0]-cells[icell].lo[0]) * 
            (cells[icell].hi[1]-cells[icell].lo[1]);
        else
          cells[icell].volume = (cells[icell].hi[0]-cells[icell].lo[0]) * 
            (cells[icell].hi[1]-cells[icell].lo[1]) *
            (cells[icell].hi[2]-cells[icell].lo[2]);
      }
      
      // set type = CELLOVERLAP for cells with surfaces
      // set cflags via all_cell_corner calls
      // some corner pts may be left unmarked as CORNERUNKNOWN
      
      for (m = 0; m < nparent; m++) {
        icell = myparent[m];
        if (cells[icell].nsurf) {
          cells[icell].type = CELLOVERLAP;
          if (domain->dimension == 2)
            surf->all_cell_corner_line(cells[icell].nsurf,csurfs[icell],
                                       cells[icell].lo,cells[icell].hi,
                                       cflags[m]);
          else if (domain->dimension == 3)
            surf->all_cell_corner_tri(cells[icell].nsurf,csurfs[icell],
                                      cells[icell].lo,cells[icell].hi,
                                      cflags[m]);
        }
      }

      assign_child();
      grid_inout();
      grid_check();

    } else {
      Cut2d *cut = new Cut2d(sparta);
      cut->surf2grid();
      surf2grid_stats();
      cut->split();
      delete cut;

      assign_child();
      grid_inout2();
      grid_check();
    }

  // if no surfs, mark owned cells and corner points as OUTSIDE
  // set volume of each cell

  } else if (!surf->exist) {
    int icell;
    for (m = 0; m < nparent; m++) {
      icell = myparent[m];
      cells[icell].type = CELLOUTSIDE;
      if (domain->dimension == 2)
        cells[icell].volume = (cells[icell].hi[0]-cells[icell].lo[0]) * 
          (cells[icell].hi[1]-cells[icell].lo[1]);
      else
        cells[icell].volume = (cells[icell].hi[0]-cells[icell].lo[0]) * 
          (cells[icell].hi[1]-cells[icell].lo[1]) *
          (cells[icell].hi[2]-cells[icell].lo[2]);
      for (i = 0; i < ncorner; i++) cflags[m][i] = CORNEROUTSIDE;
    }

    assign_child();
  }

  flow_stats();
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
  c->nsplit = 1;

  ncell++;
}

/* ---------------------------------------------------------------------- */

void Grid::add_split_cell(int icell)
{
  if (ncell+nsplit == maxcell) grow(1);
  memcpy(&cells[ncell+nsplit],&cells[icell],sizeof(OneCell));
  nsplit++;
}

/* ----------------------------------------------------------------------
   map a particle coordinate into a grid cell
   NOTE: not currently used, see loop option in CreateMolecules
   NOTE: what if particle is at upper boundary of domain
   NOTE: assumes Nx by Ny by Nz grid
------------------------------------------------------------------------- */

int Grid::which_cell(double x, double y, double z)
{
  // need to subtract off lo-corner of box to make this work
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

  assign_parent();
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

  assign_parent();
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

  assign_parent();
}

/* ----------------------------------------------------------------------
   create myparent list of parent grid cells indices I own
   also compute plocal index for owned cells
------------------------------------------------------------------------- */

void Grid::assign_parent()
{
  int me = comm->me;

  nparent = 0;
  for (int m = 0; m < ncell; m++)
    if (cells[m].proc == me) nparent++;

  memory->destroy(myparent);
  memory->create(myparent,nparent,"grid:myparent");

  nparent = 0;
  for (int m = 0; m < ncell; m++)
    if (cells[m].proc == me) {
      cells[m].plocal = nparent;
      myparent[nparent++] = m;
    }
}

/* ----------------------------------------------------------------------
   create myparent list of child grid cells indices I own
   also compute clocal index for owned cells
------------------------------------------------------------------------- */

void Grid::assign_child()
{
  int me = comm->me;

  nchild = 0;
  for (int m = 0; m < ncell+nsplit; m++)
    if (cells[m].proc == me && cells[m].nsplit <= 1) nchild++;

  memory->destroy(mychild);
  memory->create(mychild,nchild,"grid:mychild");

  nchild = 0;
  for (int m = 0; m < ncell+nsplit; m++)
    if (cells[m].proc == me && cells[m].nsplit <= 1) {
      cells[m].clocal = nchild;
      mychild[nchild++] = m;
    }
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
   each proc needs this info for all grid cells
   NOTE: no parallelism yet, but could compute on local cells and comm
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

  // tally count by double loop over surfs and grid cells within surf bbox
  // lo/hi = bounding box around surf
  // ijk lo/hi = grid index bounding box around surf
  // add epsilon to insure surf is counted in any cell it touches
  // icell = index of a grid cell within bounding box
  // NOTE: this logic is specific to regular Nx by Ny by Nz grid

  Surf::Point *pts = surf->pts;
  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  double xlo = domain->boxlo[0];
  double ylo = domain->boxlo[0];
  double zlo = domain->boxlo[0];

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

    ilo = MAX(0,static_cast<int> ((lo[0]-xlo-epsilon)*xdeltainv));
    ihi = MIN(nx-1,static_cast<int> ((hi[0]-xlo+epsilon)*xdeltainv));
    jlo = MAX(0,static_cast<int> ((lo[1]-ylo-epsilon)*ydeltainv));
    jhi = MIN(ny-1,static_cast<int> ((hi[1]-ylo+epsilon)*ydeltainv));
    klo = MAX(0,static_cast<int> ((lo[2]-zlo-epsilon)*zdeltainv));
    khi = MIN(nz-1,static_cast<int> ((hi[2]-zlo+epsilon)*zdeltainv));

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

    ilo = MAX(0,static_cast<int> ((lo[0]-xlo-epsilon)*xdeltainv));
    ihi = MIN(nx-1,static_cast<int> ((hi[0]-xlo+epsilon)*xdeltainv));
    jlo = MAX(0,static_cast<int> ((lo[1]-ylo-epsilon)*ydeltainv));
    jhi = MIN(ny-1,static_cast<int> ((hi[1]-ylo+epsilon)*ydeltainv));
    klo = MAX(0,static_cast<int> ((lo[2]-zlo-epsilon)*zdeltainv));
    khi = MIN(nz-1,static_cast<int> ((hi[2]-zlo+epsilon)*zdeltainv));

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

  for (icell = 0; icell < ncell; icell++) {
    cells[icell].nsurf = count[icell];
    //int ix = icell % nx;
    //int iy = (icell/nx) % ny;
    //int iz = icell / (nx*ny);
    //printf("AAA %d icell, %d %d %d ijk, %d count\n",
    //      icell,ix+1,iy+1,iz+1,count[icell]);
  }

  // clean up
    
  memory->destroy(count);
}

/* ----------------------------------------------------------------------
   set type and cflags of all owned cells
------------------------------------------------------------------------- */

void Grid::grid_inout()
{
  int i,j,m,n,ipt,jpt,ivalue,jvalue;
  int icell,ilocal,iface,ncorner,nface_pts;
  int mark,progress,progress_any;
  int *neigh;
  int faceflip[6] = {XHI,XLO,YHI,YLO,ZHI,ZLO};

  // corners[i][j] = J corner points of face I of a grid cell
  // works for 2d quads and 3d hexes
  
  int corners[6][4] = {{0,2,4,6}, {1,3,5,7}, {0,1,4,5}, {2,3,6,7}, 
		       {0,1,2,3}, {4,5,6,7}};

  int me = comm->me;
  int dimension = domain->dimension;
  if (dimension == 3) {
    ncorner = 8;
    nface_pts = 4;
  } else {
    ncorner = 4;
    nface_pts = 2;
  }

  // create irregular communicator for exchanging cell corner flags
  // allocate sbuf/rbuf one larger so can access buf[0][0] even if no send/recv
  // NOTE: could put this logic into Comm class to allow other cell-based
  //       comm with methods that use pack/unpack callbacks to this class

  int nsend = 0;
  for (m = 0; m < nparent; m++) {
    icell = myparent[m];
    neigh = cells[icell].neigh;
    for (j = 0; j < 6; j++) {
      if (neigh[j] < 0) continue;
      if (cells[neigh[j]].proc != me) nsend++;
    }
  }

  int *proclist;
  memory->create(proclist,nsend,"grid:proclist");

  nsend = 0;
  for (m = 0; m < nparent; m++) {
    icell = myparent[m];
    neigh = cells[icell].neigh;
    for (j = 0; j < 6; j++) {
      if (neigh[j] < 0) continue;
      if (cells[neigh[j]].proc != me) proclist[nsend++] = cells[neigh[j]].proc;
    }
  }

  Irregular *irregular = new Irregular(sparta);
  int nrecv = irregular->create(nsend,proclist);
  int **sbuf,**rbuf;
  memory->create(sbuf,nsend+1,2+ncorner,"grid:sbuf");
  memory->create(rbuf,nrecv+1,2+ncorner,"grid:rbuf");

  // loop until make no more progress marking corner flags
  // 3 ways of marking corner points of my cells:
  // a) from image points in other cells I own
  // b) via corner line/tri which shoots lines to other marked pts in cell
  // c) from image points in cells owned by other procs via irregular comm

  while (1) {
    progress = 0;

    // loop over my cells and their faces
    // skip OUTSIDE/INSIDE cells since already fully marked
    // compare each face corner pt to its image pt in owned adjacent cells
    // if neither marked, continue
    // if both marked, check consistency
    // if only one is marked, mark the other with the same value
    // if mark a pt as OUTSIDE/INSIDE and cell type is unset,
    //   mark all pts in cell and set cell type
    // repeat until can mark no more corner points

    while (1) {
      mark = 0;
      for (m = 0; m < nparent; m++) {
	icell = myparent[m];
	if (cells[icell].type == CELLOUTSIDE || 
	    cells[icell].type == CELLINSIDE) continue;

	neigh = cells[icell].neigh;
	for (j = 0; j < 6; j++) {
	  if (neigh[j] < 0) continue;
	  if (cells[neigh[j]].proc != me) continue;

	  for (i = 0; i < nface_pts; i++) {
	    ipt = corners[j][i];
	    jpt = corners[faceflip[j]][i];
	    ivalue = cflags[m][ipt];
	    jvalue = cflags[cells[neigh[j]].plocal][jpt];

	    if (ivalue != CORNERUNKNOWN && jvalue != CORNERUNKNOWN) {
	      if (ivalue != jvalue) {
		char str[128];
		sprintf(str,"Mismatched owned corner flags %d %d "
			"in two cells %d %d at corners %d %d",
			ivalue,jvalue,icell,neigh[j],ipt,jpt);
		error->one(FLERR,str);
	      }
	    } else if (ivalue != CORNERUNKNOWN) {
	      cflags[cells[neigh[j]].plocal][jpt] = ivalue;
	      if (cells[neigh[j]].type == CELLUNKNOWN && 
		  (ivalue == CORNEROUTSIDE || ivalue == CORNERINSIDE)) {
		if (flood(ivalue,ncorner,cells[neigh[j]].plocal))
		  error->one(FLERR,"Mismatched corner flags within one cell");
                if (ivalue == CORNEROUTSIDE)
                  cells[neigh[j]].type = CELLOUTSIDE;
                else if (ivalue == CORNERINSIDE)
                  cells[neigh[j]].type = CELLINSIDE;
	      }
	      mark = 1;
	    } else if (jvalue != CORNERUNKNOWN) {
	      cflags[m][ipt] = jvalue;
	      if (cells[icell].type == CELLUNKNOWN && 
		  (jvalue == CORNEROUTSIDE || jvalue == CORNERINSIDE)) {
		if (flood(jvalue,ncorner,m)) 
		  error->one(FLERR,"Mismatched corner flags within one cell");
                if (jvalue == CORNEROUTSIDE)
                  cells[icell].type = CELLOUTSIDE;
                else if (jvalue == CORNERINSIDE)
                  cells[icell].type = CELLINSIDE;
	      }
	      mark = 1;
	    }
	  }	  
	}
      }

      if (mark) progress = 1;
      if (!mark) break;
    }

    // look for unmarked corner pts in CELLOVERLAP cells
    // attempt to mark them via one_cell_corner calls

    for (m = 0; m < nparent; m++) {
      icell = myparent[m];
      if (cells[icell].type != CELLOVERLAP) continue;
      for (i = 0; i < ncorner; i++) {
	if (cflags[m][i] != CORNERUNKNOWN) continue;
	if (dimension == 2)
	  cflags[m][i] = 
	    surf->one_cell_corner_line(i,cells[icell].nsurf,csurfs[icell],
				       cells[icell].lo,cells[icell].hi,
				       cflags[m]);
	else if (dimension == 3) {
	  cflags[m][i] = 
	    surf->one_cell_corner_tri(i,cells[icell].nsurf,csurfs[icell],
				      cells[icell].lo,cells[icell].hi,
				      cflags[m]);
	}
	if (cflags[m][i] != CORNERUNKNOWN) progress = 1;
      }
    }

    // communicate cflags to neighbors of my cells owned by other procs
    // send info = icell of neighbor, iface of neighbor, cflags for my cell
    // NOTE: sending icell assumes all procs have global cell list

    nsend = 0;
    for (m = 0; m < nparent; m++) {
      neigh = cells[myparent[m]].neigh;
      for (j = 0; j < 6; j++) {
	if (neigh[j] < 0) continue;
	if (cells[neigh[j]].proc != me) {
	  sbuf[nsend][0] = neigh[j];
	  sbuf[nsend][1] = faceflip[j];
	  for (n = 0; n < ncorner; n++)
	    sbuf[nsend][2+n] = cflags[m][n];
	  nsend++;
	}
      }
    }

    // perform irregular comm of cflags info

    irregular->exchange((char *) &sbuf[0][0],(2+ncorner)*sizeof(int),
    			(char *) &rbuf[0][0]);

    // use received corner pt values to update cflags of my owned cells
    // if both marked, check consistency
    // if only other is marked, mark my pt with same value
    // if mark a pt as EXTERIOR/INTERIOR and cell type is unset,
    //   mark all pts in cell and set cell type

    for (m = 0; m < nrecv; m++) {
      icell = rbuf[m][0];
      iface = rbuf[m][1];
      ilocal = cells[icell].plocal;

      for (i = 0; i < nface_pts; i++) {
	ipt = corners[iface][i];
	jpt = corners[faceflip[iface]][i];
	ivalue = cflags[ilocal][ipt];
	jvalue = rbuf[m][2+jpt];
	if (ivalue != CORNERUNKNOWN && jvalue != CORNERUNKNOWN) {
	  if (ivalue != jvalue) {
	    char str[128];
	    sprintf(str,"Mismatched owned/ghost corner flags %d %d "
		    "in two cells %d %d at corners %d %d",ivalue,jvalue,
		    icell,cells[icell].neigh[iface],ipt,jpt);
	    error->one(FLERR,str);
	  }
	} else if (jvalue != CORNERUNKNOWN) {
	  cflags[ilocal][ipt] = jvalue;
	  if (cells[icell].type == CELLUNKNOWN && 
	      (jvalue == CORNEROUTSIDE || jvalue == CORNERINSIDE)) {
	    if (flood(jvalue,ncorner,ilocal))
	      error->one(FLERR,"Mismatched corner flags within one cell");
            if (jvalue == CORNEROUTSIDE)
              cells[icell].type = CELLOUTSIDE;
            else if (jvalue == CORNERINSIDE)
              cells[icell].type = CELLINSIDE;
	  }
	  progress = 1;
	}
      }	  
    }
    
    // if no new marks made by any processor, done

    MPI_Allreduce(&progress,&progress_any,1,MPI_INT,MPI_SUM,world);
    if (!progress_any) break;
  }

  // clean up

  delete irregular;
  memory->destroy(proclist);
  memory->destroy(sbuf);
  memory->destroy(rbuf);
}

/* ----------------------------------------------------------------------
   set type and cflags of all owned cells
------------------------------------------------------------------------- */

void Grid::grid_inout2()
{
  int i,j,m,n,ipt,jpt,ivalue,jvalue;
  int icell,ilocal,iface,ncorner,nface_pts;
  int mark,progress,progress_any;
  int *neigh;
  int faceflip[6] = {XHI,XLO,YHI,YLO,ZHI,ZLO};

  // corners[i][j] = J corner points of face I of a grid cell
  // works for 2d quads and 3d hexes
  
  int corners[6][4] = {{0,2,4,6}, {1,3,5,7}, {0,1,4,5}, {2,3,6,7}, 
		       {0,1,2,3}, {4,5,6,7}};

  int me = comm->me;
  int dimension = domain->dimension;
  if (dimension == 3) {
    ncorner = 8;
    nface_pts = 4;
  } else {
    ncorner = 4;
    nface_pts = 2;
  }

  // create irregular communicator for exchanging cell corner flags
  // allocate sbuf/rbuf one larger so can access buf[0][0] even if no send/recv
  // NOTE: could put this logic into Comm class to allow other cell-based
  //       comm with methods that use pack/unpack callbacks to this class

  int nsend = 0;
  for (m = 0; m < nparent; m++) {
    icell = myparent[m];
    neigh = cells[icell].neigh;
    for (j = 0; j < 6; j++) {
      if (neigh[j] < 0) continue;
      if (cells[neigh[j]].proc != me) nsend++;
    }
  }

  int *proclist;
  memory->create(proclist,nsend,"grid:proclist");

  nsend = 0;
  for (m = 0; m < nparent; m++) {
    icell = myparent[m];
    neigh = cells[icell].neigh;
    for (j = 0; j < 6; j++) {
      if (neigh[j] < 0) continue;
      if (cells[neigh[j]].proc != me) proclist[nsend++] = cells[neigh[j]].proc;
    }
  }

  Irregular *irregular = new Irregular(sparta);
  int nrecv = irregular->create(nsend,proclist);
  int **sbuf,**rbuf;
  memory->create(sbuf,nsend+1,2+ncorner,"grid:sbuf");
  memory->create(rbuf,nrecv+1,2+ncorner,"grid:rbuf");

  // loop until make no more progress marking corner and cell flags I own
  // 2 ways of marking corner points of my cells:
  // a) from image points in other cells I own
  // b) from image points in cells owned by other procs via irregular comm

  while (1) {
    progress = 0;

    // loop over my cells and their faces
    // skip OUTSIDE/INSIDE cells since already fully marked
    // compare each face corner pt to its image pt in owned adjacent cells
    // if both UNKNOWN marked, continue
    // if both marked, check consistency
    // if only one is marked, mark the other with the same value
    // if mark a pt as OUTSIDE/INSIDE and cell type is UNKNOWN,
    //   mark all corner pts in cell and set cell type
    // repeat until can mark no more corner points

    while (1) {
      mark = 0;
      for (m = 0; m < nparent; m++) {
	icell = myparent[m];
	if (cells[icell].type == CELLOUTSIDE || 
	    cells[icell].type == CELLINSIDE) continue;

	neigh = cells[icell].neigh;
	for (j = 0; j < 6; j++) {
	  if (neigh[j] < 0) continue;
	  if (cells[neigh[j]].proc != me) continue;

	  for (i = 0; i < nface_pts; i++) {
	    ipt = corners[j][i];
	    jpt = corners[faceflip[j]][i];
	    ivalue = cflags[m][ipt];
	    jvalue = cflags[cells[neigh[j]].plocal][jpt];

	    if (ivalue != CORNERUNKNOWN && jvalue != CORNERUNKNOWN) {
	      if (ivalue != jvalue) {
		char str[128];
		sprintf(str,"Mismatched owned corner flags %d %d "
			"in two cells %d %d at corners %d %d",
			ivalue,jvalue,icell,neigh[j],ipt,jpt);
		error->one(FLERR,str);
	      }
	    } else if (ivalue != CORNERUNKNOWN) {
	      cflags[cells[neigh[j]].plocal][jpt] = ivalue;
	      if (cells[neigh[j]].type == CELLUNKNOWN && 
		  (ivalue == CORNEROUTSIDE || ivalue == CORNERINSIDE)) {
		if (flood(ivalue,ncorner,cells[neigh[j]].plocal))
		  error->one(FLERR,"Mismatched corner flags within one cell");
                if (ivalue == CORNEROUTSIDE)
                  cells[neigh[j]].type = CELLOUTSIDE;
                else if (ivalue == CORNERINSIDE)
                  cells[neigh[j]].type = CELLINSIDE;
	      }
	      mark = 1;
	    } else if (jvalue != CORNERUNKNOWN) {
	      cflags[m][ipt] = jvalue;
	      if (cells[icell].type == CELLUNKNOWN && 
		  (jvalue == CORNEROUTSIDE || jvalue == CORNERINSIDE)) {
		if (flood(jvalue,ncorner,m)) 
		  error->one(FLERR,"Mismatched corner flags within one cell");
                if (jvalue == CORNEROUTSIDE)
                  cells[icell].type = CELLOUTSIDE;
                else if (jvalue == CORNERINSIDE)
                  cells[icell].type = CELLINSIDE;
	      }
	      mark = 1;
	    }
	  }	  
	}
      }

      if (mark) progress = 1;
      if (!mark) break;
    }

    // communicate cflags to neighbors of my cells owned by other procs
    // send info = icell of neighbor, iface of neighbor, cflags for my cell
    // NOTE: sending icell assumes all procs have global cell list

    nsend = 0;
    for (m = 0; m < nparent; m++) {
      neigh = cells[myparent[m]].neigh;
      for (j = 0; j < 6; j++) {
	if (neigh[j] < 0) continue;
	if (cells[neigh[j]].proc != me) {
	  sbuf[nsend][0] = neigh[j];
	  sbuf[nsend][1] = faceflip[j];
	  for (n = 0; n < ncorner; n++)
	    sbuf[nsend][2+n] = cflags[m][n];
	  nsend++;
	}
      }
    }

    // perform irregular comm of cflags info

    irregular->exchange((char *) &sbuf[0][0],(2+ncorner)*sizeof(int),
    			(char *) &rbuf[0][0]);

    // use received corner pt values to update cflags of my owned cells
    // if both marked, check consistency
    // if only other is marked, mark my pt with same value
    // if mark a pt as OUTSIDE/INSIDE and cell type is UNKNOWN,
    //   mark all pts in cell and set cell type

    for (m = 0; m < nrecv; m++) {
      icell = rbuf[m][0];
      iface = rbuf[m][1];
      ilocal = cells[icell].plocal;

      for (i = 0; i < nface_pts; i++) {
	ipt = corners[iface][i];
	jpt = corners[faceflip[iface]][i];
	ivalue = cflags[ilocal][ipt];
	jvalue = rbuf[m][2+jpt];
	if (ivalue != CORNERUNKNOWN && jvalue != CORNERUNKNOWN) {
	  if (ivalue != jvalue) {
	    char str[128];
	    sprintf(str,"Mismatched owned/ghost corner flags %d %d "
		    "in two cells %d %d at corners %d %d",ivalue,jvalue,
		    icell,cells[icell].neigh[iface],ipt,jpt);
	    error->one(FLERR,str);
	  }
	} else if (jvalue != CORNERUNKNOWN) {
	  cflags[ilocal][ipt] = jvalue;
	  if (cells[icell].type == CELLUNKNOWN && 
	      (jvalue == CORNEROUTSIDE || jvalue == CORNERINSIDE)) {
	    if (flood(jvalue,ncorner,ilocal))
	      error->one(FLERR,"Mismatched corner flags within one cell");
            if (jvalue == CORNEROUTSIDE)
              cells[icell].type = CELLOUTSIDE;
            else if (jvalue == CORNERINSIDE)
              cells[icell].type = CELLINSIDE;
	  }
	  progress = 1;
	}
      }	  
    }
    
    // if no new marks made by any processor, done

    MPI_Allreduce(&progress,&progress_any,1,MPI_INT,MPI_SUM,world);
    if (!progress_any) break;
  }

  // clean up

  delete irregular;
  memory->destroy(proclist);
  memory->destroy(sbuf);
  memory->destroy(rbuf);
}

/* ----------------------------------------------------------------------
   set all corner flags of local cell M to value
   return 1 as error if any flags are already set to a different value
   return 0 if success
------------------------------------------------------------------------- */

int Grid::flood(int value, int ncorner, int m)
{
  for (int i = 0; i < ncorner; i++) {
    if (cflags[m][i] != CORNERUNKNOWN && cflags[m][i] != value) return 1;
    cflags[m][i] = value;
  }
  return 0;
} 

/* ----------------------------------------------------------------------
   require all cell types be set
   check corner pts of cells with surfs
   require corner pts on global box boundary be set
   warn if interior corner pts are not set
------------------------------------------------------------------------- */

void Grid::grid_check()
{
  int i,m,icell;

  // check cell types

  int unknown = 0;
  for (m = 0; m < nparent; m++) {
    icell = myparent[m];
    if (cells[icell].type == CELLUNKNOWN) unknown++; 
  }
  int unknownall;
  MPI_Allreduce(&unknown,&unknownall,1,MPI_INT,MPI_SUM,world);

  if (unknownall) {
    char str[128];
    sprintf(str,"Grid cells marked as unknown = %d",unknownall);
    error->all(FLERR,str);
  }

  // check corner flags of cells that are CELLOVERLAP
  // warn if any interior corner flags are not set
  // error if any corner flags on global boundaries are unset

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  int dimension = domain->dimension;

  int ncorner = 4;
  if (dimension == 3) ncorner = 8;

  double x[3];
  int inside = 0;
  int outside = 0;

  for (m = 0; m < nparent; m++) {
    icell = myparent[m];
    if (cells[icell].type != CELLOVERLAP) continue;
    for (i = 0; i < ncorner; i++) {
      if (cflags[m][i] != CORNERUNKNOWN) continue;
      if (i % 2 == 0) x[0] = cells[icell].lo[0];
      else x[0] = cells[icell].hi[0];
      if ((i/2) % 2 == 0) x[1] = cells[icell].lo[1];
      else x[1] = cells[icell].hi[1];
      if (dimension == 3) {
        if (i/4 == 0) x[2] = cells[icell].lo[2];
        else x[2] = cells[icell].hi[2];
      } else x[2] = 0.0;

      if (Geometry::point_on_hex(x,boxlo,boxhi)) outside++;
      else inside++;
    }
  }

  int insideall;
  MPI_Allreduce(&inside,&insideall,1,MPI_INT,MPI_SUM,world);
  if (insideall) {
    char str[128];
    sprintf(str,"Grid cell interior corner points marked as unknown = %d",
	    insideall);
    if (comm->me == 0) error->warning(FLERR,str);
  }

  int outsideall;
  MPI_Allreduce(&outside,&outsideall,1,MPI_INT,MPI_SUM,world);
  if (outsideall) {
    char str[128];
    sprintf(str,"Grid cell boundary corner points marked as unknown = %d",
	    outsideall);
    error->all(FLERR,str);
  }
}

/* ---------------------------------------------------------------------- */

void Grid::surf2grid_stats()
{
  int i,m;
  double cmax,len,area;
  int dimension = domain->dimension;

  int icell;
  int stotal = 0;
  int smax = 0;
  double sratio = BIG;
  for (int m = 0; m < nparent; m++) {
    icell = myparent[m];
    stotal += cells[icell].nsurf;
    smax = MAX(smax,cells[icell].nsurf);
    
    cmax = MAX(cells[icell].hi[0] - cells[icell].lo[0],
	       cells[icell].hi[1] - cells[icell].lo[1]);
    if (dimension == 3) 
      cmax = MAX(cmax,cells[icell].hi[2] - cells[icell].lo[2]);
    
    if (dimension == 2) {
      for (int i = 0; i < cells[icell].nsurf; i++) {
	len = surf->line_size(csurfs[icell][i]);
	sratio = MIN(sratio,len/cmax);
      }
    } else if (dimension == 3) {
      for (int i = 0; i < cells[icell].nsurf; i++) {
	area = surf->tri_size(csurfs[icell][i],len);
	sratio = MIN(sratio,len/cmax);
      }
    }
  }
  
  int stotalall,smaxall;
  double sratioall;
  MPI_Allreduce(&stotal,&stotalall,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&smax,&smaxall,1,MPI_INT,MPI_MAX,world);
  MPI_Allreduce(&sratio,&sratioall,1,MPI_DOUBLE,MPI_MIN,world);
  
  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"  %d = total surfs in all grid cells\n",stotalall);
      fprintf(screen,"  %d = max surfs in one grid cell\n",smaxall);
      fprintf(screen,"  %g = min surf-size/cell-size ratio\n",sratioall);
    }
    if (logfile) {
      fprintf(logfile,"  %d = total surfs in all grid cells\n",stotalall);
      fprintf(logfile,"  %d = max surfs in one grid cell\n",smaxall);
      fprintf(logfile,"  %g = min surf-size/cell-size ratio\n",sratioall);
    }
  }
}

/* ---------------------------------------------------------------------- */

void Grid::flow_stats()
{
  int i,icell;

  int outside = 0;
  int inside = 0;
  int overlap = 0;
  int maxsplit = 0;
  double cellvolume = 0.0;

  for (int m = 0; m < nparent; m++) {
    icell = myparent[m];
    if (cells[icell].type == CELLOUTSIDE) outside++;
    else if (cells[icell].type == CELLINSIDE) inside++;
    else if (cells[icell].type == CELLOVERLAP) overlap++;
    maxsplit = MAX(maxsplit,cells[icell].nsplit);
  }

  for (int m = 0; m < nchild; m++) {
    icell = mychild[m];
    if (cells[icell].type != CELLINSIDE) cellvolume += cells[icell].volume;
  }

  int outall,inall,overall,maxsplitall;
  double cellvolumeall;
  MPI_Allreduce(&outside,&outall,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&inside,&inall,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&overlap,&overall,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&maxsplit,&maxsplitall,1,MPI_INT,MPI_MAX,world);
  MPI_Allreduce(&cellvolume,&cellvolumeall,1,MPI_DOUBLE,MPI_SUM,world);

  double flowvolume = flow_volume();

  int *tally = new int[maxsplitall];
  int *tallyall = new int[maxsplitall];
  for (i = 0; i < maxsplitall; i++) tally[i] = 0;
  for (int m = 0; m < nparent; m++) {
    icell = myparent[m];
    if (cells[icell].type == CELLOVERLAP) tally[cells[icell].nsplit-1]++;
  }

  MPI_Allreduce(tally,tallyall,maxsplitall,MPI_INT,MPI_SUM,world);

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"  %d %d %d = cells outside/inside/overlapping surfs\n",
	      outall,inall,overall);
      fprintf(screen," ");
      for (i = 0; i < maxsplitall; i++) fprintf(screen," %d",tallyall[i]);
      fprintf(screen," = surf cells with 1,2,etc splits\n");
      fprintf(screen,"  %g %g = cell-wise and global flow volume\n",
              cellvolumeall,flowvolume);
    }
    if (logfile) {
      fprintf(logfile,"  %d %d %d = cells outside/inside/overlapping surfs\n",
	      outall,inall,overall);
      fprintf(logfile," ");
      for (i = 0; i < maxsplitall; i++) fprintf(logfile," %d",tallyall[i]);
      fprintf(logfile," = surf cells with 1,2,etc splits\n");
      fprintf(logfile,"  %g %g = cell-wise and global flow volume\n",
              cellvolumeall,flowvolume);
    }
  }

  delete [] tally;
  delete [] tallyall;
}

/* ---------------------------------------------------------------------- */

double Grid::flow_volume()
{
  double zarea;
  double *p1,*p2,*p3;

  Surf::Point *pts = surf->pts;
  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  double vol = 0.0;

  if (domain->dimension == 2) {
    for (int i = 0; i < surf->nline; i++) {
      p1 = pts[lines[i].p1].x;
      p2 = pts[lines[i].p2].x;
      if (p1[0] < p2[0]) vol -= (0.5*(p1[1]+p2[1]) - boxlo[1]) * (p2[0]-p1[0]);
      else vol += (0.5*(p1[1]+p2[1]) - boxlo[1]) * (p1[0]-p2[0]);
    }
    if (vol <= 0.0) vol += (boxhi[0]-boxlo[0]) * (boxhi[1]-boxlo[1]); 

  } else {
    for (int i = 0; i < surf->ntri; i++) {
      p1 = pts[tris[i].p1].x;
      p2 = pts[tris[i].p2].x;
      p3 = pts[tris[i].p3].x;
      zarea = 0.5 * ((p2[0]-p1[0])*(p3[1]-p1[1]) - (p2[1]-p1[1])*(p3[0]-p1[0]));
      vol -= zarea * ((p1[2]+p2[2]+p3[2])/3.0 - boxlo[2]);
    }
    if (vol <= 0.0) vol += (boxhi[0]-boxlo[0]) * (boxhi[1]-boxlo[1]) * 
                      (boxhi[2]-boxlo[2]); 
  }

  return vol;
}

/* ----------------------------------------------------------------------
   insure cell list can hold nextra new grid cells
------------------------------------------------------------------------- */

void Grid::grow(int nextra)
{
  bigint target = (bigint) ncell + nsplit + nextra;
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
  bigint bytes = ((bigint) ncell + nsplit) * sizeof(OneCell);    // cells
  bytes += nparent*sizeof(int);                    // myparent
  bytes += nchild*sizeof(int);                     // mychild

  int ncorners = 4;                                // cflags
  if (domain->dimension == 3) ncorners = 8;
  bytes += nparent*ncorners*sizeof(int);

  if (surf->exist) {                               // csurfs
    int n = 0;
    for (int i = 0; i < ncell; i++)
      n += cells[i].nsurf;
    bytes += n*sizeof(int);
  }
  
  return bytes;
}

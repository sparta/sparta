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

#include "surf.h"
#include "math_extra.h"
#include "style_surf_collide.h"
#include "surf_collide.h"
#include "domain.h"
#include "comm.h"
#include "geometry.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

#define DELTA 4
#define EPSSQ 1.0e-12

enum{CORNERUNKNOWN,CORNEROUTSIDE,CORNERINSIDE,CORNEROVERLAP};  // same as Grid
enum{OUTSIDE,INSIDE,ONSURF2OUT,ONSURF2IN};      // same as Update

/* ---------------------------------------------------------------------- */

Surf::Surf(SPARTA *sparta) : Pointers(sparta)
{
  exist = 0;
  changed = 0;

  npoint = nline = ntri = 0;
  pts = NULL;
  lines = NULL;
  tris = NULL;

  nlocal = 0;
  mysurfs = NULL;

  nsc = maxsc = 0;
  sc = NULL;
}

/* ---------------------------------------------------------------------- */

Surf::~Surf()
{
  memory->sfree(pts);
  memory->sfree(lines);
  memory->sfree(tris);

  memory->sfree(mysurfs);

  for (int i = 0; i < nsc; i++) delete sc[i];
  memory->sfree(sc);
}

/* ---------------------------------------------------------------------- */

void Surf::init()
{
  // initialize surface collision models

  for (int i = 0; i < nsc; i++) sc[i]->init();
}

/* ----------------------------------------------------------------------
   return # of lines or triangles
------------------------------------------------------------------------- */

int Surf::nelement()
{
  if (domain->dimension == 2) return nline;
  return ntri;
}

/* ----------------------------------------------------------------------
   setup owned surf elements
   create mysurfs list of owned surfs
   compute local index for owned cells
------------------------------------------------------------------------- */

void Surf::setup_surf()
{
  int me = comm->me;
  int nprocs = comm->nprocs;

  int n = nelement();

  // assign every Pth surf element to this proc

  nlocal = n/nprocs;
  if (me < n % nprocs) nlocal++;

  memory->destroy(mysurfs);
  memory->create(mysurfs,nlocal,"surf:mysurfs");

  nlocal = 0;
  for (int m = me; m < n; m += nprocs)
    mysurfs[nlocal++] = m;
}

/* ----------------------------------------------------------------------
   compute unit outward normal vectors of N lines starting at Nstart
   outward normal = +z axis x (p2-p1)
------------------------------------------------------------------------- */

void Surf::compute_line_normal(int nstart, int n)
{
  int p1,p2;
  double z[3],delta[3],norm[3];

  z[0] = 0.0; z[1] = 0.0; z[2] = 1.0;

  int m = nstart;
  for (int i = 0; i < n; i++) {
    p1 = lines[m].p1;
    p2 = lines[m].p2;
    MathExtra::sub3(pts[p2].x,pts[p1].x,delta);
    MathExtra::cross3(z,delta,norm);
    MathExtra::norm3(norm);
    lines[m].norm[0] = norm[0];
    lines[m].norm[1] = norm[1];
    lines[m].norm[2] = 0.0;
    m++;
  }
}

/* ----------------------------------------------------------------------
   compute unit outward normal vectors of N triangles starting at Nstart
   outward normal = (p2-p1) x (p3-p1)
------------------------------------------------------------------------- */

void Surf::compute_tri_normal(int nstart, int n)
{
  int p1,p2,p3;
  double delta12[3],delta13[3];

  int m = nstart;
  for (int i = 0; i < n; i++) {
    p1 = tris[m].p1;
    p2 = tris[m].p2;
    p3 = tris[m].p3;
    MathExtra::sub3(pts[p2].x,pts[p1].x,delta12);
    MathExtra::sub3(pts[p3].x,pts[p1].x,delta13);
    MathExtra::cross3(delta12,delta13,tris[m].norm);
    MathExtra::norm3(tris[m].norm);

    m++;
  }
}

/* ----------------------------------------------------------------------
   return coords of a corner point in a 2d quad
   icorner pts 1 to 4 are ordered by x, then by y
------------------------------------------------------------------------- */

void Surf::quad_corner_point(int icorner, double *lo, double *hi, double *pt)
{
  if (icorner % 2) pt[0] = hi[0];
  else pt[0] = lo[0];
  if (icorner / 2) pt[1] = hi[1];
  else pt[1] = lo[1];
  pt[2] = 0.0;
}

/* ----------------------------------------------------------------------
   return coords of a corner point in a 3d hex
   icorner pts 1 to 8 are ordered by x, then by y, then by z
------------------------------------------------------------------------- */

void Surf::hex_corner_point(int icorner, double *lo, double *hi, double *pt)
{
  if (icorner % 2) pt[0] = hi[0];
  else pt[0] = lo[0];
  if ((icorner/2) % 2) pt[1] = hi[1];
  else pt[1] = lo[1];
  if (icorner / 4) pt[2] = hi[2];
  else pt[2] = lo[2];
}

/* ----------------------------------------------------------------------
   return length of line M
------------------------------------------------------------------------- */

double Surf::line_size(int m)
{
  double delta[3];
  MathExtra::sub3(pts[lines[m].p2].x,pts[lines[m].p1].x,delta);
  return MathExtra::len3(delta);
}

/* ----------------------------------------------------------------------
   compute side length and area of a triangle
   return len = length of shortest edge of triangle M
   return area = area of triangle M
------------------------------------------------------------------------- */

double Surf::tri_size(int m, double &len)
{
  double delta12[3],delta13[3],delta23[3],cross[3];

  MathExtra::sub3(pts[tris[m].p2].x,pts[tris[m].p1].x,delta12);
  MathExtra::sub3(pts[tris[m].p3].x,pts[tris[m].p1].x,delta13);
  MathExtra::sub3(pts[tris[m].p3].x,pts[tris[m].p2].x,delta23);
  len = MIN(MathExtra::len3(delta12),MathExtra::len3(delta13));
  len = MIN(len,MathExtra::len3(delta23));

  MathExtra::cross3(delta12,delta13,cross);
  double area = 0.5 * MathExtra::len3(cross);
  return area;
}

/* ----------------------------------------------------------------------
   add a surface collision model
------------------------------------------------------------------------- */

void Surf::add_collide(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal surf_collide command");

  // error check

  for (int i = 0; i < nsc; i++)
    if (strcmp(arg[0],sc[i]->id) == 0)
      error->all(FLERR,"Reuse of surf_collide ID");

  // extend SurfCollide list if necessary

  if (nsc == maxsc) {
    maxsc += DELTA;
    sc = (SurfCollide **)
      memory->srealloc(sc,maxsc*sizeof(SurfCollide *),"surf:sc");
  }

  // check if ID already exists

  if (0) return;

#define SURF_COLLIDE_CLASS
#define SurfCollideStyle(key,Class) \
  else if (strcmp(arg[1],#key) == 0) \
    sc[nsc] = new Class(sparta,narg,arg);
#include "style_surf_collide.h"
#undef SurfCollideStyle
#undef SURF_COLLIDE_CLASS

  else error->all(FLERR,"Invalid surf_collide style");

  nsc++;
}

/* ----------------------------------------------------------------------
   find a surface collide model by ID
   return index of surf collide model or -1 if not found
------------------------------------------------------------------------- */

int Surf::find_collide(const char *id)
{
  int isc;
  for (isc = 0; isc < nsc; isc++)
    if (strcmp(id,sc[isc]->id) == 0) break;
  if (isc == nsc) return -1;
  return isc;
}

/* ----------------------------------------------------------------------
   brute force MPI Allreduce comm of local tallies across all procs
   input values
   for vector and array
   return out = summed tallies for surfs I own
------------------------------------------------------------------------- */

void Surf::collate_vec(int nrow, int *l2g, double *in, int istride,
                       double *out, int ostride, int sumflag)
{
  int i,j,m,n;
  double *vec1,*vec2;

  int nglobal;
  if (domain->dimension == 2) nglobal = nline;
  else nglobal = ntri;
  if (nglobal == 0) return;

  double *one,*all;
  memory->create(one,nglobal,"surf:one");
  memory->create(all,nglobal,"surf:all");

  for (i = 0; i < nglobal; i++) one[i] = 0.0;

  m = 0;
  for (i = 0; i < nrow; i++) {
    one[l2g[i]] = in[m];
    m += istride;
  }

  MPI_Allreduce(one,all,nglobal,MPI_DOUBLE,MPI_SUM,world);

  if (sumflag) {
    m = 0;
    for (i = 0; i < nlocal; i++) {
      out[m] += all[mysurfs[i]];
      m += ostride;
    }
  } else {
    m = 0;
    for (i = 0; i < nlocal; i++) {
      out[m] = all[mysurfs[i]];
      m += ostride;
    }
  }

  memory->destroy(one);
  memory->destroy(all);
}

void Surf::collate_array(int nrow, int ncol, int *l2g,
                         double **in, double **out)
{
  int i,j,m,n;
  double *vec1,*vec2;

  int nglobal;
  if (domain->dimension == 2) nglobal = nline;
  else nglobal = ntri;
  if (nglobal == 0) return;

  double **one,**all;
  memory->create(one,nglobal,ncol,"surf:one");
  memory->create(all,nglobal,ncol,"surf:all");

  for (i = 0; i < nglobal; i++)
    for (j = 0; j < ncol; j++)
      one[i][j] = 0.0;

  for (i = 0; i < nrow; i++) {
    m = l2g[i];
    for (j = 0; j < ncol; j++) 
      one[m][j] = in[i][j];
  }

  MPI_Allreduce(&one[0][0],&all[0][0],nglobal*ncol,MPI_DOUBLE,MPI_SUM,world);

  for (i = 0; i < nlocal; i++) {
    m = mysurfs[i];
    for (j = 0; j < ncol; j++) 
      out[i][j] += all[m][j];
  }
  
  memory->destroy(one);
  memory->destroy(all);
}

/* ----------------------------------------------------------------------
   makr 4 corner point flags for a grid cell containing lines
     by shooting a line from each corner pt to midpt of one line in cell
   n = # of lines, list = indices of lines
   lo/hi = grid cell corner pts
   attempt to return corner flags for each corner pts (ordered by x, y)
   flag = CORNER OUTSIDE/INSIDE/OVERLAP or no change if unable to mark
   unable to mark any corner pts if no line midpt is inside cell
   unable to mark specific corner pt if:
     no collisions with any line
       can happen if start and stop points are on every infinite surf
     1st collision is within epsilon of line end pt
       avoids round-off issues with hitting 2 lines at their common end pt
------------------------------------------------------------------------- */

void Surf::all_cell_corner_line(int n, int *list, double *lo, double *hi, 
				int *corner)
{
  bool hitflag;
  int i,m,isurf,flag,cflag,side,minside;
  double param,minparam;
  double start[3],stop[3],xc[3],delta[3];
  double *x1,*x2;
  Line *line;

  // stop = any line mid point that is inside or on the cell

  for (isurf = 0; isurf < n; isurf++) {
    line = &lines[list[isurf]];
    x1 = pts[line->p1].x;
    x2 = pts[line->p2].x;
    stop[0] = 0.5 * (x1[0]+x2[0]);
    stop[1] = 0.5 * (x1[1]+x2[1]);
    stop[2] = 0.0;
    if (Geometry::point_in_hex(stop,lo,hi)) break;
  }

  // no mid points in cell, just return

  if (isurf >= n) return;

  // loop over corner points of quad

  for (int icorner = 0; icorner < 4; icorner++) {
    quad_corner_point(icorner,lo,hi,start);

    // start = stop is special case
    // don't call line_line_intersect() with zero-length line
    // simply mark corner pt as being on surface

    if (start[0] == stop[0] && start[1] == stop[1]) {
      corner[icorner] = CORNEROVERLAP;
      continue;
    }

    // find first parametric collision of start-to-stop with any line
    // enforce hit with isurf since stop pt is on that surface
    //   but no hit with isurf if start pt is on infinite surf
    // hits with other surfs detected by line_line_intersect()
    //   will not return collision if start pt is on infinite surf
    //   do not count hit if collision pt is within epsilon of line end pt
    //   this avoids round-off issue with hitting 2 surfs at end pt,
    //     one on the inside, one on the outside
    //   only test this for side = OUTSIDE/INSIDE,
    //     since is ok if start is ON another surf close to its end pt
    
    cflag = 0;
    minparam = 2.0;
    for (m = 0; m < n; m++) {
      line = &lines[list[m]];
      if (m == isurf) {
	xc[0] = stop[0];
	xc[1] = stop[1];
	xc[2] = 0.0;
	param = 1.0;
	flag = Geometry::whichside(pts[line->p1].x,line->norm,
				   start[0],start[1],start[2]);
	hitflag = 1;
	if (flag < 0) side = INSIDE;
	else if (flag > 0) side = OUTSIDE;
	else hitflag = 0;
      } else {
	line = &lines[list[m]];
	hitflag = Geometry::
	  line_line_intersect(start,stop,
			      pts[line->p1].x,pts[line->p2].x,line->norm,
			      xc,param,side);
      }
      if (hitflag && param <= minparam) {
	if (side == OUTSIDE || side == INSIDE) {
	  if (Geometry::line_fraction(xc,pts[line->p1].x,
				      pts[line->p2].x) < EPSSQ) continue;
	}
	minparam = param;
	minside = side;
	cflag = 1;
      }
    }

    // cflag = 0 means no intersection found, so do not mark corner pt
    
    if (!cflag) continue;
    else if (minside == OUTSIDE) corner[icorner] = CORNEROUTSIDE;
    else if (minside == INSIDE) corner[icorner] = CORNERINSIDE;
    else corner[icorner] = CORNEROVERLAP;
  }
}

/* ----------------------------------------------------------------------
   compute 8 corner point flags for a grid cell containing triangles
     by shooting a line from each corner pt to centroid of one tri in cell
   n = # of tris, list = indices of tris
   lo/hi = grid cell corner pts
   attempt to return corner flags for each corner pts (ordered by x, y, z)
   flag = CORNER OUTSIDE/INSIDE/OVERLAP or no change if unable to mark
   unable to mark any corner pts if no tri centroid is inside cell
   unable to mark specific corner pt if:
     no collisions with any tri
       can happen if start and stop points are on every infinite surf
     1st collision is within epsilon of tri corner pt
       avoids round-off issues with hitting 2+ tris at their common corner pt
------------------------------------------------------------------------- */

void Surf::all_cell_corner_tri(int n, int *list, double *lo, double *hi, 
			       int *corner)
{
  bool hitflag;
  int i,m,isurf,flag,cflag,side,minside;
  double param,minparam;
  double start[3],stop[3],xc[3],delta[3];
  double *x1,*x2,*x3;
  Tri *tri;

  // stop = any tri centroid that is inside or on the cell

  for (isurf = 0; isurf < n; isurf++) {
    tri = &tris[list[isurf]];
    x1 = pts[tris->p1].x;
    x2 = pts[tris->p2].x;
    x3 = pts[tris->p3].x;
    stop[0] = (x1[0]+x2[0]+x3[0])/3.0;
    stop[1] = (x1[1]+x2[1]+x3[1])/3.0;
    stop[2] = (x1[2]+x2[2]+x3[2])/3.0;
    if (Geometry::point_in_hex(stop,lo,hi)) break;
  }

  // no centroids in cell, just return

  if (isurf >= n) return;

  // loop over corner points of hex

  for (int icorner = 0; icorner < 8; icorner++) {
    hex_corner_point(icorner,lo,hi,start);

    // start = stop is special case
    // don't call line_tri_intersect() with zero-length line
    // simply mark corner pt as being on surface

    if (start[0] == stop[0] && start[1] == stop[1] && start[2] == stop[2]) {
      corner[icorner] = CORNEROVERLAP;
      continue;
    }

    // find first parametric collision of start-to-stop with any tri
    // enforce hit with isurf since stop pt is on that surface
    //   but no hit with isurf if start pt is on infinite surf
    // hits with other surfs detected by line_tri_intersect()
    //   will not return collision if start pt is on infinite surf
    //   do not count hit if collision pt is within epsilon of tri corner pt
    //   this avoids round-off issue with hitting 2+ surfs at corner pt,
    //     one on the inside, one on the outside
    //   only test this for side = OUTSIDE/INSIDE,
    //     since is ok if start is ON another surf close to its corner pt

    cflag = 0;
    minparam = 2.0;
    for (m = 0; m < n; m++) {
      tri = &tris[list[m]];
      if (m == isurf) {
	xc[0] = stop[0];
	xc[1] = stop[1];
	xc[2] = stop[2];
	param = 1.0;
	flag = Geometry::whichside(pts[tri->p1].x,tri->norm,
				   start[0],start[1],start[2]);
	hitflag = 1;
	if (flag < 0) side = INSIDE;
	else if (flag > 0) side = OUTSIDE;
	else hitflag = 0;
      } else {
	tri = &tris[list[m]];
	hitflag = Geometry::
	  line_tri_intersect(start,stop,
			     pts[tri->p1].x,pts[tri->p2].x,pts[tri->p3].x,
			     tri->norm,xc,param,side);
      }
      if (hitflag && param <= minparam) {
	if (side == OUTSIDE || side == INSIDE) {
	  if (Geometry::tri_fraction(xc,pts[tri->p1].x,pts[tri->p2].x,
				     pts[tri->p2].x) < EPSSQ) continue;
	}
	minparam = param;
	minside = side;
	cflag = 1;
      }
    }

    // cflag = 0 means no intersection found, so do not mark corner pt
    
    if (!cflag) continue;
    else if (minside == OUTSIDE) corner[icorner] = CORNEROUTSIDE;
    else if (minside == INSIDE) corner[icorner] = CORNERINSIDE;
    else corner[icorner] = CORNEROVERLAP;
  }
}

/* ----------------------------------------------------------------------
   compute 1 corner point flag for a grid cell containing lines
     by shooting lines from corner pt to other marked corners pts of cell
   ic = corner point index (0 to 3)
   n = # of lines, list = indices of lines
   lo/hi = grid cell corner pts
   corner = all corner flags for cell
   return CORNER OUTSIDE/INSIDE/OVERLAP for corner pt if can mark it
   return CORNERUNKNOWN if cannot mark it
------------------------------------------------------------------------- */

int Surf::one_cell_corner_line(int ic, int n, int *list, 
			       double *lo, double *hi, int *corner)
{
  int i,m;
  int hitflag,cflag,side,minside;
  double param,minparam;
  double start[3],stop[3],xc[3],delta[3];
  Line *line;

  // start = unset corner pt

  quad_corner_point(ic,lo,hi,start);
  
  // loop over all corner pts that are OUTSIDE/INSIDE
  // shoot line from unmarked corner pt to marked corner pt
  // if no hit, return mark of marked corner pt
  // if find first collision, return marking it induces

  for (i = 0; i < 4; i++) {
    if (corner[i] != CORNEROUTSIDE && corner[i] != CORNERINSIDE) continue;
    quad_corner_point(i,lo,hi,stop);

    // find first parametric collision of start-to-stop with any line
    // set cflag = 2 if collision pt is within epsilon of line end pt
    // this is b/c could hit different sides of 2 lines within round-off

    cflag = 0;
    minparam = 2.0;
    for (m = 0; m < n; m++) {
      line = &lines[list[m]];
      hitflag = Geometry::
	line_line_intersect(start,stop,
			    pts[line->p1].x,pts[line->p2].x,line->norm,
			    xc,param,side);
      if (hitflag && param <= minparam) {
	minparam = param;
	minside = side;
	cflag = 1;
	if (Geometry::line_fraction(xc,pts[line->p1].x,
				    pts[line->p2].x) < EPSSQ) cflag = 2;
      }
    }

    // if no collisions, return stop point mark
    // if first collision near line end pt, don't count it
    // else return mark based on first collision

    if (!cflag) return corner[i];
    else if (cflag == 2) continue;
    else if (minside == OUTSIDE) return CORNEROUTSIDE;
    else if (minside == INSIDE) return CORNERINSIDE;
    else return CORNEROVERLAP;
  }

  // no collisions at all, return unmarked

  return CORNERUNKNOWN;
}

/* ----------------------------------------------------------------------
   compute 1 corner point flag for a grid cell containing triangles
     by shooting lines from corner pt to other marked corners pts of cell
   ic = corner point index (0 to 7)
   n = # of tris, list = indices of tris
   lo/hi = grid cell corner pts
   corner = all corner flags for cell
   return CORNER OUTSIDE/INSIDE/OVERLAP for corner pt if can mark it
   return CORNERUNKNOWN if cannot mark it
------------------------------------------------------------------------- */

int Surf::one_cell_corner_tri(int ic, int n, int *list, 
			      double *lo, double *hi, int *corner)
{
  int i,m;
  int hitflag,cflag,side,minside;
  double param,minparam;
  double start[3],stop[3],xc[3],delta[3];
  Tri *tri;

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  // start = unset corner pt

  hex_corner_point(ic,lo,hi,start);
  
  // loop over all corner pts that are OUTSIDE/INSIDE
  // shoot line from unmarked corner pt to marked corner pt
  // if no hit, return mark of marked corner pt
  // if find first collision, return marking it induces

  for (i = 0; i < 8; i++) {
    if (corner[i] != CORNEROUTSIDE && corner[i] != CORNERINSIDE) continue;
    hex_corner_point(i,lo,hi,stop);

    // find first parametric collision of start-to-stop with any tri
    // set cflag = 2 if collision pt is within epsilon of tri corner pt
    // this is b/c could hit different sides of 2 triangles within round-off

    cflag = 0;
    minparam = 2.0;
    for (m = 0; m < n; m++) {
      tri = &tris[list[m]];
      hitflag = Geometry::
	line_tri_intersect(start,stop,
			   pts[tri->p1].x,pts[tri->p2].x,pts[tri->p3].x,
			   tri->norm,xc,param,side);
      if (hitflag && param <= minparam) {
	minparam = param;
	minside = side;
	cflag = 1;
	if (Geometry::tri_fraction(xc,pts[tri->p1].x,pts[tri->p2].x,
				   pts[tri->p2].x) < EPSSQ) cflag = 2;
      }
    }
    
    // if no collisions, return stop point mark
    // if first collision near tri corner pt, don't count it
    // else return mark based on first collision

    if (!cflag) return corner[i];
    else if (cflag == 2) continue;
    else if (minside == OUTSIDE) return CORNEROUTSIDE;
    else if (minside == INSIDE) return CORNERINSIDE;
    else return CORNEROVERLAP;
  }

  // no collisions at all, return unmarked

  return CORNERUNKNOWN;
}

/* ---------------------------------------------------------------------- */

bigint Surf::memory_usage()
{
  bigint bytes = 0;
  bytes += (bigint) npoint * sizeof(Point);
  bytes += (bigint) nline * sizeof(Line);
  bytes += (bigint) ntri * sizeof(Tri);
  return bytes;
}

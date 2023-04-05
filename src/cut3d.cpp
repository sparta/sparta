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

#include "math.h"
#include "string.h"
#include "cut3d.h"
#include "cut2d.h"
#include "surf.h"
#include "domain.h"
#include "grid.h"
#include "comm.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};     // several files
enum{CTRI,CTRIFACE,FACEPGON,FACE};
enum{EXTERIOR,INTERIOR,BORDER};
enum{ENTRY,EXIT,TWO,CORNER};              // same as Cut2d

#define EPSEDGE 1.0e-9    // minimum edge length (fraction of cell size)
#define SHRINK 1.0e-8     // shrink grid cell by this fraction when split fails

// cell ID for 2d or 3d cell

//#define VERBOSE
//#define VERBOSE_ID 61275

/* ---------------------------------------------------------------------- */

Cut3d::Cut3d(SPARTA *sparta) : Pointers(sparta)
{
  implicit = surf->implicit;

  ntiny = 0;
  nshrink = 0;

  cut2d = new Cut2d(sparta,0);
  memory->create(path1,12,3,"cut3d:path1");
  memory->create(path2,12,3,"cut3d:path2");
}

/* ---------------------------------------------------------------------- */

Cut3d::~Cut3d()
{
  delete cut2d;
  memory->destroy(path1);
  memory->destroy(path2);
}

/* ----------------------------------------------------------------------
   compute intersections of a grid cell with all surfs
   csurfs = indices into global surf list
   return nsurf = # of surfs
   return -1 if nsurf > max
------------------------------------------------------------------------- */

int Cut3d::surf2grid(cellint id_caller, double *lo_caller, double *hi_caller,
                     surfint *surfs_caller, int max)
{
  id = id_caller;
  lo[0] = lo_caller[0]; lo[1] = lo_caller[1]; lo[2] = lo_caller[2];
  hi[0] = hi_caller[0]; hi[1] = hi_caller[1]; hi[2] = hi_caller[2];
  surfs = surfs_caller;

  Surf::Tri *tris = surf->tris;
  int ntotal = surf->nsurf;

  double value;
  double *x1,*x2,*x3;

  nsurf = 0;
  for (int m = 0; m < ntotal; m++) {
    x1 = tris[m].p1;
    x2 = tris[m].p2;
    x3 = tris[m].p3;

    value = MAX(x1[0],x2[0]);
    if (MAX(value,x3[0]) < lo[0]) continue;
    value = MIN(x1[0],x2[0]);
    if (MIN(value,x3[0]) > hi[0]) continue;

    value = MAX(x1[1],x2[1]);
    if (MAX(value,x3[1]) < lo[1]) continue;
    value = MIN(x1[1],x2[1]);
    if (MIN(value,x3[1]) > hi[1]) continue;

    value = MAX(x1[2],x2[2]);
    if (MAX(value,x3[2]) < lo[2]) continue;
    value = MIN(x1[2],x2[2]);
    if (MIN(value,x3[2]) > hi[2]) continue;

    // 3 versions of this:
    // 1 = tri_hex_intersect with geometric line_tri_intersect,
    //     devel/cut3d.old1.py
    // 2 = tri_hex_intersect with tetvol line_tri_intersect, here
    // 3 = Sutherland-Hodgman clip algorithm, here and in devel/cut3d.py

    //if (tri_hex_intersect(x1,x2,x3,tris[m].norm)) {
    //  if (nsurf == max) return -1;
    //  surfs[nsurf++] = m;
    // }

    if (clip(x1,x2,x3)) {
      if (nsurf < max) surfs[nsurf] = m;
      nsurf++;
    }
  }

  return nsurf;
}

/* ----------------------------------------------------------------------
   compute intersections of a grid cell with a provided list of surfs
   csurfs = indices into global surf list
   nlist, list = vector of surf indices of length nlist
   return nsurf = # of surfs
   return -1 if nsurf > max
   called by AdaptGrid via Grid::surf2grid_one
------------------------------------------------------------------------- */

int Cut3d::surf2grid_list(cellint id_caller,
                          double *lo_caller, double *hi_caller,
                          int nlist, surfint *list,
                          surfint *surfs_caller, int max)
{
  id = id_caller;
  lo[0] = lo_caller[0]; lo[1] = lo_caller[1]; lo[2] = lo_caller[2];
  hi[0] = hi_caller[0]; hi[1] = hi_caller[1]; hi[2] = hi_caller[2];
  surfs = surfs_caller;

  Surf::Tri *tris = surf->tris;

  int m;
  double value;
  double *x1,*x2,*x3;

  nsurf = 0;
  for (int i = 0; i < nlist; i++) {
    m = list[i];
    x1 = tris[m].p1;
    x2 = tris[m].p2;
    x3 = tris[m].p3;

    value = MAX(x1[0],x2[0]);
    if (MAX(value,x3[0]) < lo[0]) continue;
    value = MIN(x1[0],x2[0]);
    if (MIN(value,x3[0]) > hi[0]) continue;

    value = MAX(x1[1],x2[1]);
    if (MAX(value,x3[1]) < lo[1]) continue;
    value = MIN(x1[1],x2[1]);
    if (MIN(value,x3[1]) > hi[1]) continue;

    value = MAX(x1[2],x2[2]);
    if (MAX(value,x3[2]) < lo[2]) continue;
    value = MIN(x1[2],x2[2]);
    if (MIN(value,x3[2]) > hi[2]) continue;

    // 3 versions of this:
    // 1 = tri_hex_intersect with geometric line_tri_intersect,
    //     devel/cut3d.old1.py
    // 2 = tri_hex_intersect with tetvol line_tri_intersect, here
    // 3 = Sutherland-Hodgman clip algorithm, here and in devel/cut3d.py

    //if (tri_hex_intersect(x1,x2,x3,tris[m].norm)) {
    //  if (nsurf == max) return -1;
    //  surfs[nsurf++] = m;
    // }

    if (clip(x1,x2,x3)) {
      if (nsurf < max) surfs[nsurf] = m;
      nsurf++;
    }
  }

  return nsurf;
}

/* ----------------------------------------------------------------------
   compute intersections of a grid cell with a single surf
   p012 = corner pts of surf
   lo,hi = grid cell corner points
   return 1 if intersects, 0 if not
   called by Grid::surf2grid2
------------------------------------------------------------------------- */

int Cut3d::surf2grid_one(double *p0, double *p1, double *p2,
                         double *lo_caller, double *hi_caller)
{
  lo[0] = lo_caller[0]; lo[1] = lo_caller[1]; lo[2] = lo_caller[2];
  hi[0] = hi_caller[0]; hi[1] = hi_caller[1]; hi[2] = hi_caller[2];
  return clip(p0,p1,p2);
}

/* ----------------------------------------------------------------------
   calculate cut volume of a grid cell that contains nsurf tris
   also calculate if cell is split into distinct sub-volumes by tris
   return nsplit = # of splits, 1 for no split
   return vols = ptr to vector of vols = one vol per split
     if nsplit = 1, cut vol
     if nsplit > 1, one vol per split cell
   return corners = UNKNOWN/INSIDE/OUTSIDE for each of 8 corner pts
   if nsplit > 1, also return:
     surfmap = which sub-cell (0 to Nsurfs-1) each surf is in
             = -1 if not in any sub-cell, discarded by add_tris
     xsplit = coords of a point in one of the split cells
     xsub = which sub-cell (0 to Nsplit-1) xsplit is in
   work is done by split_try()
   call it once more with shrunk grid cell if first attempt fails
------------------------------------------------------------------------- */

int Cut3d::split(cellint id_caller, double *lo_caller, double *hi_caller,
                 int nsurf_caller, surfint *surfs_caller,
                 double *&vols_caller, int *surfmap,
                 int *corners, int &xsub, double *xsplit)
{
  lo[0] = lo_caller[0]; lo[1] = lo_caller[1]; lo[2] = lo_caller[2];
  hi[0] = hi_caller[0]; hi[1] = hi_caller[1]; hi[2] = hi_caller[2];

  // perform split with full-size grid cell

  int nsplit;
  int errflag =
    split_try(id_caller,nsurf_caller,surfs_caller,vols_caller,surfmap,
              corners,xsub,xsplit,nsplit);

  // error return
  // try again after shrinking grid cell by SHRINK factor
  // this gets rid of pesky errors due to tri pts/edges on cell faces

  if (errflag) {
    nshrink++;

    double newlo = lo[0] + SHRINK*(hi[0]-lo[0]);
    double newhi = hi[0] - SHRINK*(hi[0]-lo[0]);
    lo[0] = newlo;
    hi[0] = newhi;

    newlo = lo[1] + SHRINK*(hi[1]-lo[1]);
    newhi = hi[1] - SHRINK*(hi[1]-lo[1]);
    lo[1] = newlo;
    hi[1] = newhi;

    newlo = lo[2] + SHRINK*(hi[2]-lo[2]);
    newhi = hi[2] - SHRINK*(hi[2]-lo[2]);
    lo[2] = newlo;
    hi[2] = newhi;

    errflag =
      split_try(id_caller,nsurf_caller,surfs_caller,vols_caller,surfmap,
                corners,xsub,xsplit,nsplit);
  }

  // could not perform cut/split -> fatal error
  // print info about cell and final error message
  // NOTE: could unshrink cell first

  if (errflag) {
    failed_cell();
    split_error(errflag);
  }

  return nsplit;
}

/* ----------------------------------------------------------------------
   clip triangle P0 P1 P2 against cell with corners CLO,CHI
   tri may or may not intersect cell (due to rounding)
   return # of clipped points, can be 0,1,2,3 up to 8 (I think)
   return clipped points in cpath as series of x,y,z triplets
   called externally, depends on no class variables
   duplicate points in cpath are deleted
   uses Sutherland-Hodgman clipping algorithm
------------------------------------------------------------------------- */

int Cut3d::clip_external(double *p0, double *p1, double *p2,
                         double *clo, double *chi, double *cpath)
{
  int i,npath,nnew;
  double value;
  double *s,*e;
  double **path,**newpath;

  // initial path = tri vertices

  nnew = 3;
  memcpy(path1[0],p0,3*sizeof(double));
  memcpy(path1[1],p1,3*sizeof(double));
  memcpy(path1[2],p2,3*sizeof(double));

  // clip tri against each of 6 grid face planes

  for (int dim = 0; dim < 3; dim++) {
    path = path1;
    newpath = path2;
    npath = nnew;
    nnew = 0;

    value = clo[dim];
    s = path[npath-1];
    for (i = 0; i < npath; i++) {
      e = path[i];
      if (e[dim] >= value) {
        if (s[dim] < value) between(s,e,dim,value,newpath[nnew++]);
        memcpy(newpath[nnew++],e,3*sizeof(double));
      } else if (s[dim] >= value) between(e,s,dim,value,newpath[nnew++]);
      s = e;
    }
    if (!nnew) return 0;

    path = path2;
    newpath = path1;
    npath = nnew;
    nnew = 0;

    value = chi[dim];
    s = path[npath-1];
    for (i = 0; i < npath; i++) {
      e = path[i];
      if (e[dim] <= value) {
        if (s[dim] > value) between(s,e,dim,value,newpath[nnew++]);
        memcpy(newpath[nnew++],e,3*sizeof(double));
      } else if (s[dim] <= value) between(e,s,dim,value,newpath[nnew++]);
      s = e;
    }
    if (!nnew) return 0;
  }

  // copy path points to cpath
  // delete any duplicate points while copying
  // inner-loop check is to not add a point that duplicates previous point
  // post-loop check is for duplicate due to 1st point = last point

  int m = 0;
  int n = nnew;
  nnew = 0;
  for (i = 0; i < n; i++) {
    if (m) {
      if (path1[i][0] == cpath[m-3] && path1[i][1] == cpath[m-2] &&
          path1[i][2] == cpath[m-1]) continue;
    }
    cpath[m++] = path1[i][0];
    cpath[m++] = path1[i][1];
    cpath[m++] = path1[i][2];
    nnew++;
  }
  if (nnew > 1)
    if (cpath[0] == cpath[m-3] && cpath[1] == cpath[m-2] &&
        cpath[2] == cpath[m-1]) nnew--;

  return nnew;
}

/* ----------------------------------------------------------------------
   check if pts A,B,C are on same face of cell
   return 1,2,3,4,5,6 = both on left,right,lower,upper,bottom,top face
   return 0 if not, including inside
------------------------------------------------------------------------- */

int Cut3d::sameface(double *a, double *b, double *c)
{
  if (a[0] == lo[0] and b[0] == lo[0] and c[0] == lo[0]) return 1;
  if (a[0] == hi[0] and b[0] == hi[0] and c[0] == hi[0]) return 2;
  if (a[1] == lo[1] and b[1] == lo[1] and c[1] == lo[1]) return 3;
  if (a[1] == hi[1] and b[1] == hi[1] and c[1] == hi[1]) return 4;
  if (a[2] == lo[2] and b[2] == lo[2] and c[2] == lo[2]) return 5;
  if (a[2] == hi[2] and b[2] == hi[2] and c[2] == hi[2]) return 6;
  return 0;
}

/* ----------------------------------------------------------------------
   check if pts A,B,C are on same face of cell
   return 1,2,3,4,5,6 = both on left,right,lower,upper,bottom,top face
   return 0 if not, including inside
   called externally, depends on no class variables
------------------------------------------------------------------------- */

int Cut3d::sameface_external(double *a, double *b, double *c, double *clo, double *chi)
{
  if (a[0] == clo[0] and b[0] == clo[0] and c[0] == clo[0]) return 1;
  if (a[0] == chi[0] and b[0] == chi[0] and c[0] == chi[0]) return 2;
  if (a[1] == clo[1] and b[1] == clo[1] and c[1] == clo[1]) return 3;
  if (a[1] == chi[1] and b[1] == chi[1] and c[1] == chi[1]) return 4;
  if (a[2] == clo[2] and b[2] == clo[2] and c[2] == clo[2]) return 5;
  if (a[2] == chi[2] and b[2] == chi[2] and c[2] == chi[2]) return 6;
  return 0;
}

// ----------------------------------------------------------------------
// internal methods
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   Sutherland-Hodgman clipping algorithm
   don't need to delete duplicate points since touching counts as intersection
------------------------------------------------------------------------- */

int Cut3d::clip(double *p0, double *p1, double *p2)
{
  int i,npath,nnew;
  double value;
  double *s,*e;
  double **path,**newpath;

  // initial path = tri vertices

  nnew = 3;
  memcpy(path1[0],p0,3*sizeof(double));
  memcpy(path1[1],p1,3*sizeof(double));
  memcpy(path1[2],p2,3*sizeof(double));

  // intersect if any of tri vertices is within grid cell

  if (p0[0] >= lo[0] && p0[0] <= hi[0] &&
      p0[1] >= lo[1] && p0[1] <= hi[1] &&
      p0[2] >= lo[2] && p0[2] <= hi[2] &&
      p1[0] >= lo[0] && p1[0] <= hi[0] &&
      p1[1] >= lo[1] && p1[1] <= hi[1] &&
      p1[2] >= lo[2] && p1[2] <= hi[2] &&
      p2[0] >= lo[0] && p2[0] <= hi[0] &&
      p2[1] >= lo[1] && p2[1] <= hi[1] &&
      p2[2] >= lo[2] && p2[2] <= hi[2]) return 1;

  // clip tri against each of 6 grid face planes

  for (int dim = 0; dim < 3; dim++) {
    path = path1;
    newpath = path2;
    npath = nnew;
    nnew = 0;

    value = lo[dim];
    s = path[npath-1];
    for (i = 0; i < npath; i++) {
      e = path[i];
      if (e[dim] >= value) {
        if (s[dim] < value) between(s,e,dim,value,newpath[nnew++]);
        memcpy(newpath[nnew++],e,3*sizeof(double));
      } else if (s[dim] >= value) between(e,s,dim,value,newpath[nnew++]);
      s = e;
    }
    if (!nnew) return 0;

    path = path2;
    newpath = path1;
    npath = nnew;
    nnew = 0;

    value = hi[dim];
    s = path[npath-1];
    for (i = 0; i < npath; i++) {
      e = path[i];
      if (e[dim] <= value) {
        if (s[dim] > value) between(s,e,dim,value,newpath[nnew++]);
        memcpy(newpath[nnew++],e,3*sizeof(double));
      } else if (s[dim] <= value) between(e,s,dim,value,newpath[nnew++]);
      s = e;
    }
    if (!nnew) return 0;
  }

  return nnew;
}

/* ----------------------------------------------------------------------
   attempt to split cell
   return 0 if successful, otherwise return an error flag
   called multiple times by split() with slightly different grid cell sizes
------------------------------------------------------------------------- */

int Cut3d::split_try(cellint id_caller,
                     int nsurf_caller, surfint *surfs_caller,
                     double *&vols_caller, int *surfmap,
                     int *corners, int &xsub, double *xsplit, int &nsplit)
{
  id = id_caller;
  nsurf = nsurf_caller;
  surfs = surfs_caller;

  int errflag = add_tris();
  if (errflag) return errflag;

#ifdef VERBOSE
  if (id == VERBOSE_ID) print_bpg("BPG after added tris");
#endif

  clip_tris();
  clip_adjust();

#ifdef VERBOSE
  if (id == VERBOSE_ID) print_bpg("BPG after clipped tris");
#endif

  // all triangles just touched cell surface
  // mark corner points based on grazecount or touchmark value
  // return vol = 0.0 for UNKNOWN/INSIDE, full cell vol for OUTSIDE
  // vol is changed in Grid::set_inout() if OVERLAP cell corners are marked

  if (empty) {
    int mark = UNKNOWN;
    if (grazecount || touchmark == INSIDE) mark = INSIDE;
    else if (touchmark == OUTSIDE) mark = OUTSIDE;
    corners[0] = corners[1] = corners[2] = corners[3] =
      corners[4] = corners[5] = corners[6] = corners[7] = mark;

    double vol = 0.0;
    if (mark == OUTSIDE) vol = (hi[0]-lo[0]) * (hi[1]-lo[1]) * (hi[2]-lo[2]);

    vols.grow(1);
    vols[0] = vol;
    vols_caller = &vols[0];
    nsplit = 1;
    return 0;
  }

  ctri_volume();
  errflag = edge2face();
  if (errflag) return errflag;

  double lo2d[2],hi2d[2];

  for (int iface = 0; iface < 6; iface++) {
    if (facelist[iface].n) {
      face_from_cell(iface,lo2d,hi2d);
      edge2clines(iface);
      errflag = cut2d->split_face(id,iface,lo2d,hi2d);
      if (errflag) return errflag;
      errflag = add_face_pgons(iface);
      if (errflag) return errflag;
    } else {
      face_from_cell(iface,lo2d,hi2d);
      errflag = add_face(iface,lo2d,hi2d);
      if (errflag) return errflag;
    }
  }

  remove_faces();

#ifdef VERBOSE
  if (id == VERBOSE_ID) print_bpg("BPG after faces");
#endif

  errflag = check();
  if (errflag) return errflag;

  walk();

#ifdef VERBOSE
  if (id == VERBOSE_ID) print_loops();
#endif

  errflag = loop2ph();

  // loop2ph detected no positive-volume loops, cell is inside the surf

  if (errflag == 4) {
    corners[0] = corners[1] = corners[2] = corners[3] =
      corners[4] = corners[5] = corners[6] = corners[7] = INSIDE;
    double volume = 0.0;
    vols.grow(1);
    vols[0] = volume;
    vols_caller = &vols[0];
    nsplit = 1;
    return 0;
  }

  // other error returns from loop2ph

  if (errflag) return errflag;

  // if multiple splits, find a split point

  nsplit = phs.n;
  if (nsplit > 1) {
    create_surfmap(surfmap);
    if (implicit) errflag = split_point_implicit(surfmap,xsplit,xsub);
    else errflag = split_point_explicit(surfmap,xsplit,xsub);
  }
  if (errflag) return errflag;

  // successful cut/split
  // set corners = OUTSIDE if corner pt is in list of edge points
  // else set corners = INSIDE

  int icorner;
  double *p1,*p2;

  corners[0] = corners[1] = corners[2] = corners[3] =
    corners[4] = corners[5] = corners[6] = corners[7] = INSIDE;

  int nedge = edges.n;
  for (int iedge = 0; iedge < nedge; iedge++) {
    if (!edges[iedge].active) continue;
    p1 = edges[iedge].p1;
    p2 = edges[iedge].p2;
    icorner = corner(p1);
    if (icorner >= 0) corners[icorner] = OUTSIDE;
    icorner = corner(p2);
    if (icorner >= 0) corners[icorner] = OUTSIDE;
  }

  // store volumes in vector so can return ptr to it

  vols.grow(nsplit);
  for (int i = 0; i < nsplit; i++) vols[i] = phs[i].volume;
  vols_caller = &vols[0];

  // successful exit

  return 0;
}

/* ----------------------------------------------------------------------
   2-letter prefix is which method encountered error
------------------------------------------------------------------------- */

void Cut3d::split_error(int errflag)
{
  if (errflag == 1)
    error->one(FLERR,"FE: Found edge in same direction");
  if (errflag == 2)
    error->one(FLERR,"EF: Singlet BPG edge not on cell face");
  if (errflag == 3)
    error->one(FLERR,"EF: BPG edge on more than 2 faces");
  if (errflag == 4)
    error->one(FLERR,"LP: No positive volumes in cell");
  if (errflag == 5)
    error->one(FLERR,"LP: More than one positive volume with "
               "a negative volume");
  if (errflag == 6)
    error->one(FLERR,"LP: Single volume is negative, inverse donut");
  if (errflag == 7)
    error->one(FLERR,"SP: Could not find split point in split cell");

  if (errflag == 11)
    error->one(FLERR,"CH: Vertex has less than 3 edges");
  if (errflag == 12)
    error->one(FLERR,"CH: Vertex contains invalid edge");
  if (errflag == 13)
    error->one(FLERR,"CH: Vertex contains edge that doesn't point to it");
  if (errflag == 14)
    error->one(FLERR,"CH: Vertex contains duplicate edge");
  if (errflag == 15)
    error->one(FLERR,"CH: Vertex pointers to last edge are invalid");
  if (errflag == 16)
    error->one(FLERR,"CH: Edge not part of 2 vertices");
  if (errflag == 17)
    error->one(FLERR,"CH: Edge part of same vertex twice");
  if (errflag == 18)
    error->one(FLERR,"CH: Edge part of invalid vertex");
  if (errflag == 19)
    error->one(FLERR,"CH: Edge part of invalid vertex");

  if (errflag == 21)
    error->one(FLERR,"WB: Point appears first in more than one CLINE");
  if (errflag == 22)
    error->one(FLERR,"WB: Point appears last in more than one CLINE");
  if (errflag == 23)
    error->one(FLERR,"WB: Singlet CLINES point not on cell border");
  if (errflag == 24)
    error->one(FLERR,"LP: No positive areas in cell");
  if (errflag == 25)
    error->one(FLERR,"LP: More than one positive area with a negative area");
  if (errflag == 26)
    error->one(FLERR,"LP: Single area is negative, inverse donut");
}

/* ----------------------------------------------------------------------
   add each triangle as vertex and edges to BPG
   add full edge even if outside cell, clipping comes later
   skip transparent surfs
------------------------------------------------------------------------- */

int Cut3d::add_tris()
{
  int i,m;
  int e1,e2,e3,dir1,dir2,dir3;
  double p1[3],p2[3],p3[3];
  Surf::Tri *tri;
  Vertex *vert;
  Edge *edge;

  Surf::Tri *tris = surf->tris;

  verts.grow(nsurf);
  edges.grow(3*nsurf);
  verts.n = 0;
  edges.n = 0;

  int nvert = 0;
  for (i = 0; i < nsurf; i++) {
    m = surfs[i];
    tri = &tris[m];
    if (tri->transparent) continue;

    memcpy(p1,tri->p1,3*sizeof(double));
    memcpy(p2,tri->p2,3*sizeof(double));
    memcpy(p3,tri->p3,3*sizeof(double));

    vert = &verts[nvert];
    vert->active = 1;
    vert->style = CTRI;
    vert->label = i;
    vert->nedge = 0;
    vert->volume = 0.0;
    vert->norm = tri->norm;

    // look for each edge of tri
    // add to edges in forward dir if doesn't yet exist
    // add to edges in reverse dir if already exists

    e1 = findedge(p1,p2,0,dir1);
    if (e1 == -2) return 1;

    if (e1 < 0) {
      e1 = edges.n++;
      dir1 = 0;
      edge = &edges[e1];
      edge->style = CTRI;
      edge->nvert = 0;
      edge->verts[0] = edge->verts[1] = -1;
      memcpy(edge->p1,p1,3*sizeof(double));
      memcpy(edge->p2,p2,3*sizeof(double));
    }
    edge_insert(e1,dir1,nvert,-1,-1,-1,-1);

    e2 = findedge(p2,p3,0,dir2);
    if (e2 == -2) return 1;

    if (e2 < 0) {
      e2 = edges.n++;
      dir2 = 0;
      edge = &edges[e2];
      edge->style = CTRI;
      edge->nvert = 0;
      edge->verts[0] = edge->verts[1] = -1;
      memcpy(edge->p1,p2,3*sizeof(double));
      memcpy(edge->p2,p3,3*sizeof(double));
    }
    edge_insert(e2,dir2,nvert,e1,dir1,-1,-1);

    e3 = findedge(p3,p1,0,dir3);
    if (e3 == -2) return 1;

    if (e3 < 0) {
      e3 = edges.n++;
      dir3 = 0;
      edge = &edges[e3];
      edge->style = CTRI;
      edge->nvert = 0;
      edge->verts[0] = edge->verts[1] = -1;
      memcpy(edge->p1,p3,3*sizeof(double));
      memcpy(edge->p2,p1,3*sizeof(double));
    }
    edge_insert(e3,dir3,nvert,e2,dir2,-1,-1);

    nvert++;
  }

  verts.n = nvert;
  return 0;
}

/* ----------------------------------------------------------------------
   clip collection of tris that overlap cell by 6 faces of cell
   loop over faces, loop over tris, loop over edges in tri
   edges fully outside the cell are removed
   shared edges that intersect the cell are clipped consistently
------------------------------------------------------------------------- */

void Cut3d::clip_tris()
{
  int i,n,dim,lohi,ivert,iedge,jedge,idir,jdir,nedge;
  int p1flag,p2flag;
  double value;
  double *p1,*p2,*p3;
  Edge *edge,*newedge;

  // loop over all 6 faces of cell

  int nvert = verts.n;

  for (int iface = 0; iface < 6; iface++) {
    dim = iface / 2;
    lohi = iface % 2;
    if (lohi == 0) value = lo[dim];
    else value = hi[dim];

    // mark all edges as unclipped
    // some may have been clipped and not cleared on previous face

    nedge = edges.n;
    for (iedge = 0; iedge < nedge; iedge++)
      if (edges[iedge].active) edges[iedge].clipped = 0;

    // loop over vertices, clip each of its edges to face
    // if edge already clipped, unset clip flag and keep edge as-is

    for (ivert = 0; ivert < nvert; ivert++) {
      iedge = verts[ivert].first;
      idir = verts[ivert].dirfirst;
      nedge = verts[ivert].nedge;

      for (i = 0; i < nedge; i++) {
        edge = &edges[iedge];

        if (edge->clipped) {
          edge->clipped = 0;
          iedge = edge->next[idir];
          idir = edge->dirnext[idir];
          continue;
        }

        // p1/p2 are pts in order of traversal

        if (idir == 0) {
          p1 = edge->p1;
          p2 = edge->p2;
        } else {
          p1 = edge->p2;
          p2 = edge->p1;
        }

        // p1/p2 flag = OUTSIDE/ON/INSIDE for edge pts

        if (lohi == 0) {
          if (p1[dim] < value) p1flag = OUTSIDE;
          else if (p1[dim] > value) p1flag = INSIDE;
          else p1flag = OVERLAP;
          if (p2[dim] < value) p2flag = OUTSIDE;
          else if (p2[dim] > value) p2flag = INSIDE;
          else p2flag = OVERLAP;
        } else {
          if (p1[dim] < value) p1flag = INSIDE;
          else if (p1[dim] > value) p1flag = OUTSIDE;
          else p1flag = OVERLAP;
          if (p2[dim] < value) p2flag = INSIDE;
          else if (p2[dim] > value) p2flag = OUTSIDE;
          else p2flag = OVERLAP;
        }

        // if both OUTSIDE or one OUTSIDE and other ON, delete edge
        // if both INSIDE or one INSIDE and other ON or both ON, keep as-is
        // if one INSIDE and one OUTSIDE, replace OUTSIDE pt with clip pt

#ifdef VERBOSE
        /*
        if (id == VERBOSE_ID && ivert == 1) {
          printf("EDGE %d %d: %d %d: %d %d: %d\n",iedge,idir,p1flag,p2flag,
                 edge->verts[0],edge->verts[1],edge->prev[0]);
        }
        */
#endif

        if (p1flag == OUTSIDE) {
          if (p2flag == OUTSIDE || p2flag == OVERLAP) edge_remove(edge,idir);
          else {
            if (idir == 0) between(p1,p2,dim,value,edge->p1);
            else between(p1,p2,dim,value,edge->p2);
            edge->clipped = 1;
          }
        } else if (p1flag == INSIDE) {
          if (p2flag == OUTSIDE) {
            if (idir == 0) between(p1,p2,dim,value,edge->p2);
            else between(p1,p2,dim,value,edge->p1);
            edge->clipped = 1;
          }
        } else {
          if (p2flag == OUTSIDE) edge_remove(edge,idir);
        }

        iedge = edge->next[idir];
        idir = edge->dirnext[idir];
      }

#ifdef VERBOSE
      /*
      if (id == VERBOSE_ID) {
        char str[24];
        sprintf(str,"Partial FACE %d %d\n",iface,ivert);
        print_bpg(str);
      }
      */
#endif

      // loop over edges in vertex again
      // iedge = this edge, jedge = next edge
      // p1 = last pt in iedge, pt = first pt in jedge
      // if p1 != p2, add edge between them

      edges.grow(edges.n + verts[ivert].nedge);
      iedge = verts[ivert].first;
      idir = verts[ivert].dirfirst;

      for (i = 0; i < verts[ivert].nedge; i++) {
        edge = &edges[iedge];
        jedge = edge->next[idir];
        jdir = edge->dirnext[idir];
        if (jedge < 0) {
          jedge = verts[ivert].first;
          jdir = verts[ivert].dirfirst;
        }

        if (idir == 0) p1 = edge->p2;
        else p1 = edge->p1;
        if (jdir == 0) p2 = edges[jedge].p1;
        else p2 = edges[jedge].p2;

        if (!samepoint(p1,p2)) {
          n = edges.n++;
          newedge = &edges[n];
          newedge->style = CTRI;
          newedge->nvert = 0;
          newedge->verts[0] = newedge->verts[1] = -1;
          memcpy(newedge->p1,p1,3*sizeof(double));
          memcpy(newedge->p2,p2,3*sizeof(double));
          // convert jedge back to -1 for last vertex
          if (jedge == verts[ivert].first) jedge = -1;
          edge_insert(n,0,ivert,iedge,idir,jedge,jdir);
          i++;
        }

        iedge = jedge;
        idir = jdir;
      }
    }

#ifdef VERBOSE
    /*
    if (id == VERBOSE_ID) {
      char str[24];
      sprintf(str,"After FACE %d\n",iface);
      print_bpg(str);
    }
    */
#endif
  }
}

/* ----------------------------------------------------------------------
   adjust the collection of clipped triangles (vertices)
   discard if clipped tri is a single point, increment touchcount
     touchmark = corner point marking inferred from touching tri orientations
   discard if grazes cell with outward normal, increment grazecount
   if all clipped tris are discarded
     set and return empty, touchcount, grazecount, touchmark
------------------------------------------------------------------------- */

void Cut3d::clip_adjust()
{
  int nvert,nedge,nface1,nface2;
  int faces1[6],faces2[6];
  double pboth[3],move1[3],move2[3];
  double *p1,*p2,*p3;
  Edge *edge;

#ifdef VERBOSE
  if (id == VERBOSE_ID) print_bpg("BPG after initial clipping");
#endif

  // epsilon = EPSEDGE fraction of largest cell dimension

  epsilon = EPSEDGE*(hi[0]-lo[0]);
  epsilon = MAX(epsilon,EPSEDGE*(hi[1]-lo[1]));
  epsilon = MAX(epsilon,EPSEDGE*(hi[2]-lo[2]));

  // collapse edges shorter than epsilon to a single point (so will be removed)
  // one or both of the points should be on cell faces
  // tiny edges can occur for several reasons:
  // (1) pre-clip end point is within epsilon of face
  // (2) a tri barely grazes cell edge or corner
  // (3) clipping round-off produces tiny edges near cell faces

  nedge = edges.n;

  for (int iedge = 0; iedge < nedge; iedge++) {
    if (!edges[iedge].active) continue;

    edge = &edges[iedge];
    p1 = edge->p1;
    p2 = edge->p2;
    double dx = p1[0]-p2[0];
    double dy = p1[1]-p2[1];
    double dz = p1[2]-p2[2];
    double edgelen = sqrt(dx*dx+dy*dy+dz*dz);

    if (edgelen < epsilon) {
      //printf("TINY EDGE id %ld nsurf %d i/nedge %d %d len %g eps %g\n",
      //       id,nsurf,iedge,nedge,edgelen,epsilon);
      ntiny++;

      nface1 = on_faces(p1,faces1);
      nface2 = on_faces(p2,faces2);

      // set both p1 and p2 to same pboth
      // if both pts are interior (should not happen), pboth = p1
      // if only one pt X is on a face, pboth = X
      // if both pts are on one or more faces:
      //   push both to face(s), recalculate on_faces()
      //   if pt X is on more faces, pboth = X, else pboth = p1

      if (!nface1 && !nface2) {
        memcpy(pboth,p1,3*sizeof(double));
        //printf("INTERIOR EDGE %ld %d %d %g %g\n",
        //       id,iedge,nedge,edgelen,epsilon);
      } else if (nface1 && !nface2) {
        memcpy(pboth,p1,3*sizeof(double));
      } else if (nface2 && !nface1) {
        memcpy(pboth,p2,3*sizeof(double));
      } else {
        memcpy(move1,p1,3*sizeof(double));
        memcpy(move2,p2,3*sizeof(double));
        move_to_faces(move1);
        move_to_faces(move2);
        nface1 = on_faces(move1,faces1);
        nface2 = on_faces(move2,faces2);
        if (nface2 > nface1) memcpy(pboth,move2,3*sizeof(double));
        else memcpy(pboth,move1,3*sizeof(double));
      }

      // set all points that are same as old p1 or p2 to pboth
      // reset first for all jedge != iedge, then reset iedge

      for (int jedge = 0; jedge < nedge; jedge++) {
        if (!edges[jedge].active) continue;
        if (jedge == iedge) continue;

        if (samepoint(edges[jedge].p1,p1))
          memcpy(edges[jedge].p1,pboth,3*sizeof(double));
        if (samepoint(edges[jedge].p2,p1))
          memcpy(edges[jedge].p2,pboth,3*sizeof(double));

        if (samepoint(edges[jedge].p1,p2))
          memcpy(edges[jedge].p1,pboth,3*sizeof(double));
        if (samepoint(edges[jedge].p2,p2))
          memcpy(edges[jedge].p2,pboth,3*sizeof(double));
      }

      memcpy(edges[iedge].p1,pboth,3*sizeof(double));
      memcpy(edges[iedge].p2,pboth,3*sizeof(double));
    }
  }

#ifdef VERBOSE
  if (id == VERBOSE_ID) print_bpg("BPG after tiny edge collapse");
#endif

  // remove zero-length edges

  nedge = edges.n;

  for (int iedge = 0; iedge < nedge; iedge++) {
    if (!edges[iedge].active) continue;
    edge = &edges[iedge];
    if (samepoint(edge->p1,edge->p2)) edge_remove(edge);
  }

#ifdef VERBOSE
  if (id == VERBOSE_ID) print_bpg("BPG after remove zero-length edges");
#endif

  // remove vertices (triangles) which now have less than 3 edges
  // do this after deleting zero-length edges so vertices are updated
  // removals should have 2 or 0 edges, no verts should have 1 edge
  // if all triangles are removed, BPG will be empty,
  //   which will result in cell corner pts being left UNKNOWN in split()
  // to try and avoid this, tally inside/outside for all removed tris
  //   tri is "outside" if it implies cell is outside the surf (in the flow)
  //   tri is "inside" if it implies cell is inside the surf (not in the flow)
  //   some tris may not follow this rule, but most should
  // outside = tri norm from tri ctr points towards cell ctr
  // inside = tri norm from tri ctr points away from cell ctr
  // cbox = cell center pt, ctri = triangle center pt, t2b = cbox-ctri

  touchcount = 0;
  grazecount = 0;

  int noutside = 0;
  int ninside = 0;

  double cbox[3],ctri[3],t2b[3];

  nvert = verts.n;

  for (int ivert = 0; ivert < nvert; ivert++)
    if (verts[ivert].nedge <= 2) {
      touchcount++;
      cbox[0] = 0.5*(lo[0]+hi[0]);
      cbox[1] = 0.5*(lo[1]+hi[1]);
      cbox[2] = 0.5*(lo[2]+hi[2]);
      int itri = surfs[verts[ivert].label];
      p1 = surf->tris[itri].p1;
      p2 = surf->tris[itri].p2;
      p3 = surf->tris[itri].p3;
      ctri[0] = (p1[0]+p2[0]+p3[0])/3.0;
      ctri[1] = (p1[1]+p2[1]+p3[1])/3.0;
      ctri[2] = (p1[2]+p2[2]+p3[2])/3.0;
      MathExtra::sub3(cbox,ctri,t2b);
      double dot = MathExtra::dot3(verts[ivert].norm,t2b);
      if (dot > 0.0) noutside++;
      if (dot < 0.0) ninside++;
      vertex_remove(&verts[ivert]);
    }

  // discard clipped tri if lies on a cell face w/ normal out of cell
  // increment grazecount in this case

  for (int ivert = 0; ivert < nvert; ivert++) {
    if (!verts[ivert].active) continue;
    if (grazing(&verts[ivert])) {
      vertex_remove(&verts[ivert]);
      grazecount++;
    }
  }

  // remove edges which now have no vertices

  for (int iedge = 0; iedge < nedge; iedge++) {
    if (!edges[iedge].active) continue;
    if (edges[iedge].nvert == 0) edges[iedge].active = 0;
  }

  // set BPG empty flag if no active vertices

  empty = 1;
  for (int ivert = 0; ivert < nvert; ivert++)
    if (verts[ivert].active) {
      empty = 0;
      break;
    }

  // if no lines, set touchmark which will be used to mark corner points
  // only set touchmark if all deleted tris had same orientation
  // NOTE: setting touchmark based on orientation of just majority
  //   triggered later corner marking error for cone_test/in.cone problem

  touchmark = UNKNOWN;
  if (empty) {
    if (ninside && noutside == 0) touchmark = INSIDE;
    else if (noutside && ninside == 0) touchmark = OUTSIDE;
  }
}

/* ----------------------------------------------------------------------
   compute volume of vertices
   when called, only clipped triangles exist
------------------------------------------------------------------------- */

void Cut3d::ctri_volume()
{
  int i,iedge,idir,nedge;
  double zarea,volume;
  double *p0,*p1,*p2;
  Edge *edge;

  int nvert = verts.n;
  for (int ivert = 0; ivert < nvert; ivert++) {
    if (!verts[ivert].active) continue;
    iedge = verts[ivert].first;
    idir = verts[ivert].dirfirst;
    nedge = verts[ivert].nedge;

    if (idir == 0) p0 = edges[iedge].p1;
    else p0 = edges[iedge].p2;

    volume = 0.0;

    for (i = 0; i < nedge; i++) {
      edge = &edges[iedge];

      // compute projected volume of a convex polygon to zlo face
      // split polygon into triangles
      // each tri makes a tri-capped volume with zlo face
      // zarea = area of oriented tri projected into z plane
      // volume based on height of z midpt of tri above zlo face

      if (idir == 0) {
        p1 = edge->p1;
        p2 = edge->p2;
      } else {
        p1 = edge->p2;
        p2 = edge->p1;
      }
      zarea = 0.5 * ((p1[0]-p0[0])*(p2[1]-p0[1]) - (p1[1]-p0[1])*(p2[0]-p0[0]));
      volume -= zarea * ((p0[2]+p1[2]+p2[2])/3.0 - lo[2]);

      iedge = edge->next[idir];
      idir = edge->dirnext[idir];
    }

    verts[ivert].volume = volume;
  }
}

/* ----------------------------------------------------------------------
   assign all singlet edges to faces (0-5)
   singlet edge must be on one or two faces, two if on cell edge
   if along cell edge, assign to one of two faces based on
     which has larger dot product of its inward face norm
     and the norm of the tri containing the edge
------------------------------------------------------------------------- */

int Cut3d::edge2face()
{
  int n,iface,nface,ivert;
  int faces[6];
  double dot0,dot1;
  double norm_inward[3];
  double *trinorm;
  Edge *edge;

  // insure each facelist has sufficient length

  int nedge = edges.n;
  for (int i = 0; i < 6; i++) {
    facelist[i].grow(nedge);
    facelist[i].n = 0;
  }

  // loop over edges, assign singlets to exactly one face

  for (int iedge = 0; iedge < nedge; iedge++) {
    if (!edges[iedge].active) continue;
    if (edges[iedge].nvert == 3) continue;
    edge = &edges[iedge];

    nface = which_faces(edge->p1,edge->p2,faces);
    if (nface == 0) {
      //printf("ERROR RETURN id %ld nedge %d iedge %d pt1 "
      //       "%20.16g %20.16g %20.16g %20.16g %20.16g %20.16g\n",
      //       id,nedge,iedge,
      //       edge->p1[0],edge->p1[1],edge->p1[2],
      //       edge->p2[0],edge->p2[1],edge->p2[2]);
      return 2;
    }

    else if (nface == 1) iface = faces[0];

    else if (nface == 2) {
      if (edge->nvert == 1) ivert = edge->verts[0];
      else ivert = edge->verts[1];
      trinorm = verts[ivert].norm;

      iface = faces[0];
      norm_inward[0] = norm_inward[1] = norm_inward[2] = 0.0;
      if (iface % 2) norm_inward[iface/2] = -1.0;
      else norm_inward[iface/2] = 1.0;
      dot0 = norm_inward[0]*trinorm[0] + norm_inward[1]*trinorm[1] +
        norm_inward[2]*trinorm[2];
      if (dot0 > 0.0) iface = faces[1];

    } else return 3;

    n = facelist[iface].n;
    facelist[iface][n++] = iedge;
    facelist[iface].n = n;
  }

  return 0;
}

/* ----------------------------------------------------------------------
   build a 2d CLINES data structure
   from all singlet edges assigned to iface (0-5)
   order pts in edge for tri traversing edge in forward order
   flip edge if in a flip face = faces 0,3,4
   edge label in clines = edge index in BPG
------------------------------------------------------------------------- */

void Cut3d::edge2clines(int iface)
{
  int iedge;
  double *p1,*p2;
  Edge *edge;
  Cut2d::Cline *cline;

  MyVec<Cut2d::Cline> *clines = &cut2d->clines;

  int flip = 0;
  if (iface == 0 || iface == 3 || iface == 4) flip = 1;

  int nline = facelist[iface].n;
  clines->n = 0;
  clines->grow(nline);

  for (int i = 0; i < nline; i++) {
    iedge = facelist[iface][i];
    edge = &edges[iedge];
    if (edge->nvert == 1) {
      p1 = edge->p1;
      p2 = edge->p2;
    } else {
      p1 = edge->p2;
      p2 = edge->p1;
    }
    cline = &(*clines)[i];
    cline->line = iedge;
    if (flip) {
      compress2d(iface,p1,cline->y);
      compress2d(iface,p2,cline->x);
    } else {
      compress2d(iface,p1,cline->x);
      compress2d(iface,p2,cline->y);
    }
  }

  clines->n = nline;
}

/* ----------------------------------------------------------------------
   add one or more face polygons as vertices to BPG
   have to convert pts computed by cut2d back into 3d pts on face
------------------------------------------------------------------------- */

int Cut3d::add_face_pgons(int iface)
{
  int iloop,mloop,nloop,ipt,mpt,npt;
  int iedge,dir,prev,dirprev;
  double p1[3],p2[3];
  Vertex *vert;
  Edge *edge;
  Cut2d::PG *pg;
  Cut2d::Loop *loop;
  Cut2d::Point *p12d,*p22d;

  MyVec<Cut2d::PG> *pgs = &cut2d->pgs;
  MyVec<Cut2d::Loop> *loops = &cut2d->loops;
  MyVec<Cut2d::Point> *points = &cut2d->points;

  int flip = 0;
  if (iface == 0 || iface == 3 || iface == 4) flip = 1;

  double value;
  int dim = iface / 2;
  int lohi = iface % 2;
  if (lohi == 0) value = lo[dim];
  else value = hi[dim];

  int npg = pgs->n;
  int nvert = verts.n;
  verts.grow(nvert+npg);

  for (int ipg = 0; ipg < npg; ipg++) {
    pg = &(*pgs)[ipg];

    vert = &verts[nvert];
    vert->active = 1;
    vert->style = FACEPGON;
    vert->label = iface;
    if (iface == 5) vert->volume = pg->area * (hi[2]-lo[2]);
    else vert->volume = 0.0;
    vert->nedge = 0;
    vert->norm = NULL;

    prev = -1;
    dirprev = -1;

    nloop = pg->n;
    mloop = pg->first;
    for (iloop = 0; iloop < nloop; iloop++) {
      loop = &(*loops)[mloop];
      npt = loop->n;
      mpt = loop->first;
      edges.grow(edges.n + npt);

      for (ipt = 0; ipt < npt; ipt++) {
        p12d = &(*points)[mpt];
        mpt = p12d->next;
        p22d = &(*points)[mpt];
        expand2d(iface,value,p12d->x,p1);
        expand2d(iface,value,p22d->x,p2);

        // edge was from a CTRI vertex
        // match in opposite order that CTRI vertex matched it

        if (p12d->type == ENTRY || p12d->type == TWO) {
          iedge = p12d->line;
          edge = &edges[iedge];
          edge->style = CTRIFACE;
          if (edge->nvert == 1) dir = 1;
          else dir = 0;
          edge_insert(iedge,dir,nvert,prev,dirprev,-1,-1);
          prev = iedge;
          dirprev = dir;
          continue;
        }

        // face edge not from a CTRI
        // unflip edge if in a flip face

        if (flip) iedge = findedge(p2,p1,0,dir);
        else iedge = findedge(p1,p2,0,dir);
        if (iedge == -2) return 1;

        if (iedge >= 0) {
          edge_insert(iedge,dir,nvert,prev,dirprev,-1,-1);
          prev = iedge;
          dirprev = 1;
          continue;
        }

        iedge = edges.n++;
        edge = &edges[iedge];
        edge->style = FACEPGON;
        edge->nvert = 0;
        edge->verts[0] = edge->verts[1] = -1;
        if (flip) {
          memcpy(edge->p1,p2,3*sizeof(double));
          memcpy(edge->p2,p1,3*sizeof(double));
        } else {
          memcpy(edge->p1,p1,3*sizeof(double));
          memcpy(edge->p2,p2,3*sizeof(double));
        }
        dir = 0;
        edge_insert(iedge,dir,nvert,prev,dirprev,-1,-1);
        prev = iedge;
        dirprev = 0;
      }
      mloop = loop->next;
    }

    nvert++;
  }

  verts.n = nvert;
  return 0;
}

/* ----------------------------------------------------------------------
   add an entire cell face as vertex to BPG
   if outerflag2d = 0, create new vertex
   else face polygon already exists, so add edges to it
   caller sets outerflag2d if cut2d requires adding perimeter face edges
------------------------------------------------------------------------- */

int Cut3d::add_face(int iface, double *lo2d, double *hi2d)
{
  int i,j,iedge,dir,prev,dirprev;
  double p1[3],p2[3];
  Vertex *vert;
  Edge *edge;

  int nvert = verts.n++;
  verts.grow(nvert + 1);
  vert = &verts[nvert];
  vert->active = 1;
  vert->style = FACE;
  vert->label = iface;
  if (iface == 5)
    vert->volume = (hi[0]-lo[0]) * (hi[1]-lo[1]) * (hi[2]-lo[2]);
  else vert->volume = 0.0;
  vert->nedge = 0;
  vert->norm = NULL;

  double value;
  int dim = iface / 2;
  int lohi = iface % 2;
  if (lohi == 0) value = lo[dim];
  else value = hi[dim];

  // usual ordering of points in face as LL,LR,UR,UL
  // flip order if in a flip face

  int flip = 0;
  if (iface == 0 || iface == 3 || iface == 4) flip = 1;

  double cpts[4][2];

  if (flip) {
    cpts[0][0] = lo2d[0]; cpts[0][1] = lo2d[1];
    cpts[1][0] = lo2d[0]; cpts[1][1] = hi2d[1];
    cpts[2][0] = hi2d[0]; cpts[2][1] = hi2d[1];
    cpts[3][0] = hi2d[0]; cpts[3][1] = lo2d[1];
  } else {
    cpts[0][0] = lo2d[0]; cpts[0][1] = lo2d[1];
    cpts[1][0] = hi2d[0]; cpts[1][1] = lo2d[1];
    cpts[2][0] = hi2d[0]; cpts[2][1] = hi2d[1];
    cpts[3][0] = lo2d[0]; cpts[3][1] = hi2d[1];
  }

  if (vert->nedge) {
    prev = vert->last;
    dirprev = vert->dirlast;
  } else {
    prev = -1;
    dirprev = -1;
  }

  edges.grow(edges.n + 4);

  for (i = 0; i < 4; i++) {
    j = i+1;
    if (j == 4) j = 0;
    expand2d(iface,value,&cpts[i][0],p1);
    expand2d(iface,value,&cpts[j][0],p2);
    iedge = findedge(p1,p2,1,dir);
    if (iedge == -2) return 1;

    if (iedge >= 0) {
      edge_insert(iedge,dir,nvert,prev,dirprev,-1,-1);
      prev = iedge;
      dirprev = 1;
      continue;
    }

    iedge = edges.n++;
    edge = &edges[iedge];
    edge->style = vert->style;
    edge->nvert = 0;
    edge->verts[0] = edge->verts[1] = -1;
    memcpy(edge->p1,p1,3*sizeof(double));
    memcpy(edge->p2,p2,3*sizeof(double));
    dir = 0;
    edge_insert(iedge,dir,nvert,prev,dirprev,-1,-1);
    prev = iedge;
    dirprev = 0;
  }

  return 0;
}

/* ----------------------------------------------------------------------
   remove any FACE vertices with one or more unconnected edges
   unconnected means the edge is only part of this vertex
   mark vertex as inactive, decrement their edges,
     also set edges inactive if they are no longer attached to vertices
   iterate twice since another face may become unconnected
------------------------------------------------------------------------- */

void Cut3d::remove_faces()
{
  int i,ivert,iedge,dir;
  Vertex *vert;
  Edge *edge;

  int nvert = verts.n;

  for (int iter = 0; iter < 2; iter++)
    for (ivert = 0; ivert < nvert; ivert++) {
      if (!verts[ivert].active) continue;
      if (verts[ivert].style != FACE) continue;
      vert = &verts[ivert];

      iedge = vert->first;
      dir = vert->dirfirst;
      for (i = 0; i < 4; i++) {
        edge = &edges[iedge];
        if (edge->nvert == 1 || edge->nvert == 2) break;
        iedge = edge->next[dir];
        dir = edge->dirnext[dir];
      }
      if (i < 4) vertex_remove(vert);
    }
}

/* ----------------------------------------------------------------------
   check BPG for consistency
   vertices have 3 or more unique edges that point back to it
   edges have 2 unique vertices
------------------------------------------------------------------------- */

int Cut3d::check()
{
  int i,iedge,dir,nedge,last,dirlast;
  Vertex *vert;
  Edge *edge;

  // mark all edges as unclipped
  // use for detecting duplicate edges in same vertex

  nedge = edges.n;
  for (iedge = 0; iedge < nedge; iedge++)
    if (edges[iedge].active) edges[iedge].clipped = 0;

  // check vertices
  // for each vertex: mark edges as see them, unmark all edges at end

  int nvert = verts.n;
  for (int ivert = 0; ivert < nvert; ivert++) {
    if (!verts[ivert].active) continue;
    vert = &verts[ivert];
    if (vert->nedge < 3) return 11;

    nedge = vert->nedge;
    iedge = vert->first;
    dir = vert->dirfirst;

    for (i = 0; i < nedge; i++) {
      edge = &edges[iedge];
      if (!edge->active) return 12;
      if (edge->verts[dir] != ivert) return 13;
      if (edge->clipped) return 14;
      edge->clipped = 1;
      last = iedge;
      dirlast = dir;
      iedge = edge->next[dir];
      dir = edge->dirnext[dir];
    }

    if (last != vert->last || dirlast != vert->dirlast) return 15;

    iedge = vert->first;
    dir = vert->dirfirst;
    for (i = 0; i < nedge; i++) {
      edge = &edges[iedge];
      edge->clipped = 0;
      iedge = edge->next[dir];
      dir = edge->dirnext[dir];
    }
  }

  // check edges

  nedge = edges.n;
  for (int iedge = 0; iedge < nedge; iedge++) {
    if (!edges[iedge].active) continue;
    edge = &edges[iedge];

    if (edge->nvert != 3) return 16;
    if (edge->verts[0] == edge->verts[1]) return 17;
    if (edge->verts[0] >= nvert || !verts[edge->verts[0]].active) return 18;
    if (edge->verts[1] >= nvert || !verts[edge->verts[1]].active) return 19;
  }

  return 0;
}

/* ----------------------------------------------------------------------
   convert BPG into simple closed polyhedra, not nested
   walk BPG from any unused vertex, flagging vertices as used
   stack is list of new vertices to process
   loop over edges of pgon, add its unused neighbors to stack
   when stack is empty, loop is closed
   accumulate volume of polyhedra as walk it from volume of each vertex
------------------------------------------------------------------------- */

void Cut3d::walk()
{
  int flag,ncount,ivert,firstvert,iedge,dir,nedge,prev;
  double volume;
  Vertex *vert;
  Edge *edge;

  // used = 0/1 flag for whether a vertex is already part of a loop
  // only active vertices are eligible

  int nvert = verts.n;
  used.grow(nvert);
  for (int ivert = 0; ivert < nvert; ivert++) {
    if (verts[ivert].active) used[ivert] = 0;
    else used[ivert] = 1;
  }
  used.n = nvert;

  // stack = list of vertex indices to process
  // max size = # of vertices

  stack.grow(nvert);
  int nstack = 0;

  // iterate over all vertices
  // start a loop at any unused vertex
  // add more vertices to loop via stack
  // check all neighbor vertices via shared edges
  // if neighbor vertex is unused, add to stack
  // stop when stack is empty

  int nloop = 0;

  for (int i = 0; i < nvert; i++) {
    if (used[i]) continue;
    volume = 0.0;
    flag = INTERIOR;
    ncount = 0;

    stack[0] = firstvert = i;
    nstack = 1;
    used[i] = 1;
    prev = -1;

    while (nstack) {
      nstack--;
      ivert = stack[nstack];
      ncount++;

      vert = &verts[ivert];
      if (vert->style != CTRI) flag = BORDER;
      volume += vert->volume;

      nedge = vert->nedge;
      iedge = vert->first;
      dir = vert->dirfirst;

      for (i = 0; i < nedge; i++) {
        edge = &edges[iedge];
        if (!used[edge->verts[0]]) {
          stack[nstack++] = edge->verts[0];
          used[edge->verts[0]] = 1;
        }
        if (!used[edge->verts[1]]) {
          stack[nstack++] = edge->verts[1];
          used[edge->verts[1]] = 1;
        }
        iedge = edge->next[dir];
        dir = edge->dirnext[dir];
      }

      if (prev >= 0) verts[prev].next = ivert;
      prev = ivert;
    }
    vert->next = -1;

    loops.grow(nloop+1);
    loops[nloop].volume = volume;
    loops[nloop].flag = flag;
    loops[nloop].n = ncount;
    loops[nloop].first = firstvert;
    nloop++;
  }

  loops.n = nloop;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

int Cut3d::loop2ph()
{
  int positive = 0;
  int negative = 0;

  int nloop = loops.n;

  for (int i = 0; i < nloop; i++) {
    if (loops[i].volume > 0.0) positive++;
    else negative++;
  }

  // if no positive vols, cell is entirely inside the surf, caller handles it
  // this can happen due to epsilon size polyhedron(s)
  // e.g. when a tri barely cuts off a cell corner

  if (positive == 0) return 4;

  // do not allow mulitple positive with one or more negative
  // too difficult to figure out which positive each negative is inside of

  if (positive > 1 && negative) return 5;

  // positive = 1 means 1 PH with vol = sum of all pos/neg loops
  // positive > 1 means each loop is a PH

  phs.grow(positive);

  if (positive == 1) {
    double volume = 0.0;
    for (int i = 0; i < nloop; i++) {
      volume += loops[i].volume;
      loops[i].next = i+1;
    }
    loops[nloop-1].next = -1;

    if (volume < 0.0) return 6;

    phs[0].volume = volume;
    phs[0].n = nloop;
    phs[0].first = 0;

  } else {
    for (int i = 0; i < nloop; i++) {
      phs[i].volume = loops[i].volume;
      phs[i].n = 1;
      phs[i].first = i;
      loops[i].next = -1;
    }
  }

  phs.n = positive;
  return 0;
}

/* ----------------------------------------------------------------------
   assign each tri index in list to one of the split cells in PH
   return surfmap[i] = which PH the Ith tri index is assigned to
   set surfmap[i] = -1 if the tri did not end up in a PH
     could have been discarded in clip_tris()
     due to touching cell or lying along a cell edge or face
------------------------------------------------------------------------- */

void Cut3d::create_surfmap(int *surfmap)
{
  for (int i = 0; i < nsurf; i++) surfmap[i] = -1;

  int iloop,nloop,mloop,ivert,nvert,mvert;

  int nph = phs.n;
  for (int iph = 0; iph < nph; iph++) {
    nloop = phs[iph].n;
    mloop = phs[iph].first;
    for (iloop = 0; iloop < nloop; iloop++) {
      nvert = loops[mloop].n;
      mvert = loops[mloop].first;
      for (ivert = 0; ivert < nvert; ivert++) {
        if (verts[mvert].style == CTRI || verts[mvert].style == CTRIFACE)
          surfmap[verts[mvert].label] = iph;
        mvert = verts[mvert].next;
      }
      mloop = loops[mloop].next;
    }
  }
}

/* ----------------------------------------------------------------------
   find a surf point that is inside or on the boundary of the current cell
   for external surfs and cells already been flagged as a split cell
   surfmap = sub-cell index each surf is part of (-1 if not eligible)
   return xsplit = coords of point
   return xsub = sub-cell index the chosen surf is in
------------------------------------------------------------------------- */

int Cut3d::split_point_explicit(int *surfmap, double *xsplit, int &xsub)
{
  int itri;
  double *x1,*x2,*x3;

  Surf::Tri *tris = surf->tris;

  // if end pt of any line with non-negative surfmap is in/on cell, return

  for (int i = 0; i < nsurf; i++) {
    if (surfmap[i] < 0) continue;
    itri = surfs[i];
    x1 = tris[itri].p1;
    x2 = tris[itri].p2;
    x3 = tris[itri].p3;
    if (ptflag(x1) != EXTERIOR) {
      xsplit[0] = x1[0]; xsplit[1] = x1[1]; xsplit[2] = x1[2];
      xsub = surfmap[i];
      return 0;
    }
    if (ptflag(x2) != EXTERIOR) {
      xsplit[0] = x2[0]; xsplit[1] = x2[1]; xsplit[2] = x2[2];
      xsub = surfmap[i];
      return 0;
    }
    if (ptflag(x3) != EXTERIOR) {
      xsplit[0] = x3[0]; xsplit[1] = x3[1]; xsplit[2] = x3[2];
      xsub = surfmap[i];
      return 0;
    }
  }

  // clip 1st line with non-negative surfmap to cell, and return clip point

  for (int i = 0; i < nsurf; i++) {
    if (surfmap[i] < 0) continue;
    itri = surfs[i];
    x1 = tris[itri].p1;
    x2 = tris[itri].p2;
    x3 = tris[itri].p3;
    clip(x1,x2,x3);
    xsplit[0] = path1[0][0]; xsplit[1] = path1[0][1]; xsplit[2] = path1[0][2];
    xsub = surfmap[i];
    return 0;
  }

  // error return

  return 7;
}

/* ----------------------------------------------------------------------
   find a surf point that is inside or on the boundary of the current cell
   for implicit surfs and cells already been flagged as a split cell
   surfmap = sub-cell index each surf is part of (-1 if not eligible)
   return xsplit = coords of point
   return xsub = sub-cell index the chosen surf is in
------------------------------------------------------------------------- */

int Cut3d::split_point_implicit(int *surfmap, double *xsplit, int &xsub)
{
  Surf::Tri *tris = surf->tris;

  // i = 1st surf with non-negative surfmap

  int i = 0;
  while (surfmap[i] < 0 && i < nsurf) i++;
  if (i == nsurf) return 7;

  // xsplit = center point of triangle wholly contained in cell

  int itri = surfs[i];
  double onethird = 1.0/3.0;
  xsplit[0] = onethird * (tris[itri].p1[0] + tris[itri].p2[0] + tris[itri].p3[0]);
  xsplit[1] = onethird * (tris[itri].p1[1] + tris[itri].p2[1] + tris[itri].p3[1]);
  xsplit[2] = onethird * (tris[itri].p1[2] + tris[itri].p2[2] + tris[itri].p3[2]);

  xsub = surfmap[i];

  return 0;
}

/* ----------------------------------------------------------------------
   insert edge IEDGE in DIR for ivert
   also update vertex info for added edge
------------------------------------------------------------------------- */

void Cut3d::edge_insert(int iedge, int dir, int ivert,
                        int iprev, int dirprev, int inext, int dirnext)
{
  Edge *edge = &edges[iedge];

  if (dir == 0) {
    edge->nvert += 1;
    edge->verts[0] = ivert;
  } else {
    edge->nvert += 2;
    edge->verts[1] = ivert;
  }

  edge->active = 1;
  edge->clipped = 0;

  // set prev/next pointers for doubly linked list of edges

  edge->next[dir] = inext;
  edge->prev[dir] = iprev;

  if (inext >= 0) {
    edge->dirnext[dir] = dirnext;
    Edge *next = &edges[inext];
    next->prev[dirnext] = iedge;
    next->dirprev[dirnext] = dir;
  } else edge->dirnext[dir] = -1;

  if (iprev >= 0) {
    edge->dirprev[dir] = dirprev;
    Edge *prev = &edges[iprev];
    prev->next[dirprev] = iedge;
    prev->dirnext[dirprev] = dir;
  } else edge->dirprev[dir] = -1;

  // add edge info to owning vertex

  verts[ivert].nedge++;
  if (iprev < 0) {
    verts[ivert].first = iedge;
    verts[ivert].dirfirst = dir;
  }
  if (inext < 0) {
    verts[ivert].last = iedge;
    verts[ivert].dirlast = dir;
  }
}

/* ----------------------------------------------------------------------
   complete edge removal in both dirs
   will leave edge marked inactive
------------------------------------------------------------------------- */

void Cut3d::edge_remove(Edge *edge)
{
  if (edge->verts[0] >= 0) edge_remove(edge,0);
  if (edge->verts[1] >= 0) edge_remove(edge,1);
}

/* ----------------------------------------------------------------------
   edge removal in DIR
   also update vertex info for removed edge
   mark edge inactive if its nvert -> 0
------------------------------------------------------------------------- */

void Cut3d::edge_remove(Edge *edge, int dir)
{
  int ivert = edge->verts[dir];
  edge->verts[dir] = -1;
  if (dir == 0) edge->nvert--;
  else edge->nvert -= 2;
  if (edge->nvert == 0) edge->active = 0;

  // reset prev/next pointers for doubly linked list to skip this edge

  if (edge->prev[dir] >= 0) {
    Edge *prev = &edges[edge->prev[dir]];
    int dirprev = edge->dirprev[dir];
    prev->next[dirprev] = edge->next[dir];
    prev->dirnext[dirprev] = edge->dirnext[dir];
  }

  if (edge->next[dir] >= 0) {
    Edge *next = &edges[edge->next[dir]];
    int dirnext = edge->dirnext[dir];
    next->prev[dirnext] = edge->prev[dir];
    next->dirprev[dirnext] = edge->dirprev[dir];
  }

  // update vertex for removal of this edge

  verts[ivert].nedge--;
  if (edge->prev[dir] < 0) {
    verts[ivert].first = edge->next[dir];
    verts[ivert].dirfirst = edge->dirnext[dir];
  }
  if (edge->next[dir] < 0) {
    verts[ivert].last = edge->prev[dir];
    verts[ivert].dirlast = edge->dirprev[dir];
  }
}

/* ----------------------------------------------------------------------
   remove a vertex and all edges it includes
------------------------------------------------------------------------- */

void Cut3d::vertex_remove(Vertex *vert)
{
  Edge *edge;

  vert->active = 0;

  int iedge = vert->first;
  int dir = vert->dirfirst;
  int nedge = vert->nedge;

  for (int i = 0; i < nedge; i++) {
    edge = &edges[iedge];
    if (dir == 0) edge->nvert--;
    else edge->nvert -= 2;
    if (edge->nvert == 0) edge->active = 0;
    edge->verts[dir] = -1;
    iedge = edge->next[dir];
    dir = edge->dirnext[dir];
  }
}

/* ----------------------------------------------------------------------
   a planar polygon is grazing if it lies entirely in plane of any face of cell
   and its normal is outward with respect to cell
   return 1 if grazing else 0
------------------------------------------------------------------------- */

int Cut3d::grazing(Vertex *vert)
{
  int count[6];
  double *p;
  Edge *edge;

  int iedge = vert->first;
  int idir = vert->dirfirst;
  int nedge = vert->nedge;

  count[0] = count[1] = count[2] = count[3] = count[4] = count[5] = 0;

  for (int i = 0; i < nedge; i++) {
    edge = &edges[iedge];
    if (idir == 0) p = edge->p1;
    else p = edge->p2;

    if (p[0] == lo[0]) count[0]++;
    if (p[0] == hi[0]) count[1]++;
    if (p[1] == lo[1]) count[2]++;
    if (p[1] == hi[1]) count[3]++;
    if (p[2] == lo[2]) count[4]++;
    if (p[2] == hi[2]) count[5]++;

    iedge = edge->next[idir];
    idir = edge->dirnext[idir];
  }

  double *norm = vert->norm;
  if (count[0] == nedge && norm[0] < 0.0) return 1;
  if (count[1] == nedge && norm[0] > 0.0) return 1;
  if (count[2] == nedge && norm[1] < 0.0) return 1;
  if (count[3] == nedge && norm[1] > 0.0) return 1;
  if (count[4] == nedge && norm[2] < 0.0) return 1;
  if (count[5] == nedge && norm[2] > 0.0) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   identify which cell faces point P is on
   return list of face IDs (0-5)
   list length can be 0,1,2
------------------------------------------------------------------------- */

int Cut3d::on_faces(double *p, int *faces)
{
  int n = 0;
  if (p[0] == lo[0]) faces[n++] = 0;
  if (p[0] == hi[0]) faces[n++] = 1;
  if (p[1] == lo[1]) faces[n++] = 2;
  if (p[1] == hi[1]) faces[n++] = 3;
  if (p[2] == lo[2]) faces[n++] = 4;
  if (p[2] == hi[2]) faces[n++] = 5;
  return n;
}

/* ----------------------------------------------------------------------
   identify which cell faces edge between p1,p2 is on
   p1,p2 assumed to be on surface or interior of cell
   return list of face IDs (0-5)
   list length can be 0,1,2
------------------------------------------------------------------------- */

int Cut3d::which_faces(double *p1, double *p2, int *faces)
{
  int n = 0;
  if (p1[0] == lo[0] && p2[0] == lo[0]) faces[n++] = 0;
  if (p1[0] == hi[0] && p2[0] == hi[0]) faces[n++] = 1;
  if (p1[1] == lo[1] && p2[1] == lo[1]) faces[n++] = 2;
  if (p1[1] == hi[1] && p2[1] == hi[1]) faces[n++] = 3;
  if (p1[2] == lo[2] && p2[2] == lo[2]) faces[n++] = 4;
  if (p1[2] == hi[2] && p2[2] == hi[2]) faces[n++] = 5;
  return n;
}

/* ----------------------------------------------------------------------
# extract 2d cell from iface (0-5) of 3d cell
# return lo2d/hi2d = xlo,xhi,ylo,yhi
# for XLO/XHI, keep (y,z) -> (x,y), look at face from inside/outside 3d cell
# for YLO/YHI, keep (x,z) -> (x,y), look at face from outside/inside 3d cell
# for ZLO/ZHI, keep (x,y) -> (x,y), look at face from inside/outside 3d cell
------------------------------------------------------------------------- */

void Cut3d::face_from_cell(int iface, double *lo2d, double *hi2d)
{
  if (iface < 2) {
    lo2d[0] = lo[1]; hi2d[0] = hi[1];
    lo2d[1] = lo[2]; hi2d[1] = hi[2];
  } else if (iface < 4) {
    lo2d[0] = lo[0]; hi2d[0] = hi[0];
    lo2d[1] = lo[2]; hi2d[1] = hi[2];
  } else {
    lo2d[0] = lo[0]; hi2d[0] = hi[0];
    lo2d[1] = lo[1]; hi2d[1] = hi[1];
  }
}

/* ----------------------------------------------------------------------
   compress a 3d pt into a 2d pt on iface
------------------------------------------------------------------------- */

void Cut3d::compress2d(int iface, double *p3, double *p2)
{
  if (iface < 2) {
    p2[0] = p3[1]; p2[1] = p3[2];
  } else if (iface < 4) {
    p2[0] = p3[0]; p2[1] = p3[2];
  } else {
    p2[0] = p3[0]; p2[1] = p3[1];
  }
}

/* ----------------------------------------------------------------------
   expand a 2d pt into 3d pt on iface with extra coord = value
------------------------------------------------------------------------- */

void Cut3d::expand2d(int iface, double value, double *p2, double *p3)
{
  if (iface < 2) {
    p3[0] = value; p3[1] = p2[0]; p3[2] = p2[1];
  } else if (iface < 4) {
    p3[0] = p2[0]; p3[1] = value; p3[2] = p2[1];
  } else {
    p3[0] = p2[0]; p3[1] = p2[1]; p3[2] = value;
  }
}

/* ----------------------------------------------------------------------
   look for edge (x,y) in list of edges
   match as (x,y) or (y,x)
   if flag, do not match edges that are part of a CTRI (style = CTRI,CTRIFACE)
     this is used by add_face() when adding edges of an entire face
     this avoids matching an on-face CTRI with norm into cell
   error if find edge and it is already part of a vertex in that dir
   return = index (0 to nedge-1) if find it, -1 if do not find it
   return dir = 0 if matches as (x,y), 1 if matches as (y,x), -1 if no match
   return -2 as error if edge already exists in same dir as this one
------------------------------------------------------------------------- */

int Cut3d::findedge(double *x, double *y, int flag, int &dir)
{
  double *p1,*p2;

  int nedge = edges.n;

  for (int i = 0; i < nedge; i++) {
    if (!edges[i].active) continue;
    if (flag && (edges[i].style == CTRI || edges[i].style == CTRIFACE))
      continue;
    p1 = edges[i].p1;
    p2 = edges[i].p2;
    if (samepoint(x,p1) && samepoint(y,p2)) {
      if (edges[i].nvert % 2 == 1) return -2;
      dir = 0;
      return i;
    }
    if (samepoint(x,p2) && samepoint(y,p1)) {
      if (edges[i].nvert / 2 == 1) return -2;
      dir = 1;
      return i;
    }
  }

  dir = -1;
  return -1;
}

/* ----------------------------------------------------------------------
   return intersection pt C of line segment A,B in dim with coord value
   guaranteed to intersect by caller
   C can be same as A or B, will just overwrite
------------------------------------------------------------------------- */

void Cut3d::between(double *a, double *b, int dim, double value, double *c)
{
  if (dim == 0) {
    c[1] = a[1] + (value-a[dim])/(b[dim]-a[dim]) * (b[1]-a[1]);
    c[2] = a[2] + (value-a[dim])/(b[dim]-a[dim]) * (b[2]-a[2]);
    c[0] = value;
  } else if (dim == 1) {
    c[0] = a[0] + (value-a[dim])/(b[dim]-a[dim]) * (b[0]-a[0]);
    c[2] = a[2] + (value-a[dim])/(b[dim]-a[dim]) * (b[2]-a[2]);
    c[1] = value;
  } else {
    c[0] = a[0] + (value-a[dim])/(b[dim]-a[dim]) * (b[0]-a[0]);
    c[1] = a[1] + (value-a[dim])/(b[dim]-a[dim]) * (b[1]-a[1]);
    c[2] = value;
  }
}

/* ----------------------------------------------------------------------
   return 1 if x,y are same point, else 0
------------------------------------------------------------------------- */

int Cut3d::samepoint(double *x, double *y)
{
  if (x[0] == y[0] && x[1] == y[1] && x[2] == y[2]) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   return 0-7 if pt is a corner pt of grid cell
   else return -1
------------------------------------------------------------------------- */

int Cut3d::corner(double *pt)
{
  if (pt[2] == lo[2]) {
    if (pt[1] == lo[1]) {
      if (pt[0] == lo[0]) return 0;
      else if (pt[0] == hi[0]) return 1;
    } else if (pt[1] == hi[1]) {
      if (pt[0] == lo[0]) return 2;
      else if (pt[0] == hi[0]) return 3;
    }
  } else if (pt[2] == hi[2]) {
    if (pt[1] == lo[1]) {
      if (pt[0] == lo[0]) return 4;
      else if (pt[0] == hi[0]) return 5;
    } else if (pt[1] == hi[1]) {
      if (pt[0] == lo[0]) return 6;
      else if (pt[0] == hi[0]) return 7;
    }
  }

  return -1;
}

/* ----------------------------------------------------------------------
   move point within epsilon of any cell face to be on cell faces
------------------------------------------------------------------------- */

void Cut3d::move_to_faces(double *pt)
{
  if (fabs(pt[0]-lo[0]) < epsilon) pt[0] = lo[0];
  if (fabs(pt[0]-hi[0]) < epsilon) pt[0] = hi[0];
  if (fabs(pt[1]-lo[1]) < epsilon) pt[1] = lo[1];
  if (fabs(pt[1]-hi[1]) < epsilon) pt[1] = hi[1];
  if (fabs(pt[2]-lo[2]) < epsilon) pt[2] = lo[2];
  if (fabs(pt[2]-hi[2]) < epsilon) pt[2] = hi[2];
}

/* ----------------------------------------------------------------------
   check if pt is inside or outside or on cell border
   return EXTERIOR,BORDER,INTERIOR
------------------------------------------------------------------------- */

int Cut3d::ptflag(double *pt)
{
  double x = pt[0];
  double y = pt[1];
  double z = pt[2];
  if (x < lo[0] || x > hi[0] || y < lo[1] || y > hi[1] ||
      z < lo[2] || z > hi[2]) return EXTERIOR;
  if (x > lo[0] && x < hi[0] && y > lo[1] && y < hi[1] &&
      z > lo[2] && z < hi[2]) return INTERIOR;
  return BORDER;
}

/* ----------------------------------------------------------------------
   print out cell info for cell which failed at cut/split operation
------------------------------------------------------------------------- */

void Cut3d::failed_cell()
{
  printf("Cut3d failed on proc %d in cell ID: " CELLINT_FORMAT "\n",comm->me,id);
  Surf::Tri *tris = surf->tris;
  printf("  lo corner %g %g %g\n",lo[0],lo[1],lo[2]);
  printf("  hi corner %g %g %g\n",hi[0],hi[1],hi[2]);
  printf("  # of surfs = %d out of " BIGINT_FORMAT "\n",nsurf,surf->nsurf);
  for (int i = 0; i < nsurf; i++) {
    printf("  surf " SURFINT_FORMAT ":\n",tris[surfs[i]].id);
    printf("     p1: %g %g %g\n",
           tris[surfs[i]].p1[0],tris[surfs[i]].p1[1],tris[surfs[i]].p1[2]);
    printf("     p2: %g %g %g\n",
           tris[surfs[i]].p2[0],tris[surfs[i]].p2[1],tris[surfs[i]].p2[2]);
    printf("     p3: %g %g %g\n",
           tris[surfs[i]].p3[0],tris[surfs[i]].p3[1],tris[surfs[i]].p3[2]);
  }
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void Cut3d::print_bpg(const char *str)
{
  int iedge,dir,newedge,newdir,prevedge,prevdir;
  double *p1,*p2;

  printf("%s " CELLINT_FORMAT "\n",str,id);
  printf("  Sizes: %d %d\n",verts.n,edges.n);

  printf("  Verts:\n");
  for (int i = 0; i < verts.n; i++) {
    if (verts[i].active == 0) continue;

    printf("   %d %d %d %d:",i,
           verts[i].active,verts[i].style,verts[i].label);

    printf(" edges [");
    iedge = verts[i].first;
    dir = verts[i].dirfirst;
    for (int j = 0; j < verts[i].nedge; j++) {
      printf("%d",iedge);
      if (j < verts[i].nedge-1) printf(" ");
      newedge = edges[iedge].next[dir];
      newdir = edges[iedge].dirnext[dir];
      iedge = newedge;
      dir = newdir;
    }
    printf("]");

    printf(" dirs [");
    iedge = verts[i].first;
    dir = verts[i].dirfirst;
    for (int j = 0; j < verts[i].nedge; j++) {
      printf("%d",dir);
      if (j < verts[i].nedge-1) printf(" ");
      newedge = edges[iedge].next[dir];
      newdir = edges[iedge].dirnext[dir];
      iedge = newedge;
      dir = newdir;
    }
    printf("]");

    printf(" edgelens [");
    iedge = verts[i].first;
    dir = verts[i].dirfirst;
    for (int j = 0; j < verts[i].nedge; j++) {
      p1 = edges[iedge].p1;
      p2 = edges[iedge].p2;
      double dx = p1[0]-p2[0];
      double dy = p1[1]-p2[1];
      double dz = p1[2]-p2[2];
      double delta = sqrt(dx*dx+dy*dy+dz*dz);
      printf("%g",delta);
      if (j < verts[i].nedge-1) printf(" ");
      newedge = edges[iedge].next[dir];
      newdir = edges[iedge].dirnext[dir];
      iedge = newedge;
      dir = newdir;
    }
    printf("]");

    printf(" samept [");
    iedge = verts[i].first;
    dir = verts[i].dirfirst;
    for (int j = 0; j < verts[i].nedge; j++) {
      if (j > 0) {
        prevedge = edges[iedge].prev[dir];
        prevdir = edges[iedge].dirprev[dir];
      } else {
        prevedge = verts[i].last;
        prevdir = verts[i].dirlast;
      }
      if (dir == 0) p1 = edges[iedge].p1;
      else p1 = edges[iedge].p2;
      if (prevdir == 0) p2 = edges[prevedge].p2;
      else p2 = edges[prevedge].p1;
      printf("%d",samepoint(p1,p2));
      if (j < verts[i].nedge-1) printf(" ");
      newedge = edges[iedge].next[dir];
      newdir = edges[iedge].dirnext[dir];
      iedge = newedge;
      dir = newdir;
    }
    printf("]");

    if (verts[i].norm) {
      printf(" norm [%g %g %g]\n",
             verts[i].norm[0],verts[i].norm[1],verts[i].norm[2]);
    } else printf(" [NULL]\n");
  }

  printf("  Edges:\n");
  for (int i = 0; i < edges.n; i++) {
    if (edges[i].active == 0) continue;

    printf("   %d %d %d",i,edges[i].active,edges[i].style);
    printf(" (%g %g %g)",edges[i].p1[0],edges[i].p1[1],edges[i].p1[2]);
    printf(" (%g %g %g)",edges[i].p2[0],edges[i].p2[1],edges[i].p2[2]);
    if (edges[i].nvert == 0) printf(" [-1]");
    if (edges[i].nvert == 1) {
      printf(" [%d]",edges[i].verts[0]);
      printf(" p1: [%d %d]",edges[i].prev[0],edges[i].dirprev[0]);
      printf(" n1: [%d %d]",edges[i].next[0],edges[i].dirnext[0]);
    }
    if (edges[i].nvert == 2) {
      printf(" [%d]",edges[i].verts[1]);
      printf(" p1: [%d %d]",edges[i].prev[1],edges[i].dirprev[1]);
      printf(" n1: [%d %d]",edges[i].next[1],edges[i].dirnext[1]);
    }
    if (edges[i].nvert == 3) {
      printf(" [%d %d]",edges[i].verts[0],edges[i].verts[1]);
      printf(" p1: [%d %d]",edges[i].prev[0],edges[i].dirprev[0]);
      printf(" n1: [%d %d]",edges[i].next[0],edges[i].dirnext[0]);
      printf(" p2: [%d %d]",edges[i].prev[1],edges[i].dirprev[1]);
      printf(" n2: [%d %d]",edges[i].next[1],edges[i].dirnext[1]);
    }
    if (edges[i].nvert > 3) printf(" [BIG %d]",edges[i].nvert);
    printf("\n");
  }
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void Cut3d::print_loops()
{
  printf("LOOP " CELLINT_FORMAT "\n",id);
  printf("  loops %d\n",loops.n);
  for (int i = 0; i < loops.n; i++) {
    printf("  loop %d\n",i);
    printf("    flag %d\n",loops[i].flag);
    printf("    volume %g\n",loops[i].volume);
    printf("    nverts %d\n",loops[i].n);
    printf("    verts: [");
    int ivert = loops[i].first;
    for (int j = 0; j < loops[i].n; j++) {
      printf("%d ",ivert);
      ivert = verts[ivert].next;
    }
    printf("]\n");
  }
}

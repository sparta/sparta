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

#include "geometry.h"
#include "math_extra.h"

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

#define EPS 1.0e-8
#define EPSSQ 1.0e-16
#define EPSSQNEG -1.0e-16

enum{OUTSIDE,INSIDE,ONSURF2OUT,ONSURF2IN};    // same as Update

namespace Geometry {

/* ----------------------------------------------------------------------
   compute whether line intersects an orthogonal 2d quad cell
   intersection is defined as
     any line pt (interior, vertex) in common with
     any quad pt (interior, edge, vertex)
   v0,v1 and norm = 2 vertices of line and unit normal vec
   lo,hi = opposite corner pts of quad
   return 1 if intersection, else 0
------------------------------------------------------------------------- */

int line_quad_intersect(double *v0, double *v1, double *norm,
			double *lo, double *hi)
{
  int sum,side;
  double xlo,xhi,ylo,yhi,param;
  double b[3],e[3],point[3];

  xlo = lo[0];
  xhi = hi[0];
  ylo = lo[1];
  yhi = hi[1];

  // if either of line vertices are inside quad, intersection
  // use <= and >= so touching quad surface is same as inside it
  // important to do this test first, b/c whichside test can be epsilon off

  if (v0[0] >= xlo && v0[0] <= xhi && v0[1] >= ylo && v0[1] <= yhi) return 1;
  if (v1[0] >= xlo && v1[0] <= xhi && v1[1] >= ylo && v1[1] <= yhi) return 1;

  // if all 4 quad pts are on same side of line, no intersection

  sum = whichside(v0,norm,xlo,ylo,0.0);
  sum += whichside(v0,norm,xhi,ylo,0.0);
  sum += whichside(v0,norm,xlo,yhi,0.0);
  sum += whichside(v0,norm,xhi,yhi,0.0);
  
  if (sum == 4 || sum == -4) return 0;
	
  // test 4 quad edges for intersection with line
  // b,e = begin/end of quad edge line segment

  b[0] = xlo;   b[1] = ylo;   b[2] = 0.0;
  e[0] = xhi;   e[1] = ylo;   e[2] = 0.0;
  if (line_line_intersect(b,e,v0,v1,norm,point,param,side)) return 1;

  b[0] = xhi;   b[1] = ylo;   b[2] = 0.0;
  e[0] = xhi;   e[1] = yhi;   e[2] = 0.0;
  if (line_line_intersect(b,e,v0,v1,norm,point,param,side)) return 1;

  b[0] = xhi;   b[1] = yhi;   b[2] = 0.0;
  e[0] = xlo;   e[1] = yhi;   e[2] = 0.0;
  if (line_line_intersect(b,e,v0,v1,norm,point,param,side)) return 1;

  b[0] = xlo;   b[1] = yhi;   b[2] = 0.0;
  e[0] = xlo;   e[1] = ylo;   e[2] = 0.0;
  if (line_line_intersect(b,e,v0,v1,norm,point,param,side)) return 1;

  return 0;
}

/* ----------------------------------------------------------------------
   compute any intersection of edges of orthogonal 2d quad cell with a line
   line interior to quad cell has no intersection
   v0,v1 and norm = 2 vertices of line and unit normal vec
   lo,hi = opposite corner pts of quad
   return 1 if intersection, else 0
   return xc = intersection point if there is one
------------------------------------------------------------------------- */

int quad_line_intersect_point(double *v0, double *v1, double *norm,
			      double *lo, double *hi, double *xc)
{
  int side;
  double xlo,xhi,ylo,yhi,param;
  double b[3],e[3];

  xlo = lo[0];
  xhi = hi[0];
  ylo = lo[1];
  yhi = hi[1];

  // test 4 quad edges for intersection with line
  // b,e = begin/end of quad edge line segment

  b[0] = xlo;   b[1] = ylo;   b[2] = 0.0;
  e[0] = xhi;   e[1] = ylo;   e[2] = 0.0;
  if (line_line_intersect(b,e,v0,v1,norm,xc,param,side)) return 1;

  b[0] = xhi;   b[1] = ylo;   b[2] = 0.0;
  e[0] = xhi;   e[1] = yhi;   e[2] = 0.0;
  if (line_line_intersect(b,e,v0,v1,norm,xc,param,side)) return 1;

  b[0] = xhi;   b[1] = yhi;   b[2] = 0.0;
  e[0] = xlo;   e[1] = yhi;   e[2] = 0.0;
  if (line_line_intersect(b,e,v0,v1,norm,xc,param,side)) return 1;

  b[0] = xlo;   b[1] = yhi;   b[2] = 0.0;
  e[0] = xlo;   e[1] = ylo;   e[2] = 0.0;
  if (line_line_intersect(b,e,v0,v1,norm,xc,param,side)) return 1;

  return 0;
}

/* ----------------------------------------------------------------------
   compute whether line touches iface of orthogonal 2d quad cell
   touch is defined as
     line end pt = any face pt (edge, vertex)
   v0,v1 = 2 vertices of line
   iface = 0 to 3 = XLO,XHI,YLO,YHI
   lo,hi = opposite corner pts of quad
   return 1 if touches, else 0
------------------------------------------------------------------------- */

int line_quad_face_touch(double *v0, double *v1, int iface,
			 double *lo, double *hi)
{
  // value = position of face

  int dim = iface / 2;
  int other = dim ? 0 : 1;
  double value = iface % 2 ? hi[dim] : lo[dim];

  // check if either line vertex is within face

  if (v0[dim] == value) {
    if (v0[other] >= lo[other] && v0[other] <= hi[other]) return 1;
  }
  if (v1[dim] == value) {
    if (v1[other] >= lo[other] && v1[other] <= hi[other]) return 1;
  }

  return 0;
}

/* ----------------------------------------------------------------------
   compute whether triangle intersects an orthogonal 3d hex cell
   intersection is defined as
     any triangle pt (interior, edge, vertex) in common with
     any hex pt (interior, face, edge, vertex)
   v0,v1,v2 and norm = 3 vertices of triangle and unit normal vec
   lo,hi = opposite corner pts of hex
   return 1 if intersection, else 0
------------------------------------------------------------------------- */

int tri_hex_intersect(double *v0, double *v1, double *v2, double *norm,
		      double *lo, double *hi)
{
  int sum,side;
  double xlo,xhi,ylo,yhi,zlo,zhi,param;
  double b[3],e[3],h0[3],h1[3],h2[3],h3[3],n[3],point[3];

  xlo = lo[0];
  xhi = hi[0];
  ylo = lo[1];
  yhi = hi[1];
  zlo = lo[2];
  zhi = hi[2];

  // if any of 3 tri vertices are inside hex, intersection
  // use <= and >= so touching hex surface is same as inside it
  // important to do this test first, b/c whichside test can be epsilon off

  if (v0[0] >= xlo && v0[0] <= xhi && v0[1] >= ylo && v0[1] <= yhi &&
      v0[2] >= zlo && v0[2] <= zhi) return 1;

  if (v1[0] >= xlo && v1[0] <= xhi && v1[1] >= ylo && v1[1] <= yhi &&
      v1[2] >= zlo && v1[2] <= zhi) return 1;

  if (v2[0] >= xlo && v2[0] <= xhi && v2[1] >= ylo && v2[1] <= yhi &&
      v2[2] >= zlo && v2[2] <= zhi) return 1;
  
  // if all 8 hex pts are on same side of tri plane, no intersection

  sum = whichside(v0,norm,xlo,ylo,zlo);
  sum += whichside(v0,norm,xhi,ylo,zlo);
  sum += whichside(v0,norm,xlo,yhi,zlo);
  sum += whichside(v0,norm,xhi,yhi,zlo);
  sum += whichside(v0,norm,xlo,ylo,zhi);
  sum += whichside(v0,norm,xhi,ylo,zhi);
  sum += whichside(v0,norm,xlo,yhi,zhi);
  sum += whichside(v0,norm,xhi,yhi,zhi);
  if (sum == 8 || sum == -8) return 0;
  
  // test 12 hex edges for intersection with tri
  // b,e = begin/end of hex edge line segment

  b[0] = xlo;   b[1] = ylo;   b[2] = zlo;
  e[0] = xhi;   e[1] = ylo;   e[2] = zlo;
  if (line_tri_intersect(b,e,v0,v1,v2,norm,point,param,side)) return 1;

  b[0] = xlo;   b[1] = yhi;   b[2] = zlo;
  e[0] = xhi;   e[1] = yhi;   e[2] = zlo;
  if (line_tri_intersect(b,e,v0,v1,v2,norm,point,param,side)) return 1;

  b[0] = xlo;   b[1] = ylo;   b[2] = zhi;
  e[0] = xhi;   e[1] = ylo;   e[2] = zhi;
  if (line_tri_intersect(b,e,v0,v1,v2,norm,point,param,side)) return 1;

  b[0] = xlo;   b[1] = yhi;   b[2] = zhi;
  e[0] = xhi;   e[1] = yhi;   e[2] = zhi;
  if (line_tri_intersect(b,e,v0,v1,v2,norm,point,param,side)) return 1;

  b[0] = xlo;   b[1] = ylo;   b[2] = zlo;
  e[0] = xlo;   e[1] = yhi;   e[2] = zlo;
  if (line_tri_intersect(b,e,v0,v1,v2,norm,point,param,side)) return 1;

  b[0] = xhi;   b[1] = ylo;   b[2] = zlo;
  e[0] = xhi;   e[1] = yhi;   e[2] = zlo;
  if (line_tri_intersect(b,e,v0,v1,v2,norm,point,param,side)) return 1;

  b[0] = xlo;   b[1] = ylo;   b[2] = zhi;
  e[0] = xlo;   e[1] = yhi;   e[2] = zhi;
  if (line_tri_intersect(b,e,v0,v1,v2,norm,point,param,side)) return 1;

  b[0] = xhi;   b[1] = ylo;   b[2] = zhi;
  e[0] = xhi;   e[1] = yhi;   e[2] = zhi;
  if (line_tri_intersect(b,e,v0,v1,v2,norm,point,param,side)) return 1;

  b[0] = xlo;   b[1] = ylo;   b[2] = zlo;
  e[0] = xlo;   e[1] = ylo;   e[2] = zhi;
  if (line_tri_intersect(b,e,v0,v1,v2,norm,point,param,side)) return 1;

  b[0] = xhi;   b[1] = ylo;   b[2] = zlo;
  e[0] = xhi;   e[1] = ylo;   e[2] = zhi;
  if (line_tri_intersect(b,e,v0,v1,v2,norm,point,param,side)) return 1;

  b[0] = xlo;   b[1] = yhi;   b[2] = zlo;
  e[0] = xlo;   e[1] = yhi;   e[2] = zhi;
  if (line_tri_intersect(b,e,v0,v1,v2,norm,point,param,side)) return 1;

  b[0] = xhi;   b[1] = yhi;   b[2] = zlo;
  e[0] = xhi;   e[1] = yhi;   e[2] = zhi;
  if (line_tri_intersect(b,e,v0,v1,v2,norm,point,param,side)) return 1;

  // test 3 tri edges for intersection with 6 faces of hex
  // h0,h1,h2,h3 = 4 corner pts of hex face
  // n = normal to xyz faces, depends on vertex ordering
  // each face is treated as 2 triangles -> 6 tests per face
  
  h0[0] = xlo;  h0[1] = ylo;  h0[2] = zlo;
  h1[0] = xlo;  h1[1] = yhi;  h1[2] = zlo;
  h2[0] = xlo;  h2[1] = yhi;  h2[2] = zhi;
  h3[0] = xlo;  h3[1] = ylo;  h3[2] = zhi;
  n[0]  = 1.0;  n[1]  = 0.0;  n[2]  = 0.0;
  
  if (line_tri_intersect(v0,v1,h0,h1,h2,n,point,param,side) ||
      line_tri_intersect(v1,v2,h0,h1,h2,n,point,param,side) ||
      line_tri_intersect(v2,v0,h0,h1,h2,n,point,param,side) ||
      line_tri_intersect(v0,v1,h0,h2,h3,n,point,param,side) ||
      line_tri_intersect(v1,v2,h0,h2,h3,n,point,param,side) ||
      line_tri_intersect(v2,v0,h0,h2,h3,n,point,param,side)) return 1;

  h0[0] = h1[0] = h2[0] = h3[0] = xhi;

  if (line_tri_intersect(v0,v1,h0,h1,h2,n,point,param,side) ||
      line_tri_intersect(v1,v2,h0,h1,h2,n,point,param,side) ||
      line_tri_intersect(v2,v0,h0,h1,h2,n,point,param,side) ||
      line_tri_intersect(v0,v1,h0,h2,h3,n,point,param,side) ||
      line_tri_intersect(v1,v2,h0,h2,h3,n,point,param,side) ||
      line_tri_intersect(v2,v0,h0,h2,h3,n,point,param,side)) return 1;

  h0[0] = xlo;  h0[1] = ylo;  h0[2] = zlo;
  h1[0] = xhi;  h1[1] = ylo;  h1[2] = zlo;
  h2[0] = xhi;  h2[1] = ylo;  h2[2] = zhi;
  h3[0] = xlo;  h3[1] = ylo;  h3[2] = zhi;
  n[0]  = 0.0;  n[1]  = -1.0;  n[2]  = 0.0;
  
  if (line_tri_intersect(v0,v1,h0,h1,h2,n,point,param,side) ||
      line_tri_intersect(v1,v2,h0,h1,h2,n,point,param,side) ||
      line_tri_intersect(v2,v0,h0,h1,h2,n,point,param,side) ||
      line_tri_intersect(v0,v1,h0,h2,h3,n,point,param,side) ||
      line_tri_intersect(v1,v2,h0,h2,h3,n,point,param,side) ||
      line_tri_intersect(v2,v0,h0,h2,h3,n,point,param,side)) return 1;

  h0[1] = h1[1] = h2[1] = h3[1] = yhi;

  if (line_tri_intersect(v0,v1,h0,h1,h2,n,point,param,side) ||
      line_tri_intersect(v1,v2,h0,h1,h2,n,point,param,side) ||
      line_tri_intersect(v2,v0,h0,h1,h2,n,point,param,side) ||
      line_tri_intersect(v0,v1,h0,h2,h3,n,point,param,side) ||
      line_tri_intersect(v1,v2,h0,h2,h3,n,point,param,side) ||
      line_tri_intersect(v2,v0,h0,h2,h3,n,point,param,side)) return 1;

  h0[0] = xlo;  h0[1] = ylo;  h0[2] = zlo;
  h1[0] = xhi;  h1[1] = ylo;  h1[2] = zlo;
  h2[0] = xhi;  h2[1] = yhi;  h2[2] = zlo;
  h3[0] = xlo;  h3[1] = yhi;  h3[2] = zlo;
  n[0]  = 0.0;  n[1]  = 0.0;  n[2]  = 1.0;
  
  if (line_tri_intersect(v0,v1,h0,h1,h2,n,point,param,side) ||
      line_tri_intersect(v1,v2,h0,h1,h2,n,point,param,side) ||
      line_tri_intersect(v2,v0,h0,h1,h2,n,point,param,side) ||
      line_tri_intersect(v0,v1,h0,h2,h3,n,point,param,side) ||
      line_tri_intersect(v1,v2,h0,h2,h3,n,point,param,side) ||
      line_tri_intersect(v2,v0,h0,h2,h3,n,point,param,side)) return 1;

  h0[2] = h1[2] = h2[2] = h3[2] = zhi;
  
  if (line_tri_intersect(v0,v1,h0,h1,h2,n,point,param,side) ||
      line_tri_intersect(v1,v2,h0,h1,h2,n,point,param,side) ||
      line_tri_intersect(v2,v0,h0,h1,h2,n,point,param,side) ||
      line_tri_intersect(v0,v1,h0,h2,h3,n,point,param,side) ||
      line_tri_intersect(v1,v2,h0,h2,h3,n,point,param,side) ||
      line_tri_intersect(v2,v0,h0,h2,h3,n,point,param,side)) return 1;

  return 0;
}

/* ----------------------------------------------------------------------
   compute any intersection of edges/faces of orthogonal 3d hex cell with a tri
   tri interior to quad cell has no intersection
   v0,v1,v2 and norm = 3 vertices of triangle and unit normal vec
   lo,hi = opposite corner pts of hex
   return 1 if intersection, else 0
   return xc = intersection point if there is one
------------------------------------------------------------------------- */

int hex_tri_intersect_point(double *v0, double *v1, double *v2, double *norm,
			    double *lo, double *hi, double *xc)
{
  int side;
  double xlo,xhi,ylo,yhi,zlo,zhi,param;
  double b[3],e[3],h0[3],h1[3],h2[3],h3[3],n[3];

  xlo = lo[0];
  xhi = hi[0];
  ylo = lo[1];
  yhi = hi[1];
  zlo = lo[2];
  zhi = hi[2];

  // test 12 hex edges for intersection with tri
  // b,e = begin/end of hex edge line segment

  b[0] = xlo;   b[1] = ylo;   b[2] = zlo;
  e[0] = xhi;   e[1] = ylo;   e[2] = zlo;
  if (line_tri_intersect(b,e,v0,v1,v2,norm,xc,param,side)) return 1;

  b[0] = xlo;   b[1] = yhi;   b[2] = zlo;
  e[0] = xhi;   e[1] = yhi;   e[2] = zlo;
  if (line_tri_intersect(b,e,v0,v1,v2,norm,xc,param,side)) return 1;

  b[0] = xlo;   b[1] = ylo;   b[2] = zhi;
  e[0] = xhi;   e[1] = ylo;   e[2] = zhi;
  if (line_tri_intersect(b,e,v0,v1,v2,norm,xc,param,side)) return 1;

  b[0] = xlo;   b[1] = yhi;   b[2] = zhi;
  e[0] = xhi;   e[1] = yhi;   e[2] = zhi;
  if (line_tri_intersect(b,e,v0,v1,v2,norm,xc,param,side)) return 1;

  b[0] = xlo;   b[1] = ylo;   b[2] = zlo;
  e[0] = xlo;   e[1] = yhi;   e[2] = zlo;
  if (line_tri_intersect(b,e,v0,v1,v2,norm,xc,param,side)) return 1;

  b[0] = xhi;   b[1] = ylo;   b[2] = zlo;
  e[0] = xhi;   e[1] = yhi;   e[2] = zlo;
  if (line_tri_intersect(b,e,v0,v1,v2,norm,xc,param,side)) return 1;

  b[0] = xlo;   b[1] = ylo;   b[2] = zhi;
  e[0] = xlo;   e[1] = yhi;   e[2] = zhi;
  if (line_tri_intersect(b,e,v0,v1,v2,norm,xc,param,side)) return 1;

  b[0] = xhi;   b[1] = ylo;   b[2] = zhi;
  e[0] = xhi;   e[1] = yhi;   e[2] = zhi;
  if (line_tri_intersect(b,e,v0,v1,v2,norm,xc,param,side)) return 1;

  b[0] = xlo;   b[1] = ylo;   b[2] = zlo;
  e[0] = xlo;   e[1] = ylo;   e[2] = zhi;
  if (line_tri_intersect(b,e,v0,v1,v2,norm,xc,param,side)) return 1;

  b[0] = xhi;   b[1] = ylo;   b[2] = zlo;
  e[0] = xhi;   e[1] = ylo;   e[2] = zhi;
  if (line_tri_intersect(b,e,v0,v1,v2,norm,xc,param,side)) return 1;

  b[0] = xlo;   b[1] = yhi;   b[2] = zlo;
  e[0] = xlo;   e[1] = yhi;   e[2] = zhi;
  if (line_tri_intersect(b,e,v0,v1,v2,norm,xc,param,side)) return 1;

  b[0] = xhi;   b[1] = yhi;   b[2] = zlo;
  e[0] = xhi;   e[1] = yhi;   e[2] = zhi;
  if (line_tri_intersect(b,e,v0,v1,v2,norm,xc,param,side)) return 1;

  // test 3 tri edges for intersection with 6 faces of hex
  // h0,h1,h2,h3 = 4 corner pts of hex face
  // n = normal to xyz faces, depends on vertex ordering
  // each face is treated as 2 triangles -> 6 tests per face
  
  h0[0] = xlo;  h0[1] = ylo;  h0[2] = zlo;
  h1[0] = xlo;  h1[1] = yhi;  h1[2] = zlo;
  h2[0] = xlo;  h2[1] = yhi;  h2[2] = zhi;
  h3[0] = xlo;  h3[1] = ylo;  h3[2] = zhi;
  n[0]  = 1.0;  n[1]  = 0.0;  n[2]  = 0.0;
  
  if (line_tri_intersect(v0,v1,h0,h1,h2,n,xc,param,side) ||
      line_tri_intersect(v1,v2,h0,h1,h2,n,xc,param,side) ||
      line_tri_intersect(v2,v0,h0,h1,h2,n,xc,param,side) ||
      line_tri_intersect(v0,v1,h0,h2,h3,n,xc,param,side) ||
      line_tri_intersect(v1,v2,h0,h2,h3,n,xc,param,side) ||
      line_tri_intersect(v2,v0,h0,h2,h3,n,xc,param,side)) return 1;

  h0[0] = h1[0] = h2[0] = h3[0] = xhi;

  if (line_tri_intersect(v0,v1,h0,h1,h2,n,xc,param,side) ||
      line_tri_intersect(v1,v2,h0,h1,h2,n,xc,param,side) ||
      line_tri_intersect(v2,v0,h0,h1,h2,n,xc,param,side) ||
      line_tri_intersect(v0,v1,h0,h2,h3,n,xc,param,side) ||
      line_tri_intersect(v1,v2,h0,h2,h3,n,xc,param,side) ||
      line_tri_intersect(v2,v0,h0,h2,h3,n,xc,param,side)) return 1;

  h0[0] = xlo;  h0[1] = ylo;  h0[2] = zlo;
  h1[0] = xhi;  h1[1] = ylo;  h1[2] = zlo;
  h2[0] = xhi;  h2[1] = ylo;  h2[2] = zhi;
  h3[0] = xlo;  h3[1] = ylo;  h3[2] = zhi;
  n[0]  = 0.0;  n[1]  = -1.0;  n[2]  = 0.0;
  
  if (line_tri_intersect(v0,v1,h0,h1,h2,n,xc,param,side) ||
      line_tri_intersect(v1,v2,h0,h1,h2,n,xc,param,side) ||
      line_tri_intersect(v2,v0,h0,h1,h2,n,xc,param,side) ||
      line_tri_intersect(v0,v1,h0,h2,h3,n,xc,param,side) ||
      line_tri_intersect(v1,v2,h0,h2,h3,n,xc,param,side) ||
      line_tri_intersect(v2,v0,h0,h2,h3,n,xc,param,side)) return 1;

  h0[1] = h1[1] = h2[1] = h3[1] = yhi;

  if (line_tri_intersect(v0,v1,h0,h1,h2,n,xc,param,side) ||
      line_tri_intersect(v1,v2,h0,h1,h2,n,xc,param,side) ||
      line_tri_intersect(v2,v0,h0,h1,h2,n,xc,param,side) ||
      line_tri_intersect(v0,v1,h0,h2,h3,n,xc,param,side) ||
      line_tri_intersect(v1,v2,h0,h2,h3,n,xc,param,side) ||
      line_tri_intersect(v2,v0,h0,h2,h3,n,xc,param,side)) return 1;

  h0[0] = xlo;  h0[1] = ylo;  h0[2] = zlo;
  h1[0] = xhi;  h1[1] = ylo;  h1[2] = zlo;
  h2[0] = xhi;  h2[1] = yhi;  h2[2] = zlo;
  h3[0] = xlo;  h3[1] = yhi;  h3[2] = zlo;
  n[0]  = 0.0;  n[1]  = 0.0;  n[2]  = 1.0;
  
  if (line_tri_intersect(v0,v1,h0,h1,h2,n,xc,param,side) ||
      line_tri_intersect(v1,v2,h0,h1,h2,n,xc,param,side) ||
      line_tri_intersect(v2,v0,h0,h1,h2,n,xc,param,side) ||
      line_tri_intersect(v0,v1,h0,h2,h3,n,xc,param,side) ||
      line_tri_intersect(v1,v2,h0,h2,h3,n,xc,param,side) ||
      line_tri_intersect(v2,v0,h0,h2,h3,n,xc,param,side)) return 1;

  h0[2] = h1[2] = h2[2] = h3[2] = zhi;
  
  if (line_tri_intersect(v0,v1,h0,h1,h2,n,xc,param,side) ||
      line_tri_intersect(v1,v2,h0,h1,h2,n,xc,param,side) ||
      line_tri_intersect(v2,v0,h0,h1,h2,n,xc,param,side) ||
      line_tri_intersect(v0,v1,h0,h2,h3,n,xc,param,side) ||
      line_tri_intersect(v1,v2,h0,h2,h3,n,xc,param,side) ||
      line_tri_intersect(v2,v0,h0,h2,h3,n,xc,param,side)) return 1;

  return 0;
}

/* ----------------------------------------------------------------------
   compute whether triangle touches iface of orthogonal 3d hex cell
   touch is defined as
     triangle corner pt = any face pt (interior, edge, vertex)
   v0,v1,v2 = 3 vertices of triangle
   iface = 0 to 5 = XLO,XHI,YLO,YHI,ZLO,ZHI
   lo,hi = opposite corner pts of quad
   return 1 if touches, else 0
------------------------------------------------------------------------- */

int tri_hex_face_touch(double *v0, double *v1, double *v2, int iface,
		       double *lo, double *hi)
{
  // value = position of face

  int dim = iface / 2;
  int other1,other2;
  if (dim == 0) {
    other1 = 1; other2 = 2;
  } else if (dim == 1) {
    other1 = 0; other2 = 2;
  } else if (dim == 2) {
    other1 = 0; other2 = 1;
  }
  double value = iface % 2 ? hi[dim] : lo[dim];

  // check if any triangle vertex is within face

  if (v0[dim] == value) {
    if (v0[other1] >= lo[other1] && v0[other1] <= hi[other1] &&
	v0[other2] >= lo[other2] && v0[other2] <= hi[other2]) return 1;
  }
  if (v1[dim] == value) {
    if (v1[other1] >= lo[other1] && v1[other1] <= hi[other1] &&
	v1[other2] >= lo[other2] && v1[other2] <= hi[other2]) return 1;
  }
  if (v2[dim] == value) {
    if (v2[other1] >= lo[other1] && v2[other1] <= hi[other1] &&
	v2[other2] >= lo[other2] && v2[other2] <= hi[other2]) return 1;
  }

  return 0;
}

/* ----------------------------------------------------------------------
   detect intersection between a directed line segment A and line segment B
   intersection is defined as any A pt (including end pts)
     in common with any B pt (interior,vertex)
   one exception is if both A end pts are on infinite line B,
     then is NOT an intersection
   start,stop = end points of directed line segment A
     A must have non-zero length
   v0,v1 = 2 vertices of line segment B
   norm = unit vector normal to line segment B
     pointing OUTSIDE via right-hand rule
   return TRUE if there is an intersection, else FALSE
   if TRUE also return:
     point = pt of intersection
     param = intersection pt is this fraction along line A (0-1 inclusive)
     side = side of B that was hit = OUTSIDE,INSIDE,ONSURF2OUT,ONSURF2IN
------------------------------------------------------------------------- */

bool line_line_intersect(double *start, double *stop,
			 double *v0, double *v1, double *norm, 
			 double *point, double &param, int &side, int id)
{
  double vec[3],start2stop[3],edge[3],pvec[3];

  // if start,stop are on same side of line B, no intersection
  // if start,stop are both on infinite line B, no intersection

  MathExtra::sub3(start,v0,vec);
  double dotstart = MathExtra::dot3(norm,vec);
  MathExtra::sub3(stop,v0,vec);
  double dotstop = MathExtra::dot3(norm,vec);

  if (dotstart < 0.0 && dotstop < 0.0) return false;
  if (dotstart > 0.0 && dotstop > 0.0) return false;
  if (dotstart == 0.0 && dotstop == 0.0) return false;

  // param = parametric distance from start to stop
  //   at which line B is intersected
  // param = must be 0.0 to 1.0 inclusive, else no intersection

  MathExtra::sub3(v0,start,vec);
  MathExtra::sub3(stop,start,start2stop);
  param = MathExtra::dot3(norm,vec) / MathExtra::dot3(norm,start2stop);

  if (param < 0.0 || param > 1.0) return false;

  // point = intersection pt with line B

  point[0] = start[0] + param * start2stop[0];
  point[1] = start[1] + param * start2stop[1];
  point[2] = 0.0;

  // test if intersection pt is inside line B
  // edge = line B vector from v0 to v1
  // pvec = vector from either line B vertex to intersection point
  // if dot product of edge with pvec < or > 0.0 for each pvec,
  //   intersection point is outside line B
  // use EPSSQ and EPSSQNEG instead of 0.0 for following case:
  //   intersection pt is on line end pt where 2 lines come together
  //   want it to detect collision with at least one of lines
  //   point can be epsilon away from end pt
  //   this leads to pvec being epsilon vec in opposite dirs for 2 lines
  //   this can lead to dot3() being negative espilon^2 for both lines,
  //     depending on direction of 2 lines
  //   thus this can lead to no collision with either line
  //   typical observed dot values were 1.0e-18, so use EPSSQ = 1.0e-16

  MathExtra::sub3(v1,v0,edge);
  MathExtra::sub3(point,v0,pvec);
  if (MathExtra::dot3(edge,pvec) < EPSSQNEG) return false;
  MathExtra::sub3(point,v1,pvec);
  if (MathExtra::dot3(edge,pvec) > EPSSQ) return false;

  // there is a valid intersection with line B
  // set side to ONSUFR, OUTSIDE, or INSIDE
  // if start point is inside or outside then side = same
  // if particle started on line B, side = ONSURF OUT/IN based on dotstop

  if (dotstart < 0.0) side = INSIDE;
  else if (dotstart > 0.0) side = OUTSIDE;
  else if (dotstop > 0.0) side = ONSURF2OUT;
  else side = ONSURF2IN;

  return true;
}

/* ----------------------------------------------------------------------
   check for axisymmetric move crossing general line segment in (x,r) space
   line segment from v1 to v2 in axisymmetry plane
   3d move starting at x with v for tdelta
   solve quadratic eq to determine if non-linear trajectory intersects line seg
   return 1 if yes, 0 if no
   if yes, also return:
     param = fraction of tdelta at collision pt
     xc = x,y position of collision pt in axisymmetry plane
     vc = vy,vz velocity at collision pt in axisymmetry plane
     side = side of line that was hit = OUTSIDE,INSIDE,ONSURF2OUT,ONSURF2IN
------------------------------------------------------------------------- */

bool axi_line_intersect(double tdelta, double *x, double *v,
                        int outface, double *lo, double *hi,
                        double *v1, double *v2, double *norm,
                        double *xc, double *vc, 
                        double &param, int &side)
{
  double edge[3],pvec[3];
  double tc,frac;
        
  if (v1[0] == v2[0]) {    // vertical line segment
    if (outface == 0 && v1[0] == lo[0]) tc = tdelta;
    else if (outface == 1 && v1[0] == hi[0]) tc = tdelta;
    else {
      if (v[0] == 0.0) return false;
      tc = (v1[0] - x[0]) / v[0];
    }

  } else if (v1[1] == v2[1]) {    // horizontal line segment
    if (outface == 2 && v1[1] == lo[1]) tc = tdelta;
    else if (outface == 3 && v1[1] == hi[1]) tc = tdelta;
    else {
      if (axi_horizontal_line(tdelta,x,v,v1[1],frac)) tc = frac*tdelta;
      else return false;
    }

  } else {
    double x21 = v2[0] - v1[0];
    double y21 = v2[1] - v1[1];
    double x21sq = x21*x21;
    double y21sq = y21*y21;
    double dconst = x21*v1[1] - y21*v1[0];

    double a = x21sq*(v[1]*v[1] + v[2]*v[2]) - y21sq*v[0]*v[0];
    if (a == 0.0) return false;
    double b = x21sq*x[1]*v[1] - y21sq*x[0]*v[0] - y21*v[0]*dconst;
    double c = x21sq*x[1]*x[1] - y21sq*x[0]*x[0] - 
      2.0*y21*x[0]*dconst - dconst*dconst;

    double arg = b*b - a*c;
    if (arg < 0.0) return false;
    double sarg = sqrt(arg);

    double t1 = (-b + sarg) / a;
    double t2 = (-b - sarg) / a;

    tc = MIN(t1,t2);
    if (tc < 0.0) tc = MAX(t1,t2);
  }

  if (tc < 0.0) return false;
  if (tc > tdelta) return false;

  // set xc[0,1] and vc[1,2] to values in axisymmetric plane
  // insure xc is on vertical or horizontal surf

  xc[0] = x[0] + tc*v[0];
  if (v1[0] == v2[0]) xc[0] = v1[0];
  double ynew = x[1] + tc*v[1];
  double znew = x[2] + tc*v[2];
  xc[1] = sqrt(ynew*ynew + znew*znew);
  if (v1[1] == v2[1]) xc[1] = v1[1];
  xc[2] = 0.0;

  double rn = ynew / xc[1];
  double wn = znew / xc[1];
  vc[0] = v[0];
  vc[1] = v[1]*rn + v[2]*wn;
  vc[2] = -v[1]*wn + v[2]*rn;

  // test if intersection pt is inside line segment
  // edge = line vector from v1 to v2
  // pvec = vector from either line vertex to intersection point
  // if dot product of edge with pvec < or > 0.0 for each pvec,
  //   intersection point is outside line segment
  // use EPSSQ and EPSSQNEG instead of 0.0 for following case:
  //   intersection pt is on line end pt where 2 lines come together
  //   want it to detect collision with at least one of lines
  //   point can be epsilon away from end pt
  //   this leads to pvec being epsilon vec in opposite dirs for 2 lines
  //   this can lead to dot3() being negative espilon^2 for both lines,
  //     depending on direction of 2 lines
  //   thus this can lead to no collision with either line
  //   typical observed dot values were 1.0e-18, so use EPSSQ = 1.0e-16

  if (v1[0] <= v2[0]) {
    if (xc[0] < v1[0] || xc[0] > v2[0]) return false;
  } else {
    if (xc[0] < v2[0] || xc[0] > v1[0]) return false;
  }

  /*
  MathExtra::sub3(v2,v1,edge);
  MathExtra::sub3(xc,v1,pvec);

  if (MathExtra::dot3(edge,pvec) < EPSSQNEG) return false;
  MathExtra::sub3(xc,v2,pvec);
  if (MathExtra::dot3(edge,pvec) > EPSSQ) return false;
  */

  // there is a valid intersection with line segment
  // set side to OUTSIDE or INSIDE
  // no ONSURF case b/c surface is curved
  //  what if tc = param = 0.0 ??
  //  don't think particle should start on surf
  // in axisymmetric plane, particle path is curved
  // regardless of where particle starts, it can hit front or back or surf
  // use velocity vector at collision pt to determine side

  double dot = MathExtra::dot3(norm,vc);
  if (dot < 0.0) side = OUTSIDE;
  else side = INSIDE;

  param = tc/tdelta;
  return true;
}

/* ----------------------------------------------------------------------
   check for axisymmetric move crossing horizontal line in (x,r) space
   horizontal line is at yhoriz
   move starting at x with v for tdelta
   return 1 if crosses, 0 if not
   if crosses, also return frac = fraction of tdelta at crossing
------------------------------------------------------------------------- */

bool axi_horizontal_line(double tdelta, double *x, double *v, 
                         double yhoriz, double &frac)
{
  double a = v[1]*v[1] + v[2]*v[2];
  if (a == 0.0) return false;
  double b = -v[1]*x[1];
  double arg = yhoriz*yhoriz*a - v[2]*v[2]*x[1]*x[1];
  if (arg < 0.0) return false;
  double sarg = sqrt(arg);

  double t1 = (b + sarg) / a;
  double t2 = (b - sarg) / a;

  // if particle starts on line (e.g. due to cell crossing):
  // discard crossing at time = 0.0 or epsilon (due to round-off)

  double tc;
  if (x[1] == yhoriz) {
    if (fabs(t1) < fabs(t2)) tc = t2;
    else tc = t1;
  } else {
    tc = MIN(t1,t2);
    if (tc < 0.0) tc = MAX(t1,t2);
  }

  if (tc < 0.0 || tc > tdelta) return false;

  frac = tc/tdelta;
  return true;
}

/* ----------------------------------------------------------------------
   detect intersection between a directed line segment and a triangle
   intersection is defined as any line segment pt (including end pts)
     in common with any triangle pt (interior, edge, vertex)
   one exception is if both line end pts are in plane of triangle,
     then is NOT an intersection
   start,stop = end points of directed line segment, can have zero length
   v0,v1,v2 = 3 vertices of triangle
   norm = unit vector normal to triangle plane
     pointing OUTSIDE via right-hand rule
   return TRUE if there is an intersection, else FALSE
   if TRUE also return:
     point = pt of intersection
     param = intersection pt is this fraction along line (0-1 inclusive)
     side = side of B that was hit = OUTSIDE,INSIDE,ONSURF2OUT,ONSURF2IN
------------------------------------------------------------------------- */

bool line_tri_intersect(double *start, double *stop,
			double *v0, double *v1, double *v2, double *norm, 
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
  // force param to be 0.0 to 1.0 inclusive

  MathExtra::sub3(v0,start,vec);
  MathExtra::sub3(stop,start,start2stop);
  param = MathExtra::dot3(norm,vec) / MathExtra::dot3(norm,start2stop);
  param = MAX(param,0.0);
  param = MIN(param,1.0);

  // point = intersection pt with plane of triangle

  point[0] = start[0] + param * start2stop[0];
  point[1] = start[1] + param * start2stop[1];
  point[2] = start[2] + param * start2stop[2];

  // test if intersection pt is inside triangle
  // edge = edge vector of triangle
  // pvec = vector from triangle vertex to intersection point
  // xproduct = cross product of edge with pvec
  // if dot product of xproduct with norm < 0.0 for any of 3 edges,
  //   intersection point is outside tri
  // use EPSSQNEG instead of 0.0 for following case:
  //   intersection pt is on tri edge where 2 tris come together
  //   want it to detect collision with at least one of tris
  //   point can be epsilon away from edge
  //   this leads to xproduct being epsilon vec in opposite dirs for 2 tris
  //   this can lead to dot3() being negative espilon^2 for both tris,
  //     depending on direction of 2 tri norms
  //   thus this can lead to no collision with either tri
  //   typical observed dot values were -1.0e-18, so use EPSSQNEG = -1.0e-16

  MathExtra::sub3(v1,v0,edge);
  MathExtra::sub3(point,v0,pvec);
  MathExtra::cross3(edge,pvec,xproduct);
  if (MathExtra::dot3(xproduct,norm) < EPSSQNEG) return false;

  MathExtra::sub3(v2,v1,edge);
  MathExtra::sub3(point,v1,pvec);
  MathExtra::cross3(edge,pvec,xproduct);
  if (MathExtra::dot3(xproduct,norm) < EPSSQNEG) return false;

  MathExtra::sub3(v0,v2,edge);
  MathExtra::sub3(point,v2,pvec);
  MathExtra::cross3(edge,pvec,xproduct);
  if (MathExtra::dot3(xproduct,norm) < EPSSQNEG) return false;

  // there is a valid intersection with triangle
  // set side to ONSUFR, OUTSIDE, or INSIDE
  // if start point is inside or outside then side = same
  // if particle started on triangle, side = ONSURF OUT/IN based on dotstop

  if (dotstart < 0.0) side = INSIDE;
  else if (dotstart > 0.0) side = OUTSIDE;
  else if (dotstop > 0.0) side = ONSURF2OUT;
  else side = ONSURF2IN;

  return true;
}

/* ----------------------------------------------------------------------
   determine which side of plane the point x,y,z is on
   plane is defined by vertex pt v and unit normal vec
   return -1,0,1 for below,on,above plane
------------------------------------------------------------------------- */

int whichside(double *v, double *norm, double x, double y, double z)
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
   determine if point x lies on surface of hex defined by lo and hi
   return 1 if it does, 0 if not
------------------------------------------------------------------------- */

int point_on_hex(double *x, double *lo, double *hi)
{
  if ((x[0] == lo[0] || x[0] == hi[0]) && 
      x[1] >= lo[1] && x[1] <= hi[1] && x[2] >= lo[2] && x[2] <= hi[2])
    return 1;
  if ((x[1] == lo[1] || x[1] == hi[1]) && 
      x[0] >= lo[0] && x[0] <= hi[0] && x[2] >= lo[2] && x[2] <= hi[2])
    return 1;
  if ((x[2] == lo[2] || x[2] == hi[2]) && 
      x[0] >= lo[0] && x[0] <= hi[0] && x[1] >= lo[1] && x[1] <= hi[1])
    return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   determine if point x lies inside or on surface of hex defined by lo and hi
   return 1 if it does, 0 if not
------------------------------------------------------------------------- */

int point_in_hex(double *x, double *lo, double *hi)
{
  if (x[0] >= lo[0] && x[0] <= hi[0] && 
      x[1] >= lo[1] && x[1] <= hi[1] && 
      x[2] >= lo[2] && x[2] <= hi[2]) return 1;
  return 0;
}


/* ----------------------------------------------------------------------
   determine if point x lies on edge or inside of tri defined by p1,p2,p3
   return 1 if it does, 0 if not
------------------------------------------------------------------------- */

int point_in_tri(double *x, double *p1, double *p2, double *p3, double *norm)
{
  // if not in plane of tri, then not inside tri

  if (whichside(p1,norm,x[0],x[1],x[2])) return 0;

  // enorm123 = 3 vecs normal to each edge of tri
  // are in plane of tri, pointing towards center of tri
  // enorms are NOT unit vectors

  double enorm1[3],enorm2[3],enorm3[3];

  double diff[3];
  MathExtra::sub3(p2,p1,diff);
  MathExtra::cross3(norm,diff,enorm1);
  MathExtra::sub3(p3,p2,diff);
  MathExtra::cross3(norm,diff,enorm2);
  MathExtra::sub3(p1,p3,diff);
  MathExtra::cross3(norm,diff,enorm3);

  // if (pt - vertex) dotted into tri edge normal < 0, then outside tri

  MathExtra::sub3(p1,x,diff);
  if (MathExtra::dot3(diff,enorm1) < 0.0) return 0;
  MathExtra::sub3(p2,x,diff);
  if (MathExtra::dot3(diff,enorm2) < 0.0) return 0;
  MathExtra::sub3(p3,x,diff);
  if (MathExtra::dot3(diff,enorm3) < 0.0) return 0;
  return 1;
}

/* ----------------------------------------------------------------------
   compute distance bewteen a point X and line segment (P1,P2)
   distance = nearest distance to any point on line segment
------------------------------------------------------------------------- */

double distsq_point_line(double *x, double *p1, double *p2)
{
  // A = vector from P1 to X
  // B = vector from P1 to P2

  double a[3],b[3];
  MathExtra::sub3(x,p1,a);
  MathExtra::sub3(p2,p1,b);

  // let P = projected point on infinite P1 to P2 line that is closest to X
  // alpha = fraction of distance from P1 to P2 that P is at
  // alpha can be < 0, or between 0 to 1, or > 1
 
  double alpha = MathExtra::dot3(a,b)/MathExtra::lensq3(b);

  // reset A = vector from point on P1P2 line to X
  // if alpha < 0.0, point on line is P1
  // if alpha > 1.0, point on line is P2
  // else point on line is P1 + alpha*(P2-P1)

  if (alpha >= 1.0) MathExtra::sub3(x,p2,a);
  else if (alpha > 0.0) {
    a[0] = p1[0] + alpha*b[0];
    a[1] = p1[1] + alpha*b[1];
    a[2] = p1[2] + alpha*b[2];
  }

  // return length of A

  return MathExtra::lensq3(a);
}

/* ----------------------------------------------------------------------
   compute distance bewteen a point X and triangle (P1,P2,P3) with NORM
   distance = nearest distance to any point within 2d triangle surface
------------------------------------------------------------------------- */

double distsq_point_tri(double *x, double *p1, double *p2, double *p3,
                        double *norm)
{
  double a[3],point[3],edge[3],pvec[3],xproduct[3];

  // A = vector from P1 to X

  MathExtra::sub3(x,p1,a);

  // point = projected point on infinite triangle plane
  // pdistsq = projected distance to plane

  double alpha = MathExtra::dot3(a,norm);
  point[0] = x[0] - alpha*norm[0];
  point[1] = x[1] - alpha*norm[1];
  point[2] = x[2] - alpha*norm[2];

  MathExtra::sub3(x,point,a);
  double pdistsq = MathExtra::lensq3(a);

  // test if projected point is inside triangle
  // edge = edge vector of triangle
  // pvec = vector from triangle vertex to projected point
  // xproduct = cross product of edge with pvec
  // if dot product of xproduct with norm < 0.0 for any of 3 edges,
  //   projected point is outside tri
  // if inside, return projected distance to plane

  int inside = 1;

  MathExtra::sub3(p2,p1,edge);
  MathExtra::sub3(point,p1,pvec);
  MathExtra::cross3(edge,pvec,xproduct);
  if (MathExtra::dot3(xproduct,norm) < 0.0) inside = 0;

  MathExtra::sub3(p3,p2,edge);
  MathExtra::sub3(point,p2,pvec);
  MathExtra::cross3(edge,pvec,xproduct);
  if (MathExtra::dot3(xproduct,norm) < 0.0) inside = 0;

  MathExtra::sub3(p1,p3,edge);
  MathExtra::sub3(point,p3,pvec);
  MathExtra::cross3(edge,pvec,xproduct);
  if (MathExtra::dot3(xproduct,norm) < 0.0) inside = 0;

  if (inside) return pdistsq;

  // projected point is outside triangle
  // compute minimum distance to any of 3 triangle edges
  // return sum of min distance and projected distance

  double rsq = distsq_point_line(point,p1,p2);
  rsq = MIN(rsq,distsq_point_line(point,p2,p3));
  rsq = MIN(rsq,distsq_point_line(point,p3,p1));
  return rsq + pdistsq;
}

/* ----------------------------------------------------------------------
   return minimum fractional distance that X is from line pts V0 or V1
   lensq = length of line
   frac = (length of X-V0 or X-V1) / lensq
   return fracsq = square of frac
------------------------------------------------------------------------- */

double line_fraction(double *x, double *v0, double *v1)
{
  double segment[3];

  MathExtra::sub3(v0,v1,segment);
  double lensq = MathExtra::lensq3(segment);

  MathExtra::sub3(x,v0,segment);
  double fracsq = MathExtra::lensq3(segment)/lensq;
  MathExtra::sub3(x,v1,segment);
  fracsq = MIN(fracsq,MathExtra::lensq3(segment)/lensq);

  return fracsq;
}

/* ----------------------------------------------------------------------
   return minimum fractional distance that X is from triangle pts V0, V1, V2
   lensq = min length of any of 3 triangle edges
   frac = (length of X-V0 or X-V1 or X-V2) / lensq
   return fracsq = square of frac
------------------------------------------------------------------------- */

double tri_fraction(double *x, double *v0, double *v1, double *v2)
{
  double segment[3];

  MathExtra::sub3(v0,v1,segment);
  double lensq = MathExtra::lensq3(segment);
  MathExtra::sub3(v1,v2,segment);
  lensq = MIN(lensq,MathExtra::lensq3(segment));
  MathExtra::sub3(v0,v2,segment);
  lensq = MIN(lensq,MathExtra::lensq3(segment));

  MathExtra::sub3(x,v0,segment);
  double fracsq = MathExtra::lensq3(segment)/lensq;
  MathExtra::sub3(x,v1,segment);
  fracsq = MIN(fracsq,MathExtra::lensq3(segment)/lensq);
  MathExtra::sub3(x,v2,segment);
  fracsq = MIN(fracsq,MathExtra::lensq3(segment)/lensq);

  return fracsq;
}

/* ---------------------------------------------------------------------- */

}

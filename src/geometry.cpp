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

#include "geometry.h"
#include "math_extra.h"

enum{NONE,OUTSIDE,INSIDE,BOTH};

namespace Geom {

/* ----------------------------------------------------------------------
   compute intersection of a line with a quad cell
   intersection is defined as
     any line pt (interior, vertex) in common with
     any rectangle pt (interior, edge, vertex)
   v0,v1 and norm = 2 vertices of line and unit normal vec
   lo,hi = opposite corner pts of rectangle
   return 1 if intersection, else 0
------------------------------------------------------------------------- */

int line_quad_intersect(double *v0, double *v1, double *norm,
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

int tri_hex_intersect(double *v0, double *v1, double *v2, double *norm,
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

bool line_line_intersect(double *v0, double *v1,
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

bool tri_line_intersect(double *v0, double *v1, double *v2,
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

int lineside(double *v, double *norm, double x, double y)
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

int triside(double *v, double *norm, double x, double y, double z)
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

/* ---------------------------------------------------------------------- */

}

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

#ifndef SPARTA_GEOMETRY_H
#define SPARTA_GEOMETRY_H

namespace Geometry {
  int line_quad_intersect(double *, double *, double *,
			  double *, double *);
  int quad_line_intersect_point(double *, double *, double *,
				double *, double *, double *);
  int line_touch_quad_face(double *, double *, int, double *, double *);

  int tri_hex_intersect(double *, double *, double *, double *,
			double *, double *);
  int hex_tri_intersect_point(double *, double *, double *, double *,
			      double *, double *, double *);
  int tri_touch_hex_face(double *, double *, double *, int, double *, double *);
  int tri_on_hex_face(double *, double *, double *, double *, double *);
  int edge_on_hex_face(double *, double *, double *, double *);

  bool line_line_intersect(double *, double *,
			   double *, double *, double *,
			   double *, double &param, int &, int=0);

  bool axi_line_intersect(double, double *, double *, int, double *, double *,
                          double *, double *, double *, int,
                          double *, double *, double &, int &);
  bool axi_horizontal_line(double, double *, double *, double,
                           int &, double &, double &);

  bool line_tri_intersect(double *, double *,
			  double *, double *, double *, double *,
			  double *, double &param, int &);
  int whichside(double *, double *, double, double, double);
  int point_on_hex(double *, double *, double *);
  int point_in_hex(double *, double *, double *);
  int point_in_tri(double *, double *, double *, double *, double *);

  double distsq_point_line(double *, double *, double *);
  double distsq_point_tri(double *, double *, double *, double *, double *);

  double dist_line_quad(double *, double *, double *, double *);
  double dist_tri_hex(double *, double *, double *, double *,
                      double *, double *);

  double line_fraction(double *, double *, double *);
  double tri_fraction(double *, double *, double *, double *);
}


// legal banner and code for SNL-modified PQP (Proximity Query Package) follow:

namespace Geometry_PQP {  
/*************************************************************************\

  Copyright 1999 The University of North Carolina at Chapel Hill.
  All Rights Reserved.

  Permission to use, copy, modify and distribute this software and its
  documentation for educational, research and non-profit purposes, without
  fee, and without a written agreement is hereby granted, provided that the
  above copyright notice and the following three paragraphs appear in all
  copies.

  IN NO EVENT SHALL THE UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL BE
  LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
  CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE
  USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY
  OF NORTH CAROLINA HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH
  DAMAGES.

  THE UNIVERSITY OF NORTH CAROLINA SPECIFICALLY DISCLAIM ANY
  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE
  PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
  NORTH CAROLINA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT,
  UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

  The authors may be contacted via:

  US Mail:             E. Larsen
                       Department of Computer Science
                       Sitterson Hall, CB #3175
                       University of N. Carolina
                       Chapel Hill, NC 27599-3175

  Phone:               (919)962-1749

  EMail:               geom@cs.unc.edu


\**************************************************************************/
typedef double PQP_REAL;
 
void
SegPoints(PQP_REAL VEC[3], 
	  PQP_REAL X[3], PQP_REAL Y[3],
          const PQP_REAL P[3], const PQP_REAL A[3],
          const PQP_REAL Q[3], const PQP_REAL B[3]);
  
PQP_REAL 
TriDist(PQP_REAL p[3], PQP_REAL q[3], 
	const PQP_REAL s[3][3], const PQP_REAL t[3][3]);

inline
void
VmV(PQP_REAL Vr[3], const PQP_REAL V1[3], const PQP_REAL V2[3])
{
  Vr[0] = V1[0] - V2[0];
  Vr[1] = V1[1] - V2[1];
  Vr[2] = V1[2] - V2[2];
}

inline
void
VpV(PQP_REAL Vr[3], const PQP_REAL V1[3], const PQP_REAL V2[3])
{
  Vr[0] = V1[0] + V2[0];
  Vr[1] = V1[1] + V2[1];
  Vr[2] = V1[2] + V2[2];
}
 
inline
PQP_REAL
VdotV(const PQP_REAL V1[3], const PQP_REAL V2[3])
{
  return (V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2]);
}

inline
void
VcV(PQP_REAL Vr[3], const PQP_REAL V[3])
{
  Vr[0] = V[0];  Vr[1] = V[1];  Vr[2] = V[2];
}

inline
void
VcrossV(PQP_REAL Vr[3], const PQP_REAL V1[3], const PQP_REAL V2[3])
{
  Vr[0] = V1[1]*V2[2] - V1[2]*V2[1];
  Vr[1] = V1[2]*V2[0] - V1[0]*V2[2];
  Vr[2] = V1[0]*V2[1] - V1[1]*V2[0];
}

inline
void
VpVxS(PQP_REAL Vr[3], const PQP_REAL V1[3], const PQP_REAL V2[3], PQP_REAL s)
{
  Vr[0] = V1[0] + V2[0] * s;
  Vr[1] = V1[1] + V2[1] * s;
  Vr[2] = V1[2] + V2[2] * s;
}

inline
void
VxS(PQP_REAL Vr[3], const PQP_REAL V[3], PQP_REAL s)
{
  Vr[0] = V[0] * s;
  Vr[1] = V[1] * s;
  Vr[2] = V[2] * s;
}

inline
PQP_REAL
VdistV2(const PQP_REAL V1[3], const PQP_REAL V2[3])
{
  return ( (V1[0]-V2[0]) * (V1[0]-V2[0]) + 
	   (V1[1]-V2[1]) * (V1[1]-V2[1]) + 
	   (V1[2]-V2[2]) * (V1[2]-V2[2]));
} 
  
}

#endif



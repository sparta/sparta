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

#ifndef DSMC_GEOMETRY_H
#define DSMC_GEOMETRY_H

namespace Geometry {
  int line_quad_intersect(double *, double *, double *,
			  double *, double *);
  int quad_line_intersect_point(double *, double *, double *,
				double *, double *, double *);
  int line_quad_face_touch(double *, double *, int, double *, double *);

  int tri_hex_intersect(double *, double *, double *, double *,
			double *, double *);
  int hex_tri_intersect_point(double *, double *, double *, double *,
			      double *, double *, double *);
  int tri_hex_face_touch(double *, double *, double *, int, double *, double *);

  bool line_line_intersect(double *, double *, 
			   double *, double *, double *,
			   double *, double &param, int &);
  bool line_tri_intersect(double *, double *, 
			  double *, double *, double *, double *,
			  double *, double &param, int &);
  int whichside(double *, double *, double, double, double);
  int point_on_hex(double *, double *, double *);
  int point_in_hex(double *, double *, double *);

  double line_fraction(double *, double *, double *);
  double tri_fraction(double *, double *, double *, double *);
}

#endif

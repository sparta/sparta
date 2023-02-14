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

#endif

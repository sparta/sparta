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

#ifndef DSMC_SURF_H
#define DSMC_SURF_H

#include "pointers.h"

namespace DSMC_NS {

class Surf : protected Pointers {
 public:
  int surf_exist;

  struct Point {
    double x[3];
  };

  struct Line {
    int id;                // index of surface model it belongs to ???
    int p1,p2;                // indices of points in line segment
                              // rhand rule: z x (p2-p1) = outward normal
    double norm[3];           // outward normal to line segment
  };

  struct Tri {
    int id;                // index of surface model it belongs to
    int p1,p2,p3;             // indices of points in triangle
                              // rhand rule: (p2-p1) x (p3-p1) = outward normal
    double norm[3];           // outward normal to triangle
  };

  Point *pts;
  Line *lines;
  Tri *tris;
  int npoint;
  int nline;
  int ntri;

  class SurfCollide **sc;      // list of surface collision models
  int nsc;                     // # of surface collision models
  int maxsc;

  Surf(class DSMC *);
  ~Surf();
  void init() {}
  int add_id(char *);
  void compute_line_normal(int, int);
  void compute_tri_normal(int, int);
  void add_collide(int, char **);
  int find_collide(const char *);
  bigint memory_usage();

 private:
  int nsurf;                    // # of read-in surfaces, each with ID
  char **ids;                   // read-in surface IDs

};

}

#endif

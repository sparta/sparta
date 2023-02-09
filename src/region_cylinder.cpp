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
#include "stdlib.h"
#include "string.h"
#include "region_cylinder.h"
#include "error.h"

using namespace SPARTA_NS;

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

RegCylinder::RegCylinder(SPARTA *sparta, int narg, char **arg) :
  Region(sparta, narg, arg)
{
  options(narg-8,&arg[8]);

  if (strcmp(arg[2],"x") && strcmp(arg[2],"y") && strcmp(arg[2],"z"))
    error->all(FLERR,"Illegal region cylinder command");

  axis = arg[2][0];
  c1 = atof(arg[3]);
  c2 = atof(arg[4]);
  radius = atof(arg[5]);

  if (strcmp(arg[6],"INF") == 0) lo = -BIG;
  else lo = atof(arg[6]);
  if (strcmp(arg[7],"INF") == 0) hi = BIG;
  else hi = atof(arg[7]);

  // error check

  if (radius <= 0.0) error->all(FLERR,"Illegal region cylinder command");

  // extent of cylinder
  // for variable radius, uses initial radius

  if (interior) {
    bboxflag = 1;
    if (axis == 'x') {
      extent_xlo = lo;
      extent_xhi = hi;
      extent_ylo = c1 - radius;
      extent_yhi = c1 + radius;
      extent_zlo = c2 - radius;
      extent_zhi = c2 + radius;
    }
    if (axis == 'y') {
      extent_xlo = c1 - radius;
      extent_xhi = c1 + radius;
      extent_ylo = lo;
      extent_yhi = hi;
      extent_zlo = c2 - radius;
      extent_zhi = c2 + radius;
    }
    if (axis == 'z') {
      extent_xlo = c1 - radius;
      extent_xhi = c1 + radius;
      extent_ylo = c2 - radius;
      extent_yhi = c2 + radius;
      extent_zlo = lo;
      extent_zhi = hi;
    }
  } else bboxflag = 0;
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int RegCylinder::inside(double *x)
{
  double del1,del2,dist;
  int inside;

  if (axis == 'x') {
    del1 = x[1] - c1;
    del2 = x[2] - c2;
    dist = sqrt(del1*del1 + del2*del2);
    if (dist <= radius && x[0] >= lo && x[0] <= hi) inside = 1;
    else inside = 0;
  } else if (axis == 'y') {
    del1 = x[0] - c1;
    del2 = x[2] - c2;
    dist = sqrt(del1*del1 + del2*del2);
    if (dist <= radius && x[1] >= lo && x[1] <= hi) inside = 1;
    else inside = 0;
  } else {
    del1 = x[0] - c1;
    del2 = x[1] - c2;
    dist = sqrt(del1*del1 + del2*del2);
    if (dist <= radius && x[2] >= lo && x[2] <= hi) inside = 1;
    else inside = 0;
  }

  return inside;
}

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
#include "region_plane.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

RegPlane::RegPlane(SPARTA *sparta, int narg, char **arg) :
  Region(sparta, narg, arg)
{
  options(narg-8,&arg[8]);

  xp = atof(arg[2]);
  yp = atof(arg[3]);
  zp = atof(arg[4]);
  normal[0] = atof(arg[5]);
  normal[1] = atof(arg[6]);
  normal[2] = atof(arg[7]);

  // enforce unit normal

  double rsq = normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2];
  if (rsq == 0.0) error->all(FLERR,"Illegal region plane command");
  normal[0] /= sqrt(rsq);
  normal[1] /= sqrt(rsq);
  normal[2] /= sqrt(rsq);

  // plane has no bounding box

  bboxflag = 0;
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is on normal side of plane or on plane
   inside = 0 if x,y,z is on non-normal side of plane and not on plane
   x,y,z is inside if (x-xp) dot normal >= 0
------------------------------------------------------------------------- */

int RegPlane::inside(double *x)
{
  double dot = (x[0]-xp)*normal[0] + (x[1]-yp)*normal[1] + (x[2]-zp)*normal[2];

  if (dot >= 0.0) return 1;
  return 0;
}

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

#ifdef SURF_COLLIDE_DIFFUSE_CLASS

SurfCollideStyle(diffuse,SurfCollideDiffuse)

#else

#ifndef DSMC_SURF_COLLIDE_DIFFUSE_H
#define DSMC_SURF_COLLIDE_DIFFUSE_H

#include "surf_collide.h"

namespace DSMC_NS {

class SurfCollideDiffuse : public SurfCollide {
 public:
  SurfCollideDiffuse(class DSMC *, int, char **);
  ~SurfCollideDiffuse();
  void init();
  void collide(Particle::OnePart *, double *);

  void dynamic();

 private:
  double twall;              // surface temperature
  double acc;                // surface accomodation coeff
  double vx,vy,vz;           // translational velocity of surface
  double wx,wy,wz;           // angular velocity of surface
  double px,py,pz;           // point to rotate surface around
  int tflag,rflag;           // flags for translation and rotation
  int trflag;                // 1 if either tflag or rflag is set

  char *tstr;                // temperature variable name (NULL if constant)
  int tvar;                  // index of equal-style variable

  class RanPark *random;     // RNG for particle reflection
};

}

#endif
#endif

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

#ifdef SURF_COLLIDE_CLASS

SurfCollideStyle(diffuse,SurfCollideDiffuse)

#else

#ifndef SPARTA_SURF_COLLIDE_DIFFUSE_H
#define SPARTA_SURF_COLLIDE_DIFFUSE_H

#include "surf_collide.h"
#include "surf.h"

namespace SPARTA_NS {

class SurfCollideDiffuse : public SurfCollide {
 public:
  SurfCollideDiffuse(class SPARTA *, int, char **);
  SurfCollideDiffuse(class SPARTA *sparta) : SurfCollide(sparta) {}
  ~SurfCollideDiffuse();
  void init();
  Particle::OnePart *collide(Particle::OnePart *&, double *, double &,
                             int, int &, int);

  void dynamic();

 protected:
  double twall;              // surface temperature
  double acc;                // surface accomodation coeff
  double vx,vy,vz;           // translational velocity of surface
  double wx,wy,wz;           // angular velocity of surface
  double px,py,pz;           // point to rotate surface around
  int tflag,rflag;           // flags for translation and rotation
  int trflag;                // 1 if either tflag or rflag is set

  char *tstr;                // temperature variable name (NULL if constant)
  int tvar;                  // index of equal-style variable

  Surf::Line *lines;
  Surf::Tri *tris;

  int distributed,implicit;  // Surf settings
  int nsurf;                 // nlocal or nown

  double vstream[3];
  class RanKnuth *random;     // RNG for particle reflection

  void diffuse(Particle::OnePart *, double *, int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Surf_collide diffuse rotation invalid for 2d

Specified rotation vector must be in z-direction.

E: Surf_collide diffuse variable name does not exist

Self-explanatory.

E: Surf_collide diffuse variable is invalid style

It must be an equal-style variable.

*/

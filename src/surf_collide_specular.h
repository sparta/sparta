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

#ifdef SURF_COLLIDE_CLASS

SurfCollideStyle(specular,SurfCollideSpecular)

#else

#ifndef DSMC_SURF_COLLIDE_SPECULAR_H
#define DSMC_SURF_COLLIDE_SPECULAR_H

#include "surf_collide.h"

namespace DSMC_NS {

class SurfCollideSpecular : public SurfCollide {
 public:
  SurfCollideSpecular(class DSMC *, int, char **);
  ~SurfCollideSpecular() {}
  void init() {}
  void collide(Particle::OnePart *, double *);
};

}

#endif
#endif

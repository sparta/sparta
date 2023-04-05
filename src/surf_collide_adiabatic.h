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

#ifdef SURF_COLLIDE_CLASS

SurfCollideStyle(adiabatic,SurfCollideAdiabatic)

#else

#ifndef SPARTA_SURF_COLLIDE_ADIABATIC_H
#define SPARTA_SURF_COLLIDE_ADIABATIC_H

#include "surf_collide.h"

namespace SPARTA_NS {

class SurfCollideAdiabatic : public SurfCollide {
 public:
  SurfCollideAdiabatic(class SPARTA *, int, char **);
  SurfCollideAdiabatic(class SPARTA *sparta) : SurfCollide(sparta) {}
  ~SurfCollideAdiabatic();
  Particle::OnePart *collide(Particle::OnePart *&, double &,
                             int, double *, int, int &);
  void wrapper(Particle::OnePart *, double *, int *, double*);
  void flags_and_coeffs(int *, double *) {}

 protected:
  class RanKnuth *random; // RNG for particle reflection

  void scatter_isotropic(Particle::OnePart *, double *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

*/

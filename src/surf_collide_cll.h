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

SurfCollideStyle(cll,SurfCollideCLL)

#else

#ifndef SPARTA_SURF_COLLIDE_CLL_H
#define SPARTA_SURF_COLLIDE_CLL_H

#include "pointers.h"
#include "surf_collide.h"

namespace SPARTA_NS {

class SurfCollideCLL : public SurfCollide {
 public:
  SurfCollideCLL(class SPARTA *, int, char **);
  ~SurfCollideCLL();
  void init();
  Particle::OnePart *collide(Particle::OnePart *&, double &,
                             int, double *, int, int &);
  void wrapper(Particle::OnePart *, double *, int *, double*);
  void flags_and_coeffs(int *, double *);

  void dynamic();

 private:
  double twall;                         // surface temperature
  double acc_n,acc_t,acc_rot,acc_vib;   // surface accomodation coeffs
  double vx,vy,vz;                      // translational velocity of surface
  double wx,wy,wz;                      // angular velocity of surface
  double px,py,pz;                      // point to rotate surface around
  double eccen;                         // 1 if fully diffuse scattering
                                        // < 1 if partial diffuse scattering

  int tflag,rflag;           // flags for translation and rotation
  int trflag;                // 1 if either tflag or rflag is set
  int pflag;                 // 1 if partially energy accommodation
                             // with partial/fully diffuse scattering

  int tmode;                 // Twall is NUMERIC,VARIABLE,CUSTOM
  char *tstr;                // temperature variable name (NULL if constant)
  int tvar;                  // index of equal-style variable
  double *tvector;           // custom per-surf temperature vector

  double vstream[3];
  class RanKnuth *random;     // RNG for particle reflection

  void cll(Particle::OnePart *, double *);
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

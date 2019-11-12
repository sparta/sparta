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

SurfCollideStyle(impulsive,SurfCollideImpulsive)

#else

#ifndef SPARTA_SURF_COLLIDE_IMPULSIVE_H
#define SPARTA_SURF_COLLIDE_IMPULSIVE_H

#include "surf_collide.h"

namespace SPARTA_NS {

class SurfCollideImpulsive : public SurfCollide {
 public:
  SurfCollideImpulsive(class SPARTA *, int, char **);
  ~SurfCollideImpulsive();
  void init();
  Particle::OnePart *collide(Particle::OnePart *&, double *, double &, int);
  
  void dynamic();

 private:
  double twall;                   // surface temperature
  double eng_ratio,eff_mass;      // energy ratio and effective mass
                                  // of the surface for soft-sphere model
  double var_alpha,var_alpha_sq;  // alpha value related to variance
                                  // from Rettner's expression 
  double theta_peak,cos_theta_pow;   // cosine power law varaition for theta
  double cos_phi_pow;                // cosine power law varaition for phi
  double step_size,cos_theta_pow_2;  // step_size and double cosine power
  double slow_frac,slow_u0_a,slow_u0_b;  // slow IS model
  double slow_u0,slow_alpha,slow_alpha_sq,slow_barrier;
  double vx,vy,vz;                 // translational velocity of surface
  double wx,wy,wz;                 // angular velocity of surface
  double px,py,pz;                 // point to rotate surface around
  int tflag,rflag;                 // flags for translation and rotation
  int trflag;                      // 1 if either tflag or rflag is set
  int step_flag,double_flag,slow_flag;   // optional model flags

  char *tstr;                // temperature variable name (NULL if constant)
  int tvar;                  // index of equal-style variable

  double vstream[3];
  class RanPark *random;     // RNG for particle reflection

  void impulsive(Particle::OnePart *, double *);
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

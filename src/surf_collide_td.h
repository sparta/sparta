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

SurfCollideStyle(td,SurfCollideTD)

#else

#ifndef SPARTA_SURF_COLLIDE_TD_H
#define SPARTA_SURF_COLLIDE_TD_H

#include "surf_collide.h"

namespace SPARTA_NS {

class SurfCollideTD : public SurfCollide {
 public:
  SurfCollideTD(class SPARTA *, int, char **);
  ~SurfCollideTD();
  void init();
  Particle::OnePart *collide(Particle::OnePart *&, double &,
                             int, double *, int, int &);
  void wrapper(Particle::OnePart *, double *, int *, double*);
  void flags_and_coeffs(int *, double *);
  void dynamic();

 private:
  double twall;              // surface temperature

  double barrier_val;
  double initen_trans, initen_rot, initen_vib;
  double bond_trans, bond_rot, bond_vib;

  double vx,vy,vz;           // translational velocity of surface
  double wx,wy,wz;           // angular velocity of surface
  double px,py,pz;           // point to rotate surface around

  int barrier_flag, initen_flag,bond_flag; // optional flags

  int tmode;                 // Twall is NUMERIC,VARIABLE,CUSTOM
  char *tstr;                // temperature variable name (NULL if constant)
  int tvar;                  // index of equal-style variable
  double *tvector;           // custom per-surf temperature vector

  double vstream[3];
  class RanKnuth *random;     // RNG for particle reflection

  void td(Particle::OnePart *, double *);
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

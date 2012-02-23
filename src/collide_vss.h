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

#ifdef COLLIDE_CLASS

CollideStyle(vss,CollideVSS)

#else

#ifndef DSMC_COLLIDE_VSS_H
#define DSMC_COLLIDE_VSS_H

#include "collide.h"
#include "particle.h"

namespace DSMC_NS {

class CollideVSS : public Collide {
 public:
  CollideVSS(class DSMC *, int, char **);
  ~CollideVSS();
  void init();

  double attempt_collision(int, int, int, double);
  int test_collision(int, int, int, Particle::OnePart *, Particle::OnePart *);
  void setup_collision(Particle::OnePart *, Particle::OnePart *);
  Particle::OnePart *perform_collision(Particle::OnePart *, 
				       Particle::OnePart *);

 private:
  int eng_exchange;
  double vr_indice;
  double **vrm; // static portion of max collision frequency
  double **prefactor; // static portion of collision attempt frequency
  double ***vremax;   // max relative velocity, per cell, per species pair

  void SCATTER_TwoBodyScattering(Particle::OnePart *, 
				 Particle::OnePart *);
  void EEXCHANGE_NonReactingEDisposal(Particle::OnePart *, 
				      Particle::OnePart *);

  struct State
  {
   double    vr2;
   double    vr;
   double    mr;
   double    rotdof_i;
   double    rotdof_j;
   double    vibdof_i;
   double    vibdof_j;
   double    ave_rotdof;
   double    ave_vibdof;
   double    ave_dof;
   double    etrans;
   double    erot;
   double    evib;
   double    eexchange;
   double    eint;
   double    etotal;

   double    mass_i;
   double    mass_j;
  };
 

 struct State precoln;
 struct State postcoln;

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running DSMC to see the offending line.

*/

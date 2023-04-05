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

#ifdef FIX_CLASS

FixStyle(ambipolar,FixAmbipolar)

#else

#ifndef SPARTA_FIX_AMBIPOLAR_H
#define SPARTA_FIX_AMBIPOLAR_H

#include "fix.h"
#include "particle.h"

namespace SPARTA_NS {

class FixAmbipolar : public Fix {
 public:
  int especies;               // index of electron species
  int *ions;                  // 1 if a particle species is an ionx

  FixAmbipolar(class SPARTA *, int, char **);
  FixAmbipolar(class SPARTA *sparta) : Fix(sparta) {} // needed for Kokkos
  virtual ~FixAmbipolar();
  int setmask();
  void init();
  virtual void update_custom(int, double, double, double, double *);
  void surf_react(Particle::OnePart *, int &, int &);

 protected:
  int maxion;                 // length of ions vector
  int ionindex,velindex;      // indices into particle custom data structs
  class RanKnuth *random;
};

}

#endif
#endif

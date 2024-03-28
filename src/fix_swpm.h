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

FixStyle(swpm,FixSWPM)

#else

#ifndef SPARTA_FIX_SWPM_H
#define SPARTA_FIX_SWPM_H

#include "fix.h"
#include "particle.h"

namespace SPARTA_NS {

class FixSWPM : public Fix {
 public:
  int index_swpm;

  FixSWPM(class SPARTA *, int, char **);
  FixSWPM(class SPARTA *sparta) : Fix(sparta) {} // needed for Kokkos
  virtual ~FixSWPM();
  int setmask();
  void init();
  virtual void update_custom(int, double, double, double, double *);

 protected:
  class RanKnuth *random;
};

}

#endif
#endif

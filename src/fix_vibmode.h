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

#ifdef FIX_CLASS

FixStyle(vibmode,FixVibmode)

#else

#ifndef SPARTA_FIX_VIBMODE_H
#define SPARTA_FIX_VIBMODE_H

#include "stdio.h"
#include "fix.h"

namespace SPARTA_NS {

class FixVibmode : public Fix {
 public:
  FixVibmode(class SPARTA *, int, char **);
  FixVibmode(class SPARTA *sparta) : Fix(sparta) {} // needed for Kokkos
  virtual ~FixVibmode();
  int setmask();
  void init();
  virtual void update_custom(int, double, double, double, double *);

 protected:
  int maxmode;           // max # of vibrational modes for any species
  int vibmodeindex;      // index into particle custom data structs
  class RanKnuth *random;
};

}

#endif
#endif

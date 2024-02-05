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

FixStyle(elecmode,FixElecmode)

#else

#ifndef SPARTA_FIX_ELECMODE_H
#define SPARTA_FIX_ELECMODE_H

#include "stdio.h"
#include "fix.h"

namespace SPARTA_NS {

class FixElecmode : public Fix {
 public:
  FixElecmode(class SPARTA *, int, char **);
  FixElecmode(class SPARTA *sparta) : Fix(sparta) {} // needed for Kokkos
  virtual ~FixElecmode();
  int setmask();
  void init();
  virtual void update_custom(int, double, double, double, double, double *);

 protected:
  int elecstateindex;      // index into particle custom data structs
  int eelecindex;
  class RanKnuth *random;
};

}

#endif
#endif

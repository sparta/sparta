/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(stochastic_weight,FixStochasticWeight)

#else

#ifndef SPARTA_FIX_STOCHASTIC_WEIGHT_H
#define SPARTA_FIX_STOCHASTIC_WEIGHT_H

#include "fix.h"
#include "particle.h"

namespace SPARTA_NS {

class FixStochasticWeight : public Fix {
 public:
  FixStochasticWeight(class SPARTA *, int, char **);
  FixStochasticWeight(class SPARTA *sparta) : Fix(sparta) {} // needed for Kokkos
  virtual ~FixStochasticWeight();
  int setmask();
  void init();
  virtual void update_custom(int, double, double, double, double *);

 protected:
  int stochastic_wt_index;       // index into particle custom data structs
};

}

#endif
#endif


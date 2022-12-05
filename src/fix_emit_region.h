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

FixStyle(emit/region,FixEmitRegion)

#else

#ifndef SPARTA_FIX_EMIT_REGION_H
#define SPARTA_FIX_EMIT_RECION_H

#include "fix.h"

namespace SPARTA_NS {

class FixEmitRegion : public Fix {
 public:
  FixEmitRegion(class SPARTA *, int, char **);
  virtual ~FixEmitRegion();
  int setmask();
  virtual void init();
  void start_of_step();
  double compute_vector(int);

 protected:
  int imix,groupbit;
  int mode;
  
  class RanKnuth *random;

  int nsingle,ntotal;
  bigint np;
  double voltotal;
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/

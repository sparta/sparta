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

#ifndef DSMC_RANPARK_H
#define DSMC_RANPARK_H

#include "pointers.h"

namespace DSMC_NS {

class RanPark : protected Pointers {
  friend class Set;
 public:
  RanPark(class DSMC *, int);
  double uniform();
  double gaussian();
  void reset(int);
  void reset(int, double *);
  int state();

 private:
  int seed,save;
  double second;
};

}

#endif

/* ERROR/WARNING messages:

E: Invalid seed for Park random # generator

The initial seed for this random number generator must be a positive
integer.

*/

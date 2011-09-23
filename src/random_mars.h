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

#ifndef DSMC_RANMARS_H
#define DSMC_RANMARS_H

#include "pointers.h"

namespace DSMC_NS {

class RanMars : protected Pointers {
 public:
  RanMars(class DSMC *, int);
  ~RanMars();
  double uniform();
  double gaussian();

 private:
  int seed,save;
  double second;
  double *u;
  int i97,j97;
  double c,cd,cm;
};

}

#endif

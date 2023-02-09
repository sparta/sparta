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

#ifndef SPARTA_MARCHING_SQUARES_H
#define SPARTA_MARCHING_SQUARES_H

#include "pointers.h"

namespace SPARTA_NS {

class MarchingSquares : protected Pointers {
 public:
  MarchingSquares(class SPARTA *, int, double);
  ~MarchingSquares() {}
  void invoke(double **, int *);

 private:
  int ggroup;
  double thresh;

  double interpolate(double, double, double, double);
};

}

#endif

/* ERROR/WARNING messages:

*/

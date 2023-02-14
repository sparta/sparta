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

#ifndef SPARTA_RAN_MARS_H
#define SPARTA_RAN_MARS_H

#include "pointers.h"

namespace SPARTA_NS {

class RanMars : protected Pointers {
 public:
  RanMars(class SPARTA *);
  ~RanMars();
  void init(int);
  double uniform();
  double gaussian();

 private:
  int initflag,save;
  int i97,j97;
  double c,cd,cm;
  double second;
  double *u;
};

}

#endif

/* ERROR/WARNING messages:

E: Seed command has not been used

This command should appear near the beginning of your input script,
before any random numbers are needed by other commands.

*/

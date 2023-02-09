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

#ifdef REGION_CLASS

RegionStyle(cylinder,RegCylinder)

#else

#ifndef SPARTA_REGION_CYLINDER_H
#define SPARTA_REGION_CYLINDER_H

#include "region.h"

namespace SPARTA_NS {

class RegCylinder : public Region {
 public:
  RegCylinder(class SPARTA *, int, char **);
  int inside(double *);

 private:
  char axis;
  double c1,c2;
  double radius;
  double lo,hi;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

*/

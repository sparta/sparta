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

#ifdef REGION_CLASS

RegionStyle(union,RegUnion)

#else

#ifndef SPARTA_REGION_UNION_H
#define SPARTA_REGION_UNION_H

#include "region.h"

namespace SPARTA_NS {

class RegUnion : public Region {
 public:
  RegUnion(class SPARTA *, int, char **);
  ~RegUnion();
  int inside(double *);

 private:
  int nregion;
  int *list;
};

}

#endif
#endif

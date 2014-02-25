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

#ifndef SPARTA_REACT_H
#define SPARTA_REACT_H

#include "pointers.h"
#include "particle.h"

namespace SPARTA_NS {

class React : protected Pointers {
 public:
  char *style;

  React(class SPARTA *, int, char **);
  virtual ~React();
  virtual void init() {}
  virtual int attempt(Particle::OnePart *, Particle::OnePart *, 
                      double, double, double, double &, int &) = 0;

 protected:
  class RanPark *random;
};

}

#endif

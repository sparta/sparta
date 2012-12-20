/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2012) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifndef SPARTA_MYVEC_H
#define SPARTA_MYVEC_H

#include "stdlib.h"
#include <vector>

namespace SPARTA_NS {

template<class T>
class MyVec {
 public:
  int n;

  MyVec() {
    n = nmax = 0;
  }

  inline T& operator[](int index) {return vec[index];}

  void grow(int nnew) {
    if (nnew <= nmax) return;
    vec.resize(nnew);
    nmax = nnew;
  }

 private:
  std::vector<T> vec;
  int nmax;
};

}

#endif

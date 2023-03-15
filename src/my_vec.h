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

/* ----------------------------------------------------------------------
MyVec = templated class for a vector of objects with control of resizing
  really just a wrapper on STL vector class
  vector only grows, does not shrink, so can reuse w/out reallocs
  must grow() vector to desired length before storing elements
usage:
  grow() to set size of vec
  store/retrieve values from vec
  repeat, size of vec can vary each time
inputs:
  template T = one object, e.g. int, double, struct
methods:
  vec[i] = inlined access to vector element for read/write
  void grow(int nnew) = resize vector to length nnew if needed
public variables:
  n = current # of values in vector, set and reset by caller
------------------------------------------------------------------------- */


#ifndef SPARTA_MYVEC_H
#define SPARTA_MYVEC_H

#include "stdlib.h"
#include <vector>

namespace SPARTA_NS {

template<class T>
class MyVec {
 public:
  int n;        // current length of vector

  // construct an empty vector
  // no need for destructor, since STL vec is a class member

  MyVec() {
    n = nmax = 0;
  }

  // access an element of vector
  // caller must insure index is in range 0 to Nmax-1
  // inlined so performance should be same as STL vector access

  inline T& operator[](int index) {return vec[index];}

  // resize vector to nnew, only if necessary

  void grow(int nnew) {
    if (nnew <= nmax) return;
    vec.resize(nnew);
    nmax = nnew;
  }

 private:
  std::vector<T> vec;    // underlying STL vector
  int nmax;              // length of STL vector allocation
};

}

#endif

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

#ifndef SPARTA_KK_COPY_H
#define SPARTA_KK_COPY_H

#include "grid_kokkos.h"
#include "domain_kokkos.h"

// Need a copy of classes instantiated on the stack at the class level scope.
// However, this isn't directly possible due to issues with pointers.h
//  and Kokkos allocation tracking.
// This class is a workaround, using low-level memory operations.

namespace SPARTA_NS {

template <class ClassStyle>
class KKCopy {
 public:
  ClassStyle obj;

  KKCopy(SPARTA *sparta):
  obj(sparta) {
    ptr_temp = NULL;
    obj.copy = 1;
    save();
  }

  ~KKCopy() {}

  void copy(void* orig) {
    memcpy((void*)&obj, orig, sizeof(ClassStyle));
    obj.copy = 1;
  }

  void uncopy(int copy = 0) {
    if (ptr_temp != NULL) {
      memcpy((void*)&obj, ptr_temp, sizeof(ClassStyle));
      free(ptr_temp);
      ptr_temp = NULL;
    }
    obj.copy = copy;
    obj.copymode = 0;
  }

 private:
  void* ptr_temp;

  void save() {
    ptr_temp = (ClassStyle*) malloc(sizeof(ClassStyle));
    memcpy(ptr_temp, (void*)&obj, sizeof(ClassStyle));
  }

};

}

#endif

/* ERROR/WARNING messages:

*/

/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
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

#include <cstdlib>
#include <cstring>

#include "pointers.h"

// Need a copy of classes instantiated on the stack at the class level scope.
// However, this isn't directly possible due to issues with pointers.h
//  and Kokkos allocation tracking.
// This class is a workaround, using low-level memory operations.
//
// copy() blits the live object's bytes over obj.  That deliberately bypasses
//  Kokkos View reference counting, so obj ends up holding View handles whose
//  reference count was never incremented.  Such a handle must never reach a
//  View destructor: it would decrement a count it never incremented, which
//  aborts with "SharedAllocationRecord failed decrement count = 0" and frees
//  memory the original object still owns.
//
// The constructor therefore snapshots obj's pristine, correctly reference
//  counted bytes, and the destructor puts them back before obj is destroyed.
//  Restoring in the destructor is also what makes nesting work: a wrapped
//  class that itself contains KKCopy members has those members restored by
//  ordinary C++ destruction order.
//
// A copy constructed KKCopy -- which is what happens when an enclosing class
//  is captured by value into a Kokkos functor -- must be non-owning.  Its obj
//  is properly copy constructed, so its View reference counts are already
//  correct, and the snapshot belongs to the original: freeing it here would
//  leave the original restoring from freed memory.

namespace SPARTA_NS {

template <class ClassStyle>
class KKCopy {
 public:
  ClassStyle obj;

  KKCopy(SPARTA *sparta):
  obj(sparta) {
    ptr_temp = NULL;
    save();
    obj.copy = 1;
  }

  // a copy is non-owning: it restores nothing and frees nothing

  KKCopy(const KKCopy &other) : obj(other.obj) {
    ptr_temp = NULL;
  }

  KKCopy &operator=(const KKCopy &) = delete;

  ~KKCopy() { restore(); }

  void copy(const ClassStyle* orig) {
    if (ptr_temp == NULL) save();
    memcpy((void*)&obj, (const void*)orig, sizeof(ClassStyle));
    obj.copy = 1;
  }

 private:
  void* ptr_temp;

  void save() {
    ptr_temp = (ClassStyle*) malloc(sizeof(ClassStyle));
    memcpy(ptr_temp, (void*)&obj, sizeof(ClassStyle));
  }

  void restore() {
    if (ptr_temp == NULL) return;
    memcpy((void*)&obj, ptr_temp, sizeof(ClassStyle));
    free(ptr_temp);
    ptr_temp = NULL;
  }

};

}

#endif

/* ERROR/WARNING messages:

*/

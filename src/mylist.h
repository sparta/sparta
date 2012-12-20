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

#ifndef SPARTA_MYLIST_H
#define SPARTA_MYLIST_H

#include "stdlib.h"

namespace SPARTA_NS {

template<class T>
class MyList {
 public:
  int n;
  T first,last;

  MyList() {
    reset();
  }

  void reset() {
    n = 0;
    first = last = NULL;
  }

  void append(T one) {
    if (last) last->next = one;
    one->prev = last;
    one->next = NULL;
    if (!first) first = one;
    last = one;
    n++;
  }

  void prepend(T one) {
    if (first) first->prev = one;
    one->prev = NULL;
    one->next = first;
    if (!last) last = one;
    first = one;
    n++;
  }

  void insert(T one, T prev, T next) {
    if (prev) prev->next = one;
    else first = one;
    if (next) next->prev = one;
    else last = one;
    one->prev = prev;
    one->next = next;
    n++;
  }

  void remove(T one) {
    if (one->prev) one->prev->next = one->next;
    else first = one->next;
    if (one->next) one->next->prev = one->prev;
    else last = one->prev;
    n--;
  }
};

}

#endif

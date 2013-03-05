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

/* ----------------------------------------------------------------------
MyDoubleLinkedList = templated class for a doubly linked list of datums
  T should be a ptr to a Struct that has prev/next fields
  allocates/frees no memory, that is done by caller, so can reuse
  caller can iterate over list using first,last and prev/next fields
  NOTE: could provide an iterator so prev/next could be stored internally, like FIFO
usage:
  use MyPool to provide pool of objects
  build and use a doubly linked list of object ptrs
  reset()
  repeat, size of list can vary each time
inputs:
  template T = ptr to a Struct
methods:
  all methods are O(1) operations
  void reset() = reset list to be empty
  void append(T entry)  = append entry to end of list
  void prepend(T entry) = prepend entry to beginning of list
  void insert(T entry, T prev, T next) = insert entry between prev and next, 
    either/both of which can be NULL
  void remove(T entry) = remove entry from list
public variables:
  ndatum = total # of datums in list
  first,last = ptrs to first/last elements, NULL if empty
------------------------------------------------------------------------- */

#ifndef SPARTA_MYLIST_H
#define SPARTA_MYLIST_H

#include "stdlib.h"

namespace SPARTA_NS {

template<class T>
class MyDoubleLinkedList {
 public:
  int nlist;       // # of entries in list
  T first,last;    // ptr to first and last entry of list, NULL if empty

  MyDoubleLinkedList() {reset();}

  // reset the list to be empty

  void reset() {   
    first = last = NULL;
    nlist = 0;
  }

  // append entry to end of list

  void append(T entry) {
    if (last) last->next = entry;
    entry->prev = last;
    entry->next = NULL;
    if (!first) first = entry;
    last = entry;
    nlist++;
  }

  // prepend entry to beginning of list

  void prepend(T entry) {
    if (first) first->prev = entry;
    entry->prev = NULL;
    entry->next = first;
    if (!last) last = entry;
    first = entry;
    nlist++;
  }

  // insert entry between prev and next
  // prev and/or next can be NULL

  void insert(T entry, T prev, T next) {
    if (prev) prev->next = entry;
    else first = entry;
    if (next) next->prev = entry;
    else last = entry;
    entry->prev = prev;
    entry->next = next;
    nlist++;
  }

  // remove entry from list

  void remove(T entry) {
    if (entry->prev) entry->prev->next = entry->next;
    else first = entry->next;
    if (entry->next) entry->next->prev = entry->prev;
    else last = entry->prev;
    nlist--;
  }
};

}

#endif

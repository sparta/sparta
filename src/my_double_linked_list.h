/* ----------------------------------------------------------------------
MyDoubleLinkedList = templated class for a doubly linked list of datums
  datum must be a ptr to a Struct with prev/next fields
  allocates/frees no memory, that is done by caller, so can reuse
  caller can iterate over list using first,last and prev/next fields
  NOTE: could provide an iterator so prev/next could be stored internally, 
        like FIFO
usage:
  use MyPool to provide pool of objects
  insert/remove and use a doubly linked list of object ptrs
  reset(), repeat, size of list can vary each time
inputs:
  template T = ptr to a datum = Struct with prev/next fields
methods:
  all methods are O(1) operations
  void reset() = reset list to be empty
  void append(T entry)  = append entry to end of list
  void prepend(T entry) = prepend entry to beginning of list
  void insert(T entry, T prev, T next) = insert entry between prev and next, 
    either/both of which can be NULL
  void remove(T entry) = remove entry from list
public variables:
  nlist = total # of datums in list
  first,last = ptrs to first/last elements, NULL if empty
------------------------------------------------------------------------- */

#ifndef MY_DOUBLE_LINKED_LIST_H
#define MY_DOUBLE_LINKED_LIST_H

#include "stdlib.h"

template<class T>
class MyDoubleLinkedList {
 public:
  int nlist;       // # of datums in list
  T first,last;    // ptr to first and last datum in list

  MyDoubleLinkedList() {reset();}

  // reset the list to be empty

  void reset() {   
    first = last = NULL;
    nlist = 0;
  }

  // append new datum to end of list

  void append(T datum) {
    if (last) last->next = datum;
    datum->prev = last;
    datum->next = NULL;
    if (!first) first = datum;
    last = datum;
    nlist++;
  }

  // prepend new datum to beginning of list

  void prepend(T datum) {
    if (first) first->prev = datum;
    datum->prev = NULL;
    datum->next = first;
    if (!last) last = datum;
    first = datum;
    nlist++;
  }

  // insert new datum between prev and next
  // prev and/or next can be NULL

  void insert(T datum, T prev, T next) {
    if (prev) prev->next = datum;
    else first = datum;
    if (next) next->prev = datum;
    else last = datum;
    datum->prev = prev;
    datum->next = next;
    nlist++;
  }

  // move existing datum to front of list

  void move2front(T datum) {
    if (first == datum) return;
    datum->prev->next = datum->next;
    if (datum->next) datum->next->prev = datum->prev;
    else last = datum->prev;
    datum->prev = NULL;
    datum->next = first;
    first->prev = datum;
    first = datum;
  }

  // remove existing datum from list

  void remove(T datum) {
    if (datum->prev) datum->prev->next = datum->next;
    else first = datum->next;
    if (datum->next) datum->next->prev = datum->prev;
    else last = datum->prev;
    nlist--;
  }

  // debug check if list structure is consistent
  // walk in both directions, check all pointers

  int check() {
    if (first == NULL || last == NULL) {
      if (first || last || nlist) return 1;
      return 0;
    }

    int count = 0;
    T ptr = first;
    while (ptr) {
      count++;
      ptr = ptr->next;
    }
    if (count != nlist) return 2;

    count = 0;
    ptr = last;
    while (ptr) {
      count++;
      ptr = ptr->prev;
    }
    if (count != nlist) return 3;

    if (first->prev || last->next) return 4;

    ptr = first;
    while (ptr) {
      if (ptr != first && ptr->prev == NULL) return 5;
      if (ptr != last && ptr->next == NULL) return 6;
      ptr = ptr->next;
    }

    return 0;
  }
};

#endif

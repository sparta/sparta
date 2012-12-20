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

#ifndef SPARTA_POOL_H
#define SPARTA_POOL_H

#include "stdlib.h"

namespace SPARTA_NS {

// NOTE: set this larger eventually

#define PAGESIZE 4

template<class T>
class Pool {
 public:
  Pool() {
    npage = 0;
    pages = NULL;
    allocate();
    ipage = index = 0;
    page = pages[0];
  }

  ~Pool() {
    for (int i = 0; i < npage; i++) free(pages[i]);
    free(pages);
  }

  void reset() {
    ipage = index = 0;
    page = pages[0];
  }

  T *get() {
    if (index < PAGESIZE) return &page[index++];
    index = 0;
    ipage++;
    if (ipage == npage) allocate();
    page = pages[ipage];
    return &page[index++];
  }

 private:
  T **pages;    // list of allocated pages
  int npage;    // # of allocated pages
  int ipage;    // index of current page
  T *page;      // ptr to current page
  int index;    // current index on page

  void allocate() {
    pages = (T **) realloc(pages,(npage+1)*sizeof(T *));
    //if (!pages) error;
    pages[npage] = (T *) malloc(PAGESIZE*sizeof(T));
    //if (!pages[npage]) error;
    npage++;
  }
};

}

#endif

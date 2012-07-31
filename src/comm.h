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

#ifndef SPARTA_COMM_H
#define SPARTA_COMM_H

#include "pointers.h"

namespace SPARTA_NS {

class Comm : protected Pointers {
 public:
  int me,nprocs;                    // proc info

  bigint ncomm;                     // dummy statistic for now

  Comm(class SPARTA *);
  ~Comm();
  void init() {}
  void migrate(int, int *);

 private:
  class Irregular *irregular;
  char *sbuf,*rbuf;
  int maxsend;
  int *proclist;
};

}

#endif

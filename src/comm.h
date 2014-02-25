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

#ifndef SPARTA_COMM_H
#define SPARTA_COMM_H

#include "pointers.h"

namespace SPARTA_NS {

class Comm : protected Pointers {
 public:
  int me,nprocs;                    // proc info
  bigint ncomm;                     // dummy statistic for now

  int commsortflag;                 // 1 to force sort in all irregular comms
                                    //   useful for debugging to insure
                                    //   reproducible ordering of recv datums
  int commpartstyle;                // 1 for neighbor, 0 for all
                                    //   changes how irregular comm for
                                    //   particles is performed

  Comm(class SPARTA *);
  ~Comm();
  void init() {}
  void reset_neighbors();
  int migrate_particles(int, int *);
  void migrate_cells(int);
  void ring(int, int, void *, int, void (*)(int, char *), 
            void *, int self = 1);

 private:
  class Irregular *irregular;
  class Irregular *irregular_grid;
  char *sbuf,*rbuf;
  int maxsendbuf,maxrecvbuf;
  int *pproc,*gproc,*gsize;
  int maxpproc,maxgproc;
  
  int neighflag;                    // 1 if nearest-neighbor particle comm
  int nneigh;                       // # of procs I own ghost cells of
  int *neighlist;                   // list of ghost procs
};

}

#endif

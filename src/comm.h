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
  virtual void migrate_cells(int);
  int send_cells_adapt(int, int *, char *, char **);
  int irregular_uniform_neighs(int, int *, char *, int, char **);
  int irregular_uniform(int, int *, char *, int, char **);
  void ring(int, int, void *, int, void (*)(int, char *, void *),
            void *, int, void *);
  int rendezvous(int, int, char *, int, int, int *,
                 int (*)(int, char *, int &, int *&, char *&, void *),
                 int, char *&, int, void *, int statflag=0);

 protected:
  class Irregular *iparticle,*igrid,*iuniform;
  char *sbuf,*rbuf;
  int maxsendbuf,maxrecvbuf;
  int *pproc,*gproc,*gsize;
  int maxpproc,maxgproc;
  bigint rvous_bytes;

  int neighflag;                    // 1 if nearest-neighbor particle comm
  int nneigh;                       // # of procs I own ghost cells of
  int *neighlist;                   // list of ghost procs

  int copymode;                 // 1 if copy of class (prevents deallocation of
                                // base class when child copy is destroyed)

  void migrate_cells_less_memory(int);  // small memory version of migrate_cells
  int rendezvous_irregular(int, char *, int, int, int *,
                           int (*)(int, char *, int &, int *&, char *&, void *),
                           int, char *&, int, void *, int);
  int rendezvous_all2all(int, char *, int, int, int *,
                         int (*)(int, char *, int &, int *&, char *&, void *),
                         int, char *&, int, void *, int);
  void rendezvous_stats(int, int, int, int, int, int);
};

}

#endif

/* ERROR/WARNING messages:

E: Migrate cells send buffer exceeds 2 GB

MPI does not support a communication buffer that exceeds a 4-byte
integer in size.

*/

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

#ifndef SPARTA_RCB_H
#define SPARTA_RCB_H

#include "mpi.h"
#include "pointers.h"

namespace SPARTA_NS {

class RCB : protected Pointers {
 public:
  // set by compute()

  int noriginal;              // # of dots I own before balancing
  int nfinal;                 // # of dots I own after balancing
  int nkeep;                  // how many dots of noriginal I still own
                              // will be first nkept of nfinal list
  int *recvproc;              // proc IDs of nfinal dots
  int *recvindex;             // index of nfinal dots on owning procs
                              // based on input list for compute()
  double *lo,*hi;             // final bounding box of my RCB sub-domain

  // set by invert()

  int *sendproc;              // proc to send each of my noriginal dots to
  int *sendindex;             // index of dot in receiver's nfinal list

  RCB(class SPARTA *);
  ~RCB();
  void compute(int, double **, double *, char *, int flip=0);
  void invert();
  void check();
  void stats(int);

  // RCB cut info

  struct Median {
    double totallo,totalhi;   // weight in each half of active partition
    double valuelo,valuehi;   // position of dot(s) nearest to cut
    double wtlo,wthi;         // total weight of dot(s) at that position
    int countlo,counthi;      // # of dots at that position
    int proclo,prochi;              // unique proc who owns a nearest dot
  };

  // bounding box

  struct BBox {
    double lo[3];
    double hi[3];
  };

 private:
  int me,nprocs;

  // point to balance on

  struct Dot {
    double x[3];          // coord of point
    double wt;            // weight of point
    int proc;             // owning proc
    int index;            // index on owning proc
  };

  // tree of RCB cuts

  struct Tree {
    double cut;                // position of cut
    int dim;                // dimension = 0/1/2 of cut
  };

  // inversion message

  struct Invert {
    int rindex;         // index on receiving proc
    int sproc;          // sending proc
    int sindex;         // index on sending proc
  };

  Dot *dots;        // dots on this proc
  int ndot;         // # of dots on this proc
  int maxdot;       // allocated size of dots
  int ndotorig;

  int nlist;
  int maxlist;
  int *dotlist;
  int *dotmark;

  int maxbuf;
  Dot *buf;

  int maxrecv,maxsend;

  class Irregular *irregular;

  MPI_Op box_op,med_op;
  MPI_Datatype box_type,med_type;

  int reuse;        // 1/0 to use/not use previous cuts
  int dottop;       // dots >= this index are new
  BBox rcbbox;      // bounding box of final RCB sub-domain
  Tree *tree;       // tree of RCB cuts, used by reuse()
  int counters[7];  // diagnostic counts
                    // 0 = # of median iterations
                    // 1 = # of points sent
                    // 2 = # of points received
                    // 3 = most points this proc ever owns
                    // 4 = most point memory this proc ever allocs
                    // 5 = # of times a previous cut is re-used
                    // 6 = # of reallocs of point vector
};

}

#endif

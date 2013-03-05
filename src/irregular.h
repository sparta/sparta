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

#ifndef SPARTA_IRREGULAR_H
#define SPARTA_IRREGULAR_H

#include "pointers.h"

namespace SPARTA_NS {

class Irregular : protected Pointers {
 public:
  Irregular(class SPARTA *);
  ~Irregular();
  int create(int, int *, int sort = 0);
  void exchange(char *, int, char *);

 private:
  int me,nprocs;

  // plan for irregular communication of datums

  int nsend;                 // # of messages to send
  int nrecv;                 // # of messages to recv
  int sendmax;               // # of datums in largest send message
  int nrecvdatum;            // total # of datums I recv
  int num_self;              // # of datums to copy to self
  int indexmax;              // current size of index_send
  int indexselfmax;          // current size of index_self
  int bufmax;                // current size of buf in bytes
  int *proc_send;            // procs to send to
  int *num_send;             // # of datums to send to each proc
  int *proc_recv;            // procs to recv from
  int *num_recv;             // # of datums to recv from each proc
  int *index_send;           // list of which datums to send to each proc
  int *index_self;           // list of which datums to copy to self
  int *wlist;                // work vector
  int *wcount;               // work vector
  MPI_Request *request;      // MPI requests for posted recvs
  MPI_Status *status;        // MPI statuses for WaitAll
  char *buf;                 // buffer for largest single send message
};

}

#endif

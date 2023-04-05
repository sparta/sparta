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

#ifndef SPARTA_IRREGULAR_H
#define SPARTA_IRREGULAR_H

#include "pointers.h"

namespace SPARTA_NS {

class Irregular : protected Pointers {
 public:

  // static variable across all Irregular objects, for qsort callback

  static int *proc_recv_copy;

  Irregular(class SPARTA *);
  ~Irregular();
  void create_procs(int, int *, int sort = 0);
  virtual int create_data_uniform(int, int *, int sort = 0);
  int create_data_uniform_grouped(int, int *, int sort = 0);
  int create_data_variable(int, int *, int *, int &, int sort = 0);
  virtual int augment_data_uniform(int, int *);
  void exchange_uniform(char *, int, char *);
  void exchange_variable(char *, int *, char *);
  void reverse(int, int *);

 protected:
  int me,nprocs;

  // plan for irregular communication of datums
  // same for uniform or variable sized datums

  int nsend;                 // # of messages to send, no self
  int nrecv;                 // # of messages to recv, no self
  int sendmax;               // # of datums in largest send message
  int nrecvdatum;            // total # of datums I recv
  int num_self;              // # of datums to copy to self
  int indexmax;              // current size of index_send
  int indexselfmax;          // current size of index_self
  int bufmax;                // current size of buf in bytes
  int *proc_send;            // list of procs to send to
  int *num_send;             // # of datums to send to each proc
  int *proc_recv;            // list of procs to recv from
  int *num_recv;             // # of datums to recv from each proc
  int *index_send;           // list of datum indices to send to each proc
  int *index_self;           // list of datum indices to copy to self
  int *proc2recv;            // mapping from proc IDs to recv list
  int *work1,*work2;         // work vectors
  MPI_Request *request;      // MPI requests for posted recvs
  MPI_Status *status;        // MPI statuses for WaitAll
  char *buf;                 // buffer for largest single send message
  int copymode;              // 1 if copy of class (prevents deallocation of
                             //   base class when child copy is destroyed)

  // only defined for variable sized datums

  int sendmaxbytes;          // # of bytes in largest send message
  int size_self;             // # of bytes in datums copied to self
  int offsetmax;             // current size of offset_send
  int *size_send;            // # of bytes of send to each proc
  int *size_recv;            // # of bytes to recv from each proc
  int *offset_send;          // list of byte offsets for each send datum
};

}

#endif

/* ERROR/WARNING messages:

E: Irregular comm recv buffer exceeds 2 GB

MPI does not support a communication buffer that exceeds a 4-byte
integer in size.

*/

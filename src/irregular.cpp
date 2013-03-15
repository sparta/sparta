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

#include "spatype.h"
#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "irregular.h"
#include "particle.h"
#include "domain.h"
#include "comm.h"
#include "memory.h"

using namespace SPARTA_NS;

// allocate space for static class variable
// prototype for non-class function

int *Irregular::proc_recv_copy;
int compare_standalone(const void *, const void *);

#define BUFFACTOR 1.5
#define BUFMIN 1000
#define BUFEXTRA 1000

/* ---------------------------------------------------------------------- */

Irregular::Irregular(SPARTA *sparta) : Pointers(sparta)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  // allocate fixed-length and work vectors for plan

  proc_send = new int[nprocs];
  num_send = new int[nprocs];
  proc_recv = new int[nprocs];
  num_recv = new int[nprocs];
  request = new MPI_Request[nprocs];
  status = new MPI_Status[nprocs];

  wlist = new int[nprocs];
  wcount = new int[nprocs];

  indexmax = 0;
  index_send = NULL;
  indexselfmax = 0;
  index_self = NULL;
  bufmax = 0;
  buf = NULL;
}

/* ---------------------------------------------------------------------- */

Irregular::~Irregular()
{
  delete [] proc_send;
  delete [] num_send;
  delete [] proc_recv;
  delete [] num_recv;
  delete [] index_send;
  delete [] index_self;
  delete [] wlist;
  delete [] wcount;
  delete [] request;
  delete [] status;
  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   create a communication plan
   n = # of datums to send
   proclist = proc to send each datum to (including self)
   return total # of datums I will recv (including self)
   sort = 1 if receives should be ordered by proc, default = 0
------------------------------------------------------------------------- */

int Irregular::create(int n, int *proclist, int sort)
{
  int i,m;

  // nrecv = # of messages I receive

  for (i = 0; i < nprocs; i++) {
    wlist[i] = 0;
    wcount[i] = 1;
  }
  for (i = 0; i < n; i++) wlist[proclist[i]] = 1;

  MPI_Reduce_scatter(wlist,&nrecv,wcount,MPI_INT,MPI_SUM,world);
  if (wlist[me]) nrecv--;

  // nsend = # of messages I send

  for (i = 0; i < nprocs; i++) wlist[i] = 0;
  for (i = 0; i < n; i++) wlist[proclist[i]]++;

  nsend = 0;
  for (i = 0; i < nprocs; i++)
    if (wlist[i]) nsend++;
  if (wlist[me]) nsend--;

  // reallocate send and self index lists if necessary
  // could use n-wlist[m] for length of index_send to be more precise

  if (n > indexmax) {
    indexmax = n;
    delete [] index_send;
    index_send = new int[n];
  }

  if (wlist[me] > indexselfmax) {
    indexselfmax = wlist[me];
    delete [] index_self;
    index_self = new int[wlist[me]];
  }

  // proc_send = procs I send to
  // num_send = # of datums I send to each proc
  // num_self = # of datums I copy to self
  // to balance pattern of send messages:
  //   each proc starts with iproc > me, continues until iproc = me
  // reset wlist to store which send message each proc corresponds to

  int iproc = me;
  int isend = 0;
  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    if (iproc == me) num_self = wlist[iproc];
    else if (wlist[iproc] > 0) {
      proc_send[isend] = iproc;
      num_send[isend] = wlist[iproc];
      wlist[iproc] = isend;
      isend++;
    }
  }
  wlist[me] = 0;

  // wcount = offsets into index_send for each proc I send to
  // m = ptr into index_self
  // index_send = list of which datums to send to each proc
  //   1st N1 values are datum indices for 1st proc,
  //   next N2 values are datum indices for 2nd proc, etc

  wcount[0] = 0;
  for (i = 1; i < nsend; i++) wcount[i] = wcount[i-1] + num_send[i-1];

  m = 0;
  for (i = 0; i < n; i++) {
    iproc = proclist[i];
    if (iproc == me) index_self[m++] = i;
    else {
      isend = wlist[iproc];
      index_send[wcount[isend]++] = i;
    }
  }

  // tell receivers how much data I send
  // sendmax = largest # of datums I send in a single message

  sendmax = 0;
  for (i = 0; i < nsend; i++) {
    MPI_Send(&num_send[i],1,MPI_INT,proc_send[i],0,world);
    sendmax = MAX(sendmax,num_send[i]);
  }

  // receive incoming messages
  // proc_recv = procs I recv from
  // num_recv = # of datums each proc sends me
  // nrecvdatum = total # of datums I recv

  nrecvdatum = 0;
  for (i = 0; i < nrecv; i++) {
    MPI_Recv(&num_recv[i],1,MPI_INT,MPI_ANY_SOURCE,0,world,status);
    proc_recv[i] = status->MPI_SOURCE;
    nrecvdatum += num_recv[i];
  }
  nrecvdatum += num_self;

  // sort proc_recv and num_recv by proc ID if requested
  // useful for debugging to insure reproducible behavior

  if (sort) {
    int *order = new int[nrecv];
    int *proc_recv_ordered = new int[nrecv];
    int *num_recv_ordered = new int[nrecv];

    for (i = 0; i < nrecv; i++) order[i] = i;
    proc_recv_copy = proc_recv;
    qsort(order,nrecv,sizeof(int),compare_standalone);

    int j;
    for (i = 0; i < nrecv; i++) {
      j = order[i];
      proc_recv_ordered[i] = proc_recv[j];
      num_recv_ordered[i] = num_recv[j];
    }

    memcpy(proc_recv,proc_recv_ordered,nrecv*sizeof(int));
    memcpy(num_recv,num_recv_ordered,nrecv*sizeof(int));
    delete [] order;
    delete [] proc_recv_ordered;
    delete [] num_recv_ordered;
  }

  // barrier to insure all MPI_ANY_SOURCE messages are received
  // else another proc could proceed to exchange_data() and send to me

  MPI_Barrier(world);

  // return # of datums I will receive

  return nrecvdatum;
}

/* ----------------------------------------------------------------------
   comparison function invoked by qsort()
   accesses static class member proc_recv_copy, set before call to qsort()
------------------------------------------------------------------------- */

int compare_standalone(const void *iptr, const void *jptr)
{
  int i = *((int *) iptr);
  int j = *((int *) jptr);
  int *proc_recv = Irregular::proc_recv_copy;
  if (proc_recv[i] < proc_recv[j]) return -1;
  if (proc_recv[i] > proc_recv[j]) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   communicate datums via Plan
   sendbuf = list of datums to send
   nbytes = size of each datum
   recvbuf = received datums (including copied from me)
------------------------------------------------------------------------- */

void Irregular::exchange(char *sendbuf, int nbytes, char *recvbuf)
{
  int i,m,n,offset,count;

  // post all receives, starting after self copies

  offset = num_self*nbytes;
  for (int irecv = 0; irecv < nrecv; irecv++) {
    MPI_Irecv(&recvbuf[offset],num_recv[irecv]*nbytes,MPI_CHAR,
	      proc_recv[irecv],0,world,&request[irecv]);
    offset += num_recv[irecv]*nbytes;
  }

  // reallocate buf for largest send if necessary

  if (sendmax*nbytes > bufmax) {
    memory->destroy(buf);
    bufmax = sendmax*nbytes;
    memory->create(buf,bufmax,"irregular:buf");
  }

  // send each message
  // pack buf with list of datums
  // m = index of datum in sendbuf

  n = 0;
  for (int isend = 0; isend < nsend; isend++) {
    count = num_send[isend];
    for (i = 0; i < count; i++) {
      m = index_send[n++];
      memcpy(&buf[i*nbytes],&sendbuf[m*nbytes],nbytes);
    }
    MPI_Send(buf,count*nbytes,MPI_CHAR,proc_send[isend],0,world);
  }       

  // copy datums to self, put at beginning of recvbuf

  for (i = 0; i < num_self; i++) {
    m = index_self[i];
    memcpy(&recvbuf[i*nbytes],&sendbuf[m*nbytes],nbytes);
  }

  // wait on all incoming messages

  if (nrecv) MPI_Waitall(nrecv,request,status);
}

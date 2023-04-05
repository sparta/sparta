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

#include "spatype.h"
#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "irregular_kokkos.h"
#include "particle.h"
#include "domain.h"
#include "comm.h"
#include "memory_kokkos.h"
#include "error.h"
#include "kokkos.h"



// DEBUG
#include "update.h"
#include "grid.h"

using namespace SPARTA_NS;

// allocate space for static class variable
// prototype for non-class function

//int *IrregularKokkos::proc_recv_copy;
int compare_standalone(const void *, const void *);

#define BUFFACTOR 1.5
#define BUFMIN 1000
#define BUFEXTRA 1000

/* ---------------------------------------------------------------------- */

IrregularKokkos::IrregularKokkos(SPARTA *sparta) : Irregular(sparta)
{
  k_n = DAT::tdual_int_scalar("comm:n");
}

/* ---------------------------------------------------------------------- */

IrregularKokkos::~IrregularKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_index_send,index_send);
  index_send = NULL;

  memoryKK->destroy_kokkos(k_index_self,index_self);
  index_self = NULL;
}

/* ----------------------------------------------------------------------
   create communication plan based on list of datums of uniform size
   n = # of datums to send
   proclist = proc to send each datum to, can include self
   sort = flag for sorting order of received datums by proc ID
   return total # of datums I will recv, including any to self
------------------------------------------------------------------------- */

int IrregularKokkos::create_data_uniform(int n, int *proclist, int sort)
{
  int i,m;

  // setup for collective comm
  // work1 = # of datums I send to each proc, set self to 0
  // work2 = 1 for all procs, used for ReduceScatter

  for (i = 0; i < nprocs; i++) {
    work1[i] = 0;
    work2[i] = 1;
  }
  for (i = 0; i < n; i++) work1[proclist[i]] = 1;
  work1[me] = 0;

  // nrecv = # of procs I receive messages from, not including self
  // options for performing ReduceScatter operation
  // some are more efficient on some machines at big sizes

#ifdef SPARTA_RS_ALLREDUCE_INPLACE
  MPI_Allreduce(MPI_IN_PLACE,work1,nprocs,MPI_INT,MPI_SUM,world);
  nrecv = work1[me];
#else
#ifdef SPARTA_RS_ALLREDUCE
  MPI_Allreduce(work1,work2,nprocs,MPI_INT,MPI_SUM,world);
  nrecv = work2[me];
#else
  MPI_Reduce_scatter(work1,&nrecv,work2,MPI_INT,MPI_SUM,world);
#endif
#endif

  // work1 = # of datums I send to each proc, including self
  // nsend = # of procs I send messages to, not including self

  for (i = 0; i < nprocs; i++) work1[i] = 0;
  for (i = 0; i < n; i++) work1[proclist[i]]++;

  nsend = 0;
  for (i = 0; i < nprocs; i++)
    if (work1[i]) nsend++;
  if (work1[me]) nsend--;

  // reallocate send and self index lists if necessary
  // could use n-work1[me] for length of index_send to be more precise

  if (n > indexmax) {
    indexmax = n;
    memoryKK->destroy_kokkos(k_index_send,index_send);
    memoryKK->create_kokkos(k_index_send,index_send,indexmax,"irregular:index_send");
    d_index_send = k_index_send.d_view;
  }

  if (work1[me] > indexselfmax) {
    indexselfmax = work1[me];
    memoryKK->destroy_kokkos(k_index_self,index_self);
    memoryKK->create_kokkos(k_index_self,index_self,indexselfmax,"irregular:index_self");
    d_index_self = k_index_self.d_view;
  }

  // proc_send = procs I send to
  // num_send = # of datums I send to each proc
  // num_self = # of datums I copy to self
  // to balance pattern of send messages:
  //   each proc starts with iproc > me, continues until iproc = me
  // reset work1 to store which send message each proc corresponds to

  int iproc = me;
  int isend = 0;
  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    if (iproc == me) {
      num_self = work1[iproc];
      work1[iproc] = 0;
    } else if (work1[iproc]) {
      proc_send[isend] = iproc;
      num_send[isend] = work1[iproc];
      work1[iproc] = isend;
      isend++;
    }
  }

  // work2 = offsets into index_send for each proc I send to
  // m = ptr into index_self
  // index_send = list of which datums to send to each proc
  //   1st N1 values are datum indices for 1st proc,
  //   next N2 values are datum indices for 2nd proc, etc
  // index_self = list of which datums to copy to self

  work2[0] = 0;
  for (i = 1; i < nsend; i++) work2[i] = work2[i-1] + num_send[i-1];

  m = 0;
  for (i = 0; i < n; i++) {
    iproc = proclist[i];
    if (iproc == me) index_self[m++] = i;
    else {
      isend = work1[iproc];
      index_send[work2[isend]++] = i;
    }
  }

  // tell receivers how many datums I send them
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
  // useful for debugging to insure reproducible ordering of received datums

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

  // proc2recv[I] = which recv the Ith proc ID is
  // will only be accessed by procs I actually receive from

  for (i = 0; i < nrecv; i++) proc2recv[proc_recv[i]] = i;

  // barrier to insure all MPI_ANY_SOURCE messages are received
  // else another proc could proceed to exchange_data() and send to me

  MPI_Barrier(world);

  // return # of datums I will receive

  return nrecvdatum;
}

/* ----------------------------------------------------------------------
   augment communication plan with new datums of uniform size
   called after create_procs() created initial plan
   n = # of datums to send
   proclist = proc to send each datum to, can include self
   return total # of datums I will recv
------------------------------------------------------------------------- */

int IrregularKokkos::augment_data_uniform(int n, int *proclist)
{
  int i,m,iproc,isend;

  // tally count of messages to each proc in num_send and num_self

  num_self = 0;
  for (i = 0; i < nsend; i++) work2[proc_send[i]] = 0;
  work2[me] = 0;
  for (i = 0; i < n; i++) work2[proclist[i]]++;
  for (i = 0; i < nsend; i++) num_send[i] = work2[proc_send[i]];
  num_self = work2[me];

  // reallocate send and self index lists if necessary
  // could use n-num_self for length of index_send to be more precise

  if (n > indexmax) {
    indexmax = n;
    memoryKK->destroy_kokkos(k_index_send,index_send);
    memoryKK->create_kokkos(k_index_send,index_send,indexmax,"irregular:index_send");
    d_index_send = k_index_send.d_view;
  }

  if (num_self > indexselfmax) {
    indexselfmax = num_self;
    memoryKK->destroy_kokkos(k_index_self,index_self);
    memoryKK->create_kokkos(k_index_self,index_self,indexselfmax,"irregular:index_self");
    d_index_self = k_index_self.d_view;
  }

  // work2 = offsets into index_send for each proc I send to
  // m = ptr into index_self
  // index_send = list of which datums to send to each proc
  //   1st N1 values are datum indices for 1st proc,
  //   next N2 values are datum indices for 2nd proc, etc
  // index_self = list of which datums to copy to self

  work2[0] = 0;
  for (i = 1; i < nsend; i++) work2[i] = work2[i-1] + num_send[i-1];

  if (num_self) {
    m = 0;
    for (i = 0; i < n; i++) {
      iproc = proclist[i];
      if (iproc == me) index_self[m++] = i;
      else {
        isend = work1[iproc];
        index_send[work2[isend]++] = i;
      }
    }
  } else {
    for (i = 0; i < n; i++) {
      isend = work1[proclist[i]];
      index_send[work2[isend]++] = i;
    }
  }

  // tell receivers how many datums I send them
  // sendmax = largest # of datums I send in a single message

  sendmax = 0;
  for (i = 0; i < nsend; i++) {
    MPI_Send(&num_send[i],1,MPI_INT,proc_send[i],0,world);
    sendmax = MAX(sendmax,num_send[i]);
  }

  // receive incoming messages
  // num_recv = # of datums each proc sends me
  // nrecvdatum = total # of datums I recv

  nrecvdatum = 0;
  for (i = 0; i < nrecv; i++) {
    MPI_Recv(&m,1,MPI_INT,MPI_ANY_SOURCE,0,world,status);
    iproc = status->MPI_SOURCE;
    num_recv[proc2recv[iproc]] = m;
    nrecvdatum += m;
  }
  nrecvdatum += num_self;

  // barrier to insure all MPI_ANY_SOURCE messages are received
  // else another proc could proceed to exchange_data() and send to me

  MPI_Barrier(world);

  // return # of datums I will receive

  return nrecvdatum;
}

/* ----------------------------------------------------------------------
   communicate uniform-size datums via existing plan
   sendbuf = list of datums to send
   nbytes = size of each datum
   recvbuf = received datums, including copied from me
------------------------------------------------------------------------- */

void IrregularKokkos::exchange_uniform(DAT::t_char_1d d_sendbuf_in, int nbytes_in,
                                       char* d_recvbuf_ptr_in, DAT::t_char_1d d_recvbuf_in)
{
  int offset,count;

  nbytes = nbytes_in;
  d_sendbuf = d_sendbuf_in;
  d_recvbuf_ptr = d_recvbuf_ptr_in;
  d_recvbuf = d_recvbuf_in;

  if (!sparta->kokkos->gpu_aware_flag &&
      h_recvbuf.extent(0) != d_recvbuf.extent(0)) {
    h_recvbuf = HAT::t_char_1d(Kokkos::view_alloc("irregular:d_recvbuf:mirror",Kokkos::WithoutInitializing),d_recvbuf.extent(0));
  }

  // post all receives, starting after self copies

  offset = num_self*nbytes;
  for (int irecv = 0; irecv < nrecv; irecv++) {
    if (sparta->kokkos->gpu_aware_flag) {
      MPI_Irecv(&d_recvbuf_ptr[offset],num_recv[irecv]*nbytes,MPI_CHAR,
                proc_recv[irecv],0,world,&request[irecv]);
    } else {
      MPI_Irecv(h_recvbuf.data() + offset,num_recv[irecv]*nbytes,MPI_CHAR,
                proc_recv[irecv],0,world,&request[irecv]);
    }
    offset += num_recv[irecv]*nbytes;
  }

  // reallocate buf for largest send if necessary

  if (sendmax*nbytes > bufmax) {
    bufmax = sendmax*nbytes;
    d_buf = DAT::t_char_1d("Irregular:buf",bufmax);
  } else if (d_buf.extent(0) < bufmax) {
    d_buf = DAT::t_char_1d("Irregular:buf",bufmax);
  }

  // send each message
  // pack buf with list of datums
  // m = index of datum in sendbuf

  k_index_send.modify_host();
  k_index_send.sync_device();

  k_index_self.modify_host();
  k_index_self.sync_device();

  k_n.h_view() = 0;
  k_n.modify_host();
  k_n.sync_device();
  d_n = k_n.d_view;

  for (int isend = 0; isend < nsend; isend++) {
    count = num_send[isend];

    if (!sparta->kokkos->gpu_aware_flag) {

      // allocate exact buffer size to reduce GPU <--> CPU memory transfer

      d_buf = DAT::t_char_1d(Kokkos::view_alloc("irregular:buf",Kokkos::WithoutInitializing),count*nbytes);
      h_buf = HAT::t_char_1d(Kokkos::view_alloc("irregular:buf:mirror",Kokkos::WithoutInitializing),count*nbytes);
    }

    copymode = 1;
    if (sparta->kokkos->need_atomics)
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagIrregularPackBuffer<1> >(0,count),*this);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagIrregularPackBuffer<0> >(0,count),*this);
    DeviceType().fence();
    //pack_buffer_serial(0,count);
    copymode = 0;

    if (sparta->kokkos->gpu_aware_flag)
      MPI_Send(d_buf.data(),count*nbytes,MPI_CHAR,proc_send[isend],0,world);
    else {
      Kokkos::deep_copy(h_buf,d_buf);
      MPI_Send(h_buf.data(),count*nbytes,MPI_CHAR,proc_send[isend],0,world);
    }
  }

  // copy datums to self, put at beginning of recvbuf

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagIrregularUnpackBuffer>(0,num_self),*this);
  DeviceType().fence();
  copymode = 0;

  // wait on all incoming messages

  if (nrecv) {
    MPI_Waitall(nrecv,request,status);

    if (!sparta->kokkos->gpu_aware_flag)
      Kokkos::deep_copy(d_recvbuf,h_recvbuf);
  }
}

template<int NEED_ATOMICS>
KOKKOS_INLINE_FUNCTION
void IrregularKokkos::operator()(TagIrregularPackBuffer<NEED_ATOMICS>, const int &i) const {
  int n;
  if (NEED_ATOMICS)
    n = Kokkos::atomic_fetch_add(&d_n(),1);
  else {
    n = d_n();
    d_n()++;
  }
  const int m = d_index_send[n];
  memcpy(&d_buf[i*nbytes],&d_sendbuf[m*nbytes],nbytes);
}

KOKKOS_INLINE_FUNCTION
void IrregularKokkos::operator()(TagIrregularUnpackBuffer, const int &i) const {
  const int m = d_index_self[i];
  memcpy(&d_recvbuf[i*nbytes],&d_sendbuf[m*nbytes],nbytes);
}

inline
void IrregularKokkos::pack_buffer_serial(const int start, const int end) const {
  int n = 0;
  for (int i = start; i < end; i++) {
    const int m = d_index_send[n++];
    memcpy(&d_buf[i*nbytes],&d_sendbuf[m*nbytes],nbytes);
  }
}

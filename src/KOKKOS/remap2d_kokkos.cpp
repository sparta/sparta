/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "remap2d_kokkos.h"

#include "error.h"
#include "pack2d_kokkos.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
RemapKokkos2d<DeviceType>::RemapKokkos2d(SPARTA *sparta) : Pointers(sparta)
{
  plan = nullptr;
}

template<class DeviceType>
RemapKokkos2d<DeviceType>::RemapKokkos2d(SPARTA *sparta, MPI_Comm comm,
             int in_ilo, int in_ihi, int in_jlo, int in_jhi,
             int out_ilo, int out_ihi, int out_jlo, int out_jhi,
             int nqty, int permute, int memory,
             int precision, int usecollective,
             int usegpu_aware) : Pointers(sparta)
{
  plan = remap_2d_create_plan_kokkos(comm,
                              in_ilo,in_ihi,in_jlo,in_jhi,
                              out_ilo,out_ihi,out_jlo,out_jhi,
                              nqty,permute,memory,precision,usecollective,
                              usegpu_aware);
  if (plan == nullptr) error->one(FLERR,"Could not create 2d remap plan");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
RemapKokkos2d<DeviceType>::~RemapKokkos2d()
{
  remap_2d_destroy_plan_kokkos(plan);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void RemapKokkos2d<DeviceType>::perform(typename FFT_AT::t_FFT_SCALAR_1d d_in, typename FFT_AT::t_FFT_SCALAR_1d d_out, typename FFT_AT::t_FFT_SCALAR_1d d_buf)
{
  remap_2d_kokkos(d_in,d_out,d_buf,plan);
}


/* ----------------------------------------------------------------------
   Data layout for 2d remaps:

   data set of Nfast x Nslow elements is owned by P procs
   each element = nqty contiguous datums
   on input, each proc owns a subsection of the elements
   on output, each proc will own a (presumably different) subsection
   my subsection must not overlap with any other proc's subsection,
     i.e. the union of all proc's input (or output) subsections must
     exactly tile the global Nfast x Nslow data set
   when called from C, all subsection indices are
     C-style from 0 to N-1 where N = Nfast or Nslow
   when called from F77, all subsection indices are
     F77-style from 1 to N where N = Nfast or Nslow
   a proc can own 0 elements on input or output
     by specifying hi index < lo index
   on both input and output, data is stored contiguously on a processor
     with a fast-varying and slow-varying index
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Perform 2d remap

   Arguments:
   in           starting address of input data on this proc
   out          starting address of where output data for this proc
                  will be placed (can be same as in)
   buf          extra memory required for remap
                if memory=0 was used in call to remap_2d_create_plan
                  then buf must be big enough to hold output result
                  i.e. nqty * (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1)
                if memory=1 was used in call to remap_2d_create_plan
                  then buf is not used, can just be a dummy pointer
   plan         plan returned by previous call to remap_2d_create_plan
------------------------------------------------------------------------- */

template<class DeviceType>
void RemapKokkos2d<DeviceType>::remap_2d_kokkos(typename FFT_AT::t_FFT_SCALAR_1d d_in, typename FFT_AT::t_FFT_SCALAR_1d d_out, typename FFT_AT::t_FFT_SCALAR_1d d_buf,
              struct remap_plan_2d_kokkos<DeviceType> *plan)
{
  // collective flag not yet supported

  // use point-to-point communication

  int i,isend,irecv;
  typename FFT_AT::t_FFT_SCALAR_1d d_scratch;

  if (plan->memory == 0)
    d_scratch = d_buf;
  else
    d_scratch = plan->d_scratch;

  // post all recvs into scratch space

  FFT_SCALAR* v_scratch = d_scratch.data();
  if (!plan->usegpu_aware) {
    plan->h_scratch = Kokkos::create_mirror_view(d_scratch);
    v_scratch = plan->h_scratch.data();
  }

  for (irecv = 0; irecv < plan->nrecv; irecv++) {
    FFT_SCALAR* scratch = v_scratch + plan->recv_bufloc[irecv];
    MPI_Irecv(scratch,plan->recv_size[irecv],
              MPI_FFT_SCALAR,plan->recv_proc[irecv],0,
              plan->comm,&plan->request[irecv]);
  }

  FFT_SCALAR* v_sendbuf = plan->d_sendbuf.data();
  if (!plan->usegpu_aware) {
    plan->h_sendbuf = Kokkos::create_mirror_view(plan->d_sendbuf);
    v_sendbuf = plan->h_sendbuf.data();
  }

  // send all messages to other procs

  for (isend = 0; isend < plan->nsend; isend++) {
    int in_offset = plan->send_offset[isend];
    plan->pack(d_in,in_offset,
               plan->d_sendbuf,0,&plan->packplan[isend]);

    if (!plan->usegpu_aware)
      Kokkos::deep_copy(plan->h_sendbuf,plan->d_sendbuf);

    MPI_Send(v_sendbuf,plan->send_size[isend],MPI_FFT_SCALAR,
             plan->send_proc[isend],0,plan->comm);
  }

  // copy in -> scratch -> out for self data

  if (plan->self) {
    isend = plan->nsend;
    irecv = plan->nrecv;

    int in_offset = plan->send_offset[isend];
    int scratch_offset = plan->recv_bufloc[irecv];
    int out_offset = plan->recv_offset[irecv];

    plan->pack(d_in,in_offset,
               d_scratch,scratch_offset,
               &plan->packplan[isend]);
    plan->unpack(d_scratch,scratch_offset,
                 d_out,out_offset,&plan->unpackplan[irecv]);
  }

  // unpack all messages from scratch -> out

  for (i = 0; i < plan->nrecv; i++) {
    MPI_Waitany(plan->nrecv,plan->request,&irecv,MPI_STATUS_IGNORE);

    int scratch_offset = plan->recv_bufloc[irecv];
    int out_offset = plan->recv_offset[irecv];

    if (!plan->usegpu_aware)
      Kokkos::deep_copy(d_scratch,plan->h_scratch);

    plan->unpack(d_scratch,scratch_offset,
                 d_out,out_offset,&plan->unpackplan[irecv]);
  }
}

/* ----------------------------------------------------------------------
   Create plan for performing a 2d remap

   Arguments:
   comm                 MPI communicator for the P procs which own the data
   in_ilo,in_ihi        input bounds of data I own in fast index
   in_jlo,in_jhi        input bounds of data I own in slow index
   out_ilo,out_ihi      output bounds of data I own in fast index
   out_jlo,out_jhi      output bounds of data I own in slow index
   nqty                 # of datums per element
   permute              permutation in storage order of indices on output
                          0 = no permutation
                          1 = permute = slow->fast, fast-.slow
   memory               user provides buffer memory for remap or system does
                          0 = user provides memory
                          1 = system provides memory
   precision            precision of data
                          1 = single precision (4 bytes per datum)
                          2 = double precision (8 bytes per datum)
   usecollective        whether to use collective MPI or point-to-point
   usegpu_aware         whether to use GPU-Aware MPI or not
------------------------------------------------------------------------- */

template<class DeviceType>
struct remap_plan_2d_kokkos<DeviceType>* RemapKokkos2d<DeviceType>::remap_2d_create_plan_kokkos(
  MPI_Comm comm,
  int in_ilo, int in_ihi, int in_jlo, int in_jhi,
  int out_ilo, int out_ihi, int out_jlo, int out_jhi,
  int nqty, int permute, int memory, int /*precision*/,
  int usecollective, int usegpu_aware)
{

  struct remap_plan_2d_kokkos<DeviceType> *plan;
  struct extent_2d *inarray, *outarray;
  struct extent_2d in,out,overlap;
  int i,iproc,nsend,nrecv,ibuf,size,me,nprocs;

  // query MPI info

  MPI_Comm_rank(comm,&me);
  MPI_Comm_size(comm,&nprocs);

  // allocate memory for plan data struct

  plan = new struct remap_plan_2d_kokkos<DeviceType>;
  if (plan == nullptr) return nullptr;
  plan->usecollective = usecollective;
  plan->usegpu_aware = usegpu_aware;

  // store parameters in local data structs

  in.ilo = in_ilo;
  in.ihi = in_ihi;
  in.isize = in.ihi - in.ilo + 1;

  in.jlo = in_jlo;
  in.jhi = in_jhi;
  in.jsize = in.jhi - in.jlo + 1;

  out.ilo = out_ilo;
  out.ihi = out_ihi;
  out.isize = out.ihi - out.ilo + 1;

  out.jlo = out_jlo;
  out.jhi = out_jhi;
  out.jsize = out.jhi - out.jlo + 1;

  // combine output extents across all procs

  inarray = (struct extent_2d *) malloc(nprocs*sizeof(struct extent_2d));
  if (inarray == nullptr) return nullptr;

  outarray = (struct extent_2d *) malloc(nprocs*sizeof(struct extent_2d));
  if (outarray == nullptr) return nullptr;

  MPI_Allgather(&out,sizeof(struct extent_2d),MPI_BYTE,
                outarray,sizeof(struct extent_2d),MPI_BYTE,comm);

  // count send collides, including self

  nsend = 0;
  iproc = me;
  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    nsend += remap_2d_collide(&in,&outarray[iproc],&overlap);
  }

  // malloc space for send info

  if (nsend) {
    plan->pack = PackKokkos2d<DeviceType>::pack_2d;

    plan->send_offset = (int *) malloc(nsend*sizeof(int));
    plan->send_size = (int *) malloc(nsend*sizeof(int));
    plan->send_proc = (int *) malloc(nsend*sizeof(int));
    plan->packplan = (struct pack_plan_2d *)
      malloc(nsend*sizeof(struct pack_plan_2d));

    if (plan->send_offset == nullptr || plan->send_size == nullptr ||
        plan->send_proc == nullptr || plan->packplan == nullptr) return nullptr;
  }

  // store send info, with self as last entry

  nsend = 0;
  iproc = me;
  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    if (remap_2d_collide(&in,&outarray[iproc],&overlap)) {
      plan->send_proc[nsend] = iproc;
      plan->send_offset[nsend] = nqty * ((overlap.jlo-in.jlo)*in.isize +
          (overlap.ilo-in.ilo));
      plan->packplan[nsend].nfast = nqty*overlap.isize;
      plan->packplan[nsend].nslow = overlap.jsize;
      plan->packplan[nsend].nstride = nqty*in.isize;
      plan->packplan[nsend].nqty = nqty;
      plan->send_size[nsend] = nqty*overlap.isize*overlap.jsize;
      nsend++;
    }
  }

  // plan->nsend = # of sends not including self

  if (nsend && plan->send_proc[nsend-1] == me) {
    if (plan->usecollective) // for collectives include self in nsend list
      plan->nsend = nsend;
    else
      plan->nsend = nsend - 1;
  } else
    plan->nsend = nsend;

  // combine input extents across all procs

  MPI_Allgather(&in,sizeof(struct extent_2d),MPI_BYTE,
                inarray,sizeof(struct extent_2d),MPI_BYTE,comm);

  // count recv collides, including self

  nrecv = 0;
  iproc = me;
  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    nrecv += remap_2d_collide(&out,&inarray[iproc],&overlap);
  }

  // malloc space for recv info

  if (nrecv) {
    if (permute == 0)
      plan->unpack = PackKokkos2d<DeviceType>::unpack_2d;
    else {
      if (nqty == 1)
        plan->unpack = PackKokkos2d<DeviceType>::unpack_2d_permute_1;
      else if (nqty == 2)
        plan->unpack = PackKokkos2d<DeviceType>::unpack_2d_permute_2;
      else
        plan->unpack = PackKokkos2d<DeviceType>::unpack_2d_permute_n;
    }

    plan->recv_offset = (int *) malloc(nrecv*sizeof(int));
    plan->recv_size = (int *) malloc(nrecv*sizeof(int));
    plan->recv_proc = (int *) malloc(nrecv*sizeof(int));
    plan->recv_bufloc = (int *) malloc(nrecv*sizeof(int));
    plan->request = (MPI_Request *) malloc(nrecv*sizeof(MPI_Request));
    plan->unpackplan = (struct pack_plan_2d *)
      malloc(nrecv*sizeof(struct pack_plan_2d));

    if (plan->recv_offset == nullptr || plan->recv_size == nullptr ||
        plan->recv_proc == nullptr || plan->recv_bufloc == nullptr ||
        plan->request == nullptr || plan->unpackplan == nullptr) return nullptr;
  }

  // store recv info, with self as last entry

  ibuf = 0;
  nrecv = 0;
  iproc = me;

  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    if (remap_2d_collide(&out,&inarray[iproc],&overlap)) {
      plan->recv_proc[nrecv] = iproc;
      plan->recv_bufloc[nrecv] = ibuf;

      if (permute == 0) {
        plan->recv_offset[nrecv] = nqty * ((overlap.jlo-out.jlo)*out.isize +
                  (overlap.ilo-out.ilo));
        plan->unpackplan[nrecv].nfast = nqty*overlap.isize;
        plan->unpackplan[nrecv].nslow = overlap.jsize;
        plan->unpackplan[nrecv].nstride = nqty*out.isize;
        plan->unpackplan[nrecv].nqty = nqty;
      }
      else {
        plan->recv_offset[nrecv] = nqty * ((overlap.ilo-out.ilo)*out.jsize +
                  (overlap.jlo-out.jlo));
        plan->unpackplan[nrecv].nfast = overlap.isize;
        plan->unpackplan[nrecv].nslow = overlap.jsize;
        plan->unpackplan[nrecv].nstride = nqty*out.jsize;
        plan->unpackplan[nrecv].nqty = nqty;
      }

      plan->recv_size[nrecv] = nqty*overlap.isize*overlap.jsize;
      ibuf += plan->recv_size[nrecv];
      nrecv++;
    }
  }

  // plan->nrecv = # of recvs not including self
  // for collectives include self in the nsend list

  if (nrecv && plan->recv_proc[nrecv-1] == me) {
    if (plan->usecollective) plan->nrecv = nrecv;
    else plan->nrecv = nrecv - 1;
  } else plan->nrecv = nrecv;

  // init remaining fields in remap plan

  plan->memory = memory;

  if (nrecv == plan->nrecv) plan->self = 0;
  else plan->self = 1;

  // free locally malloced space

  free(inarray);
  free(outarray);

  // find biggest send message (not including self) and malloc space for it

  size = 0;
  for (nsend = 0; nsend < plan->nsend; nsend++)
    size = MAX(size,plan->send_size[nsend]);

  if (size) {
    plan->d_sendbuf = typename FFT_AT::t_FFT_SCALAR_1d("remap2d:sendbuf",size);
    if (!plan->d_sendbuf.data()) return nullptr;
  }

  // if requested, allocate internal scratch space for recvs,
  // only need it if I will receive any data (including self)

  if (memory == 1) {
    if (nrecv > 0) {
      plan->d_scratch =
        typename FFT_AT::t_FFT_SCALAR_1d("remap2d:scratch",nqty*out.isize*out.jsize);
      if (!plan->d_scratch.data()) return nullptr;
    }
  }

  // not using collective - dup comm

  MPI_Comm_dup(comm,&plan->comm);

  // return pointer to plan

  return plan;
}

/* ----------------------------------------------------------------------
   Destroy a 2d remap plan
------------------------------------------------------------------------- */

template<class DeviceType>
void RemapKokkos2d<DeviceType>::remap_2d_destroy_plan_kokkos(struct remap_plan_2d_kokkos<DeviceType> *plan)
{
  if (plan == nullptr) return;

  // free MPI communicator

  if (!((plan->usecollective) && (plan->commringlen == 0)))
    MPI_Comm_free(&plan->comm);

  // free internal arrays

  if (plan->nsend || plan->self) {
    free(plan->send_offset);
    free(plan->send_size);
    free(plan->send_proc);
    free(plan->packplan);
  }

  if (plan->nrecv || plan->self) {
    free(plan->recv_offset);
    free(plan->recv_size);
    free(plan->recv_proc);
    free(plan->recv_bufloc);
    free(plan->request);
    free(plan->unpackplan);
  }

  // free plan itself

  delete plan;
}

namespace SPARTA_NS {
template class RemapKokkos2d<SPADeviceType>;
#ifdef SPARTA_KOKKOS_GPU
template class RemapKokkos2d<SPAHostType>;
#endif
}

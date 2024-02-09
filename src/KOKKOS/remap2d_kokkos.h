/* ----------------------------------------------------------------------
   SPARTA - Stochastic Parallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifndef SPARTA_REMAP2D_KOKKOS_H
#define SPARTA_REMAP2D_KOKKOS_H

#include "pointers.h"
#include <mpi.h>
#include "fftdata_kokkos.h"
#include "remap2d.h"

namespace SPARTA_NS {

// details of how to do a 2d remap

template<class DeviceType>
struct remap_plan_2d_kokkos {
  typedef DeviceType device_type;
  typedef FFTArrayTypes<DeviceType> FFT_AT;
  typename FFT_AT::t_FFT_SCALAR_1d d_sendbuf;                  // buffer for MPI sends
  FFT_HAT::t_FFT_SCALAR_1d h_sendbuf;                          // host buffer for MPI sends
  typename FFT_AT::t_FFT_SCALAR_1d d_scratch;                  // scratch buffer for MPI recvs
  FFT_HAT::t_FFT_SCALAR_1d h_scratch;                          // host scratch buffer for MPI recvs
  void (*pack)(typename FFT_AT::t_FFT_SCALAR_1d_um, int, typename FFT_AT::t_FFT_SCALAR_1d_um, int, struct pack_plan_2d *);
                                    // which pack function to use
  void (*unpack)(typename FFT_AT::t_FFT_SCALAR_1d_um, int, typename FFT_AT::t_FFT_SCALAR_1d_um, int, struct pack_plan_2d *);
                                    // which unpack function to use
  int *send_offset;                 // extraction loc for each send
  int *send_size;                   // size of each send message
  int *send_proc;                   // proc to send each message to
  struct pack_plan_2d *packplan;    // pack plan for each send message
  int *recv_offset;                 // insertion loc for each recv
  int *recv_size;                   // size of each recv message
  int *recv_proc;                   // proc to recv each message from
  int *recv_bufloc;                 // offset in scratch buf for each recv
  MPI_Request *request;             // MPI request for each posted recv
  struct pack_plan_2d *unpackplan;  // unpack plan for each recv message
  int nrecv;                        // # of recvs from other procs
  int nsend;                        // # of sends to other procs
  int self;                         // whether I send/recv with myself
  int memory;                       // user provides scratch space or not
  MPI_Comm comm;                    // group of procs performing remap
  int usecollective;                // use collective or point-to-point MPI
  int commringlen;                  // length of commringlist
  int *commringlist;                // ranks on communication ring of this plan
  int usegpu_aware;                 // use GPU-Aware MPI or not
};

template<class DeviceType>
class RemapKokkos2d : protected Pointers {
 public:
  typedef DeviceType device_type;
  typedef FFTArrayTypes<DeviceType> FFT_AT;
  RemapKokkos2d(class SPARTA *);
  RemapKokkos2d(class SPARTA *, MPI_Comm,
                int,int,int,int,
                int,int,int,int,
                int,int,int,int,
                int,int);
  ~RemapKokkos2d() override;
  void perform(typename FFT_AT::t_FFT_SCALAR_1d, typename FFT_AT::t_FFT_SCALAR_1d, typename FFT_AT::t_FFT_SCALAR_1d);

  struct remap_plan_2d_kokkos<DeviceType> *plan;

  void remap_2d_kokkos(typename FFT_AT::t_FFT_SCALAR_1d, typename FFT_AT::t_FFT_SCALAR_1d, typename FFT_AT::t_FFT_SCALAR_1d, struct remap_plan_2d_kokkos<DeviceType> *);
  struct remap_plan_2d_kokkos<DeviceType> *remap_2d_create_plan_kokkos(MPI_Comm,
                                             int, int, int, int,
                                             int, int, int, int,
                                             int, int, int, int, int, int);
  void remap_2d_destroy_plan_kokkos(struct remap_plan_2d_kokkos<DeviceType> *);
};

}

#endif


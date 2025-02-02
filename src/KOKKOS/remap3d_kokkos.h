/* ----------------------------------------------------------------------
   SPARTA - Stochastic Parallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifndef SPARTA_REMAP3D_KOKKOS_H
#define SPARTA_REMAP3D_KOKKOS_H

#include "pointers.h"
#include <mpi.h>
#include "fftdata_kokkos.h"
#include "remap3d.h"

namespace SPARTA_NS {

// details of how to do a 3d remap

template<class DeviceType>
struct remap_plan_3d_kokkos {
  typedef DeviceType device_type;
  typedef FFTArrayTypes<DeviceType> FFT_AT;
  typename FFT_AT::t_FFT_SCALAR_1d d_sendbuf;                  // buffer for MPI sends
  FFT_HAT::t_FFT_SCALAR_1d h_sendbuf;                          // host buffer for MPI sends
  typename FFT_AT::t_FFT_SCALAR_1d d_scratch;                  // scratch buffer for MPI recvs
  FFT_HAT::t_FFT_SCALAR_1d h_scratch;                          // host scratch buffer for MPI recvs
  void (*pack)(typename FFT_AT::t_FFT_SCALAR_1d_um, int, typename FFT_AT::t_FFT_SCALAR_1d_um, int, struct pack_plan_3d *);
                                    // which pack function to use
  void (*unpack)(typename FFT_AT::t_FFT_SCALAR_1d_um, int, typename FFT_AT::t_FFT_SCALAR_1d_um, int, struct pack_plan_3d *);
                                    // which unpack function to use
  int *send_offset;                 // extraction loc for each send
  int *send_size;                   // size of each send message
  int *send_proc;                   // proc to send each message to
  struct pack_plan_3d *packplan;    // pack plan for each send message
  int *recv_offset;                 // insertion loc for each recv
  int *recv_size;                   // size of each recv message
  int *recv_proc;                   // proc to recv each message from
  int *recv_bufloc;                 // offset in scratch buf for each recv
  int *nrecvmap;                    // maps receive index to rank index
  MPI_Request *request;             // MPI request for each posted recv
  struct pack_plan_3d *unpackplan;  // unpack plan for each recv message
  int nrecv;                        // # of recvs from other procs
  int nsend;                        // # of sends to other procs
  int self;                         // whether I send/recv with myself
  int memory;                       // user provides scratch space or not
  MPI_Comm comm;                    // group of procs performing remap
  int usecollective;                // use collective or point-to-point MPI
  int usegpu_aware;                 // use GPU-Aware MPI or not
  // variables for collective MPI only
  int commringlen;                  // length of commringlist
  int *commringlist;                // ranks on communication ring of this plan
  int *sendcnts;                    // # of elements in send buffer for each rank
  int *rcvcnts;                     // # of elements in recv buffer for each rank
  int *sdispls;                     // extraction location in send buffer for each rank
  int *rdispls;                     // extraction location in recv buffer for each rank
  int selfcommringloc;              // current proc's location in commringlist
  int selfnsendloc;                 // current proc's location in send lists
  int selfnrecvloc;                 // current proc's location in recv lists

};

template<class DeviceType>
class RemapKokkos3d : protected Pointers {
 public:
  typedef DeviceType device_type;
  typedef FFTArrayTypes<DeviceType> FFT_AT;
  RemapKokkos3d(class SPARTA *);
  RemapKokkos3d(class SPARTA *, MPI_Comm,int,int,int,int,int,int,
        int,int,int,int,int,int,int,int,int,int,int,int);
  ~RemapKokkos3d() override;
  void perform(typename FFT_AT::t_FFT_SCALAR_1d, typename FFT_AT::t_FFT_SCALAR_1d, typename FFT_AT::t_FFT_SCALAR_1d);

  struct remap_plan_3d_kokkos<DeviceType> *plan;

  void remap_3d_kokkos(typename FFT_AT::t_FFT_SCALAR_1d, typename FFT_AT::t_FFT_SCALAR_1d, typename FFT_AT::t_FFT_SCALAR_1d, struct remap_plan_3d_kokkos<DeviceType> *);
  struct remap_plan_3d_kokkos<DeviceType> *remap_3d_create_plan_kokkos(MPI_Comm,
                                             int, int, int, int, int, int,
                                             int, int, int, int, int, int,
                                             int, int, int, int, int, int);
  void remap_3d_destroy_plan_kokkos(struct remap_plan_3d_kokkos<DeviceType> *);
};

}

#endif


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

#ifndef KOKKOS_SPARTA_H
#define KOKKOS_SPARTA_H

#include "pointers.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

class KokkosSPARTA : protected Pointers {
 public:
  int kokkos_exists;
  int comm_classic;
  int atomic_reduction;
  int prewrap;
  int auto_sync;
  int nthreads,ngpus;
  int numa;
  int need_atomics;
  int gpu_direct_flag;
  int collide_retry_flag;
  double collide_extra;

  KokkosSPARTA(class SPARTA *, int, char **);
  ~KokkosSPARTA();
  void accelerator(int, char **);

  template<class DeviceType>
  int need_dup()
  {
    int value = 0;

    if (need_atomics)
      value = std::is_same<typename NeedDup<1,DeviceType>::value,Kokkos::Experimental::ScatterDuplicated>::value;

    return value;
  }
};

}

#endif

/* ERROR/WARNING messages:

E: Invalid Kokkos command-line args

Self-explanatory.  See Section ? of the manual for details.

E: Could not determine local MPI rank for multiple GPUs with Kokkos CUDA because MPI library not recognized

The local MPI rank was not found in one of four supported environment variables.

E: GPUs are requested but Kokkos has not been compiled for CUDA

Recompile Kokkos with CUDA support to use GPUs.

E: Kokkos has been compiled for CUDA but no GPUs are requested

One or more GPUs must be used when Kokkos is compiled for CUDA.
*/

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

#ifndef SPARTA_ACCELERATOR_KOKKOS_DEFS_H
#define SPARTA_ACCELERATOR_KOKKOS_DEFS_H

#ifdef SPARTA_KOKKOS

// a CMake build with Kokkos for CUDA or HIP compiles sources that do not use
// Kokkos for the host alone (see cmake/common/process/kokkos_host_only.cmake),
// which cannot use the Kokkos headers.  SPARTA_KOKKOS_GPU must still be set for
// them, since it changes the layout of structs shared with the KOKKOS package.

#ifdef SPARTA_KOKKOS_HOST_ONLY

#define SPARTA_KOKKOS_GPU

#else

#include <Kokkos_Core.hpp>

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || defined(KOKKOS_ENABLE_SYCL) || defined(KOKKOS_ENABLE_OPENMPTARGET)

#define SPARTA_KOKKOS_GPU

#endif
#endif
#endif
#endif

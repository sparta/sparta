/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefid-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#ifndef SPARTA_KOKKOS_BASE_FFT_H
#define SPARTA_KOKKOS_BASE_FFT_H

#include "fftdata_kokkos.h"

namespace SPARTA_NS {

class KokkosBaseFFT {
 public:
  KokkosBaseFFT() {}

  // Kspace
  virtual void pack_forward_grid_kokkos(int, FFT_DAT::tdual_FFT_SCALAR_1d &, int, DAT::tdual_int_2d &, int) {};
  virtual void unpack_forward_grid_kokkos(int, FFT_DAT::tdual_FFT_SCALAR_1d &, int, int, DAT::tdual_int_2d &, int) {};
  virtual void pack_reverse_grid_kokkos(int, FFT_DAT::tdual_FFT_SCALAR_1d &, int, DAT::tdual_int_2d &, int) {};
  virtual void unpack_reverse_grid_kokkos(int, FFT_DAT::tdual_FFT_SCALAR_1d &, int, int, DAT::tdual_int_2d &, int) {};
};

}

#endif


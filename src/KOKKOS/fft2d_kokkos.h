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

#ifndef SPA_FFT2D_KOKKOS_H
#define SPA_FFT2D_KOKKOS_H

#include "pointers.h"
#include "remap2d_kokkos.h"
#include "fftdata_kokkos.h"

namespace SPARTA_NS {

// -------------------------------------------------------------------------

// plan for how to perform a 2d FFT

template<class DeviceType>
struct fft_plan_2d_kokkos {
  typedef DeviceType device_type;
  typedef FFTArrayTypes<DeviceType> FFT_AT;

  struct remap_plan_2d_kokkos<DeviceType> *pre_plan;       // remap from input -> 1st FFTs
  struct remap_plan_2d_kokkos<DeviceType> *mid_plan;      // remap from 1st -> 2nd FFTs
  struct remap_plan_2d_kokkos<DeviceType> *post_plan;      // remap from 2nd FFTs -> output
  typename FFT_AT::t_FFT_DATA_1d d_copy;                   // memory for remap results (if needed)
  typename FFT_AT::t_FFT_DATA_1d d_scratch;                // scratch space for remaps
  int total1,total2;                // # of 1st and 2nd FFTs (times length)
  int length1,length2;              // length of 1st and 2nd FFTs
  int pre_target,mid_target;        // where to put remap results
  int scaled;                       // whether to scale FFT results
  int normnum;                      // # of values to rescale
  double norm;                      // normalization factor for rescaling

                                    // system specific 1d FFT info
#if defined(FFT_MKL)
  DFTI_DESCRIPTOR *handle_fast;
  DFTI_DESCRIPTOR *handle_slow;
#elif defined(FFT_FFTW3)
  FFTW_API(plan) plan_fast_forward;
  FFTW_API(plan) plan_fast_backward;
  FFTW_API(plan) plan_slow_forward;
  FFTW_API(plan) plan_slow_backward;
#elif defined(FFT_CUFFT)
  cufftHandle plan_fast;
  cufftHandle plan_slow;
#elif defined(FFT_HIPFFT)
  hipfftHandle plan_fast;
  hipfftHandle plan_slow;
#else
  kiss_fft_state_kokkos<DeviceType> cfg_fast_forward;
  kiss_fft_state_kokkos<DeviceType> cfg_fast_backward;
  kiss_fft_state_kokkos<DeviceType> cfg_slow_forward;
  kiss_fft_state_kokkos<DeviceType> cfg_slow_backward;
#endif
};

template<class DeviceType>
class FFT2dKokkos : protected Pointers {
 public:
  enum{FORWARD=1,BACKWARD=-1};
  typedef DeviceType device_type;
  typedef FFTArrayTypes<DeviceType> FFT_AT;

  FFT2dKokkos(class SPARTA *, MPI_Comm, int, int,
        int,int,int,int,
        int,int,int,int,
        int,int,int *,int,int);
  ~FFT2dKokkos() override;
  void compute(typename FFT_AT::t_FFT_SCALAR_1d, typename FFT_AT::t_FFT_SCALAR_1d, int);
  void timing1d(typename FFT_AT::t_FFT_SCALAR_1d, int, int);

 private:
  struct fft_plan_2d_kokkos<DeviceType> *plan;
  RemapKokkos2d<DeviceType> *remapKK;

#ifdef FFT_KISSFFT
  KissFFTKokkos<DeviceType> *kissfftKK;
#endif

  void fft_2d_kokkos(typename FFT_AT::t_FFT_DATA_1d, typename FFT_AT::t_FFT_DATA_1d, int, struct fft_plan_2d_kokkos<DeviceType> *);

  struct fft_plan_2d_kokkos<DeviceType> *fft_2d_create_plan_kokkos(MPI_Comm, int, int,
                                         int, int, int, int,
                                         int, int, int, int,
                                         int, int, int *, int, int, int);

  void fft_2d_destroy_plan_kokkos(struct fft_plan_2d_kokkos<DeviceType> *);

  void fft_2d_1d_only_kokkos(typename FFT_AT::t_FFT_DATA_1d, int, int, struct fft_plan_2d_kokkos<DeviceType> *);

};

}

#endif


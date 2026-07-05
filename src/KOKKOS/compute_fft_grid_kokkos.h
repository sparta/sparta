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

#ifdef COMPUTE_CLASS

ComputeStyle(fft/grid/kk,ComputeFFTGridKokkos)

#else

#ifndef SPARTA_COMPUTE_FFT_GRID_KOKKOS_H
#define SPARTA_COMPUTE_FFT_GRID_KOKKOS_H

#include "compute_fft_grid.h"
#include "irregular_kokkos.h"
#include "fftdata_kokkos.h"
#include "fft2d_kokkos.h"
#include "fft3d_kokkos.h"
#include "kokkos_base.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

class ComputeFFTGridKokkos : public ComputeFFTGrid, public KokkosBase {
 public:

  ComputeFFTGridKokkos(class SPARTA *, int, char **);

  ~ComputeFFTGridKokkos();
  void post_constructor();
  void compute_per_grid();
  void compute_per_grid_kokkos();
  void reallocate();

  void irregular_create();

 private:

  DAT::tdual_float_1d k_ingrid; // input grid values from compute,fix,variable
                                // may be NULL if ingridptr just points to c/f/v
  DAT::t_float_1d d_ingrid;

  DAT::tdual_float_1d k_fftwork; // work buf in FFT decomp, length = nfft
  DAT::t_float_1d d_fftwork;
  DAT::t_char_1d d_fftwork_char;
  DAT::t_float_1d d_gridwork;    // work buf in grid decomp, length = nglocal
  DAT::t_char_1d d_gridwork_char;

  FFT_DAT::t_FFT_SCALAR_1d d_fft; // complex buf for performing FFT, length = nfft
  DAT::t_char_1d d_fft_char;
  FFT_DAT::t_FFT_SCALAR_1d d_gridworkcomplex; // work buf in grid decomp, length = nglocal
  DAT::t_char_1d d_gridworkcomplex_char;

  DAT::t_int_1d d_map1; // mapping of received SPARTA grid values to FFT grid
                        // map1[i] = index into ordered FFT grid of
                        //           Ith value in buffer received
                        //           from SPARTA decomp via irregular comm
  DAT::t_int_1d d_map2; // mapping of received FFT grid values to SPARTA grid
                        // map2[i] = index into SPARTA grid of Ith value
                        //           in buffer received from FFT decomp via
                        //           irregular comm

  DAT::tdual_float_2d_lr k_array_grid;
  DAT::tdual_float_1d k_vector_grid;

  FFT2dKokkos<DeviceType> *fft2dKK;
  FFT3dKokkos<DeviceType> *fft3dKK;
  IrregularKokkos *irregular1KK,*irregular2KK;

  void fft_create();
  void print_FFT_info();
};

}

#endif
#endif

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

// common FFT library related defines and compilation settings

#ifndef SPARTA_FFT_SETTINGS_H
#define SPARTA_FFT_SETTINGS_H

#include "accelerator_kokkos_defs.h"

// if a user sets FFTW, it means FFTW3

#ifdef FFT_FFTW
#ifndef FFT_FFTW3
#undef FFT_FFTW
#define FFT_FFTW3
#endif
#endif

#ifdef SPARTA_KOKKOS
# ifdef FFT_KOKKOS_FFTW
#  undef FFT_KOKKOS_FFTW
#  define FFT_KOKKOS_FFTW3
# endif
# ifdef FFT_KOKKOS_FFTW_THREADS
#  if !defined(FFT_KOKKOS_FFTW3)
#   error "Must use -DFFT_KOKKOS_FFTW3 with -DFFT_KOKKOS_FFTW_THREADS"
#  endif
# endif
#endif

// set strings for library info output

#if defined(FFT_FFTW3)
#define SPARTA_FFT_LIB "FFTW3"
#elif defined(FFT_MKL)
#define SPARTA_FFT_LIB "MKL FFT"
#elif defined(FFT_CUFFT)
#define SPARTA_FFT_LIB "cuFFT"
#elif defined(FFT_HIPFFT)
#define SPARTA_FFT_LIB "hipFFT"
#else
#define SPARTA_FFT_LIB "KISS FFT"
#endif

#ifdef SPARTA_KOKKOS

// with KOKKOS in CUDA or HIP mode we can only have
//  CUFFT/HIPFFT or KISS, thus undefine all other
//  FFTs here

#ifdef KOKKOS_ENABLE_CUDA
# if defined(FFT_KOKKOS_FFTW)
#  undef FFT_KOKKOS_FFTW
# endif
# if defined(FFT_KOKKOS_FFTW3)
#  undef FFT_KOKKOS_FFTW3
# endif
# if defined(FFT_KOKKOS_MKL)
#  undef FFT_KOKKOS_MKL
# endif
# if !defined(FFT_KOKKOS_CUFFT) && !defined(FFT_KOKKOS_KISS)
#  define FFT_KOKKOS_KISS
# endif
#elif defined(KOKKOS_ENABLE_HIP)
# if defined(FFT_KOKKOS_FFTW)
#  undef FFT_KOKKOS_FFTW
# endif
# if defined(FFT_KOKKOS_FFTW3)
#  undef FFT_KOKKOS_FFTW3
# endif
# if defined(FFT_KOKKOS_MKL)
#  undef FFT_KOKKOS_MKL
# endif
# if !defined(FFT_KOKKOS_HIPFFT) && !defined(FFT_KOKKOS_KISS)
#  define FFT_KOKKOS_KISS
# endif
#else
# if defined(FFT_KOKKOS_CUFFT)
#  error "Must enable CUDA with KOKKOS to use -DFFT_KOKKOS_CUFFT"
# endif
# if defined(FFT_KOKKOS_HIPFFT)
#  error "Must enable HIP with KOKKOS to use -DFFT_KOKKOS_HIPFFT"
# endif
#endif

#if defined(FFT_KOKKOS_CUFFT)
#define SPARTA_FFT_KOKKOS_LIB "cuFFT"
#elif defined(FFT_KOKKOS_HIPFFT)
#define SPARTA_FFT_KOKKOS_LIB "hipFFT"
#elif defined(FFT_KOKKOS_FFTW3)
#define SPARTA_FFT_KOKKOS_LIB "FFTW3"
#elif defined(FFT_KOKKOS_MKL)
#define SPARTA_FFT_KOKKOS_LIB "MKL FFT"
#else
#define SPARTA_FFT_KOKKOS_LIB "KISS FFT"
#endif

#else
#define SPARTA_FFT_KOKKOS_LIB ""
#endif

#ifdef FFT_SINGLE
typedef float FFT_SCALAR;
#define FFT_PRECISION 1
#define SPARTA_FFT_PREC "single"
#define MPI_FFT_SCALAR MPI_FLOAT
#else

typedef double FFT_SCALAR;
#define FFT_PRECISION 2
#define SPARTA_FFT_PREC "double"
#define MPI_FFT_SCALAR MPI_DOUBLE
#endif

#endif

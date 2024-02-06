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

#ifndef SPA_FFT_SETTINGS_H
#define SPA_FFT_SETTINGS_H

// if user set FFTW, it means FFTW3

#ifdef FFT_FFTW
#ifndef FFT_FFTW3
#define FFT_FFTW3
#endif
#endif

// set strings for library info output

#if defined(FFT_FFTW3)
#define SPA_FFT_LIB "FFTW3"
#elif defined(FFT_MKL)
#define SPA_FFT_LIB "MKL FFT"
#elif defined(FFT_CUFFT)
#define SPA_FFT_LIB "cuFFT"
#elif defined(FFT_HIPFFT)
#define SPA_FFT_LIB "hipFFT"
#else
#define SPA_FFT_LIB "KISS FFT"
#endif

#ifdef FFT_SINGLE
typedef float FFT_SCALAR;
#define FFT_PRECISION 1
#define SPA_FFT_PREC "single"
#define MPI_FFT_SCALAR MPI_FLOAT
#else

typedef double FFT_SCALAR;
#define FFT_PRECISION 2
#define SPA_FFT_PREC "double"
#define MPI_FFT_SCALAR MPI_DOUBLE
#endif

#endif

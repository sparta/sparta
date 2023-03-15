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

// data types for 2d/3d FFTs

#ifndef FFT_DATA_H
#define FFT_DATA_H

// -------------------------------------------------------------------------

// Data types for single-precision complex

#if FFT_PRECISION == 1

#if defined(FFT_MKL)
#include "mkl_dfti.h"
typedef float _Complex FFT_DATA;
#define FFT_MKL_PREC DFTI_SINGLE

#elif defined(FFT_FFTW2)
#if defined(FFTW_SIZE)
#include "sfftw.h"
#else
#include "fftw.h"
#endif
typedef FFTW_COMPLEX FFT_DATA;

#elif defined(FFT_FFTW3)
#include "fftw3.h"
typedef fftwf_complex FFT_DATA;
#define FFTW_API(function)  fftwf_ ## function

#else

// use a stripped down version of kiss fft as default fft

#ifndef FFT_KISSFFT
#define FFT_KISSFFT
#endif
#define kiss_fft_scalar float
typedef struct {
    kiss_fft_scalar re;
    kiss_fft_scalar im;
} FFT_DATA;

struct kiss_fft_state;
typedef struct kiss_fft_state* kiss_fft_cfg;
#endif

// -------------------------------------------------------------------------

// Data types for double-precision complex

#elif FFT_PRECISION == 2

#if defined(FFT_MKL)
#include "mkl_dfti.h"
typedef double _Complex FFT_DATA;
#define FFT_MKL_PREC DFTI_DOUBLE

#elif defined(FFT_FFTW2)
#if defined(FFTW_SIZE)
#include "dfftw.h"
#else
#include "fftw.h"
#endif
typedef FFTW_COMPLEX FFT_DATA;

#elif defined(FFT_FFTW3)
#include "fftw3.h"
typedef fftw_complex FFT_DATA;
#define FFTW_API(function)  fftw_ ## function

#else

// use a stripped down version of kiss fft as default fft

#ifndef FFT_KISSFFT
#define FFT_KISSFFT
#endif
#define kiss_fft_scalar double
typedef struct {
    kiss_fft_scalar re;
    kiss_fft_scalar im;
} FFT_DATA;

struct kiss_fft_state;
typedef struct kiss_fft_state* kiss_fft_cfg;
#endif

// -------------------------------------------------------------------------

#else
#error "FFT_PRECISION needs to be either 1 (=single) or 2 (=double)"
#endif

// -------------------------------------------------------------------------

#endif

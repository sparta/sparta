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

// data types for 2d/3d FFTs

#ifndef FFT_DATA_H
#define FFT_DATA_H

// -------------------------------------------------------------------------

// Data types for single-precision complex

#if FFT_PRECISION == 1

#if defined(FFT_SGI)
#include "fft.h"
typedef complex FFT_DATA;
#define FFT_1D cfft1d
#define FFT_1D_INIT cfft1di
extern "C" {
  int cfft1d(int, int, FFT_DATA *, int, FFT_DATA *);
  FFT_DATA *cfft1di(int, FFT_DATA *);
}

#elif defined(FFT_SCSL)
#include <scsl_fft.h>
typedef scsl_complex FFT_DATA;
typedef float FFT_PREC;
#define FFT_1D ccfft
#define FFT_1D_INIT ccfft
extern "C" {
  int ccfft(int, int, FFT_PREC, FFT_DATA *, FFT_DATA *,
                      FFT_PREC *, FFT_PREC *, int *);
}

#elif defined(FFT_ACML)
typedef struct {
  float re;
  float im;
} FFT_DATA;
#define FFT_1D cfft1m_
extern "C" {
  void cfft1m_(int *, int *, int *, FFT_DATA *, FFT_DATA *, int *);
}

#elif defined(FFT_INTEL)
typedef struct {
  float re;
  float im;
} FFT_DATA;
#define FFT_1D cfft1d_
#define FFT_1D_INIT cfft1d_
extern "C" {
  void cfft1d_(FFT_DATA *, int *, int *, FFT_DATA *);
}

#elif defined(FFT_MKL)
#include "mkl_dfti.h"
typedef float _Complex FFT_DATA;
#define FFT_MKL_PREC DFTI_SINGLE

#elif defined(FFT_DEC)
typedef struct {
  float re;
  float im;
} FFT_DATA;
#define FFT_1D cfft_
extern "C" {
  void cfft_(char *, char *, char *, FFT_DATA *, FFT_DATA *, int *, int *);
}

#elif defined(FFT_T3E)
#include <complex.h>
typedef complex single FFT_DATA;
#define FFT_1D GGFFT
#define FFT_1D_INIT GGFFT
extern "C" {
  void GGFFT(int *, int *, double *, FFT_DATA *, FFT_DATA *,
             double *, double *, int *);
}

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

#if defined(FFT_SGI)
#include "fft.h"
typedef zomplex FFT_DATA;
#define FFT_1D zfft1d
#define FFT_1D_INIT zfft1di
extern "C" {
  int zfft1d(int, int, FFT_DATA *, int, FFT_DATA *);
  FFT_DATA *zfft1di(int, FFT_DATA *);
}

#elif defined(FFT_SCSL)
#include <scsl_fft.h>
typedef scsl_zomplex FFT_DATA;
typedef double FFT_PREC;
#define FFT_1D zzfft
#define FFT_1D_INIT zzfft
extern "C" {
  int zzfft(int, int, FFT_PREC, FFT_DATA *, FFT_DATA *,
                      FFT_PREC *, FFT_PREC *, int *);
}

#elif defined(FFT_ACML)
typedef struct {
  double re;
  double im;
} FFT_DATA;
#define FFT_1D zfft1m_
extern "C" {
  void zfft1m_(int *, int *, int *, FFT_DATA *, FFT_DATA *, int *);
}

#elif defined(FFT_INTEL)
typedef struct {
  double re;
  double im;
} FFT_DATA;
#define FFT_1D zfft1d_
#define FFT_1D_INIT zfft1d_
extern "C" {
  void zfft1d_(FFT_DATA *, int *, int *, FFT_DATA *);
}

#elif defined(FFT_MKL)
#include "mkl_dfti.h"
typedef double _Complex FFT_DATA;
#define FFT_MKL_PREC DFTI_DOUBLE

#elif defined(FFT_DEC)
typedef struct {
  double re;
  double im;
} FFT_DATA;
#define FFT_1D zfft_
extern "C" {
  void zfft_(char *, char *, char *, FFT_DATA *, FFT_DATA *, int *, int *);
}

#elif defined(FFT_T3E)
#include <complex.h>
typedef complex double FFT_DATA;
#define FFT_1D CCFFT
#define FFT_1D_INIT CCFFT
extern "C" {
  void CCFFT(int *, int *, double *, FFT_DATA *, FFT_DATA *,
             double *, double *, int *);
}

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

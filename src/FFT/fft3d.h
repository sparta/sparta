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

// User-settable FFT precision

// FFT_PRECISION = 1 is single-precision complex (4-byte real, 4-byte imag)
// FFT_PRECISION = 2 is double-precision complex (8-byte real, 8-byte imag)

#ifdef FFT_SINGLE
#define FFT_PRECISION 1
typedef float FFT_SCALAR;
#else
#define FFT_PRECISION 2
typedef double FFT_SCALAR;
#endif

#ifdef FFT_FFTW
#define FFT_FFTW2
#endif

#include "fftdata.h"

// -------------------------------------------------------------------------

// plan for how to perform a 3d FFT

struct fft_plan_3d {
  struct remap_plan_3d *pre_plan;       // remap from input -> 1st FFTs
  struct remap_plan_3d *mid1_plan;      // remap from 1st -> 2nd FFTs
  struct remap_plan_3d *mid2_plan;      // remap from 2nd -> 3rd FFTs
  struct remap_plan_3d *post_plan;      // remap from 3rd FFTs -> output
  FFT_DATA *copy;                   // memory for remap results (if needed)
  FFT_DATA *scratch;                // scratch space for remaps
  int total1,total2,total3;         // # of 1st,2nd,3rd FFTs (times length)
  int length1,length2,length3;      // length of 1st,2nd,3rd FFTs
  int pre_target;                   // where to put remap results
  int mid1_target,mid2_target;
  int scaled;                       // whether to scale FFT results
  int normnum;                      // # of values to rescale
  double norm;                      // normalization factor for rescaling

                                    // system specific 1d FFT info
#if defined(FFT_SGI)
  FFT_DATA *coeff1;
  FFT_DATA *coeff2;
  FFT_DATA *coeff3;
#elif defined(FFT_SCSL)
  FFT_PREC *coeff1;
  FFT_PREC *coeff2;
  FFT_PREC *coeff3;
  FFT_PREC *work1;
  FFT_PREC *work2;
  FFT_PREC *work3;
#elif defined(FFT_ACML)
  FFT_DATA *coeff1;
  FFT_DATA *coeff2;
  FFT_DATA *coeff3;
#elif defined(FFT_INTEL)
  FFT_DATA *coeff1;
  FFT_DATA *coeff2;
  FFT_DATA *coeff3;
#elif defined(FFT_MKL)
  DFTI_DESCRIPTOR *handle_fast;
  DFTI_DESCRIPTOR *handle_mid;
  DFTI_DESCRIPTOR *handle_slow;
#elif defined(FFT_T3E)
  double *coeff1;
  double *coeff2;
  double *coeff3;
  double *work1;
  double *work2;
  double *work3;
#elif defined(FFT_FFTW2)
  fftw_plan plan_fast_forward;
  fftw_plan plan_fast_backward;
  fftw_plan plan_mid_forward;
  fftw_plan plan_mid_backward;
  fftw_plan plan_slow_forward;
  fftw_plan plan_slow_backward;
#elif defined(FFT_FFTW3)
  FFTW_API(plan) plan_fast_forward;
  FFTW_API(plan) plan_fast_backward;
  FFTW_API(plan) plan_mid_forward;
  FFTW_API(plan) plan_mid_backward;
  FFTW_API(plan) plan_slow_forward;
  FFTW_API(plan) plan_slow_backward;
#elif defined(FFT_KISSFFT)
  kiss_fft_cfg cfg_fast_forward;
  kiss_fft_cfg cfg_fast_backward;
  kiss_fft_cfg cfg_mid_forward;
  kiss_fft_cfg cfg_mid_backward;
  kiss_fft_cfg cfg_slow_forward;
  kiss_fft_cfg cfg_slow_backward;
#endif
};

// function prototypes

extern "C" {
  void fft_3d(FFT_DATA *, FFT_DATA *, int, struct fft_plan_3d *);
  struct fft_plan_3d *fft_3d_create_plan(MPI_Comm, int, int, int,
                                         int, int, int, int, int,
                                         int, int, int, int, int, int, int,
                                         int, int, int *, int);
  void fft_3d_destroy_plan(struct fft_plan_3d *);
  void factor3d(int, int *, int *);
  void bifactor(int, int *, int *);
  void fft_3d_1d_only(FFT_DATA *, int, int, struct fft_plan_3d *);
}

/* ERROR/WARNING messages:

*/

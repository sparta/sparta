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

#include <mpi.h>

#include "spafftsettings.h"
#include "fftdata.h"

// -------------------------------------------------------------------------

// plan for how to perform a 2d FFT

struct fft_plan_2d {
  struct remap_plan_2d *pre_plan;       // remap from input -> 1st FFTs
  struct remap_plan_2d *mid_plan;       // remap from 1st -> 2nd FFTs
  struct remap_plan_2d *post_plan;      // remap from 2nd FFTs -> output
  FFT_DATA *copy;                   // memory for remap results (if needed)
  FFT_DATA *scratch;                // scratch space for remaps
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
#elif defined(FFT_KISS)
  kiss_fft_cfg cfg_fast_forward;
  kiss_fft_cfg cfg_fast_backward;
  kiss_fft_cfg cfg_slow_forward;
  kiss_fft_cfg cfg_slow_backward;
#endif
};

// function prototypes

extern "C" {
  void fft_2d(FFT_DATA *, FFT_DATA *, int, struct fft_plan_2d *);
  struct fft_plan_2d *fft_2d_create_plan(MPI_Comm, int, int,
                                         int, int, int, int, int, int, int, int,
                                         int, int, int *, int);
  void fft_2d_destroy_plan(struct fft_plan_2d *);
  void factor_2d(int, int *, int *);
  void fft_2d_1d_only(FFT_DATA *, int, int, struct fft_plan_2d *);
}

/* ERROR/WARNING messages:

*/

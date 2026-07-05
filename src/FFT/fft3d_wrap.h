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

#ifndef SPARTA_FFT3D_WRAP_H
#define SPARTA_FFT3D_WRAP_H

#include "fft3d.h"
#include "pointers.h"

#ifdef FFT_HEFFTE
#include "heffte.h"
// select the backend
#if defined(FFT_HEFFTE_FFTW)
using heffte_backend = heffte::backend::fftw;
#elif defined(FFT_HEFFTE_MKL)
using heffte_backend = heffte::backend::mkl;
#else
using heffte_backend = heffte::backend::stock;
#endif

#endif // FFT_HEFFTE

namespace SPARTA_NS {

class FFT3d : protected Pointers {
 public:
  enum { FORWARD = 1, BACKWARD = -1 };

  FFT3d(class SPARTA *, MPI_Comm, int, int, int, int, int, int, int, int, int, int, int, int, int,
        int, int, int, int, int *, int);
  ~FFT3d() override;
  void compute(FFT_SCALAR *, FFT_SCALAR *, int);
  void timing1d(FFT_SCALAR *, int, int);

 private:
  #ifdef FFT_HEFFTE
  // the heFFTe plan supersedes the internal fft_plan_3d
  std::unique_ptr<heffte::fft3d<heffte_backend>> heffte_plan;
  std::vector<std::complex<FFT_SCALAR>> heffte_workspace;
  heffte::scale hscale;
  #else
  struct fft_plan_3d *plan;
  #endif
};

}    // namespace SPARTA_NS

#endif

/* ERROR/WARNING messages:

E: Could not create 3d FFT plan

The FFT setup for the PPPM solver failed, typically due
to lack of memory.  This is an unusual error.  Check the
size of the FFT grid you are requesting.

*/

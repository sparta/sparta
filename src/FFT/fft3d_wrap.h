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

#ifndef SPA_FFT3D_WRAP_H
#define SPA_FFT3D_WRAP_H

#include "pointers.h"
#include "fft3d.h"

namespace SPARTA_NS {

class FFT3D : protected Pointers {
 public:
  FFT3D(class SPARTA *, MPI_Comm,
        int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,
        int,int,int *,int);
  ~FFT3D();
  void compute(FFT_SCALAR *, FFT_SCALAR *, int);
  void timing1d(FFT_SCALAR *, int, int);

 private:
  struct fft_plan_3d *plan;
};

}

#endif

/* ERROR/WARNING messages:

E: Could not create 3d FFT plan

The FFT setup for the PPPM solver failed, typically due
to lack of memory.  This is an unusual error.  Check the
size of the FFT grid you are requesting.

*/

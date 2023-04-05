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

#include <mpi.h>
#include "fft2d_wrap.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

FFT2D::FFT2D(SPARTA *spa, MPI_Comm comm, int nfast, int nslow,
             int in_ilo, int in_ihi, int in_jlo, int in_jhi,
             int out_ilo, int out_ihi, int out_jlo, int out_jhi,
             int scaled, int permute, int *nbuf, int usecollective) :
  Pointers(spa)
{
  plan = fft_2d_create_plan(comm,nfast,nslow,
                            in_ilo,in_ihi,in_jlo,in_jhi,
                            out_ilo,out_ihi,out_jlo,out_jhi,
                            scaled,permute,nbuf,usecollective);
  if (plan == NULL) error->one(FLERR,"Could not create 2d FFT plan");
}

/* ---------------------------------------------------------------------- */

FFT2D::~FFT2D()
{
  fft_2d_destroy_plan(plan);
}

/* ---------------------------------------------------------------------- */

void FFT2D::compute(FFT_SCALAR *in, FFT_SCALAR *out, int flag)
{
  fft_2d((FFT_DATA *) in,(FFT_DATA *) out,flag,plan);
}

/* ---------------------------------------------------------------------- */

void FFT2D::timing1d(FFT_SCALAR *in, int nsize, int flag)
{
  fft_2d_1d_only((FFT_DATA *) in,nsize,flag,plan);
}

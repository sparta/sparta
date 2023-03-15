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
#include "fft3d_wrap.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

FFT3D::FFT3D(SPARTA *spa, MPI_Comm comm, int nfast, int nmid, int nslow,
             int in_ilo, int in_ihi, int in_jlo, int in_jhi,
             int in_klo, int in_khi,
             int out_ilo, int out_ihi, int out_jlo, int out_jhi,
             int out_klo, int out_khi,
             int scaled, int permute, int *nbuf, int usecollective) :
  Pointers(spa)
{
  plan = fft_3d_create_plan(comm,nfast,nmid,nslow,
                            in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
                            out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,
                            scaled,permute,nbuf,usecollective);
  if (plan == NULL) error->one(FLERR,"Could not create 3d FFT plan");
}

/* ---------------------------------------------------------------------- */

FFT3D::~FFT3D()
{
  fft_3d_destroy_plan(plan);
}

/* ---------------------------------------------------------------------- */

void FFT3D::compute(FFT_SCALAR *in, FFT_SCALAR *out, int flag)
{
  fft_3d((FFT_DATA *) in,(FFT_DATA *) out,flag,plan);
}

/* ---------------------------------------------------------------------- */

void FFT3D::timing1d(FFT_SCALAR *in, int nsize, int flag)
{
  fft_3d_1d_only((FFT_DATA *) in,nsize,flag,plan);
}

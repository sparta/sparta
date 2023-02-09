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

#include "comm.h"
#include "rand_pool_wrap.h"
#include "sparta.h"
#include "kokkos.h"
#include "random_knuth.h"
#include "random_mars.h"
#include "update.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

RandPoolWrap::RandPoolWrap(int, SPARTA *sparta) : Pointers(sparta)
{
  random_thr =  NULL;
  nthreads = sparta->kokkos->nthreads;
}

/* ---------------------------------------------------------------------- */

RandPoolWrap::~RandPoolWrap()
{

}

void RandPoolWrap::destroy()
{
  if (random_thr) {
    for (int i=1; i < nthreads; ++i)
      delete random_thr[i];

    delete[] random_thr;
    random_thr = NULL;
  }
}

void RandPoolWrap::init(RanKnuth* random)
{
  // deallocate pool of RNGs
  if (random_thr) {
    for (int i=1; i < this->nthreads; ++i)
      delete random_thr[i];

    delete[] random_thr;
  }

  // allocate pool of RNGs
  // generate a random number generator instance for
  // all threads != 0. make sure we use unique seeds.
  nthreads = sparta->kokkos->nthreads;
  random_thr = new RanKnuth*[nthreads];
  for (int tid = 1; tid < nthreads; ++tid) {
    random_thr[tid] = new RanKnuth(update->ranmaster->uniform());
    double seed = update->ranmaster->uniform();
    random_thr[tid]->reset(seed,comm->me + comm->nprocs*tid,100);
  }

  // to ensure full compatibility with the serial style
  // we use the serial random number generator instance for thread 0
  random_thr[0] = random;
}
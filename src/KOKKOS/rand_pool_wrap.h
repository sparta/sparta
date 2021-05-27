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

#ifndef RAND_POOL_WRAP_H
#define RAND_POOL_WRAP_H

#include "pointers.h"
#include "kokkos_type.h"
#include "random_park.h"
#include "error.h"

namespace SPARTA_NS {

struct RandWrap {
  class RanPark* rng;

  KOKKOS_INLINE_FUNCTION
  RandWrap() {
    rng = NULL;
  }

  KOKKOS_INLINE_FUNCTION
  double drand() {
    return rng->uniform();
  }

  KOKKOS_INLINE_FUNCTION
  double normal() {
    return rng->gaussian();
  }
};

class RandPoolWrap : protected Pointers {
 public:
  RandPoolWrap(int, class SPARTA *);
  ~RandPoolWrap();
  void destroy();
  void init(RanPark*);

  KOKKOS_INLINE_FUNCTION
  RandWrap get_state() const
  {
#ifdef SPARTA_KOKKOS_GPU
    error->all(FLERR,"Cannot use Park RNG with GPUs");
#endif

    RandWrap rand_wrap;

    typedef Kokkos::Experimental::UniqueToken<
      DeviceType, Kokkos::Experimental::UniqueTokenScope::Global> unique_token_type;

    unique_token_type unique_token;
    int tid = unique_token.acquire();
    rand_wrap.rng = random_thr[tid];
    unique_token.release(tid);

    return rand_wrap;
  }

  KOKKOS_INLINE_FUNCTION
  void free_state(RandWrap) const
  {

  }

 private:
  class RanPark **random_thr;
  int nthreads;
};

}

#endif

/* ERROR/WARNING messages:

*/

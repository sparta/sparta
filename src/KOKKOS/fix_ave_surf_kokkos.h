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

#ifdef FIX_CLASS

FixStyle(ave/surf/kk,FixAveSurfKokkos)

#else

#ifndef SPARTA_FIX_AVE_SURF_KOKKOS_H
#define SPARTA_FIX_AVE_SURF_KOKKOS_H

#include "fix_ave_surf.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

struct TagFixAveSurf_Add_tally{};

class FixAveSurfKokkos : public FixAveSurf {
 public:
  FixAveSurfKokkos(class SPARTA *, int, char **);
  ~FixAveSurfKokkos();
  void init();
  void setup();
  void end_of_step();

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixAveSurf_Add_tally, const int&) const;

 private:
  int kokkosable;        // 1 if the all-tally path runs on device
                         // 0 if delegating entirely to the host base class
  int nstally;           // # of local surf rows in a compute tally
                         //   = surf->nlocal + surf->nghost
  int acc_m,acc_col;     // value index / compute-tally column for current kernel

  DAT::t_float_2d_lr d_acc;     // per-interval tally accumulator [nstally][nvalues]
  DAT::t_float_2d_lr d_tally;   // current compute device tally (set per value)

  surfint *tally2surf_all;      // surfID of each local surf row (host)
  double *acc_local_vec;        // host copy of d_acc for collate (nvalues == 1)
  double **acc_local;           // host copy of d_acc for collate (nvalues > 1)

  void reallocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/

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

#ifndef KOKKOS_BASE_H
#define KOKKOS_BASE_H

#include "kokkos_type.h"

namespace SPARTA_NS {

class KokkosBase {
 public:
  KokkosBase() {}

  // Compute
  virtual void compute_per_grid_kokkos() {}
  virtual int query_tally_grid_kokkos(DAT::t_float_2d_lr&) {return 0;}
  virtual void post_process_grid_kokkos(int, int, DAT::t_float_2d_lr, int *,
                                   DAT::t_float_1d_strided) {}

  DAT::t_float_1d d_vector;          // Kokkos version of computed vector
  DAT::t_float_2d_lr d_array_grid;   // Kokkos version of computed per-grid array
  DAT::t_float_2d_lr d_array_particle;   // Kokkos version of computed per-particle array

  DAT::tdual_float_2d_lr k_array;    // Kokkos version of computed array
  DAT::t_float_2d_lr d_array;        // Kokkos version of computed array
};

}

#endif

/* ERROR/WARNING messages:

*/

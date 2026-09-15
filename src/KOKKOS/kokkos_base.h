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

#ifndef KOKKOS_BASE_H
#define KOKKOS_BASE_H

#include "kokkos_type.h"
#include "region_prim_kokkos.h"

namespace SPARTA_NS {

class KokkosBase {
 public:
  KokkosBase() {}

  // Compute
  virtual void compute_per_grid_kokkos() {}

  // Bring this style's per-grid output up to date on the device.
  //
  // d_vector_grid and d_array_grid below are read by other Kokkos styles --
  // compute lambda/grid, compute dt/grid, fix ave/histo and the rest all take
  // the handle straight out of the producer and launch a kernel on it.  A
  // compute is asked to recompute first, so its device copy is always fresh; a
  // fix is not, and fix ave/grid keeps its per-grid arrays on the host while
  // grid adaptation moves cells between them, leaving the device holding the
  // pre-migration values.  A style that can be in that position overrides this
  // and the reader calls it before it reads.
  virtual void sync_pergrid_device_kokkos() {}
  virtual int query_tally_grid_kokkos(DAT::t_float_2d_lr&) {return 0;}
  virtual void post_process_grid_kokkos(int, int, DAT::t_float_2d_lr, int *,
                                   DAT::t_float_1d_strided) {}

  //DAT::t_float_1d d_vector;        // Kokkos version of global vector
  DAT::t_float_2d_lr d_array;        // Kokkos version of global array
  DAT::t_float_1d d_vector_grid;     // Kokkos version of per-grid vector
  DAT::t_float_2d_lr d_array_grid;   // Kokkos version of per-grid array
  DAT::t_float_1d d_vector_particle;     // Kokkos version of per-particle vector
  DAT::t_float_2d_lr d_array_particle;   // Kokkos version of per-particle array

  DAT::tdual_float_2d_lr k_array;    // Kokkos DualView of global array

  // Region
  virtual void match_all_kokkos(DAT::tdual_int_1d) {}

  // flatten this region into a device-resident postfix (RPN) token stream,
  //   so a kernel can test a point against it without virtual dispatch and
  //   without the caller holding a typed copy of every region style.
  //   see region_prim_kokkos.h for the encoding.
  //   fills the DualView (already synced to device) and returns the number
  //   of tokens in it.  returns 0 if this region cannot be flattened.
  //   a composite region ignores the DualView it is handed and returns its
  //   own buffer, which stays valid until that same region flattens again;
  //   a primitive fills the DualView it is handed, growing it if needed.
  virtual int flatten_region_kokkos(tdual_region_token_1d &) {return 0;}

  KOKKOS_INLINE_FUNCTION
  int match_kokkos(double x, double y, double z) const {return 0;}
};

}

#endif

/* ERROR/WARNING messages:

*/

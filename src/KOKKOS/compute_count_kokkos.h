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

#ifdef COMPUTE_CLASS

ComputeStyle(count/kk,ComputeCountKokkos)

#else

#ifndef SPARTA_COMPUTE_COUNT_KOKKOS_H
#define SPARTA_COMPUTE_COUNT_KOKKOS_H

#include "compute_count.h"
#include "kokkos_base.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

template<int NEED_ATOMICS>
struct TagComputeCount_per_species_tally_atomic{};

class ComputeCountKokkos : public ComputeCount, public KokkosBase {
 public:
  ComputeCountKokkos(class SPARTA *, int, char **);
  ~ComputeCountKokkos();

  void init();
  double compute_scalar();
  void compute_vector();

  template<int NEED_ATOMICS>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeCount_per_species_tally_atomic<NEED_ATOMICS>, const int&) const;

 private:

  DAT::tdual_int_1d k_count;
  DAT::t_int_1d d_count;
  int need_dup;
  Kokkos::Experimental::ScatterView<int*, typename DAT::t_int_1d::array_layout,DeviceType,typename Kokkos::Experimental::ScatterSum,typename Kokkos::Experimental::ScatterDuplicated> dup_count;
  Kokkos::Experimental::ScatterView<int*, typename DAT::t_int_1d::array_layout,DeviceType,typename Kokkos::Experimental::ScatterSum,typename Kokkos::Experimental::ScatterNonDuplicated> ndup_count;

  t_particle_1d d_particles;

  void per_species_tally_kokkos();
};


}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

*/

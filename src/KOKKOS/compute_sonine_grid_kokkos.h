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

ComputeStyle(sonine/grid/kk,ComputeSonineGridKokkos)

#else

#ifndef SPARTA_COMPUTE_SONINE_GRID_KOKKOS_H
#define SPARTA_COMPUTE_SONINE_GRID_KOKKOS_H

#include "compute_sonine_grid.h"
#include "kokkos_base.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

// vcom calculations
template<int NEED_ATOMICS>
struct TagComputeSonineGrid_compute_vcom_init_atomic{};
struct TagComputeSonineGrid_compute_vcom{};
struct TagComputeSonineGrid_normalize_vcom{};

// compute_per_grid calculations
template<int NEED_ATOMICS>
struct TagComputeSonineGrid_compute_per_grid_atomic{};
struct TagComputeSonineGrid_compute_per_grid{};

struct TagComputeSonineGrid_post_process_grid{};


class ComputeSonineGridKokkos : public ComputeSonineGrid, public KokkosBase {
 public:
  ComputeSonineGridKokkos(class SPARTA *, int, char **);
  ~ComputeSonineGridKokkos();
  void compute_per_grid();
  void compute_per_grid_kokkos();
  int query_tally_grid_kokkos(DAT::t_float_2d_lr &);
  void post_process_grid_kokkos(int, int, DAT::t_float_2d_lr, int *,
                                  DAT::t_float_1d_strided);
  void reallocate();

  template<int NEED_ATOMICS>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeSonineGrid_compute_vcom_init_atomic<NEED_ATOMICS>, const int&) const;
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeSonineGrid_normalize_vcom, const int&) const;
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeSonineGrid_compute_vcom, const int&) const;

  template<int NEED_ATOMICS>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeSonineGrid_compute_per_grid_atomic<NEED_ATOMICS>, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeSonineGrid_compute_per_grid, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeSonineGrid_post_process_grid, const int&) const;

  DAT::tdual_float_1d k_vector_grid;

 private:
  DAT::t_float_3d d_vcom;
  Kokkos::Experimental::ScatterView<F_FLOAT***, typename DAT::t_float_3d::array_layout,DeviceType,typename Kokkos::Experimental::ScatterSum,typename Kokkos::Experimental::ScatterDuplicated> dup_vcom_tally;
  Kokkos::Experimental::ScatterView<F_FLOAT***, typename DAT::t_float_3d::array_layout,DeviceType,typename Kokkos::Experimental::ScatterSum,typename Kokkos::Experimental::ScatterNonDuplicated> ndup_vcom_tally;

  DAT::tdual_float_2d_lr k_tally;
  DAT::t_float_2d_lr d_tally;
  int need_dup;
  Kokkos::Experimental::ScatterView<F_FLOAT**, typename DAT::t_float_2d_lr::array_layout,DeviceType,typename Kokkos::Experimental::ScatterSum,typename Kokkos::Experimental::ScatterDuplicated> dup_tally;
  Kokkos::Experimental::ScatterView<F_FLOAT**, typename DAT::t_float_2d_lr::array_layout,DeviceType,typename Kokkos::Experimental::ScatterSum,typename Kokkos::Experimental::ScatterNonDuplicated> ndup_tally;

  DAT::t_float_2d_lr d_etally;
  DAT::t_float_1d_strided d_vec;

  t_cinfo_1d d_cinfo;
  t_particle_1d d_particles;
  t_species_1d d_species;
  DAT::t_int_2d d_s2g;

  DAT::t_int_1d d_cellcount;
  DAT::t_int_2d d_plist;

  DAT::tdual_int_1d k_which;
  DAT::tdual_int_1d k_moment;
  DAT::tdual_int_1d k_order;
  DAT::t_int_1d d_which;
  DAT::t_int_1d d_moment;
  DAT::t_int_1d d_order;

  int mass,numerator;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Compute sonine/grid mixture ID does not exist

Self-explanatory.

E: Number of groups in compute sonine/grid mixture has changed

This mixture property cannot be changed after this compute command is
issued.

E: Invalid call to ComputeSonineGrid::post_process_grid()

This indicates a coding error.  Please report the issue to the SPARTA
developers.

*/

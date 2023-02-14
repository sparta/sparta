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

#ifdef COMPUTE_CLASS

ComputeStyle(tvib/grid/kk,ComputeTvibGridKokkos)

#else

#ifndef SPARTA_COMPUTE_TVIB_GRID_KOKKOS_H
#define SPARTA_COMPUTE_TVIB_GRID_KOKKOS_H

#include "compute_tvib_grid.h"
#include "kokkos_base.h"
#include "kokkos_type.h"
#include "particle_kokkos.h"

namespace SPARTA_NS {

template<int NEED_ATOMICS>
struct TagComputeTvibGrid_compute_per_grid_atomic{};
struct TagComputeTvibGrid_compute_per_grid{};
struct TagComputeTvibGrid_post_process_grid{};

class ComputeTvibGridKokkos : public ComputeTvibGrid, public KokkosBase {
 public:
  ComputeTvibGridKokkos(class SPARTA *, int, char **);
  ~ComputeTvibGridKokkos();
  void compute_per_grid();
  void compute_per_grid_kokkos();
  int query_tally_grid_kokkos(DAT::t_float_2d_lr &);
  void post_process_grid_kokkos(int, int, DAT::t_float_2d_lr, int *,
                                DAT::t_float_1d_strided);
  void reallocate();

  template<int NEED_ATOMICS>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTvibGrid_compute_per_grid_atomic<NEED_ATOMICS>, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTvibGrid_compute_per_grid, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeTvibGrid_post_process_grid, const int&) const;

  DAT::tdual_float_1d k_vector_grid;

 private:
  int nstride,count,evib,nsp,imode;
  double boltz;

  DAT::tdual_float_2d_lr k_tally;
  DAT::t_float_2d_lr d_tally;

  int need_dup;
  Kokkos::Experimental::ScatterView<F_FLOAT**, typename DAT::t_float_2d_lr::array_layout,DeviceType,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterDuplicated> dup_tally;
  Kokkos::Experimental::ScatterView<F_FLOAT**, typename DAT::t_float_2d_lr::array_layout,DeviceType,Kokkos::Experimental::ScatterSum,Kokkos::Experimental::ScatterNonDuplicated> ndup_tally;

  DAT::t_float_2d_lr d_etally;
  DAT::t_float_1d_strided d_vec;

  t_cinfo_1d d_cinfo;
  t_particle_1d d_particles;
  t_species_1d d_species;
  DAT::t_int_2d d_s2g;

  DAT::t_int_1d d_cellcount;
  DAT::t_int_2d d_plist;

  DAT::t_float_1d d_tspecies;
  DAT::t_float_2d_lr d_tspecies_mode;

  DAT::tdual_int_1d k_s2t;
  DAT::tdual_int_1d k_t2s;
  DAT::tdual_int_1d k_t2s_mode;
  DAT::tdual_int_2d k_s2t_mode;

  DAT::t_int_1d d_s2t;
  DAT::t_int_1d d_t2s;
  DAT::t_int_1d d_t2s_mode;
  DAT::t_int_2d d_s2t_mode;

  DAT::t_int_1d d_ewhich;
  tdual_struct_tdual_int_2d_1d k_eiarray;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Compute grid mixture ID does not exist

Self-explanatory.

E: Number of groups in compute tvib/grid mixture has changed

This mixture property cannot be changed after this compute command is
issued.

E: Number of species in compute tvib/grid mixture has changed

This mixture property cannot be changed after this compute command is
issued.

*/

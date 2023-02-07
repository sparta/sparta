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

ComputeStyle(grid/kk,ComputeGridKokkos)

#else

#ifndef SPARTA_COMPUTE_GRID_KOKKOS_H
#define SPARTA_COMPUTE_GRID_KOKKOS_H

#include "compute_grid.h"
#include "kokkos_base.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

template<int NEED_ATOMICS>
struct TagComputeGrid_compute_per_grid_atomic{};

struct TagComputeGrid_compute_per_grid{};
struct TagComputeGrid_NUM{};
struct TagComputeGrid_MASS{};
struct TagComputeGrid_NRHO{};
struct TagComputeGrid_MASSRHO{};
struct TagComputeGrid_NFRAC{};
struct TagComputeGrid_U{};
struct TagComputeGrid_KE{};
struct TagComputeGrid_TEMPERATURE{};
struct TagComputeGrid_EROT{};
struct TagComputeGrid_TROT{};
struct TagComputeGrid_PXRHO{};
struct TagComputeGrid_KERHO{};

class ComputeGridKokkos : public ComputeGrid, public KokkosBase {
 public:
  ComputeGridKokkos(class SPARTA *, int, char **);
  ~ComputeGridKokkos();
  void compute_per_grid();
  void compute_per_grid_kokkos();
  int query_tally_grid_kokkos(DAT::t_float_2d_lr &);
  void post_process_grid_kokkos(int, int, DAT::t_float_2d_lr, int *,
                                  DAT::t_float_1d_strided);
  void reallocate();

  template<int NEED_ATOMICS>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeGrid_compute_per_grid_atomic<NEED_ATOMICS>, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeGrid_compute_per_grid, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeGrid_NUM, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeGrid_MASS, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeGrid_NRHO, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeGrid_MASSRHO, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeGrid_NFRAC, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeGrid_U, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeGrid_KE, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeGrid_TEMPERATURE, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeGrid_EROT, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeGrid_TROT, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeGrid_PXRHO, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeGrid_KERHO, const int&) const;

  DAT::tdual_float_1d k_vector_grid;

 private:
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

  DAT::tdual_int_1d k_unique;
  DAT::t_int_1d d_unique;

  double fnum;
  int count,mass,count_or_mass,cell_count_or_mass;
  int velocity,mvsq,eng,dof,mom,ke;
  int nsample,nstride;
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

E: Number of groups in compute grid mixture has changed

This mixture property cannot be changed after this compute command is
issued.

E: Invalid call to ComputeGridKokkos::post_process_grid()

This indicates a coding error.  Please report the issue to the SPARTA
developers.

*/

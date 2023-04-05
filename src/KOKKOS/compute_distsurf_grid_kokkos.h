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

ComputeStyle(distsurf/grid/kk,ComputeDistSurfGridKokkos)

#else

#ifndef SPARTA_COMPUTE_DISTSURF_GRID_KOKKOS_H
#define SPARTA_COMPUTE_DISTSURF_GRID_KOKKOS_H

#include "compute_distsurf_grid.h"
#include "kokkos_base.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

struct TagComputeDistSurfGrid_surf_centroid{};
struct TagComputeDistSurfGrid_surf_distance{};

 class ComputeDistSurfGridKokkos : public ComputeDistSurfGrid, public KokkosBase {
 public:
  ComputeDistSurfGridKokkos(class SPARTA *, int, char **);
  ~ComputeDistSurfGridKokkos();
  void compute_per_grid();
  void compute_per_grid_kokkos();
  void reallocate();

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeDistSurfGrid_surf_centroid, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputeDistSurfGrid_surf_distance, const int&) const;

  DAT::tdual_float_1d k_vector_grid;

 private:

  DAT::tdual_int_1d k_eflag;
  DAT::tdual_int_1d k_slist;
  DAT::t_int_1d d_eflag;
  DAT::t_int_1d d_slist;
  HAT::t_int_1d h_eflag;
  HAT::t_int_1d h_slist;
  DAT::t_float_1d_3 d_sctr;
  t_cinfo_1d d_cinfo;
  t_cell_1d d_cells;
  Kokkos::Crs<int, DeviceType, void, int> d_csurfs;
  Kokkos::Crs<int, DeviceType, void, int> d_csubs;

  t_line_1d d_lines;
  t_tri_1d d_tris;
  int dim;
  int nsurf;
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

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

#ifdef FIX_CLASS

FixStyle(ave/grid/kk,FixAveGridKokkos)

#else

#ifndef LMP_FIX_AVE_GRID_KOKKOS_H
#define LMP_FIX_AVE_GRID_KOKKOS_H

#include "fix_ave_grid.h"
#include "kokkos_type.h"
#include "kokkos_base.h"

namespace SPARTA_NS {

struct TagFixAveGrid_Zero_group_vector{};
struct TagFixAveGrid_Zero_group_array{};
struct TagFixAveGrid_Zero_tally{};
struct TagFixAveGrid_Add_ctally{};
struct TagFixAveGrid_Add_compute_vector{};
struct TagFixAveGrid_Add_compute_array{};
struct TagFixAveGrid_Add_fix_vector{};
struct TagFixAveGrid_Add_fix_array{};
struct TagFixAveGrid_Norm_vector_grid{};
struct TagFixAveGrid_Norm_array_grid{};

class FixAveGridKokkos : public FixAveGrid, public KokkosBase {
 public:
  FixAveGridKokkos(class SPARTA *, int, char **);
  ~FixAveGridKokkos();
  void init();
  void setup();
  void end_of_step();

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixAveGrid_Zero_group_vector, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixAveGrid_Zero_group_array, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixAveGrid_Zero_tally, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixAveGrid_Add_ctally, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixAveGrid_Add_compute_vector, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixAveGrid_Add_compute_array, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixAveGrid_Add_fix_vector, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixAveGrid_Add_fix_array, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixAveGrid_Norm_vector_grid, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixAveGrid_Norm_array_grid, const int&) const;

  DAT::tdual_float_1d k_vector_grid;

  DAT::tdual_float_2d_lr k_array_grid;

 private:
  DAT::tdual_float_2d_lr k_tally;
  DAT::t_float_2d_lr d_tally;

  DAT::t_float_2d_lr d_ctally;

  DAT::t_float_1d d_compute_vector;
  DAT::t_float_2d_lr d_compute_array;

  DAT::t_float_1d d_fix_vector;
  DAT::t_float_2d_lr d_fix_array;

  t_cinfo_1d d_cinfo;

  DAT::tdual_float_1d k_numap;
  DAT::t_float_1d d_numap;

  DAT::tdual_float_2d k_umap,k_uomap;
  DAT::t_float_2d d_umap,d_uomap;

  int j,k,kk,jm1,m,ntally;

  void grow_percell(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Compute ID for fix ave/grid does not exist

Self-explanatory.

E: Fix ID for fix ave/grid does not exist

Self-explanatory.

E: Fix ave/grid compute does not calculate per-grid values

Self-explanatory.

E: Fix ave/grid compute does not calculate a per-grid vector

Self-explanatory.

E: Fix ave/grid compute does not calculate a per-grid array

Self-explanatory.

E: Fix ave/grid compute array is accessed out-of-range

Self-explanatory.

E: Fix ave/grid fix does not calculate per-grid values

Self-explanatory.

E: Fix ave/grid fix does not calculate a per-grid vector

Self-explanatory.

E: Fix ave/grid fix does not calculate a per-grid array

Self-explanatory.

E: Fix ave/grid fix array is accessed out-of-range

Self-explanatory.

E: Fix for fix ave/grid not computed at compatible time

Fixes generate values on specific timesteps.  Fix ave/grid is
requesting a value on a non-allowed timestep.

E: Variable name for fix ave/grid does not exist

Self-explanatory.

E: Fix ave/grid variable is not grid-style variable

Self-explanatory.

*/

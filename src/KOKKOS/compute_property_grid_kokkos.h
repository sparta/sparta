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

ComputeStyle(property/grid/kk,ComputePropertyGridKokkos)

#else

#ifndef SPARTA_COMPUTE_PROPERTY_GRID_KOKKOS_H
#define SPARTA_COMPUTE_PROPERTY_GRID_KOKKOS_H

#include "compute_property_grid.h"
#include "kokkos_base.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

  struct TagComputePropertyGrid_ComputePerGrid_vector{};
  struct TagComputePropertyGrid_ComputePerGrid_array{};

 class ComputePropertyGridKokkos : public ComputePropertyGrid, public KokkosBase {
 public:
  ComputePropertyGridKokkos(class SPARTA *, int, char **);
  ~ComputePropertyGridKokkos();

  void compute_per_grid();
  void compute_per_grid_kokkos();
  void reallocate();

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputePropertyGrid_ComputePerGrid_vector, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagComputePropertyGrid_ComputePerGrid_array, const int&) const;

  DAT::tdual_float_1d k_vector_grid;
  DAT::tdual_float_2d_lr k_array_grid;

 private:
  t_cell_1d d_cells;
  t_cinfo_1d d_cinfo;
  DAT::tdual_int_1d k_index;
  DAT::t_int_1d d_index;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Invalid compute property/grid field for 2d simulation

Fields that reference z-dimension properties cannot be used
in a 2d simulation.

E: Invalid keyword in compute property/grid command

Self-explanatory.

*/

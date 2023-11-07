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

ComputeStyle(lambda/grid/kk,ComputeLambdaGridKokkos)

#else

#ifndef SPARTA_COMPUTE_LAMBDA_GRID_KOKKOS_H
#define SPARTA_COMPUTE_LAMBDA_GRID_KOKKOS_H

#include "compute_lambda_grid.h"
#include "kokkos_base.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

  struct TagComputeLambdaGrid_LoadNrhoVecFromArray{};
  struct TagComputeLambdaGrid_LoadTempVecFromArray{};
  struct TagComputeLambdaGrid_ComputePerGrid{};

  class ComputeLambdaGridKokkos : public ComputeLambdaGrid, public KokkosBase {
  public:
    ComputeLambdaGridKokkos(class SPARTA *, int, char **);
    ~ComputeLambdaGridKokkos();

    void compute_per_grid();
    void compute_per_grid_kokkos();
    void reallocate();

    KOKKOS_INLINE_FUNCTION
    void operator()(TagComputeLambdaGrid_LoadNrhoVecFromArray, const int&) const;

    KOKKOS_INLINE_FUNCTION
    void operator()(TagComputeLambdaGrid_LoadTempVecFromArray, const int&) const;

    KOKKOS_INLINE_FUNCTION
    void operator()(TagComputeLambdaGrid_ComputePerGrid, const int&) const;

    DAT::tdual_float_1d k_vector_grid;
    DAT::tdual_float_2d_lr k_array_grid;

  private:
    DAT::t_float_1d d_nrho_vector;
    DAT::t_float_1d d_temp_vector;

    int dimension;
    t_cell_1d d_cells;

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

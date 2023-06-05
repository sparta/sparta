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

ComputeStyle(dt/grid/kk,ComputeDtGridKokkos)

#else

#ifndef SPARTA_COMPUTE_DT_GRID_KOKKOS_H
#define SPARTA_COMPUTE_DT_GRID_KOKKOS_H

#include "compute_dt_grid.h"
#include "kokkos_base.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

  struct TagComputeDtGrid_LoadLambdaVecFromArray{};
  struct TagComputeDtGrid_LoadTempVecFromArray{};
  struct TagComputeDtGrid_LoadUsqVecFromArray{};
  struct TagComputeDtGrid_LoadVsqVecFromArray{};
  struct TagComputeDtGrid_LoadWsqVecFromArray{};
  struct TagComputeDtGrid_ComputePerGrid{};

  class ComputeDtGridKokkos : public ComputeDtGrid, public KokkosBase {
  public:
    ComputeDtGridKokkos(class SPARTA *, int, char **);
    ~ComputeDtGridKokkos();

    void compute_per_grid();
    void compute_per_grid_kokkos();
    void reallocate();

    KOKKOS_INLINE_FUNCTION void
    operator()(TagComputeDtGrid_LoadLambdaVecFromArray, const int&) const;

    KOKKOS_INLINE_FUNCTION void
    operator()(TagComputeDtGrid_LoadTempVecFromArray, const int&) const;

    KOKKOS_INLINE_FUNCTION void
    operator()(TagComputeDtGrid_LoadUsqVecFromArray, const int&) const;

    KOKKOS_INLINE_FUNCTION void
    operator()(TagComputeDtGrid_LoadVsqVecFromArray, const int&) const;

    KOKKOS_INLINE_FUNCTION void
    operator()(TagComputeDtGrid_LoadWsqVecFromArray, const int&) const;

    KOKKOS_INLINE_FUNCTION
    void operator()(TagComputeDtGrid_ComputePerGrid, const int&) const;

    DAT::tdual_float_1d k_vector_grid;

  private:
    t_cell_1d d_cells;
    t_cinfo_1d d_cinfo;

    DAT::t_float_1d d_lambda_vector;
    DAT::t_float_1d d_temp_vector;
    DAT::t_float_1d d_usq_vector;
    DAT::t_float_1d d_vsq_vector;
    DAT::t_float_1d d_wsq_vector;

    int dimension;
    double boltz;

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

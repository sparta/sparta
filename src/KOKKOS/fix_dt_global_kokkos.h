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

#ifdef FIX_CLASS

FixStyle(dt_global/kk,FixDtGlobalKokkos)

#else

#ifndef SPARTA_FIX_DT_GLOBAL_KOKKOS_H
#define SPARTA_FIX_DT_GLOBAL_KOKKOS_H

#include "fix_dt_global.h"
#include "grid_kokkos.h"
#include "kokkos_base.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

  struct TagFixDtGlobal_LoadLambdaVecFromArray{};
  struct TagFixDtGlobal_LoadTempVecFromArray{};
  struct TagFixDtGlobal_LoadUsqVecFromArray{};
  struct TagFixDtGlobal_LoadVsqVecFromArray{};
  struct TagFixDtGlobal_LoadWsqVecFromArray{};
  struct TagFixDtGlobal_LimitCellDtDesired{};

  template<int ATOMIC_REDUCTION>
  struct TagFixDtGlobal_SetCellDtDesired{};

  class FixDtGlobalKokkos : public FixDtGlobal, public KokkosBase {
  public:

    FixDtGlobalKokkos(class SPARTA *, int, char **);
    ~FixDtGlobalKokkos();
    //  void init();
    void end_of_step() override;
    void reallocate();

    KOKKOS_INLINE_FUNCTION void
    operator()(TagFixDtGlobal_LoadLambdaVecFromArray, const int&) const;

    KOKKOS_INLINE_FUNCTION void
    operator()(TagFixDtGlobal_LoadTempVecFromArray, const int&) const;

    KOKKOS_INLINE_FUNCTION void
    operator()(TagFixDtGlobal_LoadUsqVecFromArray, const int&) const;

    KOKKOS_INLINE_FUNCTION void
    operator()(TagFixDtGlobal_LoadVsqVecFromArray, const int&) const;

    KOKKOS_INLINE_FUNCTION void
    operator()(TagFixDtGlobal_LoadWsqVecFromArray, const int&) const;

    KOKKOS_INLINE_FUNCTION void
    operator()(TagFixDtGlobal_LimitCellDtDesired, const int&) const;

    template<int ATOMIC_REDUCTION>
    KOKKOS_INLINE_FUNCTION void
    operator()(TagFixDtGlobal_SetCellDtDesired<ATOMIC_REDUCTION>, const int&) const;

  protected:

    GridKokkos* grid_kk;

  private:
    t_cell_1d d_cells;

    DAT::tdual_float_scalar k_dtsum;
    DAT::t_float_scalar d_dtsum;
    HAT::t_float_scalar h_dtsum;

    DAT::tdual_float_scalar k_dtmin;
    DAT::t_float_scalar d_dtmin;
    HAT::t_float_scalar h_dtmin;

    DAT::tdual_float_scalar k_dtmax;
    DAT::t_float_scalar d_dtmax;
    HAT::t_float_scalar h_dtmax;

    DAT::t_float_1d d_lambda_vector;
    DAT::t_float_1d d_temp_vector;
    DAT::t_float_1d d_usq_vector;
    DAT::t_float_1d d_vsq_vector;
    DAT::t_float_1d d_wsq_vector;
    double dt_global;
    double boltz;
    //    double min_species_mass;
    int dimension;

};

}

#endif
#endif

/* ERROR/WARNING messages:
 */


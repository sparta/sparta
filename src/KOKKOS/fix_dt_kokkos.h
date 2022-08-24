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

FixStyle(dt/kk,FixDtKokkos)

#else

#ifndef SPARTA_FIX_DT_KOKKOS_H
#define SPARTA_FIX_DT_KOKKOS_H

#include "fix_dt.h"
#include "grid_kokkos.h"
#include "kokkos_base.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

  struct TagFixDt_LoadLambdaVecFromArray{};
  struct TagFixDt_LoadTempVecFromArray{};
  struct TagFixDt_LoadUsqVecFromArray{};
  struct TagFixDt_LoadVsqVecFromArray{};
  struct TagFixDt_LoadWsqVecFromArray{};

  template<int ATOMIC_REDUCTION>
  struct TagFixDt_SetCellDtDesired{};

  class FixDtKokkos : public FixDt, public KokkosBase {
  public:

    FixDtKokkos(class SPARTA *, int, char **);
    ~FixDtKokkos();
    //  void init();
    void end_of_step() override;
    void reallocate();

    KOKKOS_INLINE_FUNCTION void
    operator()(TagFixDt_LoadLambdaVecFromArray, const int&) const;

    KOKKOS_INLINE_FUNCTION void
    operator()(TagFixDt_LoadTempVecFromArray, const int&) const;

    KOKKOS_INLINE_FUNCTION void
    operator()(TagFixDt_LoadUsqVecFromArray, const int&) const;

    KOKKOS_INLINE_FUNCTION void
    operator()(TagFixDt_LoadVsqVecFromArray, const int&) const;

    KOKKOS_INLINE_FUNCTION void
    operator()(TagFixDt_LoadWsqVecFromArray, const int&) const;

    template<int ATOMIC_REDUCTION>
    KOKKOS_INLINE_FUNCTION void
    operator()(TagFixDt_SetCellDtDesired<ATOMIC_REDUCTION>, const int&) const;

  protected:

    GridKokkos* grid_kk;

  private:
    t_cell_1d d_cells;

    DAT::tdual_int_scalar k_ncells_with_a_particle;
    DAT::t_int_scalar d_ncells_with_a_particle;
    HAT::t_int_scalar h_ncells_with_a_particle;

    DAT::tdual_float_scalar k_dtsum;
    DAT::t_float_scalar d_dtsum;
    HAT::t_float_scalar h_dtsum;

    DAT::tdual_float_scalar k_dtmin;
    DAT::t_float_scalar d_dtmin;
    HAT::t_float_scalar h_dtmin;

    DAT::t_int_1d d_cellcount;

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


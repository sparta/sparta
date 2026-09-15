/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(dt/reset/kk,FixDtResetKokkos)

#else

#ifndef SPARTA_FIX_DT_RESET_KOKKOS_H
#define SPARTA_FIX_DT_RESET_KOKKOS_H

#include "fix_dt_reset.h"
#include "grid_kokkos.h"
#include "kokkos_base.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

  struct TagFixDtReset_LoadGridstepVecFromArray{};
  struct TagFixDtReset_MinMaxSumReduction{};

  class FixDtResetKokkos : public FixDtReset, public KokkosBase {
  public:

    FixDtResetKokkos(class SPARTA *, int, char **);
    ~FixDtResetKokkos();
    void end_of_step() override;

    KOKKOS_INLINE_FUNCTION
    void operator()(TagFixDtReset_LoadGridstepVecFromArray, const int&) const;

    KOKKOS_INLINE_FUNCTION
    void operator()(TagFixDtReset_MinMaxSumReduction, const int&, double&, double&, double&, int&) const;


  private:

    DAT::t_float_1d d_gridstep;

    // copy the first N cells of a source per-grid vector into d_gridstep
    // d_gridstep is sized to a high-water mark which the current cell count
    //   can be below, so the two views cannot be deep_copied whole

    void copy_gridstep(DAT::t_float_1d, int);

};

}

#endif
#endif

/* ERROR/WARNING messages:
 */


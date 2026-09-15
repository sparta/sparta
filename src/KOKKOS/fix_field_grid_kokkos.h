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

FixStyle(field/grid/kk,FixFieldGridKokkos)

#else

#ifndef SPARTA_FIX_FIELD_GRID_KOKKOS_H
#define SPARTA_FIX_FIELD_GRID_KOKKOS_H

#include "fix_field_grid.h"
#include "kokkos_base.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

// the grid-style variables this fix evaluates are host-only (see
//   VariableKokkos), so the evaluation itself stays on the host and only the
//   result is published to the device.  that is cheap here: UpdateKokkos
//   calls compute_field() once per fieldfreq steps, or once per run when
//   fieldfreq is 0, not once per timestep

class FixFieldGridKokkos : public FixFieldGrid, public KokkosBase {
 public:
  FixFieldGridKokkos(class SPARTA *, int, char **);
  ~FixFieldGridKokkos() override;
  void compute_field() override;

 private:
  DAT::tdual_float_2d_lr k_array_grid;
};

}

#endif
#endif

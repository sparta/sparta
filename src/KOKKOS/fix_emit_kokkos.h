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

#ifndef SPARTA_FIX_EMIT_KOKKOS_H
#define SPARTA_FIX_EMIT_KOKKOS_H

#include "kokkos_type.h"
#include "math_const.h"

namespace SPARTA_NS {

/* ----------------------------------------------------------------------
   device version of FixEmit::mol_inflow()
   calculate flux of particles of a species with vscale/fraction
     entering a grid cell
   shared by the Kokkos emit styles (fix emit/face/kk, fix emit/surf/kk)
   see comments for FixEmit::mol_inflow() in fix_emit.cpp
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
double mol_inflow_kokkos(double indot, double vscale, double fraction)
{
  const double scosine = indot / vscale;
  if (scosine < -3.0) return 0.0;
  const double inward_number_flux = vscale*fraction *
    (exp(-scosine*scosine) + MathConst::MY_PIS*scosine*(1.0 + erf(scosine))) /
    (2*MathConst::MY_PIS);
  return inward_number_flux;
}

}

#endif

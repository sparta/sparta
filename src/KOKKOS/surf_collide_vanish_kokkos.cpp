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

#include "math.h"
#include "surf_collide_vanish_kokkos.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

SurfCollideVanishKokkos::SurfCollideVanishKokkos(SPARTA *sparta, int narg, char **arg) :
  SurfCollideVanish(sparta, narg, arg)
{
  k_nsingle = DAT::tdual_int_scalar("SurfCollide:nsingle");
  d_nsingle = k_nsingle.d_view;
  h_nsingle = k_nsingle.h_view;

  allowreact = 0;
}

SurfCollideVanishKokkos::SurfCollideVanishKokkos(SPARTA *sparta) :
  SurfCollideVanish(sparta)
{
  id = NULL;
  style = NULL;
}

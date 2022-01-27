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

/* ----------------------------------------------------------------------
   Contributing author: Alan Stagg (Sandia)
------------------------------------------------------------------------- */

#include "surf_collide_piston_kokkos.h"
#include "update.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

SurfCollidePistonKokkos::SurfCollidePistonKokkos(SPARTA *sparta, int narg, char **arg) :
  SurfCollidePiston(sparta, narg, arg)
{
  k_nsingle = DAT::tdual_int_scalar("SurfCollide:nsingle");
  d_nsingle = k_nsingle.d_view;
  h_nsingle = k_nsingle.h_view;

  allowreact = 0;
}

/* ---------------------------------------------------------------------- */

SurfCollidePistonKokkos::SurfCollidePistonKokkos(SPARTA *sparta) :
  SurfCollidePiston(sparta)
{
  id = NULL;
  style = NULL;  
}

/* ---------------------------------------------------------------------- */

void SurfCollidePistonKokkos::init()
{
  SurfCollidePiston::init();
}

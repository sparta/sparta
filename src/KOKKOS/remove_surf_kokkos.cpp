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

#include "remove_surf_kokkos.h"
#include "surf_kokkos.h"
#include "sparta_masks.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

RemoveSurfKokkos::RemoveSurfKokkos(SPARTA *sparta) : RemoveSurf(sparta) {}

/* ---------------------------------------------------------------------- */

void RemoveSurfKokkos::command(int narg, char **arg)
{
  SurfKokkos* surf_kk = (SurfKokkos*) surf;
  surf_kk->sync(Host,ALL_MASK);

  RemoveSurf::command(narg,arg);

  surf_kk->modify(Host,ALL_MASK);
}

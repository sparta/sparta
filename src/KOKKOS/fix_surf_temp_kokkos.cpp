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

/* ----------------------------------------------------------------------
   Contributing author: Arnaud Borner (NASA Ames)
------------------------------------------------------------------------- */

#include "fix_surf_temp_kokkos.h"
#include "surf_kokkos.h"
#include "sparta_masks.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

FixSurfTempKokkos::FixSurfTempKokkos(SPARTA *sparta, int narg, char **arg) :
  FixSurfTemp(sparta, narg, arg)
{
  kokkos_flag = 1;
  execution_space = Host;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

FixSurfTempKokkos::~FixSurfTempKokkos()
{
}

/* ---------------------------------------------------------------------- */

void FixSurfTempKokkos::init()
{
  SurfKokkos* surf_kk = (SurfKokkos*) surf;
  surf_kk->sync(Host,SURF_CUSTOM_MASK);

  FixSurfTemp::init();

  surf_kk->modify(Host,SURF_CUSTOM_MASK);
}

/* ----------------------------------------------------------------------
   compute new surface element temperatures based on heat flux
   only invoked once every Nevery steps
------------------------------------------------------------------------- */

void FixSurfTempKokkos::end_of_step()
{
  SurfKokkos* surf_kk = (SurfKokkos*) surf;
  surf_kk->sync(Host,LINE_MASK|TRI_MASK|SURF_CUSTOM_MASK);

  FixSurfTemp::end_of_step();

  surf_kk->modify(Host,SURF_CUSTOM_MASK);
}

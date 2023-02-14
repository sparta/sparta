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

#include "string.h"
#include "fix_move_surf_kokkos.h"
#include "grid_kokkos.h"
#include "particle_kokkos.h"
#include "move_surf.h"
#include "update.h"
#include "comm.h"
#include "input.h"
#include "grid.h"
#include "surf_kokkos.h"
#include "domain.h"
#include "memory_kokkos.h"
#include "error.h"
#include "sparta_masks.h"

using namespace SPARTA_NS;

enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};           // several files

/* ---------------------------------------------------------------------- */

FixMoveSurfKokkos::FixMoveSurfKokkos(SPARTA *sparta, int narg, char **arg) :
  FixMoveSurf(sparta, narg, arg)
{
  kokkos_flag = 0; // need auto sync
  execution_space = Host;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ----------------------------------------------------------------------
   move surface points incrementally
------------------------------------------------------------------------- */

void FixMoveSurfKokkos::end_of_step()
{
  GridKokkos* grid_kk = (GridKokkos*) grid;
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  SurfKokkos* surf_kk = (SurfKokkos*) surf;

  grid_kk->sync(Host,ALL_MASK);
  particle_kk->sync(Host,ALL_MASK);
  surf_kk->sync(Host,ALL_MASK);

  FixMoveSurf::end_of_step();

  grid_kk->modify(Host,ALL_MASK);
  particle_kk->modify(Host,ALL_MASK);
  surf_kk->modify(Host,ALL_MASK);
  particle_kk->sorted_kk = 0;

  grid_kk->wrap_kokkos_graphs();
  grid_kk->update_hash();
}

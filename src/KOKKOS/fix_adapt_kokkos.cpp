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
#include "fix_adapt_kokkos.h"
#include "adapt_grid.h"
#include "grid_kokkos.h"
#include "particle_kokkos.h"
#include "surf.h"
#include "comm.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "input.h"
#include "output.h"
#include "dump.h"
#include "error.h"
#include "sparta_masks.h"

using namespace SPARTA_NS;

enum{NONE,REFINE,COARSEN};              // also in AdaptGrid

/* ---------------------------------------------------------------------- */

FixAdaptKokkos::FixAdaptKokkos(SPARTA *sparta, int narg, char **arg) :
  FixAdapt(sparta, narg, arg)
{
  kokkos_flag = 0; // need auto sync
  execution_space = Host;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

}

/* ----------------------------------------------------------------------
 *    perform dynamic load balancing
 *    ------------------------------------------------------------------------- */

void FixAdaptKokkos::end_of_step()
{
  GridKokkos* grid_kk = (GridKokkos*) grid;
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;

  grid_kk->sync(Host,ALL_MASK);
  particle_kk->sync(Host,PARTICLE_MASK);

  FixAdapt::end_of_step();

  grid_kk->modify(Host,ALL_MASK);
  particle_kk->modify(Host,PARTICLE_MASK);
  particle_kk->sorted_kk = 0;

  grid_kk->wrap_kokkos_graphs();
  grid_kk->update_hash();
}


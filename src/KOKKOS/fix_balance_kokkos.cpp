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
#include "stdlib.h"
#include "fix_balance_kokkos.h"
#include "balance_grid.h"
#include "update.h"
#include "grid_kokkos.h"
#include "particle_kokkos.h"
#include "surf_kokkos.h"
#include "comm.h"
#include "rcb.h"
#include "modify.h"
#include "compute.h"
#include "output.h"
#include "dump.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "memory_kokkos.h"
#include "error.h"
#include "sparta_masks.h"

using namespace SPARTA_NS;

enum{RANDOM,PROC,BISECTION};
enum{CELL,PARTICLE};

#define ZEROPARTICLE 0.1

/* ---------------------------------------------------------------------- */

FixBalanceKokkos::FixBalanceKokkos(SPARTA *sparta, int narg, char **arg) :
  FixBalance(sparta, narg, arg)
{
  kokkos_flag = 0; // need auto sync
  execution_space = Host;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ----------------------------------------------------------------------
   perform dynamic load balancing
------------------------------------------------------------------------- */

void FixBalanceKokkos::end_of_step()
{
  GridKokkos* grid_kk = (GridKokkos*) grid;
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  SurfKokkos* surf_kk = (SurfKokkos*) surf;

  grid_kk->sync(Host,CELL_MASK|CINFO_MASK|SINFO_MASK);
  particle_kk->sync(Host,PARTICLE_MASK);
  surf_kk->sync(Host,ALL_MASK);

  FixBalance::end_of_step();

  grid_kk->modify(Host,CELL_MASK|CINFO_MASK|SINFO_MASK);
  particle_kk->modify(Host,PARTICLE_MASK);
  surf_kk->modify(Host,ALL_MASK);
  particle_kk->sorted_kk = 0;

  grid_kk->wrap_kokkos_graphs();
  grid_kk->update_hash();
}


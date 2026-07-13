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

#include "compute_gas_collision_grid_kokkos.h"
#include "particle_kokkos.h"
#include "grid_kokkos.h"
#include "memory_kokkos.h"
#include "sparta_masks.h"
#include "kokkos.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

ComputeGasCollisionGridKokkos::ComputeGasCollisionGridKokkos(SPARTA *sparta, int narg, char **arg) :
  ComputeGasCollisionGrid(sparta, narg, arg)
{
  kokkos_flag = 1;
}

/* ---------------------------------------------------------------------- */

ComputeGasCollisionGridKokkos::ComputeGasCollisionGridKokkos(SPARTA *sparta) :
  ComputeGasCollisionGrid(sparta)
{
  copy = 1;
  uncopy = 0;
  vector_grid = NULL;
  nglocal = 0;
}

/* ---------------------------------------------------------------------- */

ComputeGasCollisionGridKokkos::~ComputeGasCollisionGridKokkos()
{
  if (copy || copymode) return;

  memoryKK->destroy_kokkos(k_vector_grid,vector_grid);
  vector_grid = NULL;
}

/* ----------------------------------------------------------------------
   zero the tally array on device
   called by Update at beginning of timesteps gas tallying is done
------------------------------------------------------------------------- */

void ComputeGasCollisionGridKokkos::clear()
{
  Kokkos::deep_copy(d_vector_grid,0.0);
}

/* ----------------------------------------------------------------------
   setup device views and scatter view before gas tallying
   called by Collide before the collision kernel
------------------------------------------------------------------------- */

void ComputeGasCollisionGridKokkos::pre_gas_tally()
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  d_s2g = particle_kk->k_species2group.view_device();

  GridKokkos* grid_kk = (GridKokkos*) grid;
  grid_kk->sync(Device,CINFO_MASK);
  d_cinfo = grid_kk->k_cinfo.view_device();
}

/* ----------------------------------------------------------------------
   finalize gas tallying, sync result to host
   called by Collide after the collision kernel
------------------------------------------------------------------------- */

void ComputeGasCollisionGridKokkos::post_gas_tally()
{
  k_vector_grid.modify_device();
  k_vector_grid.sync_host();
}

/* ----------------------------------------------------------------------
   reallocate data storage if nglocal has changed
   called by init() and whenever grid changes
------------------------------------------------------------------------- */

void ComputeGasCollisionGridKokkos::reallocate()
{
  if (grid->nlocal == nglocal) return;

  memoryKK->destroy_kokkos(k_vector_grid,vector_grid);
  nglocal = grid->nlocal;
  memoryKK->create_kokkos(k_vector_grid,vector_grid,nglocal,"gas/collision/grid:vector_grid");
  d_vector_grid = k_vector_grid.view_device();

  // clear counts b/c may be accessed before tallying is done

  Kokkos::deep_copy(d_vector_grid,0.0);
  k_vector_grid.modify_device();
  k_vector_grid.sync_host();
}

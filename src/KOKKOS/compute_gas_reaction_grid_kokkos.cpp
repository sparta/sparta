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

#include "compute_gas_reaction_grid_kokkos.h"
#include "particle_kokkos.h"
#include "grid_kokkos.h"
#include "react.h"
#include "memory_kokkos.h"
#include "sparta_masks.h"
#include "kokkos.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

ComputeGasReactionGridKokkos::ComputeGasReactionGridKokkos(SPARTA *sparta, int narg, char **arg) :
  ComputeGasReactionGrid(sparta, narg, arg)
{
  kokkos_flag = 1;
}

/* ---------------------------------------------------------------------- */

ComputeGasReactionGridKokkos::ComputeGasReactionGridKokkos(SPARTA *sparta) :
  ComputeGasReactionGrid(sparta)
{
  copy = 1;
  uncopy = 0;
  vector_grid = NULL;
  array_grid = NULL;
  ncol = 0;
  nglocal = 0;
}

/* ---------------------------------------------------------------------- */

ComputeGasReactionGridKokkos::~ComputeGasReactionGridKokkos()
{
  if (copy || copymode) return;

  if (ncol == 0) memoryKK->destroy_kokkos(k_vector_grid,vector_grid);
  else memoryKK->destroy_kokkos(k_array_grid,array_grid);
  vector_grid = NULL;
  array_grid = NULL;
}

/* ---------------------------------------------------------------------- */

void ComputeGasReactionGridKokkos::init()
{
  ComputeGasReactionGrid::init();

  // device copy of reaction -> column map for SELECT mode

  if (mode == SELECT) {
    int n = react->nlist + 1;
    d_reaction2col = DAT::t_int_1d("gas/reaction/grid:reaction2col",n);
    auto h_reaction2col = Kokkos::create_mirror_view(d_reaction2col);
    for (int i = 0; i < n; i++) h_reaction2col(i) = reaction2col[i];
    Kokkos::deep_copy(d_reaction2col,h_reaction2col);
  }
}

/* ----------------------------------------------------------------------
   zero the tally array on device
   called by Update at beginning of timesteps gas tallying is done
------------------------------------------------------------------------- */

void ComputeGasReactionGridKokkos::clear()
{
  if (ncol == 0) Kokkos::deep_copy(d_vector_grid,0.0);
  else Kokkos::deep_copy(d_array_grid,0.0);
}

/* ----------------------------------------------------------------------
   setup device views and scatter view before gas tallying
   called by Collide before the collision kernel
------------------------------------------------------------------------- */

void ComputeGasReactionGridKokkos::pre_gas_tally()
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

void ComputeGasReactionGridKokkos::post_gas_tally()
{
  if (ncol == 0) {
    k_vector_grid.modify_device();
    k_vector_grid.sync_host();
  } else {
    k_array_grid.modify_device();
    k_array_grid.sync_host();
  }
}

/* ----------------------------------------------------------------------
   reallocate data storage if nglocal has changed
   called by init() and whenever grid changes
------------------------------------------------------------------------- */

void ComputeGasReactionGridKokkos::reallocate()
{
  if (grid->nlocal == nglocal) return;

  if (ncol == 0) memoryKK->destroy_kokkos(k_vector_grid,vector_grid);
  else memoryKK->destroy_kokkos(k_array_grid,array_grid);
  nglocal = grid->nlocal;
  if (ncol == 0) {
    memoryKK->create_kokkos(k_vector_grid,vector_grid,nglocal,"gas/reaction/grid:vector_grid");
    d_vector_grid = k_vector_grid.view_device();
    Kokkos::deep_copy(d_vector_grid,0.0);
    k_vector_grid.modify_device();
    k_vector_grid.sync_host();
  } else {
    memoryKK->create_kokkos(k_array_grid,array_grid,nglocal,ncol,"gas/reaction/grid:array_grid");
    d_array_grid = k_array_grid.view_device();
    Kokkos::deep_copy(d_array_grid,0.0);
    k_array_grid.modify_device();
    k_array_grid.sync_host();
  }
}

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
   Contributing author: Alan Stagg (Sandia)
------------------------------------------------------------------------- */

#include "string.h"
#include "compute_lambda_grid_kokkos.h"
#include "update.h"
#include "grid_kokkos.h"
#include "domain.h"
#include "collide.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "math_const.h"
#include "memory.h"
#include "memory_kokkos.h"
#include "error.h"
#include "kokkos.h"
#include "sparta_masks.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{NONE,COMPUTE,FIX};
enum{KNONE,KALL,KX,KY,KZ};

#define INVOKED_PER_GRID 16
#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

ComputeLambdaGridKokkos::ComputeLambdaGridKokkos(SPARTA *sparta, int narg, char **arg) :
  ComputeLambdaGrid(sparta, narg, arg)
{
  kokkos_flag = 1;
}

/* ---------------------------------------------------------------------- */

ComputeLambdaGridKokkos::~ComputeLambdaGridKokkos()
{
  if (copymode) return;

  if (kflag == KNONE)
    memoryKK->destroy_kokkos(k_vector_grid,vector_grid);
  else
    memoryKK->destroy_kokkos(k_array_grid,array_grid);
  vector_grid = NULL;
  array_grid = NULL;
}

/* ---------------------------------------------------------------------- */

void ComputeLambdaGridKokkos::compute_per_grid()
{
  if (sparta->kokkos->prewrap)
    ComputeLambdaGrid::compute_per_grid();
  else
    compute_per_grid_kokkos();
}

/* ---------------------------------------------------------------------- */

void ComputeLambdaGridKokkos::compute_per_grid_kokkos()
{
  invoked_per_grid = update->ntimestep;

  if (nrhowhich == FIX && update->ntimestep % fnrho->per_grid_freq)
    error->all(FLERR,"Compute lambda/grid fix not computed at compatible time");
  if (tempwhich == FIX && update->ntimestep % ftemp->per_grid_freq)
    error->all(FLERR,"Compute lambda/grid fix not computed at compatible time");

  // grab nrho and temp values from compute or fix
  // invoke nrho and temp computes as needed

  if (nrhowhich == COMPUTE) {
    if (!cnrho->kokkos_flag)
      error->all(FLERR,"Cannot (yet) use non-Kokkos computes with compute lambda/grid/kk");
    KokkosBase* computeKKBase = dynamic_cast<KokkosBase*>(cnrho);
    if (!(cnrho->invoked_flag & INVOKED_PER_GRID)) {
      computeKKBase->compute_per_grid_kokkos();
      cnrho->invoked_flag |= INVOKED_PER_GRID;
    }

    if (cnrho->post_process_grid_flag)
      computeKKBase->post_process_grid_kokkos(nrhoindex,1,DAT::t_float_2d_lr(),NULL,DAT::t_float_1d_strided());

    if (nrhoindex == 0 || cnrho->post_process_grid_flag)
      Kokkos::deep_copy(d_nrho_vector, computeKKBase->d_vector);
    else {
      d_array = computeKKBase->d_array_grid;
      copymode = 1;
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeLambdaGrid_LoadNrhoVecFromArray>(0,nglocal),*this);
      copymode = 0;
    }
  } else if (nrhowhich == FIX) {
    if (!fnrho->kokkos_flag)
      error->all(FLERR,"Cannot (yet) use non-Kokkos fixes with compute lambda/grid/kk");
    KokkosBase* computeKKBase = dynamic_cast<KokkosBase*>(fnrho);
    if (nrhoindex == 0)
      d_nrho_vector = computeKKBase->d_vector;
    else {
      d_array = computeKKBase->d_array_grid;
      copymode = 1;
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeLambdaGrid_LoadNrhoVecFromArray>(0,nglocal),*this);
      copymode = 0;
    }
  }

  if (tempwhich == COMPUTE) {
    if (!ctemp->kokkos_flag)
      error->all(FLERR,"Cannot (yet) use non-Kokkos computes with compute lambda/grid/kk");
    KokkosBase* computeKKBase = dynamic_cast<KokkosBase*>(ctemp);
    if (!(ctemp->invoked_flag & INVOKED_PER_GRID)) {
      computeKKBase->compute_per_grid_kokkos();
      ctemp->invoked_flag |= INVOKED_PER_GRID;
    }

    if (ctemp->post_process_grid_flag)
      computeKKBase->post_process_grid_kokkos(tempindex,1,DAT::t_float_2d_lr(),NULL,DAT::t_float_1d_strided());

    if (tempindex == 0 || ctemp->post_process_grid_flag)
      Kokkos::deep_copy(d_temp_vector, computeKKBase->d_vector);
    else{
      d_array = computeKKBase->d_array_grid;
      copymode = 1;
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeLambdaGrid_LoadTempVecFromArray>(0,nglocal),*this);
      copymode = 0;
    }
  } else if (tempwhich == FIX) {
    if (!ftemp->kokkos_flag)
      error->all(FLERR,"Cannot (yet) use non-Kokkos fixes with compute lambda/grid/kk");
    KokkosBase* computeKKBase = dynamic_cast<KokkosBase*>(ftemp);
    if (tempindex == 0)
      d_temp_vector = computeKKBase->d_vector;
    else {
      d_array = computeKKBase->d_array_grid;
      copymode = 1;
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeLambdaGrid_LoadTempVecFromArray>(0,nglocal),*this);
      copymode = 0;
    }
  }

  GridKokkos* grid_kk = ((GridKokkos*)grid);
  grid_kk->sync(Device,CELL_MASK);
  d_cells = grid_kk->k_cells.d_view;
  dimension = domain->dimension;
  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeLambdaGrid_ComputePerGrid>(0,nglocal),*this);
  copymode = 0;

  if (kflag == KNONE) {
    k_vector_grid.modify_device();
    k_vector_grid.sync_host();
  } else {
    k_array_grid.modify_device();
    k_array_grid.sync_host();
  }
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeLambdaGridKokkos::operator()(TagComputeLambdaGrid_LoadNrhoVecFromArray, const int &i) const {
  d_nrho_vector(i) = d_array(i,nrhoindex-1);
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeLambdaGridKokkos::operator()(TagComputeLambdaGrid_LoadTempVecFromArray, const int &i) const {
  d_temp_vector(i) = d_array(i,tempindex-1);
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeLambdaGridKokkos::operator()(TagComputeLambdaGrid_ComputePerGrid, const int &i) const {
  double lambda;

  if (d_nrho_vector(i) == 0.0) lambda = BIG;
  else if (tempwhich == NONE || d_temp_vector(i) == 0.0)
    lambda = 1.0 / (prefactor * d_nrho_vector(i));
  else
    lambda = 1.0 / (prefactor * d_nrho_vector(i) * pow(tref/d_temp_vector(i),omega-0.5));

  if (kflag == KNONE) d_vector(i) = lambda;
  else d_array_grid(i,0) = lambda;

  // calculate per-cell Knudsen number

  if (kflag == KNONE) return;

  if (kflag == KALL) {
    double size;
    size =  (d_cells(i).hi[0] - d_cells(i).lo[0]);
    size += (d_cells(i).hi[1] - d_cells(i).lo[1]);
    if (dimension == 2) size *= 0.5;
    else {
      size += (d_cells(i).hi[2] - d_cells(i).lo[2]);
      size /= 3.0;
    }
    d_array_grid(i,1) = d_array_grid(i,0) / size;
  } else if (kflag == KX) {
    d_array_grid(i,1) = d_array_grid(i,0) / (d_cells(i).hi[0] - d_cells(i).lo[0]);
  } else if (kflag == KY) {
    d_array_grid(i,1) = d_array_grid(i,0) / (d_cells(i).hi[1] - d_cells(i).lo[1]);
  } else if (kflag == KZ) {
    d_array_grid(i,1) = d_array_grid(i,0) / (d_cells(i).hi[2] - d_cells(i).lo[2]);
  }
}

/* ----------------------------------------------------------------------
   reallocate arrays if nglocal has changed
   called by init() and load balancer
------------------------------------------------------------------------- */

void ComputeLambdaGridKokkos::reallocate()
{
  if (grid->nlocal == nglocal) return;

  nglocal = grid->nlocal;
  if (kflag == KNONE) {
    memoryKK->destroy_kokkos(k_vector_grid,vector_grid);
    memoryKK->create_kokkos(k_vector_grid,vector_grid,nglocal,"lambda/grid:vector_grid");
    d_vector = k_vector_grid.d_view;
  } else {
    memoryKK->destroy_kokkos(k_array_grid,array_grid);
    memoryKK->create_kokkos(k_array_grid,array_grid,nglocal,2,"lambda/grid:array_grid");
    d_array_grid = k_array_grid.d_view;
  }
  d_nrho_vector = DAT::t_float_1d ("d_nrho_vector", nglocal);
  if (tempwhich != NONE)
    d_temp_vector = DAT::t_float_1d ("d_temp_vector", nglocal);
}

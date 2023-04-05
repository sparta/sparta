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
#include "compute_property_grid_kokkos.h"
#include "grid_kokkos.h"
#include "update.h"
#include "domain.h"
#include "memory_kokkos.h"
#include "error.h"
#include "kokkos.h"
#include "sparta_masks.h"

using namespace SPARTA_NS;

enum{ID,PROC,XLO,YLO,ZLO,XHI,YHI,ZHI,XC,YC,ZC,VOL};

/* ---------------------------------------------------------------------- */

ComputePropertyGridKokkos::ComputePropertyGridKokkos(SPARTA *sparta, int narg, char **arg) :
  ComputePropertyGrid(sparta, narg, arg)
{
  kokkos_flag = 1;
  k_index = DAT::tdual_int_1d("compute/property/grid:index",nvalues);
  for (int n = 0; n < nvalues; n++){
    k_index.h_view(n) = index[n];
  }
  k_index.modify_host();
  k_index.sync_device();
  d_index = k_index.d_view;
}

/* ---------------------------------------------------------------------- */

ComputePropertyGridKokkos::~ComputePropertyGridKokkos()
{
  if (copymode) return;
  if (nvalues == 1)
    memoryKK->destroy_kokkos(k_vector_grid,vector_grid);
  else
    memoryKK->destroy_kokkos(k_array_grid,array_grid);
  vector_grid = NULL;
  array_grid = NULL;
}

/* ---------------------------------------------------------------------- */

void ComputePropertyGridKokkos::compute_per_grid()
{
  if (sparta->kokkos->prewrap) {
    ComputePropertyGrid::compute_per_grid();
  } else {
    compute_per_grid_kokkos();
    if (nvalues == 1) {
      k_vector_grid.modify_device();
      k_vector_grid.sync_host();
    } else {
      k_array_grid.modify_device();
      k_array_grid.sync_host();
    }
  }
}

/* ---------------------------------------------------------------------- */
void ComputePropertyGridKokkos::compute_per_grid_kokkos()
{
  GridKokkos* grid_kk = ((GridKokkos*)grid);
  d_cells = grid_kk->k_cells.d_view;
  d_cinfo = grid_kk->k_cinfo.d_view;
  grid_kk->sync(Device,CELL_MASK|CINFO_MASK);

  copymode = 1;
  if (nvalues == 1)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputePropertyGrid_ComputePerGrid_vector>(0,nglocal),*this);
  else
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputePropertyGrid_ComputePerGrid_array>(0,nglocal),*this);
  copymode = 0;
}

/* ---------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void ComputePropertyGridKokkos::operator()(TagComputePropertyGrid_ComputePerGrid_vector, const int &i) const
{
  switch (d_index[0]) {
  case ID:
    if (d_cinfo[i].mask & groupbit) d_vector[i] = d_cells[i].id;
    else d_vector[i] = 0.0;
    break;
  case PROC:
    if (d_cinfo[i].mask & groupbit) d_vector[i] = d_cells[i].proc;
    else d_vector[i] = 0.0;
    break;
  case XLO:
    if (d_cinfo[i].mask & groupbit) d_vector[i] = d_cells[i].lo[0];
    else d_vector[i] = 0.0;
    break;
  case YLO:
    if (d_cinfo[i].mask & groupbit) d_vector[i] = d_cells[i].lo[1];
    else d_vector[i] = 0.0;
    break;
  case ZLO:
    if (d_cinfo[i].mask & groupbit) d_vector[i] = d_cells[i].lo[2];
    else d_vector[i] = 0.0;
    break;
  case XHI:
    if (d_cinfo[i].mask & groupbit) d_vector[i] = d_cells[i].hi[0];
    else d_vector[i] = 0.0;
    break;
  case YHI:
    if (d_cinfo[i].mask & groupbit) d_vector[i] = d_cells[i].hi[1];
    else d_vector[i] = 0.0;
    break;
  case ZHI:
    if (d_cinfo[i].mask & groupbit) d_vector[i] = d_cells[i].hi[2];
    else d_vector[i] = 0.0;
    break;
  case XC:
    if (d_cinfo[i].mask & groupbit)
      d_vector[i] = 0.5 * (d_cells[i].lo[0] + d_cells[i].hi[0]);
    else d_vector[i] = 0.0;
    break;
  case YC:
    if (d_cinfo[i].mask & groupbit)
      d_vector[i] = 0.5 * (d_cells[i].lo[1] + d_cells[i].hi[1]);
    else d_vector[i] = 0.0;
    break;
  case ZC:
    if (d_cinfo[i].mask & groupbit)
      d_vector[i] = 0.5 * (d_cells[i].lo[2] + d_cells[i].hi[2]);
    else d_vector[i] = 0.0;
    break;
  case VOL:
    if (d_cinfo[i].mask & groupbit) d_vector[i] = d_cinfo[i].volume;
    else d_vector[i] = 0.0;
    break;
  }
}

/* ---------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void ComputePropertyGridKokkos::operator()(TagComputePropertyGrid_ComputePerGrid_array, const int &i) const
{
  for (int n=0; n<nvalues; ++n) {
    switch (d_index[n]) {
    case ID:
      if (d_cinfo[i].mask & groupbit) d_array_grid(i,n) = d_cells[i].id;
      else d_array_grid(i,n) = 0.0;
      break;
    case PROC:
      if (d_cinfo[i].mask & groupbit) d_array_grid(i,n) = d_cells[i].proc;
      else d_array_grid(i,n) = 0.0;
      break;
    case XLO:
      if (d_cinfo[i].mask & groupbit) d_array_grid(i,n) = d_cells[i].lo[0];
      else d_array_grid(i,n) = 0.0;
      break;
    case YLO:
      if (d_cinfo[i].mask & groupbit) d_array_grid(i,n) = d_cells[i].lo[1];
      else d_array_grid(i,n) = 0.0;
      break;
    case ZLO:
      if (d_cinfo[i].mask & groupbit) d_array_grid(i,n) = d_cells[i].lo[2];
      else d_array_grid(i,n) = 0.0;
      break;
    case XHI:
      if (d_cinfo[i].mask & groupbit) d_array_grid(i,n) = d_cells[i].hi[0];
      else d_array_grid(i,n) = 0.0;
      break;
    case YHI:
      if (d_cinfo[i].mask & groupbit) d_array_grid(i,n) = d_cells[i].hi[1];
      else d_array_grid(i,n) = 0.0;
      break;
    case ZHI:
      if (d_cinfo[i].mask & groupbit) d_array_grid(i,n) = d_cells[i].hi[2];
      else d_array_grid(i,n) = 0.0;
      break;
    case XC:
      if (d_cinfo[i].mask & groupbit)
        d_array_grid(i,n) = 0.5 * (d_cells[i].lo[0] + d_cells[i].hi[0]);
      else d_array_grid(i,n) = 0.0;
      break;
    case YC:
      if (d_cinfo[i].mask & groupbit)
        d_array_grid(i,n) = 0.5 * (d_cells[i].lo[1] + d_cells[i].hi[1]);
      else d_array_grid(i,n) = 0.0;
      break;
    case ZC:
      if (d_cinfo[i].mask & groupbit)
        d_array_grid(i,n) = 0.5 * (d_cells[i].lo[2] + d_cells[i].hi[2]);
      else d_array_grid(i,n) = 0.0;
      break;
    case VOL:
      if (d_cinfo[i].mask & groupbit) d_array_grid(i,n) = d_cinfo[i].volume;
      else d_array_grid(i,n) = 0.0;
      break;
    }
  }
}

/* ----------------------------------------------------------------------
   reallocate arrays if nglocal has changed
   called by init() and load balancer
------------------------------------------------------------------------- */

void ComputePropertyGridKokkos::reallocate()
{
  if (grid->nlocal == nglocal) return;
  nglocal = grid->nlocal;
  if (nvalues == 1) {
    memoryKK->destroy_kokkos(k_vector_grid,vector_grid);
    memoryKK->create_kokkos(k_vector_grid,vector_grid,nglocal,"property/grid:vector_grid");
    d_vector = k_vector_grid.d_view;
  } else {
    memoryKK->destroy_kokkos(k_array_grid,array_grid);
    memoryKK->create_kokkos(k_array_grid,array_grid,nglocal,nvalues,"property/grid:array_grid");
    d_array_grid = k_array_grid.d_view;
  }
}

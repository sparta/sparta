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

/* ----------------------------------------------------------------------
   Contributing author: Alan Stagg (Sandia)
------------------------------------------------------------------------- */

#include "compute_dt_grid_kokkos.h"
#include "update.h"
#include "grid_kokkos.h"
#include "domain.h"
#include "fix.h"
#include "compute.h"
#include "memory.h"
#include "memory_kokkos.h"
#include "error.h"
#include "kokkos.h"
#include "sparta_masks.h"

using namespace SPARTA_NS;

enum{NONE,COMPUTE,FIX};

#define INVOKED_PER_GRID 16
#define BIG 1.0e20
#define EPSILON 1.0e-12

/* ---------------------------------------------------------------------- */

ComputeDtGridKokkos::ComputeDtGridKokkos(SPARTA *sparta, int narg, char **arg) :
  ComputeDtGrid(sparta, narg, arg)
{
  kokkos_flag = 1;
  dimension = domain->dimension;
  boltz = update->boltz;
}

/* ---------------------------------------------------------------------- */

ComputeDtGridKokkos::~ComputeDtGridKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_vector_grid,vector_grid);
  vector_grid = NULL;
}

/* ---------------------------------------------------------------------- */

void ComputeDtGridKokkos::compute_per_grid()
{
  if (sparta->kokkos->prewrap)
    ComputeDtGrid::compute_per_grid();
  else
    compute_per_grid_kokkos();
}

/* ---------------------------------------------------------------------- */

void ComputeDtGridKokkos::compute_per_grid_kokkos()
{
  invoked_per_grid = update->ntimestep;

  if (tau_which == FIX && update->ntimestep % ftau->per_grid_freq)
    error->all(FLERR,"Compute dt/grid tau fix not computed at compatible time");
  if (temp_which == FIX && update->ntimestep % ftemp->per_grid_freq)
    error->all(FLERR,"Compute dt/grid temp fix not computed at compatible time");
  if (usq_which == FIX && update->ntimestep % fusq->per_grid_freq)
    error->all(FLERR,"Compute dt/grid usq fix not computed at compatible time");
  if (vsq_which == FIX && update->ntimestep % fvsq->per_grid_freq)
    error->all(FLERR,"Compute dt/grid vsq fix not computed at compatible time");
  if (wsq_which == FIX && update->ntimestep % fwsq->per_grid_freq)
    error->all(FLERR,"Compute dt/grid wsq fix not computed at compatible time");

  // grab per grid cell tau from compute or fix, invoke compute if needed
  // ditto for temp,us,vsq,wsq
  if (tau_which == COMPUTE) {
    if (!ctau->kokkos_flag)
      error->all(FLERR,"Cannot (yet) use non-Kokkos computes with compute dt/grid/kk");
    KokkosBase* computeKKBase = dynamic_cast<KokkosBase*>(ctau);
    if (!(ctau->invoked_flag & INVOKED_PER_GRID)) {
      computeKKBase->compute_per_grid_kokkos();
      ctau->invoked_flag |= INVOKED_PER_GRID;
    }
    if (ctau->post_process_grid_flag)
      computeKKBase->post_process_grid_kokkos(tau_index,1,DAT::t_float_2d_lr(),NULL,DAT::t_float_1d_strided());

    if (tau_index == 0 || ctau->post_process_grid_flag)
      Kokkos::deep_copy(d_tau_vector, computeKKBase->d_vector_grid);
    else {
      d_array = computeKKBase->d_array_grid;
      copymode = 1;
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeDtGrid_LoadTauVecFromArray>(0,nglocal),*this);
      copymode = 0;
    }
  } else if (tau_which == FIX) {
    if (!ftau->kokkos_flag)
      error->all(FLERR,"Cannot (yet) use non-Kokkos fixes with compute dt/grid/kk");
    KokkosBase* computeKKBase = dynamic_cast<KokkosBase*>(ftau);
    if (tau_index == 0)
      d_tau_vector = computeKKBase->d_vector_grid;
    else {
      d_array = computeKKBase->d_array_grid;
      copymode = 1;
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeDtGrid_LoadTauVecFromArray>(0,nglocal),*this);
      copymode = 0;
    }
  }

  if (temp_which == COMPUTE) {
    if (!ctemp->kokkos_flag)
      error->all(FLERR,"Cannot (yet) use non-Kokkos computes with compute dt/grid/kk");
    KokkosBase* computeKKBase = dynamic_cast<KokkosBase*>(ctemp);
    if (!(ctemp->invoked_flag & INVOKED_PER_GRID)) {
      computeKKBase->compute_per_grid_kokkos();
      ctemp->invoked_flag |= INVOKED_PER_GRID;
    }
    if (ctemp->post_process_grid_flag)
      computeKKBase->post_process_grid_kokkos(temp_index,1,DAT::t_float_2d_lr(),NULL,DAT::t_float_1d_strided());

    if (temp_index == 0 || ctemp->post_process_grid_flag)
      Kokkos::deep_copy(d_temp_vector, computeKKBase->d_vector_grid);
    else {
      d_array = computeKKBase->d_array_grid;
      copymode = 1;
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeDtGrid_LoadTempVecFromArray>(0,nglocal),*this);
      copymode = 0;
    }
  } else if (temp_which == FIX) {
    if (!ftemp->kokkos_flag)
      error->all(FLERR,"Cannot (yet) use non-Kokkos fixes with compute dt/grid/kk");
    KokkosBase* computeKKBase = dynamic_cast<KokkosBase*>(ftemp);
    if (temp_index == 0)
      d_temp_vector = computeKKBase->d_vector_grid;
    else {
      d_array = computeKKBase->d_array_grid;
      copymode = 1;
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeDtGrid_LoadTempVecFromArray>(0,nglocal),*this);
      copymode = 0;
    }
  }

  if (usq_which == COMPUTE) {
    if (!cusq->kokkos_flag)
      error->all(FLERR,"Cannot (yet) use non-Kokkos computes with compute dt/grid/kk");
    KokkosBase* computeKKBase = dynamic_cast<KokkosBase*>(cusq);
    if (!(cusq->invoked_flag & INVOKED_PER_GRID)) {
      computeKKBase->compute_per_grid_kokkos();
      cusq->invoked_flag |= INVOKED_PER_GRID;
    }
    if (cusq->post_process_grid_flag)
      computeKKBase->post_process_grid_kokkos(usq_index,1,DAT::t_float_2d_lr(),NULL,DAT::t_float_1d_strided());

    if (usq_index == 0 || cusq->post_process_grid_flag)
      Kokkos::deep_copy(d_usq_vector, computeKKBase->d_vector_grid);
    else {
      d_array = computeKKBase->d_array_grid;
      copymode = 1;
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeDtGrid_LoadUsqVecFromArray>(0,nglocal),*this);
      copymode = 0;
    }
  } else if (usq_which == FIX) {
    if (!fusq->kokkos_flag)
      error->all(FLERR,"Cannot (yet) use non-Kokkos fixes with compute dt/grid/kk");
    KokkosBase* computeKKBase = dynamic_cast<KokkosBase*>(fusq);
    if (usq_index == 0)
      d_usq_vector = computeKKBase->d_vector_grid;
    else {
      d_array = computeKKBase->d_array_grid;
      copymode = 1;
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeDtGrid_LoadUsqVecFromArray>(0,nglocal),*this);
      copymode = 0;
    }
  }

  if (vsq_which == COMPUTE) {
    if (!cvsq->kokkos_flag)
      error->all(FLERR,"Cannot (yet) use non-Kokkos computes with compute dt/grid/kk");
    KokkosBase* computeKKBase = dynamic_cast<KokkosBase*>(cvsq);
    if (!(cvsq->invoked_flag & INVOKED_PER_GRID)) {
      computeKKBase->compute_per_grid_kokkos();
      cvsq->invoked_flag |= INVOKED_PER_GRID;
    }
    if (cvsq->post_process_grid_flag)
      computeKKBase->post_process_grid_kokkos(vsq_index,1,DAT::t_float_2d_lr(),NULL,DAT::t_float_1d_strided());

    if (vsq_index == 0 || cvsq->post_process_grid_flag)
      Kokkos::deep_copy(d_vsq_vector, computeKKBase->d_vector_grid);
    else {
      d_array = computeKKBase->d_array_grid;
      copymode = 1;
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeDtGrid_LoadVsqVecFromArray>(0,nglocal),*this);
      copymode = 0;
    }
  } else if (vsq_which == FIX) {
    if (!fvsq->kokkos_flag)
      error->all(FLERR,"Cannot (yet) use non-Kokkos fixes with compute dt/grid/kk");
    KokkosBase* computeKKBase = dynamic_cast<KokkosBase*>(fvsq);
    if (vsq_index == 0)
      d_vsq_vector = computeKKBase->d_vector_grid;
    else {
      d_array = computeKKBase->d_array_grid;
      copymode = 1;
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeDtGrid_LoadVsqVecFromArray>(0,nglocal),*this);
      copymode = 0;
    }
  }

  if (wsq_which == COMPUTE) {
    if (!cwsq->kokkos_flag)
      error->all(FLERR,"Cannot (yet) use non-Kokkos computes with compute dt/grid/kk");
    KokkosBase* computeKKBase = dynamic_cast<KokkosBase*>(cwsq);
    if (!(cwsq->invoked_flag & INVOKED_PER_GRID)) {
      computeKKBase->compute_per_grid_kokkos();
      cwsq->invoked_flag |= INVOKED_PER_GRID;
    }
    if (cwsq->post_process_grid_flag)
      computeKKBase->post_process_grid_kokkos(wsq_index,1,DAT::t_float_2d_lr(),NULL,DAT::t_float_1d_strided());

    if (wsq_index == 0 || cwsq->post_process_grid_flag)
      Kokkos::deep_copy(d_wsq_vector, computeKKBase->d_vector_grid);
    else {
      d_array = computeKKBase->d_array_grid;
      copymode = 1;
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeDtGrid_LoadWsqVecFromArray>(0,nglocal),*this);
      copymode = 0;
    }
  } else if (wsq_which == FIX) {
    if (!fwsq->kokkos_flag)
      error->all(FLERR,"Cannot (yet) use non-Kokkos fixes with compute dt/grid/kk");
    KokkosBase* computeKKBase = dynamic_cast<KokkosBase*>(fwsq);
    if (wsq_index == 0)
      d_wsq_vector = computeKKBase->d_vector_grid;
    else {
      d_array = computeKKBase->d_array_grid;
      copymode = 1;
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeDtGrid_LoadWsqVecFromArray>(0,nglocal),*this);
      copymode = 0;
    }
  }

  // calculate per grid cell timestep for cells in group
  GridKokkos* grid_kk = ((GridKokkos*)grid);
  grid_kk->sync(Device,CELL_MASK);
  d_cinfo = grid_kk->k_cinfo.d_view;
  d_cells = grid_kk->k_cells.d_view;
  d_cellcount = grid_kk->d_cellcount;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeDtGrid_ComputePerGrid>(0,nglocal),*this);
  copymode = 0;

  k_vector_grid.modify_device();
  k_vector_grid.sync_host();
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeDtGridKokkos::operator()(TagComputeDtGrid_LoadTauVecFromArray, const int &i) const {
  d_tau_vector(i) = d_array(i,tau_index-1);
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeDtGridKokkos::operator()(TagComputeDtGrid_LoadTempVecFromArray, const int &i) const {
  d_temp_vector(i) = d_array(i,temp_index-1);
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeDtGridKokkos::operator()(TagComputeDtGrid_LoadUsqVecFromArray, const int &i) const {
  d_usq_vector(i) = d_array(i,usq_index-1);
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeDtGridKokkos::operator()(TagComputeDtGrid_LoadVsqVecFromArray, const int &i) const {
  d_vsq_vector(i) = d_array(i,vsq_index-1);
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeDtGridKokkos::operator()(TagComputeDtGrid_LoadWsqVecFromArray, const int &i) const {
  d_wsq_vector(i) = d_array(i,wsq_index-1);
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeDtGridKokkos::operator()(TagComputeDtGrid_ComputePerGrid, const int &i) const {

  d_vector_grid(i) = 0.;

  // exclude cells not in the specified group
  if (!(d_cinfo[i].mask & groupbit)) return;

  // exclude cells that have no particles
  //  (includes split cells and unsplit cells interior to surface objects)
  if (d_cellcount[i] == 0) return;

  // check sufficiency of cell data to calculate cell dt
  double velmag2 = d_usq_vector(i) + d_vsq_vector(i) + d_wsq_vector(i);
  if ( !( (d_tau_vector(i) > EPSILON) && (velmag2 > EPSILON)) ) return;

  double cell_dt_desired = collision_fraction*d_tau_vector(i);

  // cell size
  double dx = d_cells[i].hi[0] - d_cells[i].lo[0];
  double dy = d_cells[i].hi[1] - d_cells[i].lo[1];
  double dz = d_cells[i].hi[2] - d_cells[i].lo[2];

  // cell dt based on transit time using average velocities
  double dt_candidate;
  double umag = sqrt(d_usq_vector(i));
  if (umag > 0.) {
    dt_candidate = transit_fraction*dx/umag;
    cell_dt_desired = MIN(dt_candidate,cell_dt_desired);
  }
  double vmag = sqrt(d_vsq_vector(i));
  if (vmag > 0.) {
    dt_candidate = transit_fraction*dy/vmag;
    cell_dt_desired = MIN(dt_candidate,cell_dt_desired);
  }
  if (dimension == 3) {
    double wmag = sqrt(d_wsq_vector(i));
    if (wmag > 0.) {
      dt_candidate = transit_fraction*dz/wmag;
      cell_dt_desired = MIN(dt_candidate,cell_dt_desired);
    }
  }

  // cell dt based on transit time using maximum most probable speed
  double vrm_max = sqrt(2.0*boltz * d_temp_vector(i) / min_species_mass);
  dt_candidate = transit_fraction*dx/vrm_max;
  cell_dt_desired = MIN(dt_candidate,cell_dt_desired);
  dt_candidate = transit_fraction*dy/vrm_max;
  cell_dt_desired = MIN(dt_candidate,cell_dt_desired);
  if (dimension == 3) {
    dt_candidate = transit_fraction*dz/vrm_max;
    cell_dt_desired = MIN(dt_candidate,cell_dt_desired);
  }

  d_vector_grid(i) = cell_dt_desired;
}

/* ----------------------------------------------------------------------
   reallocate arrays
------------------------------------------------------------------------- */
void ComputeDtGridKokkos::reallocate()
{
  if (grid->nlocal == nglocal) return;

  nglocal = grid->nlocal;

  memoryKK->destroy_kokkos(k_vector_grid,vector_grid);
  memoryKK->create_kokkos(k_vector_grid,vector_grid,nglocal,"ComputeDtGridKokkos:vector_grid");
  d_vector_grid = k_vector_grid.d_view;

  d_tau_vector = DAT::t_float_1d ("d_tau_vector", nglocal);
  d_temp_vector = DAT::t_float_1d ("d_temp_vector", nglocal);
  d_usq_vector = DAT::t_float_1d ("d_usq_vector", nglocal);
  d_vsq_vector = DAT::t_float_1d ("d_vsq_vector", nglocal);
  d_wsq_vector = DAT::t_float_1d ("d_wsq_vector", nglocal);
}

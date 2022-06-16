/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"
#include "fix_dt_kokkos.h"

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
#include "particle.h"
#include "mixture.h"
#include "Kokkos_Atomic.hpp"

enum{NONE,COMPUTE,FIX};

#define INVOKED_PER_GRID 16

using namespace SPARTA_NS;

/* ----------------------------------------------------------------------
   set global timestep (as a function of cell desired timesteps)
------------------------------------------------------------------------- */

FixDtKokkos::FixDtKokkos(SPARTA *sparta, int narg, char **arg) :
  FixDt(sparta, narg, arg)
{
  kokkos_flag = 1;
  execution_space = Device;

  nglocal = 0;
  reallocate();

  grid_kk = (GridKokkos*) grid;
  dimension = domain->dimension;
  boltz = update->boltz;
  dt_global = grid_kk->dt_global;

  k_dtsum = DAT::tdual_float_scalar("FixDtKokkos:dtsum");
  d_dtsum = k_dtsum.d_view;
  h_dtsum = k_dtsum.h_view;

  k_dtmin = DAT::tdual_float_scalar("FixDtKokkos:dtmin");
  d_dtmin = k_dtmin.d_view;
  h_dtmin = k_dtmin.h_view;

  k_dtmax = DAT::tdual_float_scalar("FixDtKokkos:dtmax");
  d_dtmax = k_dtmax.d_view;
  h_dtmax = k_dtmax.h_view;
}

/* ---------------------------------------------------------------------- */

FixDtKokkos::~FixDtKokkos()
{
  if (copymode) return;
}

/* ---------------------------------------------------------------------- */
void FixDtKokkos::end_of_step()
{
  reallocate();

  if (lambdawhich == FIX && update->ntimestep % flambda->per_grid_freq)
    error->all(FLERR,"Compute lambda/grid fix not computed at compatible time");
  if (tempwhich == FIX && update->ntimestep % ftemp->per_grid_freq)
    error->all(FLERR,"Compute lambda/grid fix not computed at compatible time");

  // grab lambda values from compute or fix
  // invoke lambda compute as needed
  if (lambdawhich == COMPUTE) {
    if (!clambda->kokkos_flag)
      error->all(FLERR,"Cannot (yet) use non-Kokkos computes with fix dt/kk");
    KokkosBase* computeKKBase = dynamic_cast<KokkosBase*>(clambda);
    if (!(clambda->invoked_flag & INVOKED_PER_GRID)) {
      computeKKBase->compute_per_grid_kokkos();
      clambda->invoked_flag |= INVOKED_PER_GRID;
    }

    if (clambda->post_process_grid_flag)
      computeKKBase->post_process_grid_kokkos(lambdaindex,1,DAT::t_float_2d_lr(),NULL,DAT::t_float_1d_strided());

    if (lambdaindex == 0 || clambda->post_process_grid_flag)
      Kokkos::deep_copy(d_lambda_vector, computeKKBase->d_vector);
    else {
      d_array = computeKKBase->d_array_grid;
      copymode = 1;
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixDt_LoadLambdaVecFromArray>(0,nglocal),*this);
      DeviceType().fence();
      copymode = 0;
    }
  } else if (lambdawhich == FIX) {
    if (!flambda->kokkos_flag)
      error->all(FLERR,"Cannot (yet) use non-Kokkos fixes with fix dt/kk");
    KokkosBase* computeKKBase = dynamic_cast<KokkosBase*>(flambda);
    if (lambdaindex == 0)
      d_lambda_vector = computeKKBase->d_vector;
    else {
      d_array = computeKKBase->d_array_grid;
      copymode = 1;
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixDt_LoadLambdaVecFromArray>(0,nglocal),*this);
      DeviceType().fence();
      copymode = 0;
    }
  }

  // grab temperature values from fix
  if (!ftemp->kokkos_flag)
    error->all(FLERR,"Cannot (yet) use non-Kokkos fixes with fix dt/kk");
  KokkosBase* computeKKBase = dynamic_cast<KokkosBase*>(ftemp);
  if (tempindex == 0)
    d_temp_vector = computeKKBase->d_vector;
  else {
    d_array = computeKKBase->d_array_grid;
    copymode = 1;
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixDt_LoadTempVecFromArray>(0,nglocal),*this);
    DeviceType().fence();
    copymode = 0;
  }

  // grab usq values from fix
  if (!fusq->kokkos_flag)
    error->all(FLERR,"Cannot (yet) use non-Kokkos fixes with fix dt/kk");
  computeKKBase = dynamic_cast<KokkosBase*>(fusq);
  if (usqindex == 0)
    d_usq_vector = computeKKBase->d_vector;
  else {
    d_array = computeKKBase->d_array_grid;
    copymode = 1;
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixDt_LoadUsqVecFromArray>(0,nglocal),*this);
    DeviceType().fence();
    copymode = 0;
  }

  // grab vsq values from fix
  if (!fvsq->kokkos_flag)
    error->all(FLERR,"Cannot (yet) use non-Kokkos fixes with fix dt/kk");
  computeKKBase = dynamic_cast<KokkosBase*>(fvsq);
  if (vsqindex == 0)
    d_vsq_vector = computeKKBase->d_vector;
  else {
    d_array = computeKKBase->d_array_grid;
    copymode = 1;
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixDt_LoadVsqVecFromArray>(0,nglocal),*this);
    DeviceType().fence();
    copymode = 0;
  }

  // grab wsq values from fix
  if (!fwsq->kokkos_flag)
    error->all(FLERR,"Cannot (yet) use non-Kokkos fixes with fix dt/kk");
  computeKKBase = dynamic_cast<KokkosBase*>(fwsq);
  if (wsqindex == 0)
    d_wsq_vector = computeKKBase->d_vector;
  else {
    d_array = computeKKBase->d_array_grid;
    copymode = 1;
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixDt_LoadWsqVecFromArray>(0,nglocal),*this);
    DeviceType().fence();
    copymode = 0;
  }

  // allow dt growth and reduce later if necessary
  grid_kk->dt_global *= 2;
  dt_global = grid_kk->dt_global;

  // compute cell desired timestep
  d_cells = grid_kk->k_cells.d_view;
  Kokkos::deep_copy(d_dtsum, 0.0);
  copymode = 1;
  if (!sparta->kokkos->need_atomics)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixDt_SetCellDtDesired<0> >(0,nglocal),*this);
  else if (sparta->kokkos->atomic_reduction)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixDt_SetCellDtDesired<1> >(0,nglocal),*this);
  else {
    error->all(FLERR,"Parallel_reduce not implemented for fix dt/kk");
    exit(3);
  }
  DeviceType().fence();
  copymode = 0;

  // limit cell desired timestep
  k_dtsum.sync_host();
  double dt_sum_global = 0.;
  double dt_sum = h_dtsum();
  MPI_Allreduce(&dt_sum,&dt_sum_global,1,MPI_DOUBLE,MPI_SUM,world);
  double avg_cell_dt = dt_sum_global/grid_kk->ncell;
  std::cout << "*********** avg_cell_dt = " << avg_cell_dt << std::endl;
  double dt_factor = 1000.0;
  h_dtmin() = avg_cell_dt/dt_factor;
  h_dtmax() = avg_cell_dt*dt_factor;
  k_dtmin.sync_device();
  k_dtmax.sync_device();
  grid_kk->dt_global = 2.*h_dtmin();  // reset global timestep
  dt_global = grid_kk->dt_global;
  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixDt_LimitCellDtDesired>(0,nglocal),*this);
  DeviceType().fence();
  copymode = 0;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixDtKokkos::operator()(TagFixDt_LoadLambdaVecFromArray, const int &i) const {
  d_lambda_vector(i) = d_array(i,lambdaindex-1);
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixDtKokkos::operator()(TagFixDt_LoadTempVecFromArray, const int &i) const {
  d_temp_vector(i) = d_array(i,tempindex-1);
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixDtKokkos::operator()(TagFixDt_LoadUsqVecFromArray, const int &i) const {
  d_usq_vector(i) = d_array(i,tempindex-1);
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixDtKokkos::operator()(TagFixDt_LoadVsqVecFromArray, const int &i) const {
  d_vsq_vector(i) = d_array(i,tempindex-1);
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixDtKokkos::operator()(TagFixDt_LoadWsqVecFromArray, const int &i) const {
  d_wsq_vector(i) = d_array(i,tempindex-1);
}

/* ---------------------------------------------------------------------- */

template<int ATOMIC_REDUCTION>
KOKKOS_INLINE_FUNCTION
void FixDtKokkos::operator()(TagFixDt_SetCellDtDesired<ATOMIC_REDUCTION>, const int &i) const {

  // cell dt based on mean collision time
  double transit_fraction = 0.25;
  double collision_fraction = 0.1;
  double mean_collision_time = d_lambda_vector(i)/sqrt(d_usq_vector(i) +
                                                       d_vsq_vector(i) +
                                                       d_wsq_vector(i));
  if (mean_collision_time > 0.)
    d_cells[i].dt_desired = collision_fraction*mean_collision_time;

  // cell size
  double dx = d_cells[i].hi[0] - d_cells[i].lo[0];
  double dy = d_cells[i].hi[1] - d_cells[i].lo[1];
  double dz = 0.;

  // cell dt based on transit time using average velocities
  double dt_candidate = transit_fraction*dx/sqrt(d_usq_vector(i));
  d_cells[i].dt_desired = MIN(dt_candidate,d_cells[i].dt_desired);

  dt_candidate = transit_fraction*dy/sqrt(d_vsq_vector(i));
  d_cells[i].dt_desired = MIN(dt_candidate,d_cells[i].dt_desired);

  if (dimension == 3) {
    dz = d_cells[i].hi[2] - d_cells[i].lo[2];
    dt_candidate = transit_fraction*dz/sqrt(d_wsq_vector(i));
    d_cells[i].dt_desired = MIN(dt_candidate,d_cells[i].dt_desired);
  }

  // cell dt based on transit time using maximum most probable speed
  double vrm_max = sqrt(2.0*boltz * d_temp_vector(i) / min_species_mass);

  dt_candidate = transit_fraction*dx/vrm_max;
  d_cells[i].dt_desired = MIN(dt_candidate,d_cells[i].dt_desired);

  dt_candidate = transit_fraction*dy/vrm_max;
  d_cells[i].dt_desired = MIN(dt_candidate,d_cells[i].dt_desired);

  if (dimension == 3) {
    dt_candidate = transit_fraction*dz/vrm_max;
    d_cells[i].dt_desired = MIN(dt_candidate,d_cells[i].dt_desired);
  }

  if (ATOMIC_REDUCTION == 1)
    Kokkos::atomic_add(&d_dtsum(), d_cells[i].dt_desired);
  else if (ATOMIC_REDUCTION == 0)
    d_dtsum() += d_cells[i].dt_desired;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixDtKokkos::operator()(TagFixDt_LimitCellDtDesired, const int &i) const {

  if (d_cells[i].dt_desired < d_dtmin())
    d_cells[i].dt_desired = d_dtmin();
  if (d_cells[i].dt_desired > d_dtmax())
    d_cells[i].dt_desired = d_dtmax();
  if (d_cells[i].dt_desired < 0.51*dt_global)
    d_cells[i].dt_desired = 0.51*dt_global;
}

/* ----------------------------------------------------------------------
   reallocate arrays
------------------------------------------------------------------------- */
void FixDtKokkos::reallocate()
{
  if (grid->nlocal == nglocal) return;

  nglocal = grid->nlocal;

  d_lambda_vector = DAT::t_float_1d ("d_lambda_vector", nglocal);
  d_temp_vector = DAT::t_float_1d ("d_temp_vector", nglocal);
  d_usq_vector = DAT::t_float_1d ("d_usq_vector", nglocal);
  d_vsq_vector = DAT::t_float_1d ("d_vsq_vector", nglocal);
  d_wsq_vector = DAT::t_float_1d ("d_wsq_vector", nglocal);
}


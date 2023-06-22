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

#include "fix_dt_reset_kokkos.h"
#include "update.h"
#include "domain.h"
#include "memory_kokkos.h"
#include "error.h"
#include "kokkos.h"
#include "fix.h"
#include "compute.h"
#include "sparta_masks.h"
#include "Kokkos_Atomic.hpp"

enum{COMPUTE,FIX};

#define INVOKED_PER_GRID 16
#define BIG 1.0e20

using namespace SPARTA_NS;

/* ----------------------------------------------------------------------
   set global timestep (as a function of cell desired timesteps)
------------------------------------------------------------------------- */

FixDtResetKokkos::FixDtResetKokkos(SPARTA *sparta, int narg, char **arg) :
  FixDtReset(sparta, narg, arg)
{
  kokkos_flag = 1;
  execution_space = Device;
}

/* ---------------------------------------------------------------------- */

FixDtResetKokkos::~FixDtResetKokkos()
{
  if (copymode) return;
}

/* ---------------------------------------------------------------------- */

void FixDtResetKokkos::end_of_step()
{
  int nglocal = grid->nlocal;
  if (nglocal > maxgrid) {
    maxgrid = grid->nlocal;
    d_gridstep = DAT::t_float_1d ("d_gridstep", maxgrid);
  }

  if (step_which == FIX && update->ntimestep % fstep->per_grid_freq)
    error->all(FLERR,"Fix dt/reset fix not computed at compatible time");

  // grab dt values from compute or fix
  // invoke  compute as needed
  if (step_which == COMPUTE) {
    if (!cstep->kokkos_flag)
      error->all(FLERR,"Cannot (yet) use non-Kokkos computes with fix dt/reset/kk");
    KokkosBase* computeKKBase = dynamic_cast<KokkosBase*>(cstep);
    if (!(cstep->invoked_flag & INVOKED_PER_GRID)) {
      computeKKBase->compute_per_grid_kokkos();
      cstep->invoked_flag |= INVOKED_PER_GRID;
    }

    if (cstep->post_process_grid_flag)
      computeKKBase->post_process_grid_kokkos(step_index,1,DAT::t_float_2d_lr(),NULL,DAT::t_float_1d_strided());

    if (step_index == 0 || cstep->post_process_grid_flag)
      Kokkos::deep_copy(d_gridstep, computeKKBase->d_vector);
    else {
      d_array = computeKKBase->d_array_grid;
      copymode = 1;
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixDtReset_LoadGridstepVecFromArray>(0,nglocal),*this);
      copymode = 0;
    }
  } else if (step_which == FIX) {
    if (!fstep->kokkos_flag)
      error->all(FLERR,"Cannot (yet) use non-Kokkos fixes with fix dt/reset/kk");
    KokkosBase* computeKKBase = dynamic_cast<KokkosBase*>(fstep);
    if (step_index == 0)
      d_gridstep = computeKKBase->d_vector;
    else {
      d_array = computeKKBase->d_array_grid;
      copymode = 1;
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixDtReset_LoadGridstepVecFromArray>(0,nglocal),*this);
      copymode = 0;
    }
  }

  // calculate dtmin,dtmax,dtave
  // skip cells whose timestep is zero (e.g., no particles)
  double dtmin_me = BIG;
  double dtmax_me = -BIG;
  double dtsum_me = 0.;
  int count_me = 0;

  copymode = 1;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,TagFixDtReset_MinMaxSumReduction>(0,nglocal),
                          *this,
                          Kokkos::Min<double>(dtmin_me),
                          Kokkos::Max<double>(dtmax_me),
                          Kokkos::Sum<double>(dtsum_me),
                          Kokkos::Sum<int>(count_me) );
  copymode = 0;

  bigint bcount_me = count_me;
  bigint bcount;
  MPI_Allreduce(&bcount_me,&bcount,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  if (bcount == 0) return; // done if no cell computed a timestep

  double min_max_data_me[2];
  double min_max_data[2];
  min_max_data_me[0] = dtmin_me;
  min_max_data_me[1] = 1./dtmax_me;
  MPI_Allreduce(&min_max_data_me,&min_max_data,2,MPI_DOUBLE,MPI_MIN,world);
  dtmin = min_max_data[0];
  dtmax = 1./min_max_data[1];
  MPI_Allreduce(&dtsum_me,&dtave,1,MPI_DOUBLE,MPI_SUM,world);
  dtave /= bcount;

  // calculate new global timestep
  dtnew = (1.0-weight) * dtmin + weight * dtave;

  // reset global timestep if requested
  // also reset the global time
  if (resetflag) {
    update->time += (update->ntimestep - update->time_last_update) * update->dt;
    update->time_last_update = update->ntimestep;
    update->dt = dtnew;
  }
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixDtResetKokkos::operator()(TagFixDtReset_MinMaxSumReduction,
                                  const int &i,
                                  double& dtmin,
                                  double& dtmax,
                                  double& dtsum,
                                  int& count) const {
  if (d_gridstep(i) == 0.0) return;

  if (d_gridstep(i) < dtmin)
    dtmin = d_gridstep(i);
  if (d_gridstep(i) > dtmax)
    dtmax = d_gridstep(i);
  dtsum += d_gridstep(i);
  count++;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixDtResetKokkos::operator()(TagFixDtReset_LoadGridstepVecFromArray, const int &i) const {
  d_gridstep(i) = d_array(i,step_index-1);
}

/* ---------------------------------------------------------------------- */



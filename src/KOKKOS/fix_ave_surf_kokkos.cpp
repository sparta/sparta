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

#include "spatype.h"
#include "string.h"
#include "fix_ave_surf_kokkos.h"
#include "surf.h"
#include "domain.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "compute_surf_kokkos.h"
#include "compute_react_surf_kokkos.h"
#include "memory_kokkos.h"
#include "error.h"
#include "sparta_masks.h"

using namespace SPARTA_NS;

enum{COMPUTE,FIX,VARIABLE,CUSTOM};       // must match fix_ave_surf.cpp
enum{ONE,RUNNING};                       // must match fix_ave_surf.cpp

#define INVOKED_PER_SURF 32              // must match fix_ave_surf.cpp

/* ---------------------------------------------------------------------- */

FixAveSurfKokkos::FixAveSurfKokkos(SPARTA *sparta, int narg, char **arg) :
  FixAveSurf(sparta, narg, arg)
{
  kokkos_flag = 1;
  execution_space = Device;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  // only the all-tally path (averaging surf-tally computes) is accelerated
  //   on device.  count_tally is 0 or nvalues, enforced by the base ctor.
  // the non-tally path (fix/variable/custom inputs) runs on the host base class

  kokkosable = (count_tally && count_tally == nvalues);

  nstally = 0;
  tally2surf_all = NULL;
  acc_local_vec = NULL;
  acc_local = NULL;
}

/* ---------------------------------------------------------------------- */

FixAveSurfKokkos::~FixAveSurfKokkos()
{
  if (copymode) return;

  memory->destroy(tally2surf_all);
  memory->destroy(acc_local_vec);
  memory->destroy(acc_local);
}

/* ---------------------------------------------------------------------- */

void FixAveSurfKokkos::init()
{
  FixAveSurf::init();
}

/* ----------------------------------------------------------------------
   allocate per-local-surf device accumulator and host collate buffers
   build tally2surf_all mapping each local surf row to its surf ID
------------------------------------------------------------------------- */

void FixAveSurfKokkos::reallocate()
{
  int n = surf->nlocal + surf->nghost;
  if (n == nstally && d_acc.extent(0)) return;
  nstally = n;

  d_acc = DAT::t_float_2d_lr("ave/surf:acc",nstally,nvalues);

  memory->destroy(tally2surf_all);
  memory->destroy(acc_local_vec);
  memory->destroy(acc_local);
  acc_local_vec = NULL;
  acc_local = NULL;
  memory->create(tally2surf_all,nstally,"ave/surf:tally2surf_all");
  if (nvalues == 1) memory->create(acc_local_vec,nstally,"ave/surf:acc_local_vec");
  else memory->create(acc_local,nstally,nvalues,"ave/surf:acc_local");

  // surf ID of each local surf row, used by the host collate at output
  // matches the per-local-surf row order of the Kokkos surf-tally computes

  if (domain->dimension == 2) {
    Surf::Line *lines = surf->lines;
    for (int i = 0; i < nstally; i++) tally2surf_all[i] = lines[i].id;
  } else {
    Surf::Tri *tris = surf->tris;
    for (int i = 0; i < nstally; i++) tally2surf_all[i] = tris[i].id;
  }
}

/* ----------------------------------------------------------------------
   only does something if nvalid = current timestep
------------------------------------------------------------------------- */

void FixAveSurfKokkos::setup()
{
  if (kokkosable) reallocate();
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveSurfKokkos::end_of_step()
{
  int i,m,n;

  // non-tally path runs entirely on the host base class

  if (!kokkosable) {
    FixAveSurf::end_of_step();
    return;
  }

  // skip if not step which requires doing something

  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;

  if (nstally != surf->nlocal + surf->nghost) reallocate();

  // first sample of an averaging interval:
  //   zero the per-interval device tally accumulator (== clearing the host hash)
  //   zero the owned-surf accumulators if ave = ONE

  if (irepeat == 0) {
    Kokkos::deep_copy(d_acc,0.0);
    if (ave == ONE) {
      if (nvalues == 1)
        for (i = 0; i < nown; i++) accvec[i] = 0.0;
      else
        for (i = 0; i < nown; i++)
          for (m = 0; m < nvalues; m++) accarray[i][m] = 0.0;
    }
  }

  // accumulate this sample's compute tallies into d_acc on device
  // each value m reads a column of its compute's per-local-surf device tally
  // compute/fix/variable may invoke computes, so wrap with clear/add

  modify->clearstep_compute();

  copymode = 1;
  for (m = 0; m < nvalues; m++) {
    n = value2index[m];
    Compute *compute = modify->compute[n];

    if (!compute->kokkos_flag)
      error->all(FLERR,"Cannot (yet) use non-Kokkos computes with fix ave/surf/kk");

    if (!(compute->invoked_flag & INVOKED_PER_SURF)) {
      compute->compute_per_surf();
      compute->invoked_flag |= INVOKED_PER_SURF;
    }

    // grab the compute's per-local-surf device tally array

    if (strcmp(compute->style,"surf") == 0)
      ((ComputeSurfKokkos*) compute)->query_tally_surf_kokkos(d_tally);
    else if (strcmp(compute->style,"react/surf") == 0)
      ((ComputeReactSurfKokkos*) compute)->query_tally_surf_kokkos(d_tally);
    else
      error->all(FLERR,"Fix ave/surf/kk requires Kokkos compute surf or compute react/surf");

    acc_m = m;
    acc_col = (argindex[m] == 0) ? 0 : argindex[m] - 1;
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixAveSurf_Add_tally>(0,nstally),*this);
  }
  copymode = 0;

  // done if irepeat < nrepeat, else reset irepeat and nvalid

  nsample++;
  irepeat++;
  if (irepeat < nrepeat) {
    nvalid += nevery;
    modify->addstep_compute(nvalid);
    return;
  }

  irepeat = 0;
  nvalid = ntimestep+per_surf_freq - (nrepeat-1)*nevery;
  modify->addstep_compute(nvalid);

  // copy the device tally accumulator to the host

  auto h_acc = Kokkos::create_mirror_view(d_acc);
  Kokkos::deep_copy(h_acc,d_acc);

  // merge per-local-surf tallies to owned surfs via surf->collate (host MPI)
  // then add the collated interval sum to the owned-surf accumulators

  if (nvalues == 1) {
    for (i = 0; i < nstally; i++) acc_local_vec[i] = h_acc(i,0);
    surf->collate_vector(nstally,tally2surf_all,acc_local_vec,1,bufvec);
    for (i = 0; i < nown; i++) accvec[i] += bufvec[i];
  } else {
    for (i = 0; i < nstally; i++)
      for (m = 0; m < nvalues; m++) acc_local[i][m] = h_acc(i,m);
    surf->collate_array(nstally,nvalues,tally2surf_all,acc_local,bufarray);
    for (i = 0; i < nown; i++)
      for (m = 0; m < nvalues; m++) accarray[i][m] += bufarray[i][m];
  }

  // normalize the accumulators for output, just by # of samples

  if (ave == ONE) {
    if (nvalues == 1)
      for (i = 0; i < nown; i++) vector_surf[i] /= nsample;
    else
      for (i = 0; i < nown; i++)
        for (m = 0; m < nvalues; m++) array_surf[i][m] /= nsample;
  } else {
    if (nvalues == 1)
      for (i = 0; i < nown; i++) vector_surf[i] = accvec[i]/nsample;
    else
      for (i = 0; i < nown; i++)
        for (m = 0; m < nvalues; m++) array_surf[i][m] = accarray[i][m]/nsample;
  }

  // set values for surfs not in group to zero

  if (groupbit != 1) {
    if (nvalues == 1) {
      for (i = 0; i < nown; i++)
        if (!(masks[i] & groupbit)) vector_surf[i] = 0.0;
    } else {
      for (i = 0; i < nown; i++)
        if (!(masks[i] & groupbit))
          for (m = 0; m < nvalues; m++) array_surf[i][m] = 0.0;
    }
  }

  // reset nsample if ave = ONE

  if (ave == ONE) nsample = 0;
}

/* ----------------------------------------------------------------------
   add one value's per-local-surf compute tally column into d_acc
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixAveSurfKokkos::operator()(TagFixAveSurf_Add_tally, const int &i) const
{
  d_acc(i,acc_m) += d_tally(i,acc_col);
}

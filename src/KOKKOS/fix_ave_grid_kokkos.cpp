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
#include "stdlib.h"
#include "string.h"
#include "fix_ave_grid_kokkos.h"
#include "grid_kokkos.h"
#include "particle.h"
#include "comm.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "input.h"
#include "variable.h"
#include "memory_kokkos.h"
#include "error.h"
#include "sparta_masks.h"
#include "kokkos_base.h"


using namespace SPARTA_NS;

enum{PERGRID,PERGRIDSURF};
enum{COMPUTE,FIX,VARIABLE,CUSTOM};
enum{ONE,RUNNING};                      // multiple files
enum{INT,DOUBLE};                       // several files

#define INVOKED_PER_GRID 16
#define DELTAGRID 1024            // must be bigger than split cells per cell
#define DELTASURF 1024;

/* ---------------------------------------------------------------------- */

FixAveGridKokkos::FixAveGridKokkos(SPARTA *sparta, int narg, char **arg) :
  FixAveGrid(sparta, narg, arg)
{
  kokkos_flag = 1;
  execution_space = Device;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  // PERGRIDSURF (grid/surf inputs, e.g. compute isurf/grid) runs on the host:
  // the per-surf tally is produced on-device by the Kokkos compute and brought
  // to the host by its tallyinfo(), then collated to per-grid by the host base
  // class.  Skip all Kokkos-specific allocation and leave the host base ctor's
  // allocations intact; the overridden methods below delegate to FixAveGrid.
  // Clear kokkos_flag so device consumers of per-grid fix output do not try to
  // read the (unallocated) device views; they read the host array_grid instead.

  if (flavor == PERGRIDSURF) {
    kokkos_flag = 0;
    execution_space = Host;
    return;
  }

  nglocal = maxgrid = grid->nlocal;

  // allocate per-grid cell data storage
  // zero vector/array grid in case used by dump or load balancer

  if (nvalues == 1) {
    memory->destroy(vector_grid);
    vector_grid = NULL;
    memoryKK->grow_kokkos(k_vector_grid,vector_grid,nglocal,"ave/grid:vector_grid");
    d_vector_grid = k_vector_grid.view_device();
  } else {
    memory->destroy(array_grid);
    array_grid = NULL;
    memoryKK->grow_kokkos(k_array_grid,array_grid,nglocal,nvalues,"ave/grid:array_grid");
    d_array_grid = k_array_grid.view_device();
  }

  // allocate tally array
  // zero in case used by ave = RUNNING or accessed for immediate output

  memory->destroy(tally);
  tally = NULL;
  memoryKK->create_kokkos(k_tally,tally,nglocal,ntotal,"ave/grid:tally");
  d_tally = k_tally.view_device();

  k_numap = DAT::tdual_float_1d("ave/grid:numap",nvalues);
  k_umap = DAT::tdual_float_2d("ave/grid:umap",nvalues,tmax);
  k_uomap = DAT::tdual_float_2d("ave/grid:uomap",nvalues,tmax);

  for (int i = 0; i < nvalues; i++) {
    k_numap.view_host()(i) = numap[i];
    for (int j = 0; j < tmax; j++) {
      k_umap.view_host()(i,j) = umap[i][j];
      k_uomap.view_host()(i,j) = uomap[i][j];
    }
  }
  k_numap.modify_host();
  k_numap.sync_device();
  d_numap = k_numap.view_device();

  k_umap.modify_host();
  k_umap.sync_device();
  d_umap = k_umap.view_device();

  k_uomap.modify_host();
  k_uomap.sync_device();
  d_uomap = k_uomap.view_device();
}

/* ---------------------------------------------------------------------- */

FixAveGridKokkos::~FixAveGridKokkos()
{
  if (copymode) return;

  if (flavor == PERGRIDSURF) return;

  if (nvalues == 1) memoryKK->destroy_kokkos(k_vector_grid,vector_grid);
  else memoryKK->destroy_kokkos(k_array_grid,array_grid);
  memoryKK->destroy_kokkos(k_tally,tally);
  vector_grid = NULL;
  array_grid = tally = NULL;
}

/* ---------------------------------------------------------------------- */

void FixAveGridKokkos::init()
{
  // PERGRIDSURF path runs entirely on the host

  if (flavor == PERGRIDSURF) {
    FixAveGrid::init();
    return;
  }

  // set indices and check validity of all computes,fixes,variables,custom attributes

  for (int m = 0; m < nvalues; m++) {
    if (which[m] == COMPUTE) {
      int icompute = modify->find_compute(ids[m]);
      if (icompute < 0)
    error->all(FLERR,"Compute ID for fix ave/grid does not exist");
      value2index[m] = icompute;

    } else if (which[m] == FIX) {
      int ifix = modify->find_fix(ids[m]);
      if (ifix < 0)
    error->all(FLERR,"Fix ID for fix ave/grid does not exist");
      value2index[m] = ifix;

    } else if (which[m] == VARIABLE) {
      int ivariable = input->variable->find(ids[m]);
      if (ivariable < 0)
    error->all(FLERR,"Variable name for fix ave/grid does not exist");
      value2index[m] = ivariable;

    } else value2index[m] = -1;
  }
}

/* ----------------------------------------------------------------------
   only does something if nvalid = current timestep
------------------------------------------------------------------------- */

void FixAveGridKokkos::setup()
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveGridKokkos::end_of_step()
{
  int j,n;
  //int *itmp;

  // PERGRIDSURF path runs entirely on the host: the Kokkos compute's device
  // surf tally is brought to the host by its tallyinfo(), then the host base
  // class collates per-surf tallies to per-grid output

  if (flavor == PERGRIDSURF) {
    FixAveGrid::end_of_step();
    return;
  }

  // skip if not step which requires doing something

  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;

  // zero tally if ave = ONE and first sample
  // could do this with memset()

  copymode = 1;

  // grid cell migration (load balance / grid adaptation) reorders the
  // per-cell tally and output on the host; pull those changes onto the
  // device before accumulating so the running tally is not corrupted

  pergrid_sync(Device);

  if (ave == ONE && irepeat == 0)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixAveGrid_Zero_tally>(0,nglocal),*this);

  // accumulate results of computes,fixes,variables,custom attributes
  // compute/fix/variable may invoke computes so wrap with clear/add
  // NOTE: add more logic for fixes and variables if enable them

  modify->clearstep_compute();

  tally_on_host = 0;

  for (m = 0; m < nvalues; m++) {
    n = value2index[m];
    j = argindex[m];

    if (which[m] == COMPUTE) {
      Compute *compute = modify->compute[n];
      if (!compute->kokkos_flag)
        error->all(FLERR,"Cannot (yet) use non-Kokkos computes with fix ave/grid/kk");
    }

    // invoke compute if not previously invoked

    if (which[m] == COMPUTE) {
      Compute *compute = modify->compute[n];
      KokkosBase* computeKKBase = dynamic_cast<KokkosBase*>(compute);
      if (!(compute->invoked_flag & INVOKED_PER_GRID)) {
        computeKKBase->compute_per_grid_kokkos();
        compute->invoked_flag |= INVOKED_PER_GRID;
      }

      // accumulate one or more compute values to umap columns of tally array
      // if compute does not post-process, access its vec/array grid directly
      // else access uomap columns in its ctally array

      tally_to_device();

      if (post_process[m]) {
        ntally = numap[m];
        computeKKBase->query_tally_grid_kokkos(d_ctally);
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixAveGrid_Add_ctally>(0,nglocal),*this);
      } else {
        k = umap[m][0];
        if (j == 0) {
          d_compute_vector = computeKKBase->d_vector_grid;
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixAveGrid_Add_compute_vector>(0,nglocal),*this);
        } else {
          jm1 = j - 1;
          d_compute_array = computeKKBase->d_array_grid;
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixAveGrid_Add_compute_array>(0,nglocal),*this);
        }
      }

    // access fix fields, guaranteed to be ready
    // two paths, mirroring the host loop in FixAveGrid::end_of_step():
    //   fast path: the fix publishes its per-grid output as a device view
    //     (KokkosBase::d_vector_grid / d_array_grid), so accumulate on device
    //     with the same kernels used for computes.  do NOT gate this on
    //     fix->kokkos_flag: fix field/grid/kk deliberately clears kokkos_flag
    //     (its variable evaluation is host-only) yet still publishes a valid
    //     device array, while fix ave/grid/kk clears it in its PERGRIDSURF
    //     flavor precisely because the device views are not allocated there.
    //     the presence of a non-empty device view is the accurate test
    //   fallback: no usable device view (a non-Kokkos fix such as fix ablate,
    //     or the PERGRIDSURF flavor of fix ave/grid), so add the host
    //     per-grid output into the host tally and push it to the device,
    //     the same round trip the VARIABLE branch below performs

    } else if (which[m] == FIX) {
      Fix *ifix = modify->fix[n];
      KokkosBase *fixKKBase = dynamic_cast<KokkosBase*>(ifix);
      if (fixKKBase) fixKKBase->sync_pergrid_device_kokkos();
      k = umap[m][0];

      if (j == 0) {
        int device_ok = fixKKBase && fixKKBase->d_vector_grid.data() &&
          (int) fixKKBase->d_vector_grid.extent(0) >= nglocal;

        if (device_ok) {
          tally_to_device();
          d_fix_vector = fixKKBase->d_vector_grid;
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixAveGrid_Add_fix_vector>(0,nglocal),*this);
        } else {
          double *fix_vector = ifix->vector_grid;
          if (nglocal && !fix_vector)
            error->all(FLERR,"Fix used by fix ave/grid/kk does not produce "
                       "a per-grid vector");

          // the tally was last written on the device this step (zeroed and/or
          // accumulated by the kernels above), so mark it device-modified and
          // pull it to the host before adding to it, then push the sum back;
          // otherwise the host add would operate on stale values and the
          // sync_device() would clobber the on-device accumulation

          tally_to_host();
          for (int i = 0; i < nglocal; i++)
            tally[i][k] += fix_vector[i];
        }

      } else {
        jm1 = j - 1;
        int device_ok = fixKKBase && fixKKBase->d_array_grid.data() &&
          (int) fixKKBase->d_array_grid.extent(0) >= nglocal &&
          (int) fixKKBase->d_array_grid.extent(1) > jm1;

        if (device_ok) {
          tally_to_device();
          d_fix_array = fixKKBase->d_array_grid;
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixAveGrid_Add_fix_array>(0,nglocal),*this);
        } else {
          double **fix_array = ifix->array_grid;
          if (nglocal && !fix_array)
            error->all(FLERR,"Fix used by fix ave/grid/kk does not produce "
                       "a per-grid array");

          tally_to_host();
          for (int i = 0; i < nglocal; i++)
            tally[i][k] += fix_array[i][jm1];
        }
      }

    // evaluate grid-style variable, sum values to Kth column of tally array

    } else if (which[m] == VARIABLE) {
      k = umap[m][0];

      // compute_grid() with sumflag = 1 adds into the host tally, so the host
      //   copy has to be current first: earlier values in this same command
      //   accumulate on the device, and k_tally is not marked device-modified
      //   until after this loop.  Without the pull-down the sum would read
      //   stale values and the push-back would clobber the device work

      tally_to_host();
      input->variable->compute_grid(n,&tally[0][k],ntotal,1);

    // access custom attribute

    } else if (which[m] == CUSTOM) {

      auto gridKK = (GridKokkos*) grid;
      gridKK->sync(Device,CUSTOM_MASK);
      tally_to_device();

      k = umap[m][0];
      if (j == 0) {
        if (nvalues == 1) {
          if (grid->etype[n] == INT) {
            auto h_ewhich = gridKK->k_ewhich.view_host();
            auto h_eivec = gridKK->k_eivec.view_host();
            auto d_custom_vector = h_eivec[h_ewhich[n]].k_view.view_device();

            Kokkos::parallel_for(nglocal, SPARTA_CLASS_LAMBDA(int i) {
              d_tally(i,k) += d_custom_vector[i];
            });

          } else if (grid->etype[n] == DOUBLE) {
            auto h_ewhich = gridKK->k_ewhich.view_host();
            auto h_edvec = gridKK->k_edvec.view_host();
            auto d_custom_vector = h_edvec[h_ewhich[n]].k_view.view_device();

            Kokkos::parallel_for(nglocal, SPARTA_CLASS_LAMBDA(int i) {
              d_tally(i,k) += d_custom_vector[i];
            });

          }
        } else {
          if (grid->etype[n] == INT) {
            auto h_ewhich = gridKK->k_ewhich.view_host();
            auto h_eivec = gridKK->k_eivec.view_host();
            auto d_custom_vector = h_eivec[h_ewhich[n]].k_view.view_device();

            Kokkos::parallel_for(nglocal, SPARTA_CLASS_LAMBDA(int i) {
              d_tally(i,k) += d_custom_vector[i];
            });

          } else if (grid->etype[n] == DOUBLE) {
            auto h_ewhich = gridKK->k_ewhich.view_host();
            auto h_edvec = gridKK->k_edvec.view_host();
            auto d_custom_vector = h_edvec[h_ewhich[n]].k_view.view_device();

            Kokkos::parallel_for(nglocal, SPARTA_CLASS_LAMBDA(int i) {
              d_tally(i,k) += d_custom_vector[i];
            });

          }
        }
      } else {
        int jm1 = j - 1;
        if (nvalues == 1) {
          if (grid->etype[n] == INT) {
            auto h_ewhich = gridKK->k_ewhich.view_host();
            auto h_eiarray = gridKK->k_eiarray.view_host();
            auto d_custom_array = h_eiarray[h_ewhich[n]].k_view.view_device();

            Kokkos::parallel_for(nglocal, SPARTA_CLASS_LAMBDA(int i) {
              d_tally(i,k) += d_custom_array(i,jm1);
            });

          } else if (grid->etype[n] == DOUBLE) {
            auto h_ewhich = gridKK->k_ewhich.view_host();
            auto h_edarray = gridKK->k_edarray.view_host();
            auto d_custom_array = h_edarray[h_ewhich[n]].k_view.view_device();

            Kokkos::parallel_for(nglocal, SPARTA_CLASS_LAMBDA(int i) {
              d_tally(i,k) += d_custom_array(i,jm1);
            });

          }
        } else {
          if (grid->etype[n] == INT) {
            auto h_ewhich = gridKK->k_ewhich.view_host();
            auto h_edarray = gridKK->k_edarray.view_host();
            auto d_custom_array = h_edarray[h_ewhich[n]].k_view.view_device();

            Kokkos::parallel_for(nglocal, SPARTA_CLASS_LAMBDA(int i) {
              d_tally(i,k) += d_custom_array(i,jm1);
            });

          } else if (grid->etype[n] == DOUBLE) {
            auto h_ewhich = gridKK->k_ewhich.view_host();
            auto h_edarray = gridKK->k_edarray.view_host();
            auto d_custom_array = h_edarray[h_ewhich[n]].k_view.view_device();

            Kokkos::parallel_for(nglocal, SPARTA_CLASS_LAMBDA(int i) {
              d_tally(i,k) += d_custom_array(i,jm1);
            });

          }
        }
      }
    }
  }

  // a trailing run of host values may have left the tally on the host, so
  // push it down before claiming the device below

  tally_to_device();

  // the tally array was accumulated on the device this step
  // mark it modified so the host copy is refreshed if grid cells later
  // migrate (the migration hooks pack/unpack the host tally)

  k_tally.modify_device();

  // done if irepeat < nrepeat
  // else reset irepeat and nvalid

  nsample++;
  irepeat++;
  if (irepeat < nrepeat) {
    nvalid += nevery;
    modify->addstep_compute(nvalid);
    copymode = 0;
    return;
  }

  irepeat = 0;
  nvalid = ntimestep+per_grid_freq - (nrepeat-1)*nevery;
  modify->addstep_compute(nvalid);

  // create normalized output in vector_grid or array_grid
  // if post_process flag set, compute performs normalization via pp_grid()
  // else just divide by nsample

  if (nvalues == 1) {
    if (post_process[0]) {
      n = value2index[0];
      j = argindex[0];
      Compute *c = modify->compute[n];
      KokkosBase* cKKBase = dynamic_cast<KokkosBase*>(c);
      cKKBase->post_process_grid_kokkos(j,nsample,d_tally,map[0],d_vector_grid);
    } else {
      k = map[0][0];
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixAveGrid_Norm_vector_grid>(0,nglocal),*this);
    }

  } else {
    for (m = 0; m < nvalues; m++) {
      if (post_process[m]) {
        n = value2index[m];
        j = argindex[m];
        Compute *c = modify->compute[n];
        KokkosBase* cKKBase = dynamic_cast<KokkosBase*>(c);
        if (d_array_grid.data()) cKKBase->post_process_grid_kokkos(j,nsample,d_tally,map[m],
                             Kokkos::subview(d_array_grid,Kokkos::ALL(),m)); // need to use subview
      } else {
        k = map[m][0];
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixAveGrid_Norm_array_grid>(0,nglocal),*this);
      }
    }
  }

  if (nvalues == 1) {
    k_vector_grid.modify_device();
    k_vector_grid.sync_host();
  } else {
    k_array_grid.modify_device();
    k_array_grid.sync_host();
  }

  // set values for grid cells not in group to zero

  if (groupbit != 1) {
    GridKokkos* grid_kk = (GridKokkos*) grid;
    grid_kk->sync(Device,CINFO_MASK);
    d_cinfo = grid_kk->k_cinfo.view_device();
    if (nvalues == 1) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixAveGrid_Zero_group_vector>(0,nglocal),*this);
      k_vector_grid.modify_device();
      k_vector_grid.sync_host();
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixAveGrid_Zero_group_array>(0,nglocal),*this);
      k_array_grid.modify_device();
      k_array_grid.sync_host();
    }
  }

  // reset nsample if ave = ONE

  if (ave == ONE) nsample = 0;
  copymode = 0;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixAveGridKokkos::operator()(TagFixAveGrid_Zero_group_vector, const int &i) const {
  if (!(d_cinfo[i].mask & groupbit)) d_vector_grid(i) = 0.0;
}

KOKKOS_INLINE_FUNCTION
void FixAveGridKokkos::operator()(TagFixAveGrid_Zero_group_array, const int &i) const {
  if (!(d_cinfo[i].mask & groupbit))
    for (int m = 0; m < nvalues; m++) d_array_grid(i,m) = 0.0;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixAveGridKokkos::operator()(TagFixAveGrid_Zero_tally, const int &i) const {
  for (int j = 0; j < ntotal; j++)
    d_tally(i,j) = 0.0;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixAveGridKokkos::operator()(TagFixAveGrid_Add_ctally, const int &i) const {
  for (int itally = 0; itally < ntally; itally++) {
    const int k = d_umap(m,itally);
    const int kk = d_uomap(m,itally);
    d_tally(i,k) += d_ctally(i,kk);
  }
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixAveGridKokkos::operator()(TagFixAveGrid_Add_compute_vector, const int &i) const {
  d_tally(i,k) += d_compute_vector[i];
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixAveGridKokkos::operator()(TagFixAveGrid_Add_compute_array, const int &i) const {
  d_tally(i,k) += d_compute_array(i,jm1);
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixAveGridKokkos::operator()(TagFixAveGrid_Add_fix_vector, const int &i) const {
  d_tally(i,k) += d_fix_vector[i];
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixAveGridKokkos::operator()(TagFixAveGrid_Add_fix_array, const int &i) const {
  d_tally(i,k) += d_fix_array(i,jm1);
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixAveGridKokkos::operator()(TagFixAveGrid_Norm_vector_grid, const int &i) const {
  d_vector_grid[i] = d_tally(i,k) / nsample;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixAveGridKokkos::operator()(TagFixAveGrid_Norm_array_grid, const int &i) const {
  d_array_grid(i,m) = d_tally(i,k) / nsample;
}

/* ----------------------------------------------------------------------
   insure per-cell arrays are allocated long enough for N new cells
------------------------------------------------------------------------- */

void FixAveGridKokkos::grow_percell(int nnew)
{
  // PERGRIDSURF keeps its per-cell arrays in host memory (managed by the host
  // base class); reallocating them with Kokkos memory here would make the host
  // base destructor free a Kokkos-allocated pointer

  if (flavor == PERGRIDSURF) {
    FixAveGrid::grow_percell(nnew);
    return;
  }

  if (nglocal+nnew < maxgrid) return;
  maxgrid += DELTAGRID;
  int n = maxgrid;

  // resize with the device as the source of truth, then refresh the host
  // so both copies hold valid data and no dangling modify flag remains
  // (this is called from the grid migration hooks, which then edit the host)

  pergrid_sync(Device);

  if (nvalues == 1) {
    memoryKK->grow_kokkos(k_vector_grid,vector_grid,n,"ave/grid:vector_grid");
    d_vector_grid = k_vector_grid.view_device();
  } else {
    memoryKK->grow_kokkos(k_array_grid,array_grid,n,nvalues,"ave/grid:array_grid");
    d_array_grid = k_array_grid.view_device();
  }

  memoryKK->grow_kokkos(k_tally,tally,n,ntotal,"ave/grid:tally");
  d_tally = k_tally.view_device();

  pergrid_sync(Host);
}

/* ----------------------------------------------------------------------
   grid cell migration hooks (load balance / grid adaptation)
   the base class packs/unpacks/copies the per-cell tally and output using
   the host arrays, so bring the host up to date before it reads and mark
   the host modified after it writes; the device is refreshed lazily in
   end_of_step().  Without this the accumulated tally is only correct on the
   device and grid migration silently corrupts the averaged output when UVM
   is disabled.
------------------------------------------------------------------------- */

int FixAveGridKokkos::pack_grid_one(int icell, char *buf, int memflag)
{
  pergrid_sync(Host);
  return FixAveGrid::pack_grid_one(icell,buf,memflag);
}

/* ---------------------------------------------------------------------- */

int FixAveGridKokkos::unpack_grid_one(int icell, char *buf)
{
  pergrid_sync(Host);
  int n = FixAveGrid::unpack_grid_one(icell,buf);
  pergrid_modify(Host);
  return n;
}

/* ---------------------------------------------------------------------- */

void FixAveGridKokkos::copy_grid_one(int icell, int jcell)
{
  pergrid_sync(Host);
  FixAveGrid::copy_grid_one(icell,jcell);
  pergrid_modify(Host);
}

/* ---------------------------------------------------------------------- */

void FixAveGridKokkos::add_grid_one()
{
  pergrid_sync(Host);
  FixAveGrid::add_grid_one();
  pergrid_modify(Host);
}

/* ----------------------------------------------------------------------
   bring the per-grid output up to date on the device.

   The migration hooks above leave vector_grid/array_grid current on the host
   only, and end_of_step() is what normally refreshes the device.  Another
   Kokkos style can read this fix's output before then: fix adapt invokes
   compute lambda/grid from AdaptGrid::candidates_coarsen(), which is after
   refinement has already migrated cells, and that compute launches a kernel
   straight on d_array_grid.  Without this the kernel reads the values the
   cells held before the migration.

   PERGRIDSURF keeps its per-cell arrays in host memory and allocates no
   device views, so there is nothing to publish.
------------------------------------------------------------------------- */

void FixAveGridKokkos::sync_pergrid_device_kokkos()
{
  if (flavor == PERGRIDSURF) return;

  pergrid_sync(Device);

  // the consumer reads the cached handles, so re-point them at the views
  //   pergrid_sync() just refreshed

  if (nvalues == 1) d_vector_grid = k_vector_grid.view_device();
  else d_array_grid = k_array_grid.view_device();
}

/* ----------------------------------------------------------------------
   move k_tally to one side, only when it is not already there.  the value
   loop in end_of_step() mixes device values (computes, fixes that publish a
   device view, custom attributes) with host ones (a fix with no device view,
   a grid variable); transferring at the crossings collapses a run of host
   values into a single round trip
------------------------------------------------------------------------- */

void FixAveGridKokkos::tally_to_host()
{
  if (tally_on_host) return;
  k_tally.modify_device();
  k_tally.sync_host();
  tally_on_host = 1;
}

/* ---------------------------------------------------------------------- */

void FixAveGridKokkos::tally_to_device()
{
  if (!tally_on_host) return;
  k_tally.modify_host();
  k_tally.sync_device();
  tally_on_host = 0;
}

/* ----------------------------------------------------------------------
   sync/modify the per-grid dual views: the tally array plus whichever
   output array (vector_grid or array_grid) is in use
------------------------------------------------------------------------- */

void FixAveGridKokkos::pergrid_sync(ExecutionSpace space)
{
  if (space == Device) {
    if (nvalues == 1) k_vector_grid.sync_device();
    else k_array_grid.sync_device();
    k_tally.sync_device();
  } else {
    if (nvalues == 1) k_vector_grid.sync_host();
    else k_array_grid.sync_host();
    k_tally.sync_host();
  }
}

/* ---------------------------------------------------------------------- */

void FixAveGridKokkos::pergrid_modify(ExecutionSpace space)
{
  if (space == Device) {
    if (nvalues == 1) k_vector_grid.modify_device();
    else k_array_grid.modify_device();
    k_tally.modify_device();
  } else {
    if (nvalues == 1) k_vector_grid.modify_host();
    else k_array_grid.modify_host();
    k_tally.modify_host();
  }
}


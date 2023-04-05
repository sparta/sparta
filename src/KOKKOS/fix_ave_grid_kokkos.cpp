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
#include "input.h"
#include "variable.h"
#include "memory_kokkos.h"
#include "error.h"
#include "sparta_masks.h"
#include "kokkos_base.h"


using namespace SPARTA_NS;

enum{PERGRID,PERGRIDSURF};
enum{COMPUTE,FIX,VARIABLE};
enum{ONE,RUNNING};                // multiple files

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

  if (flavor == PERGRIDSURF)
    error->all(FLERR,"Cannot yet use Kokkos with fix ave/grid for grid/surf inputs");

  nglocal = maxgrid = grid->nlocal;

  // allocate per-grid cell data storage
  // zero vector/array grid in case used by dump or load balancer

  if (nvalues == 1) {
    memory->destroy(vector_grid);
    vector_grid = NULL;
    memoryKK->grow_kokkos(k_vector_grid,vector_grid,nglocal,"ave/grid:vector_grid");
    d_vector = k_vector_grid.d_view;
  } else {
    memory->destroy(array_grid);
    array_grid = NULL;
    memoryKK->grow_kokkos(k_array_grid,array_grid,nglocal,nvalues,"ave/grid:array_grid");
    d_array_grid = k_array_grid.d_view;
  }

  // allocate tally array
  // zero in case used by ave = RUNNING or accessed for immediate output

  memory->destroy(tally);
  tally = NULL;
  memoryKK->create_kokkos(k_tally,tally,nglocal,ntotal,"ave/grid:tally");
  d_tally = k_tally.d_view;

  k_numap = DAT::tdual_float_1d("ave/grid:numap",nvalues);
  k_umap = DAT::tdual_float_2d("ave/grid:umap",nvalues,tmax);
  k_uomap = DAT::tdual_float_2d("ave/grid:uomap",nvalues,tmax);

  for (int i = 0; i < nvalues; i++) {
    k_numap.h_view(i) = numap[i];
    for (int j = 0; j < tmax; j++) {
      k_umap.h_view(i,j) = umap[i][j];
      k_uomap.h_view(i,j) = uomap[i][j];
    }
  }
  k_numap.modify_host();
  k_numap.sync_device();
  d_numap = k_numap.d_view;

  k_umap.modify_host();
  k_umap.sync_device();
  d_umap = k_umap.d_view;

  k_uomap.modify_host();
  k_uomap.sync_device();
  d_uomap = k_uomap.d_view;
}

/* ---------------------------------------------------------------------- */

FixAveGridKokkos::~FixAveGridKokkos()
{
  if (copymode) return;

  if (nvalues == 1) memoryKK->destroy_kokkos(k_vector_grid,vector_grid);
  else memoryKK->destroy_kokkos(k_array_grid,array_grid);
  memoryKK->destroy_kokkos(k_tally,tally);
  vector_grid = NULL;
  array_grid = tally = NULL;
}

/* ---------------------------------------------------------------------- */

void FixAveGridKokkos::init()
{
  // set indices and check validity of all computes,fixes,variables

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

  // skip if not step which requires doing something

  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;

  // zero tally if ave = ONE and first sample
  // could do this with memset()

  copymode = 1;
  if (ave == ONE && irepeat == 0)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixAveGrid_Zero_tally>(0,nglocal),*this);

  // accumulate results of computes,fixes,variables
  // compute/fix/variable may invoke computes so wrap with clear/add
  // NOTE: add more logic for fixes and variables if enable them

  modify->clearstep_compute();

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

      if (post_process[m]) {
        ntally = numap[m];
        computeKKBase->query_tally_grid_kokkos(d_ctally);
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixAveGrid_Add_ctally>(0,nglocal),*this);
      } else {
        k = umap[m][0];
        if (j == 0) {
          d_compute_vector = computeKKBase->d_vector;
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixAveGrid_Add_compute_vector>(0,nglocal),*this);
        } else {
          jm1 = j - 1;
          d_compute_array = computeKKBase->d_array_grid;
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixAveGrid_Add_compute_array>(0,nglocal),*this);
        }
      }

    // access fix fields, guaranteed to be ready

    } else if (which[m] == FIX) {
      error->all(FLERR,"Cannot (yet) use fixes with fix ave/grid/kk");
      //k = umap[m][0];
      //if (j == 0) {
      //  double *d_fix_vector = modify->fix[n]->vector_grid; // need Kokkos version
      //  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixAveGrid_Add_fix_vector>(0,nglocal),*this);
      //} else {
      //  int jm1 = j - 1;
      //  double **fix_array = modify->fix[n]->array_grid; // need Kokkos version
      //  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixAveGrid_Add_fix_array>(0,nglocal),*this);
      //}

    // evaluate grid-style variable, sum values to Kth column of tally array

    } else if (which[m] == VARIABLE) {
      error->all(FLERR,"Cannot (yet) use variables with fix ave/grid/kk");
      //k = umap[m][0];
      //input->variable->compute_grid(n,&tally[0][k],ntotal,1); // need Kokkos version
    }

  }

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
      cKKBase->post_process_grid_kokkos(j,nsample,d_tally,map[0],d_vector);
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
    d_cinfo = grid_kk->k_cinfo.d_view;
    if (nvalues == 1)
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixAveGrid_Zero_group_vector>(0,nglocal),*this);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixAveGrid_Zero_group_array>(0,nglocal),*this);
  }

  // reset nsample if ave = ONE

  if (ave == ONE) nsample = 0;
  copymode = 0;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixAveGridKokkos::operator()(TagFixAveGrid_Zero_group_vector, const int &i) const {
  if (!(d_cinfo[i].mask & groupbit)) d_vector(i) = 0.0;
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
  d_vector[i] = d_tally(i,k) / nsample;
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
  if (nglocal+nnew < maxgrid) return;
  maxgrid += DELTAGRID;
  int n = maxgrid;

  if (nvalues == 1) {
    memoryKK->grow_kokkos(k_vector_grid,vector_grid,n,"ave/grid:vector_grid");
    d_vector = k_vector_grid.d_view;
    k_vector_grid.sync_host();
  } else {
    memoryKK->grow_kokkos(k_array_grid,array_grid,n,nvalues,"ave/grid:array_grid");
    d_array_grid = k_array_grid.d_view;
    k_array_grid.sync_host();
  }

  memoryKK->grow_kokkos(k_tally,tally,n,ntotal,"ave/grid:tally");
  d_tally = k_tally.d_view;
}


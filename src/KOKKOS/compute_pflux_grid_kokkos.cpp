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

/* ----------------------------------------------------------------------
   Contributing author: Alan Stagg (Sandia)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "compute_pflux_grid_kokkos.h"
#include "particle_kokkos.h"
#include "mixture.h"
#include "grid_kokkos.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "memory_kokkos.h"
#include "error.h"
#include "sparta_masks.h"
#include "kokkos.h"

using namespace SPARTA_NS;

// user keywords

enum{MOMXX,MOMYY,MOMZZ,MOMXY,MOMYZ,MOMXZ};

// internal accumulators

enum{MASSSUM,MVX,MVY,MVZ,MVXSQ,MVYSQ,MVZSQ,MVXVY,MVYVZ,MVXVZ,LASTSIZE};


/* ---------------------------------------------------------------------- */

ComputePFluxGridKokkos::ComputePFluxGridKokkos(SPARTA *sparta, int narg, char **arg) :
  ComputePFluxGrid(sparta, narg, arg)
{
  kokkos_flag = 1;

  k_unique = DAT::tdual_int_1d("compute/pflux/grid:unique",npergroup);
  for (int m = 0; m < npergroup; m++)
    k_unique.h_view(m) = unique[m];
  k_unique.modify_host();
  k_unique.sync_device();
  d_unique = k_unique.d_view;
}

/* ---------------------------------------------------------------------- */

ComputePFluxGridKokkos::~ComputePFluxGridKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_vector_grid,vector_grid);
  memoryKK->destroy_kokkos(k_tally,tally);
  vector_grid = NULL;
  tally = NULL;
}

/* ---------------------------------------------------------------------- */

void ComputePFluxGridKokkos::compute_per_grid()
{
  if (sparta->kokkos->prewrap) {
    ComputePFluxGrid::compute_per_grid();
  } else {
    compute_per_grid_kokkos();
    k_tally.modify_device();
    k_tally.sync_host();
  }
}

/* ---------------------------------------------------------------------- */

void ComputePFluxGridKokkos::compute_per_grid_kokkos()
{
  invoked_per_grid = update->ntimestep;

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,PARTICLE_MASK|SPECIES_MASK);
  d_particles = particle_kk->k_particles.d_view;
  d_species = particle_kk->k_species.d_view;

  GridKokkos* grid_kk = (GridKokkos*) grid;
  d_cellcount = grid_kk->d_cellcount;
  d_plist = grid_kk->d_plist;
  grid_kk->sync(Device,CINFO_MASK);
  d_cinfo = grid_kk->k_cinfo.d_view;

  d_s2g = particle_kk->k_species2group.d_view;
  int nlocal = particle->nlocal;

  // zero all accumulators

  Kokkos::deep_copy(d_tally,0.0);

  // loop over all particles, skip species not in mixture group

  need_dup = sparta->kokkos->need_dup<DeviceType>();
  if (particle_kk->sorted_kk && sparta->kokkos->need_atomics && !sparta->kokkos->atomic_reduction)
    need_dup = 0;

  if (need_dup)
    dup_tally = Kokkos::Experimental::create_scatter_view<typename Kokkos::Experimental::ScatterSum, typename Kokkos::Experimental::ScatterDuplicated>(d_tally);
  else
    ndup_tally = Kokkos::Experimental::create_scatter_view<typename Kokkos::Experimental::ScatterSum, typename Kokkos::Experimental::ScatterNonDuplicated>(d_tally);

  copymode = 1;
  if (particle_kk->sorted_kk && sparta->kokkos->need_atomics && !sparta->kokkos->atomic_reduction)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputePFluxGrid_compute_per_grid>(0,nglocal),*this);
  else {
    if (sparta->kokkos->need_atomics)
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputePFluxGrid_compute_per_grid_atomic<1> >(0,nlocal),*this);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputePFluxGrid_compute_per_grid_atomic<0> >(0,nlocal),*this);
  }
  DeviceType().fence();
  copymode = 0;

  if (need_dup) {
    Kokkos::Experimental::contribute(d_tally, dup_tally);
    dup_tally = decltype(dup_tally)(); // free duplicated memory
  }

  d_particles = t_particle_1d(); // destroy reference to reduce memory use
  d_plist = DAT::t_int_2d(); // destroy reference to reduce memory use
}

/* ---------------------------------------------------------------------- */

template<int NEED_ATOMICS>
KOKKOS_INLINE_FUNCTION
void ComputePFluxGridKokkos::operator()(TagComputePFluxGrid_compute_per_grid_atomic<NEED_ATOMICS>, const int &i) const {

  // The tally array is duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_tally = ScatterViewHelper<typename NeedDup<NEED_ATOMICS,DeviceType>::value,decltype(dup_tally),decltype(ndup_tally)>::get(dup_tally,ndup_tally);
  auto a_tally = v_tally.template access<typename AtomicDup<NEED_ATOMICS,DeviceType>::value>();

  const int ispecies = d_particles[i].ispecies;
  const int igroup = d_s2g(imix,ispecies);
  if (igroup < 0) return;

  const int icell = d_particles[i].icell;
  if (!(d_cinfo[icell].mask & groupbit)) return;

  const double mass = d_species[ispecies].mass;
  double *v = d_particles[i].v;

  // loop has all possible values particle needs to accumulate
  // subset defined by user values are indexed by accumulate vector

  int k = igroup*npergroup;

  for (int m = 0; m < npergroup; m++) {
    switch (d_unique[m]) {
    case MASSSUM:
      a_tally(icell,k++) += mass;
      break;
    case MVX:
      a_tally(icell,k++) += mass*v[0];
      break;
    case MVY:
      a_tally(icell,k++) += mass*v[1];
      break;
    case MVZ:
      a_tally(icell,k++) += mass*v[2];
      break;
    case MVXSQ:
      a_tally(icell,k++) += mass*v[0]*v[0];
      break;
    case MVYSQ:
      a_tally(icell,k++) += mass*v[1]*v[1];
      break;
    case MVZSQ:
      a_tally(icell,k++) += mass*v[2]*v[2];
      break;
    case MVXVY:
      a_tally(icell,k++) += mass*v[0]*v[1];
      break;
    case MVYVZ:
      a_tally(icell,k++) += mass*v[1]*v[2];
      break;
    case MVXVZ:
      a_tally(icell,k++) += mass*v[0]*v[2];
      break;
    }
  }
}

KOKKOS_INLINE_FUNCTION
void ComputePFluxGridKokkos::operator()(TagComputePFluxGrid_compute_per_grid, const int &icell) const {
  if (!(d_cinfo[icell].mask & groupbit)) return;

  const int np = d_cellcount[icell];

  for (int n = 0; n < np; n++) {

    const int i = d_plist(icell,n);

    const int ispecies = d_particles[i].ispecies;
    const int igroup = d_s2g(imix,ispecies);
    if (igroup < 0) return;

    const double mass = d_species[ispecies].mass;
    double *v = d_particles[i].v;

    int k = igroup*npergroup;

    for (int m = 0; m < npergroup; m++) {
      switch (d_unique[m]) {
      case MASSSUM:
        d_tally(icell,k++) += mass;
        break;
      case MVX:
        d_tally(icell,k++) += mass*v[0];
        break;
      case MVY:
        d_tally(icell,k++) += mass*v[1];
        break;
      case MVZ:
        d_tally(icell,k++) += mass*v[2];
        break;
      case MVXSQ:
        d_tally(icell,k++) += mass*v[0]*v[0];
        break;
      case MVYSQ:
        d_tally(icell,k++) += mass*v[1]*v[1];
        break;
      case MVZSQ:
        d_tally(icell,k++) += mass*v[2]*v[2];
        break;
      case MVXVY:
        d_tally(icell,k++) += mass*v[0]*v[1];
        break;
      case MVYVZ:
        d_tally(icell,k++) += mass*v[1]*v[2];
        break;
      case MVXVZ:
        d_tally(icell,k++) += mass*v[0]*v[2];
        break;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   query info about internal tally array for this compute
   index = which column of output (0 for vec, 1 to N for array)
   return # of tally quantities for this index
   also return array = ptr to tally array
   also return cols = ptr to list of columns in tally for this index
------------------------------------------------------------------------- */

int ComputePFluxGridKokkos::query_tally_grid_kokkos(DAT::t_float_2d_lr &d_array)
{
  d_array = d_tally;
  return 0;
}

/* ----------------------------------------------------------------------
   tally accumulated info to compute final normalized values
   index = which column of output (0 for vec, 1 to N for array)
   for etally = NULL:
     use internal tallied info for single timestep
     compute values for all grid cells
       store results in vector_grid with nstride = 1 (single col of array_grid)
   for etaylly = ptr to caller array:
     use external tallied info for many timesteps
     emap = list of etally columns to use, # of columns determined by index
     store results in caller's vec, spaced by nstride
   if norm = 0.0, set result to 0.0 directly so do not divide by 0.0
------------------------------------------------------------------------- */

void ComputePFluxGridKokkos::post_process_grid_kokkos(int index, int nsample,
                         DAT::t_float_2d_lr d_etally, int *emap, DAT::t_float_1d_strided d_vec)
{
  index--;
  int ivalue = index % nvalue;

  int lo = 0;
  int hi = nglocal;

  if (!d_etally.data()) {
    d_etally = d_tally;
    emap = map[index];
    d_vec = d_vector;
    nstride = 1;
  }

  this->d_etally = d_etally;
  this->d_vec = d_vec;
  this->nstride = nstride;

  // compute normalized final value for each grid cell
  // Vcm = Sum mv / Sum m = (Wx,Wy,Wz)
  // Wi = Sum mVi / M
  // momii = F/V Sum m(Vi - Wi)^2
  // momii = F/V [ Sum mVi^2 + M Wi^2 - 2 Wi Sum mVi ]
  // momii = F/V [ Sum mVi^2 - (Sum mVi)^2 / M ]
  // momij = F/V Sum m(Vi - Wi)(Vj - Wj)
  // momij = F/V [ Sum mViVj + M Wi Wj - Wj Sum mVi - Wi Sum mVj ]
  // momij = F/V [ Sum mViVj - Sum mVi Sum mVj / M ]

  GridKokkos* grid_kk = (GridKokkos*) grid;
  grid_kk->sync(Device,CINFO_MASK);
  d_cinfo = grid_kk->k_cinfo.d_view;
  fnum = update->fnum;

  mass = emap[0];
  switch (value[ivalue]) {

  case MOMXX:
  case MOMYY:
  case MOMZZ:
    {
      mv = emap[1];
      mvv = emap[2];
      copymode = 1;
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputePFluxGrid_post_process_grid_diag>(lo,hi),*this);
      DeviceType().fence();
      copymode = 0;
      break;
    }
  case MOMXY:
  case MOMYZ:
  case MOMXZ:
    {
      mv1 = emap[1];
      mv2 = emap[2];
      mvv = emap[3];
      copymode = 1;
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputePFluxGrid_post_process_grid_offdiag>(lo,hi),*this);
      DeviceType().fence();
      copymode = 0;
      break;
    }
  } // end switch
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputePFluxGridKokkos::operator()(TagComputePFluxGrid_post_process_grid_diag, const int &icell) const {
  double summass, summv, wt;
  summass = d_etally(icell,mass);
  if (summass == 0.0) d_vec[icell] = 0.0;
  else{
    wt = fnum * d_cinfo[icell].weight / d_cinfo[icell].volume;
    summv = d_etally(icell,mv);
    d_vec[icell] = wt * (d_etally(icell,mvv) - summv*summv/summass);
  }
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputePFluxGridKokkos::operator()(TagComputePFluxGrid_post_process_grid_offdiag, const int &icell) const {
  double summass, wt;
  summass = d_etally(icell,mass);
  if (summass == 0.0) d_vec[icell] = 0.0;
  else{
    wt = fnum * d_cinfo[icell].weight / d_cinfo[icell].volume;
    d_vec[icell] = wt * (d_etally(icell,mvv) -
                         d_etally(icell,mv1)*d_etally(icell,mv2)/summass);
  }
}

/* ----------------------------------------------------------------------
   reallocate arrays if nglocal has changed
   called by init() and whenever grid changes
------------------------------------------------------------------------- */

void ComputePFluxGridKokkos::reallocate()
{
  if (grid->nlocal == nglocal) return;

  memoryKK->destroy_kokkos(k_vector_grid,vector_grid);
  memoryKK->destroy_kokkos(k_tally,tally);
  nglocal = grid->nlocal;
  memoryKK->create_kokkos(k_vector_grid,vector_grid,nglocal,"pflux/grid:vector_grid");
  d_vector = k_vector_grid.d_view;
  memoryKK->create_kokkos(k_tally,tally,nglocal,ntotal,"pflux/grid:tally");
  d_tally = k_tally.d_view;
}

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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "compute_thermal_grid_kokkos.h"
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

enum{TEMP,PRESS};




/* ---------------------------------------------------------------------- */

ComputeThermalGridKokkos::ComputeThermalGridKokkos(SPARTA *sparta, int narg, char **arg) :
  ComputeThermalGrid(sparta, narg, arg)
{
  kokkos_flag = 1;
}

/* ---------------------------------------------------------------------- */

ComputeThermalGridKokkos::~ComputeThermalGridKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_vector_grid,vector_grid);
  memoryKK->destroy_kokkos(k_tally,tally);
  vector_grid = NULL;
  tally = NULL;
}

/* ---------------------------------------------------------------------- */

void ComputeThermalGridKokkos::compute_per_grid()
{
  if (sparta->kokkos->prewrap) {
    ComputeThermalGrid::compute_per_grid();
  } else {
    compute_per_grid_kokkos();
    k_tally.modify_device();
    k_tally.sync_host();
  }
}

/* ---------------------------------------------------------------------- */

void ComputeThermalGridKokkos::compute_per_grid_kokkos()
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
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeThermalGrid_compute_per_grid>(0,nglocal),*this);
  else {
    if (sparta->kokkos->need_atomics)
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeThermalGrid_compute_per_grid_atomic<1> >(0,nlocal),*this);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeThermalGrid_compute_per_grid_atomic<0> >(0,nlocal),*this);
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
void ComputeThermalGridKokkos::operator()(TagComputeThermalGrid_compute_per_grid_atomic<NEED_ATOMICS>, const int &i) const {

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

  // 6 tallies per particle: N, Mass, mVx, mVy, mVz, mV^2

  int k = igroup*npergroup;

  a_tally(icell,k++) += 1.0;
  a_tally(icell,k++) += mass;
  a_tally(icell,k++) += mass*v[0];
  a_tally(icell,k++) += mass*v[1];
  a_tally(icell,k++) += mass*v[2];
  a_tally(icell,k++) += mass * (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

KOKKOS_INLINE_FUNCTION
void ComputeThermalGridKokkos::operator()(TagComputeThermalGrid_compute_per_grid, const int &icell) const {
  if (!(d_cinfo[icell].mask & groupbit)) return;

  const int np = d_cellcount[icell];

  for (int n = 0; n < np; n++) {

    const int i = d_plist(icell,n);

    const int ispecies = d_particles[i].ispecies;
    const int igroup = d_s2g(imix,ispecies);
  if (igroup < 0) return;

    const int icell = d_particles[i].icell;

    const double mass = d_species[ispecies].mass;
    double *v = d_particles[i].v;

    // 6 tallies per particle: N, Mass, mVx, mVy, mVz, mV^2



    int k = igroup*npergroup;

    d_tally(icell,k++) += 1.0;
    d_tally(icell,k++) += mass;
    d_tally(icell,k++) += mass*v[0];
    d_tally(icell,k++) += mass*v[1];
    d_tally(icell,k++) += mass*v[2];
    d_tally(icell,k++) += mass * (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  }
}

/* ----------------------------------------------------------------------
   query info about internal tally array for this compute
   index = which column of output (0 for vec, 1 to N for array)
   return # of tally quantities for this index
   also return array = ptr to tally array
   also return cols = ptr to list of columns in tally for this index
------------------------------------------------------------------------- */

int ComputeThermalGridKokkos::query_tally_grid_kokkos(DAT::t_float_2d_lr &d_array)
{
  d_array = d_tally;
  return 0;
}

/* ----------------------------------------------------------------------
   tally accumulated info to compute final normalized values
   index = which column of output (0 for vec, 1 to N for array)
   for etally = NULL:
     use internal tallied info for single timestep, set nsample = 1
     compute values for all grid cells
       store results in vector_grid with nstride = 1 (single col of array_grid)
   for etaylly = ptr to caller array:
     use external tallied info for many timesteps
     nsample = additional normalization factor used by some values
     emap = list of etally columns to use, # of columns determined by index
     store results in caller's vec, spaced by nstride
   if norm = 0.0, set result to 0.0 directly so do not divide by 0.0
------------------------------------------------------------------------- */

void ComputeThermalGridKokkos::
post_process_grid_kokkos(int index, int nsample,
                  DAT::t_float_2d_lr d_etally, int *emap, DAT::t_float_1d_strided d_vec)
{
  index--;
  int ivalue = index % nvalue;

  int lo = 0;
  int hi = nglocal;
  //k = 0;

  if (!d_etally.data()) {
    nsample = 1;
    d_etally = d_tally;
    emap = map[index];
    d_vec = d_vector;
    nstride = 1;
  }

  this->d_etally = d_etally;
  this->d_vec = d_vec;
  this->nstride = nstride;
  this->nsample = nsample;

  // compute normalized final value for each grid cell
  // Vcm = Sum mv / Sum m
  // total KE = 0.5 * Sum m(v - Vcm)^2
  // KE = 0.5 * Sum(i=xyz) [(Sum mVi^2) - M Vcmx^2]
  // KE = 0.5 * Sum(i=xyz) [(Sum mVi^2) - (Sum mVx)^2 / M]
  // KE = 0.5 * Sum (mv^2) - [(Sum mVx)^2 + (Sum mVy)^2 + (Sum mVz)^2] / M
  // thermal temp = (2/(3NkB)) * KE
  // press = (2/(3V)) * KE

  tflag = 1;
  if (value[ivalue] == PRESS) tflag = 0;

  if (value[ivalue] == TEMP) prefactor = tprefactor;
  else if (value[ivalue] == PRESS) prefactor = pprefactor;

  GridKokkos* grid_kk = (GridKokkos*) grid;
  grid_kk->sync(Device,CINFO_MASK);
  d_cinfo = grid_kk->k_cinfo.d_view;

  n = emap[0];

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeThermalGrid_post_process_grid>(lo,hi),*this);
  DeviceType().fence();
  copymode = 0;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeThermalGridKokkos::operator()(TagComputeThermalGrid_post_process_grid, const int &icell) const {
  const double ncount = d_etally(icell,n);
  if (ncount == 0.0) d_vec[icell] = 0.0;
  else {
    const double mass = d_etally(icell,n+1);
    const double mvx = d_etally(icell,n+2);
    const double mvy = d_etally(icell,n+3);
    const double mvz = d_etally(icell,n+4);
    const double mvsq = d_etally(icell,n+5);
    d_vec[icell] = mvsq - (mvx*mvx + mvy*mvy + mvz*mvz)/mass;
    d_vec[icell] *= prefactor;
    if (tflag) d_vec[icell] /= ncount;
    else d_vec[icell] *= d_cinfo[icell].weight / d_cinfo[icell].volume / nsample;
  }
}

/* ----------------------------------------------------------------------
   reallocate arrays if nglocal has changed
   called by init() and whenever grid changes
------------------------------------------------------------------------- */

void ComputeThermalGridKokkos::reallocate()
{
  if (grid->nlocal == nglocal) return;

  memoryKK->destroy_kokkos(k_vector_grid,vector_grid);
  memoryKK->destroy_kokkos(k_tally,tally);
  nglocal = grid->nlocal;
  memoryKK->create_kokkos(k_vector_grid,vector_grid,nglocal,"thermal/grid:vector_grid");
  d_vector = k_vector_grid.d_view;
  memoryKK->create_kokkos(k_tally,tally,nglocal,ntotal,"thermal/grid:tally");
  d_tally = k_tally.d_view;
}

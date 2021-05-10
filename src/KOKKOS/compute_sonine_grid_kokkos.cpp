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
#include "compute_sonine_grid_kokkos.h"
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

enum{AMOM,BMOM};
enum{X,Y,Z};

/* ---------------------------------------------------------------------- */

ComputeSonineGridKokkos::ComputeSonineGridKokkos(SPARTA *sparta, int narg, char **arg) :
  ComputeSonineGrid(sparta, narg, arg)
{
  kokkos_flag = 1;

  k_which = DAT::tdual_int_1d("compute/sonine/grid:which",nvalue);
  k_moment = DAT::tdual_int_1d("compute/sonine/grid:moment",nvalue);
  k_order = DAT::tdual_int_1d("compute/sonine/grid:order",nvalue);
  for (int n = 0; n < nvalue; n++){
    k_which.h_view(n) = which[n];
    k_moment.h_view(n) = moment[n];
    k_order.h_view(n) = order[n];
  }
  k_which.modify_host();
  k_moment.modify_host();
  k_order.modify_host();

  k_which.sync_device();
  k_moment.sync_device();
  k_order.sync_device();

  d_which = k_which.d_view;
  d_moment = k_moment.d_view;
  d_order = k_order.d_view;
}

/* ---------------------------------------------------------------------- */

ComputeSonineGridKokkos::~ComputeSonineGridKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_vector_grid,vector_grid);
  memoryKK->destroy_kokkos(k_tally,tally);
  vector_grid = NULL;
  tally = NULL;
}

/* ---------------------------------------------------------------------- */

void ComputeSonineGridKokkos::compute_per_grid()
{
  if (sparta->kokkos->prewrap) {
    ComputeSonineGrid::compute_per_grid();
  } else {
    compute_per_grid_kokkos();
    k_tally.modify_device();
    k_tally.sync_host();
  }
}

/* ---------------------------------------------------------------------- */
void ComputeSonineGridKokkos::compute_per_grid_kokkos()
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
  Kokkos::deep_copy(d_vcom,0.0);

  need_dup = sparta->kokkos->need_dup<DeviceType>();
  if (particle_kk->sorted_kk && sparta->kokkos->need_atomics && !sparta->kokkos->atomic_reduction)
    need_dup = 0;

  if (need_dup) {
    dup_tally = Kokkos::Experimental::create_scatter_view<typename Kokkos::Experimental::ScatterSum, typename Kokkos::Experimental::ScatterDuplicated>(d_tally);
    dup_vcom_tally = Kokkos::Experimental::create_scatter_view<typename Kokkos::Experimental::ScatterSum, typename Kokkos::Experimental::ScatterDuplicated>(d_vcom);
  }
  else {
    ndup_tally = Kokkos::Experimental::create_scatter_view<typename Kokkos::Experimental::ScatterSum, typename Kokkos::Experimental::ScatterNonDuplicated>(d_tally);
    ndup_vcom_tally = Kokkos::Experimental::create_scatter_view<typename Kokkos::Experimental::ScatterSum, typename Kokkos::Experimental::ScatterNonDuplicated>(d_vcom);
  }

  // calculate center of mass velocity for each cell and group
  copymode = 1;
  if (particle_kk->sorted_kk && sparta->kokkos->need_atomics && !sparta->kokkos->atomic_reduction)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeSonineGrid_compute_vcom>(0,nglocal),*this);
  else {
    if (sparta->kokkos->need_atomics)
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeSonineGrid_compute_vcom_init_atomic<1> >(0,nlocal),*this);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeSonineGrid_compute_vcom_init_atomic<0> >(0,nlocal),*this);
    DeviceType().fence();
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeSonineGrid_normalize_vcom>(0,nglocal),*this);
  }
  DeviceType().fence();

  if (need_dup) {
    Kokkos::Experimental::contribute(d_vcom, dup_vcom_tally);
    dup_vcom_tally = decltype(dup_vcom_tally)(); // free duplicated memory
  }

  // tally sonine moments
  if (particle_kk->sorted_kk && sparta->kokkos->need_atomics && !sparta->kokkos->atomic_reduction)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeSonineGrid_compute_per_grid>(0,nglocal),*this);
  else {
    if (sparta->kokkos->need_atomics)
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeSonineGrid_compute_per_grid_atomic<1> >(0,nlocal),*this);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeSonineGrid_compute_per_grid_atomic<0> >(0,nlocal),*this);
  }
  DeviceType().fence();
  copymode = 0;

  if (need_dup) {
    Kokkos::Experimental::contribute(d_tally, dup_tally);
    dup_tally = decltype(dup_tally)(); // free duplicated memory
  }
}

/* ---------------------------------------------------------------------- */

template<int NEED_ATOMICS>
KOKKOS_INLINE_FUNCTION
void ComputeSonineGridKokkos::operator()(TagComputeSonineGrid_compute_vcom_init_atomic<NEED_ATOMICS>, const int &i) const {
  // compute COM velocity on this timestep for each cell and group

  // The tally array is duplicated for OpenMP, atomic for CUDA, and neither for Serial
  auto v_com_tally = ScatterViewHelper<typename NeedDup<NEED_ATOMICS,DeviceType>::value,decltype(dup_vcom_tally),decltype(ndup_vcom_tally)>::get(dup_vcom_tally,ndup_vcom_tally);
  auto a_vcom_tally = v_com_tally.template access<typename AtomicDup<NEED_ATOMICS,DeviceType>::value>();

  const int ispecies = d_particles[i].ispecies;
  const int igroup = d_s2g(imix,ispecies);
  if (igroup < 0) return;
  const int icell = d_particles[i].icell;
  if (!(d_cinfo[icell].mask & groupbit)) return;

  const double mass = d_species[ispecies].mass;
  double *v = d_particles[i].v;

  a_vcom_tally(icell,igroup,0) += mass * v[0];
  a_vcom_tally(icell,igroup,1) += mass * v[1];
  a_vcom_tally(icell,igroup,2) += mass * v[2];
  a_vcom_tally(icell,igroup,3) += mass;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeSonineGridKokkos::operator()(TagComputeSonineGrid_compute_vcom, const int &icell) const {
  // compute COM velocity on this timestep for each cell and group
  const int np = d_cellcount[icell];

  for (int n = 0; n < np; n++) {

    const int i = d_plist(icell,n);

    const int ispecies = d_particles[i].ispecies;
    const int igroup = d_s2g(imix,ispecies);
    if (igroup < 0) continue;
    if (!(d_cinfo[icell].mask & groupbit)) continue;

    const double mass = d_species[ispecies].mass;
    double *v = d_particles[i].v;

    d_vcom(icell,igroup,0) += mass * v[0];
    d_vcom(icell,igroup,1) += mass * v[1];
    d_vcom(icell,igroup,2) += mass * v[2];
    d_vcom(icell,igroup,3) += mass;
  }

  double norm;
  for (int j=0; j<ngroup; j++) {
    norm = d_vcom(icell,j,3);
    if (norm == 0.0) continue;
    d_vcom(icell,j,0) /= norm;
    d_vcom(icell,j,1) /= norm;
    d_vcom(icell,j,2) /= norm;
  }
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeSonineGridKokkos::operator()(TagComputeSonineGrid_normalize_vcom, const int &icell) const {
  const int ispecies = d_particles[icell].ispecies;
  const int igroup = d_s2g(imix,ispecies);

  double norm;
  for (int j=0; j<ngroup; j++) {
    norm = d_vcom(icell,j,3);
    if (norm == 0.0) continue;
    d_vcom(icell,j,0) /= norm;
    d_vcom(icell,j,1) /= norm;
    d_vcom(icell,j,2) /= norm;
  }
}

/* ---------------------------------------------------------------------- */

template<int NEED_ATOMICS>
KOKKOS_INLINE_FUNCTION
void ComputeSonineGridKokkos::operator()(TagComputeSonineGrid_compute_per_grid_atomic<NEED_ATOMICS>, const int &i) const {

  // The tally array is duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_tally = ScatterViewHelper<typename NeedDup<NEED_ATOMICS,DeviceType>::value,decltype(dup_tally),decltype(ndup_tally)>::get(dup_tally,ndup_tally);
  auto a_tally = v_tally.template access<typename AtomicDup<NEED_ATOMICS,DeviceType>::value>();

  const int ispecies = d_particles[i].ispecies;
  const int igroup = d_s2g(imix,ispecies);
  if (igroup < 0) return;
  const int icell = d_particles[i].icell;
  if (!(d_cinfo[icell].mask & groupbit)) return;

  int k = igroup*npergroup;

  const double mass = d_species[ispecies].mass;
  double *v = d_particles[i].v;
  a_tally(icell,k++) += mass;

  double vthermal[3];
  double csq;
  vthermal[0] = v[0] - d_vcom(icell,igroup,0);
  vthermal[1] = v[1] - d_vcom(icell,igroup,1);
  vthermal[2] = v[2] - d_vcom(icell,igroup,2);
  csq = vthermal[0]*vthermal[0] + vthermal[1]*vthermal[1] + vthermal[2]*vthermal[2];

  double value;
  for (int m=0; m<nvalue; m++) {
    if (d_which[m] == AMOM) {
      value = mass*vthermal[d_moment[m]] * csq;
      a_tally(icell,k++) += value;
      for (int n=1; n<d_order[m]; n++) {
        value *= csq;
        a_tally(icell,k++) += value;
      }
    } else {
      value = mass * vthermal[d_moment[m]/3] * vthermal[d_moment[m]%3] * csq;
      a_tally(icell,k++) += value;
      for (int n=1; n<d_order[m]; n++) {
        value *= csq;
        a_tally(icell,k++) += value;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeSonineGridKokkos::operator()(TagComputeSonineGrid_compute_per_grid, const int &icell) const {

  const int np = d_cellcount[icell];

  for (int n = 0; n < np; n++) {

    const int i = d_plist(icell,n);

    const int ispecies = d_particles[i].ispecies;
    const int igroup = d_s2g(imix,ispecies);
    if (igroup < 0) continue;
    if (!(d_cinfo[icell].mask & groupbit)) continue;

    const double mass = d_species[ispecies].mass;
    double *v = d_particles[i].v;

    int k = igroup*npergroup;

    double vthermal[3];
    double csq;
    vthermal[0] = v[0] - d_vcom(icell,igroup,0);
    vthermal[1] = v[1] - d_vcom(icell,igroup,1);
    vthermal[2] = v[2] - d_vcom(icell,igroup,2);
    csq = vthermal[0]*vthermal[0] + vthermal[1]*vthermal[1] + vthermal[2]*vthermal[2];

    double value;
    for (int m=0; m<nvalue; m++) {
      if (d_which[m] == AMOM) {
        value = mass*vthermal[d_moment[m]] * csq;
        d_tally(icell,k++) += value;
        for (int n=1; n<d_order[m]; n++) {
          value *= csq;
          d_tally(icell,k++) += value;
        }
      } else {
        value = mass * vthermal[d_moment[m]/3] * vthermal[d_moment[m]%3] * csq;
        d_tally(icell,k++) += value;
        for (int n=1; n<d_order[m]; n++) {
          value *= csq;
          d_tally(icell,k++) += value;
        }
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

int ComputeSonineGridKokkos::query_tally_grid_kokkos(DAT::t_float_2d_lr &d_array)
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

void ComputeSonineGridKokkos::post_process_grid_kokkos(int index,
                                                         int nsample,
                                                         DAT::t_float_2d_lr d_etally,
                                                         int *emap,
                                                         DAT::t_float_1d_strided d_vec)
{
  index--;

  int lo = 0;
  int hi = nglocal;

  if (!d_etally.data()) {
    d_etally = d_tally;
    emap = map[index];
    d_vec = d_vector;
  }

  this->d_etally = d_etally;
  this->d_vec = d_vec;
  numerator = emap[0];
  mass = emap[1];

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeSonineGrid_post_process_grid>(lo,hi),*this);
  DeviceType().fence();
  copymode = 0;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeSonineGridKokkos::operator()(TagComputeSonineGrid_post_process_grid, const int &icell) const
{
  double norm = d_etally(icell,mass);
  if (norm == 0.0) d_vec[icell] = 0.0;
  else d_vec[icell] = d_etally(icell,numerator)/norm;
}

/* ----------------------------------------------------------------------
   reallocate arrays if nglocal has changed
   called by init() and load balancer
------------------------------------------------------------------------- */

void ComputeSonineGridKokkos::reallocate()
{
  if (grid->nlocal == nglocal) return;

  memoryKK->destroy_kokkos(k_vector_grid,vector_grid);
  memoryKK->destroy_kokkos(k_tally,tally);

  nglocal = grid->nlocal;
  memoryKK->create_kokkos(k_vector_grid,vector_grid,nglocal,"sonine/grid:vector_grid");
  d_vector = k_vector_grid.d_view;
  memoryKK->create_kokkos(k_tally,tally,nglocal,ntotal,"sonine/grid:tally");
  d_tally = k_tally.d_view;
  d_vcom = DAT::t_float_3d ("d_vcom",nglocal,ngroup,4);
}

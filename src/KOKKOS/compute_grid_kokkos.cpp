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

#include "string.h"
#include "compute_grid_kokkos.h"
#include "particle_kokkos.h"
#include "mixture.h"
#include "grid_kokkos.h"
#include "update.h"
#include "modify.h"
#include "memory_kokkos.h"
#include "error.h"
#include "sparta_masks.h"
#include "kokkos.h"

using namespace SPARTA_NS;

// user keywords

enum{NUM,NRHO,NFRAC,MASS,MASSRHO,MASSFRAC,
     U,V,W,USQ,VSQ,WSQ,KE,TEMPERATURE,EROT,TROT,EVIB,TVIB,
     PXRHO,PYRHO,PZRHO,KERHO};

// internal accumulators

enum{COUNT,MASSSUM,MVX,MVY,MVZ,MVXSQ,MVYSQ,MVZSQ,MVSQ,
     ENGROT,ENGVIB,DOFROT,DOFVIB,CELLCOUNT,CELLMASS,LASTSIZE};

// max # of quantities to accumulate for any user value

#define MAXACCUMULATE 2

/* ---------------------------------------------------------------------- */

ComputeGridKokkos::ComputeGridKokkos(SPARTA *sparta, int narg, char **arg) :
  ComputeGrid(sparta, narg, arg)
{
  kokkos_flag = 1;

  k_unique = DAT::tdual_int_1d("compute/grid:unique",npergroup);
  for (int m = 0; m < npergroup; m++)
    k_unique.h_view(m) = unique[m];
  k_unique.modify_host();
  k_unique.sync_device();
  d_unique = k_unique.d_view;
}

/* ---------------------------------------------------------------------- */

ComputeGridKokkos::~ComputeGridKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_vector_grid,vector_grid);
  memoryKK->destroy_kokkos(k_tally,tally);
  vector_grid = NULL;
  tally = NULL;
}

/* ---------------------------------------------------------------------- */

void ComputeGridKokkos::compute_per_grid()
{
  if (sparta->kokkos->prewrap) {
    ComputeGrid::compute_per_grid();
  } else {
    compute_per_grid_kokkos();
    k_tally.modify_device();
    k_tally.sync_host();
  }
}

/* ---------------------------------------------------------------------- */

void ComputeGridKokkos::compute_per_grid_kokkos()
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
  // skip cells not in grid group
  // perform all tallies needed for each particle
  // depends on its species group and the user-requested values

  need_dup = sparta->kokkos->need_dup<DeviceType>();
  if (particle_kk->sorted_kk && sparta->kokkos->need_atomics && !sparta->kokkos->atomic_reduction)
    need_dup = 0;

  if (need_dup)
    dup_tally = Kokkos::Experimental::create_scatter_view<typename Kokkos::Experimental::ScatterSum, typename Kokkos::Experimental::ScatterDuplicated>(d_tally);
  else
    ndup_tally = Kokkos::Experimental::create_scatter_view<typename Kokkos::Experimental::ScatterSum, typename Kokkos::Experimental::ScatterNonDuplicated>(d_tally);

  copymode = 1;
  if (particle_kk->sorted_kk && sparta->kokkos->need_atomics && !sparta->kokkos->atomic_reduction)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeGrid_compute_per_grid>(0,nglocal),*this);
  else {
    if (sparta->kokkos->need_atomics)
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeGrid_compute_per_grid_atomic<1> >(0,nlocal),*this);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeGrid_compute_per_grid_atomic<0> >(0,nlocal),*this);
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
void ComputeGridKokkos::operator()(TagComputeGrid_compute_per_grid_atomic<NEED_ATOMICS>, const int &i) const {

  // The tally array is duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_tally = ScatterViewHelper<typename NeedDup<NEED_ATOMICS,DeviceType>::value,decltype(dup_tally),decltype(ndup_tally)>::get(dup_tally,ndup_tally);
  auto a_tally = v_tally.template access<typename AtomicDup<NEED_ATOMICS,DeviceType>::value>();

  const int ispecies = d_particles[i].ispecies;
  const int igroup = d_s2g(imix,ispecies);
  if (igroup < 0) return;

  const int icell = d_particles[i].icell;
  if (!(d_cinfo[icell].mask & groupbit)) return;

  const double mass = d_species[ispecies].mass;
  double* v = d_particles[i].v;

  if (cellmass) a_tally(icell,cellmass) += mass;
  if (cellcount) a_tally(icell,cellcount) += 1.0;

  // loop has all possible values particle needs to accumulate
  // subset defined by user values are indexed by accumulate vector

  int k = igroup*npergroup;

  for (int m = 0; m < npergroup; m++) {
    switch (d_unique[m]) {
    case COUNT:
      a_tally(icell,k++) += 1.0;
      break;
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
    case MVSQ:
      a_tally(icell,k++) += mass * (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
      break;
    case ENGROT:
      a_tally(icell,k++) += d_particles[i].erot;
      break;
    case ENGVIB:
      a_tally(icell,k++) += d_particles[i].evib;
      break;
    case DOFROT:
      a_tally(icell,k++) += d_species[ispecies].rotdof;
      break;
    case DOFVIB:
      a_tally(icell,k++) += d_species[ispecies].vibdof;
      break;
    }
  }
}

KOKKOS_INLINE_FUNCTION
void ComputeGridKokkos::operator()(TagComputeGrid_compute_per_grid, const int &icell) const {
  if (!(d_cinfo[icell].mask & groupbit)) return;

  const int np = d_cellcount[icell];

  for (int n = 0; n < np; n++) {

    const int i = d_plist(icell,n);

    const int ispecies = d_particles[i].ispecies;
    const int igroup = d_s2g(imix,ispecies);
    if (igroup < 0) return;

    const double mass = d_species[ispecies].mass;
    double* v = d_particles[i].v;

    if (cellmass) d_tally(icell,cellmass) += mass;
    if (cellcount) d_tally(icell,cellcount) += 1.0;

    // loop has all possible values particle needs to accumulate
    // subset defined by user values are indexed by accumulate vector

    int k = igroup*npergroup;

    for (int m = 0; m < npergroup; m++) {
      switch (d_unique[m]) {
      case COUNT:
        d_tally(icell,k++) += 1.0;
        break;
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
      case MVSQ:
        d_tally(icell,k++) += mass * (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
        break;
      case ENGROT:
        d_tally(icell,k++) += d_particles[i].erot;
        break;
      case ENGVIB:
        d_tally(icell,k++) += d_particles[i].evib;
        break;
      case DOFROT:
        d_tally(icell,k++) += d_species[ispecies].rotdof;
        break;
      case DOFVIB:
        d_tally(icell,k++) += d_species[ispecies].vibdof;
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

int ComputeGridKokkos::query_tally_grid_kokkos(DAT::t_float_2d_lr &d_array)
{
  d_array = d_tally;
  return 0;
}

/* ----------------------------------------------------------------------
   tally accumulated info to compute final normalized values
   index = which column of output (0 for vec, 1 to N for array)
   for etally = NULL:
     use internal tallied info for single timestep, set nsample = 1
     if onecell = -1, compute values for all grid cells
       store results in vector_grid with nstride = 1 (single col of array_grid)
     if onecell >= 0, compute single value for onecell and return it
   for etally = ptr to caller array:
     use external tallied info for many timesteps
     nsample = additional normalization factor used by some values
     emap = list of etally columns to use, # of columns determined by index
     store results in caller's vec, spaced by nstride
   if norm = 0.0, set result to 0.0 directly so do not divide by 0.0
------------------------------------------------------------------------- */

void ComputeGridKokkos::post_process_grid_kokkos(int index, int nsample,
                                      DAT::t_float_2d_lr d_etally, int *emap,
                                      DAT::t_float_1d_strided d_vec)
{
  index--;
  int ivalue = index % nvalue;

  int lo = 0;
  int hi = nglocal;

  if (!d_etally.data()) {
    nsample = 1;
    d_etally = d_tally;
    emap = map[index];
    d_vec = d_vector;
    nstride = 1;
  }

  this->nsample = nsample;
  this->d_etally = d_etally;
  this->d_vec = d_vec;
  this->nstride = nstride;

  GridKokkos* grid_kk = (GridKokkos*) grid;
  grid_kk->sync(Device,CINFO_MASK);
  d_cinfo = grid_kk->k_cinfo.d_view;

  fnum = update->fnum;

  copymode = 1;

  // compute normalized final value for each grid cell

  switch (value[ivalue]) {

  case NUM:
    {
      count = emap[0];
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeGrid_NUM>(lo,hi),*this);
      DeviceType().fence();
      break;
    }

  case MASS:
    {
      mass = emap[0];
      count = emap[1];
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeGrid_MASS>(lo,hi),*this);
      DeviceType().fence();
      break;
    }

  case NRHO:
    {
      count = emap[0];
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeGrid_NRHO>(lo,hi),*this);
      DeviceType().fence();
      break;
    }

  case MASSRHO:
    {
      mass = emap[0];
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeGrid_MASSRHO>(lo,hi),*this);
      DeviceType().fence();
      break;
    }

  case NFRAC:
  case MASSFRAC:
    {
      count_or_mass = emap[0];
      cell_count_or_mass = emap[1];
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeGrid_NFRAC>(lo,hi),*this);
      DeviceType().fence();
      break;
    }

  case U:
  case V:
  case W:
  case USQ:
  case VSQ:
  case WSQ:
    {
      velocity = emap[0];
      mass = emap[1];
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeGrid_U>(lo,hi),*this);
      DeviceType().fence();
      break;
    }

  case KE:
    {
      mvsq = emap[0];
      count = emap[1];
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeGrid_KE>(lo,hi),*this);
      DeviceType().fence();
      break;
    }

  case TEMPERATURE:
    {
      mvsq = emap[0];
      count = emap[1];
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeGrid_TEMPERATURE>(lo,hi),*this);
      DeviceType().fence();
      break;
    }

  case EROT:
  case EVIB:
    {
      eng = emap[0];
      count = emap[1];
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeGrid_EROT>(lo,hi),*this);
      DeviceType().fence();
      break;
    }

  case TROT:
  case TVIB:
    {
      eng = emap[0];
      dof = emap[1];
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeGrid_TROT>(lo,hi),*this);
      DeviceType().fence();
      break;
    }

  case PXRHO:
  case PYRHO:
  case PZRHO:
    {
      mom = emap[0];
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeGrid_PXRHO>(lo,hi),*this);
      DeviceType().fence();
      break;
    }

  case KERHO:
    {
      ke = emap[0];
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeGrid_KERHO>(lo,hi),*this);
      DeviceType().fence();
      break;
    }
  }
  copymode = 0;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeGridKokkos::operator()(TagComputeGrid_NUM, const int &icell) const {
  d_vec[icell] = d_etally(icell,count) / nsample;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeGridKokkos::operator()(TagComputeGrid_MASS, const int &icell) const {
  const double norm = d_etally(icell,count);
  if (norm == 0.0) d_vec[icell] = 0.0;
  else d_vec[icell] = d_etally(icell,mass) / norm;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeGridKokkos::operator()(TagComputeGrid_NRHO, const int &icell) const {
  const double norm = d_cinfo[icell].volume;
  if (norm == 0.0) d_vec[icell] = 0.0;
  else {
    const double wt = fnum * d_cinfo[icell].weight / norm;
    d_vec[icell] = wt * d_etally(icell,count) / nsample;
  }
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeGridKokkos::operator()(TagComputeGrid_MASSRHO, const int &icell) const {
  const double norm = d_cinfo[icell].volume;
  if (norm == 0.0) d_vec[icell] = 0.0;
  else {
    const double wt = fnum * d_cinfo[icell].weight / norm;
    d_vec[icell] = wt * d_etally(icell,mass) / nsample;
  }
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeGridKokkos::operator()(TagComputeGrid_NFRAC, const int &icell) const {
  const double norm = d_etally(icell,cell_count_or_mass);
  if (norm == 0.0) d_vec[icell] = 0.0;
  else d_vec[icell] = d_etally(icell,count_or_mass) / norm;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeGridKokkos::operator()(TagComputeGrid_U, const int &icell) const {
  const double norm = d_etally(icell,mass);
  if (norm == 0.0) d_vec[icell] = 0.0;
  else d_vec[icell] = d_etally(icell,velocity) / norm;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeGridKokkos::operator()(TagComputeGrid_KE, const int &icell) const {
  const double norm = d_etally(icell,count);
  if (norm == 0.0) d_vec[icell] = 0.0;
  else d_vec[icell] = eprefactor * d_etally(icell,mvsq) / norm;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeGridKokkos::operator()(TagComputeGrid_TEMPERATURE, const int &icell) const {
  const double norm = d_etally(icell,count);
  if (norm == 0.0) d_vec[icell] = 0.0;
  else d_vec[icell] = tprefactor * d_etally(icell,mvsq) / norm;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeGridKokkos::operator()(TagComputeGrid_EROT, const int &icell) const {
  const double norm = d_etally(icell,count);
  if (norm == 0.0) d_vec[icell] = 0.0;
  else d_vec[icell] = d_etally(icell,eng) / norm;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeGridKokkos::operator()(TagComputeGrid_TROT, const int &icell) const {
  const double norm = d_etally(icell,dof);
  if (norm == 0.0) d_vec[icell] = 0.0;
  else d_vec[icell] = rvprefactor * d_etally(icell,eng) / norm;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeGridKokkos::operator()(TagComputeGrid_PXRHO, const int &icell) const {
  const double norm = d_cinfo[icell].volume;
  if (norm == 0.0) d_vec[icell] = 0.0;
  else {
    const double wt = fnum * d_cinfo[icell].weight / norm;
    d_vec[icell] = wt * d_etally(icell,mom) / nsample;
  }
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeGridKokkos::operator()(TagComputeGrid_KERHO, const int &icell) const {
  const double norm = d_cinfo[icell].volume;
  if (norm == 0.0) d_vec[icell] = 0.0;
  else {
    const double wt = fnum * d_cinfo[icell].weight / norm;
    d_vec[icell] = eprefactor * wt * d_etally(icell,ke) / nsample;
  }
}

/* ----------------------------------------------------------------------
   reallocate data storage if nglocal has changed
   called by init() and whenever grid changes
------------------------------------------------------------------------- */

void ComputeGridKokkos::reallocate()
{
  if (grid->nlocal == nglocal) return;

  memoryKK->destroy_kokkos(k_vector_grid,vector_grid);
  memoryKK->destroy_kokkos(k_tally,tally);
  nglocal = grid->nlocal;
  memoryKK->create_kokkos(k_vector_grid,vector_grid,nglocal,"compute_grid:vector_grid");
  d_vector = k_vector_grid.d_view;
  memoryKK->create_kokkos(k_tally,tally,nglocal,ntotal,"compute_grid:tally");
  d_tally = k_tally.d_view;
}

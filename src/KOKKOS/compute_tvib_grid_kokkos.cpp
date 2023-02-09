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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "compute_tvib_grid_kokkos.h"
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

enum{NONE,COUNT,MASSWT,DOF};

/* ---------------------------------------------------------------------- */

ComputeTvibGridKokkos::ComputeTvibGridKokkos(SPARTA *sparta, int narg, char **arg) :
  ComputeTvibGrid(sparta, narg, arg)
{
  kokkos_flag = 1;

  if (modeflag == 0) {
    k_t2s = DAT::tdual_int_1d("compute/tvib/grid:t2s",ntally);
    k_s2t = DAT::tdual_int_1d("compute/tvib/grid:s2t",nspecies);
    d_tspecies = DAT::t_float_1d("d_tspecies",nspecies);

    for (int n = 0; n < nspecies; n++)
      k_s2t.h_view(n) = s2t[n];

    for (int n = 0; n < ntally; n++)
      k_t2s.h_view(n) = t2s[n];

    k_s2t.modify_host();
    k_t2s.modify_host();

    k_s2t.sync_device();
    k_t2s.sync_device();

    d_s2t = k_s2t.d_view;
    d_t2s = k_t2s.d_view;
  } else {
    k_t2s_mode = DAT::tdual_int_1d("compute/tvib/grid:t2s_mode",ntally);
    k_s2t_mode = DAT::tdual_int_2d("compute/tvib/grids2t_mode",nspecies,maxmode);
    d_tspecies_mode = DAT::t_float_2d_lr("d_tspecies_mode",nspecies,maxmode);

    for (int n = 0; n < nspecies; n++)
      for (int m = 0; m < maxmode; m++)
        k_s2t_mode.h_view(n,m) = s2t_mode[n][m];

    for (int n = 0; n < ntally; n++)
      k_t2s_mode.h_view(n) = t2s_mode[n];

    k_s2t_mode.modify_host();
    k_t2s_mode.modify_host();

    k_s2t_mode.sync_device();
    k_t2s_mode.sync_device();

    d_s2t_mode = k_s2t_mode.d_view;
    d_t2s_mode = k_t2s_mode.d_view;
  }

  boltz = update->boltz;
}

/* ---------------------------------------------------------------------- */

ComputeTvibGridKokkos::~ComputeTvibGridKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_vector_grid,vector_grid);
  memoryKK->destroy_kokkos(k_tally,tally);
  vector_grid = NULL;
  tally = NULL;
}

/* ---------------------------------------------------------------------- */

void ComputeTvibGridKokkos::compute_per_grid()
{
  if (sparta->kokkos->prewrap) {
    ComputeTvibGrid::compute_per_grid();
  } else {
    compute_per_grid_kokkos();
    k_tally.modify_device();
    k_tally.sync_host();
  }
}

/* ---------------------------------------------------------------------- */

void ComputeTvibGridKokkos::compute_per_grid_kokkos()
{
  invoked_per_grid = update->ntimestep;

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,PARTICLE_MASK|SPECIES_MASK|CUSTOM_MASK);
  d_particles = particle_kk->k_particles.d_view;
  d_species = particle_kk->k_species.d_view;
  d_ewhich = particle_kk->k_ewhich.d_view;
  k_eiarray = particle_kk->k_eiarray;

  GridKokkos* grid_kk = (GridKokkos*) grid;
  d_cellcount = grid_kk->d_cellcount;
  d_plist = grid_kk->d_plist;
  grid_kk->sync(Device,CINFO_MASK);
  d_cinfo = grid_kk->k_cinfo.d_view;

  d_s2g = particle_kk->k_species2group.view<DeviceType>();
  int nlocal = particle->nlocal;

  // zero all accumulators

  Kokkos::deep_copy(d_tally,0.0);

  // loop over all particles, skip species not in mixture group
  // mode = 0: tally vib eng and count for each species
  // mode >= 1: tally vib level and count for each species and each vib mode

  need_dup = sparta->kokkos->need_dup<DeviceType>();
  if (particle_kk->sorted_kk && sparta->kokkos->need_atomics && !sparta->kokkos->atomic_reduction)
    need_dup = 0;

  if (need_dup)
    dup_tally = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_tally);
  else
    ndup_tally = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_tally);

  copymode = 1;
  if (particle_kk->sorted_kk && sparta->kokkos->need_atomics && !sparta->kokkos->atomic_reduction)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeTvibGrid_compute_per_grid>(0,nglocal),*this);
  else {
    if (sparta->kokkos->need_atomics)
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeTvibGrid_compute_per_grid_atomic<1> >(0,nlocal),*this);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeTvibGrid_compute_per_grid_atomic<0> >(0,nlocal),*this);
  }
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
void ComputeTvibGridKokkos::operator()(TagComputeTvibGrid_compute_per_grid_atomic<NEED_ATOMICS>, const int &i) const {
  // The tally array is duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_tally = ScatterViewHelper<typename NeedDup<NEED_ATOMICS,DeviceType>::value,decltype(dup_tally),decltype(ndup_tally)>::get(dup_tally,ndup_tally);
  auto a_tally = v_tally.template access<typename AtomicDup<NEED_ATOMICS,DeviceType>::value>();

  const int ispecies = d_particles[i].ispecies;
  if (!d_species[ispecies].vibdof) return;
  const int igroup = d_s2g(imix,ispecies);
  if (igroup < 0) return;
  const int icell = d_particles[i].icell;
  if (!(d_cinfo[icell].mask & groupbit)) return;

  if (modeflag == 0) {
    const int j = d_s2t[ispecies];
    a_tally(icell,j) += d_particles[i].evib;
    a_tally(icell,j+1) += 1.0;
  } else if (modeflag >= 1) {
    auto &d_vibmode = k_eiarray.d_view[d_ewhich[index_vibmode]].k_view.d_view;

    // tally only the modes this species has

    const int nmode = d_species[ispecies].nvibmode;
    for (int imode = 0; imode < nmode; imode++) {
      const int j = d_s2t_mode(ispecies,imode);
      if (nmode > 1) a_tally(icell,j) += d_vibmode(i,imode);
      else a_tally(icell,j) +=
             d_particles[i].evib / (boltz*d_species[ispecies].vibtemp[0]);
      a_tally(icell,j+1) += 1.0;
    }
  }
}

KOKKOS_INLINE_FUNCTION
void ComputeTvibGridKokkos::operator()(TagComputeTvibGrid_compute_per_grid, const int &icell) const {
  if (!(d_cinfo[icell].mask & groupbit)) return;

  const int np = d_cellcount[icell];

  for (int n = 0; n < np; n++) {
    const int i = d_plist(icell,n);
    const int ispecies = d_particles[i].ispecies;
    if (!d_species[ispecies].vibdof) continue;
    const int igroup = d_s2g(imix,ispecies);
    if (igroup < 0) continue;

    if (modeflag == 0) {
      const int j = d_s2t[ispecies];
      d_tally(icell,j) += d_particles[i].evib;
      d_tally(icell,j+1) += 1.0;
    } else if (modeflag >= 1) {
      auto &d_vibmode = k_eiarray.d_view[d_ewhich[index_vibmode]].k_view.d_view;

      // tally only the modes this species has

      const int nmode = d_species[ispecies].nvibmode;
      for (int imode = 0; imode < nmode; imode++) {
        const int j = d_s2t_mode(ispecies,imode);
        if (nmode > 1) d_tally(icell,j) += d_vibmode(i,imode);
        else d_tally(icell,j) +=
               d_particles[i].evib / (boltz*d_species[ispecies].vibtemp[0]);
        d_tally(icell,j+1) += 1.0;
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

int ComputeTvibGridKokkos::query_tally_grid_kokkos(DAT::t_float_2d_lr &d_array)
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
   for etaylly = ptr to caller array:
     use external tallied info for many timesteps
     nsample = additional normalization factor used by some values
     emap = list of etally columns to use, # of columns determined by index
     store results in caller's vec, spaced by nstride
   if norm = 0.0, set result to 0.0 directly so do not divide by 0.0
------------------------------------------------------------------------- */

void ComputeTvibGridKokkos::
post_process_grid_kokkos(int index, int nsample,
                         DAT::t_float_2d_lr d_etally, int *emap,
                         DAT::t_float_1d_strided d_vec)
{
  index--;

  int lo = 0;
  int hi = nglocal;

  if (!d_etally.data()) {
    nsample = 1;
    d_etally = d_tally;
    emap = map[index];
    d_vec = d_vector;
    nstride = 1;
  }

  this->d_etally = d_etally;
  this->d_vec = d_vec;

  // compute normalized single Tgroup value for each grid cell
  // compute Ibar and Tspecies for each species in the requested group
  // ditto for individual vibrational modes if modeflag = 1 or 2
  // loop over species/modes in group to compute normalized Tgroup

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,SPECIES_MASK|CUSTOM_MASK);
  d_species = particle_kk->k_species.d_view;

  if (modeflag == 0) {
    nsp = nmap[index] / 2;
    evib = emap[0];
    count = emap[1];
  } else if (modeflag == 1) {
    nsp = nmap[index] / maxmode / 2;
    evib = emap[0];
    count = emap[1];
  } else if (modeflag == 2) {
    nsp = nmap[index] / maxmode / 2;
    imode = index % maxmode;
    evib = emap[2*imode];
    count = emap[2*imode+1];
  }

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeTvibGrid_post_process_grid>(lo,hi),*this);
  copymode = 0;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeTvibGridKokkos::operator()(TagComputeTvibGrid_post_process_grid, const int &icell) const {
  int evb = evib;
  int cnt = evib+1;

  // modeflag = 0, no vib modes exist
  // nsp = # of species in the group
  // inputs: 2*nsp tallies
  // output: Tgroup = weighted sum over all Tsp for species in group

  if (modeflag == 0) {

    for (int isp = 0; isp < nsp; ++isp) {
      const int ispecies = d_t2s[evb-evib];
      const double theta = d_species[ispecies].vibtemp[0];
      if (theta == 0.0 || d_etally(icell,cnt) == 0.0) {
        d_tspecies[isp] = 0.0;
        evb += 2;
        cnt = evb+1;
        continue;
      }
      const double ibar = d_etally(icell,evb) / (d_etally(icell,cnt) * boltz * theta);
      if (ibar == 0.0) {
        d_tspecies[isp] = 0.0;
        evb += 2;
        cnt = evb+1;
        continue;
      }
      d_tspecies[isp] = theta / (log(1.0 + 1.0/ibar));
      //denom = boltz * etally[icell][count] * ibar * log(1.0 + 1.0/ibar);
      //tspecies[isp] = etally[icell][evb] / denom;
      evb += 2;
      cnt = evb+1;
    }

    double numer = 0.0;
    double denom = 0.0;
    cnt = count;
    for (int isp = 0; isp < nsp; isp++) {
      numer += d_tspecies[isp]*d_etally(icell,cnt);
      denom += d_etally(icell,cnt);
      cnt += 2;
    }

    if (denom == 0.0) d_vec[icell] = 0.0;
    else d_vec[icell] = numer/denom;

  // modeflag = 1, vib modes exist
  // Tgroup = weighted sum over all Tsp and modes for species in group
  // nsp = # of species in the group
  // maxmode = max # of modes for any species (unused values are zero)
  // inputs: 2*nsp*maxmode tallies

  } else if (modeflag == 1) {
    auto &d_vibmode = k_eiarray.d_view[d_ewhich[index_vibmode]].k_view.d_view;

    for (int isp = 0; isp < nsp; isp++) {
      const int ispecies = d_t2s_mode[evb-evib];
      for (int imode = 0; imode < maxmode; imode++) {
        const double theta = d_species[ispecies].vibtemp[imode];
        if (theta == 0.0 || d_etally(icell,cnt) == 0.0) {
          d_tspecies_mode(isp,imode) = 0.0;
          evb += 2;
          cnt = evb+1;
          continue;
        }
        const double ibar = d_etally(icell,evb) / d_etally(icell,cnt);
        if (ibar == 0.0) {
          d_tspecies_mode(isp,imode) = 0.0;
          evb += 2;
          cnt = evb+1;
          continue;
        }
        d_tspecies_mode(isp,imode) = theta / (log(1.0 + 1.0/ibar));
        //denom = boltz * etally[icell][count] * ibar * log(1.0 + 1.0/ibar);
        //tspecies_mode[isp][imode] = etally[icell][evib] / denom;
        evb += 2;
        cnt = evb+1;
      }
    }

    // loop over species in group and all their modes
    // to accumulate numerator & denominator

    double numer = 0.0;
    double denom = 0.0;
    cnt = count;
    for (int isp = 0; isp < nsp; isp++) {
      for (int imode = 0; imode < maxmode; imode++) {
        numer += d_tspecies_mode(isp,imode)*d_etally(icell,cnt);
        denom += d_etally(icell,cnt);
        cnt += 2;
      }
    }
    if (denom == 0.0) d_vec[icell] = 0.0;
    else d_vec[icell] = numer/denom;

  // modeflag = 2, vib modes exist
  // Tgroup = weighted sum over all Tsp and single mode for species in group
  // nsp = # of species in the group
  // imode = single mode correpsonding to caller index
  // inputs: 2*nsp tallies strided by maxmode

  } else if (modeflag == 2) {
    auto &d_vibmode = k_eiarray.d_view[d_ewhich[index_vibmode]].k_view.d_view;

    for (int isp = 0; isp < nsp; isp++) {
      const int ispecies = d_t2s_mode[evb-evib];
      const double theta = d_species[ispecies].vibtemp[imode];
      if (theta == 0.0 || d_etally(icell,cnt) == 0.0) {
        d_tspecies_mode(isp,imode) = 0.0;
        evb += 2*maxmode;
        cnt = evb+1;
        continue;
      }
      const double ibar = d_etally(icell,evb) / d_etally(icell,cnt);
      if (ibar == 0.0) {
        d_tspecies_mode(isp,imode) = 0.0;
        evb += 2*maxmode;
        cnt = evb+1;
        continue;
      }
      d_tspecies_mode(isp,imode) = theta / (log(1.0 + 1.0/ibar));
      //denom = boltz * etally[icell][count] * ibar * log(1.0 + 1.0/ibar);
      //tspecies_mode[isp][imode] = etally[icell][evib] / denom;
      evb += 2*maxmode;
      cnt = evb+1;
    }

    // loop over species in group and single mode for each species
    // to accumulate numerator & denominator

    double numer = 0.0;
    double denom = 0.0;
    cnt = count;
    for (int isp = 0; isp < nsp; isp++) {
      numer += d_tspecies_mode(isp,imode)*d_etally(icell,cnt);
      denom += d_etally(icell,cnt);
      cnt += 2*maxmode;
    }

    if (denom == 0.0) d_vec[icell] = 0.0;
    else d_vec[icell] = numer/denom;
  }

}

/* ----------------------------------------------------------------------
   reallocate arrays if nglocal has changed
   called by init() and whenever grid changes
------------------------------------------------------------------------- */

void ComputeTvibGridKokkos::reallocate()
{
  if (grid->nlocal == nglocal) return;

  memoryKK->destroy_kokkos(k_vector_grid,vector_grid);
  memoryKK->destroy_kokkos(k_tally,tally);
  nglocal = grid->nlocal;
  memoryKK->create_kokkos(k_vector_grid,vector_grid,nglocal,"tvib/grid:vector_grid");
  d_vector = k_vector_grid.d_view;
  memoryKK->create_kokkos(k_tally,tally,nglocal,ntally,"tvib/grid:tally");
  d_tally = k_tally.d_view;
}

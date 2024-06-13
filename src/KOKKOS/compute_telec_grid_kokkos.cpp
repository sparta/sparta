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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "compute_telec_grid_kokkos.h"
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

/* ---------------------------------------------------------------------- */

ComputeTelecGridKokkos::ComputeTelecGridKokkos(SPARTA *sparta, int narg, char **arg) :
  ComputeTelecGrid(sparta, narg, arg)
{
  kokkos_flag = 1;

  int **groupspecies = particle->mixture[imix]->groupspecies;
  int *groupsize = particle->mixture[imix]->groupsize;
  int mixspecies = particle->mixture[imix]->nspecies;

  d_groupspecies = DAT::t_int_2d("compute/telec/groupspecies",ngroup,mixspecies);
  auto h_groupspecies = Kokkos::create_mirror_view(d_groupspecies);

  for (int igroup = 0; igroup < ngroup; igroup++)
    for (int n = 0; n < groupsize[igroup]; n++)
      h_groupspecies(igroup,n) = groupspecies[igroup][n];

  Kokkos::deep_copy(d_groupspecies,h_groupspecies);

  k_s2t = DAT::tdual_int_1d("compute/telec/grid:s2t",nspecies);
  d_tspecies = DAT::t_float_1d("d_tspecies",nspecies);

  for (int n = 0; n < nspecies; n++)
    k_s2t.h_view(n) = s2t[n];

  k_s2t.modify_host();
  k_s2t.sync_device();

  d_s2t = k_s2t.d_view;

  boltz = update->boltz;
}

/* ---------------------------------------------------------------------- */

ComputeTelecGridKokkos::~ComputeTelecGridKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_vector_grid,vector_grid);
  memoryKK->destroy_kokkos(k_tally,tally);
  vector_grid = NULL;
  tally = NULL;
}

/* ---------------------------------------------------------------------- */

void ComputeTelecGridKokkos::compute_per_grid()
{
  if (sparta->kokkos->prewrap) {
    ComputeTelecGrid::compute_per_grid();
  } else {
    compute_per_grid_kokkos();
    k_tally.modify_device();
    k_tally.sync_host();
  }
}

/* ---------------------------------------------------------------------- */

void ComputeTelecGridKokkos::compute_per_grid_kokkos()
{
  invoked_per_grid = update->ntimestep;

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,PARTICLE_MASK|SPECIES_MASK|CUSTOM_MASK);
  d_particles = particle_kk->k_particles.d_view;
  d_ewhich = particle_kk->k_ewhich.d_view;
  k_edvec = particle_kk->k_edvec;

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
  // tally electronic eng and count for species

  need_dup = sparta->kokkos->need_dup<DeviceType>();
  if (particle_kk->sorted_kk && sparta->kokkos->need_atomics && !sparta->kokkos->atomic_reduction)
    need_dup = 0;

  if (need_dup)
    dup_tally = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_tally);
  else
    ndup_tally = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_tally);

  copymode = 1;
  if (particle_kk->sorted_kk && sparta->kokkos->need_atomics && !sparta->kokkos->atomic_reduction)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeTelecGrid_compute_per_grid>(0,nglocal),*this);
  else {
    if (sparta->kokkos->need_atomics)
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeTelecGrid_compute_per_grid_atomic<1> >(0,nlocal),*this);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeTelecGrid_compute_per_grid_atomic<0> >(0,nlocal),*this);
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
void ComputeTelecGridKokkos::operator()(TagComputeTelecGrid_compute_per_grid_atomic<NEED_ATOMICS>, const int &i) const {
  // The tally array is duplicated for OpenMP, atomic for GPUs, and neither for Serial

  auto v_tally = ScatterViewHelper<typename NeedDup<NEED_ATOMICS,DeviceType>::value,decltype(dup_tally),decltype(ndup_tally)>::get(dup_tally,ndup_tally);
  auto a_tally = v_tally.template access<typename AtomicDup<NEED_ATOMICS,DeviceType>::value>();

  const int ispecies = d_particles[i].ispecies;
  const int igroup = d_s2g(imix,ispecies);
  if (igroup < 0) return;
  const int icell = d_particles[i].icell;
  if (!(d_cinfo[icell].mask & groupbit)) return;

  const int j = d_s2t[ispecies];
  auto &d_eelecs = k_edvec.d_view[d_ewhich[index_eelec]].k_view.d_view;
  a_tally(icell,j) += d_eelecs[i];
  a_tally(icell,j+1) += 1.0;
}

KOKKOS_INLINE_FUNCTION
void ComputeTelecGridKokkos::operator()(TagComputeTelecGrid_compute_per_grid, const int &icell) const {
  if (!(d_cinfo[icell].mask & groupbit)) return;

  const int np = d_cellcount[icell];

  for (int n = 0; n < np; n++) {
    const int i = d_plist(icell,n);
    const int ispecies = d_particles[i].ispecies;
    const int igroup = d_s2g(imix,ispecies);
    if (igroup < 0) continue;

    const int j = d_s2t[ispecies];
    auto &d_eelecs = k_edvec.d_view[d_ewhich[index_eelec]].k_view.d_view;
    d_tally(icell,j) += d_eelecs[i];
    d_tally(icell,j+1) += 1.0;
  }
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
double ComputeTelecGridKokkos::elec_energy(int icell, int isp, double temp_elec) const
{
  // copied from Particle::electronic_distribution_func()

  auto &d_distribution = d_cumulative_probabilities;
  double partition_function = 0.0;

  for (int i = 0; i < d_nelecstates[isp]; ++i) {

    // Calculate boltzmann fractions

    d_distribution(icell,i) = d_elecstates(isp,i).degen*exp(-d_elecstates(isp,i).temp/temp_elec);

    // Calculate partition function

    partition_function += d_distribution(icell,i);
  }

  for (int i = 0; i < d_nelecstates[isp]; ++i)
    d_distribution(icell,i) /= partition_function;

  auto &d_state_probabilities = d_cumulative_probabilities;

  double total_energy = 0.0;
  for (int i = 0; i < d_nelecstates[isp]; ++i)
    total_energy += d_state_probabilities(icell,i)*d_elecstates(isp,i).temp*boltz;

  return total_energy;
}

/* ----------------------------------------------------------------------
   query info about internal tally array for this compute
   index = which column of output (0 for vec, 1 to N for array)
   return # of tally quantities for this index
   also return array = ptr to tally array
   also return cols = ptr to list of columns in tally for this index
------------------------------------------------------------------------- */

int ComputeTelecGridKokkos::query_tally_grid_kokkos(DAT::t_float_2d_lr &d_array)
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
   for etally = ptr to caller array:
     use external tallied info for many timesteps
     nsample = additional normalization factor used by some values
     emap = list of etally columns to use, # of columns determined by index
     store results in caller's vec, spaced by nstride
   if norm = 0.0, set result to 0.0 directly so do not divide by 0.0
------------------------------------------------------------------------- */

void ComputeTelecGridKokkos::
post_process_grid_kokkos(int index, int /*nsample*/,
                         DAT::t_float_2d_lr d_etally, int *emap,
                         DAT::t_float_1d_strided d_vec)
{
  index--;
  this->index = index;

  int lo = 0;
  int hi = nglocal;

  if (!d_etally.data()) {
    d_etally = d_tally;
    emap = map[index];
    d_vec = d_vector_grid;
    nstride = 1;
  }

  this->d_etally = d_etally;
  this->d_vec = d_vec;

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,SPECIES_MASK|CUSTOM_MASK);
  d_species = particle_kk->k_species.d_view;
  d_nelecstates = particle_kk->d_nelecstates;
  d_elecstates = particle_kk->d_elecstates;

  if (grid->maxlocal > (int)d_cumulative_probabilities.extent(0))
    MemKK::realloc_kokkos(d_cumulative_probabilities,"collide:cumulative_probabilities",grid->maxlocal,particle->maxelecstate);

  nsp = nmap[index] / 2;
  eelec = emap[0];
  count = emap[1];

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeTelecGrid_post_process_grid>(lo,hi),*this);
  copymode = 0;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeTelecGridKokkos::operator()(TagComputeTelecGrid_post_process_grid, const int &icell) const {
  int eelc = eelec;
  int cnt = eelec+1;

  auto &d_eelecs = k_edvec.d_view[d_ewhich[index_eelec]].k_view.d_view;

  for (int isp = 0; isp < nsp; isp++) {
    const int ispecies = d_groupspecies(index,isp);
    if (d_nelecstates[isp] == 0 ||
        d_etally(icell,eelc) == 0.0) {
      d_tspecies[isp] = 0.0;
      eelc += 2;
      cnt = eelc+1;
      return; //////
    }

    // We calculate a first guess at the temp assuming
    //   all the electronic energy is stored in the first
    //   excited state

    const double first_elec_eng = d_species[ispecies].elecdat->states[1].temp*boltz;
    const double degen0 = d_elecstates(ispecies,0).degen;
    const double degen1 = d_elecstates(ispecies,1).degen;
    double t_elec = first_elec_eng / (boltz*(
        - log(d_etally(icell,eelc)*degen0 /
            (d_etally(icell,cnt)*first_elec_eng*degen1)
           ))
      );

    // If the electronic excitation is very high, the above calculation will
    //   give negative numbers, including -inf. Negative temperatures are physically
    //   meaningful if over 50% of particles are in excited states. However, in the
    //   case of a truly broken value (-inf) we initialize the guess at something
    //   more sane

    if (isinf(t_elec)) t_elec = -d_elecstates(ispecies,1).temp;

    // Bisection method to find T accurate to 1%

    double target_energy_per_part = d_etally(icell,eelc)/d_etally(icell,cnt);

    // Find initial bounds based on our first guess

    double T_low = 0.9*t_elec;
    while (elec_energy(icell, isp, T_low) > target_energy_per_part) {
      T_low /= 2.0;
    }

    double T_high = t_elec*1.1;
    while (elec_energy(icell, isp, T_high) < target_energy_per_part) {
      T_high *= 2.0;
    }

    // Bisect

    double T_mid = t_elec;
    double e_mid = elec_energy(icell, isp, T_mid);
    while ((T_high - T_low) > 0.01) {
      if (e_mid > target_energy_per_part) {
        T_high = T_mid;
      } else {
        T_low = T_mid;
      }
      T_mid = (T_high - T_low)/2.0 + T_low;
      e_mid = elec_energy(icell, isp, T_mid);
    }

    d_tspecies[isp] = T_mid;
    eelc += 2;
    cnt = eelc+1;
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
}

/* ----------------------------------------------------------------------
   reallocate arrays if nglocal has changed
   called by init() and whenever grid changes
------------------------------------------------------------------------- */

void ComputeTelecGridKokkos::reallocate()
{
  if (grid->nlocal == nglocal) return;

  memoryKK->destroy_kokkos(k_vector_grid,vector_grid);
  memoryKK->destroy_kokkos(k_tally,tally);
  nglocal = grid->nlocal;
  memoryKK->create_kokkos(k_vector_grid,vector_grid,nglocal,"telec/grid:vector_grid");
  d_vector_grid = k_vector_grid.d_view;
  memoryKK->create_kokkos(k_tally,tally,nglocal,ntally,"telec/grid:tally");
  d_tally = k_tally.d_view;
}

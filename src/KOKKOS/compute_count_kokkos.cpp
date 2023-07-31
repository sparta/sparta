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

#include "mpi.h"
#include "string.h"
#include "compute_count_kokkos.h"
#include "update.h"
#include "particle.h"
#include "mixture.h"
#include "comm.h"
#include "error.h"
#include "memory_kokkos.h"
#include "particle_kokkos.h"
#include "sparta_masks.h"
#include "kokkos.h"

using namespace SPARTA_NS;

enum{SPECIES,MIXTURE};


/* ---------------------------------------------------------------------- */

ComputeCountKokkos::ComputeCountKokkos(SPARTA *sparta, int narg, char **arg) :
  ComputeCount(sparta, narg, arg)
{
  kokkos_flag = 1;
}

/* ---------------------------------------------------------------------- */

ComputeCountKokkos::~ComputeCountKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_count, count);
  count = NULL;
}


/* ---------------------------------------------------------------------- */

void ComputeCountKokkos::init()
{
  // insure count vector is long enough for species count
  if (maxspecies < particle->nspecies) {
    maxspecies = particle->nspecies;
    memoryKK->destroy_kokkos(k_count,count);
    memoryKK->create_kokkos(k_count,count,maxspecies,"compute/count:count");
    d_count = k_count.d_view;
  }

  // check if the group count in any accessed mixtures has changed
  int warn = 0;
  int err = 0;
  for (int i = 0; i < nvalues; i++) {
    if (spmix[i] == SPECIES) continue;
    int imix = index[i];
    int igroup = indexgroup[i];
    int ngroup = mixgroups[i];
    if (ngroup != particle->mixture[imix]->ngroup) warn = 1;
    if (igroup >= particle->mixture[imix]->ngroup) err = 1;
  }

  if (err)
    error->all(FLERR,
               "Group in mixture used by compute count no longer valid");
  if (warn && comm->me == 0)
    error->warning(FLERR,
                   "Group count in mixture used by compute count has changed");
}

/* ---------------------------------------------------------------------- */

double ComputeCountKokkos::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  per_species_tally_kokkos();
  k_count.modify_device();
  k_count.sync_host();
  for (int m=0; m<maxspecies; ++m)
    count[m] = k_count.h_view(m);

  bigint one;
  if (spmix[0] == SPECIES) one = count[index[0]];
  else {
    int nspecies = particle->mixture[index[0]]->groupsize[indexgroup[0]];
    one = 0;
    for (int m = 0; m < nspecies; m++) {
      int isp = particle->mixture[index[0]]->groupspecies[indexgroup[0]][m];
      one += count[isp];
    }
  }

  bigint sum;
  MPI_Allreduce(&one,&sum,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  scalar = sum;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeCountKokkos::compute_vector()
{
  int i,m;
  invoked_scalar = update->ntimestep;

  per_species_tally_kokkos();
  k_count.modify_device();
  k_count.sync_host();
  for (int m=0; m<maxspecies; ++m)
    count[m] = k_count.h_view(m);

  for (i = 0; i < nvalues; i++) {
    onevec[i] = 0;
    if (spmix[i] == SPECIES) onevec[i] = count[index[i]];
    else {
      int nspecies = particle->mixture[index[i]]->groupsize[indexgroup[i]];
      for (m = 0; m < nspecies; m++) {
        int isp = particle->mixture[index[i]]->groupspecies[indexgroup[i]][m];
        onevec[i] += count[isp];
      }
    }
  }

  MPI_Allreduce(onevec,sumvec,nvalues,MPI_SPARTA_BIGINT,MPI_SUM,world);
  for (i = 0; i < nvalues; i++) vector[i] = sumvec[i];
}

/* ---------------------------------------------------------------------- */

void ComputeCountKokkos::per_species_tally_kokkos()
{
  if (lasttally == update->ntimestep) return;
  lasttally = update->ntimestep;

  Kokkos::deep_copy(d_count,0);

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,PARTICLE_MASK);
  d_particles = particle_kk->k_particles.d_view;

  need_dup = sparta->kokkos->need_dup<DeviceType>();
  if (need_dup)
    dup_count = Kokkos::Experimental::create_scatter_view<typename Kokkos::Experimental::ScatterSum, typename Kokkos::Experimental::ScatterDuplicated>(d_count);
  else
    ndup_count = Kokkos::Experimental::create_scatter_view<typename Kokkos::Experimental::ScatterSum, typename Kokkos::Experimental::ScatterNonDuplicated>(d_count);

  // loop over all particles
  int nlocal = particle->nlocal;
  copymode = 1;
  if (sparta->kokkos->need_atomics)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeCount_per_species_tally_atomic<1> >(0,nlocal),*this);
  else
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeCount_per_species_tally_atomic<0> >(0,nlocal),*this);
  DeviceType().fence();
  copymode = 0;

  if (need_dup) {
    Kokkos::Experimental::contribute(d_count, dup_count);
    dup_count = decltype(dup_count)(); // free duplicated memory
  }

  d_particles = t_particle_1d(); // destroy reference to reduce memory use
}

/* ---------------------------------------------------------------------- */

template<int NEED_ATOMICS>
KOKKOS_INLINE_FUNCTION
void ComputeCountKokkos::operator()(TagComputeCount_per_species_tally_atomic<NEED_ATOMICS>, const int& i) const {

  // The tally (count) array is duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_count = ScatterViewHelper<typename NeedDup<NEED_ATOMICS,DeviceType>::value,decltype(dup_count),decltype(ndup_count)>::get(dup_count,ndup_count);
  auto a_count = v_count.template access<typename AtomicDup<NEED_ATOMICS,DeviceType>::value>();

  a_count[d_particles[i].ispecies] += 1;
}

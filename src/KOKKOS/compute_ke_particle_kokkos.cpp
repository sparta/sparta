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

#include "string.h"
#include "compute_ke_particle_kokkos.h"
#include "particle_kokkos.h"
#include "memory_kokkos.h"
#include "kokkos.h"
#include "update.h"
#include "sparta_masks.h"
#include "modify.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

ComputeKEParticleKokkos::ComputeKEParticleKokkos(SPARTA *sparta, int narg, char **arg) :
  ComputeKEParticle(sparta, narg, arg)
{
  kokkos_flag = 1;
}

/* ---------------------------------------------------------------------- */

ComputeKEParticleKokkos::~ComputeKEParticleKokkos()
{
  if (copymode) return;
  memoryKK->destroy_kokkos(k_vector_particle,vector_particle);
  vector_particle = NULL;
}

/* ---------------------------------------------------------------------- */

void ComputeKEParticleKokkos::compute_per_particle()
{
  if (sparta->kokkos->prewrap) {
    ComputeKEParticle::compute_per_particle();
  }  else {
    compute_per_particle_kokkos();
  }
}

/* ---------------------------------------------------------------------- */

void ComputeKEParticleKokkos::compute_per_particle_kokkos()
{
  // grow ke array (d_vector) if necessary
  if (particle->nlocal > nmax) {
    memoryKK->destroy_kokkos(k_vector_particle,vector_particle);
    nmax = particle->maxlocal;
    memoryKK->create_kokkos(k_vector_particle,vector_particle,nmax,"ke/particle:vector_particle");
    d_vector = k_vector_particle.d_view;
  }

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,PARTICLE_MASK|SPECIES_MASK);
  d_particles = particle_kk->k_particles.d_view;
  d_species = particle_kk->k_species.d_view;
  int nlocal = particle->nlocal;

  // compute kinetic energy for each atom in group
  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,nlocal),*this);
  copymode = 0;

  d_particles = t_particle_1d(); // destroy reference to reduce memory use
}

KOKKOS_INLINE_FUNCTION
void ComputeKEParticleKokkos::operator()(const int &i) const {
  const int ispecies = d_particles[i].ispecies;
  const double mass = d_species[ispecies].mass;
  const double mvv2e = update->mvv2e;
  double *v = d_particles[i].v;
  d_vector[i] = 0.5 * mvv2e * mass * (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

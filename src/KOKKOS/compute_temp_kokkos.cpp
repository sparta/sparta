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
   Contributing author: Tim Fuller (Sandia)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "string.h"
#include "compute_temp_kokkos.h"
#include "update.h"
#include "particle_kokkos.h"
#include "domain.h"
#include "error.h"
#include "kokkos.h"
#include "sparta_masks.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

ComputeTempKokkos::ComputeTempKokkos(SPARTA *sparta, int narg, char **arg) :
  ComputeTemp(sparta, narg, arg)
{
  if (narg != 2) error->all(FLERR,"Illegal compute temp command");
  scalar_flag = 1;
  kokkos_flag = 1;
}

/* ---------------------------------------------------------------------- */

double ComputeTempKokkos::compute_scalar()
{
  if (sparta->kokkos->prewrap) {
    return ComputeTemp::compute_scalar();
  } else {
    return compute_scalar_kokkos();
  }
}

KOKKOS_INLINE_FUNCTION
void ComputeTempKokkos::operator()(const int& i, double& lsum) const {
  double* v = d_particles[i].v;
  const int ispecies = d_particles[i].ispecies;
  const double mass = d_species[ispecies].mass;
  lsum += (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) * mass;
}

double ComputeTempKokkos::compute_scalar_kokkos()
{
  copymode = 1;
  invoked_scalar = update->ntimestep;
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device, PARTICLE_MASK|SPECIES_MASK);
  d_particles = particle_kk->k_particles.d_view;
  d_species = particle_kk->k_species.d_view;
  int nlocal = particle->nlocal;

  double t = 0.0;
  auto range_policy = Kokkos::RangePolicy<DeviceType>(0, nlocal);
  Kokkos::parallel_reduce(range_policy, *this, t);
  copymode = 0;

  d_particles = t_particle_1d(); // destroy reference to reduce memory use

  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);

  bigint n = particle->nlocal;
  MPI_Allreduce(&n,&particle->nglobal,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  if (particle->nglobal == 0) return 0.0;

  // normalize with 3 instead of dim since even 2d has 3 velocity components

  double factor = update->mvv2e / (3.0 * particle->nglobal * update->boltz);
  scalar *= factor;
  return scalar;
}

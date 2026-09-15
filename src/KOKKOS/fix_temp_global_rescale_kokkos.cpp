/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "fix_temp_global_rescale_kokkos.h"
#include "update.h"
#include "particle_kokkos.h"
#include "sparta_masks.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

FixTempGlobalRescaleKokkos::FixTempGlobalRescaleKokkos(SPARTA *sparta, int narg, char **arg) :
  FixTempGlobalRescale(sparta, narg, arg)
{
  kokkos_flag = 1;
  execution_space = Device;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

void FixTempGlobalRescaleKokkos::end_of_step()
{
  if (update->ntimestep % nevery) return;

  // set current t_target

  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;
  double t_target = tstart + delta * (tstop-tstart);

  // t_current = global temperature
  // just return if no particles or t_current = 0.0

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,PARTICLE_MASK|SPECIES_MASK);
  d_particles = particle_kk->k_particles.view_device();
  d_species = particle_kk->k_species.view_device();

  int nlocal = particle->nlocal;

  // 1st pass: t = sum over my particles of mass*(v.v)

  double t = 0.0;

  copymode = 1;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixTempGlobalRescale_reduce>(0,nlocal),*this,t);
  copymode = 0;

  double t_current;
  MPI_Allreduce(&t,&t_current,1,MPI_DOUBLE,MPI_SUM,world);

  bigint n = particle->nlocal;
  MPI_Allreduce(&n,&particle->nglobal,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  if (particle->nglobal == 0 || t_current == 0.0) return;

  double tscale = update->mvv2e / (3.0 * particle->nglobal * update->boltz);
  t_current *= tscale;

  // rescale all particle velocities

  t_target = t_current - fraction*(t_current-t_target);
  vscale = sqrt(t_target/t_current);

  // 2nd pass: rescale velocities of all my particles

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixTempGlobalRescale_scale>(0,nlocal),*this);
  copymode = 0;

  particle_kk->modify(Device,PARTICLE_MASK);
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixTempGlobalRescaleKokkos::operator()(TagFixTempGlobalRescale_reduce,
                                            const int &i, double &t) const {
  const double *v = d_particles[i].v;
  t += (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) *
    d_species[d_particles[i].ispecies].mass;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixTempGlobalRescaleKokkos::operator()(TagFixTempGlobalRescale_scale,
                                            const int &i) const {
  double *v = d_particles[i].v;
  v[0] *= vscale;
  v[1] *= vscale;
  v[2] *= vscale;
}

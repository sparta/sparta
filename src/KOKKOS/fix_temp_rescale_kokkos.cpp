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

#include "fix_temp_rescale_kokkos.h"
#include "update.h"
#include "particle_kokkos.h"
#include "grid_kokkos.h"
#include "sparta_masks.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

FixTempRescaleKokkos::FixTempRescaleKokkos(SPARTA *sparta, int narg, char **arg) :
  FixTempRescale(sparta, narg, arg)
{
  kokkos_flag = 1;
  execution_space = Device;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

void FixTempRescaleKokkos::end_of_step()
{
  if (update->ntimestep % nevery) return;

  // set current t_target

  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;
  t_target = tstart + delta * (tstop-tstart);

  // sort particles by grid cell if needed

  if (!particle->sorted) particle->sort();

  // loop over grid cells and twice over particles in each cell
  // 1st pass: calc thermal temp via same logic as in ComputeThermalGrid
  // 2nd pass: rescale thermal velocity components

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,PARTICLE_MASK|SPECIES_MASK);
  d_particles = particle_kk->k_particles.d_view;
  d_species = particle_kk->k_species.d_view;

  GridKokkos* grid_kk = (GridKokkos*) grid;
  d_cellcount = grid_kk->d_cellcount;
  d_plist = grid_kk->d_plist;
  grid_kk->sync(Device,CINFO_MASK);

  int nglocal = grid->nlocal;

  copymode = 1;

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixTempRescale_end_of_step>(0,nglocal),*this);

  copymode = 0;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixTempRescaleKokkos::operator()(TagFixTempRescale_end_of_step, const int &icell) const {

  // loop over grid cells with more than 1 particle

  int ip,ispecies;
  double mass;
  double count,totmass,mvx,mvy,mvz,mvsq;
  double invtotmass,vscale;
  double vxcom,vycom,vzcom;
  double t_current;
  double *v;

  const int np = d_cellcount[icell];
  if (np <= 1) return;

  count = 0.0;
  totmass = 0.0;
  mvx = mvy = mvz = 0.0;
  mvsq = 0.0;

  // 1st pass: loop over particles in cell
  // 6 tallies per particle: N, Mass, mVx, mVy, mVz, mV^2

  for (int n = 0; n < np; n++) {
    const int ip = d_plist(icell,n);

    ispecies = d_particles[ip].ispecies;
    mass = d_species[ispecies].mass;
    v = d_particles[ip].v;

    count += 1.0;
    totmass += mass;
    mvx += mass*v[0];
    mvy += mass*v[1];
    mvz += mass*v[2];
    mvsq += mass * (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  }

  // COM velocity of particles in grid cell

  invtotmass = 1.0/totmass;
  vxcom = mvx * invtotmass;
  vycom = mvy * invtotmass;
  vzcom = mvz * invtotmass;

  // t_current = thermal T of particles in grid cell

  t_current = mvsq - (mvx*mvx + mvy*mvy + mvz*mvz)*invtotmass;
  t_current *= tprefactor/count;

  // vscale = scale factor for thermal velocity components

  vscale = sqrt(t_target/t_current);

  // 2nd pass: loop over particles in cell
  // rescale thermal velocity components

  for (int n = 0; n < np; n++) {
    const int ip = d_plist(icell,n);
    ispecies = d_particles[ip].ispecies;
    mass = d_species[ispecies].mass;
    v = d_particles[ip].v;

    v[0] = vscale*(v[0]-vxcom) + vxcom;
    v[1] = vscale*(v[1]-vycom) + vycom;
    v[2] = vscale*(v[2]-vzcom) + vzcom;
  }
}

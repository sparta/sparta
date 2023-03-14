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

/* ----------------------------------------------------------------------
   current thermal temperature is calculated on a per-cell basis
---------------------------------------------------------------------- */

void FixTempRescaleKokkos::end_of_step_no_average(double t_target_in)
{
  t_target = t_target_in;

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

  int nglocal = grid->nlocal;

  copymode = 1;

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixTempRescale_end_of_step_no_average>(0,nglocal),*this);

  copymode = 0;

  d_plist = decltype(d_plist)();
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixTempRescaleKokkos::operator()(TagFixTempRescale_end_of_step_no_average, const int &icell) const {

  // loop over grid cells with more than 1 particle

  double totmass,mvx,mvy,mvz,mvsq;
  double *v;

  const int count = d_cellcount[icell];
  if (count <= 1) return;

  totmass = 0.0;
  mvx = mvy = mvz = 0.0;
  mvsq = 0.0;

  // 1st pass: loop over particles in cell
  // 5 tallies per particle: Mass, mVx, mVy, mVz, mV^2

  for (int n = 0; n < count; n++) {
    const int ip = d_plist(icell,n);

    const int ispecies = d_particles[ip].ispecies;
    const double mass = d_species[ispecies].mass;
    v = d_particles[ip].v;

    totmass += mass;
    mvx += mass*v[0];
    mvy += mass*v[1];
    mvz += mass*v[2];
    mvsq += mass * (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  }

  // COM velocity of particles in grid cell

  const double invtotmass = 1.0/totmass;
  const double vxcom = mvx * invtotmass;
  const double vycom = mvy * invtotmass;
  const double vzcom = mvz * invtotmass;

  // t_current = thermal T of particles in grid cell

  double t_current = mvsq - (mvx*mvx + mvy*mvy + mvz*mvz)*invtotmass;
  t_current *= tprefactor/count;

  // vscale = scale factor for thermal velocity components

  const double vscale = sqrt(t_target/t_current);

  // 2nd pass: loop over particles in cell
  // rescale thermal velocity components

  for (int n = 0; n < count; n++) {
    const int ip = d_plist(icell,n);
    const int ispecies = d_particles[ip].ispecies;
    const double mass = d_species[ispecies].mass;
    v = d_particles[ip].v;

    v[0] = vscale*(v[0]-vxcom) + vxcom;
    v[1] = vscale*(v[1]-vycom) + vycom;
    v[2] = vscale*(v[2]-vzcom) + vzcom;
  }
}

/* ----------------------------------------------------------------------
   current thermal temperature is averaged over all per-cell temperatures
---------------------------------------------------------------------- */

void FixTempRescaleKokkos::end_of_step_average(double t_target_in)
{
  t_target = t_target_in;

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
  d_cells = grid_kk->k_cells.d_view;
  grid_kk->sync(Device,CELL_MASK);

  int nglocal = grid->nlocal;

  // resize d_vcom if needed

  if (nglocal > maxgrid) {
    d_vcom = DAT::t_float_1d_3();
    maxgrid = nglocal + grid->nghost;
    d_vcom = DAT::t_float_1d_3("temp/rescale:d_vcom",maxgrid);
  }

  // loop over grid cells to compute thermal T of each
  // only unsplit and sub cells, skip split cells

  REDUCE current_mine;

  copymode = 1;

  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixTempRescale_end_of_step_average1>(0,nglocal),*this,current_mine);

  double t_current_mine = current_mine.t;
  bigint n_current_mine = current_mine.n;

  // t_current = average of cellwise thermal T across all cells
  // n_current = total # of cells contributing to t_current

  double t_current;
  MPI_Allreduce(&t_current_mine,&t_current,1,MPI_DOUBLE,MPI_SUM,world);

  bigint n_current;
  MPI_Allreduce(&n_current_mine,&n_current,1,MPI_SPARTA_BIGINT,MPI_SUM,world);

  // t_current = cellwise averaged thermal T
  // scale all particles in all cells by vscale

  t_current /= n_current;
  vscale = sqrt(t_target/t_current);

  // loop over grid cells to rescale velocity of particles in each
  // single-particle cells are also rescaled, their d_vcom = 0.0

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixTempRescale_end_of_step_average2>(0,nglocal),*this);

  copymode = 0;
  d_plist = decltype(d_plist)();
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixTempRescaleKokkos::operator()(TagFixTempRescale_end_of_step_average1, const int &icell, REDUCE &current_mine) const {

  double totmass,mvx,mvy,mvz,mvsq;
  double t_one;
  double *v;

  if (d_cells[icell].nsplit > 1) return;

  const int count = d_cellcount[icell];
  totmass = 0.0;
  mvx = mvy = mvz = 0.0;
  mvsq = 0.0;

  // loop over particles in cell
  // 5 tallies per particle: Mass, mVx, mVy, mVz, mV^2

  for (int n = 0; n < count; n++) {
    const int ip = d_plist(icell,n);
    const int ispecies = d_particles[ip].ispecies;
    const double mass = d_species[ispecies].mass;
    v = d_particles[ip].v;

    totmass += mass;
    mvx += mass*v[0];
    mvy += mass*v[1];
    mvz += mass*v[2];
    mvsq += mass * (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  }

  // t_one = thermal T of particles in the cell
  // d_vcom = COM velocity of particle in the cell
  // if cell has <= 1 particle: set t_one = t_target, d_vcom = zero
  // likewise if t_one = 0.0: set t_one = t_target, d_vcom = zero
  //   corner case when all particles have same velocity

  if (count > 1.0) {
    const double invtotmass = 1.0/totmass;
    t_one = mvsq - (mvx*mvx + mvy*mvy + mvz*mvz)*invtotmass;
    t_one *= tprefactor/count;
    d_vcom(icell,0) = mvx * invtotmass;
    d_vcom(icell,1) = mvy * invtotmass;
    d_vcom(icell,2) = mvz * invtotmass;

  } else {
    t_one = t_target;
    d_vcom(icell,0) = d_vcom(icell,1) = d_vcom(icell,2) = 0.0;
  }

  if (t_one == 0.0) {
    t_one = t_target;
    d_vcom(icell,0) = d_vcom(icell,1) = d_vcom(icell,2) = 0.0;
  }

  // accumulate thermal T over my cells

  current_mine.t += t_one;
  current_mine.n++;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixTempRescaleKokkos::operator()(TagFixTempRescale_end_of_step_average2, const int &icell) const {

  const int count = d_cellcount[icell];
  if (count == 0) return;

  // loop over particles in cell
  // rescale thermal velocity components

  for (int n = 0; n < count; n++) {
    const int ip = d_plist(icell,n);
    const int ispecies = d_particles[ip].ispecies;
    const double mass = d_species[ispecies].mass;
    double* v = d_particles[ip].v;

    v[0] = vscale*(v[0]-d_vcom(icell,0)) + d_vcom(icell,0);
    v[1] = vscale*(v[1]-d_vcom(icell,1)) + d_vcom(icell,1);
    v[2] = vscale*(v[2]-d_vcom(icell,2)) + d_vcom(icell,2);
  }
}

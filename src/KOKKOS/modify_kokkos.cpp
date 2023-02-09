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

#include "stdio.h"
#include "string.h"
#include "modify_kokkos.h"
#include "domain.h"
#include "update.h"
#include "compute.h"
#include "fix.h"
#include "style_compute.h"
#include "style_fix.h"
#include "memory_kokkos.h"
#include "error.h"
#include "particle_kokkos.h"
#include "grid_kokkos.h"
#include "kokkos.h"

using namespace SPARTA_NS;

#define DELTA 4

// mask settings - same as in fix.cpp

#define START_OF_STEP  1
#define END_OF_STEP    2

/* ---------------------------------------------------------------------- */

ModifyKokkos::ModifyKokkos(SPARTA *sparta) : Modify(sparta)
{
  particle_kk = (ParticleKokkos*) particle;
  grid_kk = (GridKokkos*) grid;
}

/* ---------------------------------------------------------------------- */

ModifyKokkos::~ModifyKokkos()
{

}

/* ----------------------------------------------------------------------
   start-of-timestep call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::start_of_step()
{
  for (int i = 0; i < n_start_of_step; i++) {
    int j = list_start_of_step[i];
    particle_kk->sync(fix[j]->execution_space,fix[j]->datamask_read);
    int prev_auto_sync = sparta->kokkos->auto_sync;
    if (!fix[j]->kokkos_flag) sparta->kokkos->auto_sync = 1;

    fix[list_start_of_step[i]]->start_of_step();

    sparta->kokkos->auto_sync = prev_auto_sync;
    particle_kk->modify(fix[j]->execution_space,fix[j]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   end-of-timestep call, only for relevant fixes
   only call fix->end_of_step() on timesteps that are multiples of nevery
------------------------------------------------------------------------- */

void ModifyKokkos::end_of_step()
{
  for (int i = 0; i < n_end_of_step; i++)
    if (update->ntimestep % end_of_step_every[i] == 0) {
      int j = list_end_of_step[i];
      particle_kk->sync(fix[j]->execution_space,fix[j]->datamask_read);
      int prev_auto_sync = sparta->kokkos->auto_sync;
      if (!fix[j]->kokkos_flag) sparta->kokkos->auto_sync = 1;

      fix[list_end_of_step[i]]->end_of_step();

      sparta->kokkos->auto_sync = prev_auto_sync;
      particle_kk->modify(fix[j]->execution_space,fix[j]->datamask_modify);
    }
}

/* ----------------------------------------------------------------------
   pack_grid_one call, only for relevant fixes
   invoked by load balancer when grid cells migrate
------------------------------------------------------------------------- */

int ModifyKokkos::pack_grid_one(int icell, char *buf, int memflag)
{
  char *ptr = buf;
  for (int i = 0; i < n_pergrid; i++) {
    int j = list_pergrid[i];
    particle_kk->sync(fix[j]->execution_space,fix[j]->datamask_read);
    int prev_auto_sync = sparta->kokkos->auto_sync;
    if (!fix[j]->kokkos_flag) sparta->kokkos->auto_sync = 1;

    ptr += fix[list_pergrid[i]]->pack_grid_one(icell,ptr,memflag);

    sparta->kokkos->auto_sync = prev_auto_sync;
    particle_kk->modify(fix[j]->execution_space,fix[j]->datamask_modify);
  }
  return ptr-buf;
}

/* ----------------------------------------------------------------------
   unpack_grid_one call, only for relevant fixes
   invoked by load balancer when grid cells migrate
------------------------------------------------------------------------- */

int ModifyKokkos::unpack_grid_one(int icell, char *buf)
{
  char *ptr = buf;
  for (int i = 0; i < n_pergrid; i++) {
    int j = list_pergrid[i];
    particle_kk->sync(fix[j]->execution_space,fix[j]->datamask_read);
    int prev_auto_sync = sparta->kokkos->auto_sync;
    if (!fix[j]->kokkos_flag) sparta->kokkos->auto_sync = 1;

    ptr += fix[list_pergrid[i]]->unpack_grid_one(icell,ptr);

    sparta->kokkos->auto_sync = prev_auto_sync;
    particle_kk->modify(fix[j]->execution_space,fix[j]->datamask_modify);
  }
  return ptr-buf;
}

/* ----------------------------------------------------------------------
   copy_grid call, only for relevant fixes
   invoked when a grod cell is removed
------------------------------------------------------------------------- */

void ModifyKokkos::copy_grid_one(int icell, int jcell)
{
  for (int i = 0; i < n_pergrid; i++) {
    int j = list_pergrid[i];
    particle_kk->sync(fix[j]->execution_space,fix[j]->datamask_read);
    int prev_auto_sync = sparta->kokkos->auto_sync;
    if (!fix[j]->kokkos_flag) sparta->kokkos->auto_sync = 1;

    fix[j]->copy_grid_one(icell,jcell);

    sparta->kokkos->auto_sync = prev_auto_sync;
    particle_kk->modify(fix[j]->execution_space,fix[j]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   add_grid_one call, only for relevant fixes
   invoked by adapt_grid and fix adapt when new child cells are created
------------------------------------------------------------------------- */

void ModifyKokkos::add_grid_one()
{
  for (int i = 0; i < n_pergrid; i++) {
    int j = list_pergrid[i];
    particle_kk->sync(fix[j]->execution_space,fix[j]->datamask_read);
    int prev_auto_sync = sparta->kokkos->auto_sync;
    if (!fix[j]->kokkos_flag) sparta->kokkos->auto_sync = 1;

    fix[j]->add_grid_one();

    sparta->kokkos->auto_sync = prev_auto_sync;
    particle_kk->modify(fix[j]->execution_space,fix[j]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   reset_grid call, only for relevant fixes
   invoked after all grid cell removals
------------------------------------------------------------------------- */

void ModifyKokkos::reset_grid_count(int nlocal)
{
  for (int i = 0; i < n_pergrid; i++) {
    int j = list_pergrid[i];
    particle_kk->sync(fix[j]->execution_space,fix[j]->datamask_read);
    int prev_auto_sync = sparta->kokkos->auto_sync;
    if (!fix[j]->kokkos_flag) sparta->kokkos->auto_sync = 1;

    fix[j]->reset_grid_count(nlocal);

    sparta->kokkos->auto_sync = prev_auto_sync;
    particle_kk->modify(fix[j]->execution_space,fix[j]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   grid_changed call, only for relevant fixes
   invoked after per-processor list of grid cells has changed
------------------------------------------------------------------------- */

void ModifyKokkos::grid_changed()
{
  for (int i = 0; i < n_pergrid; i++) {
    int j = list_pergrid[i];
    particle_kk->sync(fix[j]->execution_space,fix[j]->datamask_read);
    int prev_auto_sync = sparta->kokkos->auto_sync;
    if (!fix[j]->kokkos_flag) sparta->kokkos->auto_sync = 1;

    fix[j]->grid_changed();

    sparta->kokkos->auto_sync = prev_auto_sync;
    particle_kk->modify(fix[j]->execution_space,fix[j]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   invoke update_custom() method, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::update_custom(int index, double temp_thermal,
                                 double temp_rot, double temp_vib, double *vstream)
{
  for (int i = 0; i < n_update_custom; i++) {
    int j = list_update_custom[i];
    particle_kk->sync(fix[j]->execution_space,fix[j]->datamask_read);
    int prev_auto_sync = sparta->kokkos->auto_sync;
    if (!fix[j]->kokkos_flag) sparta->kokkos->auto_sync = 1;

    fix[list_update_custom[i]]->update_custom(index,temp_thermal,temp_rot,
                                            temp_vib,vstream);

    sparta->kokkos->auto_sync = prev_auto_sync;
    particle_kk->modify(fix[j]->execution_space,fix[j]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   invoke gas_react() method, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::gas_react(int index)
{
  for (int i = 0; i < n_gas_react; i++) {
    int j = list_gas_react[i];
    particle_kk->sync(fix[j]->execution_space,fix[j]->datamask_read);
    int prev_auto_sync = sparta->kokkos->auto_sync;
    if (!fix[j]->kokkos_flag) sparta->kokkos->auto_sync = 1;

    fix[list_gas_react[i]]->gas_react(index);

    sparta->kokkos->auto_sync = prev_auto_sync;
    particle_kk->modify(fix[j]->execution_space,fix[j]->datamask_modify);
  }
}

/* ----------------------------------------------------------------------
   invoke surf_react() method, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::surf_react(Particle::OnePart *iorig, int &i, int &)
{
  for (int m = 0; m < n_surf_react; m++) {
    int j = list_surf_react[m];
    particle_kk->sync(fix[j]->execution_space,fix[j]->datamask_read);
    int prev_auto_sync = sparta->kokkos->auto_sync;
    if (!fix[j]->kokkos_flag) sparta->kokkos->auto_sync = 1;

    fix[list_surf_react[m]]->surf_react(iorig,i,j);

    sparta->kokkos->auto_sync = prev_auto_sync;
    particle_kk->modify(fix[j]->execution_space,fix[j]->datamask_modify);
  }
}

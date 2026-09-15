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

#include "string.h"
#include "fix_rigid_kokkos.h"
#include "update.h"
#include "grid_kokkos.h"
#include "particle_kokkos.h"
#include "surf_kokkos.h"
#include "compute_surf.h"
#include "error.h"
#include "sparta_masks.h"

using namespace SPARTA_NS;

/* ----------------------------------------------------------------------
   KOKKOS version of fix rigid
   the body time integration, force/torque reduction, swept collision
     lists, grid re-map and deletion of particles inside the body run on
     the host, exactly as in FixRigid; the moving-surf collision tests
     run in the KOKKOS particle mover (UpdateKokkos::move)
   this class brackets the host work with the device<->host transfers
     it needs, and rebuilds the device per-cell surf graphs whenever the
     host work changed the per-cell surf lists
------------------------------------------------------------------------- */

FixRigidKokkos::FixRigidKokkos(SPARTA *sparta, int narg, char **arg) :
  FixRigid(sparta, narg, arg)
{
  kokkos_flag = 0; // need auto sync
  execution_space = Host;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
  kokkosable = 1;
}

/* ---------------------------------------------------------------------- */

void FixRigidKokkos::init()
{
  ((GridKokkos*) grid)->sync(Host,ALL_MASK);
  ((SurfKokkos*) surf)->sync(Host,ALL_MASK);

  FixRigid::init();

  // force/torque tallies must come from the KOKKOS mover, which only
  //   tallies into the KOKKOS variant of compute surf

  if (!csurf->kokkos_flag)
    error->all(FLERR,"Fix rigid/kk requires compute surf/kk");
}

/* ----------------------------------------------------------------------
   1 if this is the last-defined rigid fix, which performs the all-body
     grid operations (swept assignment, re-map) each step
------------------------------------------------------------------------- */

int FixRigidKokkos::last_body()
{
  return (update->fixrigidlist[update->nfixrigid-1] == this);
}

/* ----------------------------------------------------------------------
   bring grid, particles and surfs to the host before host-side work
------------------------------------------------------------------------- */

void FixRigidKokkos::host_begin()
{
  ((GridKokkos*) grid)->sync(Host,ALL_MASK);
  ((ParticleKokkos*) particle)->sync(Host,PARTICLE_MASK|CUSTOM_MASK);
  ((SurfKokkos*) surf)->sync(Host,ALL_MASK);
}

/* ----------------------------------------------------------------------
   after host-side work: the host copies are authoritative
   if the grid was rebuilt (full re-map, flagged via Grid::changed),
     GridKokkos::resync_after_host_change() re-establishes the device
     grid, hash and per-cell surf graphs; otherwise only the per-cell
     surf lists changed (swept lists, incremental re-cut), so rewrap
     the surf graphs
------------------------------------------------------------------------- */

void FixRigidKokkos::host_end()
{
  GridKokkos *grid_kk = (GridKokkos*) grid;
  ParticleKokkos *particle_kk = (ParticleKokkos*) particle;
  SurfKokkos *surf_kk = (SurfKokkos*) surf;

  // particles: compress_rebalance() moves particles and their custom data
  // surfs: the body geometry was regenerated on the host
  // per-cell surf lists: only rewrap the device graphs if a list changed

  grid_kk->modify(Host,ALL_MASK);
  particle_kk->modify(Host,PARTICLE_MASK|CUSTOM_MASK);
  surf_kk->modify(Host,ALL_MASK);
  particle_kk->sorted_kk = 0;

  if (grid->changed) grid_kk->resync_after_host_change();
  else if (any_lists_changed()) grid_kk->wrap_kokkos_graphs();
  clear_lists_changed();
}

/* ---------------------------------------------------------------------- */

void FixRigidKokkos::setup()
{
  host_begin();
  FixRigid::setup();
  host_end();
}

/* ----------------------------------------------------------------------
   integrate the body and install the swept collision lists (host),
     then rebuild the device surf graphs the mover reads
   only the last-defined fix changes the per-cell lists, so only it
     rewraps; surfs and cells are not otherwise changed here
------------------------------------------------------------------------- */

void FixRigidKokkos::start_of_step()
{
  GridKokkos *grid_kk = (GridKokkos*) grid;
  SurfKokkos *surf_kk = (SurfKokkos*) surf;

  grid_kk->sync(Host,ALL_MASK);
  surf_kk->sync(Host,ALL_MASK);

  FixRigid::start_of_step();

  if (last_body()) {
    grid_kk->modify(Host,CELL_MASK);
    if (any_lists_changed()) grid_kk->wrap_kokkos_graphs();
    clear_lists_changed();
  }
}

/* ----------------------------------------------------------------------
   restore swept lists, reduce force/torque, move the body, re-map the
     grid and delete particles inside the body (host), then re-establish
     the device state
   the first-defined fix restores the swept lists and the last-defined
     fix re-maps, so the device is re-established after the last one
------------------------------------------------------------------------- */

void FixRigidKokkos::end_of_step()
{
  host_begin();
  FixRigid::end_of_step();
  if (last_body()) host_end();
  else {
    ((SurfKokkos*) surf)->modify(Host,ALL_MASK);
    ((GridKokkos*) grid)->modify(Host,CELL_MASK);
  }
}

/* ----------------------------------------------------------------------
   grid cells were rebuilt, adapted, or migrated: FixRigid re-establishes
     its local body-surf copies on the host, which may append surfs and
     re-index ghost cell lists, so flag the host surfs and cells as
     modified; the caller (grid rebuild, fix balance/kk, fix adapt/kk)
     rewraps the device surf graphs afterward
------------------------------------------------------------------------- */

void FixRigidKokkos::grid_changed()
{
  FixRigid::grid_changed();
  if (surf->distributed) ((SurfKokkos*) surf)->modify(Host,ALL_MASK);
  ((GridKokkos*) grid)->modify(Host,CELL_MASK);
}

/* ----------------------------------------------------------------------
   1 if any rigid fix changed a per-cell surf list on the host since the
     device graphs were last rewrapped: swept lists installed or
     restored, or cells re-cut incrementally
------------------------------------------------------------------------- */

int FixRigidKokkos::any_lists_changed()
{
  for (int m = 0; m < update->nfixrigid; m++)
    if (update->fixrigidlist[m]->listschanged) return 1;
  return 0;
}

void FixRigidKokkos::clear_lists_changed()
{
  for (int m = 0; m < update->nfixrigid; m++)
    update->fixrigidlist[m]->listschanged = 0;
}

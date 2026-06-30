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
#include "compute_isurf_grid_kokkos.h"
#include "particle_kokkos.h"
#include "mixture.h"
#include "surf_kokkos.h"
#include "grid.h"
#include "update.h"
#include "memory_kokkos.h"
#include "error.h"
#include "sparta_masks.h"
#include "kokkos.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

ComputeISurfGridKokkos::ComputeISurfGridKokkos(SPARTA *sparta, int narg, char **arg) :
  ComputeISurfGrid(sparta, narg, arg)
{
  kokkos_flag = 1;

  // hash is allocated/used only on the host; not needed for device tally

  d_which = DAT::t_int_1d("isurf/grid:which",nvalue);
}

ComputeISurfGridKokkos::ComputeISurfGridKokkos(SPARTA *sparta) :
  ComputeISurfGrid(sparta)
{
  copy = 1;
  uncopy = 0;
}

/* ---------------------------------------------------------------------- */

ComputeISurfGridKokkos::~ComputeISurfGridKokkos()
{
  if (copy) return;

  memoryKK->destroy_kokkos(k_tally2surf,tally2surf);
  memoryKK->destroy_kokkos(k_array_surf_tally,array_surf_tally);
}

/* ---------------------------------------------------------------------- */

void ComputeISurfGridKokkos::init()
{
  ComputeISurfGrid::init();

  auto h_which = Kokkos::create_mirror_view(d_which);
  for (int n=0; n<nvalue; n++)
    h_which(n) = which[n];
  Kokkos::deep_copy(d_which,h_which);
}

/* ----------------------------------------------------------------------
   set normflux for all surfs I store, mirror to device
------------------------------------------------------------------------- */

void ComputeISurfGridKokkos::init_normflux()
{
  ComputeISurfGrid::init_normflux();

  int nsurf = surf->nlocal + surf->nghost;

  d_normflux = DAT::t_float_1d("isurf/grid:normflux",nsurf);
  auto h_normflux = Kokkos::create_mirror_view(d_normflux);
  for (int n=0; n<nsurf; n++)
    h_normflux(n) = normflux[n];
  Kokkos::deep_copy(d_normflux,h_normflux);

  // Cannot realloc inside a Kokkos parallel region, so size tally2surf as nsurf

  memoryKK->grow_kokkos(k_tally2surf,tally2surf,nsurf,"isurf/grid:tally2surf");
  d_tally2surf = k_tally2surf.view_device();
  d_surf2tally = DAT::t_int_1d("isurf/grid:surf2tally",nsurf);
  Kokkos::deep_copy(d_surf2tally,-1);

  memoryKK->grow_kokkos(k_array_surf_tally,array_surf_tally,nsurf,ntotal,"isurf/grid:array_surf_tally");
  d_array_surf_tally = k_array_surf_tally.view_device();
}

/* ---------------------------------------------------------------------- */

void ComputeISurfGridKokkos::clear()
{
  // reset all set surf2tally values to -1
  // called by Update at beginning of timesteps surf tallying is done

  Kokkos::deep_copy(d_array_surf_tally,0);
  Kokkos::deep_copy(d_surf2tally,-1);

  ntally = 0;
  combined = 0;
}

/* ---------------------------------------------------------------------- */

void ComputeISurfGridKokkos::pre_surf_tally()
{
  mvv2e = update->mvv2e;

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,SPECIES_MASK);
  d_species = particle_kk->k_species.view_device();
  d_s2g = particle_kk->k_species2group.view_device();

  SurfKokkos* surf_kk = (SurfKokkos*) surf;
  surf_kk->sync(Device,ALL_MASK);
  d_lines = surf_kk->k_lines.view_device();
  d_tris = surf_kk->k_tris.view_device();

  need_dup = sparta->kokkos->need_dup<DeviceType>();
  if (need_dup)
    dup_array_surf_tally = Kokkos::Experimental::create_scatter_view<typename Kokkos::Experimental::ScatterSum, typename Kokkos::Experimental::ScatterDuplicated>(d_array_surf_tally);
  else
    ndup_array_surf_tally = Kokkos::Experimental::create_scatter_view<typename Kokkos::Experimental::ScatterSum, typename Kokkos::Experimental::ScatterNonDuplicated>(d_array_surf_tally);
}

/* ---------------------------------------------------------------------- */

void ComputeISurfGridKokkos::post_surf_tally()
{
  if (need_dup) {
    Kokkos::Experimental::contribute(d_array_surf_tally, dup_array_surf_tally);
    dup_array_surf_tally = {}; // free duplicated memory
  }

  k_tally2surf.modify_device();
  k_array_surf_tally.modify_device();
}

/* ----------------------------------------------------------------------
   sync device tallies to host and compress to dense list (ntally tallies)
   matches ComputeSurfKokkos::tallyinfo();
   host fix ave/grid (PERGRIDSURF) consumes array_surf_tally + tally2surf
------------------------------------------------------------------------- */

int ComputeISurfGridKokkos::tallyinfo(surfint *&ptr)
{
  k_tally2surf.sync_host();
  ptr = tally2surf;

  k_array_surf_tally.sync_host();
  auto h_surf2tally = Kokkos::create_mirror_view(d_surf2tally);
  Kokkos::deep_copy(h_surf2tally,d_surf2tally);

  // compress array_surf_tally

  int nsurf = surf->nlocal + surf->nghost;
  int istart = 0;
  int iend = nsurf-1;

  while (1) {
    while (h_surf2tally[istart] != -1 && istart < nsurf-2) istart++;
    while (h_surf2tally[iend] == -1 && iend > 0) iend--;
    if (istart >= iend) {
      ntally = istart;
      break;
    }
    for (int k = 0; k < ntotal; k++) {
      array_surf_tally[istart][k] = array_surf_tally[iend][k];
    }
    h_surf2tally[istart] = h_surf2tally[iend];
    h_surf2tally[iend] = -1;
    tally2surf[istart] = tally2surf[iend];
  }

  return ntally;
}

/* ---------------------------------------------------------------------- */

void ComputeISurfGridKokkos::grow_tally()
{
  // Cannot realloc inside a Kokkos parallel region, so size as nsurf

  int nsurf = surf->nlocal + surf->nghost;

  memoryKK->grow_kokkos(k_tally2surf,tally2surf,nsurf,"isurf/grid:tally2surf");
  d_tally2surf = k_tally2surf.view_device();

  memoryKK->grow_kokkos(k_array_surf_tally,array_surf_tally,nsurf,ntotal,"isurf/grid:array_surf_tally");
  d_array_surf_tally = k_array_surf_tally.view_device();
}

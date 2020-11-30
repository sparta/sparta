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

#include "string.h"
#include "compute_surf_kokkos.h"
#include "particle_kokkos.h"
#include "mixture.h"
#include "surf_kokkos.h"
#include "grid.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "memory_kokkos.h"
#include "error.h"
#include "sparta_masks.h"
#include "kokkos.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

ComputeSurfKokkos::ComputeSurfKokkos(SPARTA *sparta, int narg, char **arg) :
  ComputeSurf(sparta, narg, arg)
{
  kokkos_flag = 1;
  d_which = DAT::t_int_1d("surf:which",nvalue);

  d_ntally = DAT::t_int_scalar("surf:ntally");
}

ComputeSurfKokkos::ComputeSurfKokkos(SPARTA *sparta) :
  ComputeSurf(sparta)
{
  hash = NULL;
  which = NULL;
  array_surf_tally = NULL;
  tally2surf = NULL;
  array_surf = NULL;
  vector_surf = NULL;
  normflux = NULL;
  id = NULL;
  style = NULL;
  tlist = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeSurfKokkos::~ComputeSurfKokkos()
{
  if (copy || copymode) return;

  memoryKK->destroy_kokkos(k_tally2surf,tally2surf);
  memoryKK->destroy_kokkos(k_array_surf_tally,array_surf_tally);
}

/* ---------------------------------------------------------------------- */

void ComputeSurfKokkos::init()
{
  ComputeSurf::init();

  auto h_which = Kokkos::create_mirror_view(d_which);
  for (int n=0; n<nvalue; n++)
    h_which(n) = which[n];
  Kokkos::deep_copy(d_which,h_which);
}

/* ----------------------------------------------------------------------
   set normflux for all surfs I store
   all: just nlocal
   distributed: nlocal + nghost
   called by init before each run (in case dt or fnum has changed)
   called whenever grid changes
------------------------------------------------------------------------- */

void ComputeSurfKokkos::init_normflux()
{
  ComputeSurf::init_normflux();

  int nsurf = surf->nlocal + surf->nghost;

  d_normflux = DAT::t_float_1d("surf:normflux",nsurf);
  auto h_normflux = Kokkos::create_mirror_view(d_normflux);
  for (int n=0; n<nsurf; n++)
    h_normflux(n) = normflux[n];
  Kokkos::deep_copy(d_normflux,h_normflux);

  // Cannot realloc inside a Kokkos parallel region, so size tally2surf as nsurf
  memoryKK->grow_kokkos(k_tally2surf,tally2surf,nsurf,"surf:tally2surf");
  d_tally2surf = k_tally2surf.d_view;
  d_surf2tally = DAT::t_int_1d("surf:surf2tally",nsurf);

  memoryKK->grow_kokkos(k_array_surf_tally,array_surf_tally,nsurf,ntotal,"surf:array_surf_tally");
  d_array_surf_tally = k_array_surf_tally.d_view;
}

/* ---------------------------------------------------------------------- */

void ComputeSurfKokkos::clear()
{
  // reset all set surf2tally values to -1
  // called by Update at beginning of timesteps surf tallying is done

  ntally = 0;
  combined = 0;
  Kokkos::deep_copy(d_ntally,0);
  Kokkos::deep_copy(d_array_surf_tally,0);

  Kokkos::deep_copy(d_surf2tally,-1);
}

/* ---------------------------------------------------------------------- */

void ComputeSurfKokkos::pre_surf_tally()
{
  mvv2e = update->mvv2e;

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,SPECIES_MASK);
  d_species = particle_kk->k_species.d_view;
  d_s2g = particle_kk->k_species2group.d_view;

  SurfKokkos* surf_kk = (SurfKokkos*) surf;
  surf_kk->sync(Device,ALL_MASK);
  d_lines = surf_kk->k_lines.d_view;
  d_tris = surf_kk->k_tris.d_view;

  need_dup = sparta->kokkos->need_dup<DeviceType>();
  if (need_dup)
    dup_array_surf_tally = Kokkos::Experimental::create_scatter_view<typename Kokkos::Experimental::ScatterSum, typename Kokkos::Experimental::ScatterDuplicated>(d_array_surf_tally);
  else
    ndup_array_surf_tally = Kokkos::Experimental::create_scatter_view<typename Kokkos::Experimental::ScatterSum, typename Kokkos::Experimental::ScatterNonDuplicated>(d_array_surf_tally);
}

/* ---------------------------------------------------------------------- */

void ComputeSurfKokkos::post_surf_tally()
{
  if (need_dup) {
    Kokkos::Experimental::contribute(d_array_surf_tally, dup_array_surf_tally);
    dup_array_surf_tally = decltype(dup_array_surf_tally)(); // free duplicated memory
  }

  k_tally2surf.modify_device();
  k_array_surf_tally.modify_device();
}

/* ----------------------------------------------------------------------
   return ptr to norm vector used by column N
------------------------------------------------------------------------- */

int ComputeSurfKokkos::tallyinfo(surfint *&ptr)
{
  k_tally2surf.sync_host();
  ptr = tally2surf;

  k_array_surf_tally.sync_host();

  auto h_ntally = Kokkos::create_mirror_view(d_ntally);
  Kokkos::deep_copy(h_ntally,d_ntally);
  return h_ntally();
}

/* ---------------------------------------------------------------------- */

void ComputeSurfKokkos::grow_tally()
{
  // Cannot realloc inside a Kokkos parallel region, so size tally2surf the
  //  same as surf2tally

  int nsurf = surf->nlocal + surf->nghost;

  memoryKK->grow_kokkos(k_tally2surf,tally2surf,nsurf,"surf:tally2surf");
  d_tally2surf = k_tally2surf.d_view;

  memoryKK->grow_kokkos(k_array_surf_tally,array_surf_tally,nsurf,ntotal,"surf:array_surf_tally");
  d_array_surf_tally = k_array_surf_tally.d_view;
}


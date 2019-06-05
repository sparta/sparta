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
  which = NULL;
  surf2tally = NULL;
  tally2surf = NULL;
  array_surf_tally = NULL;
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

  d_surf2tally = DAT::t_int_1d("surf:surf2tally",nsurf);
  Kokkos::deep_copy(d_surf2tally,-1);

  d_normflux = DAT::t_float_1d("surf:normflux",nsurf);
  auto h_normflux = Kokkos::create_mirror_view(d_normflux);
  for (int n=0; n<nsurf; n++)
    h_normflux(n) = normflux[n];
  Kokkos::deep_copy(d_normflux,h_normflux);

  // Cannot realloc inside a Kokkos parallel region, so size tally2surf the 
  //  same as surf2tally 

  memoryKK->destroy_kokkos(k_tally2surf,tally2surf);
  memoryKK->create_kokkos(k_tally2surf,tally2surf,nsurf,"surf:tally2surf");
  d_tally2surf = k_tally2surf.d_view;

  memoryKK->destroy_kokkos(k_array_surf_tally,array_surf_tally);
  memoryKK->create_kokkos(k_array_surf_tally,array_surf_tally,nsurf,ntotal,"surf:array_surf_tally");
  d_array_surf_tally = k_array_surf_tally.d_view;
}

/* ---------------------------------------------------------------------- */

void ComputeSurfKokkos::clear()
{
  // reset all set surf2tally values to -1
  // called by Update at beginning of timesteps surf tallying is done

  auto h_ntally = Kokkos::create_mirror_view(d_ntally);
  Kokkos::deep_copy(h_ntally,d_ntally);
  ntally = h_ntally();

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeSurf_clear>(0,ntally),*this);
  DeviceType::fence();
  copymode = 0;

  ntally = 0;
  Kokkos::deep_copy(d_ntally,0);
  Kokkos::deep_copy(d_array_surf_tally,0);
}

KOKKOS_INLINE_FUNCTION
void ComputeSurfKokkos::operator()(TagComputeSurf_clear, const int &i) const {
  d_surf2tally[d_tally2surf[i]] = -1;
}

/* ---------------------------------------------------------------------- */

void ComputeSurfKokkos::pre_surf_tally()
{
  mvv2e = update->mvv2e;

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,SPECIES_MASK);
  d_species = particle_kk->k_species.d_view;
  d_s2g = particle_kk->k_species2group.view<DeviceType>();

  SurfKokkos* surf_kk = (SurfKokkos*) surf;
  surf_kk->sync(Device,ALL_MASK);
  d_lines = surf_kk->k_lines.d_view;
  d_tris = surf_kk->k_tris.d_view;
}

/* ---------------------------------------------------------------------- */

void ComputeSurfKokkos::post_surf_tally()
{
  k_array_surf_tally.modify<DeviceType>();
  k_array_surf_tally.sync<SPAHostType>();
}

/* ----------------------------------------------------------------------
   return ptr to norm vector used by column N
------------------------------------------------------------------------- */

int ComputeSurfKokkos::tallyinfo(int *&locptr)
{
  k_tally2surf.modify<DeviceType>();
  k_tally2surf.sync<SPAHostType>();
  locptr = tally2surf;

  auto h_ntally = Kokkos::create_mirror_view(d_ntally);
  Kokkos::deep_copy(h_ntally,d_ntally);
  return h_ntally();
}

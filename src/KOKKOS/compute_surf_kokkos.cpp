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

#define VAL_1(X) X
#define VAL_2(X) VAL_1(X), VAL_1(X)

/* ---------------------------------------------------------------------- */

ComputeSurfKokkos::ComputeSurfKokkos(SPARTA *sparta, int narg, char **arg) :
  ComputeSurf(sparta, narg, arg)
#ifdef SPARTA_KOKKOS_FIXED_LISTS
  , sr_kk_global_copy{VAL_2(KKCopy<SurfReactGlobalKokkos>(sparta))}
  , sr_kk_prob_copy{VAL_2(KKCopy<SurfReactProbKokkos>(sparta))}
#endif
{
  kokkos_flag = 1;
  compressed = 0;
  d_which = DAT::t_int_1d("surf:which",nvalue);
}

ComputeSurfKokkos::ComputeSurfKokkos(SPARTA *sparta) :
  ComputeSurf(sparta)
#ifdef SPARTA_KOKKOS_FIXED_LISTS
  , sr_kk_global_copy{VAL_2(KKCopy<SurfReactGlobalKokkos>(sparta))}
  , sr_kk_prob_copy{VAL_2(KKCopy<SurfReactProbKokkos>(sparta))}
#endif
{
  copy = 1;
  compressed = 0;
}

/* ---------------------------------------------------------------------- */

ComputeSurfKokkos::~ComputeSurfKokkos()
{
  if (copy) return;

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
  d_tally2surf = k_tally2surf.view_device();
  d_surf2tally = DAT::t_int_1d("surf:surf2tally",nsurf);
  Kokkos::deep_copy(d_surf2tally,-1);

  memoryKK->grow_kokkos(k_array_surf_tally,array_surf_tally,nsurf,ntotal,"surf:array_surf_tally");
  d_array_surf_tally = k_array_surf_tally.view_device();
}

/* ---------------------------------------------------------------------- */

void ComputeSurfKokkos::clear()
{
  // tallyinfo() compresses the tally list in place on the host and claims the
  // host for it.  A new cycle starts here: the device copy below is about to be
  // zeroed and refilled by the tally kernels, so the two are deliberately
  // parted and neither owes the other a copy.  Without this the modify_device()
  // in post_surf_tally() would meet the outstanding host claim and abort.

  k_tally2surf.clear_sync_state();
  k_array_surf_tally.clear_sync_state();

  // reset all set surf2tally values to -1
  // called by Update at beginning of timesteps surf tallying is done

  Kokkos::deep_copy(d_array_surf_tally,0);

  Kokkos::deep_copy(d_surf2tally,-1);

  ntally = 0;
  combined = 0;
  compressed = 0;
}

/* ---------------------------------------------------------------------- */

void ComputeSurfKokkos::pre_surf_tally()
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

#ifdef SPARTA_KOKKOS_FIXED_LISTS
  if (surf->nsr > KOKKOS_MAX_TOT_SURF_REACT)
    error->all(FLERR,"Kokkos currently supports a limited number of surface reaction methods");
#else

  // the buffers must be sized before anything is blitted into them.  surf->nsr
  //   bounds the count of every individual style, so sizing both of them to it
  //   needs no counting pass, and the loop below still runs pre_react() in
  //   surf react list order

  sr_idx_resize(k_sr_type_list,d_sr_type_list,surf->nsr);
  sr_idx_resize(k_sr_map,d_sr_map,surf->nsr);
  sr_buf_resize<SurfReactGlobalKokkos>(k_sr_global,d_sr_global,surf->nsr);
  sr_buf_resize<SurfReactProbKokkos>(k_sr_prob,d_sr_prob,surf->nsr);
#endif

  if (surf->nsr > 0) {
    int nglob,nprob;
    nglob = nprob = 0;
    for (int n = 0; n < surf->nsr; n++) {
      if (!surf->sr[n]->kokkosable)
        error->all(FLERR,"Must use Kokkos-enabled surface reaction method with Kokkos");
      if (strcmp(surf->sr[n]->style,"global") == 0) {
#ifdef SPARTA_KOKKOS_FIXED_LISTS
        if (nglob >= KOKKOS_MAX_SURF_REACT_PER_TYPE)
          error->all(FLERR,"Kokkos currently supports two instances of each surface reaction method");
        sr_kk_global_copy[nglob].copy((SurfReactGlobalKokkos*)(surf->sr[n]));
#else
        sr_buf_blit(k_sr_global,nglob,(SurfReactGlobalKokkos*)(surf->sr[n]));
#endif
        KK_SR_H_GLOBAL(nglob).pre_react();
        KK_SR_H_TYPE(n) = 0;
        KK_SR_H_MAP(n) = nglob;
        nglob++;
      } else if (strcmp(surf->sr[n]->style,"prob") == 0) {
#ifdef SPARTA_KOKKOS_FIXED_LISTS
        if (nprob >= KOKKOS_MAX_SURF_REACT_PER_TYPE)
          error->all(FLERR,"Kokkos currently supports two instances of each surface reaction method");
        sr_kk_prob_copy[nprob].copy((SurfReactProbKokkos*)(surf->sr[n]));
#else
        sr_buf_blit(k_sr_prob,nprob,(SurfReactProbKokkos*)(surf->sr[n]));
#endif
        KK_SR_H_PROB(nprob).pre_react();
        KK_SR_H_TYPE(n) = 1;
        KK_SR_H_MAP(n) = nprob;
        nprob++;
      } else {
        error->all(FLERR,"Unknown Kokkos surface reaction method");
      }
    }

#ifndef SPARTA_KOKKOS_FIXED_LISTS

    // the models were blitted into the host image of the buffers and their
    //   pre_react() ran there; push the result to the device

    sr_buf_sync(k_sr_global,d_sr_global);
    sr_buf_sync(k_sr_prob,d_sr_prob);
    sr_idx_sync(k_sr_type_list,d_sr_type_list);
    sr_idx_sync(k_sr_map,d_sr_map);
#endif
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSurfKokkos::post_surf_tally()
{
  if (need_dup) {
    Kokkos::Experimental::contribute(d_array_surf_tally, dup_array_surf_tally);
    dup_array_surf_tally = {}; // free duplicated memory
  }

  k_tally2surf.modify_device();
  k_array_surf_tally.modify_device();

  // pre_surf_tally() has each active surf react model retain a reference to
  //  the particle list.  Release it before the next pre_surf_tally() blits
  //  over the member: a blit does not release, so the reference would be
  //  orphaned and its allocation never freed

  for (int n = 0; n < surf->nsr; n++) {
    if (KK_SR_H_TYPE(n) == 0) KK_SR_H_GLOBAL(KK_SR_H_MAP(n)).post_react();
    else KK_SR_H_PROB(KK_SR_H_MAP(n)).post_react();
  }
}

/* ----------------------------------------------------------------------
   return ptr to norm vector used by column N
------------------------------------------------------------------------- */

int ComputeSurfKokkos::tallyinfo(surfint *&ptr)
{
  // compressing below is destructive, so only do it once per clear() cycle

  if (compressed) {
    ptr = tally2surf;
    return ntally;
  }
  compressed = 1;

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
    while (istart < nsurf && h_surf2tally[istart] != -1) istart++;
    while (iend > 0 && h_surf2tally[iend] == -1) iend--;
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

  // the compression above rewrote both arrays on the host.  Claim it: the dense
  // list is what every consumer reads, and an unclaimed write is discarded by
  // the next sync_host() or by the realloc in init_normflux(), which resizes
  // whichever side the counters call newer.

  k_tally2surf.modify_host();
  k_array_surf_tally.modify_host();

  return ntally;
}

/* ----------------------------------------------------------------------
   sum tally values to owning surfs
   ComputeSurf::post_process_surf() collates ntally rows, but ntally is only
     computed by tallyinfo(), which is also what copies the device tallies to
     the host.  Only fix ave/surf calls tallyinfo(); dump surf, compute reduce
     and surf-style variables call post_process_surf() directly, and would
     otherwise collate the ntally = 0 left by clear() and report all zeroes.
------------------------------------------------------------------------- */

void ComputeSurfKokkos::post_process_surf()
{
  if (combined) return;

  surfint *ptr;
  tallyinfo(ptr);

  ComputeSurf::post_process_surf();
}

/* ---------------------------------------------------------------------- */

void ComputeSurfKokkos::grow_tally()
{
  // Cannot realloc inside a Kokkos parallel region, so size tally2surf the
  //  same as surf2tally

  int nsurf = surf->nlocal + surf->nghost;

  memoryKK->grow_kokkos(k_tally2surf,tally2surf,nsurf,"surf:tally2surf");
  d_tally2surf = k_tally2surf.view_device();

  memoryKK->grow_kokkos(k_array_surf_tally,array_surf_tally,nsurf,ntotal,"surf:array_surf_tally");
  d_array_surf_tally = k_array_surf_tally.view_device();
}


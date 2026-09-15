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

#include "stdio.h"
#include "region_union_kokkos.h"
#include "domain.h"
#include "particle_kokkos.h"
#include "error.h"
#include "sparta_masks.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

RegUnionKokkos::RegUnionKokkos(SPARTA *sparta, int narg, char **arg) :
  RegUnion(sparta, narg, arg)
{
  kokkos_flag = 1;
  ntoken = 0;
}

/* ---------------------------------------------------------------------- */

RegUnionKokkos::~RegUnionKokkos()
{
}

/* ----------------------------------------------------------------------
   flatten this region and its sub-regions into one device-resident postfix
     (RPN) token stream: sub-stream(0) sub-stream(1) OP sub-stream(2) OP ...
     followed by a NOT token when this composite is an exterior region
   a sub-region contributes a whole sub-stream, so a sub-region may itself be
     a region union or region intersect -- nesting is supported to any depth
     whose evaluation fits the RKK_MAX_STACK boolean stack, which is checked
     here on the host and errors out rather than truncating
------------------------------------------------------------------------- */

int RegUnionKokkos::flatten_region_kokkos(tdual_region_token_1d &k_tokens_out)
{
  Region **regions = domain->regions;

  ntoken = 0;

  for (int i = 0; i < nregion; i++) {
    Region *r = regions[list[i]];
    KokkosBase *rkk = dynamic_cast<KokkosBase*>(r);
    if (!rkk || !r->kokkos_flag)
      error->all(FLERR,"KOKKOS package does not (yet) support the region style "
                 "used inside region union");

    // k_sub is declared inside the loop on purpose: a composite sub-region
    //   hands back its own buffer, and a later primitive sub-region handed
    //   that same handle would write into it

    tdual_region_token_1d k_sub;
    const int nsub = rkk->flatten_region_kokkos(k_sub);
    if (nsub <= 0)
      error->all(FLERR,"KOKKOS package does not (yet) support the region style "
                 "used inside region union");

    // append the sub-region's stream, then the op that folds it into the
    //   running result (every sub-region past the first)

    region_token_grow(k_tokens,ntoken+nsub+2);
    for (int j = 0; j < nsub; j++)
      k_tokens.view_host()[ntoken++] = k_sub.view_host()[j];
    if (i) k_tokens.view_host()[ntoken++].type = RKK_TOK_UNION;
  }

  // this composite's own interior/exterior sense: !(hit ^ interior)

  if (!interior) {
    region_token_grow(k_tokens,ntoken+1);
    k_tokens.view_host()[ntoken++].type = RKK_TOK_NOT;
  }

  // bound the boolean stack the kernel will need, on the host, before any
  //   particle is tested against the stream

  const int depth = region_token_depth(k_tokens,ntoken);
  if (depth < 0)
    error->all(FLERR,"Internal error flattening region union for the KOKKOS package");
  if (depth > RKK_MAX_STACK) {
    char str[128];
    snprintf(str,sizeof(str),"Region union is nested too deeply for the KOKKOS package "
            "(needs a boolean stack of %d, max is %d)",depth,RKK_MAX_STACK);
    error->all(FLERR,str);
  }

  k_tokens.modify_host();
  k_tokens.sync_device();

  k_tokens_out = k_tokens;
  return ntoken;
}

/* ---------------------------------------------------------------------- */

void RegUnionKokkos::match_all_kokkos(DAT::tdual_int_1d k_match_in)
{
  tdual_region_token_1d k_tokens_local;
  const int ntoken_local = flatten_region_kokkos(k_tokens_local);

  d_match = k_match_in.view_device();
  ParticleKokkos* particleKK = (ParticleKokkos*) particle;
  particleKK->sync(Device, PARTICLE_MASK);
  d_particles = particleKK->k_particles.view_device();
  const int nlocal = particle->nlocal;

  auto l_tokens = k_tokens_local.view_device();
  auto l_match = d_match;
  auto l_particles = d_particles;
  const int l_ntoken = ntoken_local;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,nlocal),
    KOKKOS_LAMBDA(const int &i) {
      const double x = l_particles[i].x[0];
      const double y = l_particles[i].x[1];
      const double z = l_particles[i].x[2];
      l_match[i] = region_match_kk(l_tokens,l_ntoken,x,y,z);
    });
  copymode = 0;
  k_match_in.modify_device();
}

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
#include "compute_react_isurf_grid_kokkos.h"
#include "surf_kokkos.h"
#include "surf_react.h"
#include "grid.h"
#include "update.h"
#include "memory_kokkos.h"
#include "error.h"
#include "sparta_masks.h"
#include "kokkos.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

ComputeReactISurfGridKokkos::ComputeReactISurfGridKokkos(SPARTA *sparta, int narg, char **arg) :
  ComputeReactISurfGrid(sparta, narg, arg)
{
  kokkos_flag = 1;
}

ComputeReactISurfGridKokkos::ComputeReactISurfGridKokkos(SPARTA *sparta) :
  ComputeReactISurfGrid(sparta)
{
  copy = 1;
  uncopy = 0;
}

/* ---------------------------------------------------------------------- */

ComputeReactISurfGridKokkos::~ComputeReactISurfGridKokkos()
{
  if (copy) return;

  memoryKK->destroy_kokkos(k_tally2surf,tally2surf);
  memoryKK->destroy_kokkos(k_array_surf_tally,array_surf_tally);
}

/* ---------------------------------------------------------------------- */

void ComputeReactISurfGridKokkos::init()
{
  ComputeReactISurfGrid::init();

  // flatten reaction2col to device (only used when rpflag)

  if (rpflag) {
    int nreaction = surf->sr[isr]->nlist;
    d_reaction2col = DAT::t_int_2d("react/isurf/grid:reaction2col",nreaction,ntotal);
    auto h_r2c = Kokkos::create_mirror_view(d_reaction2col);
    for (int i = 0; i < nreaction; i++)
      for (int j = 0; j < ntotal; j++)
        h_r2c(i,j) = reaction2col[i][j];
    Kokkos::deep_copy(d_reaction2col,h_r2c);
  }

  // size per-surf tally storage to nsurf (implicit-surf ablation only shrinks)

  resize_device(surf->nlocal + surf->nghost);
}

/* ---------------------------------------------------------------------- */

void ComputeReactISurfGridKokkos::resize_device(int nsurf)
{
  if (nsurf < 1) nsurf = 1;

  memoryKK->grow_kokkos(k_tally2surf,tally2surf,nsurf,"react/isurf/grid:tally2surf");
  d_tally2surf = k_tally2surf.view_device();

  d_surf2tally = DAT::t_int_1d("react/isurf/grid:surf2tally",nsurf);
  Kokkos::deep_copy(d_surf2tally,-1);

  memoryKK->grow_kokkos(k_array_surf_tally,array_surf_tally,nsurf,ntotal,
                        "react/isurf/grid:array_surf_tally");
  d_array_surf_tally = k_array_surf_tally.view_device();

  nsurf_tally_alloc = nsurf;
}

/* ---------------------------------------------------------------------- */

void ComputeReactISurfGridKokkos::clear()
{
  // called by Update at beginning of timesteps surf tallying is done

  int nsurf = surf->nlocal + surf->nghost;
  if (nsurf > nsurf_tally_alloc) resize_device(nsurf);

  Kokkos::deep_copy(d_array_surf_tally,0);
  Kokkos::deep_copy(d_surf2tally,-1);

  ntally = 0;
  combined = 0;
}

/* ---------------------------------------------------------------------- */

void ComputeReactISurfGridKokkos::pre_surf_tally()
{
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

void ComputeReactISurfGridKokkos::post_surf_tally()
{
  if (need_dup) {
    Kokkos::Experimental::contribute(d_array_surf_tally, dup_array_surf_tally);
    dup_array_surf_tally = {};
  }

  k_tally2surf.modify_device();
  k_array_surf_tally.modify_device();
}

/* ----------------------------------------------------------------------
   sync device tallies to host and compress to dense list (ntally tallies)
   matches ComputeISurfGridKokkos::tallyinfo()
------------------------------------------------------------------------- */

int ComputeReactISurfGridKokkos::tallyinfo(surfint *&ptr)
{
  k_tally2surf.sync_host();
  ptr = tally2surf;

  k_array_surf_tally.sync_host();
  auto h_surf2tally = Kokkos::create_mirror_view(d_surf2tally);
  Kokkos::deep_copy(h_surf2tally,d_surf2tally);

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
    for (int k = 0; k < ntotal; k++)
      array_surf_tally[istart][k] = array_surf_tally[iend][k];
    h_surf2tally[istart] = h_surf2tally[iend];
    h_surf2tally[iend] = -1;
    tally2surf[istart] = tally2surf[iend];
  }

  return ntally;
}

/* ----------------------------------------------------------------------
   sync the device per-surf tally to the host (tallyinfo) before the host
   base class collates it to per-grid; consumers (e.g. fix ablate) read the
   compute directly via post_process_isurf_grid() rather than tallyinfo()
------------------------------------------------------------------------- */

void ComputeReactISurfGridKokkos::post_process_isurf_grid()
{
  if (combined) return;
  surfint *dummy;
  tallyinfo(dummy);
  ComputeReactISurfGrid::post_process_isurf_grid();
}

/* ---------------------------------------------------------------------- */

void ComputeReactISurfGridKokkos::grow_tally()
{
  resize_device(surf->nlocal + surf->nghost);
}

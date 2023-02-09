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

/* ----------------------------------------------------------------------
   Contributing author: Alan Stagg (Sandia)
------------------------------------------------------------------------- */

#include "compute_distsurf_grid_kokkos.h"
#include "update.h"
#include "grid_kokkos.h"
#include "surf_kokkos.h"
#include "domain.h"
#include "geometry_kokkos.h"
#include "math_extra_kokkos.h"
#include "memory_kokkos.h"
#include "kokkos.h"
#include "sparta_masks.h"

using namespace SPARTA_NS;

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

ComputeDistSurfGridKokkos::
ComputeDistSurfGridKokkos(SPARTA *sparta, int narg, char **arg) :
  ComputeDistSurfGrid(sparta, narg, arg)
{
  kokkos_flag = 1;
}

/* ---------------------------------------------------------------------- */

ComputeDistSurfGridKokkos::~ComputeDistSurfGridKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_vector_grid,vector_grid);
  vector_grid = NULL;
}

/* ---------------------------------------------------------------------- */

void ComputeDistSurfGridKokkos::compute_per_grid()
{
  if (sparta->kokkos->prewrap) {
    ComputeDistSurfGrid::compute_per_grid();
  } else {
    compute_per_grid_kokkos();
    k_vector_grid.modify_device();
    k_vector_grid.sync_host();
  }
}

/* ---------------------------------------------------------------------- */

void ComputeDistSurfGridKokkos::compute_per_grid_kokkos()
{
  SurfKokkos* surf_kk = (SurfKokkos*) surf;
  surf_kk->sync(Device,ALL_MASK);
  d_lines = surf_kk->k_lines.d_view;
  d_tris = surf_kk->k_tris.d_view;
  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  int ntotal = surf_kk->nlocal;

  k_eflag = DAT::tdual_int_1d("compute/distsurf/grid:eflag",ntotal);
  k_slist = DAT::tdual_int_1d("compute/distsurf/grid:slist",ntotal);
  h_eflag = k_eflag.h_view;
  h_slist = k_slist.h_view;
  d_slist = k_slist.d_view;
  d_eflag = k_eflag.d_view;

  invoked_per_grid = update->ntimestep;
  dim = domain->dimension;

  nsurf = 0;
  if (dim == 2) {
    for (int i = 0; i < ntotal; i++) {
      h_eflag[i] = 0;
      if (!(lines[i].mask & sgroupbit)) continue;
      if (MathExtraKokkos::dot3(lines[i].norm,sdir) <= 0.0) {
        h_eflag[i] = 1;
        h_slist[nsurf++] = i;
      }
    }
  } else {
    for (int i = 0; i < ntotal; i++) {
      h_eflag[i] = 0;
      if (!(tris[i].mask & sgroupbit)) continue;
      if (MathExtraKokkos::dot3(tris[i].norm,sdir) <= 0.0) {
        h_eflag[i] = 1;
        h_slist[nsurf++] = i;
      }
    }
  }
  k_eflag.modify_host();
  k_slist.modify_host();

  k_eflag.sync_device();
  k_slist.sync_device();

  d_sctr = DAT::t_float_1d_3("compute/distsurf/grid:sctr",nsurf);

  // pre-compute center point of each eligible surf
  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeDistSurfGrid_surf_centroid>(0,nsurf),*this);
  copymode = 0;

  // loop over my unsplit/split grid cells
  // if surfs in cell and any are eligible, dist = 0.0
  // else loop over eligible surfs in slist:
  //   if vector from cell center to surf center is against surf norm, skip it
  //   compute distance from cell to surf via Geometry method
  //   dist = minimum distance to any eligible surf
  // if assign dist to split cell, also assign dist to all its sub cells
  GridKokkos* grid_kk = ((GridKokkos*)grid);
  d_cells = grid_kk->k_cells.d_view;
  d_cinfo = grid_kk->k_cinfo.d_view;
  d_csurfs = grid_kk->d_csurfs;
  d_csubs = grid_kk->d_csubs;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeDistSurfGrid_surf_distance>(0,nglocal),*this);
  copymode = 0;

  memoryKK->destroy_kokkos(k_eflag);
  memoryKK->destroy_kokkos(k_slist);
  d_sctr = DAT::t_float_1d_3();
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeDistSurfGridKokkos::operator()(TagComputeDistSurfGrid_surf_centroid, const int &i) const {

  // compute surf centroids
  double invthird = 1.0/3.0;
  double *p1,*p2,*p3;

  int m = d_slist[i];
  if (dim == 2) {
    p1 = d_lines[m].p1;
    p2 = d_lines[m].p2;
    d_sctr(i,0) = 0.5 * (p1[0] + p2[0]);
    d_sctr(i,1) = 0.5 * (p1[1] + p2[1]);
    d_sctr(i,2) = 0.0;
  } else {
    p1 = d_tris[m].p1;
    p2 = d_tris[m].p2;
    p3 = d_tris[m].p3;
    d_sctr(i,0) = invthird * (p1[0] + p2[0] + p3[0]);
    d_sctr(i,1) = invthird * (p1[1] + p2[1] + p3[1]);
    d_sctr(i,2) = invthird * (p1[2] + p2[2] + p3[2]);
  }
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeDistSurfGridKokkos::operator()(TagComputeDistSurfGrid_surf_distance, const int &icell) const {
  int i,m,n;
  double dist,mindist;
  double *lo,*hi;
  double cctr[3],cell2surf[3];

  if (!(d_cinfo[icell].mask & groupbit)) return;
  if (d_cells[icell].nsplit < 1) return;

  if (d_cells[icell].nsurf) {
    n = d_cells[icell].nsurf;
    auto csurfs_begin = d_csurfs.row_map(icell);
    for (i = 0; i < n; i++) {
      m = d_csurfs.entries(csurfs_begin + i);
      if (d_eflag[m]) break;
    }

    // cell is overlapped, set dist = 0.0 and return
    // if split cell, also set vector_grid = 0.0 for sub-cells

    if (i < n) {
      d_vector[icell] = 0.0;
      if (d_cells[icell].nsplit > 1) {
        n = d_cells[icell].nsplit;
        int isplit = d_cells[icell].isplit;
        auto csubs_begin = d_csubs.row_map(isplit);
        for (i = 0; i < n; i++) {
          m = d_csubs.entries(csubs_begin + i);
          d_vector[m] = 0.0;
        }
      }
      return;
    }
  }

  lo = d_cells[icell].lo;
  hi = d_cells[icell].hi;
  cctr[0] = 0.5 * (lo[0]+hi[0]);
  cctr[1] = 0.5 * (lo[1]+hi[1]);
  if (dim == 3) cctr[2] = 0.5 * (lo[2]+hi[2]);
  else cctr[2] = 0.0;

  mindist = BIG;
  for (i = 0; i < nsurf; i++) {
    m = d_slist[i];

    cell2surf[0] = d_sctr(i,0) - cctr[0];
    cell2surf[1] = d_sctr(i,1) - cctr[1];
    cell2surf[2] = d_sctr(i,2) - cctr[2];

    if (dim == 2) {
      if (MathExtraKokkos::dot3(cell2surf,d_lines[m].norm) > 0.0) continue;
      dist = GeometryKokkos::dist_line_quad(d_lines[m].p1,d_lines[m].p2,lo,hi);
    } else {
      if (MathExtraKokkos::dot3(cell2surf,d_tris[m].norm) > 0.0) continue;
      dist = GeometryKokkos::dist_tri_hex(d_tris[m].p1,d_tris[m].p2,d_tris[m].p3,
                                          d_tris[m].norm,lo,hi);
    }
    mindist = MIN(mindist,dist);
  }

  d_vector[icell] = mindist;
}

/* ----------------------------------------------------------------------
   reallocate vector if nglocal has changed
   called by init() and load balancer
------------------------------------------------------------------------- */

void ComputeDistSurfGridKokkos::reallocate()
{
  if (grid->nlocal == nglocal) return;

  memoryKK->destroy_kokkos(k_vector_grid,vector_grid);

  memory->destroy(vector_grid);
  nglocal = grid->nlocal;
  memoryKK->create_kokkos(k_vector_grid,vector_grid,nglocal,"distsurf/grid:vector_grid");
  d_vector = k_vector_grid.d_view;
}

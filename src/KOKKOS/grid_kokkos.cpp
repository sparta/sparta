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

#include "string.h"
#include "grid_kokkos.h"
#include "geometry.h"
#include "domain.h"
#include "comm.h"
#include "irregular.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "error.h"
#include "kokkos.h"
#include "sparta_masks.h"

using namespace SPARTA_NS;
using namespace MathConst;

#define DELTA 8192
#define DELTAPARENT 1024
#define BIG 1.0e20
#define MAXGROUP 32
#define MAXLEVEL 32

// default value, can be overridden by global command

#define MAXSURFPERCELL 100

enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};         // same as Domain
enum{PERIODIC,OUTFLOW,REFLECT,SURFACE,AXISYM};  // same as Domain
enum{REGION_ALL,REGION_ONE,REGION_CENTER};      // same as Surf

// cell is entirely outside/inside surfs or has any overlap with surfs
// corner pt is outside/inside surfs or is on a surf

enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};           // several files
enum{NCHILD,NPARENT,NUNKNOWN,NPBCHILD,NPBPARENT,NPBUNKNOWN,NBOUND};  // Update
enum{NOWEIGHT,VOLWEIGHT,RADWEIGHT,RADUNWEIGHT};

// allocate space for static class variable

//Grid *Grid::gptr;

// corners[i][j] = J corner points of face I of a grid cell
// works for 2d quads and 3d hexes

//int corners[6][4] = {{0,2,4,6}, {1,3,5,7}, {0,1,4,5}, {2,3,6,7},
//                     {0,1,2,3}, {4,5,6,7}};

/* ---------------------------------------------------------------------- */

GridKokkos::GridKokkos(SPARTA *sparta) : Grid(sparta)
{
  delete [] plevels;
  memoryKK->create_kokkos(k_plevels,plevels,MAXLEVEL,"grid:plevels");
}

GridKokkos::~GridKokkos()
{
  if (copy || copymode) return;

  cells = NULL;
  cinfo = NULL;
  sinfo = NULL;
  pcells = NULL;
  plevels = NULL;
}

///////////////////////////////////////////////////////////////////////////
// grow cell list data structures
///////////////////////////////////////////////////////////////////////////

/* ----------------------------------------------------------------------
   insure cells and cinfo can hold N and M new cells respectively
------------------------------------------------------------------------- */

void GridKokkos::grow_cells(int n, int m)
{
  if (sparta->kokkos->prewrap) {
    Grid::grow_cells(n,m);
  } else {

    if (nlocal+nghost+n >= maxcell) {
      while (maxcell < nlocal+nghost+n) maxcell += DELTA;
      if (cells == NULL)
          k_cells = tdual_cell_1d("grid:cells",maxcell);
      else {
        this->sync(Device,CELL_MASK); // force resize on device
        k_cells.resize(maxcell);
        this->modify(Device,CELL_MASK); // needed for auto sync
      }
      cells = k_cells.h_view.data();
    }

    if (nlocal+m >= maxlocal) {
      while (maxlocal < nlocal+m) maxlocal += DELTA;
      if (cinfo == NULL)
          k_cinfo = tdual_cinfo_1d("grid:cinfo",maxlocal);
      else {
        this->sync(Device,CINFO_MASK); // force resize on device
        k_cinfo.resize(maxlocal);
        this->modify(Device,CINFO_MASK); // needed for auto sync
      }
      cinfo = k_cinfo.h_view.data();
    }
  }
}

/* ----------------------------------------------------------------------
   grow pcells
------------------------------------------------------------------------- */

void GridKokkos::grow_pcells()
{
  if (sparta->kokkos->prewrap) {
    Grid::grow_pcells();
  } else {

    maxparent += DELTA;
    if (pcells == NULL)
        k_pcells = tdual_pcell_1d("grid:pcells",maxparent);
    else {
      this->sync(Device,PCELL_MASK); // force resize on device
      k_pcells.resize(maxparent);
      this->modify(Device,PCELL_MASK); // needed for auto sync
    }
    pcells = k_pcells.h_view.data();
  }
}

/* ----------------------------------------------------------------------
   insure sinfo can hold N new split cells
------------------------------------------------------------------------- */

void GridKokkos::grow_sinfo(int n)
{
  if (sparta->kokkos->prewrap) {
    Grid::grow_sinfo(n);
  } else {

    if (nsplitlocal+nsplitghost+n >= maxsplit) {
      while (maxsplit < nsplitlocal+nsplitghost+n) maxsplit += DELTA;
      if (sinfo == NULL)
          k_sinfo = tdual_sinfo_1d("grid:sinfo",maxsplit);
      else {
        this->sync(Device,SINFO_MASK); // force resize on device
        k_sinfo.resize(maxsplit);
        this->modify(Device,SINFO_MASK); // needed for auto sync
      }
      sinfo = k_sinfo.h_view.data();
    }
  }
}

/* ---------------------------------------------------------------------- */

void GridKokkos::wrap_kokkos_graphs()
{
  if (!surf->exist) return;

  // csurfs

  Kokkos::Crs<int, SPAHostType, void, int> h_csurfs;
  auto csurfs_lambda = [=](int icell, int* fill) {
    int nsurf = cells[icell].nsurf;
    if (nsurf < 0) nsurf = 0;
    else if (fill) {
      surfint* csurfs = cells[icell].csurfs;
      // d_csurfs doesn't need to be surfint because at this point there are only
      //   local (not global) ids stored in csurfs
      for (int j = 0; j < nsurf; ++j) fill[j] = (int) csurfs[j];
    }
    return nsurf;
  };
  Kokkos::count_and_fill_crs(h_csurfs, nlocal+nghost, csurfs_lambda);
  d_csurfs.row_map = decltype(d_csurfs.row_map)(
      "csurfs.row_map", h_csurfs.row_map.size());
  d_csurfs.entries = decltype(d_csurfs.entries)(
      "csurfs.entries", h_csurfs.entries.size());
  Kokkos::deep_copy(d_csurfs.row_map, h_csurfs.row_map);
  Kokkos::deep_copy(d_csurfs.entries, h_csurfs.entries);

  if (sinfo != NULL) {

    Kokkos::Crs<int, SPAHostType, void, int> h_csplits;
    auto csplits_lambda = [=](int isplit, int* fill) {
      int icell = sinfo[isplit].icell;
      int nsurf = cells[icell].nsurf;
      int nsplit = cells[icell].nsplit;
      if (nsurf < 0 || nsplit <= 1) nsurf = 0;
      else if (fill) {
        int* csplits = sinfo[isplit].csplits;
        for (int j = 0; j < nsurf; ++j) fill[j] = csplits[j];
      }
      return nsurf;
    };
    Kokkos::count_and_fill_crs(h_csplits, nsplitlocal+nsplitghost, csplits_lambda);
    d_csplits.row_map = decltype(d_csplits.row_map)(
        "csplits.row_map", h_csplits.row_map.size());
    d_csplits.entries = decltype(d_csplits.entries)(
        "csplits.entries", h_csplits.entries.size());
    Kokkos::deep_copy(d_csplits.row_map, h_csplits.row_map);
    Kokkos::deep_copy(d_csplits.entries, h_csplits.entries);

    Kokkos::Crs<int, SPAHostType, void, int> h_csubs;
    auto csubs_lambda = [=](int isplit, int* fill) {
      int icell = sinfo[isplit].icell;
      int nsurf = cells[icell].nsurf;
      int nsplit = cells[icell].nsplit;
      if (nsurf < 0 || nsplit <= 1) nsplit = 0;
      else if (fill) {
        int* csubs = sinfo[isplit].csubs;
        for (int j = 0; j < nsplit; ++j) fill[j] = csubs[j];
      }
      return nsplit;
    };
    Kokkos::count_and_fill_crs(h_csubs, nsplitlocal+nsplitghost, csubs_lambda);
    d_csubs.row_map = decltype(d_csubs.row_map)(
        "csubs.row_map", h_csubs.row_map.size());
    d_csubs.entries = decltype(d_csubs.entries)(
        "csubs.entries", h_csubs.entries.size());
    Kokkos::deep_copy(d_csubs.row_map, h_csubs.row_map);
    Kokkos::deep_copy(d_csubs.entries, h_csubs.entries);

  }
}

/* ---------------------------------------------------------------------- */

void GridKokkos::wrap_kokkos()
{
  // cells

  if (cells != k_cells.h_view.data()) {
    memoryKK->wrap_kokkos(k_cells,cells,maxcell,"grid:cells");
    k_cells.modify_host();
    k_cells.sync_device();
    memory->sfree(cells);
    cells = k_cells.h_view.data();
  }

  // cinfo

  if (cinfo != k_cinfo.h_view.data()) {
    memoryKK->wrap_kokkos(k_cinfo,cinfo,maxlocal,"grid:cinfo");
    k_cinfo.modify_host();
    k_cinfo.sync_device();
    memory->sfree(cinfo);
    cinfo = k_cinfo.h_view.data();
  }

  // sinfo

  if (sinfo != k_sinfo.h_view.data()) {
    memoryKK->wrap_kokkos(k_sinfo,sinfo,maxsplit,"grid:sinfo");
    k_sinfo.modify_host();
    k_sinfo.sync_device();
    memory->sfree(sinfo);
    sinfo = k_sinfo.h_view.data();
  }

  wrap_kokkos_graphs();

  // pcells

  if (pcells != k_pcells.h_view.data()) {
    memoryKK->wrap_kokkos(k_pcells,pcells,maxparent,"grid:pcells");
    k_pcells.modify_host();
    k_pcells.sync_device();
    memory->sfree(pcells);
    pcells = k_pcells.h_view.data();
  }

  // plevels doesn't need wrap but was modified on host

  k_plevels.modify_host();
  k_plevels.sync_device();
}

/* ---------------------------------------------------------------------- */

void GridKokkos::sync(ExecutionSpace space, unsigned int mask)
{
  if (sparta->kokkos->prewrap) {
    if (space == Device)
      error->one(FLERR,"Sync Device before wrap");
    else
      return;
  }

  if (space == Device) {
    if (sparta->kokkos->auto_sync)
      modify(Host,mask);
    if (mask & CELL_MASK) k_cells.sync_device();
    if (mask & CINFO_MASK) k_cinfo.sync_device();
    if (mask & PCELL_MASK) k_pcells.sync_device();
    if (mask & SINFO_MASK) k_sinfo.sync_device();
    if (mask & PLEVEL_MASK) k_plevels.sync_device();
  } else {
    if (mask & CELL_MASK) k_cells.sync_host();
    if (mask & CINFO_MASK) k_cinfo.sync_host();
    if (mask & PCELL_MASK) k_pcells.sync_host();
    if (mask & SINFO_MASK) k_sinfo.sync_host();
    if (mask & PLEVEL_MASK) k_plevels.sync_host();
  }
}

/* ---------------------------------------------------------------------- */

void GridKokkos::modify(ExecutionSpace space, unsigned int mask)
{
  if (sparta->kokkos->prewrap) {
    if (space == Device)
      error->one(FLERR,"Modify Device before wrap");
    else
      return;
  }

  if (space == Device) {
    if (mask & CELL_MASK) k_cells.modify_device();
    if (mask & CINFO_MASK) k_cinfo.modify_device();
    if (mask & PCELL_MASK) k_pcells.modify_device();
    if (mask & SINFO_MASK) k_sinfo.modify_device();
    if (mask & PLEVEL_MASK) k_plevels.modify_device();
    if (sparta->kokkos->auto_sync)
      sync(Host,mask);
  } else {
    if (mask & CELL_MASK) k_cells.modify_host();
    if (mask & CINFO_MASK) k_cinfo.modify_host();
    if (mask & PCELL_MASK) k_pcells.modify_host();
    if (mask & SINFO_MASK) k_sinfo.modify_host();
    if (mask & PLEVEL_MASK) k_plevels.modify_host();
  }
}

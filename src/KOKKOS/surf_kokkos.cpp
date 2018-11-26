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

#include "ctype.h"
#include "surf_kokkos.h"
#include "style_surf_collide.h"
#include "style_surf_react.h"
#include "surf_collide.h"
#include "surf_react.h"
#include "domain.h"
#include "region.h"
#include "comm.h"
#include "geometry.h"
#include "input.h"
#include "math_extra.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "error.h"
#include "sparta_masks.h"
#include "kokkos.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{TALLYAUTO,TALLYREDUCE,TALLYLOCAL};         // same as Update
enum{REGION_ALL,REGION_ONE,REGION_CENTER};      // same as Grid
enum{TYPE,MOLECULE,ID};
enum{LT,LE,GT,GE,EQ,NEQ,BETWEEN};

#define DELTA 4
#define EPSSQ 1.0e-12
#define BIG 1.0e20
#define MAXGROUP 32

/* ---------------------------------------------------------------------- */

SurfKokkos::SurfKokkos(SPARTA *sparta) : Surf(sparta)
{


}

/* ---------------------------------------------------------------------- */

SurfKokkos::~SurfKokkos()
{
  lines = NULL;
  tris = NULL;
}

/* ---------------------------------------------------------------------- */

void SurfKokkos::wrap_kokkos()
{
  if (lines != NULL && nline > 0) {
    if (lines != k_lines.h_view.data()) {
      memoryKK->wrap_kokkos(k_lines,lines,nline,"surf:lines");
      k_lines.modify<SPAHostType>();
      k_lines.sync<DeviceType>();
      memory->sfree(lines);
      lines = k_lines.h_view.data();
    }
  }

  if (tris != NULL && ntri > 0) {
    if (tris != k_tris.h_view.data()) {
      memoryKK->wrap_kokkos(k_tris,tris,ntri,"surf:tris");
      k_tris.modify<SPAHostType>();
      k_tris.sync<DeviceType>();
      memory->sfree(tris);
      tris = k_tris.h_view.data();
    }
  }
}

/* ----------------------------------------------------------------------
   grow surface data structures
------------------------------------------------------------------------- */
void SurfKokkos::grow()
{
  if (sparta->kokkos->prewrap) {
    Surf::grow();
  } else {
    SurfKokkos* surf_kk = (SurfKokkos*) surf;

    if (lines == NULL)
        surf_kk->k_lines = tdual_line_1d("surf:lines",nline);
    else {
      surf_kk->sync(Host,LINE_MASK);
      surf_kk->modify(Host,LINE_MASK); // force resize on host
      surf_kk->k_lines.resize(nline);
    }
    lines = surf_kk->k_lines.h_view.data();

    if (tris == NULL)
        surf_kk->k_tris = tdual_tri_1d("surf:tris",ntri);
    else {
      surf_kk->sync(Host,TRI_MASK);
      surf_kk->modify(Host,TRI_MASK); // force resize on host
      surf_kk->k_tris.resize(ntri);
    }
    tris = surf_kk->k_tris.h_view.data();
  }
}

/* ---------------------------------------------------------------------- */

void SurfKokkos::sync(ExecutionSpace space, unsigned int mask)
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
    if (mask & LINE_MASK) k_lines.sync<SPADeviceType>();
    if (mask & TRI_MASK) k_tris.sync<SPADeviceType>();
  } else {
    if (mask & LINE_MASK) k_lines.sync<SPAHostType>();
    if (mask & TRI_MASK) k_tris.sync<SPAHostType>();
  }
}

/* ---------------------------------------------------------------------- */

void SurfKokkos::modify(ExecutionSpace space, unsigned int mask)
{
  if (sparta->kokkos->prewrap) {
    if (space == Device)
      error->one(FLERR,"Modify Device before wrap");
    else
      return;
  }

  if (space == Device) {
    if (mask & LINE_MASK) k_lines.modify<SPADeviceType>();
    if (mask & TRI_MASK) k_tris.modify<SPADeviceType>();
    if (sparta->kokkos->auto_sync)
      sync(Host,mask);
  } else {
    if (mask & LINE_MASK) k_lines.modify<SPAHostType>();
    if (mask & TRI_MASK) k_tris.modify<SPAHostType>();
  }
}

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

  mylines = NULL;
  mytris = NULL;
}

/* ---------------------------------------------------------------------- */

void SurfKokkos::wrap_kokkos()
{
  if (domain->dimension == 2) {
    if (lines != NULL && nmax > 0) {
      if (lines != k_lines.h_view.data()) {
        memoryKK->wrap_kokkos(k_lines,lines,nmax,"surf:lines");
        k_lines.modify<SPAHostType>();
        k_lines.sync<DeviceType>();
        memory->sfree(lines);
        lines = k_lines.h_view.data();
      }
    }
    if (mylines != NULL && nown > 0) {
      if (mylines != k_mylines.h_view.data()) {
        memoryKK->wrap_kokkos(k_mylines,mylines,nown,"surf:lines");
        k_mylines.modify<SPAHostType>();
        k_mylines.sync<DeviceType>();
        memory->sfree(mylines);
        mylines = k_mylines.h_view.data();
      }
    }
  } else {
    if (tris != NULL && nmax > 0) {
      if (tris != k_tris.h_view.data()) {
        memoryKK->wrap_kokkos(k_tris,tris,nmax,"surf:tris");
        k_tris.modify<SPAHostType>();
        k_tris.sync<DeviceType>();
        memory->sfree(tris);
        tris = k_tris.h_view.data();
      }
    }
    if (mytris != NULL && nown > 0) {
      if (mytris != k_mytris.h_view.data()) {
        memoryKK->wrap_kokkos(k_mytris,mytris,nown,"surf:mytris");
        k_mytris.modify<SPAHostType>();
        k_mytris.sync<DeviceType>();
        memory->sfree(mytris);
        mytris = k_mytris.h_view.data();
      }
    }
  }
}

/* ----------------------------------------------------------------------
   grow surface data structures
------------------------------------------------------------------------- */
void SurfKokkos::grow(int old)
{
  if (sparta->kokkos->prewrap) {
    Surf::grow(old);
  } else {
    SurfKokkos* surf_kk = (SurfKokkos*) surf;

    if (domain->dimension == 2) {
      if (lines == NULL)
          surf_kk->k_lines = tdual_line_1d("surf:lines",nmax);
      else {
        surf_kk->sync(Host,LINE_MASK);
        surf_kk->modify(Host,LINE_MASK); // force resize on host
        surf_kk->k_lines.resize(nmax);
      }
      lines = surf_kk->k_lines.h_view.data();
    } else {
      if (tris == NULL)
          surf_kk->k_tris = tdual_tri_1d("surf:tris",nmax);
      else {
        surf_kk->sync(Host,TRI_MASK);
        surf_kk->modify(Host,TRI_MASK); // force resize on host
        surf_kk->k_tris.resize(nmax);
      }
      tris = surf_kk->k_tris.h_view.data();
    }
  }
}


/* ----------------------------------------------------------------------
   grow surface data structures
------------------------------------------------------------------------- */
void SurfKokkos::grow_own(int old)
{
  if (sparta->kokkos->prewrap) {
    Surf::grow_own(old);
  } else {
    SurfKokkos* surf_kk = (SurfKokkos*) surf;

    if (domain->dimension == 2) {
      if (mylines == NULL)
          surf_kk->k_mylines = tdual_line_1d("surf:mylines",nown);
      else {
        surf_kk->sync(Host,LINE_MASK);
        surf_kk->modify(Host,LINE_MASK); // force resize on host
        surf_kk->k_mylines.resize(nown);
      }
      mylines = surf_kk->k_mylines.h_view.data();
    } else {
      if (mytris == NULL)
          surf_kk->k_mytris = tdual_tri_1d("surf:mytris",nown);
      else {
        surf_kk->sync(Host,TRI_MASK);
        surf_kk->modify(Host,TRI_MASK); // force resize on host
        surf_kk->k_mytris.resize(nown);
      }
      mytris = surf_kk->k_mytris.h_view.data();
    }
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

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

  // deallocate views of views in serial to prevent race condition in profiling tools

  for (int i = 0; i < k_eivec.extent(0); i++)
    k_eivec.h_view(i).k_view = {};

  for (int i = 0; i < k_eiarray.extent(0); i++)
    k_eiarray.h_view(i).k_view = {};

  for (int i = 0; i < k_edvec.extent(0); i++)
    k_edvec.h_view(i).k_view = {};

  for (int i = 0; i < k_edarray.extent(0); i++)
    k_edarray.h_view(i).k_view = {};

  for (int i = 0; i < k_eivec_local.extent(0); i++)
    k_eivec_local.h_view(i).k_view = {};

  for (int i = 0; i < k_eiarray_local.extent(0); i++)
    k_eiarray_local.h_view(i).k_view = {};

  for (int i = 0; i < k_edvec_local.extent(0); i++)
    k_edvec_local.h_view(i).k_view = {};

  for (int i = 0; i < k_edarray_local.extent(0); i++)
    k_edarray_local.h_view(i).k_view = {};

  eivec = NULL;
  eiarray = NULL;
  edvec = NULL;
  edarray = NULL;

  ewhich = NULL;
  eicol = NULL;
  edcol = NULL;

  ncustom_ivec = ncustom_iarray = 0;
  ncustom_dvec = ncustom_darray = 0;
}

/* ---------------------------------------------------------------------- */

void SurfKokkos::clear_explicit()
{
  nsurf = 0;
  nlocal = nghost = nmax = 0;
  nown = maxown = 0;

  k_lines = {};
  k_tris = {};
  k_mylines = {};
  k_mytris = {};

  lines = NULL;
  tris = NULL;
  mylines = NULL;
  mytris = NULL;

  hash->clear();
  hashfilled = 0;
}

/* ---------------------------------------------------------------------- */

void SurfKokkos::wrap_kokkos()
{
  if (domain->dimension == 2) {
    if (lines != NULL && nmax > 0) {
      if (lines != k_lines.h_view.data()) {
        memoryKK->wrap_kokkos(k_lines,lines,nmax,"surf:lines");
        k_lines.modify_host();
        k_lines.sync_device();
        memory->sfree(lines);
        lines = k_lines.h_view.data();
      }
    }
    if (mylines != NULL && nown > 0) {
      if (mylines != k_mylines.h_view.data()) {
        memoryKK->wrap_kokkos(k_mylines,mylines,nown,"surf:lines");
        k_mylines.modify_host();
        k_mylines.sync_device();
        memory->sfree(mylines);
        mylines = k_mylines.h_view.data();
      }
    }
  } else {
    if (tris != NULL && nmax > 0) {
      if (tris != k_tris.h_view.data()) {
        memoryKK->wrap_kokkos(k_tris,tris,nmax,"surf:tris");
        k_tris.modify_host();
        k_tris.sync_device();
        memory->sfree(tris);
        tris = k_tris.h_view.data();
      }
    }
    if (mytris != NULL && nown > 0) {
      if (mytris != k_mytris.h_view.data()) {
        memoryKK->wrap_kokkos(k_mytris,mytris,nown,"surf:mytris");
        k_mytris.modify_host();
        k_mytris.sync_device();
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
  if (sparta->kokkos->prewrap)
    if (space == Device)
      error->one(FLERR,"Sync Device before wrap");

  if (space == Device) {
    if (sparta->kokkos->auto_sync)
      modify(Host,mask);
    if (mask & LINE_MASK) k_lines.sync_device();
    if (mask & TRI_MASK) k_tris.sync_device();
    if (mask & CUSTOM_MASK) {
      if (ncustom) {
        if (ncustom_ivec)
          for (int i = 0; i < ncustom_ivec; i++)
            k_eivec.h_view[i].k_view.sync_device();

        if (ncustom_iarray)
          for (int i = 0; i < ncustom_iarray; i++)
            k_eiarray.h_view[i].k_view.sync_device();

        if (ncustom_dvec)
          for (int i = 0; i < ncustom_dvec; i++)
            k_edvec.h_view[i].k_view.sync_device();

        if (ncustom_darray)
          for (int i = 0; i < ncustom_darray; i++)
            k_edarray.h_view[i].k_view.sync_device();
      }
    }
  } else {
    if (mask & LINE_MASK) k_lines.sync_host();
    if (mask & TRI_MASK) k_tris.sync_host();
    if (mask & CUSTOM_MASK) {
      if (ncustom_ivec)
        for (int i = 0; i < ncustom_ivec; i++)
          k_eivec.h_view[i].k_view.sync_host();

      if (ncustom_iarray)
        for (int i = 0; i < ncustom_iarray; i++)
          k_eiarray.h_view[i].k_view.sync_host();

      if (ncustom_dvec)
        for (int i = 0; i < ncustom_dvec; i++)
          k_edvec.h_view[i].k_view.sync_host();

      if (ncustom_darray)
        for (int i = 0; i < ncustom_darray; i++)
          k_edarray.h_view[i].k_view.sync_host();
    }
  }
}

/* ---------------------------------------------------------------------- */

void SurfKokkos::modify(ExecutionSpace space, unsigned int mask)
{
  if (sparta->kokkos->prewrap)
    if (space == Device)
      error->one(FLERR,"Modify Device before wrap");

  if (space == Device) {
    if (mask & LINE_MASK) k_lines.modify_device();
    if (mask & TRI_MASK) k_tris.modify_device();
    if (mask & CUSTOM_MASK) {
      if (ncustom) {
        if (ncustom_ivec)
          for (int i = 0; i < ncustom_ivec; i++)
            k_eivec.h_view[i].k_view.modify_device();

        if (ncustom_iarray)
          for (int i = 0; i < ncustom_iarray; i++)
            k_eiarray.h_view[i].k_view.modify_device();

        if (ncustom_dvec)
          for (int i = 0; i < ncustom_dvec; i++)
            k_edvec.h_view[i].k_view.modify_device();

        if (ncustom_darray)
          for (int i = 0; i < ncustom_darray; i++)
            k_edarray.h_view[i].k_view.modify_device();
      }
    }

    if (sparta->kokkos->auto_sync)
      sync(Host,mask);
  } else {
    if (mask & LINE_MASK) k_lines.modify_host();
    if (mask & TRI_MASK) k_tris.modify_host();
    if (mask & CUSTOM_MASK) {
      if (ncustom) {
        if (ncustom_ivec)
          for (int i = 0; i < ncustom_ivec; i++)
            k_eivec.h_view[i].k_view.modify_host();

        if (ncustom_iarray)
          for (int i = 0; i < ncustom_iarray; i++)
            k_eiarray.h_view[i].k_view.modify_host();

        if (ncustom_dvec)
          for (int i = 0; i < ncustom_dvec; i++)
            k_edvec.h_view[i].k_view.modify_host();

        if (ncustom_darray)
          for (int i = 0; i < ncustom_darray; i++)
            k_edarray.h_view[i].k_view.modify_host();
      }
    }
  }
}

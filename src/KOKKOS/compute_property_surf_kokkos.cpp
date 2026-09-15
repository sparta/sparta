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
#include "compute_property_surf_kokkos.h"
#include "surf_kokkos.h"
#include "domain.h"
#include "update.h"
#include "memory_kokkos.h"
#include "error.h"
#include "sparta_masks.h"
#include "kokkos.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

ComputePropertySurfKokkos::ComputePropertySurfKokkos(SPARTA *sparta, int narg, char **arg) :
  ComputePropertySurf(sparta, narg, arg)
{
  kokkos_flag = 1;

  // map field keywords to device index enum (must match base parse order)

  d_index = DAT::t_int_1d("property/surf:index",nvalues);
  auto h_index = Kokkos::create_mirror_view(d_index);
  for (int i = 0; i < nvalues; i++) {
    char *a = arg[3+i];
    int idx = -1;
    if (strcmp(a,"id") == 0) idx = ID;
    else if (strcmp(a,"v1x") == 0) idx = V1X;
    else if (strcmp(a,"v1y") == 0) idx = V1Y;
    else if (strcmp(a,"v1z") == 0) idx = V1Z;
    else if (strcmp(a,"v2x") == 0) idx = V2X;
    else if (strcmp(a,"v2y") == 0) idx = V2Y;
    else if (strcmp(a,"v2z") == 0) idx = V2Z;
    else if (strcmp(a,"v3x") == 0) idx = V3X;
    else if (strcmp(a,"v3y") == 0) idx = V3Y;
    else if (strcmp(a,"v3z") == 0) idx = V3Z;
    else if (strcmp(a,"xc") == 0) idx = XC;
    else if (strcmp(a,"yc") == 0) idx = YC;
    else if (strcmp(a,"zc") == 0) idx = ZC;
    else if (strcmp(a,"area") == 0) idx = AREA;
    else if (strcmp(a,"normx") == 0) idx = NORMX;
    else if (strcmp(a,"normy") == 0) idx = NORMY;
    else if (strcmp(a,"normz") == 0) idx = NORMZ;
    h_index(i) = idx;
  }
  Kokkos::deep_copy(d_index,h_index);
}

/* ---------------------------------------------------------------------- */

ComputePropertySurfKokkos::~ComputePropertySurfKokkos()
{
}

/* ---------------------------------------------------------------------- */

void ComputePropertySurfKokkos::init()
{
  ComputePropertySurf::init();

  // copy cglobal (owned-in-group surf indices) to device

  d_cglobal = DAT::t_int_1d("property/surf:cglobal",MAX(nchoose,1));
  auto h_cglobal = Kokkos::create_mirror_view(d_cglobal);
  for (int i = 0; i < nchoose; i++) h_cglobal(i) = cglobal[i];
  Kokkos::deep_copy(d_cglobal,h_cglobal);

  // device output storage (sized nsown to match host vector_surf/array_surf)

  int n = MAX(nsown,1);
  if (nvalues == 1) {
    d_vector_surf = DAT::t_float_1d("property/surf:vector_surf",n);
    k_vector_surf = DAT::tdual_float_1d("property/surf:vector_surf",n);
  } else {
    d_array_surf = DAT::t_float_2d_lr("property/surf:array_surf",n,nvalues);
    k_array_surf = DAT::tdual_float_2d_lr("property/surf:array_surf",n,nvalues);
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertySurfKokkos::compute_per_surf()
{
  if (sparta->kokkos->prewrap) {
    ComputePropertySurf::compute_per_surf();
  } else {
    compute_per_surf_kokkos();
    if (nvalues == 1) {
      Kokkos::deep_copy(k_vector_surf.view_device(),d_vector_surf);
      k_vector_surf.modify_device();
      k_vector_surf.sync_host();
      auto h = k_vector_surf.view_host();
      for (int i = 0; i < nsown; i++) vector_surf[i] = h(i);
    } else {
      Kokkos::deep_copy(k_array_surf.view_device(),d_array_surf);
      k_array_surf.modify_device();
      k_array_surf.sync_host();
      auto h = k_array_surf.view_host();
      for (int i = 0; i < nsown; i++)
        for (int n = 0; n < nvalues; n++) array_surf[i][n] = h(i,n);
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertySurfKokkos::compute_per_surf_kokkos()
{
  invoked_per_surf = update->ntimestep;

  dim = domain->dimension;

  SurfKokkos* surf_kk = (SurfKokkos*) surf;
  surf_kk->sync(Device,ALL_MASK);
  if (distributed) {
    d_lines = surf_kk->k_mylines.view_device();
    d_tris = surf_kk->k_mytris.view_device();
  } else {
    d_lines = surf_kk->k_lines.view_device();
    d_tris = surf_kk->k_tris.view_device();
  }

  if (nvalues == 1) Kokkos::deep_copy(d_vector_surf,0.0);
  else Kokkos::deep_copy(d_array_surf,0.0);

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,nchoose),*this);
  copymode = 0;
}

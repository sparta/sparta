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
#include "compute_boundary_kokkos.h"
#include "particle_kokkos.h"
#include "mixture.h"
#include "grid.h"
#include "surf.h"
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "memory_kokkos.h"
#include "error.h"
#include "sparta_masks.h"
#include "kokkos.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

ComputeBoundaryKokkos::ComputeBoundaryKokkos(SPARTA *sparta, int narg, char **arg) :
  ComputeBoundary(sparta, narg, arg)
{
  kokkos_flag = 1;

  memory->destroy(array);
  memory->destroy(myarray);
  memoryKK->create_kokkos(k_array,array,size_array_rows,size_array_cols,"boundary:array");
  memoryKK->create_kokkos(k_myarray,myarray,size_array_rows,size_array_cols,"boundary:myarray");
  d_myarray = k_myarray.d_view;
  d_array = k_array.d_view;
  d_which = DAT::t_int_1d("boundary:which",nvalue);
}

ComputeBoundaryKokkos::ComputeBoundaryKokkos(SPARTA *sparta) :
  ComputeBoundary(sparta)
{
  array = NULL;
  myarray = NULL;
  which = NULL;
  id = NULL;
  style = NULL;
  tlist = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeBoundaryKokkos::~ComputeBoundaryKokkos()
{
  if (copy || copymode) return;

  memoryKK->destroy_kokkos(k_array,array);
  memoryKK->destroy_kokkos(k_myarray,myarray);
}

/* ---------------------------------------------------------------------- */

void ComputeBoundaryKokkos::init()
{
  ComputeBoundary::init();

  auto h_which = Kokkos::create_mirror_view(d_which);
  for (int n=0; n<nvalue; n++)
    h_which(n) = which[n];
  Kokkos::deep_copy(d_which,h_which);
}

/* ---------------------------------------------------------------------- */

void ComputeBoundaryKokkos::compute_array()
{
  invoked_array = update->ntimestep;

  // sum tallies across processors

  if (sparta->kokkos->gpu_direct_flag) {
    MPI_Allreduce(d_myarray.data(),d_array.data(),nrow*ntotal,
                  MPI_DOUBLE,MPI_SUM,world);
    k_array.modify_device();
    k_array.sync_host();
  } else {
    k_myarray.modify_device();
    k_myarray.sync_host();
    MPI_Allreduce(k_myarray.h_view.data(),k_array.h_view.data(),nrow*ntotal,
                  MPI_DOUBLE,MPI_SUM,world);
  }

  // normalize tallies

  int m;
  for (int j = 0; j < ntotal; j++) {
    m = j % nvalue;
    if (which[m] != NUM && which[m] != NUMWT)
      for (int i = 0; i < size_array_rows; i++)
        array[i][j] /= normflux[i];
  }
}

/* ---------------------------------------------------------------------- */

void ComputeBoundaryKokkos::clear()
{
  // reset tally values to zero
  // called by Update at beginning of timesteps boundary tallying is done

  Kokkos::deep_copy(d_myarray,0.0);
}

void ComputeBoundaryKokkos::pre_boundary_tally()
{
  mvv2e = update->mvv2e;

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,PARTICLE_MASK|SPECIES_MASK);
  d_species = particle_kk->k_species.d_view;
  d_s2g = particle_kk->k_species2group.d_view;

  need_dup = sparta->kokkos->need_dup<DeviceType>();
  if (need_dup)
    dup_myarray = Kokkos::Experimental::create_scatter_view<typename Kokkos::Experimental::ScatterSum, typename Kokkos::Experimental::ScatterDuplicated>(d_myarray);
  else
    ndup_myarray = Kokkos::Experimental::create_scatter_view<typename Kokkos::Experimental::ScatterSum, typename Kokkos::Experimental::ScatterNonDuplicated>(d_myarray);
}

void ComputeBoundaryKokkos::post_boundary_tally()
{
  if (need_dup) {
    Kokkos::Experimental::contribute(d_myarray, dup_myarray);
    dup_myarray = decltype(dup_myarray)(); // free duplicated memory
  }
}

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
#include "compute_react_boundary_kokkos.h"
#include "surf.h"
#include "surf_react.h"
#include "update.h"
#include "domain.h"
#include "comm.h"
#include "memory_kokkos.h"
#include "error.h"
#include "kokkos.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

ComputeReactBoundaryKokkos::
ComputeReactBoundaryKokkos(SPARTA *sparta, int narg, char **arg) :
  ComputeReactBoundary(sparta, narg, arg)
{
  kokkos_flag = 1;

  memory->destroy(myarray);
  memoryKK->create_kokkos(k_myarray,myarray,size_array_rows,size_array_cols,
                          "react/boundary:myarray");
  d_myarray = k_myarray.view_device();

  // reaction2col never changes after construction, so flatten it once
  //   rpflag = 0 leaves it unallocated on the host; the kernel does not read
  //   it in that case, but give it one element so the view is always valid

  if (rpflag) {
    DAT::tdual_int_2d k_r2c("react/boundary:reaction2col",
                            surf->sr[isr]->nlist,ntotal);
    for (int i = 0; i < surf->sr[isr]->nlist; i++)
      for (int j = 0; j < ntotal; j++)
        k_r2c.view_host()(i,j) = reaction2col[i][j];
    k_r2c.modify_host();
    k_r2c.sync_device();
    d_reaction2col = k_r2c.view_device();
  } else {
    d_reaction2col = DAT::t_int_2d("react/boundary:reaction2col",1,1);
  }

  d_surf_react = DAT::t_int_1d("react/boundary:surf_react",6);
}

/* ---------------------------------------------------------------------- */

ComputeReactBoundaryKokkos::ComputeReactBoundaryKokkos(SPARTA *sparta) :
  ComputeReactBoundary(sparta)
{
  copy = 1;
}

/* ---------------------------------------------------------------------- */

ComputeReactBoundaryKokkos::~ComputeReactBoundaryKokkos()
{
  if (copy || copymode) return;

  memoryKK->destroy_kokkos(k_myarray,myarray);
  myarray = NULL;
}

/* ---------------------------------------------------------------------- */

void ComputeReactBoundaryKokkos::init()
{
  if (!domain->surfreactany && comm->me == 0)
    error->warning(FLERR,"Using compute react/boundary "
                   "when no box faces are assigned a reaction model");

  clear();
}

/* ----------------------------------------------------------------------
   called by Update at the start of any timestep boundary tallying is done
------------------------------------------------------------------------- */

void ComputeReactBoundaryKokkos::clear()
{
  Kokkos::deep_copy(d_myarray,0.0);
}

/* ----------------------------------------------------------------------
   called by UpdateKokkos before the move kernel
------------------------------------------------------------------------- */

void ComputeReactBoundaryKokkos::pre_boundary_tally()
{
  // domain->surf_react is a small host array indexed by box face, and
  //   bound_modify can change it between runs, so refresh it each step

  auto h_surf_react = Kokkos::create_mirror_view(d_surf_react);
  for (int i = 0; i < 6; i++) h_surf_react(i) = domain->surf_react[i];
  Kokkos::deep_copy(d_surf_react,h_surf_react);

  need_dup = sparta->kokkos->need_dup<DeviceType>();
  if (need_dup)
    dup_myarray = Kokkos::Experimental::create_scatter_view<typename Kokkos::Experimental::ScatterSum, typename Kokkos::Experimental::ScatterDuplicated>(d_myarray);
  else
    ndup_myarray = Kokkos::Experimental::create_scatter_view<typename Kokkos::Experimental::ScatterSum, typename Kokkos::Experimental::ScatterNonDuplicated>(d_myarray);
}

/* ----------------------------------------------------------------------
   called by UpdateKokkos after the move kernel
------------------------------------------------------------------------- */

void ComputeReactBoundaryKokkos::post_boundary_tally()
{
  if (need_dup) {
    Kokkos::Experimental::contribute(d_myarray, dup_myarray);
    dup_myarray = {}; // free duplicated memory
  }
}

/* ----------------------------------------------------------------------
   sum tallies across processors, as the host version does
------------------------------------------------------------------------- */

void ComputeReactBoundaryKokkos::compute_array()
{
  invoked_array = update->ntimestep;

  // the reduction is a host MPI call, so bring the device tallies back first

  k_myarray.modify_device();
  k_myarray.sync_host();

  MPI_Allreduce(&myarray[0][0],&array[0][0],nrow*ntotal,
                MPI_DOUBLE,MPI_SUM,world);
}

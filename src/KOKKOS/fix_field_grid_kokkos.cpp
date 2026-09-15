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

#include "fix_field_grid_kokkos.h"
#include "grid.h"
#include "memory_kokkos.h"
#include <type_traits>
#include "sparta_masks.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

FixFieldGridKokkos::FixFieldGridKokkos(SPARTA *sparta, int narg, char **arg) :
  FixFieldGrid(sparta, narg, arg)
{
  // the variable evaluation in compute_field() runs on the host, so this fix
  //   is a host fix.  it is not invoked through Modify -- UpdateKokkos calls
  //   compute_field() directly -- so the datamasks only describe the fact
  //   that no particle data is read or written here

  kokkos_flag = 0;
  execution_space = Host;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

FixFieldGridKokkos::~FixFieldGridKokkos()
{
}

/* ----------------------------------------------------------------------
   evaluate the per-grid-cell field on the host, then publish it to the
     device view UpdateKokkos::field_per_grid() reads
------------------------------------------------------------------------- */

void FixFieldGridKokkos::compute_field()
{
  FixFieldGrid::compute_field();

  const int nglocal = grid->nlocal;
  if (!nglocal) return;
  const int ncols = size_per_grid_cols;

  if ((int) k_array_grid.extent(0) < nglocal ||
      (int) k_array_grid.extent(1) != ncols)
    MemKK::realloc_kokkos(k_array_grid,"field/grid/kk:array_grid",nglocal,ncols);

  // array_grid is one contiguous row-major block: Memory::create(TYPE**&,n1,n2)
  //   (memory.h:114-127) does a single allocation of n1*n2 and points each
  //   row pointer into it.  tdual_float_2d_lr is LayoutRight, so the two have
  //   the same element order and the host side needs no copy at all -- wrap
  //   the existing buffer in an unmanaged View and hand it straight to
  //   deep_copy.  On a host backend the DualView's two views are the same
  //   memory and this is a no-op; on a GPU it is the one H2D transfer that
  //   has to happen either way.  The element-wise loop this replaces was
  //   pure overhead on every step.

  static_assert(std::is_same<F_FLOAT,double>::value,
                "wrapping array_grid (double**) in an F_FLOAT view assumes "
                "F_FLOAT is double; use a converting deep_copy if that changes");

  Kokkos::View<F_FLOAT**,Kokkos::LayoutRight,Kokkos::HostSpace,
               Kokkos::MemoryTraits<Kokkos::Unmanaged> >
    h_array_grid(array_grid[0],nglocal,ncols);

  // both sides can be longer than nglocal: the host array is sized to
  //   maxparticle/maxgrid, and the DualView guard above is grow-only
  //   (extent(0) < nglocal), so it keeps a high-water-mark row count after the
  //   count drops.  Kokkos::deep_copy throws on any extent mismatch
  //   (Kokkos_CopyViews.hpp:1163), so copy the leading nglocal rows rather
  //   than the whole allocation.  A leading row range of a LayoutRight view
  //   is still contiguous, so the subview costs nothing at runtime

  auto d_rows = Kokkos::subview(k_array_grid.view_device(),
                                Kokkos::make_pair(0,nglocal),Kokkos::ALL());
  Kokkos::deep_copy(d_rows,h_array_grid);
  k_array_grid.modify_device();

  d_array_grid = k_array_grid.view_device();
}

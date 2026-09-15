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
#include "error.h"
#include "grid_kokkos.h"
#include "domain.h"
#include "memory.h"

using namespace SPARTA_NS;

void GridKokkos::update_hash()
{
  typedef hash_type::size_type size_type;    // uint32_t
  typedef hash_type::key_type key_type;      // cellint
  typedef hash_type::value_type value_type;  // int
  typedef hash_type::host_mirror_type host_hash_type;

  size_type failed_count = 0;

  // Copy the keys:values from hash to Kokkos::UnorderedMap that lives on host
  if (2*(bigint) hash->size() > 4294967295LL)
    error->one(FLERR,"Grid cell hash too large for Kokkos::UnorderedMap");
  host_hash_type hash_h(2*hash->size()); // double hash capacity to prevent insertion failure
  hash_kk = hash_type(2*hash->size());
  for (volatile auto it : *hash) { // volatile keyword works around a suspected compiler bug
    key_type key = static_cast<key_type>(it.first);
    value_type val = static_cast<value_type>(it.second);
    auto insert_result = hash_h.insert(key, val);
    failed_count += insert_result.failed() ? 1 : 0;
  }
  if (failed_count) {
    error->one(FLERR, "Kokkos::UnorderedMap insertion failed");
  }

  Kokkos::deep_copy(hash_kk, hash_h);

  update_halo_index();
}

/* ----------------------------------------------------------------------
   upload Grid::halo_index, the dense cell-position -> local-index map used by
     the uniform-grid fast path, to the device
   the map itself is built by Grid::update_halo_index(), called from rehash()
     just above, so host and device agree by construction
   leaves d_halo_index empty (extent 0) when the host map is unavailable, in
     which case callers use hash_kk instead
------------------------------------------------------------------------- */

void GridKokkos::update_halo_index()
{
  d_halo_index = DAT::t_int_1d();

  if (!halo_index) return;

  const size_t nindex = (size_t) halo_nx*halo_ny*halo_nz;
  if (!nindex) return;

  HAT::t_int_1d h_halo_index(
    Kokkos::view_alloc("grid:halo_index",Kokkos::WithoutInitializing),nindex);
  for (size_t i = 0; i < nindex; i++) h_halo_index[i] = halo_index[i];

  d_halo_index = DAT::t_int_1d(
    Kokkos::view_alloc("grid:halo_index",Kokkos::WithoutInitializing),nindex);
  Kokkos::deep_copy(d_halo_index,h_halo_index);
}

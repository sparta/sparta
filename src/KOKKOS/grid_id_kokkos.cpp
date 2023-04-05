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
#include "error.h"
#include "grid_kokkos.h"

using namespace SPARTA_NS;

void GridKokkos::update_hash()
{
  typedef hash_type::size_type size_type;    // uint32_t
  typedef hash_type::key_type key_type;      // cellint
  typedef hash_type::value_type value_type;  // int
  typedef hash_type::HostMirror host_hash_type;

  size_type failed_count = 0;

  // Copy the keys:values from hash to Kokkos::UnorderedMap that lives on host
  host_hash_type hash_h(2*hash->size()); // double hash capacity to prevent insertion failure
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
}


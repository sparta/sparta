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

#ifndef SPARTA_GRID_KOKKOS_H
#define SPARTA_GRID_KOKKOS_H

#include "grid.h"
#include "kokkos_type.h"
#include <Kokkos_UnorderedMap.hpp>

namespace SPARTA_NS {

class GridKokkos : public Grid {
 public:

  // make into a view
  //ChildCell *cells;           // list of owned and ghost child cells

  // methods

  GridKokkos(class SPARTA *);
  ~GridKokkos();
  void wrap_kokkos();
  void wrap_kokkos_graphs();
  void sync(ExecutionSpace, unsigned int);
  void modify(ExecutionSpace, unsigned int);

// operations with grid cell IDs
  void update_hash();

/* ----------------------------------------------------------------------
   find child cell within iparent cell assumed to be contain pt x
   recurse down thru parents until reach a child cell
   pt can be on any boundary of parent cell
   if I don't store child cell as owned or ghost, return -1 for unknown
   else return local index of child cell
   NOTE: replace recursive with while loop
------------------------------------------------------------------------- */
  KOKKOS_INLINE_FUNCTION
  int id_find_child(int iparent, double *x) const
  {
    typedef hash_type::size_type size_type;    // uint32_t
    typedef hash_type::key_type key_type;      // cellint
    typedef hash_type::value_type value_type;  // int

    ParentCell *p = &k_pcells.d_view[iparent];
    double *lo = p->lo;
    double *hi = p->hi;
    int nx = p->nx;
    int ny = p->ny;
    int nz = p->nz;
    int ix = static_cast<int> ((x[0]-lo[0]) * nx/(hi[0]-lo[0]));
    int iy = static_cast<int> ((x[1]-lo[1]) * ny/(hi[1]-lo[1]));
    int iz = static_cast<int> ((x[2]-lo[2]) * nz/(hi[2]-lo[2]));
    if (ix == nx) ix--;
    if (iy == ny) iy--;
    if (iz == nz) iz--;

    cellint ichild = iz*nx*ny + iy*nx + ix + 1;
    cellint idchild = p->id | (ichild << p->nbits);

    size_type h_index = hash_kk.find(static_cast<key_type>(idchild));
    if (!hash_kk.valid_at(h_index)) return -1;
    int index = static_cast<int>(hash_kk.value_at(h_index));
    if (index > 0) return index-1;
    return id_find_child(-index-1,x);
  }

  // extract/return neighbor flag for iface from per-cell nmask
  // inlined for efficiency

  KOKKOS_INLINE_FUNCTION
  int neigh_decode(int nmask, int iface) const {
    return (nmask & neighmask[iface]) >> neighshift[iface];
  }

  // overwrite neighbor flag for iface in per-cell nmask
  // first line zeroes the iface bits via one's complement of mask
  // inlined for efficiency
  // return updated nmask

  KOKKOS_INLINE_FUNCTION
  int neigh_encode(int flag, int nmask, int iface) const {
    nmask &= ~neighmask[iface];
    nmask |= flag << neighshift[iface];
    return nmask;
  }

  tdual_cell_1d k_cells;
  tdual_cinfo_1d k_cinfo;
  tdual_sinfo_1d k_sinfo;
  tdual_pcell_1d k_pcells;

  Kokkos::Crs<int, DeviceType, void, int> d_csurfs;
  Kokkos::Crs<int, DeviceType, void, int> d_csplits;
  Kokkos::Crs<int, DeviceType, void, int> d_csubs;

  DAT::t_int_1d d_cellcount;
  DAT::t_int_2d d_plist;

 private:
  void grow_cells(int, int);
  void grow_pcells(int);
  void grow_sinfo(int);

  // hash for all cell IDs (owned,ghost,parent).  The _d postfix refers to the
  // fact that this hash lives on "device"

  typedef Kokkos::UnorderedMap<cellint,int> hash_type;
  hash_type hash_kk; 

};

}

#endif

/* ERROR/WARNING messages:

*/

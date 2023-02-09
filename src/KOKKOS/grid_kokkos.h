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

#ifndef SPARTA_GRID_KOKKOS_H
#define SPARTA_GRID_KOKKOS_H

#include "grid.h"
#include "kokkos_type.h"
#include <Kokkos_UnorderedMap.hpp>

namespace SPARTA_NS {

class GridKokkos : public Grid {
 public:
  typedef ArrayTypes<DeviceType> AT;

  typedef Kokkos::UnorderedMap<cellint,int> hash_type;
  typedef hash_type::size_type size_type;    // uint32_t
  typedef hash_type::key_type key_type;      // cellint
  typedef hash_type::value_type value_type;  // int

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
     compute lo/hi extent of a specific child cell within a parent cell
     plevel = level of parent
     plo/phi = parent cell corner points
     ichild ranges from 1 to Nx*Ny*Nz within parent cell
     return clo/chi corner points, caller must allocate them
  ------------------------------------------------------------------------- */

  KOKKOS_INLINE_FUNCTION
  void id_child_lohi(int plevel, double *plo, double *phi,
                     cellint ichild, double *clo, double *chi) const
  {
    int nx = k_plevels.d_view[plevel].nx;
    int ny = k_plevels.d_view[plevel].ny;
    int nz = k_plevels.d_view[plevel].nz;

    ichild--;
    int ix = ichild % nx;
    int iy = (ichild/nx) % ny;
    int iz = ichild / ((cellint) nx*ny);

    clo[0] = plo[0] + ix*(phi[0]-plo[0])/nx;
    clo[1] = plo[1] + iy*(phi[1]-plo[1])/ny;
    clo[2] = plo[2] + iz*(phi[2]-plo[2])/nz;

    chi[0] = plo[0] + (ix+1)*(phi[0]-plo[0])/nx;
    chi[1] = plo[1] + (iy+1)*(phi[1]-plo[1])/ny;
    chi[2] = plo[2] + (iz+1)*(phi[2]-plo[2])/nz;

    if (ix == nx-1) chi[0] = phi[0];
    if (iy == ny-1) chi[1] = phi[1];
    if (iz == nz-1) chi[2] = phi[2];
  }

  /* ----------------------------------------------------------------------
     find child cell within parentID which contains pt X
     level = level of parent cell
     oplo/ophi = original parent cell corner pts
     pt X can be inside or on any boundary of parent cell
     recurse from parent downward until find a child cell or reach maxlevel
     if find child cell this proc stores (owned or ghost), return its local index
     else return -1 for unknown
  ------------------------------------------------------------------------- */

  KOKKOS_INLINE_FUNCTION
  int id_find_child(cellint parentID, int plevel,
                    double *oplo, double *ophi, double *x) const
  {
    int ix,iy,iz,nx,ny,nz;
    double plo[3],phi[3],clo[3],chi[3];
    cellint childID,ichild;

    cellint id = parentID;
    int level = plevel;
    double *lo = oplo;
    double *hi = ophi;

    while (level < maxlevel) {
      nx = k_plevels.d_view[level].nx;
      ny = k_plevels.d_view[level].ny;
      nz = k_plevels.d_view[level].nz;
      ix = static_cast<int> ((x[0]-lo[0]) * nx/(hi[0]-lo[0]));
      iy = static_cast<int> ((x[1]-lo[1]) * ny/(hi[1]-lo[1]));
      iz = static_cast<int> ((x[2]-lo[2]) * nz/(hi[2]-lo[2]));
      if (ix == nx) ix--;
      if (iy == ny) iy--;
      if (iz == nz) iz--;

      ichild = (cellint) iz*nx*ny + (cellint) iy*nx + ix + 1;
      childID = (ichild << k_plevels.d_view[level].nbits) | id;

      size_type h_index = hash_kk.find(static_cast<key_type>(childID));
      if (hash_kk.valid_at(h_index)) return static_cast<int>(hash_kk.value_at(h_index));

      id = childID;
      id_child_lohi(level,lo,hi,ichild,clo,chi);
      plo[0] = clo[0]; plo[1] = clo[1]; plo[2] = clo[2];
      phi[0] = chi[0]; phi[1] = chi[1]; phi[2] = chi[2];
      lo = plo; hi = phi;
      level++;
    }

    return -1;
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
  tdual_plevel_1d k_plevels;

  Kokkos::Crs<int, DeviceType, void, int> d_csurfs;
  Kokkos::Crs<int, DeviceType, void, int> d_csplits;
  Kokkos::Crs<int, DeviceType, void, int> d_csubs;

  DAT::t_int_1d d_cellcount;
  DAT::t_int_2d d_plist;

  // hash for all cell IDs (owned,ghost,parent).  The _d postfix refers to the
  // fact that this hash lives on "device"
  hash_type hash_kk;

 private:
  void grow_cells(int, int);
  void grow_sinfo(int);
  void grow_pcells();
};

}

#endif

/* ERROR/WARNING messages:

*/

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

#ifndef SPARTA_COMM_KOKKOS_H
#define SPARTA_COMM_KOKKOS_H

#include "comm.h"
#include "grid.h"
#include "particle.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

template<int NEED_ATOMICS>
struct TagCommMigrateParticles{};

template<int NEED_ATOMICS>
struct TagCommMigrateUnpackParticles{};

class CommKokkos : public Comm {
 public:
  typedef ArrayTypes<DeviceType> AT;

  CommKokkos(class SPARTA *);
  ~CommKokkos();
  int migrate_particles(int, int*, DAT::t_int_1d);
  void migrate_cells(int);

  template<int NEED_ATOMICS>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagCommMigrateParticles<NEED_ATOMICS>, const int&) const;

  template<int NEED_ATOMICS>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagCommMigrateUnpackParticles<NEED_ATOMICS>, const int&) const;

  inline
  void pack_serial(const int, const int) const;

 private:
  DAT::tdual_int_scalar k_nsend;
  DAT::t_int_scalar d_nsend;
  HAT::t_int_scalar h_nsend;
  DAT::t_int_scalar d_nlocal;

  typedef Kokkos::
    DualView<Grid::ChildCell*, Kokkos::LayoutRight, DeviceType> tdual_cell_1d;
  typedef tdual_cell_1d::t_dev t_cell_1d;
  t_cell_1d d_cells;

  typedef Kokkos::
    DualView<Particle::OnePart*, Kokkos::LayoutRight, DeviceType> tdual_particle_1d;
  typedef tdual_particle_1d::t_dev t_particle_1d;
  t_particle_1d d_particles;

  DAT::t_int_1d d_plist;
  DAT::tdual_int_1d k_pproc;
  DAT::t_int_1d d_pproc;
  DAT::t_char_1d d_sbuf;
  DAT::t_char_1d d_rbuf;

  int nbytes_particle;
};

}

#endif

/* ERROR/WARNING messages:

*/

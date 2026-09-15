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

#ifdef COMPUTE_CLASS

ComputeStyle(gas/collision/grid/kk,ComputeGasCollisionGridKokkos)

#else

#ifndef SPARTA_COMPUTE_GAS_COLLISION_GRID_KOKKOS_H
#define SPARTA_COMPUTE_GAS_COLLISION_GRID_KOKKOS_H

#include "compute_gas_collision_grid.h"
#include "kokkos_base.h"
#include "kokkos_type.h"
#include "particle.h"

namespace SPARTA_NS {

class ComputeGasCollisionGridKokkos : public ComputeGasCollisionGrid, public KokkosBase {
 public:
  ComputeGasCollisionGridKokkos(class SPARTA *, int, char **);
  ComputeGasCollisionGridKokkos(class SPARTA *);
  ~ComputeGasCollisionGridKokkos();
  void compute_per_grid_kokkos() {}   // tallying happens in Collide, not here
  void clear();
  void pre_gas_tally();
  void post_gas_tally();
  void reallocate();

  // tally a single gas collision in icell on device
  // reaction = 0 for a collision that did not induce a reaction
  // this compute only tallies non-reacting collisions
  // Collide parallelizes over grid cells (one icell per thread), so the
  //   per-cell tally has no write contention and needs no atomics/duplication
  // ATOMIC_REDUCTION template arg is unused, kept for a uniform call interface

  template<int ATOMIC_REDUCTION>
  KOKKOS_INLINE_FUNCTION
  void gas_tally_kk(int icell, int reaction,
                    Particle::OnePart *iorig, Particle::OnePart *jorig,
                    Particle::OnePart * /*ip*/, Particle::OnePart * /*jp*/,
                    Particle::OnePart * /*kp*/) const
  {
    // skip if a reaction (reactions tallied by compute gas/reaction/grid)

    if (reaction) return;

    // skip if icell not in grid group

    if (!(d_cinfo[icell].mask & groupbit)) return;

    // skip if either particle species not in mixture group

    int igroup = d_s2g(imix,iorig->ispecies);
    int jgroup = d_s2g(imix,jorig->ispecies);
    if (igroup < 0 || jgroup < 0) return;

    // tally the collision to its grid cell

    d_vector_grid(icell) += 1.0;
  }

 private:
  DAT::tdual_float_1d k_vector_grid;
  // d_vector_grid is inherited from KokkosBase (read by fix ave/grid/kk)

  t_cinfo_1d d_cinfo;
  DAT::t_int_2d d_s2g;
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/

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

ComputeStyle(react/boundary/kk,ComputeReactBoundaryKokkos)

#else

#ifndef SPARTA_REACT_BOUNDARY_KOKKOS_H
#define SPARTA_REACT_BOUNDARY_KOKKOS_H

#include "compute_react_boundary.h"
#include "kokkos_base.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

class ComputeReactBoundaryKokkos : public ComputeReactBoundary, public KokkosBase {
 public:
  ComputeReactBoundaryKokkos(class SPARTA *, int, char **);
  ComputeReactBoundaryKokkos(class SPARTA *);
  ~ComputeReactBoundaryKokkos() override;

  void init() override;
  void compute_array() override;
  void clear() override;

  // called by UpdateKokkos around the move kernel

  void pre_boundary_tally();
  void post_boundary_tally();

  /* ----------------------------------------------------------------------
     tally a surface reaction on box face iface
     mirrors ComputeReactBoundary::boundary_tally()
       (compute_react_boundary.cpp:139-166) exactly
     norm is unused here; it is in the signature so UpdateKokkos can dispatch
       every boundary tally compute through one call
  ------------------------------------------------------------------------- */

  template<int ATOMIC_REDUCTION>
  KOKKOS_INLINE_FUNCTION
  void boundary_tally_kk(double /*dtremain*/, int iface, int /*istyle*/,
                         int reaction, Particle::OnePart * /*iorig*/,
                         Particle::OnePart * /*ip*/, Particle::OnePart * /*jp*/,
                         const double * /*norm*/) const
  {
    // skip if no reaction

    if (reaction == 0) return;
    reaction--;

    // skip if this face's reaction model is not a match

    if (d_surf_react[iface] != isr) return;

    auto v_myarray = ScatterViewHelper<typename NeedDup<ATOMIC_REDUCTION,DeviceType>::value,decltype(dup_myarray),decltype(ndup_myarray)>::get(dup_myarray,ndup_myarray);
    auto a_myarray = v_myarray.template access<typename AtomicDup<ATOMIC_REDUCTION,DeviceType>::value>();

    // for rpflag, tally each column whose reaction2col entry is set
    // for rpflag = 0, tally the reaction directly

    if (rpflag) {
      for (int i = 0; i < ntotal; i++)
        if (d_reaction2col(reaction,i)) a_myarray(iface,i) += 1.0;
    } else a_myarray(iface,reaction) += 1.0;
  }

 private:
  DAT::tdual_float_2d_lr k_myarray;     // local accumulator array
  DAT::t_float_2d_lr d_myarray;

  DAT::t_int_1d d_surf_react;           // per box face: its surf react model
  DAT::t_int_2d d_reaction2col;         // 1 if ireaction tallies into icol

  int need_dup;
  Kokkos::Experimental::ScatterView<F_FLOAT**, typename DAT::t_float_2d_lr::array_layout,DeviceType,typename Kokkos::Experimental::ScatterSum,typename Kokkos::Experimental::ScatterDuplicated> dup_myarray;
  Kokkos::Experimental::ScatterView<F_FLOAT**, typename DAT::t_float_2d_lr::array_layout,DeviceType,typename Kokkos::Experimental::ScatterSum,typename Kokkos::Experimental::ScatterNonDuplicated> ndup_myarray;
};

}

#endif
#endif

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

ComputeStyle(react/surf/kk,ComputeReactSurfKokkos)

#else

#ifndef SPARTA_COMPUTE_REACT_SURF_KOKKOS_H
#define SPARTA_COMPUTE_REACT_SURF_KOKKOS_H

#include "compute_react_surf.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

class ComputeReactSurfKokkos : public ComputeReactSurf {
 public:
  ComputeReactSurfKokkos(class SPARTA *, int, char **);
  ComputeReactSurfKokkos(class SPARTA *);
  ~ComputeReactSurfKokkos();
  void init();
  void clear();
  int tallyinfo(surfint *&);
  void post_process_surf();
  void pre_surf_tally();
  void post_surf_tally();

/* ----------------------------------------------------------------------
   tally a surface reaction for particle colliding with surf element isurf
   mirrors ComputeReactSurf::surf_tally(); per-surf tally compressed to host
------------------------------------------------------------------------- */

  template <int ATOMIC_REDUCTION>
  KOKKOS_INLINE_FUNCTION
  void surf_tally_kk(double /*dtremain*/, int isurf, int /*icell*/, int reaction,
                     Particle::OnePart * /*iorig*/,
                     Particle::OnePart * /*ip*/, Particle::OnePart * /*jp*/) const
  {
    if (reaction == 0) return;
    reaction--;

    surfint surfID;
    if (dim == 2) {
      if (!(d_lines[isurf].mask & groupbit)) return;
      if (d_lines[isurf].isr != isr) return;
      surfID = d_lines[isurf].id;
    } else {
      if (!(d_tris[isurf].mask & groupbit)) return;
      if (d_tris[isurf].isr != isr) return;
      surfID = d_tris[isurf].id;
    }

    int itally = isurf;
    d_tally2surf(itally) = surfID;
    d_surf2tally(isurf) = isurf;

    auto v_array_surf_tally = ScatterViewHelper<typename NeedDup<ATOMIC_REDUCTION,DeviceType>::value,decltype(dup_array_surf_tally),decltype(ndup_array_surf_tally)>::get(dup_array_surf_tally,ndup_array_surf_tally);
    auto a_array_surf_tally = v_array_surf_tally.template access<typename AtomicDup<ATOMIC_REDUCTION,DeviceType>::value>();

    if (rpflag) {
      for (int i = 0; i < ntotal; i++)
        if (d_reaction2col(reaction,i)) a_array_surf_tally(itally,i) += 1.0;
    } else a_array_surf_tally(itally,reaction) += 1.0;
  }

 private:
  DAT::t_int_2d d_reaction2col;

  DAT::tdual_float_2d_lr k_array_surf_tally;
  DAT::t_float_2d_lr d_array_surf_tally;

  int need_dup;
  Kokkos::Experimental::ScatterView<F_FLOAT**, typename DAT::t_float_2d_lr::array_layout,DeviceType,typename Kokkos::Experimental::ScatterSum,typename Kokkos::Experimental::ScatterDuplicated> dup_array_surf_tally;
  Kokkos::Experimental::ScatterView<F_FLOAT**, typename DAT::t_float_2d_lr::array_layout,DeviceType,typename Kokkos::Experimental::ScatterSum,typename Kokkos::Experimental::ScatterNonDuplicated> ndup_array_surf_tally;

  DAT::t_surfint_1d d_tally2surf;
  DAT::tdual_surfint_1d k_tally2surf;
  DAT::t_int_1d d_surf2tally;

  t_line_1d d_lines;
  t_tri_1d d_tris;

  int nsurf_tally_alloc;

  void grow_tally();
  void resize_device(int);
};

}

#endif
#endif

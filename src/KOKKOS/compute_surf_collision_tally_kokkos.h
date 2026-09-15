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

ComputeStyle(surf/collision/tally/kk,ComputeSurfCollisionTallyKokkos)

#else

#ifndef SPARTA_COMPUTE_SURF_COLLISION_TALLY_KOKKOS_H
#define SPARTA_COMPUTE_SURF_COLLISION_TALLY_KOKKOS_H

#include "compute_surf_collision_tally.h"
#include "kokkos_base.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

// unlike the per-grid and per-surf tally computes, this one appends one row
//   per collision event, so the row count is not known before the move
//   kernel runs.  rows are claimed with an atomic counter; if a claim lands
//   past the end of the buffer the row is dropped and d_tally_overflow is
//   raised, which makes UpdateKokkos grow every such compute and repeat the
//   move.  truncating instead would silently corrupt dump tally output.

class ComputeSurfCollisionTallyKokkos : public ComputeSurfCollisionTally, public KokkosBase {
 public:
  ComputeSurfCollisionTallyKokkos(class SPARTA *, int, char **);
  ComputeSurfCollisionTallyKokkos(class SPARTA *);
  ~ComputeSurfCollisionTallyKokkos() override;

  void clear() override;

  // called by UpdateKokkos around the move kernel

  void pre_surf_tally();
  void post_surf_tally();
  void grow_tally_kokkos(int);

  // grow to what the overflowed attempt actually needed; the device counter
  //   kept counting past the end of the buffer, so it is that number

  void grow_after_overflow()
  {
    Kokkos::deep_copy(h_ntally,d_ntally);
    grow_tally_kokkos(h_ntally());
  }

  // every attempt of the retry loop re-runs the whole pass, so the rows the
  //   aborted attempt appended have to be taken back before the next one.
  //   Not by zeroing: the move kernel runs once per migration iteration and
  //   the tally accumulates across all of them, so an attempt rewinds to
  //   where its own pass started, which mark_ntally() records.

  void mark_ntally()
  {
    Kokkos::deep_copy(h_ntally,d_ntally);
    ntally_mark = h_ntally();
  }

  void rewind_ntally() { Kokkos::deep_copy(d_ntally,ntally_mark); }

  DAT::t_int_scalar d_overflow;         // set by UpdateKokkos each step

  KOKKOS_INLINE_FUNCTION
  void surf_tally_kk(double dtremain, int isurf, int icell, int reaction,
                     Particle::OnePart *iorig,
                     Particle::OnePart *ip, Particle::OnePart *jp) const
  {
    // this compute tallies only collisions that induce no reaction;
    //   reactions belong to compute surf/reaction/tally

    if (reaction) return;

    if (dim == 2) {
      if (!(d_lines(isurf).mask & groupbit)) return;
    } else {
      if (!(d_tris(isurf).mask & groupbit)) return;
    }

    const int origspecies = iorig->ispecies;
    if (d_s2g(imix,origspecies) < 0) return;

    const int itally = Kokkos::atomic_fetch_add(&d_ntally(),1);
    if (itally >= (int) d_array_tally.extent(0)) {
      d_overflow() = 1;
      return;
    }

    for (int m = 0; m < nvalue; m++) {
      switch (d_which[m]) {
      case IDSURF:
        if (dim == 2) d_array_tally(itally,m) = d_ubuf(d_lines(isurf).id).d;
        else d_array_tally(itally,m) = d_ubuf(d_tris(isurf).id).d;
        break;
      case ID:     d_array_tally(itally,m) = d_ubuf(ip->id).d; break;
      case TYPE:   d_array_tally(itally,m) = d_ubuf(ip->ispecies+1).d; break;
      case XC:     d_array_tally(itally,m) = iorig->x[0]; break;
      case YC:     d_array_tally(itally,m) = iorig->x[1]; break;
      case ZC:     d_array_tally(itally,m) = iorig->x[2]; break;
      case TIME:   d_array_tally(itally,m) = dt - dtremain; break;
      case VXPRE:  d_array_tally(itally,m) = iorig->v[0]; break;
      case VYPRE:  d_array_tally(itally,m) = iorig->v[1]; break;
      case VZPRE:  d_array_tally(itally,m) = iorig->v[2]; break;
      case VXPOST: d_array_tally(itally,m) = ip->v[0]; break;
      case VYPOST: d_array_tally(itally,m) = ip->v[1]; break;
      case VZPOST: d_array_tally(itally,m) = ip->v[2]; break;
      }
    }
  }

 private:
  // must match the enum in compute_surf_collision_tally.cpp

  enum{IDSURF,ID,TYPE,TIME,XC,YC,ZC,VXPRE,VYPRE,VZPRE,VXPOST,VYPOST,VZPOST};

  int maxtally_host;                 // rows array_tally is allocated for
  DAT::tdual_float_2d_lr k_array_tally;
  DAT::t_float_2d_lr d_array_tally;
  DAT::t_int_scalar d_ntally;
  int ntally_mark;
  HAT::t_int_scalar h_ntally;

  DAT::t_int_1d d_which;
  DAT::t_int_2d d_s2g;

  t_line_1d d_lines;
  t_tri_1d d_tris;

  double dt;
};

}

#endif
#endif

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

ComputeStyle(gas/reaction/tally/kk,ComputeGasReactionTallyKokkos)

#else

#ifndef SPARTA_COMPUTE_GAS_REACTION_TALLY_KOKKOS_H
#define SPARTA_COMPUTE_GAS_REACTION_TALLY_KOKKOS_H

#include "compute_gas_reaction_tally.h"
#include "kokkos_base.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

// unlike the per-grid and per-surf tally computes, this one appends one row
//   per collision event, so the row count is not known before the collide
//   kernel runs.  rows are claimed with an atomic counter; if a claim lands
//   past the end of the buffer the row is dropped and d_tally_overflow is
//   raised, which makes CollideVSSKokkos grow every such compute and repeat the
//   collision pass.  truncating instead would silently corrupt dump tally output.

class ComputeGasReactionTallyKokkos : public ComputeGasReactionTally, public KokkosBase {
 public:
  ComputeGasReactionTallyKokkos(class SPARTA *, int, char **);
  ComputeGasReactionTallyKokkos(class SPARTA *);
  ~ComputeGasReactionTallyKokkos() override;

  void clear() override;

  // called by CollideVSSKokkos around the collision kernel

  void pre_gas_tally();
  void post_gas_tally();
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

  DAT::t_int_scalar d_overflow;         // set by CollideVSSKokkos each step

  template<int ATOMIC_REDUCTION>
  KOKKOS_INLINE_FUNCTION
  void gas_tally_kk(int icell, int reaction,
                    Particle::OnePart *iorig, Particle::OnePart *jorig,
                    Particle::OnePart *ip, Particle::OnePart *jp,
                    Particle::OnePart *kp) const
  {
    // this compute tallies only collisions that induce a reaction;
    //   plain collisions belong to compute gas/collision/tally

    if (!reaction) return;

    if (!(d_cinfo[icell].mask & groupbit)) return;
    if (d_s2g(imix,iorig->ispecies) < 0) return;
    if (d_s2g(imix,jorig->ispecies) < 0) return;

    const int itally = Kokkos::atomic_fetch_add(&d_ntally(),1);
    if (itally >= (int) d_array_tally.extent(0)) {
      d_overflow() = 1;
      return;
    }

    for (int m = 0; m < nvalue; m++) {
      switch (d_which[m]) {
      case REACTION:   d_array_tally(itally,m) = d_ubuf(reaction).d; break;

      case IDCELL:     d_array_tally(itally,m) = d_ubuf(d_cells[icell].id).d; break;
      case ID1PRE:     d_array_tally(itally,m) = d_ubuf(iorig->id).d; break;
      case ID2PRE:     d_array_tally(itally,m) = d_ubuf(jorig->id).d; break;
      case ID1POST:    d_array_tally(itally,m) = d_ubuf(ip->id).d; break;
      case ID2POST:
        d_array_tally(itally,m) = (jp == NULL) ? d_ubuf(0).d : d_ubuf(jp->id).d; break;
      case ID3POST:
        d_array_tally(itally,m) = (kp == NULL) ? d_ubuf(0).d : d_ubuf(kp->id).d; break;

      case TYPE1PRE:   d_array_tally(itally,m) = d_ubuf(iorig->ispecies+1).d; break;
      case TYPE2PRE:   d_array_tally(itally,m) = d_ubuf(jorig->ispecies+1).d; break;
      case TYPE1POST:  d_array_tally(itally,m) = d_ubuf(ip->ispecies+1).d; break;
      case TYPE2POST:
        d_array_tally(itally,m) = (jp == NULL) ? d_ubuf(0).d : d_ubuf(jp->ispecies+1).d; break;
      case TYPE3POST:
        d_array_tally(itally,m) = (kp == NULL) ? d_ubuf(0).d : d_ubuf(kp->ispecies+1).d; break;

      case VX1PRE:     d_array_tally(itally,m) = iorig->v[0]; break;
      case VY1PRE:     d_array_tally(itally,m) = iorig->v[1]; break;
      case VZ1PRE:     d_array_tally(itally,m) = iorig->v[2]; break;
      case VX2PRE:     d_array_tally(itally,m) = jorig->v[0]; break;
      case VY2PRE:     d_array_tally(itally,m) = jorig->v[1]; break;
      case VZ2PRE:     d_array_tally(itally,m) = jorig->v[2]; break;

      case VX1POST:    d_array_tally(itally,m) = ip->v[0]; break;
      case VY1POST:    d_array_tally(itally,m) = ip->v[1]; break;
      case VZ1POST:    d_array_tally(itally,m) = ip->v[2]; break;
      case VX2POST:    d_array_tally(itally,m) = (jp == NULL) ? 0.0 : jp->v[0]; break;
      case VY2POST:    d_array_tally(itally,m) = (jp == NULL) ? 0.0 : jp->v[1]; break;
      case VZ2POST:    d_array_tally(itally,m) = (jp == NULL) ? 0.0 : jp->v[2]; break;
      case VX3POST:    d_array_tally(itally,m) = (kp == NULL) ? 0.0 : kp->v[0]; break;
      case VY3POST:    d_array_tally(itally,m) = (kp == NULL) ? 0.0 : kp->v[1]; break;
      case VZ3POST:    d_array_tally(itally,m) = (kp == NULL) ? 0.0 : kp->v[2]; break;
      }
    }
  }

 private:
  // must match the enum in compute_gas_reaction_tally.cpp

  enum{REACTION,IDCELL,ID1PRE,ID2PRE,ID1POST,ID2POST,ID3POST,TYPE1PRE,TYPE2PRE,
       TYPE1POST,TYPE2POST,TYPE3POST,VX1PRE,VY1PRE,VZ1PRE,VX2PRE,VY2PRE,VZ2PRE,
       VX1POST,VY1POST,VZ1POST,VX2POST,VY2POST,VZ2POST,VX3POST,VY3POST,VZ3POST};

  int maxtally_host;                 // rows array_tally is allocated for
  DAT::tdual_float_2d_lr k_array_tally;
  DAT::t_float_2d_lr d_array_tally;
  DAT::t_int_scalar d_ntally;
  int ntally_mark;
  HAT::t_int_scalar h_ntally;

  DAT::t_int_1d d_which;
  DAT::t_int_2d d_s2g;

  t_cell_1d d_cells;
  t_cinfo_1d d_cinfo;
};

}

#endif
#endif

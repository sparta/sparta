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

#ifdef REGION_CLASS

RegionStyle(intersect/kk,RegIntersectKokkos)

#else

#ifndef SPARTA_REGION_INTERSECT_KOKKOS_H
#define SPARTA_REGION_INTERSECT_KOKKOS_H

#include "region_intersect.h"

#include "kokkos_base.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

// a composite region cannot dispatch to its sub-regions on the device, so it
//   flattens them into a postfix (RPN) token stream instead; see
//   region_prim_kokkos.h.  each sub-region contributes its own whole
//   sub-stream, so a sub-region may itself be a region union or region
//   intersect: composites nest to arbitrary depth, bounded only by the
//   RKK_MAX_STACK boolean stack depth checked here at flatten time.

class RegIntersectKokkos : public RegIntersect, public KokkosBase {

 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  RegIntersectKokkos(class SPARTA *, int, char **);
  ~RegIntersectKokkos() override;

  void match_all_kokkos(DAT::tdual_int_1d) override;
  int flatten_region_kokkos(tdual_region_token_1d &) override;

 private:
  int groupbit;
  typename AT::t_int_1d d_match;
  t_particle_1d d_particles;

  tdual_region_token_1d k_tokens;   // this region's own token stream
  int ntoken;
};

}

#endif
#endif

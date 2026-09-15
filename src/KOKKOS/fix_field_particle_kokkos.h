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

#ifdef FIX_CLASS

FixStyle(field/particle/kk,FixFieldParticleKokkos)

#else

#ifndef SPARTA_FIX_FIELD_PARTICLE_KOKKOS_H
#define SPARTA_FIX_FIELD_PARTICLE_KOKKOS_H

#include "fix_field_particle.h"
#include "kokkos_base.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

// the particle-style variables this fix evaluates are host-only (see
//   VariableKokkos), so the evaluation itself stays on the host and only the
//   result is published to the device.  unlike field/grid this runs every
//   timestep (update.cpp:672), so it carries a per-step host round trip --
//   correctness first; a device variable evaluator would remove it

class FixFieldParticleKokkos : public FixFieldParticle, public KokkosBase {
 public:
  FixFieldParticleKokkos(class SPARTA *, int, char **);
  ~FixFieldParticleKokkos() override;
  void compute_field() override;

 private:
  DAT::tdual_float_2d_lr k_array_particle;
};

}

#endif
#endif

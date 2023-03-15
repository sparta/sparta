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

#ifdef COMPUTE_CLASS

ComputeStyle(ke/particle/kk,ComputeKEParticleKokkos)

#else

#ifndef SPARTA_COMPUTE_KE_PARTICLE_KOKKOS_H
#define SPARTA_COMPUTE_KE_PARTICLE_KOKKOS_H

#include "compute_ke_particle.h"
#include "kokkos_type.h"
#include "kokkos_base.h"

namespace SPARTA_NS {

class ComputeKEParticleKokkos : public ComputeKEParticle, public KokkosBase {
 public:
  ComputeKEParticleKokkos(class SPARTA *, int, char **);
  ~ComputeKEParticleKokkos();
  void compute_per_particle();
  void compute_per_particle_kokkos();

  KOKKOS_INLINE_FUNCTION
  void operator()(const int&) const;

  DAT::tdual_float_1d k_vector_particle;

 private:
  t_particle_1d d_particles;
  t_species_1d d_species;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

W: More than one compute ke/particle

This may be inefficient since each such compute stores a vector
of length equal to the number of particles.

*/

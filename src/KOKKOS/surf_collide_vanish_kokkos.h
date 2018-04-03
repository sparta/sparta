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

#ifdef SURF_COLLIDE_CLASS

SurfCollideStyle(vanish/kk,SurfCollideVanishKokkos)

#else

#ifndef SPARTA_SURF_COLLIDE_VANISH_KOKKOS_H
#define SPARTA_SURF_COLLIDE_VANISH_KOKKOS_H

#include "surf_collide_vanish.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

class SurfCollideVanishKokkos : public SurfCollideVanish {
 public:
  typedef ArrayTypes<DeviceType> AT;

  SurfCollideVanishKokkos(class SPARTA *, int, char **);
  SurfCollideVanishKokkos(class SPARTA *);
  ~SurfCollideVanishKokkos() {}
  Particle::OnePart *collide(Particle::OnePart *&, double *, double &, int) { return NULL; }

  DAT::tdual_int_scalar k_nsingle;
  typename AT::t_int_scalar d_nsingle;
  HAT::t_int_scalar h_nsingle;

  /* ----------------------------------------------------------------------
     particle collision with surface with optional chemistry
     ip = particle with current x = collision pt, current v = incident v
     norm = surface normal unit vector
     simply return ip = NULL to delete particle
  ------------------------------------------------------------------------- */
  
  KOKKOS_INLINE_FUNCTION
  Particle::OnePart*
  collide_kokkos(Particle::OnePart *&ip, const double *, double &, int) const
  {
    Kokkos::atomic_fetch_add(&d_nsingle(),1);

    ip = NULL;
    return NULL;
  }
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

*/

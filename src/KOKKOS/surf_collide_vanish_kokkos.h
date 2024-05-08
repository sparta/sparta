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

  SurfCollideVanishKokkos(class SPARTA *, int, char **);
  SurfCollideVanishKokkos(class SPARTA *);
  ~SurfCollideVanishKokkos() {}
  void pre_collide();
  void post_collide();

 private:

  DAT::tdual_int_scalar k_nsingle;
  DAT::t_int_scalar d_nsingle;
  HAT::t_int_scalar h_nsingle;

 public:

  /* ----------------------------------------------------------------------
     particle collision with surface with optional chemistry
     ip = particle with current x = collision pt, current v = incident v
     norm = surface normal unit vector
     simply return ip = NULL to delete particle
     return reaction = 0 = no reaction took place
  ------------------------------------------------------------------------- */

  template<int REACT, int ATOMIC_REDUCTION>
  KOKKOS_INLINE_FUNCTION
  Particle::OnePart* collide_kokkos(Particle::OnePart *&ip, double &,
                                    int, const double *, int, int &,
                                    const DAT::t_int_scalar &, const DAT::t_int_scalar &) const
  {
    if (ATOMIC_REDUCTION == 0)
      d_nsingle()++;
    else
      Kokkos::atomic_increment(&d_nsingle());

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

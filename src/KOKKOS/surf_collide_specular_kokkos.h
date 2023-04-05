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

SurfCollideStyle(specular/kk,SurfCollideSpecularKokkos)

#else

#ifndef SPARTA_SURF_COLLIDE_SPECULAR_KOKKOS_H
#define SPARTA_SURF_COLLIDE_SPECULAR_KOKKOS_H

#include "surf_collide_specular.h"
#include "kokkos_type.h"
#include "math_extra_kokkos.h"

namespace SPARTA_NS {

class SurfCollideSpecularKokkos : public SurfCollideSpecular {
 public:

  SurfCollideSpecularKokkos(class SPARTA *, int, char **);
  SurfCollideSpecularKokkos(class SPARTA *);
  ~SurfCollideSpecularKokkos() {}

  /* ----------------------------------------------------------------------
     particle collision with surface with optional chemistry
     ip = particle with current x = collision pt, current v = incident v
     isurf = index of surface element
     norm = surface normal unit vector
     isr = index of reaction model if >= 0, -1 for no chemistry
     ip = set to NULL if destroyed by chemistry
     return jp = new particle if created by chemistry
     return reaction = index of reaction (1 to N) that took place, 0 = no reaction
     resets particle(s) to post-collision outward velocity
  ------------------------------------------------------------------------- */

  KOKKOS_INLINE_FUNCTION
  Particle::OnePart* collide_kokkos(Particle::OnePart *&ip, double &,
                                    int, const double *norm, int, int &) const
  {
    Kokkos::atomic_increment(&d_nsingle());

    // if surface chemistry defined, attempt reaction
    // reaction > 0 if reaction took place

    //Particle::OnePart iorig;
    Particle::OnePart *jp = NULL;
    //reaction = 0;
    int velreset = 0;

    //if (isr >= 0) {
    //  if (modify->n_surf_react) memcpy(&iorig,ip,sizeof(Particle::OnePart));
    //  reaction = surf->sr[isr]->react(ip,isurf,norm,jp,velreset);
    //  if (reaction) surf->nreact_one++;
    //}

    // specular or noslip reflection for each particle
    // only if SurfReact did not already reset velocities
    // also both partiticles need to trigger any fixes
    //   to update per-particle properties which depend on
    //   temperature of the particle, e.g. fix vibmode and fix ambipolar
    // NOTE: not doing this for this specular model,
    //   since temperature does not change, would need to add a twall arg

    if (ip) {
      if (!velreset) {
        if (noslip_flag) {

          // noslip reflection
          // reflect incident v, all three components.

          MathExtraKokkos::negate3(ip->v);
        } else {

          // specular reflection

          MathExtraKokkos::reflect3(ip->v,norm);
        }
      }

      //if (modify->n_update_custom) {
      //  int i = ip - particle->particles;
      //  modify->update_custom(i,twall,twall,twall,vstream);
      //}
    }

    //if (jp) {
    //  if (!velreset) {
    //    if (noslip_flag) {

          // noslip reflection
          // reflect incident v, all three components.

    //      MathExtra::negate3(jp->v);
    //    } else {

          // specular reflection

    //      MathExtra::reflect3(jp->v,norm);
    //    }
    //  }

    //  if (modify->n_update_custom) {
    //    int j = jp - particle->particles;
    //    modify->update_custom(j,twall,twall,twall,vstream);
    //}
    //}

    // call any fixes with a surf_react() method
    // they may reset j to -1, e.g. fix ambipolar
    //   in which case newly created j is deleted

    //if (reaction && modify->n_surf_react) {
    //  int i = -1;
    //  if (ip) i = ip - particle->particles;
    //  int j = -1;
    //  if (jp) j = jp - particle->particles;
    //  modify->surf_react(&iorig,i,j);
    //  if (jp && j < 0) {
    //  jp = NULL;
    //  particle->nlocal--;
    //  }
    //}

    return jp;
  };

  DAT::tdual_int_scalar k_nsingle;
  DAT::t_int_scalar d_nsingle;
  HAT::t_int_scalar h_nsingle;
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

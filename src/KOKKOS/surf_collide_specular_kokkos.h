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
  typedef ArrayTypes<DeviceType> AT;

  SurfCollideSpecularKokkos(class SPARTA *, int, char **);
  SurfCollideSpecularKokkos(class SPARTA *);
  ~SurfCollideSpecularKokkos() {}
  Particle::OnePart *collide(Particle::OnePart *&, double *, double &, int) { return NULL; }

  /* ----------------------------------------------------------------------
     particle collision with surface with optional chemistry
     ip = particle with current x = collision pt, current v = incident v
     norm = surface normal unit vector
     isr = index of reaction model if >= 0, -1 for no chemistry
     ip = set to NULL if destroyed by chemistry
     return jp = new particle if created by chemistry
     return reaction = index of reaction (1 to N) that took place, 0 = no reaction
     resets particle(s) to post-collision outward velocity
  ------------------------------------------------------------------------- */
  
  KOKKOS_INLINE_FUNCTION
  Particle::OnePart* collide_kokkos(Particle::OnePart *&ip, const double *norm, double &, int, int &) const
  {
    Kokkos::atomic_fetch_add(&d_nsingle(),1);
  
    // if surface chemistry defined, attempt reaction
    // reaction > 0 if reaction took place
  
    //Particle::OnePart iorig;
    Particle::OnePart *jp = NULL;
    //reaction = 0;
  
    //if (isr >= 0) {
    //  if (modify->n_surf_react) memcpy(&iorig,ip,sizeof(Particle::OnePart));
    //  reaction = surf->sr[isr]->react(ip,norm,jp);
    //  if (reaction) surf->nreact_one++;
    //}
  
    // specular reflection for each particle
    // reflect incident v around norm
  
    if (ip) MathExtraKokkos::reflect3(ip->v,norm);
    //if (jp) MathExtra::reflect3(jp->v,norm);
  
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
    //    jp = NULL;
    //    particle->nlocal--;
    //  }
    //}
      
    return jp;
  };

  DAT::tdual_int_scalar k_nsingle;
  typename AT::t_int_scalar d_nsingle;
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

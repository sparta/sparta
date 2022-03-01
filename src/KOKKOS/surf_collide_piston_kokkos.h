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

SurfCollideStyle(piston/kk,SurfCollidePistonKokkos)

#else

#ifndef SPARTA_SURF_COLLIDE_PISTON_KOKKOS_H
#define SPARTA_SURF_COLLIDE_PISTON_KOKKOS_H

#include "surf_collide_piston.h"
#include "kokkos_type.h"
#include "error.h"

namespace SPARTA_NS {

class SurfCollidePistonKokkos : public SurfCollidePiston {
 public:

  SurfCollidePistonKokkos(class SPARTA *, int, char **);
  SurfCollidePistonKokkos(class SPARTA *);
  ~SurfCollidePistonKokkos() {}

  void init();

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
  Particle::OnePart* collide_kokkos(Particle::OnePart *&ip, const double *norm, double &dtremain, int, int &) const
  {
    Kokkos::atomic_increment(&d_nsingle());

    // if surface chemistry defined, attempt reaction
    // reaction > 0 if reaction took place

    //  Particle::OnePart iorig;
    Particle::OnePart *jp = NULL;
    // reaction = 0;

    //    if (isr >= 0) {
    //      if (modify->n_surf_react) memcpy(&iorig,ip,sizeof(Particle::OnePart));
    //      reaction = surf->sr[isr]->react(ip,norm,jp);
    //      if (reaction) surf->nreact_one++;
    //    }

    // norm will be in single coordinate direction
    // dir = 0,1,2 for wall (or surface) with norm parallel to x,y,z
    // which = 0/1 for wall (or surface) with +/- normal (lo/hi wall)

    int dim,which;

    if (norm[0] != 0.0) {
      dim = 0;
      if (norm[0] < 0.0) which = 1;
      else which = 0;
    } else if (norm[1] != 0.0) {
      dim = 1;
      if (norm[1] < 0.0) which = 1;
      else which = 0;
    } else {
      dim = 2;
      if (norm[2] < 0.0) which = 1;
      else which = 0;
    }

    double *x = ip->x;
    double *v = ip->v;
    double xwall = x[dim];
    double xorig = xwall - v[dim]*(dt - dtremain);
    double vorig = v[dim];

    // piston reflection: see eqs 12.30 and 12.31 in Bird 1994, p 288
    // uprime = post-collision velocity component
    // xprime = post-collision coordinate component
    // delete particle and return if xprime is not inside box
    // formula for dtremain works for both which = 0/1
    //   since numerator and denominator are always same sign

    double uprime,xprime;

    if (which == 0) {
      uprime = -2.0*vwall - vorig;
      xprime = 2.0*xwall - xorig + uprime*dt;
      if (xprime <= xwall) {
        ip = NULL;
        return NULL;
      }
    } else {
      uprime = 2.0*vwall - vorig;
      xprime = 2.0*xwall - xorig + uprime*dt;
      if (xprime >= xwall) {
        ip = NULL;
        return NULL;
      }
    }

    dtremain = (xprime - xwall) / uprime;
    v[dim] = uprime;

    // call any fixes with a surf_react() method
    // they may reset j to -1, e.g. fix ambipolar
    //   in which case newly created j is deleted

    //    if (reaction && modify->n_surf_react) {
    //      int i = -1;
    //      if (ip) i = ip - particle->particles;
    //      int j = -1;
    //      if (jp) j = jp - particle->particles;
    //      modify->surf_react(&iorig,i,j);
    //      if (jp && j < 0) {
    //        jp = NULL;
    //        particle->nlocal--;
    //      }
    //    }

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

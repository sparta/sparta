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

SurfCollideStyle(piston/kk,SurfCollidePistonKokkos)

#else

#ifndef SPARTA_SURF_COLLIDE_PISTON_KOKKOS_H
#define SPARTA_SURF_COLLIDE_PISTON_KOKKOS_H

#include "surf_collide_piston.h"
#include "kokkos_type.h"
#include "error.h"
#include "Kokkos_Random.hpp"
#include "rand_pool_wrap.h"
#include "kokkos_copy.h"
#include "fix_ambipolar_kokkos.h"
#include "fix_vibmode_kokkos.h"
#include "surf_react_global_kokkos.h"
#include "surf_react_prob_kokkos.h"

namespace SPARTA_NS {

#define KOKKOS_MAX_SURF_REACT_PER_TYPE 2
#define KOKKOS_MAX_TOT_SURF_REACT 4

class SurfCollidePistonKokkos : public SurfCollidePiston {
 public:

  enum{PKEEP,PINSERT,PDONE,PDISCARD,PENTRY,PEXIT,PSURF};   // several files

  SurfCollidePistonKokkos(class SPARTA *, int, char **);
  SurfCollidePistonKokkos(class SPARTA *);
  ~SurfCollidePistonKokkos();
  void init();
  void pre_collide();
  void post_collide();
  void backup();
  void restore();

 private:

  typedef Kokkos::DualView<int[2], DeviceType::array_layout, DeviceType> tdual_int_2;
  typedef tdual_int_2::t_dev t_int_2;
  typedef tdual_int_2::t_host t_host_int_2;
  t_int_2 d_scalars;
  t_host_int_2 h_scalars;

  DAT::t_int_scalar d_nsingle;
  DAT::t_int_scalar d_nreact_one;

  HAT::t_int_scalar h_nsingle;
  HAT::t_int_scalar h_nreact_one;

  t_particle_1d d_particles;

  int ambi_flag,vibmode_flag;
  FixAmbipolarKokkos* afix_kk;
  FixVibmodeKokkos* vfix_kk;
  KKCopy<FixAmbipolarKokkos> fix_ambi_kk_copy;
  KKCopy<FixVibmodeKokkos> fix_vibmode_kk_copy;

  int sr_type_list[KOKKOS_MAX_TOT_SURF_REACT];
  int sr_map[KOKKOS_MAX_TOT_SURF_REACT];
  KKCopy<SurfReactGlobalKokkos> sr_kk_global_copy[KOKKOS_MAX_SURF_REACT_PER_TYPE];
  KKCopy<SurfReactProbKokkos> sr_kk_prob_copy[KOKKOS_MAX_SURF_REACT_PER_TYPE];

 public:

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

  template<int REACT, int ATOMIC_REDUCTION>
  KOKKOS_INLINE_FUNCTION
  Particle::OnePart* collide_kokkos(Particle::OnePart *&ip, double &dtremain,
                                    int isurf, const double *norm, int isr, int &reaction,
                                    const DAT::t_int_scalar &d_retry, const DAT::t_int_scalar &d_nlocal) const
  {
    if (ATOMIC_REDUCTION == 0)
      d_nsingle()++;
    else
      Kokkos::atomic_increment(&d_nsingle());

    // if surface chemistry defined, attempt reaction
    // reaction = 1 to N for which reaction took place, 0 for none
    // velreset = 1 if reaction reset post-collision velocity, else 0

    Particle::OnePart iorig;
    Particle::OnePart *jp = NULL;
    reaction = 0;
    int velreset = 0;

    if (REACT) {
      if (ambi_flag || vibmode_flag) memcpy(&iorig,ip,sizeof(Particle::OnePart));

      int sr_type = sr_type_list[isr];
      int m = sr_map[isr];

      if (sr_type == 0) {
        reaction = sr_kk_global_copy[m].obj.
          react_kokkos<ATOMIC_REDUCTION>(ip,isurf,norm,jp,velreset,d_retry,d_nlocal);
      } else if (sr_type == 1) {
        reaction = sr_kk_prob_copy[m].obj.
          react_kokkos<ATOMIC_REDUCTION>(ip,isurf,norm,jp,velreset,d_retry,d_nlocal);
      }

      if (reaction) {
        if (ATOMIC_REDUCTION == 0)
          d_nreact_one()++;
        else
          Kokkos::atomic_increment(&d_nreact_one());
      }
    }

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

    // xwall = initial position of wall (collision pt)
    // xorig = initial coordinate component
    // vorig = initial velocity component
    // vwall = user-specified wall velocity (always >= 0)

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

    if (REACT && reaction && ambi_flag) {
      int i = -1;
      if (ip) i = ip - d_particles.data();
      int j = -1;
      if (jp) j = jp - d_particles.data();
      int j_orig = j;
      fix_ambi_kk_copy.obj.surf_react_kokkos(&iorig,i,j);
      if (jp && j < 0) {
        d_particles[j_orig].flag = PDISCARD;
        jp = NULL;
      }
    }

    return jp;
  };
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

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

class SurfCollideSpecularKokkos : public SurfCollideSpecular {
 public:

  enum{PKEEP,PINSERT,PDONE,PDISCARD,PENTRY,PEXIT,PSURF};   // several files

  SurfCollideSpecularKokkos(class SPARTA *, int, char **);
  SurfCollideSpecularKokkos(class SPARTA *);
  ~SurfCollideSpecularKokkos();
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
  Particle::OnePart* collide_kokkos(Particle::OnePart *&ip, double &,
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

    // specular or noslip reflection for each particle
    // only if SurfReact did not already reset velocities
    // also both particles need to trigger any fixes
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

    if (REACT && jp) {
      if (!velreset) {
        if (noslip_flag) {

          // noslip reflection
          // reflect incident v, all three components.

          MathExtraKokkos::negate3(jp->v);
        } else {

          // specular reflection

          MathExtraKokkos::reflect3(jp->v,norm);
        }
      }

    //  if (modify->n_update_custom) {
    //    int j = jp - particle->particles;
    //    modify->update_custom(j,twall,twall,twall,vstream);
    //}
    }

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

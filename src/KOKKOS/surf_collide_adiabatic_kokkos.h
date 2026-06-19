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

#ifdef SURF_COLLIDE_CLASS

SurfCollideStyle(adiabatic/kk,SurfCollideAdiabaticKokkos)

#else

#ifndef SPARTA_SURF_COLLIDE_ADIABATIC_KOKKOS_H
#define SPARTA_SURF_COLLIDE_ADIABATIC_KOKKOS_H

#include "surf_collide_adiabatic.h"
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

class SurfCollideAdiabaticKokkos : public SurfCollideAdiabatic {
 public:

  enum{PKEEP,PINSERT,PDONE,PDISCARD,PENTRY,PEXIT,PSURF};   // several files

  SurfCollideAdiabaticKokkos(class SPARTA *, int, char **);
  SurfCollideAdiabaticKokkos(class SPARTA *);
  ~SurfCollideAdiabaticKokkos();
  void init();
  void pre_collide();
  void post_collide();
  void backup();
  void restore();

 private:

#ifndef SPARTA_KOKKOS_EXACT
  Kokkos::Random_XorShift64_Pool<DeviceType> rand_pool;
  typedef typename Kokkos::Random_XorShift64_Pool<DeviceType>::generator_type rand_type;
#else
  RandPoolWrap rand_pool;
  typedef RandWrap rand_type;
#endif

  RanKnuth* random_backup;

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
  t_species_1d d_species;

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

     note that the adiabatic condition (i.e. no energy transfer of flow to
     surf) only applies to particle collisions.  Chemistry (e.g. particle
     adsorptions) can lead to energy transfer in both directions, so the
     velocities reset by SurfReact are kept as-is.
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
      Kokkos::atomic_inc(&d_nsingle());

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
          Kokkos::atomic_inc(&d_nreact_one());
      }
    }

    // isotropic scattering conserving velocity magnitude (kinetic energy)
    //   of each particle
    // only if SurfReact did not already reset velocities
    // cannot trigger fixes that require temperature of particle here
    //   because temperature of wall is not known

    if (ip) {
      if (!velreset) scatter_isotropic(ip,norm);
    }
    if (REACT && jp) {
      if (!velreset) scatter_isotropic(jp,norm);
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

 private:

  /* ----------------------------------------------------------------------
     particle collision with adiabatic surface
     p = particle with current x = collision pt, current v = incident v
     norm = surface normal unit vector
     particle is scattered isotropically while conserving its velocity
     magnitude (i.e. no energy transfer to surf)
  ------------------------------------------------------------------------- */

  KOKKOS_INLINE_FUNCTION
  void scatter_isotropic(Particle::OnePart *p, const double *norm) const
  {
    rand_type rand_gen = rand_pool.get_state();

    double *v = p->v;
    double dot = MathExtraKokkos::dot3(v,norm);

    // tangent1/2 = surface tangential unit vectors

    double tangent1[3],tangent2[3];
    tangent1[0] = v[0] - dot*norm[0];
    tangent1[1] = v[1] - dot*norm[1];
    tangent1[2] = v[2] - dot*norm[2];

    // if mag(tangent1) == 0, normal collision: choose a random tangent vector

    if (MathExtraKokkos::lensq3(tangent1) == 0.0) {
      tangent2[0] = rand_gen.drand();
      tangent2[1] = rand_gen.drand();
      tangent2[2] = rand_gen.drand();
      MathExtraKokkos::cross3(norm,tangent2,tangent1);
    }

    MathExtraKokkos::norm3(tangent1);
    MathExtraKokkos::cross3(norm,tangent1,tangent2);

    // isotropic scattering
    // vmag = magnitude of incident particle velocity vector
    // vperp = velocity component perpendicular to surface along norm
    // vtan1/2 = 2 remaining velocity components tangential to surface

    double vmag = MathExtraKokkos::len3(v);

    double theta = MathConst::MY_2PI * rand_gen.drand();
    double f_phi = rand_gen.drand();
    double sqrt_f_phi = sqrt(f_phi);

    double vperp = vmag * sqrt(1.0 - f_phi);
    double vtan1 = vmag * sqrt_f_phi * sin(theta);
    double vtan2 = vmag * sqrt_f_phi * cos(theta);

    v[0] = vperp*norm[0] + vtan1*tangent1[0] + vtan2*tangent2[0];
    v[1] = vperp*norm[1] + vtan1*tangent1[1] + vtan2*tangent2[1];
    v[2] = vperp*norm[2] + vtan1*tangent1[2] + vtan2*tangent2[2];

    // p->erot and p->evib stay identical

    rand_pool.free_state(rand_gen);
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

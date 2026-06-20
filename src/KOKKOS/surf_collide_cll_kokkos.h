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

SurfCollideStyle(cll/kk,SurfCollideCLLKokkos)

#else

#ifndef SPARTA_SURF_COLLIDE_CLL_KOKKOS_H
#define SPARTA_SURF_COLLIDE_CLL_KOKKOS_H

#include "surf_collide_cll.h"
#include "kokkos_type.h"
#include "math_extra_kokkos.h"
#include "Kokkos_Random.hpp"
#include "rand_pool_wrap.h"
#include "kokkos_copy.h"
#include "fix_ambipolar_kokkos.h"
#include "fix_vibmode_kokkos.h"
#include "surf_react_global_kokkos.h"
#include "surf_react_prob_kokkos.h"
#include "surf_react_adsorb_kokkos.h"

namespace SPARTA_NS {

class SurfCollideCLLKokkos : public SurfCollideCLL {
 public:

  enum{NONE,DISCRETE,SMOOTH};                              // several files
  enum{PKEEP,PINSERT,PDONE,PDISCARD,PENTRY,PEXIT,PSURF};   // several files

  SurfCollideCLLKokkos(class SPARTA *, int, char **);
  SurfCollideCLLKokkos(class SPARTA *);
  ~SurfCollideCLLKokkos();
  void init();
  void dynamic();
  void pre_collide();
  void post_collide();
  void backup();
  void restore();

 private:
  double boltz;
  int rotstyle, vibstyle;

#ifndef SPARTA_KOKKOS_EXACT
  Kokkos::Random_XorShift64_Pool<DeviceType> rand_pool;
  typedef typename Kokkos::Random_XorShift64_Pool<DeviceType>::generator_type rand_type;
#else
  RandPoolWrap rand_pool;
  typedef RandWrap rand_type;
#endif

  RanKnuth* random_backup;

  DAT::t_float_1d d_t_persurf;

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
  KKCopy<SurfReactAdsorbKokkos> sr_kk_adsorb_copy[KOKKOS_MAX_SURF_REACT_PER_TYPE];

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
      Kokkos::atomic_inc(&d_nsingle());

    // if surface chemistry defined, attempt reaction
    // reaction = 1 to N for which reaction took place, 0 for none
    // velreset = 1 if reaction reset post-collision velocity, else 0

    Particle::OnePart iorig;
    Particle::OnePart *jp = NULL;
    reaction = 0;
    int velreset = 0;

    if (REACT && isr >= 0) {
      if (ambi_flag || vibmode_flag) memcpy(&iorig,ip,sizeof(Particle::OnePart));

      int sr_type = sr_type_list[isr];
      int m = sr_map[isr];

      if (sr_type == 0) {
        reaction = sr_kk_global_copy[m].obj.
          react_kokkos<ATOMIC_REDUCTION>(ip,isurf,norm,jp,velreset,d_retry,d_nlocal);
      } else if (sr_type == 1) {
        reaction = sr_kk_prob_copy[m].obj.
          react_kokkos<ATOMIC_REDUCTION>(ip,isurf,norm,jp,velreset,d_retry,d_nlocal);
      } else if (sr_type == 2) {
        reaction = sr_kk_adsorb_copy[m].obj.
          react_kokkos<ATOMIC_REDUCTION>(ip,isurf,norm,jp,velreset,d_retry,d_nlocal);
      }

      if (reaction) {
        if (ATOMIC_REDUCTION == 0)
          d_nreact_one()++;
        else
          Kokkos::atomic_inc(&d_nreact_one());
      }
    }

    // set temperature of isurf if VARSURF or CUSTOM

    double tsurf_local = tsurf;
    if (persurf_temperature) {
      tsurf_local = d_t_persurf[isurf];
      if (tsurf_local <= 0.0) Kokkos::abort("Surf_collide tsurf <= 0.0");
    }

    // CLL reflection for each particle
    // only if SurfReact did not already reset velocities
    // also both particles need to trigger any fixes
    //   to update per-particle properties which depend on
    //   temperature of the particle, e.g. fix vibmode and fix ambipolar

    if (ip) {
      if (!velreset) cll(ip,norm,tsurf_local);
      int i = ip - d_particles.data();
      if (ambi_flag)
        fix_ambi_kk_copy.obj.update_custom_kokkos(i,tsurf_local,tsurf_local,tsurf_local,vstream);
      if (vibmode_flag)
        fix_vibmode_kk_copy.obj.update_custom_kokkos(i,tsurf_local,tsurf_local,tsurf_local,vstream);
    }
    if (REACT && jp) {
      if (!velreset) cll(jp,norm,tsurf_local);
      int j = jp - d_particles.data();
      if (ambi_flag)
        fix_ambi_kk_copy.obj.update_custom_kokkos(j,tsurf_local,tsurf_local,tsurf_local,vstream);
      if (vibmode_flag)
        fix_vibmode_kk_copy.obj.update_custom_kokkos(j,tsurf_local,tsurf_local,tsurf_local,vstream);
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
     cll reflection
  ------------------------------------------------------------------------- */

  KOKKOS_INLINE_FUNCTION
  void cll(Particle::OnePart *p, const double *norm, const double twall) const
  {
    rand_type rand_gen = rand_pool.get_state();

    double tangent1[3],tangent2[3];
    int ispecies = p->ispecies;
    double beta_un,normalized_distbn_fn;

    double *v = p->v;
    double dot = MathExtraKokkos::dot3(v,norm);
    double vrm, vperp, vtan1, vtan2;

    tangent1[0] = v[0] - dot*norm[0];
    tangent1[1] = v[1] - dot*norm[1];
    tangent1[2] = v[2] - dot*norm[2];

    if (MathExtraKokkos::lensq3(tangent1) == 0.0) {
      tangent2[0] = rand_gen.drand();
      tangent2[1] = rand_gen.drand();
      tangent2[2] = rand_gen.drand();
      MathExtraKokkos::cross3(norm,tangent2,tangent1);
    }

    MathExtraKokkos::norm3(tangent1);
    MathExtraKokkos::cross3(norm,tangent1,tangent2);

    double tan1 = MathExtraKokkos::dot3(v,tangent1);

    vrm = sqrt(2.0*boltz * twall / d_species[ispecies].mass);

    // CLL model normal velocity

    double r_1 = sqrt(-acc_n*log(rand_gen.drand()));
    double theta_1 = MathConst::MY_2PI * rand_gen.drand();
    double dot_norm = dot/vrm * sqrt(1-acc_n);
    vperp = vrm * sqrt(r_1*r_1 + dot_norm*dot_norm + 2*r_1*dot_norm*cos(theta_1));

    // CLL model tangential velocities

    double r_2 = sqrt(-acc_t*log(rand_gen.drand()));
    double theta_2 = MathConst::MY_2PI * rand_gen.drand();
    double vtangent = tan1/vrm * sqrt(1-acc_t);
    vtan1 = vrm * (vtangent + r_2*cos(theta_2));
    vtan2 = vrm * r_2 * sin(theta_2);

    // partial keyword
    // incomplete energy accommodation with partial/fully diffuse scattering
    // adjust the final angle of the particle while keeping
    //   the velocity magnitude or speed according to CLL scattering

    if (pflag) {
      double tan2 = MathExtraKokkos::dot3(v,tangent2);
      double phi_i, psi_i, theta_f, phi_f, psi_f, cos_beta;

      psi_i = acos(dot*dot/MathExtraKokkos::lensq3(v));
      phi_i = atan2(tan2,tan1);

      double v_mag = sqrt(vperp*vperp + vtan1*vtan1 + vtan2*vtan2);

      double P = 0;
      while (rand_gen.drand() > P) {
        phi_f = MathConst::MY_2PI*rand_gen.drand();
        psi_f = acos(1-rand_gen.drand());
        cos_beta =  cos(psi_i)*cos(psi_f) +
          sin(psi_i)*sin(psi_f)*cos(phi_i - phi_f);
        P = (1-eccen)/(1-eccen*cos_beta);
      }

      theta_f = acos(sqrt(cos(psi_f)));

      vperp = v_mag * cos(theta_f);
      vtan1 = v_mag * sin(theta_f) * cos(phi_f);
      vtan2 = v_mag * sin(theta_f) * sin(phi_f);
    }

    // add in translation or rotation vector if specified
    // only keep portion of vector tangential to surface element

    if (trflag) {
      double vxdelta,vydelta,vzdelta;
      if (tflag) {
        vxdelta = vx; vydelta = vy; vzdelta = vz;
        double dot = vxdelta*norm[0] + vydelta*norm[1] + vzdelta*norm[2];

        if (fabs(dot) > 0.001) {
          dot /= vrm;
          do {
            do {
              beta_un = (6.0*rand_gen.normal() - 3.0);
            } while (beta_un + dot < 0.0);
            normalized_distbn_fn = 2.0 * (beta_un + dot) /
              (dot + sqrt(dot*dot + 2.0)) *
              exp(0.5 + (0.5*dot)*(dot-sqrt(dot*dot + 2.0)) - beta_un*beta_un);
          } while (normalized_distbn_fn < rand_gen.drand());
          vperp = beta_un*vrm;
        }

      } else {
        double *x = p->x;
        vxdelta = wy*(x[2]-pz) - wz*(x[1]-py);
        vydelta = wz*(x[0]-px) - wx*(x[2]-pz);
        vzdelta = wx*(x[1]-py) - wy*(x[0]-px);
        double dot = vxdelta*norm[0] + vydelta*norm[1] + vzdelta*norm[2];
        vxdelta -= dot*norm[0];
        vydelta -= dot*norm[1];
        vzdelta -= dot*norm[2];
      }

      v[0] = vperp*norm[0] + vtan1*tangent1[0] + vtan2*tangent2[0] + vxdelta;
      v[1] = vperp*norm[1] + vtan1*tangent1[1] + vtan2*tangent2[1] + vydelta;
      v[2] = vperp*norm[2] + vtan1*tangent1[2] + vtan2*tangent2[2] + vzdelta;

    // no translation or rotation

    } else {
      v[0] = vperp*norm[0] + vtan1*tangent1[0] + vtan2*tangent2[0];
      v[1] = vperp*norm[1] + vtan1*tangent1[1] + vtan2*tangent2[1];
      v[2] = vperp*norm[2] + vtan1*tangent1[2] + vtan2*tangent2[2];
    }

    // rotational component

    if (rotstyle == NONE || d_species[ispecies].rotdof < 2) p->erot = 0.0;

    else {
      double erot_mag = sqrt(p->erot*(1-acc_rot)/(boltz*twall));

      double r_rot,cos_theta_rot,A_rot,X_rot;
      if (d_species[ispecies].rotdof == 2) {
        r_rot = sqrt(-acc_rot*log(rand_gen.drand()));
        cos_theta_rot = cos(MathConst::MY_2PI*rand_gen.drand());
      }
      else if (d_species[ispecies].rotdof > 2) {
        A_rot = 0;
        while (A_rot < rand_gen.drand()) {
          X_rot = 4*rand_gen.drand();
          A_rot = 2.71828182845904523536028747*X_rot*X_rot*exp(-X_rot*X_rot);
        }
        r_rot = sqrt(acc_rot)*X_rot;
        cos_theta_rot = 2*rand_gen.drand() - 1;
      }

      p->erot = boltz * twall *
        (r_rot*r_rot + erot_mag*erot_mag + 2*r_rot*erot_mag*cos_theta_rot);
      }

    // vibrational component

    int vibdof = d_species[ispecies].vibdof;
    double r_vib, cos_theta_vib, A_vib, X_vib, evib_mag, evib_val;

    if (vibstyle == NONE || vibdof < 2)
      p->evib = 0.0;

    else if (vibstyle == DISCRETE && vibdof == 2) {
      double evib_star =
        -log(1 - rand_gen.drand() *
             (1 - exp(-boltz*d_species[ispecies].vibtemp[0])));
      evib_val = p->evib + evib_star;
      evib_mag = sqrt(evib_val*(1-acc_vib)/(boltz*twall));
      r_vib = sqrt(-acc_vib*log(rand_gen.drand()));
      cos_theta_vib = cos(MathConst::MY_2PI*rand_gen.drand());
      evib_val = boltz * twall *
        (r_vib*r_vib + evib_mag*evib_mag + 2*r_vib*evib_mag*cos_theta_vib);
      int ivib =  evib_val / (boltz*d_species[ispecies].vibtemp[0]);
      p->evib = ivib * boltz * d_species[ispecies].vibtemp[0];
    }

    else if (vibstyle == SMOOTH || vibdof >= 2) {
      evib_mag = sqrt(p->evib*(1-acc_vib)/(boltz*twall));
      if (vibdof == 2) {
        r_vib = sqrt(-acc_vib*log(rand_gen.drand()));
        cos_theta_vib = cos(MathConst::MY_2PI*rand_gen.drand());
      } else if (vibdof > 2) {
        A_vib = 0;
        while (A_vib < rand_gen.drand()) {
          X_vib = 4*rand_gen.drand();
          A_vib = 2.71828182845904523536028747*X_vib*X_vib*exp(-X_vib*X_vib);
        }
        r_vib = sqrt(acc_vib)*X_vib;
        cos_theta_vib = 2*rand_gen.drand() - 1;
      }

      p->evib = boltz * twall *
        (r_vib*r_vib + evib_mag*evib_mag + 2*r_vib*evib_mag*cos_theta_vib);
    }

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

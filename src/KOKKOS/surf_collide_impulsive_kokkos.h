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

SurfCollideStyle(impulsive/kk,SurfCollideImpulsiveKokkos)

#else

#ifndef SPARTA_SURF_COLLIDE_IMPULSIVE_KOKKOS_H
#define SPARTA_SURF_COLLIDE_IMPULSIVE_KOKKOS_H

#include "surf_collide_impulsive.h"
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

class SurfCollideImpulsiveKokkos : public SurfCollideImpulsive {
 public:

  enum{NONE,DISCRETE,SMOOTH};                              // several files
  enum{PKEEP,PINSERT,PDONE,PDISCARD,PENTRY,PEXIT,PSURF};   // several files

  SurfCollideImpulsiveKokkos(class SPARTA *, int, char **);
  SurfCollideImpulsiveKokkos(class SPARTA *);
  ~SurfCollideImpulsiveKokkos();
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

    // set temperature of isurf if VARSURF or CUSTOM

    double tsurf_local = tsurf;
    if (persurf_temperature) {
      tsurf_local = d_t_persurf[isurf];
      if (tsurf_local <= 0.0) Kokkos::abort("Surf_collide tsurf <= 0.0");
    }

    // impulsive reflection for each particle
    // only if SurfReact did not already reset velocities
    // also both particles need to trigger any fixes
    //   to update per-particle properties which depend on
    //   temperature of the particle, e.g. fix vibmode and fix ambipolar

    if (ip) {
      if (!velreset) impulsive(ip,norm,tsurf_local);
      int i = ip - d_particles.data();
      if (ambi_flag)
        fix_ambi_kk_copy.obj.update_custom_kokkos(i,tsurf_local,tsurf_local,tsurf_local,vstream);
      if (vibmode_flag)
        fix_vibmode_kk_copy.obj.update_custom_kokkos(i,tsurf_local,tsurf_local,tsurf_local,vstream);
    }
    if (REACT && jp) {
      if (!velreset) impulsive(jp,norm,tsurf_local);
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
     impulsive reflection
  ------------------------------------------------------------------------- */

  KOKKOS_INLINE_FUNCTION
  void impulsive(Particle::OnePart *p, const double *norm, const double twall) const
  {
    rand_type rand_gen = rand_pool.get_state();

    double tangent1[3],tangent2[3];
    int ispecies = p->ispecies;

    double vperp, vtan1, vtan2;
    double mass = d_species[ispecies].mass;

    double *v = p->v;
    double dot = MathExtraKokkos::dot3(v,norm);

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

    // compute final polar (theta) and azimuthal (phi) angles

    double tan1 = MathExtraKokkos::dot3(v,tangent1);
    double tan2 = MathExtraKokkos::dot3(v,tangent2);

    double v_i_mag_sq = MathExtraKokkos::lensq3(v);
    double E_i = 0.5 * mass * v_i_mag_sq;
    double theta_i = acos(-dot/sqrt(v_i_mag_sq));
    double phi_i = atan2(tan2,tan1);
    double phi_peak = MathConst::MY_2PI - phi_i;

    double theta_f, phi_f;
    double P = 0.0;

    // theta_f calculation

    while (rand_gen.drand() > P) {
      theta_f = MathConst::MY_PI2 * rand_gen.drand();
      P = pow(cos( theta_f - theta_peak ),cos_theta_pow) * sin(theta_f);
      if (double_flag) {
        if (theta_f > theta_peak)
          P = pow(cos( theta_f - theta_peak ),cos_theta_pow_2) * sin(theta_f);
      }

      if (step_flag) {
        double func_step = 0.0;
        double tan_theta = tan(theta_f);
        double cotangent = 1.0/tan_theta;
        if (cotangent > step_size) func_step = 1 - step_size*tan_theta;
        P *= func_step;
      }
    }

    // phi_f calculations

    P = 0.0;
    while (rand_gen.drand() > P) {
      phi_f = phi_peak + MathConst::MY_PI * (2*rand_gen.drand() - 1);
      P = pow(cos( 0.5*(phi_f - phi_peak) ),cos_phi_pow);
    }

    if (phi_f > MathConst::MY_PI) phi_f -= MathConst::MY_2PI;
    else if (phi_f < -MathConst::MY_PI) phi_f += MathConst::MY_2PI;

    double v_f_avg = 0.0;
    if (softsphere_flag) {
      double mu = d_species[ispecies].molwt/eff_mass;
      double cos_khi = cos(MathConst::MY_PI - theta_i - theta_f);
      double sin_khi_sq = 1 - cos_khi*cos_khi;
      double dE, E_f_avg;

      dE = 2*mu/((mu+1)*(mu+1)) *
        (1 + mu*sin_khi_sq + eng_ratio*(mu+1)/(2*mu) -
         cos_khi*sqrt(1 - mu*mu*sin_khi_sq - eng_ratio*(mu + 1)));
      E_f_avg = E_i * (1 - dE);
      v_f_avg = var_alpha_sq * sqrt(mass/(2*E_f_avg)) *
        (2*E_f_avg/(mass*var_alpha_sq) - 1);
    } else {
      v_f_avg = u0_a*twall + u0_b;
    }

    double v_f_max = 0.5 * (v_f_avg + sqrt(v_f_avg*v_f_avg + 6*var_alpha_sq));
    double f_max = v_f_max*v_f_max*v_f_max *
      exp(-(v_f_max - v_f_avg) * (v_f_max - v_f_avg)/(var_alpha_sq));

    double v_f_mag;
    P = 0.0;
    while (rand_gen.drand() > P) {
      v_f_mag = v_f_max + 3 * var_alpha * ( 2 * rand_gen.drand() - 1 );
      P = v_f_mag*v_f_mag*v_f_mag/(f_max) *
        exp(-(v_f_mag - v_f_avg)*(v_f_mag - v_f_avg)/(var_alpha_sq));
    }

    vperp = v_f_mag * cos(theta_f);
    vtan1 = v_f_mag * sin(theta_f) * cos(phi_f);
    vtan2 = v_f_mag * sin(theta_f) * sin(phi_f);

    v[0] = vperp*norm[0] + vtan1*tangent1[0] + vtan2*tangent2[0];
    v[1] = vperp*norm[1] + vtan1*tangent1[1] + vtan2*tangent2[1];
    v[2] = vperp*norm[2] + vtan1*tangent1[2] + vtan2*tangent2[2];

    if (intenergy_flag) {
      double E_f = 0.5 * mass * v_f_mag * v_f_mag;
      double extra_energy = E_i - E_f;

      // rotational component

      if (rotstyle == NONE || d_species[ispecies].rotdof < 2) p->erot = 0.0;
      else p->erot += rot_frac*extra_energy;

      // vibrational component

      int vibdof = d_species[ispecies].vibdof;

      if (vibstyle == NONE || vibdof < 2) {
        p->evib = 0.0;
      } else {
        double *vibtemp = d_species[ispecies].vibtemp;
        double evib_val = p->evib + vib_frac*extra_energy;

        if (vibstyle == SMOOTH) p->evib = evib_val;
        if (vibstyle == DISCRETE && vibdof==2) {
          int ivib =  evib_val / (boltz*vibtemp[0]);
          p->evib = ivib * boltz * vibtemp[0];
        } else {
          int nvibmode = d_species[ispecies].nvibmode;
          int *vibdegen = d_species[ispecies].vibdegen;
          double tot_temp=0.0;
          double evib_sum = 0.0;

          for (int imode=0; imode<nvibmode; imode++)
            tot_temp += vibtemp[imode]*vibdegen[imode];

          for (int imode=0; imode<nvibmode; imode++) {
            int ivib = evib_val / (boltz*tot_temp);
            evib_sum += ivib * boltz * vibtemp[imode]*vibdegen[imode];
          }

          p->evib = evib_sum;
        }
      }
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

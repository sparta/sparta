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

SurfCollideStyle(td/kk,SurfCollideTDKokkos)

#else

#ifndef SPARTA_SURF_COLLIDE_TD_KOKKOS_H
#define SPARTA_SURF_COLLIDE_TD_KOKKOS_H

#include "surf_collide_td.h"
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

class SurfCollideTDKokkos : public SurfCollideTD {
 public:

  enum{NONE,DISCRETE,SMOOTH};                              // several files
  enum{PKEEP,PINSERT,PDONE,PDISCARD,PENTRY,PEXIT,PSURF};   // several files

  SurfCollideTDKokkos(class SPARTA *, int, char **);
  SurfCollideTDKokkos(class SPARTA *);
  ~SurfCollideTDKokkos();
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

    // TD reflection for each particle
    // only if SurfReact did not already reset velocities
    // also both particles need to trigger any fixes
    //   to update per-particle properties which depend on
    //   temperature of the particle, e.g. fix vibmode and fix ambipolar

    if (ip) {
      if (!velreset) td(ip,norm,tsurf_local);
      int i = ip - d_particles.data();
      if (ambi_flag)
        fix_ambi_kk_copy.obj.update_custom_kokkos(i,tsurf_local,tsurf_local,tsurf_local,vstream);
      if (vibmode_flag)
        fix_vibmode_kk_copy.obj.update_custom_kokkos(i,tsurf_local,tsurf_local,tsurf_local,vstream);
    }
    if (REACT && jp) {
      if (!velreset) td(jp,norm,tsurf_local);
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
     TD (thermal desorption) reflection
  ------------------------------------------------------------------------- */

  KOKKOS_INLINE_FUNCTION
  void td(Particle::OnePart *p, const double *norm, const double twall) const
  {
    rand_type rand_gen = rand_pool.get_state();

    double tangent1[3],tangent2[3];
    int ispecies = p->ispecies;

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

    double mass = d_species[ispecies].mass;
    double E_i = 0.5 * mass * MathExtraKokkos::lensq3(v);

    double E_t = boltz * twall;
    if (bond_flag) E_t += boltz*bond_trans;
    if (initen_flag) E_t += E_i*initen_trans;

    double E_n = E_t;
    if (barrier_flag) E_n += boltz*barrier_val;

    double vrm_n = sqrt(2.0*E_n / mass);
    double vrm_t = sqrt(2.0*E_t / mass);
    double vperp = vrm_n * sqrt(-log(rand_gen.drand()));

    double theta = MathConst::MY_2PI * rand_gen.drand();
    double vtangent = vrm_t * sqrt(-log(rand_gen.drand()));
    double vtan1 = vtangent * sin(theta);
    double vtan2 = vtangent * cos(theta);

    v[0] = vperp*norm[0] + vtan1*tangent1[0] + vtan2*tangent2[0];
    v[1] = vperp*norm[1] + vtan1*tangent1[1] + vtan2*tangent2[1];
    v[2] = vperp*norm[2] + vtan1*tangent1[2] + vtan2*tangent2[2];

    double twall_rot = twall;
    double twall_vib = twall;

    if (bond_flag) {
      twall_rot += bond_rot;
      twall_vib += bond_vib;
    }

    if (initen_flag) {
      twall_rot += E_i*initen_rot/boltz;
      twall_vib += E_i*initen_vib/boltz;
    }

    p->erot = erot(ispecies,twall_rot,rand_gen,boltz);
    p->evib = evib(ispecies,twall_vib,rand_gen,boltz);

    rand_pool.free_state(rand_gen);
  }

  /* ----------------------------------------------------------------------
     generate random rotational energy for a particle
     only a function of species index and species properties
  ------------------------------------------------------------------------- */

  KOKKOS_INLINE_FUNCTION
  double erot(int isp, double temp_thermal, rand_type &rand_gen, double boltz) const
  {
    double eng,a,erm,b;

    if (rotstyle == NONE) return 0.0;
    if (d_species[isp].rotdof < 2) return 0.0;

    if (rotstyle == DISCRETE && d_species[isp].rotdof == 2) {
      int irot = -log(rand_gen.drand()) * temp_thermal /
        d_species[isp].rottemp[0];
      eng = irot * boltz * d_species[isp].rottemp[0];
    } else if (rotstyle == SMOOTH && d_species[isp].rotdof == 2) {
      eng = -log(rand_gen.drand()) * boltz * temp_thermal;
    } else {
      a = 0.5*d_species[isp].rotdof-1.0;
      while (1) {
        // energy cut-off at 10 kT
        erm = 10.0*rand_gen.drand();
        b = pow(erm/a,a) * exp(a-erm);
        if (b > rand_gen.drand()) break;
      }
      eng = erm * boltz * temp_thermal;
    }

   return eng;
  }

  /* ----------------------------------------------------------------------
     generate random vibrational energy for a particle
     only a function of species index and species properties
     index_vibmode = index of extra per-particle vibrational mode storage
       -1 if not defined for this model
  ------------------------------------------------------------------------- */

  KOKKOS_INLINE_FUNCTION
  double evib(int isp, double temp_thermal, rand_type &rand_gen, double boltz) const
  {
    double eng,a,erm,b;

    if (vibstyle == NONE || d_species[isp].vibdof < 2) return 0.0;

    // for DISCRETE, only need set evib for vibdof = 2
    // mode levels and evib will be set by FixVibmode::update_custom()

    eng = 0.0;

    if (vibstyle == DISCRETE && d_species[isp].vibdof == 2) {
      int ivib = -log(rand_gen.drand()) * temp_thermal /
        d_species[isp].vibtemp[0];
      eng = ivib * boltz * d_species[isp].vibtemp[0];
    } else if (vibstyle == SMOOTH || d_species[isp].vibdof >= 2) {
      if (d_species[isp].vibdof == 2)
        eng = -log(rand_gen.drand()) * boltz * temp_thermal;
      else if (d_species[isp].vibdof > 2) {
        a = 0.5*d_species[isp].vibdof-1.;
        while (1) {
          // energy cut-off at 10 kT
          erm = 10.0*rand_gen.drand();
          b = pow(erm/a,a) * exp(a-erm);
          if (b > rand_gen.drand()) break;
        }
        eng = erm * boltz * temp_thermal;
      }
    }

    return eng;
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

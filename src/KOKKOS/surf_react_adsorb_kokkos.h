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

#ifdef SURF_REACT_CLASS

SurfReactStyle(adsorb/kk,SurfReactAdsorbKokkos)

#else

#ifndef SPARTA_SURF_REACT_ADSORB_KOKKOS_H
#define SPARTA_SURF_REACT_ADSORB_KOKKOS_H

#include "surf_react_adsorb.h"
#include "kokkos_type.h"
#include "math_extra_kokkos.h"
#include "math_const.h"
#include "rand_pool_wrap.h"
#include "Kokkos_Random.hpp"
#include "particle_kokkos.h"

namespace SPARTA_NS {

// must match enums in surf_react_adsorb.cpp

namespace SRA_KK {
  enum{DISSOCIATION,EXCHANGE,RECOMBINATION,AA,DA,LH1,LH3,CD,ER,CI};
  enum{FACE,SURF};
  enum{NOMODEL,SPECULAR,DIFFUSE,ADIABATIC,CLL,TD,IMPULSIVE,MAXMODELS};
  enum{SIMPLE,ARRHENIUS};
  enum{NONE,DISCRETE,SMOOTH};   // rotstyle/vibstyle, must match collide.h
}

#define SRA_KK_MAXREACTANT 5
#define SRA_KK_MAXPRODUCT 5
#define SRA_KK_MAXCOEFF 4
#define SRA_KK_MAXPERSPECIES 16   // max GS reactions a single species can be in
#define SRA_KK_MAXMODELS 7        // = MAXMODELS
#define SRA_KK_MAXCMCOEFF 11      // max cmodel coeffs (impulsive)
#define SRA_KK_MAXCMFLAG 4        // max cmodel flags (impulsive)

class SurfReactAdsorbKokkos : public SurfReactAdsorb {
 public:
  SurfReactAdsorbKokkos(class SPARTA *, int, char **);
  SurfReactAdsorbKokkos(class SPARTA *);
  ~SurfReactAdsorbKokkos();
  void init();
  void tally_reset();
  void tally_update();

  void pre_react();
  void backup();
  void restore();

 private:
  // flattened GS reaction tables (indexed by reaction j in 0..nlist_gs)

  DAT::t_int_1d d_reactions_n;       // # of GS reactions for each species
  DAT::t_int_2d d_list;              // per-species list of reaction indices

  DAT::t_int_1d d_type;              // reaction type (DISSOCIATION,...)
  DAT::t_int_1d d_style;             // SIMPLE or ARRHENIUS
  DAT::t_float_1d d_kreact;          // precomputed rate coefficient
  DAT::t_int_1d d_kisliuk_flag;
  DAT::t_float_2d d_kisliuk;         // [j][3]
  DAT::t_int_1d d_energy_flag;
  DAT::t_float_2d d_energy;          // [j][2]
  DAT::t_float_2d d_coeff;           // [j][MAXCOEFF]
  DAT::t_int_1d d_nreactant,d_nproduct;
  DAT::t_int_1d d_nprod_g,d_nprod_g_tot;
  DAT::t_int_1d d_cmodel_ip,d_cmodel_jp;

  DAT::t_int_2d d_rstate,d_rpart,d_rstoich,d_rad;   // reactant slots [j][MAXREACTANT]
  DAT::t_int_2d d_pstate,d_ppart,d_pstoich,d_pad;   // product slots [j][MAXPRODUCT]
  DAT::t_int_2d d_products;                         // product species indices

  // per-state-slot data: FACE mode => 6 box faces; SURF mode => nlocal+nghost surfs

  int nstate_;                             // # of state slots (nface or nall)
  DAT::t_int_1d d_total_state;             // [nstate]
  DAT::t_float_1d d_area,d_weight;         // [nstate]
  DAT::t_int_2d d_species_state;           // [nstate][nspecies_surf]
  DAT::t_int_2d d_species_delta;           // [nstate][nspecies_surf] (atomic)

  DAT::tdual_int_2d k_species_delta;

  DAT::tdual_int_1d k_mark;                // [nstate] reacted-surf mark (SURF)
  DAT::t_int_1d d_mark;

  double fnum_;        // update->fnum, set in pre_react

  // post-reaction collision model (cmodel) state for bit-exact device scatter
  // per-reaction flattened cmodel coeffs/flags (ip and jp), plus collide
  //   rot/vib styles and boltz captured in pre_react

  DAT::t_int_2d d_cmip_flags,d_cmjp_flags;     // [nr][MAXCMFLAG]
  DAT::t_float_2d d_cmip_coeffs,d_cmjp_coeffs; // [nr][MAXCMCOEFF]
  double boltz_;
  int rotstyle_,vibstyle_;

  // one RNG per cmodel type, wrapping that cmodel's RanKnuth, so device
  //   scatter draws match the host wrapper bit-for-bit (EXACT serial)

  RandPoolWrap *cmodel_pool[SRA_KK_MAXMODELS];
#ifdef SPARTA_KOKKOS_EXACT
  RandWrap d_cmodel_rand[SRA_KK_MAXMODELS];
#endif

  void init_reactions_gs_kokkos();
  void init_cmodels_kokkos();

#ifndef SPARTA_KOKKOS_EXACT
  Kokkos::Random_XorShift64_Pool<DeviceType> rand_pool;
  typedef typename Kokkos::Random_XorShift64_Pool<DeviceType>::generator_type rand_type;
#else
  RandPoolWrap rand_pool;
  typedef RandWrap rand_type;
#endif

  RanKnuth* random_backup;

  DAT::t_int_1d d_scalars;
  HAT::t_int_1d h_scalars;
  DAT::t_int_scalar d_nsingle;
  DAT::t_int_1d d_tally_single;
  HAT::t_int_scalar h_nsingle;
  HAT::t_int_1d h_tally_single;

  t_particle_1d d_particles;
  t_species_1d d_species;

 public:

  /* ----------------------------------------------------------------------
     select GS surface reaction to perform for particle IP on box face
     mirrors SurfReactAdsorb::react() for mode == FACE, gsflag == 1
     return reaction 1 to N, 0 = no reaction
     only reaction types/cmodels validated at init are reachable here
  ------------------------------------------------------------------------- */

  template<int ATOMIC_REDUCTION>
  KOKKOS_INLINE_FUNCTION
  int react_kokkos(Particle::OnePart *&ip, int isurf, const double *norm,
                   Particle::OnePart *&jp, int &velreset,
                   const DAT::t_int_scalar &d_retry,
                   const DAT::t_int_scalar &d_nlocal) const
  {
    // FACE: convert negative face code to 0..5; SURF: use local surf index

    int idx;
    if (mode == SRA_KK::FACE) idx = -(isurf+1);
    else idx = isurf;

    int n = d_reactions_n[ip->ispecies];
    if (n == 0) return 0;

    double fnum = fnum_;
    long int maxstick = ceil(max_cover*d_area[idx] / (fnum*d_weight[idx]));
    double factor = fnum * d_weight[idx] / d_area[idx];
    double ms_inv = factor / max_cover;

    double prob_value[SRA_KK_MAXPERSPECIES];
    double sum_prob = 0.0;
    double scatter_prob = 0.0, correction = 1.0;
    int coeff_val = 1;

    rand_type rand_gen = rand_pool.get_state();

    for (int i = 0; i < n; i++) {
      int j = d_list(ip->ispecies,i);

      if (d_style(j) == SRA_KK::ARRHENIUS) coeff_val = 3;

      double surf_cover,S_theta,K_ads;

      switch (d_type(j)) {
      case SRA_KK::DISSOCIATION:
      case SRA_KK::EXCHANGE:
      case SRA_KK::RECOMBINATION:
        prob_value[i] = d_kreact(j);
        break;

      case SRA_KK::AA:
      case SRA_KK::DA:
      case SRA_KK::LH1:
      case SRA_KK::LH3:
      case SRA_KK::CD:
        surf_cover = d_total_state[idx] * ms_inv;
        S_theta = 0.0;
        if (d_kisliuk_flag(j)) {
          K_ads = d_kisliuk(j,0) * pow(twall,d_kisliuk(j,1)) *
            exp(-d_kisliuk(j,2)/twall);
          if (surf_cover < 1)
            S_theta = pow((1 - surf_cover) /
                          (1 - surf_cover + K_ads*surf_cover),d_coeff(j,coeff_val));
        } else {
          S_theta = pow((1-surf_cover),d_coeff(j,coeff_val));
        }
        prob_value[i] = d_kreact(j)*S_theta;
        break;

      case SRA_KK::ER:
        {
          double dot = 2.0;
          if (d_nreactant(j) == 1)
            prob_value[i] = 2.0 * d_kreact(j) *
              (maxstick - d_total_state[idx]) * ms_inv / fabs(dot);
          else
            prob_value[i] = 2.0 * d_kreact(j) / fabs(dot);
          break;
        }

      case SRA_KK::CI:
        prob_value[i] = d_kreact(j);
        if (d_energy_flag(j)) {
          double *v = ip->v;
          double dot = v[0]*norm[0]+v[1]*norm[1]+v[2]*norm[2];
          double vmag_sq = v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
          double E_i = 0.5 * d_species[ip->ispecies].mass * vmag_sq;
          double cos_theta = fabs(dot) / sqrt(vmag_sq);
          prob_value[i] *= pow(E_i,d_energy(j,0)) * pow(cos_theta,d_energy(j,1));
        }
        break;
      }

      for (int k = 1; k < d_nreactant(j); k++) {
        if (d_rstate(j,k) == 's') {
          if (d_rpart(j,k) == 0)
            prob_value[i] *= stoich_pow_kk(d_total_state[idx],d_rstoich(j,k)) *
              pow(ms_inv,d_rstoich(j,k));
          else
            prob_value[i] *= stoich_pow_kk(d_species_state(idx,d_rad(j,k)),
                                           d_rstoich(j,k)) *
              pow(ms_inv,d_rstoich(j,k));
        }
      }

      sum_prob += prob_value[i];
    }

    if (sum_prob > 1.0) correction = 1.0/sum_prob;
    else scatter_prob = 1.0 - sum_prob;

    double react_prob = scatter_prob;
    double random_prob = rand_gen.drand();

    if (react_prob > random_prob) {
      rand_pool.free_state(rand_gen);
      return 0;
    }

    for (int i = 0; i < n; i++) {
      int j = d_list(ip->ispecies,i);
      react_prob += prob_value[i] * correction;
      if (react_prob <= random_prob) continue;

      // reaction j fires

      if (ATOMIC_REDUCTION == 0) {
        d_nsingle()++;
        d_tally_single(j)++;
      } else {
        Kokkos::atomic_inc(&d_nsingle());
        Kokkos::atomic_inc(&d_tally_single(j));
      }

      // SURF mode: mark this surf element for the periodic state collate

      if (mode == SRA_KK::SURF) d_mark(idx) = 1;

      // update perspecies deltas for participating surf reactants/products

      auto a_species_delta = d_species_delta;
      for (int k = 0; k < d_nreactant(j); k++)
        if (d_rpart(j,k) == 1 && d_rstate(j,k) == 's')
          Kokkos::atomic_add(&a_species_delta(idx,d_rad(j,k)),-d_rstoich(j,k));
      for (int k = 0; k < d_nproduct(j); k++)
        if (d_ppart(j,k) == 1 && d_pstate(j,k) == 's')
          Kokkos::atomic_add(&a_species_delta(idx,d_pad(j,k)),d_pstoich(j,k));

      // post-reaction particle handling, mirrors SurfReactAdsorb::react()
      // cmodel post-reaction scatter currently supports NOMODEL and SPECULAR
      //   (validated at init); RNG-based cmodels (diffuse/cll/td/...) deferred

      switch (d_type(j)) {

      case SRA_KK::DISSOCIATION:
        {
          double x[3],v[3];
          ip->ispecies = d_products(j,0);
          int id = MAXSMALLINT*rand_gen.drand();
          memcpy(x,ip->x,3*sizeof(double));
          memcpy(v,ip->v,3*sizeof(double));
          int jp_species;
          if (d_pstoich(j,0) == 2) jp_species = d_products(j,0);
          else jp_species = d_products(j,1);

          int index;
          if (ATOMIC_REDUCTION == 0) { index = d_nlocal(); d_nlocal()++; }
          else index = Kokkos::atomic_fetch_add(&d_nlocal(),1);

          int reallocflag = ParticleKokkos::add_particle_kokkos(d_particles,index,
                              id,jp_species,ip->icell,x,v,0.0,0.0);
          if (reallocflag) {
            d_retry() = 1;
            rand_pool.free_state(rand_gen);
            return 0;
          }
          jp = &d_particles[index];
          rand_pool.free_state(rand_gen);
          return (j + 1);
        }

      case SRA_KK::EXCHANGE:
        ip->ispecies = d_products(j,0);
        rand_pool.free_state(rand_gen);
        return (j + 1);

      case SRA_KK::RECOMBINATION:
      case SRA_KK::AA:
      case SRA_KK::LH3:
      case SRA_KK::CD:
        ip = NULL;
        rand_pool.free_state(rand_gen);
        return (j + 1);

      case SRA_KK::DA:
        {
          if (d_nprod_g(j) == 0) ip = NULL;
          else {
            int nn = 1;
            for (int pj = 1; pj < d_nproduct(j); pj++) {
              if (d_pstate(j,pj) == 'g') {
                if (nn == 1) {
                  nn++;
                  ip->ispecies = d_products(j,pj);
                  scatter_cmodel(ip,norm,d_cmodel_ip(j),j,0,rand_gen);
                  if (d_pstoich(j,pj) == 2) {
                    jp = create_particle(ip,d_products(j,pj),rand_gen,d_nlocal,d_retry);
                    if (!jp) { rand_pool.free_state(rand_gen); return 0; }
                    scatter_cmodel(jp,norm,d_cmodel_ip(j),j,0,rand_gen);
                  }
                } else {
                  jp = create_particle(ip,d_products(j,pj),rand_gen,d_nlocal,d_retry);
                  if (!jp) { rand_pool.free_state(rand_gen); return 0; }
                  scatter_cmodel(jp,norm,d_cmodel_jp(j),j,1,rand_gen);
                }
              }
            }
          }
          if (d_cmodel_ip(j) != SRA_KK::NOMODEL) velreset = 1;
          rand_pool.free_state(rand_gen);
          return (j + 1);
        }

      case SRA_KK::LH1:
      case SRA_KK::ER:
        ip->ispecies = d_products(j,0);
        scatter_cmodel(ip,norm,d_cmodel_ip(j),j,0,rand_gen);
        if (d_cmodel_ip(j) != SRA_KK::NOMODEL) velreset = 1;
        rand_pool.free_state(rand_gen);
        return (j + 1);

      case SRA_KK::CI:
        {
          ip->ispecies = d_products(j,0);
          scatter_cmodel(ip,norm,d_cmodel_ip(j),j,0,rand_gen);
          if (d_nprod_g_tot(j) == 2) {
            if (d_pstoich(j,0) == 2) {
              jp = create_particle(ip,d_products(j,0),rand_gen,d_nlocal,d_retry);
              if (!jp) { rand_pool.free_state(rand_gen); return 0; }
              scatter_cmodel(jp,norm,d_cmodel_ip(j),j,0,rand_gen);
            } else {
              jp = create_particle(ip,d_products(j,1),rand_gen,d_nlocal,d_retry);
              if (!jp) { rand_pool.free_state(rand_gen); return 0; }
              scatter_cmodel(jp,norm,d_cmodel_jp(j),j,1,rand_gen);
            }
          }
          if (d_cmodel_ip(j) != SRA_KK::NOMODEL) velreset = 1;
          rand_pool.free_state(rand_gen);
          return (j + 1);
        }
      }
    }

    rand_pool.free_state(rand_gen);
    return 0;
  }

  /* ----------------------------------------------------------------------
     apply a post-reaction collision model (cmodel) scatter to particle p
     SPECULAR mirrors SurfCollideSpecular::wrapper() (reflect, no RNG)
     NOMODEL is a no-op; RNG-based cmodels are rejected at init
  ------------------------------------------------------------------------- */

  KOKKOS_INLINE_FUNCTION
  void scatter_cmodel(Particle::OnePart *p, const double *norm, int cmodel,
                      int j, int useJp, rand_type &adsorb_rand) const
  {
    if (cmodel == SRA_KK::NOMODEL) return;
    if (cmodel == SRA_KK::SPECULAR) {       // SurfCollideSpecular::wrapper
      MathExtraKokkos::reflect3(p->v,norm);
      return;
    }

    // RNG: bit-exact uses the cmodel's own RanKnuth (EXACT serial); on GPU
    //   builds fall back to the surf-react RNG (not bit-exact, not gated here)

#ifdef SPARTA_KOKKOS_EXACT
    rand_type rg = d_cmodel_rand[cmodel];
#else
    rand_type &rg = adsorb_rand;
#endif

    // gather this reaction's cmodel coeffs/flags (ip or jp slot)

    double cf[SRA_KK_MAXCMCOEFF];
    int fl[SRA_KK_MAXCMFLAG];
    for (int k = 0; k < SRA_KK_MAXCMCOEFF; k++)
      cf[k] = useJp ? d_cmjp_coeffs(j,k) : d_cmip_coeffs(j,k);
    for (int k = 0; k < SRA_KK_MAXCMFLAG; k++)
      fl[k] = useJp ? d_cmjp_flags(j,k) : d_cmip_flags(j,k);

    if (cmodel == SRA_KK::DIFFUSE) {
      diffuse_scatter(p,norm,cf[0],cf[1],rg);            // tsurf, acc
    } else if (cmodel == SRA_KK::CLL) {
      double eccen = fl[0] ? cf[5] : 0.0;
      cll_scatter(p,norm,cf[0],cf[1],cf[2],cf[3],cf[4],fl[0],eccen,rg);
    } else if (cmodel == SRA_KK::TD) {
      double tsurf = cf[0];
      int barrier_flag = fl[0], initen_flag = fl[1], bond_flag = fl[2];
      int m = 1;
      double barrier_val = 0.0;
      double initen_trans = 0.0,initen_rot = 0.0,initen_vib = 0.0;
      double bond_trans = 0.0,bond_rot = 0.0,bond_vib = 0.0;
      if (barrier_flag) barrier_val = cf[m++];
      if (initen_flag) { initen_trans = cf[m++]; initen_rot = cf[m++]; initen_vib = cf[m++]; }
      if (bond_flag) { bond_trans = cf[m++]; bond_rot = cf[m++]; bond_vib = cf[m++]; }
      td_scatter(p,norm,tsurf,barrier_flag,barrier_val,initen_flag,initen_trans,
                 initen_rot,initen_vib,bond_flag,bond_trans,bond_rot,bond_vib,rg);
    }
  }

  /* ----------------------------------------------------------------------
     replicas of the Kokkos surf-collide scatter device functions, drawing
     from the cmodel's RNG so they match the host wrapper bit-for-bit;
     cmodels never translate/rotate (trflag off)
  ------------------------------------------------------------------------- */

  KOKKOS_INLINE_FUNCTION
  double erot_kk(int isp, double temp, rand_type &rg) const
  {
    double eng,a,erm,b;
    if (rotstyle_ == SRA_KK::NONE) return 0.0;
    if (d_species[isp].rotdof < 2) return 0.0;
    if (rotstyle_ == SRA_KK::DISCRETE && d_species[isp].rotdof == 2) {
      int irot = -log(rg.drand()) * temp / d_species[isp].rottemp[0];
      eng = irot * boltz_ * d_species[isp].rottemp[0];
    } else if (rotstyle_ == SRA_KK::SMOOTH && d_species[isp].rotdof == 2) {
      eng = -log(rg.drand()) * boltz_ * temp;
    } else {
      a = 0.5*d_species[isp].rotdof - 1.0;
      while (1) {
        erm = 10.0*rg.drand();
        b = pow(erm/a,a) * exp(a-erm);
        if (b > rg.drand()) break;
      }
      eng = erm * boltz_ * temp;
    }
    return eng;
  }

  KOKKOS_INLINE_FUNCTION
  double evib_kk(int isp, double temp, rand_type &rg) const
  {
    double eng,a,erm,b;
    if (vibstyle_ == SRA_KK::NONE || d_species[isp].vibdof < 2) return 0.0;
    eng = 0.0;
    if (vibstyle_ == SRA_KK::DISCRETE && d_species[isp].vibdof == 2) {
      int ivib = -log(rg.drand()) * temp / d_species[isp].vibtemp[0];
      eng = ivib * boltz_ * d_species[isp].vibtemp[0];
    } else if (vibstyle_ == SRA_KK::SMOOTH || d_species[isp].vibdof >= 2) {
      if (d_species[isp].vibdof == 2)
        eng = -log(rg.drand()) * boltz_ * temp;
      else if (d_species[isp].vibdof > 2) {
        a = 0.5*d_species[isp].vibdof - 1.0;
        while (1) {
          erm = 10.0*rg.drand();
          b = pow(erm/a,a) * exp(a-erm);
          if (b > rg.drand()) break;
        }
        eng = erm * boltz_ * temp;
      }
    }
    return eng;
  }

  KOKKOS_INLINE_FUNCTION
  void diffuse_scatter(Particle::OnePart *p, const double *norm,
                       double twall, double acc, rand_type &rg) const
  {
    if (rg.drand() > acc) {
      MathExtraKokkos::reflect3(p->v,norm);
    } else {
      double tangent1[3],tangent2[3];
      int isp = p->ispecies;
      double vrm = sqrt(2.0*boltz_*twall / d_species[isp].mass);
      double vperp = vrm * sqrt(-log(rg.drand()));
      double theta = MathConst::MY_2PI * rg.drand();
      double vtangent = vrm * sqrt(-log(rg.drand()));
      double vtan1 = vtangent * sin(theta);
      double vtan2 = vtangent * cos(theta);
      double *v = p->v;
      double dot = MathExtraKokkos::dot3(v,norm);
      tangent1[0] = v[0] - dot*norm[0];
      tangent1[1] = v[1] - dot*norm[1];
      tangent1[2] = v[2] - dot*norm[2];
      if (MathExtraKokkos::lensq3(tangent1) == 0.0) {
        tangent2[0] = rg.drand();
        tangent2[1] = rg.drand();
        tangent2[2] = rg.drand();
        MathExtraKokkos::cross3(norm,tangent2,tangent1);
      }
      MathExtraKokkos::norm3(tangent1);
      MathExtraKokkos::cross3(norm,tangent1,tangent2);
      v[0] = vperp*norm[0] + vtan1*tangent1[0] + vtan2*tangent2[0];
      v[1] = vperp*norm[1] + vtan1*tangent1[1] + vtan2*tangent2[1];
      v[2] = vperp*norm[2] + vtan1*tangent1[2] + vtan2*tangent2[2];
      p->erot = erot_kk(isp,twall,rg);
      p->evib = evib_kk(isp,twall,rg);
    }
  }

  KOKKOS_INLINE_FUNCTION
  void cll_scatter(Particle::OnePart *p, const double *norm, double twall,
                   double acc_n, double acc_t, double acc_rot, double acc_vib,
                   int pflag, double eccen, rand_type &rg) const
  {
    double tangent1[3],tangent2[3];
    int ispecies = p->ispecies;
    double *v = p->v;
    double dot = MathExtraKokkos::dot3(v,norm);
    double vrm,vperp,vtan1,vtan2;

    tangent1[0] = v[0] - dot*norm[0];
    tangent1[1] = v[1] - dot*norm[1];
    tangent1[2] = v[2] - dot*norm[2];
    if (MathExtraKokkos::lensq3(tangent1) == 0.0) {
      tangent2[0] = rg.drand();
      tangent2[1] = rg.drand();
      tangent2[2] = rg.drand();
      MathExtraKokkos::cross3(norm,tangent2,tangent1);
    }
    MathExtraKokkos::norm3(tangent1);
    MathExtraKokkos::cross3(norm,tangent1,tangent2);
    double tan1 = MathExtraKokkos::dot3(v,tangent1);

    vrm = sqrt(2.0*boltz_ * twall / d_species[ispecies].mass);

    double r_1 = sqrt(-acc_n*log(rg.drand()));
    double theta_1 = MathConst::MY_2PI * rg.drand();
    double dot_norm = dot/vrm * sqrt(1-acc_n);
    vperp = vrm * sqrt(r_1*r_1 + dot_norm*dot_norm + 2*r_1*dot_norm*cos(theta_1));

    double r_2 = sqrt(-acc_t*log(rg.drand()));
    double theta_2 = MathConst::MY_2PI * rg.drand();
    double vtangent = tan1/vrm * sqrt(1-acc_t);
    vtan1 = vrm * (vtangent + r_2*cos(theta_2));
    vtan2 = vrm * r_2 * sin(theta_2);

    if (pflag) {
      double tan2 = MathExtraKokkos::dot3(v,tangent2);
      double phi_i,psi_i,theta_f,phi_f,psi_f,cos_beta;
      psi_i = acos(dot*dot/MathExtraKokkos::lensq3(v));
      phi_i = atan2(tan2,tan1);
      double v_mag = sqrt(vperp*vperp + vtan1*vtan1 + vtan2*vtan2);
      double P = 0;
      while (rg.drand() > P) {
        phi_f = MathConst::MY_2PI*rg.drand();
        psi_f = acos(1-rg.drand());
        cos_beta = cos(psi_i)*cos(psi_f) + sin(psi_i)*sin(psi_f)*cos(phi_i - phi_f);
        P = (1-eccen)/(1-eccen*cos_beta);
      }
      theta_f = acos(sqrt(cos(psi_f)));
      vperp = v_mag * cos(theta_f);
      vtan1 = v_mag * sin(theta_f) * cos(phi_f);
      vtan2 = v_mag * sin(theta_f) * sin(phi_f);
    }

    v[0] = vperp*norm[0] + vtan1*tangent1[0] + vtan2*tangent2[0];
    v[1] = vperp*norm[1] + vtan1*tangent1[1] + vtan2*tangent2[1];
    v[2] = vperp*norm[2] + vtan1*tangent1[2] + vtan2*tangent2[2];

    // rotational component (CLL partial accommodation)

    if (rotstyle_ == SRA_KK::NONE || d_species[ispecies].rotdof < 2) p->erot = 0.0;
    else {
      double erot_mag = sqrt(p->erot*(1-acc_rot)/(boltz_*twall));
      double r_rot,cos_theta_rot,A_rot,X_rot;
      if (d_species[ispecies].rotdof == 2) {
        r_rot = sqrt(-acc_rot*log(rg.drand()));
        cos_theta_rot = cos(MathConst::MY_2PI*rg.drand());
      } else {
        A_rot = 0;
        while (A_rot < rg.drand()) {
          X_rot = 4*rg.drand();
          A_rot = 2.71828182845904523536028747*X_rot*X_rot*exp(-X_rot*X_rot);
        }
        r_rot = sqrt(acc_rot)*X_rot;
        cos_theta_rot = 2*rg.drand() - 1;
      }
      p->erot = boltz_ * twall *
        (r_rot*r_rot + erot_mag*erot_mag + 2*r_rot*erot_mag*cos_theta_rot);
    }

    // vibrational component

    int vibdof = d_species[ispecies].vibdof;
    double r_vib,cos_theta_vib,A_vib,X_vib,evib_mag,evib_val;
    if (vibstyle_ == SRA_KK::NONE || vibdof < 2) p->evib = 0.0;
    else if (vibstyle_ == SRA_KK::DISCRETE && vibdof == 2) {
      double evib_star = -log(1 - rg.drand() *
        (1 - exp(-boltz_*d_species[ispecies].vibtemp[0])));
      evib_val = p->evib + evib_star;
      evib_mag = sqrt(evib_val*(1-acc_vib)/(boltz_*twall));
      r_vib = sqrt(-acc_vib*log(rg.drand()));
      cos_theta_vib = cos(MathConst::MY_2PI*rg.drand());
      evib_val = boltz_ * twall *
        (r_vib*r_vib + evib_mag*evib_mag + 2*r_vib*evib_mag*cos_theta_vib);
      int ivib = evib_val / (boltz_*d_species[ispecies].vibtemp[0]);
      p->evib = ivib * boltz_ * d_species[ispecies].vibtemp[0];
    }
    else if (vibstyle_ == SRA_KK::SMOOTH || vibdof >= 2) {
      evib_mag = sqrt(p->evib*(1-acc_vib)/(boltz_*twall));
      if (vibdof == 2) {
        r_vib = sqrt(-acc_vib*log(rg.drand()));
        cos_theta_vib = cos(MathConst::MY_2PI*rg.drand());
      } else {
        A_vib = 0;
        while (A_vib < rg.drand()) {
          X_vib = 4*rg.drand();
          A_vib = 2.71828182845904523536028747*X_vib*X_vib*exp(-X_vib*X_vib);
        }
        r_vib = sqrt(acc_vib)*X_vib;
        cos_theta_vib = 2*rg.drand() - 1;
      }
      p->evib = boltz_ * twall *
        (r_vib*r_vib + evib_mag*evib_mag + 2*r_vib*evib_mag*cos_theta_vib);
    }
  }

  KOKKOS_INLINE_FUNCTION
  void td_scatter(Particle::OnePart *p, const double *norm, double twall,
                  int barrier_flag, double barrier_val,
                  int initen_flag, double initen_trans, double initen_rot, double initen_vib,
                  int bond_flag, double bond_trans, double bond_rot, double bond_vib,
                  rand_type &rg) const
  {
    double tangent1[3],tangent2[3];
    int ispecies = p->ispecies;
    double *v = p->v;
    double dot = MathExtraKokkos::dot3(v,norm);

    tangent1[0] = v[0] - dot*norm[0];
    tangent1[1] = v[1] - dot*norm[1];
    tangent1[2] = v[2] - dot*norm[2];
    if (MathExtraKokkos::lensq3(tangent1) == 0.0) {
      tangent2[0] = rg.drand();
      tangent2[1] = rg.drand();
      tangent2[2] = rg.drand();
      MathExtraKokkos::cross3(norm,tangent2,tangent1);
    }
    MathExtraKokkos::norm3(tangent1);
    MathExtraKokkos::cross3(norm,tangent1,tangent2);

    double mass = d_species[ispecies].mass;
    double E_i = 0.5 * mass * MathExtraKokkos::lensq3(v);
    double E_t = boltz_ * twall;
    if (bond_flag) E_t += boltz_*bond_trans;
    if (initen_flag) E_t += E_i*initen_trans;
    double E_n = E_t;
    if (barrier_flag) E_n += boltz_*barrier_val;

    double vrm_n = sqrt(2.0*E_n / mass);
    double vrm_t = sqrt(2.0*E_t / mass);
    double vperp = vrm_n * sqrt(-log(rg.drand()));
    double theta = MathConst::MY_2PI * rg.drand();
    double vtangent = vrm_t * sqrt(-log(rg.drand()));
    double vtan1 = vtangent * sin(theta);
    double vtan2 = vtangent * cos(theta);

    v[0] = vperp*norm[0] + vtan1*tangent1[0] + vtan2*tangent2[0];
    v[1] = vperp*norm[1] + vtan1*tangent1[1] + vtan2*tangent2[1];
    v[2] = vperp*norm[2] + vtan1*tangent1[2] + vtan2*tangent2[2];

    double twall_rot = twall, twall_vib = twall;
    if (bond_flag) { twall_rot += bond_rot; twall_vib += bond_vib; }
    if (initen_flag) { twall_rot += E_i*initen_rot/boltz_; twall_vib += E_i*initen_vib/boltz_; }

    p->erot = erot_kk(ispecies,twall_rot,rg);
    p->evib = evib_kk(ispecies,twall_vib,rg);
  }

  /* ----------------------------------------------------------------------
     create a new particle (copy of ip's x,v) of species sp, mirrors the
     add_particle path in SurfReactAdsorb::react(); returns ptr or NULL on
     realloc (caller retries the whole move)
  ------------------------------------------------------------------------- */

  KOKKOS_INLINE_FUNCTION
  Particle::OnePart *create_particle(Particle::OnePart *ip, int sp,
                                     rand_type &rand_gen,
                                     const DAT::t_int_scalar &d_nlocal,
                                     const DAT::t_int_scalar &d_retry) const
  {
    double x[3],v[3];
    memcpy(x,ip->x,3*sizeof(double));
    memcpy(v,ip->v,3*sizeof(double));
    int id = MAXSMALLINT*rand_gen.drand();

    int index = Kokkos::atomic_fetch_add(&d_nlocal(),1);
    int reallocflag = ParticleKokkos::add_particle_kokkos(d_particles,index,
                        id,sp,ip->icell,x,v,0.0,0.0);
    if (reallocflag) {
      d_retry() = 1;
      return NULL;
    }
    return &d_particles[index];
  }

  KOKKOS_INLINE_FUNCTION
  double stoich_pow_kk(int base, int p) const
  {
    const double THIRD = 1.0/3.0;
    switch (p) {
    case 0: return 1.0;
    case 1: return (base >= p) ? double(base) : 0.0;
    case 2: return (base >= p) ? 0.5*base*(base-1) : 0.0;
    case 3: return (base >= p) ? 0.5*THIRD*base*(base-1)*(base-2) : 0.0;
    case 4: return (base >= p) ? 0.125*THIRD*base*(base-1)*(base-2)*(base-3) : 0.0;
    case 5: return (base >= p) ? 0.025*THIRD*base*(base-1)*(base-2)*(base-3)*(base-4) : 0.0;
    case 6: return (base >= p) ? 0.0125*THIRD*THIRD*base*(base-1)*(base-2)*(base-3)*(base-4)*(base-5) : 0.0;
    }
    return 0.0;
  }
};

}

#endif
#endif

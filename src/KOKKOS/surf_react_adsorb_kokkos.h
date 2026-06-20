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
}

#define SRA_KK_MAXREACTANT 5
#define SRA_KK_MAXPRODUCT 5
#define SRA_KK_MAXCOEFF 4
#define SRA_KK_MAXPERSPECIES 16   // max GS reactions a single species can be in

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

  // per-face state (FACE mode); small (nface <= 6)

  DAT::t_int_1d d_total_state;             // [nface]
  DAT::t_float_1d d_area,d_weight;         // [nface]
  DAT::t_int_2d d_species_state;           // [nface][nspecies_surf]
  DAT::t_int_2d d_species_delta;           // [nface][nspecies_surf] (atomic)

  DAT::tdual_int_2d k_species_delta;

  double fnum_;        // update->fnum, set in pre_react

  void init_reactions_gs_kokkos();

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
    // convert face index from negative value to 0..5 inclusive

    int iface = -(isurf+1);

    int n = d_reactions_n[ip->ispecies];
    if (n == 0) return 0;

    double fnum = fnum_;
    long int maxstick = ceil(max_cover*d_area[iface] / (fnum*d_weight[iface]));
    double factor = fnum * d_weight[iface] / d_area[iface];
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
        surf_cover = d_total_state[iface] * ms_inv;
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
              (maxstick - d_total_state[iface]) * ms_inv / fabs(dot);
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
            prob_value[i] *= stoich_pow_kk(d_total_state[iface],d_rstoich(j,k)) *
              pow(ms_inv,d_rstoich(j,k));
          else
            prob_value[i] *= stoich_pow_kk(d_species_state(iface,d_rad(j,k)),
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

      // update per-face perspecies deltas for participating surf reactants/products

      auto a_species_delta = d_species_delta;
      for (int k = 0; k < d_nreactant(j); k++)
        if (d_rpart(j,k) == 1 && d_rstate(j,k) == 's')
          Kokkos::atomic_add(&a_species_delta(iface,d_rad(j,k)),-d_rstoich(j,k));
      for (int k = 0; k < d_nproduct(j); k++)
        if (d_ppart(j,k) == 1 && d_pstate(j,k) == 's')
          Kokkos::atomic_add(&a_species_delta(iface,d_pad(j,k)),d_pstoich(j,k));

      // post-reaction particle handling
      // only types validated at init (no cmodel post-scatter) are reachable

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
      }
    }

    rand_pool.free_state(rand_gen);
    return 0;
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

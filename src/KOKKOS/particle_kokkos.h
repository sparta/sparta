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

#ifndef SPARTA_PARTICLE_KOKKOS_H
#define SPARTA_PARTICLE_KOKKOS_H

#include "stdio.h"
#include "particle.h"
#include "kokkos_type.h"
#include "Kokkos_Random.hpp"
#include "rand_pool_wrap.h"

namespace SPARTA_NS {

struct TagParticleZero_cellcount{};
struct TagParticleCompressReactions{};
struct TagCopyParticleReorderDestinations{};
struct TagFixedMemoryReorder{};
struct TagFixedMemoryReorderInit{};
struct TagSetIcellFromPlist{};
struct TagParticleReorder_COPYPARTICLELIST{};
struct TagSetDPlistNewStyle{};

template<int NEED_ATOMICS>
struct TagParticleSort{};


class ParticleKokkos : public Particle {
 public:
  typedef int value_type;

  // methods

  ParticleKokkos(class SPARTA *);
  ~ParticleKokkos();
  static KOKKOS_INLINE_FUNCTION
  int add_particle_kokkos(t_particle_1d particles, int, int, int, int,
                           double *, double *, double, double);
#ifndef SPARTA_KOKKOS_EXACT
  void compress_migrate(int, int *);
#endif
  void sort_kokkos();
  void grow(int);
  void grow_species();
  void pre_weight() override;
  void post_weight() override;
  void update_class_variables();

#ifndef SPARTA_KOKKOS_EXACT
  typedef typename Kokkos::Random_XorShift64_Pool<DeviceType>::generator_type rand_type;

  //typedef typename Kokkos::Random_XorShift1024_Pool<DeviceType>::generator_type rand_type;
#else
  typedef RandWrap rand_type;
#endif

  KOKKOS_INLINE_FUNCTION
  double erot(int, double, rand_type &) const;

  KOKKOS_INLINE_FUNCTION
  double evib(int, double, rand_type &) const;

  void wrap_kokkos();
  void sync(ExecutionSpace, unsigned int);
  void modify(ExecutionSpace, unsigned int);

  KOKKOS_INLINE_FUNCTION
  void operator()(TagParticleZero_cellcount, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagParticleCompressReactions, const int&) const;

  template<int NEED_ATOMICS>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagParticleSort<NEED_ATOMICS>, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagParticleReorder_COPYPARTICLELIST, const int, int&, const bool&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagCopyParticleReorderDestinations, const int, int&, const bool&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixedMemoryReorder, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixedMemoryReorderInit, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagSetIcellFromPlist, const int&) const;


  tdual_particle_1d k_particles;
  tdual_species_1d k_species;
  DAT::tdual_int_2d k_species2group;

  int sorted_kk;

  inline
  int get_maxcellcount() {
    return maxcellcount;
  }

  inline
  void set_maxcellcount(int in) {
    maxcellcount = in;
  }

 private:
  t_particle_1d d_particles;
  t_particle_1d d_sorted;
  t_species_1d d_species;
  int nParticlesWksp;
  DAT::tdual_int_scalar k_reorder_pass;
  DAT::t_int_scalar d_reorder_pass;
  HAT::t_int_scalar h_reorder_pass;

  int nbytes;
  int maxcellcount,ngrid;
  int collide_rot,vibstyle;
  double boltz;

  DAT::t_int_2d d_plist;
  DAT::t_int_1d d_cellcount;

  DAT::t_int_2d d_lists;
  DAT::t_int_1d d_mlist;
  DAT::t_int_1d d_slist;

  HAT::t_int_2d h_lists;
  HAT::t_int_1d h_mlist;
  HAT::t_int_1d h_slist;

  DAT::t_int_scalar d_fail_flag;
  HAT::t_int_scalar h_fail_flag;

  // work memory for reduced memory reordering
  t_particle_1d d_pswap1;
  t_particle_1d d_pswap2;
};

KOKKOS_INLINE_FUNCTION
int ParticleKokkos::add_particle_kokkos(t_particle_1d particles, int index, int id,
      int ispecies, int icell, double *x, double *v, double erot, double evib)
{
  OnePart tmp;
  tmp.id = id;
  tmp.ispecies = ispecies;
  tmp.icell = icell;
  tmp.x[0] = x[0];
  tmp.x[1] = x[1];
  tmp.x[2] = x[2];
  tmp.v[0] = v[0];
  tmp.v[1] = v[1];
  tmp.v[2] = v[2];
  tmp.erot = erot;
  tmp.evib = evib;
  enum{PKEEP,PINSERT,PDONE,PDISCARD,PENTRY,PEXIT,PSURF};  // same as .cpp file
  tmp.flag = PKEEP;

  int realloc = 0;

  if (index < particles.extent(0)) {
    particles[index] = tmp;
  } else {
    realloc = 1;
  }

  return realloc; 
}


/* ----------------------------------------------------------------------
   generate random rotational energy for a particle
   only a function of species index and species properties
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
double ParticleKokkos::erot(int isp, double temp_thermal, rand_type &erandom) const
{
 double eng,a,erm,b;

 if (!collide_rot) return 0.0;
 if (d_species[isp].rotdof < 2) return 0.0;

 if (d_species[isp].rotdof == 2)
   eng = -log(erandom.drand()) * boltz * temp_thermal;
 else {
   a = 0.5*d_species[isp].rotdof-1.0;
   while (1) {
     // energy cut-off at 10 kT
     erm = 10.0*erandom.drand();
     b = pow(erm/a,a) * exp(a-erm);
     if (b > erandom.drand()) break;
   }
   eng = erm * boltz * temp_thermal;
 }

 return eng;
}

/* ----------------------------------------------------------------------
   generate random vibrational energy for a particle
   only a function of species index and species properties
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
double ParticleKokkos::evib(int isp, double temp_thermal, rand_type &erandom) const
{
  double eng,a,erm,b;

  enum{NONE,DISCRETE,SMOOTH};            // several files
  if (vibstyle == NONE || d_species[isp].vibdof < 2) return 0.0;

  eng = 0.0;
  if (vibstyle == DISCRETE && d_species[isp].vibdof == 2) {
    int ivib = static_cast<int> (-log(erandom.drand()) * temp_thermal / 
                                 d_species[isp].vibtemp[0]);
    eng = ivib * boltz * d_species[isp].vibtemp[0];
  } else if (vibstyle == SMOOTH || d_species[isp].vibdof >= 2) {
    if (d_species[isp].vibdof == 2)
      eng = -log(erandom.drand()) * boltz * temp_thermal;
    else if (d_species[isp].vibdof > 2) {
      a = 0.5*d_species[isp].vibdof-1.;
      while (1) {
        // energy cut-off at 10 kT
        erm = 10.0*erandom.drand();
        b = pow(erm/a,a) * exp(a-erm);
        if (b > erandom.drand()) break;
      }
      eng = erm * boltz * temp_thermal;
    }
  }

  return eng;
}

}

#endif

/* ERROR/WARNING messages:

*/

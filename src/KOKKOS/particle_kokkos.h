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

#ifndef SPARTA_PARTICLE_KOKKOS_H
#define SPARTA_PARTICLE_KOKKOS_H

#include "stdio.h"
#include "particle.h"
#include "kokkos_type.h"
#include "Kokkos_Random.hpp"
#include "rand_pool_wrap.h"

namespace SPARTA_NS {

struct TagParticleCompressReactions{};
struct TagCopyParticleReorderDestinations{};
struct TagFixedMemoryReorder{};
struct TagFixedMemoryReorderInit{};
struct TagSetIcellFromPlist{};
struct TagParticleReorder_COPYPARTICLELIST1{};
struct TagParticleReorder_COPYPARTICLELIST2{};
struct TagSetDPlistNewStyle{};

template<int NEED_ATOMICS, int REORDER_FLAG>
struct TagParticleSort{};


class ParticleKokkos : public Particle {
 public:
  typedef int value_type;

  // methods

  ParticleKokkos(class SPARTA *);
  ~ParticleKokkos() override;
  static KOKKOS_INLINE_FUNCTION
  int add_particle_kokkos(t_particle_1d particles, int, int, int, int,
                           double *, double *, double, double);
#ifndef SPARTA_KOKKOS_EXACT
  void compress_migrate(int, int *) override;
#endif
  void sort_kokkos();
  void grow(int) override;
  void grow_species() override;
  void pre_weight() override;
  void post_weight() override;
  void update_class_variables();
  int add_custom(char *, int, int) override;
  void grow_custom(int, int, int) override;
  void remove_custom(int) override;
  void copy_custom(int, int) override;
  void pack_custom(int, char *) override;
  void unpack_custom(char *, int) override;

  KOKKOS_INLINE_FUNCTION
  void copy_custom_kokkos(int, int) const;

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

  KOKKOS_INLINE_FUNCTION
  void pack_custom_kokkos(int, char *) const;

  KOKKOS_INLINE_FUNCTION
  void unpack_custom_kokkos(char *, int) const;

  void wrap_kokkos();
  void sync(ExecutionSpace, unsigned int);
  void modify(ExecutionSpace, unsigned int);

  KOKKOS_INLINE_FUNCTION
  void operator()(TagParticleCompressReactions, const int&) const;

  template<int NEED_ATOMICS, int REORDER_FLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagParticleSort<NEED_ATOMICS,REORDER_FLAG>, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagParticleReorder_COPYPARTICLELIST1, const int, int&, const bool&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagParticleReorder_COPYPARTICLELIST2, const int) const;

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

  DAT::tdual_int_1d k_ewhich,k_eicol,k_edcol;

  tdual_struct_tdual_int_1d_1d k_eivec;
  tdual_struct_tdual_float_1d_1d k_edvec;

  tdual_struct_tdual_int_2d_1d k_eiarray;
  tdual_struct_tdual_float_2d_1d k_edarray;

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
  t_species_1d d_species;

  t_particle_1d d_sorted;
  DAT::t_int_1d d_sorted_id;
  DAT::t_int_1d d_offsets_part;
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

  DAT::t_int_2d_lr d_lists;
  DAT::t_int_1d d_mlist;
  DAT::t_int_1d d_slist;

  HAT::t_int_2d_lr h_lists;
  HAT::t_int_1d h_mlist;
  HAT::t_int_1d h_slist;

  DAT::t_int_scalar d_resize;
  HAT::t_int_scalar h_resize;

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

KOKKOS_INLINE_FUNCTION
void ParticleKokkos::pack_custom_kokkos(int n, char *buf) const
{
  int i,j;
  char *ptr = buf;

  if (ncustom_ivec) {
    for (i = 0; i < ncustom_ivec; i++) {
      memcpy(ptr,&(k_eivec.d_view(i).k_view.d_view(n)),sizeof(int));
      ptr += sizeof(int);
    }
  }
  if (ncustom_iarray) {
    for (i = 0; i < ncustom_iarray; i++) {
      const int ncols = k_eicol.d_view[i];
      for (j = 0; j < ncols; j++) {
        memcpy(ptr,&(k_eiarray.d_view(i).k_view.d_view(n,j)),sizeof(int));
        ptr += sizeof(int);
      }
    }
  }

  ptr = ROUNDUP(ptr);

  if (ncustom_dvec) {
    for (i = 0; i < ncustom_dvec; i++) {
      memcpy(ptr,&(k_edvec.d_view(i).k_view.d_view(n)),sizeof(double));
      ptr += sizeof(double);
    }
  }
  if (ncustom_darray) {
    for (i = 0; i < ncustom_darray; i++) {
      const int ncols = k_edcol.d_view[i];
      for (j = 0; j < ncols; j++) {
        memcpy(ptr,&(k_edarray.d_view(i).k_view.d_view(n,j)),sizeof(double));
        ptr += sizeof(double);
      }
    }
  }
}

KOKKOS_INLINE_FUNCTION
void ParticleKokkos::unpack_custom_kokkos(char *buf, int n) const
{
  int i,j;
  char *ptr = buf;

  if (ncustom_ivec) {
    for (i = 0; i < ncustom_ivec; i++) {
      memcpy(&(k_eivec.d_view(i).k_view.d_view(n)),ptr,sizeof(int));
      ptr += sizeof(int);
    }
  }
  if (ncustom_iarray) {
    for (i = 0; i < ncustom_iarray; i++) {
      const int ncols = k_eicol.d_view[i];
      for (j = 0; j < ncols; j++) {
        memcpy(&(k_eiarray.d_view(i).k_view.d_view(n,j)),ptr,sizeof(int));
        ptr += sizeof(int);
      }
    }
  }

  ptr = ROUNDUP(ptr);

  if (ncustom_dvec) {
    for (i = 0; i < ncustom_dvec; i++) {
      memcpy(&(k_edvec.d_view(i).k_view.d_view(n)),ptr,sizeof(double));
      ptr += sizeof(double);
    }
  }
  if (ncustom_darray) {
    for (i = 0; i < ncustom_darray; i++) {
      const int ncols = k_edcol.d_view[i];
      for (j = 0; j < ncols; j++) {
        memcpy(&(k_edarray.d_view(i).k_view.d_view(n,j)),ptr,sizeof(double));
        ptr += sizeof(double);
      }
    }
  }
}


}

#endif

/* ERROR/WARNING messages:

*/

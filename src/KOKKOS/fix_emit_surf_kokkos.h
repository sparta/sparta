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

#ifdef FIX_CLASS

FixStyle(emit/surf/kk,FixEmitSurfKokkos)

#else

#ifndef SPARTA_FIX_EMIT_SURF_KOKKOS_H
#define SPARTA_FIX_EMIT_SURF_KOKKOS_H

#include "fix_emit_surf.h"
#include "rand_pool_wrap.h"
#include "kokkos_copy.h"
#include "particle_kokkos.h"
#include "compute_surf_kokkos.h"
#include "kokkos_base.h"
#include "region_prim_kokkos.h"

namespace SPARTA_NS {

struct TagFixEmitSurf_ninsert{};
struct TagFixEmitSurf_perform_task{};
struct TagFixEmitSurf_subsonic_inflow{};
struct TagFixEmitSurf_subsonic_grid{};
struct TagFixEmitSurf_mflow_grid{};
struct TagFixEmitSurf_mflow_nrho{};

template<int ATOMIC_REDUCTION>
struct TagFixEmitSurf_insert_particles{};

class FixEmitSurfKokkos : public FixEmitSurf {
 public:
  typedef int value_type;

  FixEmitSurfKokkos(class SPARTA *, int, char **);
  ~FixEmitSurfKokkos() override;
  void init() override;
  void flatten_region();
  void perform_task() override;
  void perform_task_twopass() override { perform_task(); }

  void grid_changed() override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixEmitSurf_ninsert, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixEmitSurf_perform_task, const int&, int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixEmitSurf_subsonic_inflow, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixEmitSurf_subsonic_grid, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixEmitSurf_mflow_grid, const int&, double&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixEmitSurf_mflow_nrho, const int&) const;

  template<int ATOMIC_REDUCTION>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixEmitSurf_insert_particles<ATOMIC_REDUCTION>, const int&) const;

#ifndef SPARTA_KOKKOS_EXACT
  Kokkos::Random_XorShift64_Pool<DeviceType> rand_pool;
  typedef typename Kokkos::Random_XorShift64_Pool<DeviceType>::generator_type rand_type;

  //Kokkos::Random_XorShift1024_Pool<DeviceType> rand_pool;
  //typedef typename Kokkos::Random_XorShift1024_Pool<DeviceType>::generator_type rand_type;
#else
  RandPoolWrap rand_pool;
  typedef RandWrap rand_type;
#endif

 private:
  int nsurf_tally,nlocal_before,nlocal_surf,region_flag;
  double npcurrent;       // matches the non-Kokkos double (VARIABLE npmode)
  int plist_descending;   // 1 if the host walks d_plist high index -> low
  double boltz,temp_thermal_mix;
  double acoef,nrho_mflow;

  KKCopy<ParticleKokkos> particle_kk_copy;
  // the active compute surf tallies this fix drives, reached on device only
  //   through FES_SLIST().  Same two layouts as UpdateKokkos's tally lists:
  //   a fixed KKCopy array held by value in the functor, or one runtime-sized
  //   device byte buffer with no cap.  See kokkos_type.h for why both exist.

#ifdef SPARTA_KOKKOS_FIXED_LISTS
  KKCopy<ComputeSurfKokkos> slist_active_copy[KOKKOS_MAX_SLIST];
#define FES_SLIST(m) slist_active_copy[m].obj
#else
  DAT::tdual_char_1d k_slist_surf;
  DAT::t_char_1d d_slist_surf;
#define FES_SLIST(m) ((const ComputeSurfKokkos *) d_slist_surf.data())[m]
#endif
  // region flattened to a device-resident postfix token stream; replaces
  //   the per-style KKCopy members and the caps that went with them.
  //   region_flag says whether there is a region at all -- nregion_token and
  //   d_region_tokens are only meaningful when it is 1

  tdual_region_token_1d k_region_tokens;
  t_region_token_1d d_region_tokens;
  int nregion_token;

  typedef Kokkos::DualView<Task*, DeviceType::array_layout, DeviceType> tdual_task_1d;
  typedef tdual_task_1d::t_dev t_task_1d;
  tdual_task_1d k_tasks;
  t_task_1d d_tasks;

  DAT::tdual_float_2d_lr k_ntargetsp;          // # of mols to insert for each species
  DAT::tdual_float_2d_lr k_vscale;             // vscale for each species
  DAT::tdual_float_2d_lr k_path;               // path for each species
  DAT::tdual_float_2d_lr k_fracarea;           // fracarea for each species
  DAT::t_float_2d_lr d_ntargetsp;
  DAT::t_float_2d_lr d_vscale;
  DAT::t_float_2d_lr d_path;
  DAT::t_float_2d_lr d_fracarea;

  DAT::tdual_int_1d k_ninsert;
  Kokkos::View<int*, DeviceType> d_ninsert; // won't compile with DAT::t_int_1d type
  DAT::t_int_1d d_task2cand;
  DAT::t_int_1d d_cands2new;

  DAT::t_float_2d d_x;
  DAT::t_float_2d d_v;
  DAT::t_float_1d d_erot;
  DAT::t_float_1d d_evib;
  DAT::t_float_1d d_dtremain;
  DAT::t_int_1d   d_id;
  DAT::t_int_1d   d_isp;
  DAT::t_int_1d   d_task;
  Kokkos::View<int*, DeviceType> d_keep; // won't compile with DAT::t_int_1d type

  t_particle_1d d_particles;

  DAT::tdual_float_1d k_vscale_mix;
  DAT::tdual_float_1d k_cummulative_mix;
  DAT::tdual_float_2d_lr k_cummulative_custom;
  DAT::tdual_int_1d k_mspecies;

  DAT::t_float_1d d_vscale_mix;
  DAT::t_float_1d d_cummulative_mix;
  DAT::t_float_2d_lr d_cummulative_custom;
  DAT::t_int_1d d_mspecies;

  DAT::tdual_float_1d k_fraction;        // mixture fraction for each species
  DAT::t_float_1d d_fraction;

  t_line_1d d_lines;
  t_tri_1d d_tris;

  // data structs for subsonic emission

  t_particle_1d d_subsonic_particles;
  t_species_1d d_species_all;            // all particle species (mass, rotdof)
  t_cinfo_1d d_cinfo;
  DAT::t_int_2d d_plist;
  DAT::t_int_1d d_cellcount;
  DAT::t_float_scalar d_tempmax;

  void create_tasks() override;
  void grow_task() override;
  void realloc_nspecies() override;

  void subsonic_inflow() override;
  void subsonic_sort() override;
  void subsonic_grid() override;
  void mflow_grid() override;

#ifdef SPARTA_KOKKOS_FIXED_LISTS
  // unused fixed slots must not alias a compute that may be reallocated or
  //   deleted while they still reference count it
  ComputeSurfKokkos tmp_compute_surf_kk;
#endif
};

}

#endif
#endif

/* ERROR/WARNING messages:
 */


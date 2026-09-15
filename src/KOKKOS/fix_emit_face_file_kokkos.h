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

FixStyle(emit/face/file/kk,FixEmitFaceFileKokkos)

#else

#ifndef SPARTA_FIX_EMIT_FACE_FILE_KOKKOS_H
#define SPARTA_FIX_EMIT_FACE_FILE_KOKKOS_H

#include "fix_emit_face_file.h"
#include "rand_pool_wrap.h"
#include "kokkos_base.h"
#include "kokkos_copy.h"
#include "particle_kokkos.h"
#include "region_prim_kokkos.h"

namespace SPARTA_NS {

struct TagFixEmitFaceFile_ninsert{};
struct TagFixEmitFaceFile_perform_task{};
struct TagFixEmitFaceFile_subsonic_inflow{};
struct TagFixEmitFaceFile_subsonic_grid{};

class FixEmitFaceFileKokkos : public FixEmitFaceFile {
 public:
  typedef int value_type;

  FixEmitFaceFileKokkos(class SPARTA *, int, char **);
  ~FixEmitFaceFileKokkos() override;
  void init() override;

  // the Kokkos path is always two-pass: the count scan has to know every
  //   task's insertion count before candidate arrays can be sized.  so both
  //   entry points land on the same kernel pair, exactly as fix emit/face/kk
  //   does (fix_emit_face_kokkos.h:46-47)

  void flatten_region();
  void perform_task() override;
  void perform_task_twopass() override { perform_task(); }

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixEmitFaceFile_ninsert, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixEmitFaceFile_perform_task, const int&, int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixEmitFaceFile_subsonic_inflow, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixEmitFaceFile_subsonic_grid, const int&) const;

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
  int region_flag;
  int axisymmetric;      // copy of domain->axisymmetric, needed on device
                         //   FixEmitFaceFile, unlike FixEmitFace, keeps no
                         //   copy of it
  double boltz;
  double dt_step;        // update->dt for the current step, used for dtremain
                         //   kept separate from the base member dt, which the
                         //   non-Kokkos code freezes at init() and uses for
                         //   the subsonic ntarget recalculation

  KKCopy<ParticleKokkos> particle_kk_copy;

  // region flattened to a device-resident postfix token stream, so the
  //   insertion kernel needs no virtual dispatch and no typed copy per
  //   region style.  the stream carries each sub-region's interior/exterior
  //   sense and the composite's own, so nothing else needs to be passed
  //   along.  region_flag says whether there is a region at all --
  //   nregion_token and d_region_tokens are only meaningful when it is 1.
  //   see region_prim_kokkos.h

  tdual_region_token_1d k_region_tokens;
  t_region_token_1d d_region_tokens;
  int nregion_token;

  typedef Kokkos::DualView<Task*, DeviceType::array_layout, DeviceType> tdual_task_1d;
  typedef tdual_task_1d::t_dev t_task_1d;
  tdual_task_1d k_tasks;
  t_task_1d d_tasks;

  // per-task, per-species arrays.
  // unlike fix emit/face, fix emit/face/file carries a per-task fraction,
  //   cummulative and vscale ALWAYS -- they are interpolated from the file
  //   per face, not taken mixture-wide -- so all three are allocated
  //   unconditionally, and the kernels index them (i,isp) rather than (isp).
  // the host Task::fraction/cummulative/vscale/ntargetsp pointers are aimed
  //   at the host rows of these DualViews by realloc_species_views() and
  //   grow_task().  that is what lets the unmodified host interpolate() go
  //   on writing through those pointers: the file mesh setup stays entirely
  //   host-side and is flattened once, at task-build time, not per step.

  DAT::tdual_float_2d_lr k_ntargetsp;   // # of mols to insert for each species
  DAT::tdual_float_2d_lr k_vscale;      // vscale for each species
  DAT::tdual_float_2d_lr k_cummulative; // cummulative fraction for each species
  DAT::tdual_float_2d_lr k_fraction;    // fraction for each species
  DAT::t_float_2d_lr d_ntargetsp;
  DAT::t_float_2d_lr d_vscale;
  DAT::t_float_2d_lr d_cummulative;
  DAT::t_float_2d_lr d_fraction;

  Kokkos::View<int*, DeviceType> d_ninsert;
  DAT::t_int_1d d_task2cand;

  DAT::t_float_2d d_x;
  DAT::t_float_1d d_beta_un;
  DAT::t_float_1d d_theta;
  DAT::t_float_1d d_vr;
  DAT::t_float_1d d_erot;
  DAT::t_float_1d d_evib;
  DAT::t_float_1d d_dtremain;
  DAT::t_int_1d   d_id;
  DAT::t_int_1d   d_isp;
  DAT::t_int_1d   d_task;
  Kokkos::View<int*, DeviceType> d_keep; // won't compile with DAT::t_int_1d type

  DAT::tdual_int_1d k_mspecies;          // species indices of mixture
  DAT::t_int_1d d_mspecies;

  // data structs for subsonic emission

  t_particle_1d d_particles;
  t_species_1d d_species_all;            // all particle species (mass, rotdof)
  t_cinfo_1d d_cinfo;
  DAT::t_int_2d d_plist;
  DAT::t_int_1d d_cellcount;
  DAT::t_float_scalar d_tempmax;
  int plist_descending;   // 1 if the host walks d_plist high index -> low

  void create_tasks() override;
  void grow_task() override;

  void subsonic_inflow() override;
  void subsonic_sort() override;
  void subsonic_grid() override;

  // (re)size the per-task species DualViews and re-aim every Task pointer at
  //   their host rows.  FixEmitFaceFile has no realloc_nspecies() hook the
  //   way FixEmitFace does, so this is called explicitly from init(), before
  //   the base init() reaches create_tasks()

  void realloc_species_views();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Fix emit/face/file/kk requires the twopass keyword under SPARTA_KOKKOS_EXACT

The Kokkos insertion loop draws every task's insertion count before it
generates any particle.  The non-Kokkos one-pass loop interleaves the two.
Only the non-Kokkos twopass keyword consumes random numbers in the same
order, so only then can the two runs produce the same particles.

*/

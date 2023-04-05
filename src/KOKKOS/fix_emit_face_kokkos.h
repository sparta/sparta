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

#ifdef FIX_CLASS

FixStyle(emit/face/kk,FixEmitFaceKokkos)

#else

#ifndef SPARTA_FIX_EMIT_FACE_KOKKOS_H
#define SPARTA_FIX_EMIT_FACE_KOKKOS_H

#include "fix_emit_face.h"
#include "rand_pool_wrap.h"
#include "kokkos_copy.h"
#include "particle_kokkos.h"

namespace SPARTA_NS {

struct TagFixEmitFace_ninsert{};
struct TagFixEmitFace_perform_task{};

class FixEmitFaceKokkos : public FixEmitFace {
 public:
  typedef int value_type;

  FixEmitFaceKokkos(class SPARTA *, int, char **);
  ~FixEmitFaceKokkos() override;
  void init() override;
  void perform_task() override;
  void perform_task_twopass() override { perform_task(); }

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixEmitFace_ninsert, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixEmitFace_perform_task, const int&, int&) const;

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
  KKCopy<ParticleKokkos> particle_kk_copy;

  typedef Kokkos::DualView<Task*, DeviceType::array_layout, DeviceType> tdual_task_1d;
  typedef tdual_task_1d::t_dev t_task_1d;
  tdual_task_1d k_tasks;
  t_task_1d d_tasks;

  DAT::tdual_float_2d_lr k_ntargetsp;          // # of mols to insert for each species
  DAT::tdual_float_2d_lr k_vscale;             // vscale for each species
  DAT::t_float_2d_lr d_ntargetsp;
  DAT::t_float_2d_lr d_vscale;

  DAT::tdual_int_1d k_ninsert;
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

  DAT::tdual_float_1d k_mix_vscale;
  DAT::tdual_float_1d k_cummulative;
  DAT::tdual_int_1d k_species;

  DAT::t_float_1d d_mix_vscale;
  DAT::t_float_1d d_cummulative;
  DAT::t_int_1d d_species;

  void create_task(int) override;
  void grow_task() override;
  void realloc_nspecies() override;
};

}

#endif
#endif

/* ERROR/WARNING messages:
 */


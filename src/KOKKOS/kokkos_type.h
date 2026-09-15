/* -*- c++ -*- ----------------------------------------------------------
   SPARTA - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://sparta.github.io, Sandia National Laboratories
   Steve Plimpton, sjplimp@gmail.com

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifndef SPARTA_STYPE_KOKKOS_H
#define SPARTA_STYPE_KOKKOS_H

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <Kokkos_ScatterView.hpp>

#include "particle.h"
#include "grid.h"
#include "surf.h"
#include "spatype.h"
#include "accelerator_kokkos_defs.h"

#include <cstring>

// offset type for the Kokkos::Crs per-cell surf/split/sub lists.
// Under BIGBIG the total number of flattened entries on one rank can exceed
//   2^31, so the row_map offsets must be a bigint.  Under BIG they cannot, and
//   a 32-bit row_map halves the bytes touched by the per-particle surf
//   lookups in the move kernel.  All declarations of these Crs objects must
//   use this typedef so the device and host mirrors stay assignable.

#ifdef SPARTA_BIGBIG
typedef SPARTA_NS::bigint crs_size_type;
#else
typedef int crs_size_type;
#endif

#define NeighClusterSize 8

// SPARTA_KOKKOS_FIXED_LISTS restores the original fixed-size KKCopy arrays for
//   the per-type tally compute lists and for the per-style surf react lists,
//   in place of the runtime-sized device buffers that replaced them.  The
//   buffers exist to lift the instance caps below; the fixed arrays keep every
//   compute and surf react model inside the functor that is handed by value to
//   each kernel.
// Which is faster is a GPU question with no obvious answer: a smaller functor
//   can raise occupancy while costing data locality, and higher occupancy is
//   not the same thing as higher throughput.  Neither path has been measured
//   on an accelerator.  Both are kept so the two can be compared on real
//   hardware by rebuilding with -DSPARTA_KOKKOS_FIXED_LISTS, rather than by
//   reverting commits.

#ifdef SPARTA_KOKKOS_FIXED_LISTS
#define KOKKOS_MAX_SLIST 2
#define KOKKOS_MAX_BLIST 2
#define KOKKOS_MAX_GLIST 4
#define KOKKOS_MAX_SURF_REACT_PER_TYPE 2
#define KOKKOS_MAX_TOT_SURF_REACT 4
#endif

// architectures where the move kernel is dispatched with ATOMIC_REDUCTION = -1
//   (parallel_reduce) rather than atomics.  KokkosSPARTA::accelerator() clears
//   atomic_reduction for exactly these, and UpdateKokkos::move() reads the
//   per-step counters back from the reduction result only when it is clear.
//   The two decisions must agree or the counters come from the wrong place, so
//   the condition is spelled out once here instead of in both files.

#if defined(KOKKOS_ARCH_AMD_GFX940) || defined(KOKKOS_ARCH_AMD_GFX942) || \
    defined(KOKKOS_ARCH_AMD_GFX942_APU)
#define SPARTA_KOKKOS_REDUCE_ARCH 1
#else
#define SPARTA_KOKKOS_REDUCE_ARCH 0
#endif

// the active surf react models, as named by the eight classes that dispatch to
//   them: compute surf, and the seven surf collide models that support surface
//   chemistry.  Both representations are reached through these accessors, so
//   the device dispatch site in each of those classes is written exactly once
//   and the two modes cannot drift:
//
//     KK_SR_*    read on device, from the model's collide_kokkos()
//     KK_SR_H_*  the host image of the same models, used by the
//                pre_react()/post_react()/backup()/restore() lifecycle
//
//   under SPARTA_KOKKOS_FIXED_LISTS both are the same fixed KKCopy arrays, held
//   by value in the class; otherwise the device side reads the per-style device
//   buffer and the host side the host half of the same DualView.
// the accessors expand to member accesses only, and all eight classes spell
//   those members identically, so one definition here serves every one of them.

#ifdef SPARTA_KOKKOS_FIXED_LISTS

#define KK_SR_TYPE(n)     sr_type_list[n]
#define KK_SR_MAP(n)      sr_map[n]
#define KK_SR_GLOBAL(m)   sr_kk_global_copy[m].obj
#define KK_SR_PROB(m)     sr_kk_prob_copy[m].obj
#define KK_SR_ADSORB(m)   sr_kk_adsorb_copy[m].obj

#define KK_SR_H_TYPE(n)   sr_type_list[n]
#define KK_SR_H_MAP(n)    sr_map[n]
#define KK_SR_H_GLOBAL(m) sr_kk_global_copy[m].obj
#define KK_SR_H_PROB(m)   sr_kk_prob_copy[m].obj
#define KK_SR_H_ADSORB(m) sr_kk_adsorb_copy[m].obj

#else

#define KK_SR_TYPE(n)     d_sr_type_list[n]
#define KK_SR_MAP(n)      d_sr_map[n]
#define KK_SR_GLOBAL(m)   ((const SurfReactGlobalKokkos *) d_sr_global.data())[m]
#define KK_SR_PROB(m)     ((const SurfReactProbKokkos *) d_sr_prob.data())[m]
#define KK_SR_ADSORB(m)   ((const SurfReactAdsorbKokkos *) d_sr_adsorb.data())[m]

#define KK_SR_H_TYPE(n)   k_sr_type_list.view_host()[n]
#define KK_SR_H_MAP(n)    k_sr_map.view_host()[n]
#define KK_SR_H_GLOBAL(m) ((SurfReactGlobalKokkos *) k_sr_global.view_host().data())[m]
#define KK_SR_H_PROB(m)   ((SurfReactProbKokkos *) k_sr_prob.view_host().data())[m]
#define KK_SR_H_ADSORB(m) ((SurfReactAdsorbKokkos *) k_sr_adsorb.view_host().data())[m]

#endif

namespace Kokkos {
  static auto NoInit = [](std::string const& label) {
    return Kokkos::view_alloc(Kokkos::WithoutInitializing, label);
  };
}

  struct sparta_float3 {
    float x,y,z;
    KOKKOS_INLINE_FUNCTION
    sparta_float3():x(0.0f),y(0.0f),z(0.0f) {}

    KOKKOS_INLINE_FUNCTION
    void operator += (const sparta_float3& tmp) {
      x+=tmp.x;
      y+=tmp.y;
      z+=tmp.z;
    }
    KOKKOS_INLINE_FUNCTION
    void operator = (const sparta_float3& tmp) {
      x=tmp.x;
      y=tmp.y;
      z=tmp.z;
    }
  };

  struct sparta_double3 {
    double x,y,z;
    KOKKOS_INLINE_FUNCTION
    sparta_double3():x(0.0),y(0.0),z(0.0) {}

    KOKKOS_INLINE_FUNCTION
    void operator += (const sparta_double3& tmp) {
      x+=tmp.x;
      y+=tmp.y;
      z+=tmp.z;
    }
    KOKKOS_INLINE_FUNCTION
    void operator = (const sparta_double3& tmp) {
      x=tmp.x;
      y=tmp.y;
      z=tmp.z;
    }
  };

// set SPAHostype and DeviceType from Kokkos Default Types
typedef Kokkos::DefaultExecutionSpace SPADeviceType;
typedef Kokkos::HostSpace::execution_space SPAHostType;

typedef SPADeviceType DeviceType;

// set ExecutionSpace stuct with variable "space"

template<class Device>
struct ExecutionSpaceFromDevice;

template<>
struct ExecutionSpaceFromDevice<SPAHostType> {
  static const SPARTA_NS::ExecutionSpace space = SPARTA_NS::Host;
};

#ifdef KOKKOS_ENABLE_CUDA
template<>
struct ExecutionSpaceFromDevice<Kokkos::Cuda> {
  static const SPARTA_NS::ExecutionSpace space = SPARTA_NS::Device;
};
#elif defined(KOKKOS_ENABLE_HIP)
template<>
struct ExecutionSpaceFromDevice<Kokkos::Experimental::HIP> {
  static const SPARTA_NS::ExecutionSpace space = SPARTA_NS::Device;
};
#elif defined(KOKKOS_ENABLE_SYCL)
template<>
struct ExecutionSpaceFromDevice<Kokkos::Experimental::SYCL> {
  static const SPARTA_NS::ExecutionSpace space = SPARTA_NS::Device;
};
#elif defined(KOKKOS_ENABLE_OPENMPTARGET)
template<>
struct ExecutionSpaceFromDevice<Kokkos::Experimental::OpenMPTarget> {
  static const SPARTA_NS::ExecutionSpace space = SPARTA_NS::Device;
};
#endif

// set host pinned space
#if defined(KOKKOS_ENABLE_CUDA)
typedef Kokkos::CudaHostPinnedSpace SPAPinnedHostType;
#elif defined(KOKKOS_ENABLE_HIP)
typedef Kokkos::Experimental::HIPHostPinnedSpace SPAPinnedHostType;
#elif defined(KOKKOS_ENABLE_SYCL)
typedef Kokkos::Experimental::SYCLHostUSMSpace SPAPinnedHostType;
#elif defined(KOKKOS_ENABLE_OPENMPTARGET)
typedef Kokkos::Serial SPAPinnedHostType;
#endif

// Determine memory traits for atomic arrays
template<int NEED_ATOMICS>
struct AtomicView {
  enum {value = Kokkos::Unmanaged};
};

template<>
struct AtomicView<1> {
  enum {value = Kokkos::Atomic|Kokkos::Unmanaged};
};

template<>
struct AtomicView<-1> {
  enum {value = Kokkos::Atomic|Kokkos::Unmanaged};
};

// Determine memory traits for array
// Do atomic trait when running with CUDA
template<int NEED_ATOMICS, class DeviceType>
struct AtomicDup {
  using value = Kokkos::Experimental::ScatterNonAtomic;
};

#ifdef KOKKOS_ENABLE_CUDA
template<>
struct AtomicDup<1,Kokkos::Cuda> {
  using value = Kokkos::Experimental::ScatterAtomic;
};

template<>
struct AtomicDup<-1,Kokkos::Cuda> {
  using value = Kokkos::Experimental::ScatterAtomic;
};
#elif defined(KOKKOS_ENABLE_HIP)
template<>
struct AtomicDup<1,Kokkos::Experimental::HIP> {
  using value = Kokkos::Experimental::ScatterAtomic;
};

template<>
struct AtomicDup<-1,Kokkos::Experimental::HIP> {
  using value = Kokkos::Experimental::ScatterAtomic;
};
#elif defined(KOKKOS_ENABLE_SYCL)
template<>
struct AtomicDup<1,Kokkos::Experimental::SYCL> {
  using value = Kokkos::Experimental::ScatterAtomic;
};

template<>
struct AtomicDup<-1,Kokkos::Experimental::SYCL> {
  using value = Kokkos::Experimental::ScatterAtomic;
};
#elif defined(KOKKOS_ENABLE_OPENMPTARGET)
template<>
struct AtomicDup<1,Kokkos::Experimental::OpenMPTarget> {
  using value = Kokkos::Experimental::ScatterAtomic;
};

template<>
struct AtomicDup<-1,Kokkos::Experimental::OpenMPTarget> {
  using value = Kokkos::Experimental::ScatterAtomic;
};
#endif

#ifdef SPARTA_KOKKOS_USE_ATOMICS

#ifdef KOKKOS_ENABLE_OPENMP
template<>
struct AtomicDup<1,Kokkos::OpenMP> {
  using value = Kokkos::Experimental::ScatterAtomic;
};

template<>
struct AtomicDup<-1,Kokkos::OpenMP> {
  using value = Kokkos::Experimental::ScatterAtomic;
};
#endif

#ifdef KOKKOS_ENABLE_THREADS
template<>
struct AtomicDup<1,Kokkos::Threads> {
  using value = Kokkos::Experimental::ScatterAtomic;
};

template<>
struct AtomicDup<-1,Kokkos::Threads> {
  using value = Kokkos::Experimental::ScatterAtomic;
};
#endif

#endif


// Determine duplication traits for array
// Use duplication when running threaded and not using atomics
template<int NEED_ATOMICS, class DeviceType>
struct NeedDup {
  using value = Kokkos::Experimental::ScatterNonDuplicated;
};

#ifndef SPARTA_KOKKOS_USE_ATOMICS

#ifdef KOKKOS_ENABLE_OPENMP
template<>
struct NeedDup<1,Kokkos::OpenMP> {
  using value = Kokkos::Experimental::ScatterDuplicated;
};

template<>
struct NeedDup<-1,Kokkos::OpenMP> {
  using value = Kokkos::Experimental::ScatterDuplicated;
};
#endif

#ifdef KOKKOS_ENABLE_THREADS
template<>
struct NeedDup<1,Kokkos::Threads> {
  using value = Kokkos::Experimental::ScatterDuplicated;
};

template<>
struct NeedDup<-1,Kokkos::Threads> {
  using value = Kokkos::Experimental::ScatterDuplicated;
};
#endif

#endif

template<typename value, typename T1, typename T2>
class ScatterViewHelper {};

template<typename T1, typename T2>
class ScatterViewHelper<Kokkos::Experimental::ScatterDuplicated,T1,T2> {
public:
  KOKKOS_INLINE_FUNCTION
  static T1 get(const T1 &dup, const T2 & /*nondup*/) {
    return dup;
  }
};

template<typename T1, typename T2>
class ScatterViewHelper<Kokkos::Experimental::ScatterNonDuplicated,T1,T2> {
public:
  KOKKOS_INLINE_FUNCTION
  static T2 get(const T1 & /*dup*/, const T2 &nondup) {
    return nondup;
  }
};


// define precision

#ifndef SPA_PRECISION
#define SPA_PRECISION 2
#endif
#if SPA_PRECISION==1
typedef float SPARTA_FLOAT;
#else
typedef double SPARTA_FLOAT;
#endif

#ifndef PREC_FORCE
#define PREC_FORCE SPA_PRECISION
#endif

#if PREC_FORCE==1
typedef float F_FLOAT;
#else
typedef double F_FLOAT;
#endif

#ifndef PREC_ENERGY
#define PREC_ENERGY SPA_PRECISION
#endif

#if PREC_ENERGY==1
typedef float E_FLOAT;
#else
typedef double E_FLOAT;
#endif

struct s_EV_FLOAT {
  E_FLOAT evdwl;
  E_FLOAT ecoul;
  E_FLOAT v[6];
  KOKKOS_INLINE_FUNCTION
  s_EV_FLOAT() {
    evdwl = 0;
    ecoul = 0;
    v[0] = 0; v[1] = 0; v[2] = 0;
    v[3] = 0; v[4] = 0; v[5] = 0;
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(const s_EV_FLOAT &rhs) {
    evdwl += rhs.evdwl;
    ecoul += rhs.ecoul;
    v[0] += rhs.v[0];
    v[1] += rhs.v[1];
    v[2] += rhs.v[2];
    v[3] += rhs.v[3];
    v[4] += rhs.v[4];
    v[5] += rhs.v[5];
  }
};
typedef struct s_EV_FLOAT EV_FLOAT;

#ifndef PREC_POS
#define PREC_POS SPA_PRECISION
#endif

#if PREC_POS==1
typedef float X_FLOAT;
#else
typedef double X_FLOAT;
#endif

#ifndef PREC_VELOCITIES
#define PREC_VELOCITIES SPA_PRECISION
#endif

#if PREC_VELOCITIES==1
typedef float V_FLOAT;
#else
typedef double V_FLOAT;
#endif

#if PREC_KSPACE==1
typedef float K_FLOAT;
#else
typedef double K_FLOAT;
#endif

// ------------------------------------------------------------------------

// SPARTA types

namespace SPARTA_NS {

  typedef Kokkos::
    DualView<Particle::OnePart*, DeviceType::array_layout, DeviceType> tdual_particle_1d;
  typedef tdual_particle_1d::t_dev t_particle_1d;
  typedef tdual_particle_1d::t_host t_host_particle_1d;

  typedef Kokkos::
    DualView<Particle::OnePart**, DeviceType::array_layout, DeviceType> tdual_particle_2d;
  typedef tdual_particle_2d::t_dev t_particle_2d;
  typedef tdual_particle_2d::t_host t_host_particle_2d;

  typedef Kokkos::
    DualView<Particle::Species*, DeviceType::array_layout, DeviceType> tdual_species_1d;
  typedef tdual_species_1d::t_dev t_species_1d;
  typedef tdual_species_1d::t_dev_const t_species_1d_const;
  typedef tdual_species_1d::t_host t_host_species_1d;

  typedef Kokkos::
    DualView<Grid::ChildCell*, DeviceType::array_layout, DeviceType> tdual_cell_1d;
  typedef tdual_cell_1d::t_dev t_cell_1d;
  typedef tdual_cell_1d::t_host t_host_cell_1d;

  typedef Kokkos::
    DualView<Grid::ChildInfo*, DeviceType::array_layout, DeviceType> tdual_cinfo_1d;
  typedef tdual_cinfo_1d::t_dev t_cinfo_1d;
  typedef tdual_cinfo_1d::t_host t_host_cinfo_1d;

  typedef Kokkos::
    DualView<Grid::SplitInfo*, DeviceType::array_layout, DeviceType> tdual_sinfo_1d;
  typedef tdual_sinfo_1d::t_dev t_sinfo_1d;
  typedef tdual_sinfo_1d::t_host t_host_sinfo_1d;

  typedef Kokkos::
    DualView<Grid::ParentCell*, DeviceType::array_layout, DeviceType> tdual_pcell_1d;
  typedef tdual_pcell_1d::t_dev t_pcell_1d;
  typedef tdual_pcell_1d::t_host t_host_pcell_1d;

  typedef Kokkos::
    DualView<Grid::ParentLevel*, DeviceType::array_layout, DeviceType> tdual_plevel_1d;
  typedef tdual_pcell_1d::t_dev t_plevel_1d;
  typedef tdual_pcell_1d::t_host t_host_plevel_1d;

  typedef Kokkos::
    DualView<Surf::Line*, DeviceType::array_layout, DeviceType> tdual_line_1d;
  typedef tdual_line_1d::t_dev t_line_1d;
  typedef tdual_line_1d::t_host t_host_line_1d;

  typedef Kokkos::
    DualView<Surf::Tri*, DeviceType::array_layout, DeviceType> tdual_tri_1d;
  typedef tdual_tri_1d::t_dev t_tri_1d;
  typedef tdual_tri_1d::t_host t_host_tri_1d;

  // device-callable equivalent of Compute::ubuf, used by the KOKKOS tally
  //   computes to pack ints into a double buf slot from within
  //   KOKKOS_INLINE_FUNCTION methods.  Compute::ubuf's constructors are
  //   host-only, which nvcc tolerates but hipcc rejects, so the bit pattern
  //   is reproduced here as a KOKKOS_INLINE_FUNCTION union instead.

  union d_ubuf {
    double d;
    int64_t i;
    KOKKOS_INLINE_FUNCTION d_ubuf(double arg) : d(arg) {}
    KOKKOS_INLINE_FUNCTION d_ubuf(int64_t arg) : i(arg) {}
    KOKKOS_INLINE_FUNCTION d_ubuf(int arg) : i(arg) {}
    KOKKOS_INLINE_FUNCTION d_ubuf(uint32_t arg) : i(arg) {}
    KOKKOS_INLINE_FUNCTION d_ubuf(uint64_t arg) : i(arg) {}
  };
}

template <class DeviceType>
struct ArrayTypes;

template <>
struct ArrayTypes<DeviceType> {

// scalar types

typedef Kokkos::
  DualView<int, DeviceType::array_layout, DeviceType> tdual_int_scalar;
typedef tdual_int_scalar::t_dev t_int_scalar;
typedef tdual_int_scalar::t_dev_const t_int_scalar_const;
typedef tdual_int_scalar::t_dev_um t_int_scalar_um;
typedef tdual_int_scalar::t_dev_const_um t_int_scalar_const_um;

typedef Kokkos::
  DualView<SPARTA_NS::bigint, DeviceType::array_layout, DeviceType> tdual_bigint_scalar;
typedef tdual_bigint_scalar::t_dev t_bigint_scalar;
typedef tdual_bigint_scalar::t_dev_const t_bigint_scalar_const;
typedef tdual_bigint_scalar::t_dev_um t_bigint_scalar_um;
typedef tdual_bigint_scalar::t_dev_const_um t_bigint_scalar_const_um;

typedef Kokkos::
  DualView<SPARTA_FLOAT, DeviceType::array_layout, DeviceType>
  tdual_float_scalar;
typedef tdual_float_scalar::t_dev t_float_scalar;
typedef tdual_float_scalar::t_dev_const t_float_scalar_const;
typedef tdual_float_scalar::t_dev_um t_float_scalar_um;
typedef tdual_float_scalar::t_dev_const_um t_float_scalar_const_um;

// generic array types

typedef Kokkos::
  DualView<char*, DeviceType::array_layout, DeviceType> tdual_char_1d;
typedef tdual_char_1d::t_dev t_char_1d;
typedef tdual_char_1d::t_dev_const t_char_1d_const;
typedef tdual_char_1d::t_dev_um t_char_1d_um;
typedef tdual_char_1d::t_dev_const_um t_char_1d_const_um;
typedef tdual_char_1d::t_dev_const_randomread t_char_1d_randomread;

typedef Kokkos::
  DualView<int*, DeviceType::array_layout, DeviceType> tdual_int_1d;
typedef tdual_int_1d::t_dev t_int_1d;
typedef tdual_int_1d::t_dev_const t_int_1d_const;
typedef tdual_int_1d::t_dev_um t_int_1d_um;
typedef tdual_int_1d::t_dev_const_um t_int_1d_const_um;
typedef tdual_int_1d::t_dev_const_randomread t_int_1d_randomread;

typedef Kokkos::
  DualView<SPARTA_NS::bigint*, DeviceType::array_layout, DeviceType> tdual_bigint_1d;
typedef tdual_bigint_1d::t_dev t_bigint_1d;
typedef tdual_bigint_1d::t_dev_const t_bigint_1d_const;
typedef tdual_bigint_1d::t_dev_um t_bigint_1d_um;
typedef tdual_bigint_1d::t_dev_const_um t_bigint_1d_const_um;
typedef tdual_bigint_1d::t_dev_const_randomread t_bigint_1d_randomread;

typedef Kokkos::
  DualView<int*[3], DeviceType::array_layout, DeviceType> tdual_int_1d_3;
typedef tdual_int_1d_3::t_dev t_int_1d_3;
typedef tdual_int_1d_3::t_dev_const t_int_1d_3_const;
typedef tdual_int_1d_3::t_dev_um t_int_1d_3_um;
typedef tdual_int_1d_3::t_dev_const_um t_int_1d_3_const_um;
typedef tdual_int_1d_3::t_dev_const_randomread t_int_1d_3_randomread;

typedef Kokkos::
  DualView<int**, Kokkos::LayoutRight, DeviceType> tdual_int_2d_lr;
typedef tdual_int_2d_lr::t_dev t_int_2d_lr;
typedef tdual_int_2d_lr::t_dev_const t_int_2d_const_lr;
typedef tdual_int_2d_lr::t_dev_um t_int_2d_um_lr;
typedef tdual_int_2d_lr::t_dev_const_um t_int_2d_const_um_lr;
typedef tdual_int_2d_lr::t_dev_const_randomread t_int_2d_randomread_lr;

typedef Kokkos::
  DualView<int**, DeviceType::array_layout, DeviceType> tdual_int_2d;
typedef tdual_int_2d::t_dev t_int_2d;
typedef tdual_int_2d::t_dev_const t_int_2d_const;
typedef tdual_int_2d::t_dev_um t_int_2d_um;
typedef tdual_int_2d::t_dev_const_um t_int_2d_const_um;
typedef tdual_int_2d::t_dev_const_randomread t_int_2d_randomread;

typedef Kokkos::
  DualView<SPARTA_NS::cellint*, DeviceType::array_layout, DeviceType>
  tdual_cellint_1d;
typedef tdual_cellint_1d::t_dev t_cellint_1d;
typedef tdual_cellint_1d::t_dev_const t_cellint_1d_const;
typedef tdual_cellint_1d::t_dev_um t_cellint_1d_um;
typedef tdual_cellint_1d::t_dev_const_um t_cellint_1d_const_um;
typedef tdual_cellint_1d::t_dev_const_randomread t_cellint_1d_randomread;

typedef Kokkos::
  DualView<SPARTA_NS::surfint*, DeviceType::array_layout, DeviceType>
  tdual_surfint_1d;
typedef tdual_surfint_1d::t_dev t_surfint_1d;
typedef tdual_surfint_1d::t_dev_const t_surfint_1d_const;
typedef tdual_surfint_1d::t_dev_um t_surfint_1d_um;
typedef tdual_surfint_1d::t_dev_const_um t_surfint_1d_const_um;
typedef tdual_surfint_1d::t_dev_const_randomread t_surfint_1d_randomread;

// 1d float array n

typedef Kokkos::DualView<SPARTA_FLOAT*, DeviceType::array_layout, DeviceType> tdual_float_1d;
typedef tdual_float_1d::t_dev t_float_1d;
typedef tdual_float_1d::t_dev_const t_float_1d_const;
typedef tdual_float_1d::t_dev_um t_float_1d_um;
typedef tdual_float_1d::t_dev_const_um t_float_1d_const_um;
typedef tdual_float_1d::t_dev_const_randomread t_float_1d_randomread;

//1d float array strided
typedef Kokkos::DualView<SPARTA_FLOAT*, Kokkos::LayoutStride, DeviceType> tdual_float_1d_strided;
typedef tdual_float_1d_strided::t_dev t_float_1d_strided;
typedef tdual_float_1d_strided::t_dev_um t_float_1d_strided_um;

// 1d float array n[3]

typedef Kokkos::DualView<SPARTA_FLOAT*[3], DeviceType::array_layout, DeviceType> tdual_float_1d_3;
typedef tdual_float_1d_3::t_dev t_float_1d_3;
typedef tdual_float_1d_3::t_dev_const t_float_1d_3_const;
typedef tdual_float_1d_3::t_dev_um t_float_1d_3_um;
typedef tdual_float_1d_3::t_dev_const_um t_float_1d_3_const_um;
typedef tdual_float_1d_3::t_dev_const_randomread t_float_1d_3_randomread;

//2d float array n
typedef Kokkos::DualView<SPARTA_FLOAT**, DeviceType::array_layout, DeviceType> tdual_float_2d;
typedef tdual_float_2d::t_dev t_float_2d;
typedef tdual_float_2d::t_dev_const t_float_2d_const;
typedef tdual_float_2d::t_dev_um t_float_2d_um;
typedef tdual_float_2d::t_dev_const_um t_float_2d_const_um;
typedef tdual_float_2d::t_dev_const_randomread t_float_2d_randomread;

//2d float array n Kokkos::LayoutRight
typedef Kokkos::DualView<F_FLOAT**, Kokkos::LayoutRight, DeviceType> tdual_float_2d_lr;
typedef tdual_float_2d_lr::t_dev t_float_2d_lr;
typedef tdual_float_2d_lr::t_dev_const t_float_2d_lr_const;
typedef tdual_float_2d_lr::t_dev_um t_float_2d_lr_um;
typedef tdual_float_2d_lr::t_dev_const_um t_float_2d_lr_const_um;
typedef tdual_float_2d_lr::t_dev_const_randomread t_float_2d_lr_randomread;

//3d float array n
typedef Kokkos::DualView<SPARTA_FLOAT***, DeviceType::array_layout, DeviceType> tdual_float_3d;
typedef tdual_float_3d::t_dev t_float_3d;
typedef tdual_float_3d::t_dev_const t_float_3d_const;
typedef tdual_float_3d::t_dev_um t_float_3d_um;
typedef tdual_float_3d::t_dev_const_um t_float_3d_const_um;
typedef tdual_float_3d::t_dev_const_randomread t_float_3d_randomread;
};

#ifdef SPARTA_KOKKOS_GPU
template <>
struct ArrayTypes<SPAHostType> {

//Scalar Types

typedef Kokkos::DualView<int, DeviceType::array_layout, DeviceType> tdual_int_scalar;
typedef tdual_int_scalar::t_host t_int_scalar;
typedef tdual_int_scalar::t_host_const t_int_scalar_const;
typedef tdual_int_scalar::t_host_um t_int_scalar_um;
typedef tdual_int_scalar::t_host_const_um t_int_scalar_const_um;

typedef Kokkos::DualView<SPARTA_NS::bigint, DeviceType::array_layout, DeviceType> tdual_bigint_scalar;
typedef tdual_bigint_scalar::t_host t_bigint_scalar;
typedef tdual_bigint_scalar::t_host_const t_bigint_scalar_const;
typedef tdual_bigint_scalar::t_host_um t_bigint_scalar_um;
typedef tdual_bigint_scalar::t_host_const_um t_bigint_scalar_const_um;

typedef Kokkos::DualView<SPARTA_FLOAT, DeviceType::array_layout, DeviceType> tdual_float_scalar;
typedef tdual_float_scalar::t_host t_float_scalar;
typedef tdual_float_scalar::t_host_const t_float_scalar_const;
typedef tdual_float_scalar::t_host_um t_float_scalar_um;
typedef tdual_float_scalar::t_host_const_um t_float_scalar_const_um;

//Generic ArrayTypes
typedef Kokkos::
  DualView<char*, DeviceType::array_layout, DeviceType> tdual_char_1d;
typedef tdual_char_1d::t_host t_char_1d;
typedef tdual_char_1d::t_host_const t_char_1d_const;
typedef tdual_char_1d::t_host_um t_char_1d_um;
typedef tdual_char_1d::t_host_const_um t_char_1d_const_um;
typedef tdual_char_1d::t_host_const_randomread t_char_1d_randomread;

typedef Kokkos::DualView<int*, DeviceType::array_layout, DeviceType> tdual_int_1d;
typedef tdual_int_1d::t_host t_int_1d;
typedef tdual_int_1d::t_host_const t_int_1d_const;
typedef tdual_int_1d::t_host_um t_int_1d_um;
typedef tdual_int_1d::t_host_const_um t_int_1d_const_um;
typedef tdual_int_1d::t_host_const_randomread t_int_1d_randomread;

typedef Kokkos::DualView<SPARTA_NS::bigint*, DeviceType::array_layout, DeviceType> tdual_bigint_1d;
typedef tdual_bigint_1d::t_host t_bigint_1d;
typedef tdual_bigint_1d::t_host_const t_bigint_1d_const;
typedef tdual_bigint_1d::t_host_um t_bigint_1d_um;
typedef tdual_bigint_1d::t_host_const_um t_bigint_1d_const_um;
typedef tdual_bigint_1d::t_host_const_randomread t_bigint_1d_randomread;

typedef Kokkos::DualView<int*[3], DeviceType::array_layout, DeviceType> tdual_int_1d_3;
typedef tdual_int_1d_3::t_host t_int_1d_3;
typedef tdual_int_1d_3::t_host_const t_int_1d_3_const;
typedef tdual_int_1d_3::t_host_um t_int_1d_3_um;
typedef tdual_int_1d_3::t_host_const_um t_int_1d_3_const_um;
typedef tdual_int_1d_3::t_host_const_randomread t_int_1d_3_randomread;

typedef Kokkos::DualView<int**, Kokkos::LayoutRight, DeviceType> tdual_int_2d_lr;
typedef tdual_int_2d_lr::t_host t_int_2d_lr;
typedef tdual_int_2d_lr::t_host_const t_int_2d_const_lr;
typedef tdual_int_2d_lr::t_host_um t_int_2d_um_lr;
typedef tdual_int_2d_lr::t_host_const_um t_int_2d_const_um_lr;
typedef tdual_int_2d_lr::t_host_const_randomread t_int_2d_randomread_lr;

typedef Kokkos::DualView<int**, DeviceType::array_layout, DeviceType> tdual_int_2d;
typedef tdual_int_2d::t_host t_int_2d;
typedef tdual_int_2d::t_host_const t_int_2d_const;
typedef tdual_int_2d::t_host_um t_int_2d_um;
typedef tdual_int_2d::t_host_const_um t_int_2d_const_um;
typedef tdual_int_2d::t_host_const_randomread t_int_2d_randomread;

typedef Kokkos::DualView<SPARTA_NS::cellint*, DeviceType::array_layout, DeviceType> tdual_cellint_1d;
typedef tdual_cellint_1d::t_host t_cellint_1d;
typedef tdual_cellint_1d::t_host_const t_cellint_1d_const;
typedef tdual_cellint_1d::t_host_um t_cellint_1d_um;
typedef tdual_cellint_1d::t_host_const_um t_cellint_1d_const_um;
typedef tdual_cellint_1d::t_host_const_randomread t_cellint_1d_randomread;

typedef Kokkos::DualView<SPARTA_NS::surfint*, DeviceType::array_layout, DeviceType> tdual_surfint_1d;
typedef tdual_surfint_1d::t_host t_surfint_1d;
typedef tdual_surfint_1d::t_host_const t_surfint_1d_const;
typedef tdual_surfint_1d::t_host_um t_surfint_1d_um;
typedef tdual_surfint_1d::t_host_const_um t_surfint_1d_const_um;
typedef tdual_surfint_1d::t_host_const_randomread t_surfint_1d_randomread;

//1d float array
typedef Kokkos::DualView<SPARTA_FLOAT*, DeviceType::array_layout, DeviceType> tdual_float_1d;
typedef tdual_float_1d::t_host t_float_1d;
typedef tdual_float_1d::t_host_const t_float_1d_const;
typedef tdual_float_1d::t_host_um t_float_1d_um;
typedef tdual_float_1d::t_host_const_um t_float_1d_const_um;
typedef tdual_float_1d::t_host_const_randomread t_float_1d_randomread;

//1d float array strided
typedef Kokkos::DualView<SPARTA_FLOAT*, Kokkos::LayoutStride, DeviceType> tdual_float_1d_strided;
typedef tdual_float_1d_strided::t_host t_float_1d_strided;
typedef tdual_float_1d_strided::t_host_um t_float_1d_strided_um;

//1d float array n[3]
typedef Kokkos::DualView<SPARTA_FLOAT*[3], DeviceType::array_layout, DeviceType> tdual_float_1d_3;
typedef tdual_float_1d_3::t_host t_float_1d_3;
typedef tdual_float_1d_3::t_host_const t_float_1d_3_const;
typedef tdual_float_1d_3::t_host_um t_float_1d_3_um;
typedef tdual_float_1d_3::t_host_const_um t_float_1d_3_const_um;
typedef tdual_float_1d_3::t_host_const_randomread t_float_1d_3_randomread;

//2d float array
typedef Kokkos::DualView<SPARTA_FLOAT**, DeviceType::array_layout, DeviceType> tdual_float_2d;
typedef tdual_float_2d::t_host t_float_2d;
typedef tdual_float_2d::t_host_const t_float_2d_const;
typedef tdual_float_2d::t_host_um t_float_2d_um;
typedef tdual_float_2d::t_host_const_um t_float_2d_const_um;
typedef tdual_float_2d::t_host_const_randomread t_float_2d_randomread;

//2d float array LayoutRight
typedef Kokkos::DualView<F_FLOAT**, Kokkos::LayoutRight, DeviceType> tdual_float_2d_lr;
typedef tdual_float_2d_lr::t_host t_float_2d_lr;
typedef tdual_float_2d_lr::t_host_const t_float_2d_lr_const;
typedef tdual_float_2d_lr::t_host_um t_float_2d_lr_um;
typedef tdual_float_2d_lr::t_host_const_um t_float_2d_lr_const_um;
typedef tdual_float_2d_lr::t_host_const_randomread t_float_2d_lr_randomread;

//3d float array
typedef Kokkos::DualView<SPARTA_FLOAT***, DeviceType::array_layout, DeviceType> tdual_float_3d;
typedef tdual_float_3d::t_host t_float_3d;
typedef tdual_float_3d::t_host_const t_float_3d_const;
typedef tdual_float_3d::t_host_um t_float_3d_um;
typedef tdual_float_3d::t_host_const_um t_float_3d_const_um;
typedef tdual_float_3d::t_host_const_randomread t_float_3d_randomread;
};

#endif

template <typename D>
struct Graph {
  using Ints = Kokkos::View<int*, D>;
  Ints offsets;
  Ints at;
  int nedges;
  KOKKOS_INLINE_FUNCTION
  int start(int i) const { return offsets(i); }
  KOKKOS_INLINE_FUNCTION
  int end(int i) const { return offsets(i + 1); }
  KOKKOS_INLINE_FUNCTION
  int count(int i) const { return end(i) - start(i); }
  KOKKOS_INLINE_FUNCTION
  int& get(int i, int j) const { return at(start(i) + j); }
};

// default SPARTA Types
typedef struct ArrayTypes<DeviceType> DAT;
typedef struct ArrayTypes<SPAHostType> HAT;

// custom data types

namespace SPARTA_NS {

  struct struct_tdual_int_1d
  { DAT::tdual_int_1d k_view; };

  struct struct_tdual_float_1d
  { DAT::tdual_float_1d k_view; };

  struct struct_tdual_int_2d
  { DAT::tdual_int_2d_lr k_view; };

  struct struct_tdual_float_2d
  { DAT::tdual_float_2d_lr k_view; };

  typedef Kokkos::DualView<struct_tdual_int_1d*, DeviceType::array_layout, DeviceType> tdual_struct_tdual_int_1d_1d;
  typedef Kokkos::DualView<struct_tdual_float_1d*, DeviceType::array_layout, DeviceType> tdual_struct_tdual_float_1d_1d;
  typedef Kokkos::DualView<struct_tdual_int_2d*, DeviceType::array_layout, DeviceType> tdual_struct_tdual_int_2d_1d;
  typedef Kokkos::DualView<struct_tdual_float_2d*, DeviceType::array_layout, DeviceType> tdual_struct_tdual_float_2d_1d;
}

#ifndef SPARTA_KOKKOS_FIXED_LISTS

// the per-style device buffers behind the KK_SR_* accessors above, and the two
//   index lists that map a surf react index to a style and to a slot within it.
// blitting a model into a buffer is the same operation, for the same reason, as
//   KKCopy::copy() (kokkos_copy.h:71): on device the model is only read, through
//   KOKKOS_INLINE_FUNCTION members, so its vtable pointer is never used and the
//   Views it carries stay alive in the original that surf->sr holds.  The host
//   lifecycle calls are made on the host half of the same bytes -- a valid
//   object of the class as far as the host is concerned, its vtable pointer
//   copied from a live instance -- and pushed to the device by sr_buf_sync().
// shared here rather than repeated in each of the eight classes that carry
//   these lists, since all eight set them up the same way.

namespace SPARTA_NS {

  template<class T>
  void sr_buf_resize(DAT::tdual_char_1d &k, DAT::t_char_1d &d, int n)
  {
    const size_t need = (size_t) (n > 0 ? n : 1) * sizeof(T);
    if (k.view_device().extent(0) < need) {
      k = DAT::tdual_char_1d("surf_react:models",need);
      d = k.view_device();
    }
  }

  template<class T>
  void sr_buf_blit(DAT::tdual_char_1d &k, int slot, T *obj)
  {
    char *dst = k.view_host().data() + (size_t) slot*sizeof(T);
    memcpy((void*) dst, (const void*) obj, sizeof(T));
    ((T *) dst)->copy = 1;
  }

  inline void sr_buf_sync(DAT::tdual_char_1d &k, DAT::t_char_1d &d)
  {
    if (k.view_device().extent(0) == 0) return;
    k.modify_host();
    k.sync_device();
    d = k.view_device();
  }

  inline void sr_idx_resize(DAT::tdual_int_1d &k, DAT::t_int_1d &d, int n)
  {
    const size_t need = (size_t) (n > 0 ? n : 1);
    if (k.view_device().extent(0) < need) {
      k = DAT::tdual_int_1d("surf_react:index",need);
      d = k.view_device();
    }
  }

  inline void sr_idx_sync(DAT::tdual_int_1d &k, DAT::t_int_1d &d)
  {
    if (k.view_device().extent(0) == 0) return;
    k.modify_host();
    k.sync_device();
    d = k.view_device();
  }
}

#endif

template<class DeviceType, class BufferView, class DualView>
void buffer_view(BufferView &buf, DualView &view,
                 const size_t n0,
                 const size_t n1 = 0,
                 const size_t n2 = 0,
                 const size_t n3 = 0,
                 const size_t n4 = 0,
                 const size_t n5 = 0,
                 const size_t n6 = 0,
                 const size_t n7 = 0) {

  buf = BufferView(
          view.d_view.data(),
          n0,n1,n2,n3,n4,n5,n6,n7);

}

template<class DeviceType>
struct MemsetZeroFunctor {
  typedef DeviceType  execution_space ;
  void* ptr;
  KOKKOS_INLINE_FUNCTION void operator()(const int i) const {
    ((int*)ptr)[i] = 0;
  }
};

#define SPARTA_LAMBDA KOKKOS_LAMBDA
#define SPARTA_CLASS_LAMBDA KOKKOS_CLASS_LAMBDA

namespace SPARTA_NS {
template <typename Device>
Kokkos::View<int*, Device> offset_scan(Kokkos::View<int*, Device> a, int& total);
}

#ifdef SPARTA_KOKKOS_GPU
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) || defined(__SYCL_DEVICE_ONLY__)
#define SPARTA_KK_DEVICE_COMPILE
#endif
#endif

#endif

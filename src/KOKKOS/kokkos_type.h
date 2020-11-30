/* -*- c++ -*- ----------------------------------------------------------
   SPARTA - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://sparta.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

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
#include <impl/Kokkos_Timer.hpp>
#include <Kokkos_Vectorization.hpp>
#include <Kokkos_ScatterView.hpp>

#include "particle.h"
#include "grid.h"
#include "surf.h"

#define MAX_TYPES_STACKPARAMS 12
#define NeighClusterSize 8

  struct sparta_float3 {
    float x,y,z;
    KOKKOS_INLINE_FUNCTION
    sparta_float3():x(0.0f),z(0.0f),y(0.0f) {}

    KOKKOS_INLINE_FUNCTION
    void operator += (const sparta_float3& tmp) {
      x+=tmp.x;
      y+=tmp.y;
      z+=tmp.z;
    }
    KOKKOS_INLINE_FUNCTION
    void operator += (const sparta_float3& tmp) volatile {
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
    KOKKOS_INLINE_FUNCTION
    void operator = (const sparta_float3& tmp) volatile {
      x=tmp.x;
      y=tmp.y;
      z=tmp.z;
    }
  };

  struct sparta_double3 {
    double x,y,z;
    KOKKOS_INLINE_FUNCTION
    sparta_double3():x(0.0),z(0.0),y(0.0) {}

    KOKKOS_INLINE_FUNCTION
    void operator += (const sparta_double3& tmp) {
      x+=tmp.x;
      y+=tmp.y;
      z+=tmp.z;
    }
    KOKKOS_INLINE_FUNCTION
    void operator += (const sparta_double3& tmp) volatile {
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
    KOKKOS_INLINE_FUNCTION
    void operator = (const sparta_double3& tmp) volatile {
      x=tmp.x;
      y=tmp.y;
      z=tmp.z;
    }
  };

#if !defined(__CUDACC__) && !defined(__VECTOR_TYPES_H__)
  struct double2 {
    double x, y;
  };
  struct float2 {
    float x, y;
  };
  struct float4 {
    float x, y, z, w;
  };
  struct double4 {
    double x, y, z, w;
  };
#endif
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

#ifndef PRECISION
#define PRECISION 2
#endif
#if PRECISION==1
typedef float SPARTA_FLOAT;
typedef float2 SPARTA_FLOAT2;
typedef sparta_float3 SPARTA_FLOAT3;
typedef float4 SPARTA_FLOAT4;
#else
typedef double SPARTA_FLOAT;
typedef double2 SPARTA_FLOAT2;
typedef sparta_double3 SPARTA_FLOAT3;
typedef double4 SPARTA_FLOAT4;
#endif

#ifndef PREC_FORCE
#define PREC_FORCE PRECISION
#endif

#if PREC_FORCE==1
typedef float F_FLOAT;
typedef float2 F_FLOAT2;
typedef sparta_float3 F_FLOAT3;
typedef float4 F_FLOAT4;
#else
typedef double F_FLOAT;
typedef double2 F_FLOAT2;
typedef sparta_double3 F_FLOAT3;
typedef double4 F_FLOAT4;
#endif

#ifndef PREC_ENERGY
#define PREC_ENERGY PRECISION
#endif

#if PREC_ENERGY==1
typedef float E_FLOAT;
typedef float2 E_FLOAT2;
typedef float4 E_FLOAT4;
#else
typedef double E_FLOAT;
typedef double2 E_FLOAT2;
typedef double4 E_FLOAT4;
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

  KOKKOS_INLINE_FUNCTION
  void operator+=(const volatile s_EV_FLOAT &rhs) volatile {
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
#define PREC_POS PRECISION
#endif

#if PREC_POS==1
typedef float X_FLOAT;
typedef float2 X_FLOAT2;
typedef float4 X_FLOAT4;
#else
typedef double X_FLOAT;
typedef double2 X_FLOAT2;
typedef double4 X_FLOAT4;
#endif

#ifndef PREC_VELOCITIES
#define PREC_VELOCITIES PRECISION
#endif

#if PREC_VELOCITIES==1
typedef float V_FLOAT;
typedef float2 V_FLOAT2;
typedef float4 V_FLOAT4;
#else
typedef double V_FLOAT;
typedef double2 V_FLOAT2;
typedef double4 V_FLOAT4;
#endif

#if PREC_KSPACE==1
typedef float K_FLOAT;
typedef float2 K_FLOAT2;
typedef float4 K_FLOAT4;
#else
typedef double K_FLOAT;
typedef double2 K_FLOAT2;
typedef double4 K_FLOAT4;
#endif

// ------------------------------------------------------------------------

// SPARTA types

namespace SPARTA_NS {

  typedef Kokkos::
    DualView<Particle::OnePart*, DeviceType::array_layout, DeviceType> tdual_particle_1d;
  typedef tdual_particle_1d::t_dev t_particle_1d;
  typedef tdual_particle_1d::t_host t_host_particle_1d;

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
  DualView<int*[3], Kokkos::LayoutRight, DeviceType> tdual_int_1d_3;
typedef tdual_int_1d_3::t_dev t_int_1d_3;
typedef tdual_int_1d_3::t_dev_const t_int_1d_3_const;
typedef tdual_int_1d_3::t_dev_um t_int_1d_3_um;
typedef tdual_int_1d_3::t_dev_const_um t_int_1d_3_const_um;
typedef tdual_int_1d_3::t_dev_const_randomread t_int_1d_3_randomread;

typedef Kokkos::
  DualView<int**, Kokkos::LayoutRight, DeviceType> tdual_int_2d;
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
  DualView<SPARTA_NS::cellint**, Kokkos::LayoutRight, DeviceType>
  tdual_cellint_2d;
typedef tdual_cellint_2d::t_dev t_cellint_2d;
typedef tdual_cellint_2d::t_dev_const t_cellint_2d_const;
typedef tdual_cellint_2d::t_dev_um t_cellint_2d_um;
typedef tdual_cellint_2d::t_dev_const_um t_cellint_2d_const_um;
typedef tdual_cellint_2d::t_dev_const_randomread t_cellint_2d_randomread;

typedef Kokkos::
  DualView<SPARTA_NS::surfint*, DeviceType::array_layout, DeviceType>
  tdual_surfint_1d;
typedef tdual_surfint_1d::t_dev t_surfint_1d;
typedef tdual_surfint_1d::t_dev_const t_surfint_1d_const;
typedef tdual_surfint_1d::t_dev_um t_surfint_1d_um;
typedef tdual_surfint_1d::t_dev_const_um t_surfint_1d_const_um;
typedef tdual_surfint_1d::t_dev_const_randomread t_surfint_1d_randomread;

typedef Kokkos::
  DualView<SPARTA_NS::surfint**, Kokkos::LayoutRight, DeviceType>
  tdual_surfint_2d;
typedef tdual_surfint_2d::t_dev t_surfint_2d;
typedef tdual_surfint_2d::t_dev_const t_surfint_2d_const;
typedef tdual_surfint_2d::t_dev_um t_surfint_2d_um;
typedef tdual_surfint_2d::t_dev_const_um t_surfint_2d_const_um;
typedef tdual_surfint_2d::t_dev_const_randomread t_surfint_2d_randomread;

typedef Kokkos::
  DualView<double*, Kokkos::LayoutRight, DeviceType> tdual_double_1d;
typedef tdual_double_1d::t_dev t_double_1d;
typedef tdual_double_1d::t_dev_const t_double_1d_const;
typedef tdual_double_1d::t_dev_um t_double_1d_um;
typedef tdual_double_1d::t_dev_const_um t_double_1d_const_um;
typedef tdual_double_1d::t_dev_const_randomread t_double_1d_randomread;

typedef Kokkos::
  DualView<double**, Kokkos::LayoutRight, DeviceType> tdual_double_2d;
typedef tdual_double_2d::t_dev t_double_2d;
typedef tdual_double_2d::t_dev_const t_double_2d_const;
typedef tdual_double_2d::t_dev_um t_double_2d_um;
typedef tdual_double_2d::t_dev_const_um t_double_2d_const_um;
typedef tdual_double_2d::t_dev_const_randomread t_double_2d_randomread;

// 1d float array n

typedef Kokkos::DualView<SPARTA_FLOAT*, DeviceType::array_layout, DeviceType> tdual_float_1d;
typedef tdual_float_1d::t_dev t_float_1d;
typedef tdual_float_1d::t_dev_const t_float_1d_const;
typedef tdual_float_1d::t_dev_um t_float_1d_um;
typedef tdual_float_1d::t_dev_const_um t_float_1d_const_um;
typedef tdual_float_1d::t_dev_const_randomread t_float_1d_randomread;

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

//3d float array n
typedef Kokkos::DualView<SPARTA_FLOAT***, DeviceType::array_layout, DeviceType> tdual_float_3d;
typedef tdual_float_3d::t_dev t_float_3d;
typedef tdual_float_3d::t_dev_const t_float_3d_const;
typedef tdual_float_3d::t_dev_um t_float_3d_um;
typedef tdual_float_3d::t_dev_const_um t_float_3d_const_um;
typedef tdual_float_3d::t_dev_const_randomread t_float_3d_randomread;

//Position Types
//1d X_FLOAT array n
typedef Kokkos::DualView<X_FLOAT*, DeviceType::array_layout, DeviceType> tdual_xfloat_1d;
typedef tdual_xfloat_1d::t_dev t_xfloat_1d;
typedef tdual_xfloat_1d::t_dev_const t_xfloat_1d_const;
typedef tdual_xfloat_1d::t_dev_um t_xfloat_1d_um;
typedef tdual_xfloat_1d::t_dev_const_um t_xfloat_1d_const_um;
typedef tdual_xfloat_1d::t_dev_const_randomread t_xfloat_1d_randomread;

//2d X_FLOAT array n*m
typedef Kokkos::DualView<X_FLOAT**, Kokkos::LayoutRight, DeviceType> tdual_xfloat_2d;
typedef tdual_xfloat_2d::t_dev t_xfloat_2d;
typedef tdual_xfloat_2d::t_dev_const t_xfloat_2d_const;
typedef tdual_xfloat_2d::t_dev_um t_xfloat_2d_um;
typedef tdual_xfloat_2d::t_dev_const_um t_xfloat_2d_const_um;
typedef tdual_xfloat_2d::t_dev_const_randomread t_xfloat_2d_randomread;

//2d X_FLOAT array n*4
#ifdef SPARTA_KOKKOS_NO_LEGACY
typedef Kokkos::DualView<X_FLOAT*[3], Kokkos::LayoutLeft, DeviceType> tdual_x_array;
#else
typedef Kokkos::DualView<X_FLOAT*[3], Kokkos::LayoutRight, DeviceType> tdual_x_array;
#endif
typedef tdual_x_array::t_dev t_x_array;
typedef tdual_x_array::t_dev_const t_x_array_const;
typedef tdual_x_array::t_dev_um t_x_array_um;
typedef tdual_x_array::t_dev_const_um t_x_array_const_um;
typedef tdual_x_array::t_dev_const_randomread t_x_array_randomread;

//Velocity Types
//1d V_FLOAT array n
typedef Kokkos::DualView<V_FLOAT*, DeviceType::array_layout, DeviceType> tdual_vfloat_1d;
typedef tdual_vfloat_1d::t_dev t_vfloat_1d;
typedef tdual_vfloat_1d::t_dev_const t_vfloat_1d_const;
typedef tdual_vfloat_1d::t_dev_um t_vfloat_1d_um;
typedef tdual_vfloat_1d::t_dev_const_um t_vfloat_1d_const_um;
typedef tdual_vfloat_1d::t_dev_const_randomread t_vfloat_1d_randomread;

//2d V_FLOAT array n*m
typedef Kokkos::DualView<V_FLOAT**, Kokkos::LayoutRight, DeviceType> tdual_vfloat_2d;
typedef tdual_vfloat_2d::t_dev t_vfloat_2d;
typedef tdual_vfloat_2d::t_dev_const t_vfloat_2d_const;
typedef tdual_vfloat_2d::t_dev_um t_vfloat_2d_um;
typedef tdual_vfloat_2d::t_dev_const_um t_vfloat_2d_const_um;
typedef tdual_vfloat_2d::t_dev_const_randomread t_vfloat_2d_randomread;

//2d V_FLOAT array n*3
typedef Kokkos::DualView<V_FLOAT*[3], Kokkos::LayoutRight, DeviceType> tdual_v_array;
//typedef Kokkos::DualView<V_FLOAT*[3], DeviceType::array_layout, DeviceType> tdual_v_array;
typedef tdual_v_array::t_dev t_v_array;
typedef tdual_v_array::t_dev_const t_v_array_const;
typedef tdual_v_array::t_dev_um t_v_array_um;
typedef tdual_v_array::t_dev_const_um t_v_array_const_um;
typedef tdual_v_array::t_dev_const_randomread t_v_array_randomread;

//Force Types
//1d F_FLOAT array n

typedef Kokkos::DualView<F_FLOAT*, DeviceType::array_layout, DeviceType> tdual_ffloat_1d;
typedef tdual_ffloat_1d::t_dev t_ffloat_1d;
typedef tdual_ffloat_1d::t_dev_const t_ffloat_1d_const;
typedef tdual_ffloat_1d::t_dev_um t_ffloat_1d_um;
typedef tdual_ffloat_1d::t_dev_const_um t_ffloat_1d_const_um;
typedef tdual_ffloat_1d::t_dev_const_randomread t_ffloat_1d_randomread;

//2d F_FLOAT array n*m

typedef Kokkos::DualView<F_FLOAT**, Kokkos::LayoutRight, DeviceType> tdual_float_2d_lr;
typedef tdual_float_2d_lr::t_dev t_float_2d_lr;
typedef tdual_float_2d_lr::t_dev_const t_float_2d_lr_const;
typedef tdual_float_2d_lr::t_dev_um t_float_2d_lr_um;
typedef tdual_float_2d_lr::t_dev_const_um t_float_2d_lr_const_um;
typedef tdual_float_2d_lr::t_dev_const_randomread t_float_2d_lr_randomread;

//2d F_FLOAT array n*3

typedef Kokkos::DualView<F_FLOAT*[3], Kokkos::LayoutRight, DeviceType> tdual_f_array;
//typedef Kokkos::DualView<F_FLOAT*[3], DeviceType::array_layout, DeviceType> tdual_f_array;
typedef tdual_f_array::t_dev t_f_array;
typedef tdual_f_array::t_dev_const t_f_array_const;
typedef tdual_f_array::t_dev_um t_f_array_um;
typedef tdual_f_array::t_dev_const_um t_f_array_const_um;
typedef tdual_f_array::t_dev_const_randomread t_f_array_randomread;

//2d F_FLOAT array n*6 (for virial)

typedef Kokkos::DualView<F_FLOAT*[6], Kokkos::LayoutRight, DeviceType> tdual_virial_array;
typedef tdual_virial_array::t_dev t_virial_array;
typedef tdual_virial_array::t_dev_const t_virial_array_const;
typedef tdual_virial_array::t_dev_um t_virial_array_um;
typedef tdual_virial_array::t_dev_const_um t_virial_array_const_um;
typedef tdual_virial_array::t_dev_const_randomread t_virial_array_randomread;

//Energy Types
//1d E_FLOAT array n

typedef Kokkos::DualView<E_FLOAT*, DeviceType::array_layout, DeviceType> tdual_efloat_1d;
typedef tdual_efloat_1d::t_dev t_efloat_1d;
typedef tdual_efloat_1d::t_dev_const t_efloat_1d_const;
typedef tdual_efloat_1d::t_dev_um t_efloat_1d_um;
typedef tdual_efloat_1d::t_dev_const_um t_efloat_1d_const_um;
typedef tdual_efloat_1d::t_dev_const_randomread t_efloat_1d_randomread;

//2d E_FLOAT array n*m

typedef Kokkos::DualView<E_FLOAT**, Kokkos::LayoutRight, DeviceType> tdual_efloat_2d;
typedef tdual_efloat_2d::t_dev t_efloat_2d;
typedef tdual_efloat_2d::t_dev_const t_efloat_2d_const;
typedef tdual_efloat_2d::t_dev_um t_efloat_2d_um;
typedef tdual_efloat_2d::t_dev_const_um t_efloat_2d_const_um;
typedef tdual_efloat_2d::t_dev_const_randomread t_efloat_2d_randomread;

//2d E_FLOAT array n*3

typedef Kokkos::DualView<E_FLOAT*[3], Kokkos::LayoutRight, DeviceType> tdual_e_array;
typedef tdual_e_array::t_dev t_e_array;
typedef tdual_e_array::t_dev_const t_e_array_const;
typedef tdual_e_array::t_dev_um t_e_array_um;
typedef tdual_e_array::t_dev_const_um t_e_array_const_um;
typedef tdual_e_array::t_dev_const_randomread t_e_array_randomread;

//Neighbor Types

typedef Kokkos::DualView<int**, DeviceType::array_layout, DeviceType> tdual_neighbors_2d;
typedef tdual_neighbors_2d::t_dev t_neighbors_2d;
typedef tdual_neighbors_2d::t_dev_const t_neighbors_2d_const;
typedef tdual_neighbors_2d::t_dev_um t_neighbors_2d_um;
typedef tdual_neighbors_2d::t_dev_const_um t_neighbors_2d_const_um;
typedef tdual_neighbors_2d::t_dev_const_randomread t_neighbors_2d_randomread;

typedef Kokkos::DualView<double*,Kokkos::LayoutStride> tdual_float_1d_strided;
typedef tdual_float_1d_strided::t_dev t_float_1d_strided;
typedef tdual_float_1d_strided::t_dev_um t_float_1d_strided_um;
};

#ifdef KOKKOS_ENABLE_CUDA
template <>
struct ArrayTypes<SPAHostType> {

//Scalar Types

typedef Kokkos::DualView<int, DeviceType::array_layout, DeviceType> tdual_int_scalar;
typedef tdual_int_scalar::t_host t_int_scalar;
typedef tdual_int_scalar::t_host_const t_int_scalar_const;
typedef tdual_int_scalar::t_host_um t_int_scalar_um;
typedef tdual_int_scalar::t_host_const_um t_int_scalar_const_um;

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

typedef Kokkos::DualView<int*[3], Kokkos::LayoutRight, DeviceType> tdual_int_1d_3;
typedef tdual_int_1d_3::t_host t_int_1d_3;
typedef tdual_int_1d_3::t_host_const t_int_1d_3_const;
typedef tdual_int_1d_3::t_host_um t_int_1d_3_um;
typedef tdual_int_1d_3::t_host_const_um t_int_1d_3_const_um;
typedef tdual_int_1d_3::t_host_const_randomread t_int_1d_3_randomread;

typedef Kokkos::DualView<int**, Kokkos::LayoutRight, DeviceType> tdual_int_2d;
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

typedef Kokkos::
  DualView<SPARTA_NS::cellint**, Kokkos::LayoutRight, DeviceType>
  tdual_cellint_2d;
typedef tdual_cellint_2d::t_host t_cellint_2d;
typedef tdual_cellint_2d::t_host_const t_cellint_2d_const;
typedef tdual_cellint_2d::t_host_um t_cellint_2d_um;
typedef tdual_cellint_2d::t_host_const_um t_cellint_2d_const_um;
typedef tdual_cellint_2d::t_host_const_randomread t_cellint_2d_randomread;

typedef Kokkos::DualView<SPARTA_NS::surfint*, DeviceType::array_layout, DeviceType> tdual_surfint_1d;
typedef tdual_surfint_1d::t_host t_surfint_1d;
typedef tdual_surfint_1d::t_host_const t_surfint_1d_const;
typedef tdual_surfint_1d::t_host_um t_surfint_1d_um;
typedef tdual_surfint_1d::t_host_const_um t_surfint_1d_const_um;
typedef tdual_surfint_1d::t_host_const_randomread t_surfint_1d_randomread;

typedef Kokkos::
  DualView<SPARTA_NS::surfint**, Kokkos::LayoutRight, DeviceType>
  tdual_surfint_2d;
typedef tdual_surfint_2d::t_host t_surfint_2d;
typedef tdual_surfint_2d::t_host_const t_surfint_2d_const;
typedef tdual_surfint_2d::t_host_um t_surfint_2d_um;
typedef tdual_surfint_2d::t_host_const_um t_surfint_2d_const_um;
typedef tdual_surfint_2d::t_host_const_randomread t_surfint_2d_randomread;

typedef Kokkos::
  DualView<double*, Kokkos::LayoutRight, DeviceType> tdual_double_1d;
typedef tdual_double_1d::t_host t_double_1d;
typedef tdual_double_1d::t_host_const t_double_1d_const;
typedef tdual_double_1d::t_host_um t_double_1d_um;
typedef tdual_double_1d::t_host_const_um t_double_1d_const_um;
typedef tdual_double_1d::t_host_const_randomread t_double_1d_randomread;

typedef Kokkos::
  DualView<double**, Kokkos::LayoutRight, DeviceType> tdual_double_2d;
typedef tdual_double_2d::t_host t_double_2d;
typedef tdual_double_2d::t_host_const t_double_2d_const;
typedef tdual_double_2d::t_host_um t_double_2d_um;
typedef tdual_double_2d::t_host_const_um t_double_2d_const_um;
typedef tdual_double_2d::t_host_const_randomread t_double_2d_randomread;

//1d float array n
typedef Kokkos::DualView<SPARTA_FLOAT*, DeviceType::array_layout, DeviceType> tdual_float_1d;
typedef tdual_float_1d::t_host t_float_1d;
typedef tdual_float_1d::t_host_const t_float_1d_const;
typedef tdual_float_1d::t_host_um t_float_1d_um;
typedef tdual_float_1d::t_host_const_um t_float_1d_const_um;
typedef tdual_float_1d::t_host_const_randomread t_float_1d_randomread;

//1d float array n[3]
typedef Kokkos::DualView<SPARTA_FLOAT*[3], DeviceType::array_layout, DeviceType> tdual_float_1d_3;
typedef tdual_float_1d_3::t_host t_float_1d_3;
typedef tdual_float_1d_3::t_host_const t_float_1d_3_const;
typedef tdual_float_1d_3::t_host_um t_float_1d_3_um;
typedef tdual_float_1d_3::t_host_const_um t_float_1d_3_const_um;
typedef tdual_float_1d_3::t_host_const_randomread t_float_1d_3_randomread;

//2d float array n
typedef Kokkos::DualView<SPARTA_FLOAT**, DeviceType::array_layout, DeviceType> tdual_float_2d;
typedef tdual_float_2d::t_host t_float_2d;
typedef tdual_float_2d::t_host_const t_float_2d_const;
typedef tdual_float_2d::t_host_um t_float_2d_um;
typedef tdual_float_2d::t_host_const_um t_float_2d_const_um;
typedef tdual_float_2d::t_host_const_randomread t_float_2d_randomread;

//3d float array n
typedef Kokkos::DualView<SPARTA_FLOAT***, DeviceType::array_layout, DeviceType> tdual_float_3d;
typedef tdual_float_3d::t_host t_float_3d;
typedef tdual_float_3d::t_host_const t_float_3d_const;
typedef tdual_float_3d::t_host_um t_float_3d_um;
typedef tdual_float_3d::t_host_const_um t_float_3d_const_um;
typedef tdual_float_3d::t_host_const_randomread t_float_3d_randomread;


//Position Types
//1d X_FLOAT array n
typedef Kokkos::DualView<X_FLOAT*, DeviceType::array_layout, DeviceType> tdual_xfloat_1d;
typedef tdual_xfloat_1d::t_host t_xfloat_1d;
typedef tdual_xfloat_1d::t_host_const t_xfloat_1d_const;
typedef tdual_xfloat_1d::t_host_um t_xfloat_1d_um;
typedef tdual_xfloat_1d::t_host_const_um t_xfloat_1d_const_um;
typedef tdual_xfloat_1d::t_host_const_randomread t_xfloat_1d_randomread;

//2d X_FLOAT array n*m
typedef Kokkos::DualView<X_FLOAT**, Kokkos::LayoutRight, DeviceType> tdual_xfloat_2d;
typedef tdual_xfloat_2d::t_host t_xfloat_2d;
typedef tdual_xfloat_2d::t_host_const t_xfloat_2d_const;
typedef tdual_xfloat_2d::t_host_um t_xfloat_2d_um;
typedef tdual_xfloat_2d::t_host_const_um t_xfloat_2d_const_um;
typedef tdual_xfloat_2d::t_host_const_randomread t_xfloat_2d_randomread;

//2d X_FLOAT array n*3
typedef Kokkos::DualView<X_FLOAT*[3], Kokkos::LayoutRight, DeviceType> tdual_x_array;
typedef tdual_x_array::t_host t_x_array;
typedef tdual_x_array::t_host_const t_x_array_const;
typedef tdual_x_array::t_host_um t_x_array_um;
typedef tdual_x_array::t_host_const_um t_x_array_const_um;
typedef tdual_x_array::t_host_const_randomread t_x_array_randomread;

//Velocity Types
//1d V_FLOAT array n
typedef Kokkos::DualView<V_FLOAT*, DeviceType::array_layout, DeviceType> tdual_vfloat_1d;
typedef tdual_vfloat_1d::t_host t_vfloat_1d;
typedef tdual_vfloat_1d::t_host_const t_vfloat_1d_const;
typedef tdual_vfloat_1d::t_host_um t_vfloat_1d_um;
typedef tdual_vfloat_1d::t_host_const_um t_vfloat_1d_const_um;
typedef tdual_vfloat_1d::t_host_const_randomread t_vfloat_1d_randomread;

//2d V_FLOAT array n*m
typedef Kokkos::DualView<V_FLOAT**, Kokkos::LayoutRight, DeviceType> tdual_vfloat_2d;
typedef tdual_vfloat_2d::t_host t_vfloat_2d;
typedef tdual_vfloat_2d::t_host_const t_vfloat_2d_const;
typedef tdual_vfloat_2d::t_host_um t_vfloat_2d_um;
typedef tdual_vfloat_2d::t_host_const_um t_vfloat_2d_const_um;
typedef tdual_vfloat_2d::t_host_const_randomread t_vfloat_2d_randomread;

//2d V_FLOAT array n*3
typedef Kokkos::DualView<V_FLOAT*[3], Kokkos::LayoutRight, DeviceType> tdual_v_array;
//typedef Kokkos::DualView<V_FLOAT*[3], DeviceType::array_layout, DeviceType> tdual_v_array;
typedef tdual_v_array::t_host t_v_array;
typedef tdual_v_array::t_host_const t_v_array_const;
typedef tdual_v_array::t_host_um t_v_array_um;
typedef tdual_v_array::t_host_const_um t_v_array_const_um;
typedef tdual_v_array::t_host_const_randomread t_v_array_randomread;

//Force Types
//1d F_FLOAT array n
typedef Kokkos::DualView<F_FLOAT*, DeviceType::array_layout, DeviceType> tdual_ffloat_1d;
typedef tdual_ffloat_1d::t_host t_ffloat_1d;
typedef tdual_ffloat_1d::t_host_const t_ffloat_1d_const;
typedef tdual_ffloat_1d::t_host_um t_ffloat_1d_um;
typedef tdual_ffloat_1d::t_host_const_um t_ffloat_1d_const_um;
typedef tdual_ffloat_1d::t_host_const_randomread t_ffloat_1d_randomread;

//2d F_FLOAT array n*m
typedef Kokkos::DualView<F_FLOAT**, Kokkos::LayoutRight, DeviceType> tdual_float_2d_lr;
typedef tdual_float_2d_lr::t_host t_float_2d_lr;
typedef tdual_float_2d_lr::t_host_const t_float_2d_lr_const;
typedef tdual_float_2d_lr::t_host_um t_float_2d_lr_um;
typedef tdual_float_2d_lr::t_host_const_um t_float_2d_lr_const_um;
typedef tdual_float_2d_lr::t_host_const_randomread t_float_2d_lr_randomread;

//2d F_FLOAT array n*3
typedef Kokkos::DualView<F_FLOAT*[3], Kokkos::LayoutRight, DeviceType> tdual_f_array;
//typedef Kokkos::DualView<F_FLOAT*[3], DeviceType::array_layout, DeviceType> tdual_f_array;
typedef tdual_f_array::t_host t_f_array;
typedef tdual_f_array::t_host_const t_f_array_const;
typedef tdual_f_array::t_host_um t_f_array_um;
typedef tdual_f_array::t_host_const_um t_f_array_const_um;
typedef tdual_f_array::t_host_const_randomread t_f_array_randomread;

//2d F_FLOAT array n*6 (for virial)
typedef Kokkos::DualView<F_FLOAT*[6], Kokkos::LayoutRight, DeviceType> tdual_virial_array;
typedef tdual_virial_array::t_host t_virial_array;
typedef tdual_virial_array::t_host_const t_virial_array_const;
typedef tdual_virial_array::t_host_um t_virial_array_um;
typedef tdual_virial_array::t_host_const_um t_virial_array_const_um;
typedef tdual_virial_array::t_host_const_randomread t_virial_array_randomread;



//Energy Types
//1d E_FLOAT array n
typedef Kokkos::DualView<E_FLOAT*, DeviceType::array_layout, DeviceType> tdual_efloat_1d;
typedef tdual_efloat_1d::t_host t_efloat_1d;
typedef tdual_efloat_1d::t_host_const t_efloat_1d_const;
typedef tdual_efloat_1d::t_host_um t_efloat_1d_um;
typedef tdual_efloat_1d::t_host_const_um t_efloat_1d_const_um;
typedef tdual_efloat_1d::t_host_const_randomread t_efloat_1d_randomread;

//2d E_FLOAT array n*m
typedef Kokkos::DualView<E_FLOAT**, Kokkos::LayoutRight, DeviceType> tdual_efloat_2d;
typedef tdual_efloat_2d::t_host t_efloat_2d;
typedef tdual_efloat_2d::t_host_const t_efloat_2d_const;
typedef tdual_efloat_2d::t_host_um t_efloat_2d_um;
typedef tdual_efloat_2d::t_host_const_um t_efloat_2d_const_um;
typedef tdual_efloat_2d::t_host_const_randomread t_efloat_2d_randomread;

//2d E_FLOAT array n*3
typedef Kokkos::DualView<E_FLOAT*[3], Kokkos::LayoutRight, DeviceType> tdual_e_array;
typedef tdual_e_array::t_host t_e_array;
typedef tdual_e_array::t_host_const t_e_array_const;
typedef tdual_e_array::t_host_um t_e_array_um;
typedef tdual_e_array::t_host_const_um t_e_array_const_um;
typedef tdual_e_array::t_host_const_randomread t_e_array_randomread;

//Neighbor Types
typedef Kokkos::DualView<int**, DeviceType::array_layout, DeviceType> tdual_neighbors_2d;
typedef tdual_neighbors_2d::t_host t_neighbors_2d;
typedef tdual_neighbors_2d::t_host_const t_neighbors_2d_const;
typedef tdual_neighbors_2d::t_host_um t_neighbors_2d_um;
typedef tdual_neighbors_2d::t_host_const_um t_neighbors_2d_const_um;
typedef tdual_neighbors_2d::t_host_const_randomread t_neighbors_2d_randomread;

typedef Kokkos::DualView<double*,Kokkos::LayoutStride> tdual_float_1d_strided;
typedef tdual_float_1d_strided::t_host t_float_1d_strided;
typedef tdual_float_1d_strided::t_host_um t_float_1d_strided_um;
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

//default SPARTA Types
typedef struct ArrayTypes<DeviceType> DAT;
typedef struct ArrayTypes<SPAHostType> HAT;

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
          view.template d_view.data(),
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

template<class ViewType>
void memset_kokkos (ViewType &view) {
  static MemsetZeroFunctor<typename ViewType::execution_space> f;
  f.ptr = view.data();
  Kokkos::parallel_for(view.scan()*sizeof(typename ViewType::value_type)/4, f);
  ViewType::execution_space::fence();
}

#ifdef KOKKOS_ENABLE_CUDA
#define SPARTA_LAMBDA [=] __device__
#else
#define SPARTA_LAMBDA [=]
#endif

namespace SPARTA_NS {
template <typename Device>
Kokkos::View<int*, Device> offset_scan(Kokkos::View<int*, Device> a, int& total);
}

#endif

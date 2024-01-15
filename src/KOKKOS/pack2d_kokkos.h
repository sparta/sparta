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

// loop counters for doing a pack/unpack

#include "pack2d.h"

/* ----------------------------------------------------------------------
   Pack and unpack functions:

   pack routines copy strided values from data into contiguous locs in buf
   unpack routines copy contiguous values from buf into strided locs in data
   different versions of unpack depending on permutation
     and # of values/element
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   pack from data -> buf
------------------------------------------------------------------------- */

#include "fftdata_kokkos.h"

namespace SPARTA_NS {

template<class DeviceType>
class PackKokkos2d {
 public:
  typedef DeviceType device_type;
  typedef FFTArrayTypes<DeviceType> FFT_AT;

struct pack_2d_functor {
public:
  typedef DeviceType device_type;
  typedef FFTArrayTypes<DeviceType> FFT_AT;
  typename FFT_AT::t_FFT_SCALAR_1d_um d_buf,d_data;
  int buf_offset,data_offset;
  int nfast;                 // # of elements in fast index
  int nslow;                 // # of elements in slow index
  int nstride;               // stride between successive slow indices

  pack_2d_functor(typename FFT_AT::t_FFT_SCALAR_1d_um d_buf_, int buf_offset_, typename FFT_AT::t_FFT_SCALAR_1d_um d_data_, int data_offset_, struct pack_plan_2d *plan):
    d_buf(d_buf_),
    d_data(d_data_)
    {
      buf_offset    = buf_offset_        ;
      data_offset   = data_offset_       ;
      nfast         = plan->nfast        ;
      nslow         = plan->nslow        ;
      nstride       = plan->nstride      ;
    }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &index) const {
 
    const int slow = index / nfast;
    const int fast = index % nfast;
    const int in = slow*nfast + fast;
    const int out = slow*nstride + fast;
    d_buf[buf_offset + in] = d_data[data_offset + out];
  }
};

static void pack_2d(typename FFT_AT::t_FFT_SCALAR_1d_um d_data, int data_offset, typename FFT_AT::t_FFT_SCALAR_1d_um d_buf, int buf_offset, struct pack_plan_2d *plan)
{
  const int nslow = plan->nslow;
  const int nfast = plan->nfast;
  pack_2d_functor f(d_buf,buf_offset,d_data,data_offset,plan);
  Kokkos::parallel_for(nslow*nfast,f);
  DeviceType().fence();
}

/* ----------------------------------------------------------------------
   unpack from buf -> data
------------------------------------------------------------------------- */

struct unpack_2d_functor {
public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename FFT_AT::t_FFT_SCALAR_1d_um d_buf,d_data;
  int buf_offset,data_offset;
  int nfast;                 // # of elements in fast index
  int nslow;                 // # of elements in slow index
  int nstride;         // stride between successive slow indices

  unpack_2d_functor(typename FFT_AT::t_FFT_SCALAR_1d_um d_buf_, int buf_offset_, typename FFT_AT::t_FFT_SCALAR_1d_um d_data_, int data_offset_, struct pack_plan_2d *plan):
    d_buf(d_buf_),
    d_data(d_data_)
    {
      buf_offset    = buf_offset_        ;
      data_offset   = data_offset_       ;
      nfast         = plan->nfast        ;
      nslow         = plan->nslow        ;
      nstride       = plan->nstride      ;
    }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &index) const {

    const int slow = index / nfast;
    const int fast = index % nfast;
    const int out = slow*nfast + fast;
    const int in = slow*nstride + fast;
    d_buf[buf_offset + in] = d_data[data_offset + out];
  }
};

static void unpack_2d(typename FFT_AT::t_FFT_SCALAR_1d_um d_buf, int buf_offset, typename FFT_AT::t_FFT_SCALAR_1d_um d_data, int data_offset, struct pack_plan_2d *plan)
{
  const int nslow = plan->nslow;
  const int nfast = plan->nfast;
  unpack_2d_functor f(d_buf,buf_offset,d_data,data_offset,plan);
  Kokkos::parallel_for(nslow*nfast,f);
  DeviceType().fence();
}

/* ----------------------------------------------------------------------
   unpack from buf -> data, axis permutation, 1 value/element
------------------------------------------------------------------------- */


struct unpack_2d_permute_1_functor {
public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename FFT_AT::t_FFT_SCALAR_1d_um d_buf,d_data;
  int buf_offset,data_offset;
  int nfast;                 // # of elements in fast index
  int nslow;                 // # of elements in slow index
  int nstride;         // stride between successive slow indices

  unpack_2d_permute_1_functor(typename FFT_AT::t_FFT_SCALAR_1d_um d_buf_, int buf_offset_, typename FFT_AT::t_FFT_SCALAR_1d_um d_data_, int data_offset_, struct pack_plan_2d *plan):
    d_buf(d_buf_),
    d_data(d_data_)
    {
      buf_offset    = buf_offset_        ;
      data_offset   = data_offset_       ;
      nfast         = plan->nfast        ;
      nslow         = plan->nslow        ;
      nstride       = plan->nstride      ;
    }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &index) const {
    const int slow = index / nfast;
    const int fast = index % nfast;
    const int in = slow + fast*nstride;
    const int out = slow*nfast + fast;
    d_data[data_offset + in] = d_buf[buf_offset + out];
  }
};

static void unpack_2d_permute_1(typename FFT_AT::t_FFT_SCALAR_1d_um d_buf, int buf_offset, typename FFT_AT::t_FFT_SCALAR_1d_um d_data, int data_offset, struct pack_plan_2d *plan)
{
  const int nslow = plan->nslow;
  const int nfast = plan->nfast;
  unpack_2d_permute_1_functor f(d_buf,buf_offset,d_data,data_offset,plan);
  Kokkos::parallel_for(nslow*nfast,f);
  DeviceType().fence();
}
/* ----------------------------------------------------------------------
   unpack from buf -> data, axis permutation, 2 values/element
------------------------------------------------------------------------- */

struct unpack_2d_permute_2_functor {
public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename FFT_AT::t_FFT_SCALAR_1d_um d_buf,d_data;
  int buf_offset,data_offset;
  int nfast;                 // # of elements in fast index
  int nslow;                 // # of elements in slow index
  int nstride;               // stride between successive slow indices

  unpack_2d_permute_2_functor(typename FFT_AT::t_FFT_SCALAR_1d_um d_buf_, int buf_offset_, typename FFT_AT::t_FFT_SCALAR_1d_um d_data_, int data_offset_, struct pack_plan_2d *plan):
    d_buf(d_buf_),
    d_data(d_data_)
    {
      buf_offset    = buf_offset_        ;
      data_offset   = data_offset_       ;
      nfast         = plan->nfast        ;
      nslow         = plan->nslow        ;
      nstride       = plan->nstride      ;
    }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &index) const {
    const int slow = index / nfast;
    const int fast = index % nfast;
    const int in = 2*slow + fast*nstride;
    const int out = slow*nfast*2 + fast*2;
    d_data[data_offset + in] = d_buf[buf_offset + out];
    d_data[data_offset + in+1] = d_buf[buf_offset + out+1];
  }
};

static void unpack_2d_permute_2(typename FFT_AT::t_FFT_SCALAR_1d_um d_buf, int buf_offset, typename FFT_AT::t_FFT_SCALAR_1d_um d_data, int data_offset, struct pack_plan_2d *plan)
{
  const int nslow = plan->nslow;
  const int nfast = plan->nfast;
  unpack_2d_permute_2_functor f(d_buf,buf_offset,d_data,data_offset,plan);
  Kokkos::parallel_for(nslow*nfast,f);
  DeviceType().fence();
}

/* ----------------------------------------------------------------------
   unpack from buf -> data, axis permutation, nqty values/element
------------------------------------------------------------------------- */

struct unpack_2d_permute_n_functor {
public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename FFT_AT::t_FFT_SCALAR_1d_um d_buf,d_data;
  int buf_offset,data_offset;
  int nfast;                 // # of elements in fast index
  int nslow;                 // # of elements in slow index
  int nstride;               // stride between successive slow indices
  int nqty;                  // # of values/element

  unpack_2d_permute_n_functor(typename FFT_AT::t_FFT_SCALAR_1d_um d_buf_, int buf_offset_, typename FFT_AT::t_FFT_SCALAR_1d_um d_data_, int data_offset_, struct pack_plan_2d *plan):
    d_buf(d_buf_),
    d_data(d_data_)
    {
      buf_offset    = buf_offset_        ;
      data_offset   = data_offset_       ;
      nfast         = plan->nfast        ;
      nslow         = plan->nslow        ;
      nstride       = plan->nstride      ;
      nqty          = plan->nqty         ;
    }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &index) const {
    const int slow = index / nfast;
    const int fast = index % nfast;
    int in =  slow*nqty + fast*nstride;
    int out = slow*nfast*nqty + fast*nqty;
    for (int iqty = 0; iqty < nqty; iqty++)
      d_data[data_offset + in++] = d_buf[buf_offset + out++];
  }
};

static void unpack_2d_permute_n(typename FFT_AT::t_FFT_SCALAR_1d_um d_buf, int buf_offset, typename FFT_AT::t_FFT_SCALAR_1d_um d_data, int data_offset, struct pack_plan_2d *plan)
{
  const int nslow = plan->nslow;
  const int nfast = plan->nfast;
  unpack_2d_permute_n_functor f(d_buf,buf_offset,d_data,data_offset,plan);
  Kokkos::parallel_for(nslow*nfast,f);
  DeviceType().fence();
}


};

}

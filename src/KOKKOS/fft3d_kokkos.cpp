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

/* ----------------------------------------------------------------------
   Contributing authors: James Almgren-Bell (UT Austin), Sam Mish (UC Davis)
------------------------------------------------------------------------- */

#include "fft3d_kokkos.h"

#include "error.h"
#include "kokkos.h"
#include "remap3d_kokkos.h"

#include <cmath>

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FFT3dKokkos<DeviceType>::FFT3dKokkos(SPARTA *sparta, MPI_Comm comm, int nfast, int nmid, int nslow,
             int in_ilo, int in_ihi, int in_jlo, int in_jhi,
             int in_klo, int in_khi,
             int out_ilo, int out_ihi, int out_jlo, int out_jhi,
             int out_klo, int out_khi,
             int scaled, int permute, int *nbuf, int usecollective,
             int usegpu_aware) :
  Pointers(sparta)
{
  int nthreads = sparta->kokkos->nthreads;
  int ngpus = sparta->kokkos->ngpus;
  ExecutionSpace execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

#if defined(FFT_KOKKOS_MKL)
  if (ngpus > 0 && execution_space == Device)
    sparta->error->all(FLERR,"Cannot use the MKL library with Kokkos on GPUs");
#elif defined(FFT_KOKKOS_FFTW3)
  if (ngpus > 0 && execution_space == Device)
    sparta->error->all(FLERR,"Cannot use the FFTW library with Kokkos on GPUs");
#elif defined(FFT_KOKKOS_CUFFT)
  if (ngpus > 0 && execution_space == Host)
    sparta->error->all(FLERR,"Cannot use the cuFFT library with Kokkos on the host CPUs");
#elif defined(FFT_KOKKOS_HIPFFT)
  if (ngpus > 0 && execution_space == Host)
    sparta->error->all(FLERR,"Cannot use the hipFFT library with Kokkos on the host CPUs");

#elif defined(FFT_KOKKOS_KISS)
  // The compiler can't statically determine the stack size needed for
  //  recursive function calls in KISS FFT and the default per-thread
  //  stack size on GPUs needs to be increased to prevent stack overflows
  //  for reasonably sized FFTs
  #if defined (KOKKOS_ENABLE_CUDA)
    size_t stack_size;
    cudaDeviceGetLimit(&stack_size,cudaLimitStackSize);
    if (stack_size < 2048)
      cudaDeviceSetLimit(cudaLimitStackSize,2048);
  #endif
#endif

  plan = fft_3d_create_plan_kokkos(comm,nfast,nmid,nslow,
                            in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
                            out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,
                            scaled,permute,nbuf,usecollective,nthreads,usegpu_aware);
  if (plan == nullptr) error->one(FLERR,"Could not create 3d FFT plan");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FFT3dKokkos<DeviceType>::~FFT3dKokkos()
{
  fft_3d_destroy_plan_kokkos(plan);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FFT3dKokkos<DeviceType>::compute(typename FFT_AT::t_FFT_SCALAR_1d d_in, typename FFT_AT::t_FFT_SCALAR_1d d_out, int flag)
{
  typename FFT_AT::t_FFT_DATA_1d d_in_data((FFT_KOKKOS_DATA_POINTER)d_in.data(),d_in.size()/2);
  typename FFT_AT::t_FFT_DATA_1d d_out_data((FFT_KOKKOS_DATA_POINTER)d_out.data(),d_out.size()/2);

  fft_3d_kokkos(d_in_data,d_out_data,flag,plan);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FFT3dKokkos<DeviceType>::timing1d(typename FFT_AT::t_FFT_SCALAR_1d d_in, int nsize, int flag)
{
  typename FFT_AT::t_FFT_DATA_1d d_in_data((FFT_KOKKOS_DATA_POINTER)d_in.data(),d_in.size()/2);

  fft_3d_1d_only_kokkos(d_in_data,nsize,flag,plan);
}

/* ----------------------------------------------------------------------
   Data layout for 3d FFTs:

   data set of Nfast x Nmid x Nslow elements is owned by P procs
   on input, each proc owns a subsection of the elements
   on output, each proc will own a (possibly different) subsection
   my subsection must not overlap with any other proc's subsection,
     i.e. the union of all proc's input (or output) subsections must
     exactly tile the global Nfast x Nmid x Nslow data set
   when called from C, all subsection indices are
     C-style from 0 to N-1 where N = Nfast or Nmid or Nslow
   when called from F77, all subsection indices are
     F77-style from 1 to N where N = Nfast or Nmid or Nslow
   a proc can own 0 elements on input or output
     by specifying hi index < lo index
   on both input and output, data is stored contiguously on a processor
     with a fast-varying, mid-varying, and slow-varying index
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Perform 3d FFT

   Arguments:
   in           starting address of input data on this proc
   out          starting address of where output data for this proc
                  will be placed (can be same as in)
   flag         1 for forward FFT, -1 for backward FFT
   plan         plan returned by previous call to fft_3d_create_plan
------------------------------------------------------------------------- */

template<class DeviceType>
struct norm_functor {
public:
  typedef DeviceType device_type;
  typedef FFTArrayTypes<DeviceType> FFT_AT;
  typename FFT_AT::t_FFT_DATA_1d_um d_out;
  int norm;

  norm_functor(typename FFT_AT::t_FFT_DATA_1d &d_out_, int norm_):
    d_out(d_out_),norm(norm_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &i) const {
#if defined(FFT_KOKKOS_FFTW3) || defined(FFT_KOKKOS_CUFFT) || defined(FFT_KOKKOS_HIPFFT)
    FFT_SCALAR* out_ptr = (FFT_SCALAR *)(d_out.data()+i);
    *(out_ptr++) *= norm;
    *(out_ptr++) *= norm;
#elif defined(FFT_KOKKOS_MKL)
    d_out(i) *= norm;
#else // FFT_KOKKOS_KISS
    d_out(i).re *= norm;
    d_out(i).im *= norm;
#endif
  }
};

#ifdef FFT_KOKKOS_KISS
template<class DeviceType>
struct kiss_fft_functor {
public:
  typedef DeviceType device_type;
  typedef FFTArrayTypes<DeviceType> FFT_AT;
  typename FFT_AT::t_FFT_DATA_1d_um d_data,d_tmp;
  kiss_fft_state_kokkos<DeviceType> st;
  int length;

  kiss_fft_functor() = default;

  kiss_fft_functor(typename FFT_AT::t_FFT_DATA_1d &d_data_,typename FFT_AT::t_FFT_DATA_1d &d_tmp_, kiss_fft_state_kokkos<DeviceType> &st_, int length_):
    d_data(d_data_),
    d_tmp(d_tmp_),
    st(st_)
    {
      length = length_;
    }

  KOKKOS_INLINE_FUNCTION
  void operator() (const int &i) const {
    const int offset = i*length;
    KissFFTKokkos<DeviceType>::kiss_fft_kokkos(st,d_data,d_tmp,offset);
  }
};
#endif

template<class DeviceType>
void FFT3dKokkos<DeviceType>::fft_3d_kokkos(typename FFT_AT::t_FFT_DATA_1d d_in, typename FFT_AT::t_FFT_DATA_1d d_out, int flag, struct fft_plan_3d_kokkos<DeviceType> *plan)
{
  typename FFT_AT::t_FFT_DATA_1d d_data,d_copy;
  typename FFT_AT::t_FFT_SCALAR_1d d_in_scalar,d_data_scalar,d_out_scalar,d_copy_scalar,d_scratch_scalar;

  // pre-remap to prepare for 1st FFTs if needed
  // copy = loc for remap result

  if (plan->pre_plan) {
    if (plan->pre_target == 0) d_copy = d_out;
    else d_copy = plan->d_copy;

     d_in_scalar = typename FFT_AT::t_FFT_SCALAR_1d((FFT_SCALAR*)d_in.data(),d_in.size()*2);
     d_copy_scalar = typename FFT_AT::t_FFT_SCALAR_1d((FFT_SCALAR*)d_copy.data(),d_copy.size()*2);
     d_scratch_scalar = typename FFT_AT::t_FFT_SCALAR_1d((FFT_SCALAR*)plan->d_scratch.data(),plan->d_scratch.size()*2);

    remapKK->remap_3d_kokkos(d_in_scalar, d_copy_scalar,
             d_scratch_scalar, plan->pre_plan);

    d_data = d_copy;
  } else d_data = d_in;

  // 1d FFTs along fast axis

  #if defined(FFT_KOKKOS_MKL)
    if (flag == 1)
      DftiComputeForward(plan->handle_fast,d_data.data());
    else
      DftiComputeBackward(plan->handle_fast,d_data.data());
  #elif defined(FFT_KOKKOS_FFTW3)
    if (flag == 1)
      FFTW_API(execute_dft)(plan->plan_fast_forward,(FFT_KOKKOS_DATA*)d_data.data(),(FFT_KOKKOS_DATA*)d_data.data());
    else
      FFTW_API(execute_dft)(plan->plan_fast_backward,(FFT_KOKKOS_DATA*)d_data.data(),(FFT_KOKKOS_DATA*)d_data.data());
  #elif defined(FFT_KOKKOS_CUFFT)
    cufftExec(plan->plan_fast,d_data.data(),d_data.data(),-flag);
  #elif defined(FFT_KOKKOS_HIPFFT)
    hipfftExec(plan->plan_fast,d_data.data(),d_data.data(),-flag);
  #else
    int total = plan->total1;
    int length = plan->length1;

    typename FFT_AT::t_FFT_DATA_1d d_tmp =
     typename FFT_AT::t_FFT_DATA_1d(Kokkos::view_alloc("fft_3d:tmp",Kokkos::WithoutInitializing),d_data.extent(0));
    kiss_fft_functor<DeviceType> f;
    if (flag == 1)
      f = kiss_fft_functor<DeviceType>(d_data,d_tmp,plan->cfg_fast_forward,length);
    else
      f = kiss_fft_functor<DeviceType>(d_data,d_tmp,plan->cfg_fast_backward,length);
    Kokkos::parallel_for(total/length,f);
    d_data = d_tmp;
  #endif

  // 1st mid-remap to prepare for 2nd FFTs
  // copy = loc for remap result

  if (plan->mid1_target == 0) d_copy = d_out;
  else d_copy = plan->d_copy;

  d_data_scalar = typename FFT_AT::t_FFT_SCALAR_1d((FFT_SCALAR*)d_data.data(),d_data.size()*2);
  d_copy_scalar = typename FFT_AT::t_FFT_SCALAR_1d((FFT_SCALAR*)d_copy.data(),d_copy.size()*2);
  d_scratch_scalar = typename FFT_AT::t_FFT_SCALAR_1d((FFT_SCALAR*)plan->d_scratch.data(),plan->d_scratch.size()*2);

  remapKK->remap_3d_kokkos(d_data_scalar, d_copy_scalar,
           d_scratch_scalar, plan->mid1_plan);

  d_data = d_copy;

  // 1d FFTs along mid axis

  #if defined(FFT_KOKKOS_MKL)
    if (flag == 1)
      DftiComputeForward(plan->handle_mid,d_data.data());
    else
      DftiComputeBackward(plan->handle_mid,d_data.data());
  #elif defined(FFT_KOKKOS_FFTW3)
    if (flag == 1)
      FFTW_API(execute_dft)(plan->plan_mid_forward,(FFT_KOKKOS_DATA*)d_data.data(),(FFT_KOKKOS_DATA*)d_data.data());
    else
      FFTW_API(execute_dft)(plan->plan_mid_backward,(FFT_KOKKOS_DATA*)d_data.data(),(FFT_KOKKOS_DATA*)d_data.data());
  #elif defined(FFT_KOKKOS_CUFFT)
    cufftExec(plan->plan_mid,d_data.data(),d_data.data(),-flag);
  #elif defined(FFT_KOKKOS_HIPFFT)
    hipfftExec(plan->plan_mid,d_data.data(),d_data.data(),-flag);
  #else
    total = plan->total2;
    length = plan->length2;

    d_tmp = typename FFT_AT::t_FFT_DATA_1d(Kokkos::view_alloc("fft_3d:tmp",Kokkos::WithoutInitializing),d_data.extent(0));
    if (flag == 1)
      f = kiss_fft_functor<DeviceType>(d_data,d_tmp,plan->cfg_mid_forward,length);
    else
      f = kiss_fft_functor<DeviceType>(d_data,d_tmp,plan->cfg_mid_backward,length);
    Kokkos::parallel_for(total/length,f);
    d_data = d_tmp;
  #endif

  // 2nd mid-remap to prepare for 3rd FFTs
  // copy = loc for remap result

  if (plan->mid2_target == 0) d_copy = d_out;
  else d_copy = plan->d_copy;

  d_data_scalar = typename FFT_AT::t_FFT_SCALAR_1d((FFT_SCALAR*)d_data.data(),d_data.size()*2);
  d_copy_scalar = typename FFT_AT::t_FFT_SCALAR_1d((FFT_SCALAR*)d_copy.data(),d_copy.size()*2);
  d_scratch_scalar = typename FFT_AT::t_FFT_SCALAR_1d((FFT_SCALAR*)plan->d_scratch.data(),plan->d_scratch.size()*2);

  remapKK->remap_3d_kokkos(d_data_scalar, d_copy_scalar,
           d_scratch_scalar, plan->mid2_plan);

  d_data = d_copy;

  // 1d FFTs along slow axis

  #if defined(FFT_KOKKOS_MKL)
    if (flag == 1)
      DftiComputeForward(plan->handle_slow,d_data.data());
    else
      DftiComputeBackward(plan->handle_slow,d_data.data());
  #elif defined(FFT_KOKKOS_FFTW3)
    if (flag == 1)
      FFTW_API(execute_dft)(plan->plan_slow_forward,(FFT_KOKKOS_DATA*)d_data.data(),(FFT_KOKKOS_DATA*)d_data.data());
    else
      FFTW_API(execute_dft)(plan->plan_slow_backward,(FFT_KOKKOS_DATA*)d_data.data(),(FFT_KOKKOS_DATA*)d_data.data());
  #elif defined(FFT_KOKKOS_CUFFT)
    cufftExec(plan->plan_slow,d_data.data(),d_data.data(),-flag);
  #elif defined(FFT_KOKKOS_HIPFFT)
    hipfftExec(plan->plan_slow,d_data.data(),d_data.data(),-flag);
  #else
    total = plan->total3;
    length = plan->length3;

    d_tmp = typename FFT_AT::t_FFT_DATA_1d(Kokkos::view_alloc("fft_3d:tmp",Kokkos::WithoutInitializing),d_data.extent(0));
    if (flag == 1)
      f = kiss_fft_functor<DeviceType>(d_data,d_tmp,plan->cfg_slow_forward,length);
    else
      f = kiss_fft_functor<DeviceType>(d_data,d_tmp,plan->cfg_slow_backward,length);
    Kokkos::parallel_for(total/length,f);
    d_data = d_tmp;
  #endif

  // post-remap to put data in output format if needed
  // destination is always out

  if (plan->post_plan) {
    d_data_scalar = typename FFT_AT::t_FFT_SCALAR_1d((FFT_SCALAR*)d_data.data(),d_data.size()*2);
    d_out_scalar = typename FFT_AT::t_FFT_SCALAR_1d((FFT_SCALAR*)d_out.data(),d_out.size()*2);
    d_scratch_scalar = typename FFT_AT::t_FFT_SCALAR_1d((FFT_SCALAR*)plan->d_scratch.data(),plan->d_scratch.size()*2);

    remapKK->remap_3d_kokkos(d_data_scalar, d_out_scalar,
             d_scratch_scalar, plan->post_plan);
    }

  // scaling if required

  if (flag == -1 && plan->scaled) {
    FFT_SCALAR norm = plan->norm;
    int num = plan->normnum;

    norm_functor<DeviceType> f(d_out,norm);
    Kokkos::parallel_for(num,f);
  }
}

/* ----------------------------------------------------------------------
   Create plan for performing a 3d FFT

   Arguments:
   comm                 MPI communicator for the P procs which own the d_data
   nfast,nmid,nslow     size of global 3d matrix
   in_ilo,in_ihi        input bounds of data I own in fast index
   in_jlo,in_jhi        input bounds of data I own in mid index
   in_klo,in_khi        input bounds of data I own in slow index
   out_ilo,out_ihi      output bounds of data I own in fast index
   out_jlo,out_jhi      output bounds of data I own in mid index
   out_klo,out_khi      output bounds of data I own in slow index
   scaled               0 = no scaling of result, 1 = scaling
   permute              permutation in storage order of indices on output
                          0 = no permutation
                          1 = permute once = mid->fast, slow->mid, fast->slow
                          2 = permute twice = slow->fast, fast->mid, mid->slow
   nbuf                 returns size of internal storage buffers used by FFT
   usecollective        use collective MPI operations for remapping data
   usegpu_aware         use GPU-Aware MPI or not
------------------------------------------------------------------------- */

template<class DeviceType>
struct fft_plan_3d_kokkos<DeviceType>* FFT3dKokkos<DeviceType>::fft_3d_create_plan_kokkos (
       MPI_Comm comm, int nfast, int nmid, int nslow,
       int in_ilo, int in_ihi, int in_jlo, int in_jhi,
       int in_klo, int in_khi,
       int out_ilo, int out_ihi, int out_jlo, int out_jhi,
       int out_klo, int out_khi,
       int scaled, int permute, int *nbuf, int usecollective,
       int nthreads, int usegpu_aware)
{
  struct fft_plan_3d_kokkos<DeviceType> *plan;
  int me,nprocs;
  int flag,remapflag;
  int first_ilo,first_ihi,first_jlo,first_jhi,first_klo,first_khi;
  int second_ilo,second_ihi,second_jlo,second_jhi,second_klo,second_khi;
  int third_ilo,third_ihi,third_jlo,third_jhi,third_klo,third_khi;
  int out_size,first_size,second_size,third_size,copy_size,scratch_size;
  int np1,np2,ip1,ip2;

  // query MPI info

  MPI_Comm_rank(comm,&me);
  MPI_Comm_size(comm,&nprocs);

  // compute division of procs in 2 dimensions not on-processor

  bifactor(nprocs,&np1,&np2);
  ip1 = me % np1;
  ip2 = me/np1;

  // allocate memory for plan data struct

  plan = new struct fft_plan_3d_kokkos<DeviceType>;
  remapKK = new RemapKokkos3d<DeviceType>(sparta);
  if (plan == nullptr) return nullptr;

  // remap from initial distribution to layout needed for 1st set of 1d FFTs
  // not needed if all procs own entire fast axis initially
  // first indices = distribution after 1st set of FFTs

  if (in_ilo == 0 && in_ihi == nfast-1) flag = 0;
  else flag = 1;

  MPI_Allreduce(&flag,&remapflag,1,MPI_INT,MPI_MAX,comm);

  if (remapflag == 0) {
    first_ilo = in_ilo;
    first_ihi = in_ihi;
    first_jlo = in_jlo;
    first_jhi = in_jhi;
    first_klo = in_klo;
    first_khi = in_khi;
    plan->pre_plan = nullptr;
  }
  else {
    first_ilo = 0;
    first_ihi = nfast - 1;
    first_jlo = ip1*nmid/np1;
    first_jhi = (ip1+1)*nmid/np1 - 1;
    first_klo = ip2*nslow/np2;
    first_khi = (ip2+1)*nslow/np2 - 1;
    plan->pre_plan =
      remapKK->remap_3d_create_plan_kokkos(comm,in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
                           first_ilo,first_ihi,first_jlo,first_jhi,
                           first_klo,first_khi,2,0,0,FFT_PRECISION,
                           usecollective,usegpu_aware);
    if (plan->pre_plan == nullptr) return nullptr;
  }

  // 1d FFTs along fast axis

  plan->length1 = nfast;
  plan->total1 = nfast * (first_jhi-first_jlo+1) * (first_khi-first_klo+1);

  // remap from 1st to 2nd FFT
  // choose which axis is split over np1 vs np2 to minimize communication
  // second indices = distribution after 2nd set of FFTs

  second_ilo = ip1*nfast/np1;
  second_ihi = (ip1+1)*nfast/np1 - 1;
  second_jlo = 0;
  second_jhi = nmid - 1;
  second_klo = ip2*nslow/np2;
  second_khi = (ip2+1)*nslow/np2 - 1;
  plan->mid1_plan =
      remapKK->remap_3d_create_plan_kokkos(comm,
                           first_ilo,first_ihi,first_jlo,first_jhi,
                           first_klo,first_khi,
                           second_ilo,second_ihi,second_jlo,second_jhi,
                           second_klo,second_khi,2,1,0,FFT_PRECISION,
                           usecollective,usegpu_aware);
  if (plan->mid1_plan == nullptr) return nullptr;

  // 1d FFTs along mid axis

  plan->length2 = nmid;
  plan->total2 = (second_ihi-second_ilo+1) * nmid * (second_khi-second_klo+1);

  // remap from 2nd to 3rd FFT
  // if final distribution is permute=2 with all procs owning entire slow axis
  //   then this remapping goes directly to final distribution
  //  third indices = distribution after 3rd set of FFTs

  if (permute == 2 && out_klo == 0 && out_khi == nslow-1) flag = 0;
  else flag = 1;

  MPI_Allreduce(&flag,&remapflag,1,MPI_INT,MPI_MAX,comm);

  if (remapflag == 0) {
    third_ilo = out_ilo;
    third_ihi = out_ihi;
    third_jlo = out_jlo;
    third_jhi = out_jhi;
    third_klo = out_klo;
    third_khi = out_khi;
  }
  else {
    third_ilo = ip1*nfast/np1;
    third_ihi = (ip1+1)*nfast/np1 - 1;
    third_jlo = ip2*nmid/np2;
    third_jhi = (ip2+1)*nmid/np2 - 1;
    third_klo = 0;
    third_khi = nslow - 1;
  }

  plan->mid2_plan =
    remapKK->remap_3d_create_plan_kokkos(comm,
                         second_jlo,second_jhi,second_klo,second_khi,
                         second_ilo,second_ihi,
                         third_jlo,third_jhi,third_klo,third_khi,
                         third_ilo,third_ihi,2,1,0,FFT_PRECISION,
                         usecollective,usegpu_aware);
  if (plan->mid2_plan == nullptr) return nullptr;

  // 1d FFTs along slow axis

  plan->length3 = nslow;
  plan->total3 = (third_ihi-third_ilo+1) * (third_jhi-third_jlo+1) * nslow;

  // remap from 3rd FFT to final distribution
  //  not needed if permute = 2 and third indices = out indices on all procs

  if (permute == 2 &&
      out_ilo == third_ilo && out_ihi == third_ihi &&
      out_jlo == third_jlo && out_jhi == third_jhi &&
      out_klo == third_klo && out_khi == third_khi) flag = 0;
  else flag = 1;

  MPI_Allreduce(&flag,&remapflag,1,MPI_INT,MPI_MAX,comm);

  if (remapflag == 0)
    plan->post_plan = nullptr;
  else {
    plan->post_plan =
      remapKK->remap_3d_create_plan_kokkos(comm,
                           third_klo,third_khi,third_ilo,third_ihi,
                           third_jlo,third_jhi,
                           out_klo,out_khi,out_ilo,out_ihi,
                           out_jlo,out_jhi,2,(permute+1)%3,0,FFT_PRECISION,
                           usecollective,usegpu_aware);
    if (plan->post_plan == nullptr) return nullptr;
  }

  // configure plan memory pointers and allocate work space
  // out_size = amount of memory given to FFT by user
  // first/second/third_size = amount of memory needed after
  //                           pre,mid1,mid2 remaps
  // copy_size = amount needed internally for extra copy of data
  // scratch_size = amount needed internally for remap scratch space
  // for each remap:
  //   out space used for result if big enough, else require copy buffer
  //   accumulate largest required remap scratch space

  out_size = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1) * (out_khi-out_klo+1);
  first_size = (first_ihi-first_ilo+1) * (first_jhi-first_jlo+1) *
    (first_khi-first_klo+1);
  second_size = (second_ihi-second_ilo+1) * (second_jhi-second_jlo+1) *
    (second_khi-second_klo+1);
  third_size = (third_ihi-third_ilo+1) * (third_jhi-third_jlo+1) *
    (third_khi-third_klo+1);

  copy_size = 0;
  scratch_size = 0;

  if (plan->pre_plan) {
    if (first_size <= out_size)
      plan->pre_target = 0;
    else {
      plan->pre_target = 1;
      copy_size = MAX(copy_size,first_size);
    }
    scratch_size = MAX(scratch_size,first_size);
  }

  if (plan->mid1_plan) {
    if (second_size <= out_size)
      plan->mid1_target = 0;
    else {
      plan->mid1_target = 1;
      copy_size = MAX(copy_size,second_size);
    }
    scratch_size = MAX(scratch_size,second_size);
  }

  if (plan->mid2_plan) {
    if (third_size <= out_size)
      plan->mid2_target = 0;
    else {
      plan->mid2_target = 1;
      copy_size = MAX(copy_size,third_size);
    }
    scratch_size = MAX(scratch_size,third_size);
  }

  if (plan->post_plan)
    scratch_size = MAX(scratch_size,out_size);

  *nbuf = copy_size + scratch_size;

  if (copy_size) {
    plan->d_copy = typename FFT_AT::t_FFT_DATA_1d("fft3d:copy",copy_size);
  }

  if (scratch_size) {
    plan->d_scratch = typename FFT_AT::t_FFT_DATA_1d("fft3d:scratch",scratch_size);
  }

  // system specific pre-computation of 1d FFT coeffs
  // and scaling normalization

#if defined(FFT_KOKKOS_MKL)
  DftiCreateDescriptor( &(plan->handle_fast), FFT_KOKKOS_MKL_PREC, DFTI_COMPLEX, 1,
                        (MKL_LONG)nfast);
  DftiSetValue(plan->handle_fast, DFTI_NUMBER_OF_TRANSFORMS,
               (MKL_LONG)plan->total1/nfast);
  DftiSetValue(plan->handle_fast, DFTI_PLACEMENT,DFTI_INPLACE);
  DftiSetValue(plan->handle_fast, DFTI_INPUT_DISTANCE, (MKL_LONG)nfast);
  DftiSetValue(plan->handle_fast, DFTI_OUTPUT_DISTANCE, (MKL_LONG)nfast);
#if defined(FFT_KOKKOS_MKL_THREADS)
  DftiSetValue(plan->handle_fast, DFTI_NUMBER_OF_USER_THREADS, nthreads);
#endif
  DftiCommitDescriptor(plan->handle_fast);

  DftiCreateDescriptor( &(plan->handle_mid), FFT_KOKKOS_MKL_PREC, DFTI_COMPLEX, 1,
                        (MKL_LONG)nmid);
  DftiSetValue(plan->handle_mid, DFTI_NUMBER_OF_TRANSFORMS,
               (MKL_LONG)plan->total2/nmid);
  DftiSetValue(plan->handle_mid, DFTI_PLACEMENT,DFTI_INPLACE);
  DftiSetValue(plan->handle_mid, DFTI_INPUT_DISTANCE, (MKL_LONG)nmid);
  DftiSetValue(plan->handle_mid, DFTI_OUTPUT_DISTANCE, (MKL_LONG)nmid);
#if defined(FFT_KOKKOS_MKL_THREADS)
  DftiSetValue(plan->handle_mid, DFTI_NUMBER_OF_USER_THREADS, nthreads);
#endif
  DftiCommitDescriptor(plan->handle_mid);

  DftiCreateDescriptor( &(plan->handle_slow), FFT_KOKKOS_MKL_PREC, DFTI_COMPLEX, 1,
                        (MKL_LONG)nslow);
  DftiSetValue(plan->handle_slow, DFTI_NUMBER_OF_TRANSFORMS,
               (MKL_LONG)plan->total3/nslow);
  DftiSetValue(plan->handle_slow, DFTI_PLACEMENT,DFTI_INPLACE);
  DftiSetValue(plan->handle_slow, DFTI_INPUT_DISTANCE, (MKL_LONG)nslow);
  DftiSetValue(plan->handle_slow, DFTI_OUTPUT_DISTANCE, (MKL_LONG)nslow);
#if defined(FFT_KOKKOS_MKL_THREADS)
  DftiSetValue(plan->handle_slow, DFTI_NUMBER_OF_USER_THREADS, nthreads);
#endif
  DftiCommitDescriptor(plan->handle_slow);

#elif defined(FFT_KOKKOS_FFTW3)

#if defined (FFT_KOKKOS_FFTW_THREADS)
  if (nthreads > 1) {
    FFTW_API(init_threads)();
    FFTW_API(plan_with_nthreads)(nthreads);
  }
#endif

  plan->plan_fast_forward =
    FFTW_API(plan_many_dft)(1, &nfast,plan->total1/plan->length1,
                       nullptr,&nfast,1,plan->length1,
                       nullptr,&nfast,1,plan->length1,
                       FFTW_FORWARD,FFTW_ESTIMATE);

  plan->plan_fast_backward =
    FFTW_API(plan_many_dft)(1, &nfast,plan->total1/plan->length1,
                       nullptr,&nfast,1,plan->length1,
                       nullptr,&nfast,1,plan->length1,
                       FFTW_BACKWARD,FFTW_ESTIMATE);

  plan->plan_mid_forward =
    FFTW_API(plan_many_dft)(1, &nmid,plan->total2/plan->length2,
                       nullptr,&nmid,1,plan->length2,
                       nullptr,&nmid,1,plan->length2,
                       FFTW_FORWARD,FFTW_ESTIMATE);

  plan->plan_mid_backward =
    FFTW_API(plan_many_dft)(1, &nmid,plan->total2/plan->length2,
                       nullptr,&nmid,1,plan->length2,
                       nullptr,&nmid,1,plan->length2,
                       FFTW_BACKWARD,FFTW_ESTIMATE);


  plan->plan_slow_forward =
    FFTW_API(plan_many_dft)(1, &nslow,plan->total3/plan->length3,
                       nullptr,&nslow,1,plan->length3,
                       nullptr,&nslow,1,plan->length3,
                       FFTW_FORWARD,FFTW_ESTIMATE);

  plan->plan_slow_backward =
    FFTW_API(plan_many_dft)(1, &nslow,plan->total3/plan->length3,
                       nullptr,&nslow,1,plan->length3,
                       nullptr,&nslow,1,plan->length3,
                       FFTW_BACKWARD,FFTW_ESTIMATE);

#elif defined(FFT_KOKKOS_CUFFT)

  cufftPlanMany(&(plan->plan_fast), 1, &nfast,
    &nfast,1,plan->length1,
    &nfast,1,plan->length1,
    CUFFT_TYPE,plan->total1/plan->length1);

  cufftPlanMany(&(plan->plan_mid), 1, &nmid,
    &nmid,1,plan->length2,
    &nmid,1,plan->length2,
    CUFFT_TYPE,plan->total2/plan->length2);

  cufftPlanMany(&(plan->plan_slow), 1, &nslow,
    &nslow,1,plan->length3,
    &nslow,1,plan->length3,
    CUFFT_TYPE,plan->total3/plan->length3);

#elif defined(FFT_KOKKOS_HIPFFT)

  hipfftPlanMany(&(plan->plan_fast), 1, &nfast,
    &nfast,1,plan->length1,
    &nfast,1,plan->length1,
    HIPFFT_KOKKOS_TYPE,plan->total1/plan->length1);

  hipfftPlanMany(&(plan->plan_mid), 1, &nmid,
    &nmid,1,plan->length2,
    &nmid,1,plan->length2,
    HIPFFT_KOKKOS_TYPE,plan->total2/plan->length2);

  hipfftPlanMany(&(plan->plan_slow), 1, &nslow,
    &nslow,1,plan->length3,
    &nslow,1,plan->length3,
    HIPFFT_KOKKOS_TYPE,plan->total3/plan->length3);

#else  /* FFT_KOKKOS_KISS */

  kissfftKK = new KissFFTKokkos<DeviceType>();

  plan->cfg_fast_forward = KissFFTKokkos<DeviceType>::kiss_fft_alloc_kokkos(nfast,0,nullptr,nullptr);
  plan->cfg_fast_backward = KissFFTKokkos<DeviceType>::kiss_fft_alloc_kokkos(nfast,1,nullptr,nullptr);

  if (nmid == nfast) {
    plan->cfg_mid_forward = plan->cfg_fast_forward;
    plan->cfg_mid_backward = plan->cfg_fast_backward;
  }
  else {
    plan->cfg_mid_forward = KissFFTKokkos<DeviceType>::kiss_fft_alloc_kokkos(nmid,0,nullptr,nullptr);
    plan->cfg_mid_backward = KissFFTKokkos<DeviceType>::kiss_fft_alloc_kokkos(nmid,1,nullptr,nullptr);
  }

  if (nslow == nfast) {
    plan->cfg_slow_forward = plan->cfg_fast_forward;
    plan->cfg_slow_backward = plan->cfg_fast_backward;
  }
  else if (nslow == nmid) {
    plan->cfg_slow_forward = plan->cfg_mid_forward;
    plan->cfg_slow_backward = plan->cfg_mid_backward;
  }
  else {
    plan->cfg_slow_forward = KissFFTKokkos<DeviceType>::kiss_fft_alloc_kokkos(nslow,0,nullptr,nullptr);
    plan->cfg_slow_backward = KissFFTKokkos<DeviceType>::kiss_fft_alloc_kokkos(nslow,1,nullptr,nullptr);
  }

#endif

  if (scaled == 0)
    plan->scaled = 0;
  else {
    plan->scaled = 1;
    plan->norm = 1.0/(nfast*nmid*nslow);
    plan->normnum = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1) *
      (out_khi-out_klo+1);
  }

  return plan;
}

/* ----------------------------------------------------------------------
   Destroy a 3d fft plan
------------------------------------------------------------------------- */

template<class DeviceType>
void FFT3dKokkos<DeviceType>::fft_3d_destroy_plan_kokkos(struct fft_plan_3d_kokkos<DeviceType> *plan)
{
  if (plan->pre_plan) remapKK->remap_3d_destroy_plan_kokkos(plan->pre_plan);
  if (plan->mid1_plan) remapKK->remap_3d_destroy_plan_kokkos(plan->mid1_plan);
  if (plan->mid2_plan) remapKK->remap_3d_destroy_plan_kokkos(plan->mid2_plan);
  if (plan->post_plan) remapKK->remap_3d_destroy_plan_kokkos(plan->post_plan);

#if defined(FFT_KOKKOS_MKL)
  DftiFreeDescriptor(&(plan->handle_fast));
  DftiFreeDescriptor(&(plan->handle_mid));
  DftiFreeDescriptor(&(plan->handle_slow));
#elif defined(FFT_KOKKOS_FFTW3)
  FFTW_API(destroy_plan)(plan->plan_slow_forward);
  FFTW_API(destroy_plan)(plan->plan_slow_backward);
  FFTW_API(destroy_plan)(plan->plan_mid_forward);
  FFTW_API(destroy_plan)(plan->plan_mid_backward);
  FFTW_API(destroy_plan)(plan->plan_fast_forward);
  FFTW_API(destroy_plan)(plan->plan_fast_backward);

#if defined (FFT_KOKKOS_FFTW_THREADS)
  FFTW_API(cleanup_threads)();
#endif

#elif defined (FFT_KOKKOS_KISS)
  delete kissfftKK;
#endif

  delete plan;
  delete remapKK;
}

/* ----------------------------------------------------------------------
   divide n into 2 factors of as equal size as possible
------------------------------------------------------------------------- */

template<class DeviceType>
void FFT3dKokkos<DeviceType>::bifactor(int n, int *factor1, int *factor2)
{
  int n1,n2,facmax;

  facmax = static_cast<int> (sqrt((double) n));

  for (n1 = facmax; n1 > 0; n1--) {
    n2 = n/n1;
    if (n1*n2 == n) {
      *factor1 = n1;
      *factor2 = n2;
      return;
    }
  }
}

/* ----------------------------------------------------------------------
   perform just the 1d FFTs needed by a 3d FFT, no data movement
   used for timing purposes

   Arguments:
   in           starting address of input data on this proc, all set to 0.0
   nsize        size of in
   flag         1 for forward FFT, -1 for backward FFT
   plan         plan returned by previous call to fft_3d_create_plan
------------------------------------------------------------------------- */

template<class DeviceType>
void FFT3dKokkos<DeviceType>::fft_3d_1d_only_kokkos(typename FFT_AT::t_FFT_DATA_1d d_data, int nsize, int flag,
                    struct fft_plan_3d_kokkos<DeviceType> *plan)
{
  // total = size of data needed in each dim
  // length = length of 1d FFT in each dim
  // total/length = # of 1d FFTs in each dim
  // if total > nsize, limit # of 1d FFTs to available size of data

  int total1 = plan->total1;
  int length1 = plan->length1;
  int total2 = plan->total2;
  int length2 = plan->length2;
  int total3 = plan->total3;
  int length3 = plan->length3;

  // fftw3 and Dfti in MKL encode the number of transforms
  // into the plan, so we cannot operate on a smaller data set

#if defined(FFT_KOKKOS_MKL) || defined(FFT_KOKKOS_FFTW3)
  if ((total1 > nsize) || (total2 > nsize) || (total3 > nsize))
    return;
#endif
  if (total1 > nsize) total1 = (nsize/length1) * length1;
  if (total2 > nsize) total2 = (nsize/length2) * length2;
  if (total3 > nsize) total3 = (nsize/length3) * length3;

  // perform 1d FFTs in each of 3 dimensions
  // data is just an array of 0.0

#if defined(FFT_KOKKOS_MKL)
  if (flag == -1) {
    DftiComputeForward(plan->handle_fast,d_data.data());
    DftiComputeForward(plan->handle_mid,d_data.data());
    DftiComputeForward(plan->handle_slow,d_data.data());
  } else {
    DftiComputeBackward(plan->handle_fast,d_data.data());
    DftiComputeBackward(plan->handle_mid,d_data.data());
    DftiComputeBackward(plan->handle_slow,d_data.data());
  }
#elif defined(FFT_KOKKOS_FFTW3)
  if (flag == -1) {
    FFTW_API(execute_dft)(plan->plan_fast_forward,(FFT_KOKKOS_DATA*)d_data.data(),(FFT_KOKKOS_DATA*)d_data.data());
    FFTW_API(execute_dft)(plan->plan_mid_forward,(FFT_KOKKOS_DATA*)d_data.data(),(FFT_KOKKOS_DATA*)d_data.data());
    FFTW_API(execute_dft)(plan->plan_slow_forward,(FFT_KOKKOS_DATA*)d_data.data(),(FFT_KOKKOS_DATA*)d_data.data());
  } else {
    FFTW_API(execute_dft)(plan->plan_fast_backward,(FFT_KOKKOS_DATA*)d_data.data(),(FFT_KOKKOS_DATA*)d_data.data());
    FFTW_API(execute_dft)(plan->plan_mid_backward,(FFT_KOKKOS_DATA*)d_data.data(),(FFT_KOKKOS_DATA*)d_data.data());
    FFTW_API(execute_dft)(plan->plan_slow_backward,(FFT_KOKKOS_DATA*)d_data.data(),(FFT_KOKKOS_DATA*)d_data.data());
  }
#elif defined(FFT_KOKKOS_CUFFT)
  cufftExec(plan->plan_fast,d_data.data(),d_data.data(),-flag);
  cufftExec(plan->plan_mid,d_data.data(),d_data.data(),-flag);
  cufftExec(plan->plan_slow,d_data.data(),d_data.data(),-flag);
#elif defined(FFT_KOKKOS_HIPFFT)
  hipfftExec(plan->plan_fast,d_data.data(),d_data.data(),-flag);
  hipfftExec(plan->plan_mid,d_data.data(),d_data.data(),-flag);
  hipfftExec(plan->plan_slow,d_data.data(),d_data.data(),-flag);
#else
  kiss_fft_functor<DeviceType> f;
    typename FFT_AT::t_FFT_DATA_1d d_tmp =
     typename FFT_AT::t_FFT_DATA_1d(Kokkos::view_alloc("fft_3d:tmp",Kokkos::WithoutInitializing),d_data.extent(0));
  if (flag == -1) {
    f = kiss_fft_functor<DeviceType>(d_data,d_tmp,plan->cfg_fast_forward,length1);
    Kokkos::parallel_for(total1/length1,f);

    f = kiss_fft_functor<DeviceType>(d_data,d_tmp,plan->cfg_mid_forward,length2);
    Kokkos::parallel_for(total2/length2,f);

    f = kiss_fft_functor<DeviceType>(d_data,d_tmp,plan->cfg_slow_forward,length3);
    Kokkos::parallel_for(total3/length3,f);
  } else {
    f = kiss_fft_functor<DeviceType>(d_data,d_tmp,plan->cfg_fast_backward,length1);
    Kokkos::parallel_for(total1/length1,f);

    f = kiss_fft_functor<DeviceType>(d_data,d_tmp,plan->cfg_mid_backward,length2);
    Kokkos::parallel_for(total2/length2,f);

    f = kiss_fft_functor<DeviceType>(d_data,d_tmp,plan->cfg_slow_backward,length3);
    Kokkos::parallel_for(total3/length3,f);
  }
#endif

  // scaling if required
  // limit num to size of data

  if (flag == 1 && plan->scaled) {
    FFT_SCALAR norm = plan->norm;
    int num = MIN(plan->normnum,nsize);

    norm_functor<DeviceType> f(d_data,norm);
    Kokkos::parallel_for(num,f);
  }
}

namespace SPARTA_NS {
template class FFT3dKokkos<SPADeviceType>;
#ifdef SPARTA_KOKKOS_GPU
template class FFT3dKokkos<SPAHostType>;
#endif
}

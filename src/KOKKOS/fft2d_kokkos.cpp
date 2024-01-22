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

#include "fft2d_kokkos.h"

#include "error.h"
#include "kokkos.h"
#include "remap2d_kokkos.h"

#include <cmath>

using namespace SPARTA_NS;

#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FFT2dKokkos<DeviceType>::FFT2dKokkos(SPARTA *sparta, MPI_Comm comm, int nfast, int nslow,
             int in_ilo, int in_ihi, int in_jlo, int in_jhi,
             int out_ilo, int out_ihi, int out_jlo, int out_jhi,
             int scaled, int permute, int *nbuf, int usecollective,
             int usegpu_aware) :
  Pointers(sparta)
{
  int nthreads = sparta->kokkos->nthreads;
  int ngpus = sparta->kokkos->ngpus;
  ExecutionSpace execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

#if defined(FFT_MKL)
  if (ngpus > 0 && execution_space == Device)
    sparta->error->all(FLERR,"Cannot use the MKL library with Kokkos on GPUs");
#elif defined(FFT_FFTW3)
  if (ngpus > 0 && execution_space == Device)
    sparta->error->all(FLERR,"Cannot use the FFTW library with Kokkos on GPUs");
#elif defined(FFT_CUFFT)
  if (ngpus > 0 && execution_space == Host)
    sparta->error->all(FLERR,"Cannot use the cuFFT library with Kokkos on the host CPUs");
#elif defined(FFT_HIPFFT)
  if (ngpus > 0 && execution_space == Host)
    sparta->error->all(FLERR,"Cannot use the hipFFT library with Kokkos on the host CPUs");

#elif defined(FFT_KISSFFT)
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

  plan = fft_2d_create_plan_kokkos(comm,nfast,nslow,
                            in_ilo,in_ihi,in_jlo,in_jhi,
                            out_ilo,out_ihi,out_jlo,out_jhi,
                            scaled,permute,nbuf,usecollective,nthreads,usegpu_aware);
  if (plan == nullptr) error->one(FLERR,"Could not create 2d FFT plan");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FFT2dKokkos<DeviceType>::~FFT2dKokkos()
{
  fft_2d_destroy_plan_kokkos(plan);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FFT2dKokkos<DeviceType>::compute(typename FFT_AT::t_FFT_SCALAR_1d d_in, typename FFT_AT::t_FFT_SCALAR_1d d_out, int flag)
{
  typename FFT_AT::t_FFT_DATA_1d d_in_data((FFT_DATA_POINTER)d_in.data(),d_in.size()/2);
  typename FFT_AT::t_FFT_DATA_1d d_out_data((FFT_DATA_POINTER)d_out.data(),d_out.size()/2);

  fft_2d_kokkos(d_in_data,d_out_data,flag,plan);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FFT2dKokkos<DeviceType>::timing1d(typename FFT_AT::t_FFT_SCALAR_1d d_in, int nsize, int flag)
{
  typename FFT_AT::t_FFT_DATA_1d d_in_data((FFT_DATA_POINTER)d_in.data(),d_in.size()/2);

  fft_2d_1d_only_kokkos(d_in_data,nsize,flag,plan);
}

/* ----------------------------------------------------------------------
   Data layout for 2d FFTs:

   data set of Nfast x Nslow elements is owned by P procs
   on input, each proc owns a subsection of the elements
   on output, each proc will own a (possibly different) subsection
   my subsection must not overlap with any other proc's subsection,
     i.e. the union of all proc's input (or output) subsections must
     exactly tile the global Nfast x Nslow data set
   when called from C, all subsection indices are
     C-style from 0 to N-1 where N = Nfast or Nslow
   when called from F77, all subsection indices are
     F77-style from 1 to N where N = Nfast or Nslow
   a proc can own 0 elements on input or output
     by specifying hi index < lo index
   on both input and output, data is stored contiguously on a processor
     with a fast-varying and slow-varying index
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Perform 2d FFT

   Arguments:
   in           starting address of input data on this proc
   out          starting address of where output data for this proc
                  will be placed (can be same as in)
   flag         1 for forward FFT, -1 for inverse FFT
   plan         plan returned by previous call to fft_2d_create_plan
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
#if defined(FFT_FFTW3) || defined(FFT_CUFFT) || defined(FFT_HIPFFT)
    FFT_SCALAR* out_ptr = (FFT_SCALAR *)(d_out.data()+i);
    *(out_ptr++) *= norm;
    *(out_ptr++) *= norm;
#elif defined(FFT_MKL)
    d_out(i) *= norm;
#else // FFT_KISS
    d_out(i).re *= norm;
    d_out(i).im *= norm;
#endif
  }
};

#ifdef FFT_KISSFFT
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
void FFT2dKokkos<DeviceType>::fft_2d_kokkos(typename FFT_AT::t_FFT_DATA_1d d_in, typename FFT_AT::t_FFT_DATA_1d d_out, int flag, struct fft_plan_2d_kokkos<DeviceType> *plan)
{
  int total,length;
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

    remapKK->remap_2d_kokkos(d_in_scalar, d_copy_scalar,
             d_scratch_scalar, plan->pre_plan);

    d_data = d_copy;
  } else d_data = d_in;

  // 1d FFTs along fast axis

  total = plan->total1;
  length = plan->length1;

  #if defined(FFT_MKL)
    if (flag == 1)
      DftiComputeForward(plan->handle_fast,d_data.data());
    else
      DftiComputeBackward(plan->handle_fast,d_data.data());
  #elif defined(FFT_FFTW3)
    if (flag == 1)
      FFTW_API(execute_dft)(plan->plan_fast_forward,(FFT_DATA*)d_data.data(),(FFT_DATA*)d_data.data());
    else
      FFTW_API(execute_dft)(plan->plan_fast_backward,(FFT_DATA*)d_data.data(),(FFT_DATA*)d_data.data());
  #elif defined(FFT_CUFFT)
    cufftExec(plan->plan_fast,d_data.data(),d_data.data(),-flag);
  #elif defined(FFT_HIPFFT)
    hipfftExec(plan->plan_fast,d_data.data(),d_data.data(),-flag);
  #else
    typename FFT_AT::t_FFT_DATA_1d d_tmp =
     typename FFT_AT::t_FFT_DATA_1d(Kokkos::view_alloc("fft_2d:tmp",Kokkos::WithoutInitializing),d_data.extent(0));
    kiss_fft_functor<DeviceType> f;
    if (flag == 1)
      f = kiss_fft_functor<DeviceType>(d_data,d_tmp,plan->cfg_fast_forward,length);
    else
      f = kiss_fft_functor<DeviceType>(d_data,d_tmp,plan->cfg_fast_backward,length);
    Kokkos::parallel_for(total/length,f);
    d_data = d_tmp;
  #endif

  // mid-remap to prepare for 2nd FFTs
  // copy = loc for remap result

  if (plan->mid_target == 0) d_copy = d_out;
  else d_copy = plan->d_copy;

  d_data_scalar = typename FFT_AT::t_FFT_SCALAR_1d((FFT_SCALAR*)d_data.data(),d_data.size()*2);
  d_copy_scalar = typename FFT_AT::t_FFT_SCALAR_1d((FFT_SCALAR*)d_copy.data(),d_copy.size()*2);
  d_scratch_scalar = typename FFT_AT::t_FFT_SCALAR_1d((FFT_SCALAR*)plan->d_scratch.data(),plan->d_scratch.size()*2);

  remapKK->remap_2d_kokkos(d_data_scalar, d_copy_scalar,
           d_scratch_scalar, plan->mid_plan);

  d_data = d_copy;

  // 1d FFTs along slow axis

  total = plan->total2;
  length = plan->length2;

  #if defined(FFT_MKL)
    if (flag == 1)
      DftiComputeForward(plan->handle_slow,d_data.data());
    else
      DftiComputeBackward(plan->handle_slow,d_data.data());
  #elif defined(FFT_FFTW3)
    if (flag == 1)
      FFTW_API(execute_dft)(plan->plan_slow_forward,(FFT_DATA*)d_data.data(),(FFT_DATA*)d_data.data());
    else
      FFTW_API(execute_dft)(plan->plan_slow_backward,(FFT_DATA*)d_data.data(),(FFT_DATA*)d_data.data());
  #elif defined(FFT_CUFFT)
    cufftExec(plan->plan_slow,d_data.data(),d_data.data(),-flag);
  #elif defined(FFT_HIPFFT)
    hipfftExec(plan->plan_slow,d_data.data(),d_data.data(),-flag);
  #else
    d_tmp = typename FFT_AT::t_FFT_DATA_1d(Kokkos::view_alloc("fft_2d:tmp",Kokkos::WithoutInitializing),d_data.extent(0));
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

    remapKK->remap_2d_kokkos(d_data_scalar, d_out_scalar,
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
   Create plan for performing a 2d FFT

   Arguments:
   comm                 MPI communicator for the P procs which own the d_data
   nfast,nslow          size of global 2d matrix
   in_ilo,in_ihi        input bounds of data I own in fast index
   in_jlo,in_jhi        input bounds of data I own in slow index
   out_ilo,out_ihi      output bounds of data I own in fast index
   out_jlo,out_jhi      output bounds of data I own in slow index
   scaled               0 = no scaling of result, 1 = scaling
   permute              permutation in storage order of indices on output
                          0 = no permutation
                          1 = permute = slow->fast, fast->slow
   nbuf                 returns size of internal storage buffers used by FFT
   usecollective        use collective MPI operations for remapping data
   usegpu_aware        use GPU-Aware MPI or not
------------------------------------------------------------------------- */

template<class DeviceType>
struct fft_plan_2d_kokkos<DeviceType>* FFT2dKokkos<DeviceType>::fft_2d_create_plan_kokkos (
       MPI_Comm comm, int nfast, int nslow,
       int in_ilo, int in_ihi, int in_jlo, int in_jhi,
       int out_ilo, int out_ihi, int out_jlo, int out_jhi,
       int scaled, int permute, int *nbuf, int usecollective,
       int nthreads, int usegpu_aware)
{
  struct fft_plan_2d_kokkos<DeviceType> *plan;
  int me,nprocs;
  int flag,remapflag;
  int first_ilo,first_ihi,first_jlo,first_jhi;
  int second_ilo,second_ihi,second_jlo,second_jhi;
  int out_size,first_size,second_size,copy_size,scratch_size;
  int np1,np2,ip1,ip2;

  // query MPI info

  MPI_Comm_rank(comm,&me);
  MPI_Comm_size(comm,&nprocs);

  // allocate memory for plan data struct

  plan = new struct fft_plan_2d_kokkos<DeviceType>;
  remapKK = new RemapKokkos2d<DeviceType>(sparta);
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
    plan->pre_plan = nullptr;
  }
  else {
    first_ilo = 0;
    first_ihi = nfast - 1;
    first_jlo = me*nslow/nprocs;
    first_jhi = (me+1)*nslow/nprocs - 1;
    plan->pre_plan =
      remapKK->remap_2d_create_plan_kokkos(comm,in_ilo,in_ihi,in_jlo,in_jhi,
                  first_ilo,first_ihi,first_jlo,first_jhi,
                           2,0,0,FFT_PRECISION,
                           usecollective,usegpu_aware);
    if (plan->pre_plan == nullptr) return nullptr;
  }

  // 1d FFTs along fast axis

  plan->length1 = nfast;
  plan->total1 = nfast * (first_jhi-first_jlo+1);

  // remap from 1st to 2nd FFT
  // if final distribution is permute=1 with all procs owning entire slow axis
  //   then this remapping goes directly to final distribution
  // second indices = distribution after 2nd set of FFTs

  if (permute == 1 && out_jlo == 0 && out_jhi == nslow-1) flag = 0;
  else flag = 1;

  MPI_Allreduce(&flag,&remapflag,1,MPI_INT,MPI_MAX,comm);

  if (remapflag == 0) {
    second_ilo = out_ilo;
    second_ihi = out_ihi;
    second_jlo = out_jlo;
    second_jhi = out_jhi;
  }
  else {
    second_ilo = me*nfast/nprocs;
    second_ihi = (me+1)*nfast/nprocs - 1;
    second_jlo = 0;
    second_jhi = nslow - 1;
  }

  plan->mid_plan =
    remapKK->remap_2d_create_plan_kokkos(comm,first_ilo,first_ihi,first_jlo,first_jhi,
                                         second_ilo,second_ihi,second_jlo,second_jhi,
                                         FFT_PRECISION,1,0,2,
                                         usecollective,usegpu_aware);
  if (plan->mid_plan == nullptr) return nullptr;

  // 1d FFTs along slow axis
  plan->length2 = nslow;
  plan->total2 = (second_ihi-second_ilo+1) * nslow;

  // remap from 2nd FFT to final distribution
  // not needed if permute = 1 and second indices = out indices on all procs

  if (permute == 1 && out_ilo == second_ilo && out_ihi == second_ihi &&
      out_jlo == second_jlo && out_jhi == second_jhi) flag = 0;
  else flag = 1;

  MPI_Allreduce(&flag,&remapflag,1,MPI_INT,MPI_MAX,comm);

  if (remapflag == 0)
    plan->post_plan = nullptr;
  else {
    plan->post_plan =
      remapKK->remap_2d_create_plan_kokkos(comm,second_jlo,second_jhi,second_ilo,second_ihi,
                  out_jlo,out_jhi,out_ilo,out_ihi,
                  FFT_PRECISION,(permute+1)%2,0,2,
                  usecollective,usegpu_aware);
      if (plan->post_plan == nullptr) return nullptr;
  }

  // configure plan memory pointers and allocate work space
  // out_size = amount of memory given to FFT by user
  // first/second_size = amount of memory needed after
  //                           pre,mid remaps
  // copy_size = amount needed internally for extra copy of data
  // scratch_size = amount needed internally for remap scratch space
  // for each remap:
  //   out space used for result if big enough, else require copy buffer
  //   accumulate largest required remap scratch space

  out_size = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1);
  first_size = (first_ihi-first_ilo+1) * (first_jhi-first_jlo+1);
  second_size = (second_ihi-second_ilo+1) * (second_jhi-second_jlo+1);

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

  if (plan->mid_plan) {
    if (second_size <= out_size)
      plan->mid_target = 0;
    else {
      plan->mid_target = 1;
      copy_size = MAX(copy_size,second_size);
    }
    scratch_size = MAX(scratch_size,second_size);
  }

  if (plan->post_plan)
    scratch_size = MAX(scratch_size,out_size);

  *nbuf = copy_size + scratch_size;

  if (copy_size) {
    plan->d_copy = typename FFT_AT::t_FFT_DATA_1d("fft2d:copy",copy_size);
  }

  if (scratch_size) {
    plan->d_scratch = typename FFT_AT::t_FFT_DATA_1d("fft2d:scratch",scratch_size);
  }

  // system specific pre-computation of 1d FFT coeffs
  // and scaling normalization

#if defined(FFT_MKL)
  DftiCreateDescriptor( &(plan->handle_fast), FFT_MKL_PREC, DFTI_COMPLEX, 1,
                        (MKL_LONG)nfast);
  DftiSetValue(plan->handle_fast, DFTI_NUMBER_OF_TRANSFORMS,
               (MKL_LONG)plan->total1/nfast);
  DftiSetValue(plan->handle_fast, DFTI_PLACEMENT,DFTI_INPLACE);
  DftiSetValue(plan->handle_fast, DFTI_INPUT_DISTANCE, (MKL_LONG)nfast);
  DftiSetValue(plan->handle_fast, DFTI_OUTPUT_DISTANCE, (MKL_LONG)nfast);
#if defined(FFT_MKL_THREADS)
  DftiSetValue(plan->handle_fast, DFTI_NUMBER_OF_USER_THREADS, nthreads);
#endif
  DftiCommitDescriptor(plan->handle_fast);

  DftiCreateDescriptor( &(plan->handle_slow), FFT_MKL_PREC, DFTI_COMPLEX, 1,
                        (MKL_LONG)nslow);
  DftiSetValue(plan->handle_slow, DFTI_NUMBER_OF_TRANSFORMS,
               (MKL_LONG)plan->total2/nslow);
  DftiSetValue(plan->handle_slow, DFTI_PLACEMENT,DFTI_INPLACE);
  DftiSetValue(plan->handle_slow, DFTI_INPUT_DISTANCE, (MKL_LONG)nslow);
  DftiSetValue(plan->handle_slow, DFTI_OUTPUT_DISTANCE, (MKL_LONG)nslow);
#if defined(FFT_MKL_THREADS)
  DftiSetValue(plan->handle_slow, DFTI_NUMBER_OF_USER_THREADS, nthreads);
#endif
  DftiCommitDescriptor(plan->handle_slow);

#elif defined(FFT_FFTW3)

#if defined (FFT_FFTW_THREADS)
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

  plan->plan_slow_forward =
    FFTW_API(plan_many_dft)(1, &nslow,plan->total2/plan->length2,
                       nullptr,&nslow,1,plan->length2,
                       nullptr,&nslow,1,plan->length2,
                       FFTW_FORWARD,FFTW_ESTIMATE);

  plan->plan_slow_backward =
    FFTW_API(plan_many_dft)(1, &nslow,plan->total2/plan->length2,
                       nullptr,&nslow,1,plan->length2,
                       nullptr,&nslow,1,plan->length2,
                       FFTW_BACKWARD,FFTW_ESTIMATE);

#elif defined(FFT_CUFFT)

  cufftPlanMany(&(plan->plan_fast), 1, &nfast,
    &nfast,1,plan->length1,
    &nfast,1,plan->length1,
    CUFFT_TYPE,plan->total1/plan->length1);

  cufftPlanMany(&(plan->plan_slow), 1, &nslow,
    &nslow,1,plan->length2,
    &nslow,1,plan->length2,
    CUFFT_TYPE,plan->total2/plan->length2);

#elif defined(FFT_HIPFFT)

  hipfftPlanMany(&(plan->plan_fast), 1, &nfast,
    &nfast,1,plan->length1,
    &nfast,1,plan->length1,
    HIPFFT_TYPE,plan->total1/plan->length1);

  hipfftPlanMany(&(plan->plan_slow), 1, &nslow,
    &nslow,1,plan->length2,
    &nslow,1,plan->length2,
    HIPFFT_TYPE,plan->total2/plan->length2);

#else  /* FFT_KISS */

  kissfftKK = new KissFFTKokkos<DeviceType>();

  plan->cfg_fast_forward = KissFFTKokkos<DeviceType>::kiss_fft_alloc_kokkos(nfast,0,nullptr,nullptr);
  plan->cfg_fast_backward = KissFFTKokkos<DeviceType>::kiss_fft_alloc_kokkos(nfast,1,nullptr,nullptr);

  if (nslow == nfast) {
    plan->cfg_slow_forward = plan->cfg_fast_forward;
    plan->cfg_slow_backward = plan->cfg_fast_backward;
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
    plan->norm = 1.0/(nfast*nslow);
    plan->normnum = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1);
  }

  return plan;
}

/* ----------------------------------------------------------------------
   Destroy a 2d fft plan
------------------------------------------------------------------------- */

template<class DeviceType>
void FFT2dKokkos<DeviceType>::fft_2d_destroy_plan_kokkos(struct fft_plan_2d_kokkos<DeviceType> *plan)
{
  if (plan->pre_plan) remapKK->remap_2d_destroy_plan_kokkos(plan->pre_plan);
  if (plan->mid_plan) remapKK->remap_2d_destroy_plan_kokkos(plan->mid_plan);
  if (plan->post_plan) remapKK->remap_2d_destroy_plan_kokkos(plan->post_plan);

#if defined(FFT_MKL)
  DftiFreeDescriptor(&(plan->handle_fast));
  DftiFreeDescriptor(&(plan->handle_slow));
#elif defined(FFT_FFTW3)
  FFTW_API(destroy_plan)(plan->plan_slow_forward);
  FFTW_API(destroy_plan)(plan->plan_slow_backward);
  FFTW_API(destroy_plan)(plan->plan_fast_forward);
  FFTW_API(destroy_plan)(plan->plan_fast_backward);

#if defined (FFT_FFTW_THREADS)
  FFTW_API(cleanup_threads)();
#endif

#elif defined (FFT_KISSFFT)
  delete kissfftKK;
#endif

  delete plan;
  delete remapKK;
}

/* ----------------------------------------------------------------------
   perform just the 1d FFTs needed by a 2d FFT, no data movement
   used for timing purposes

   Arguments:
   in           starting address of input data on this proc, all set to 0.0
   nsize        size of in
   flag         1 for forward FFT, -1 for backward FFT
   plan         plan returned by previous call to fft_2d_create_plan
------------------------------------------------------------------------- */

template<class DeviceType>
void FFT2dKokkos<DeviceType>::fft_2d_1d_only_kokkos(typename FFT_AT::t_FFT_DATA_1d d_data, int nsize, int flag,
                    struct fft_plan_2d_kokkos<DeviceType> *plan)
{
  // total = size of data needed in each dim
  // length = length of 1d FFT in each dim
  // total/length = # of 1d FFTs in each dim
  // if total > nsize, limit # of 1d FFTs to available size of data

  int total1 = plan->total1;
  int length1 = plan->length1;
  int total2 = plan->total2;
  int length2 = plan->length2;

  // fftw3 and Dfti in MKL encode the number of transforms
  // into the plan, so we cannot operate on a smaller data set

#if defined(FFT_MKL) || defined(FFT_FFTW3)
  if ((total1 > nsize) || (total2 > nsize))
    return;
#endif
  if (total1 > nsize) total1 = (nsize/length1) * length1;
  if (total2 > nsize) total2 = (nsize/length2) * length2;

  // perform 1d FFTs in each of 2 dimensions
  // data is just an array of 0.0

#if defined(FFT_MKL)
  if (flag == -1) {
    DftiComputeForward(plan->handle_fast,d_data.data());
    DftiComputeForward(plan->handle_slow,d_data.data());
  } else {
    DftiComputeBackward(plan->handle_fast,d_data.data());
    DftiComputeBackward(plan->handle_slow,d_data.data());
  }
#elif defined(FFT_FFTW3)
  if (flag == -1) {
    FFTW_API(execute_dft)(plan->plan_fast_forward,(FFT_DATA*)d_data.data(),(FFT_DATA*)d_data.data());
    FFTW_API(execute_dft)(plan->plan_slow_forward,(FFT_DATA*)d_data.data(),(FFT_DATA*)d_data.data());
  } else {
    FFTW_API(execute_dft)(plan->plan_fast_backward,(FFT_DATA*)d_data.data(),(FFT_DATA*)d_data.data());
    FFTW_API(execute_dft)(plan->plan_slow_backward,(FFT_DATA*)d_data.data(),(FFT_DATA*)d_data.data());
  }
#elif defined(FFT_CUFFT)
  cufftExec(plan->plan_fast,d_data.data(),d_data.data(),-flag);
  cufftExec(plan->plan_slow,d_data.data(),d_data.data(),-flag);
#elif defined(FFT_HIPFFT)
  hipfftExec(plan->plan_fast,d_data.data(),d_data.data(),-flag);
  hipfftExec(plan->plan_slow,d_data.data(),d_data.data(),-flag);
#else
  kiss_fft_functor<DeviceType> f;
    typename FFT_AT::t_FFT_DATA_1d d_tmp =
     typename FFT_AT::t_FFT_DATA_1d(Kokkos::view_alloc("fft_2d:tmp",Kokkos::WithoutInitializing),d_data.extent(0));
  if (flag == -1) {
    f = kiss_fft_functor<DeviceType>(d_data,d_tmp,plan->cfg_fast_forward,length1);
    Kokkos::parallel_for(total1/length1,f);

    f = kiss_fft_functor<DeviceType>(d_data,d_tmp,plan->cfg_slow_forward,length2);
    Kokkos::parallel_for(total2/length2,f);
  } else {
    f = kiss_fft_functor<DeviceType>(d_data,d_tmp,plan->cfg_fast_backward,length1);
    Kokkos::parallel_for(total1/length1,f);

    f = kiss_fft_functor<DeviceType>(d_data,d_tmp,plan->cfg_slow_backward,length2);
    Kokkos::parallel_for(total2/length2,f);
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
template class FFT2dKokkos<SPADeviceType>;
#ifdef SPARTA_KOKKOS_GPU
template class FFT2dKokkos<SPAHostType>;
#endif
}

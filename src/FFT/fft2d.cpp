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

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "fft2d.h"
#include "remap2d.h"

// include kissfft implementation

#ifdef FFT_KISSFFT
#include "kissfft.h"
#endif

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

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

void fft_2d(FFT_DATA *in, FFT_DATA *out, int flag, struct fft_plan_2d *plan)
{
  int i,total,length,offset,num;
  double norm,*out_ptr;
  FFT_DATA *data,*copy;

  // system specific constants

#if defined(FFT_SCSL)
  int isys = 0;
  FFT_PREC scalef = 1.0;
#elif defined(FFT_DEC)
  char c = 'C';
  char f = 'F';
  char b = 'B';
  int one = 1;
#elif defined(FFT_T3E)
  int isys = 0;
  double scalef = 1.0;
#elif defined(FFT_ACML)
  int info;
#elif defined(FFT_FFTW3)
  FFTW_API(plan) theplan;
#else
  // nothing to do for other FFTs.
#endif

  // pre-remap to prepare for 1st FFTs if needed
  // copy = loc for remap result

  if (plan->pre_plan) {
    if (plan->pre_target == 0) copy = out;
    else copy = plan->copy;
    remap_2d((FFT_SCALAR *) in, (FFT_SCALAR *) copy, 
             (FFT_SCALAR *) plan->scratch, plan->pre_plan);
    data = copy;
  } else data = in;

  // 1d FFTs along fast axis

  total = plan->total1;
  length = plan->length1;

#if defined(FFT_SGI)
  for (offset = 0; offset < total; offset += length)
    FFT_1D(flag,length,&data[offset],1,plan->coeff1);
#elif defined(FFT_SCSL)
  for (offset = 0; offset < total; offset += length)
    FFT_1D(flag,length,scalef,&data[offset],&data[offset],plan->coeff1,
           plan->work1,&isys);
#elif defined(FFT_ACML)
  num=total/length;
  FFT_1D(&flag,&num,&length,data,plan->coeff1,&info);
#elif defined(FFT_INTEL)
  for (offset = 0; offset < total; offset += length)
    FFT_1D(&data[offset],&length,&flag,plan->coeff1);
#elif defined(FFT_MKL)
  if (flag == -1)
    DftiComputeForward(plan->handle_fast,data);
  else
    DftiComputeBackward(plan->handle_fast,data);
#elif defined(FFT_DEC)
  if (flag == -1)
    for (offset = 0; offset < total; offset += length)
      FFT_1D(&c,&c,&f,&data[offset],&data[offset],&length,&one);
  else
    for (offset = 0; offset < total; offset += length)
      FFT_1D(&c,&c,&b,&data[offset],&data[offset],&length,&one);
#elif defined(FFT_T3E)
  for (offset = 0; offset < total; offset += length)
    FFT_1D(&flag,&length,&scalef,&data[offset],&data[offset],plan->coeff1,
           plan->work1,&isys);
#elif defined(FFT_FFTW2)
  if (flag == -1)
    fftw(plan->plan_fast_forward,total/length,data,1,length,NULL,0,0);
  else
    fftw(plan->plan_fast_backward,total/length,data,1,length,NULL,0,0);
#elif defined(FFT_FFTW3)
  if (flag == -1)
    theplan=plan->plan_fast_forward;
  else
    theplan=plan->plan_fast_backward;
  FFTW_API(execute_dft)(theplan,data,data);
#else
  if (flag == -1)
    for (offset = 0; offset < total; offset += length)
      kiss_fft(plan->cfg_fast_forward,&data[offset],&data[offset]);
  else
    for (offset = 0; offset < total; offset += length)
      kiss_fft(plan->cfg_fast_backward,&data[offset],&data[offset]);
#endif

  // mid-remap to prepare for 2nd FFTs
  // copy = loc for remap result

  if (plan->mid_target == 0) copy = out;
  else copy = plan->copy;
  remap_2d((FFT_SCALAR *) data, (FFT_SCALAR *) copy, 
           (FFT_SCALAR *) plan->scratch, plan->mid_plan);
  data = copy;

  // 1d FFTs along slow axis

  total = plan->total2;
  length = plan->length2;

#if defined(FFT_SGI)
  for (offset = 0; offset < total; offset += length)
    FFT_1D(flag,length,&data[offset],1,plan->coeff2);
#elif defined(FFT_SCSL)
  for (offset = 0; offset < total; offset += length)
    FFT_1D(flag,length,scalef,&data[offset],&data[offset],plan->coeff2,
           plan->work2,&isys);
#elif defined(FFT_ACML)
  num=total/length;
  FFT_1D(&flag,&num,&length,data,plan->coeff2,&info);
#elif defined(FFT_INTEL)
  for (offset = 0; offset < total; offset += length)
    FFT_1D(&data[offset],&length,&flag,plan->coeff2);
#elif defined(FFT_MKL)
  if (flag == -1)
    DftiComputeForward(plan->handle_mid,data);
  else
    DftiComputeBackward(plan->handle_mid,data);
#elif defined(FFT_DEC)
  if (flag == -1)
    for (offset = 0; offset < total; offset += length)
      FFT_1D(&c,&c,&f,&data[offset],&data[offset],&length,&one);
  else
    for (offset = 0; offset < total; offset += length)
      FFT_1D(&c,&c,&b,&data[offset],&data[offset],&length,&one);
#elif defined(FFT_T3E)
  for (offset = 0; offset < total; offset += length)
    FFT_1D(&flag,&length,&scalef,&data[offset],&data[offset],plan->coeff2,
           plan->work2,&isys);
#elif defined(FFT_FFTW2)
  if (flag == -1)
    fftw(plan->plan_mid_forward,total/length,data,1,length,NULL,0,0);
  else
    fftw(plan->plan_mid_backward,total/length,data,1,length,NULL,0,0);
#elif defined(FFT_FFTW3)
  if (flag == -1)
    theplan=plan->plan_mid_forward;
  else
    theplan=plan->plan_mid_backward;
  FFTW_API(execute_dft)(theplan,data,data);
#else
  if (flag == -1)
    for (offset = 0; offset < total; offset += length)
      kiss_fft(plan->cfg_slow_forward,&data[offset],&data[offset]);
  else
    for (offset = 0; offset < total; offset += length)
      kiss_fft(plan->cfg_slow_backward,&data[offset],&data[offset]);
#endif

  // post-remap to put data in output format if needed
  // destination is always out

  if (plan->post_plan)
    remap_2d((FFT_SCALAR *) data, (FFT_SCALAR *) out, 
             (FFT_SCALAR *) plan->scratch,plan->post_plan);

  // scaling if required

#if !defined(FFT_T3E) && !defined(FFT_ACML)
  if (flag == 1 && plan->scaled) {
    norm = plan->norm;
    num = plan->normnum;
    out_ptr = (FFT_SCALAR *)out;
    for (i = 0; i < num; i++) {
#if defined(FFT_FFTW3)
      *(out_ptr++) *= norm;
      *(out_ptr++) *= norm;
#elif defined(FFT_MKL)
      out[i] *= norm;
#else
      out[i].re *= norm;
      out[i].im *= norm;
#endif
    }
  }
#endif

#ifdef FFT_T3E
  if (flag == 1 && plan->scaled) {
    norm = plan->norm;
    num = plan->normnum;
    for (i = 0; i < num; i++) out[i] *= (norm,norm);
  }
#endif

#ifdef FFT_ACML
  norm = plan->norm;
  num = plan->normnum;
  for (i = 0; i < num; i++) {
    out[i].re *= norm;
    out[i].im *= norm;
  }
#endif
}

/* ----------------------------------------------------------------------
   Create plan for performing a 2d FFT

   Arguments:
   comm                 MPI communicator for the P procs which own the data
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
------------------------------------------------------------------------- */

struct fft_plan_2d *fft_2d_create_plan(
       MPI_Comm comm, int nfast, int nslow,
       int in_ilo, int in_ihi, int in_jlo, int in_jhi,
       int out_ilo, int out_ihi, int out_jlo, int out_jhi,
       int scaled, int permute, int *nbuf, int usecollective)
{
  struct fft_plan_2d *plan;
  int me,nprocs;
  int i,num,flag,remapflag,fftflag;
  int first_ilo,first_ihi,first_jlo,first_jhi;
  int second_ilo,second_ihi,second_jlo,second_jhi;
  int out_size,first_size,second_size,copy_size,scratch_size;
  int list[50];

  // system specific variables

#ifdef FFT_SCSL
  FFT_DATA dummy_d[5];
  FFT_PREC dummy_p[5];
  int isign,isys;
  FFT_PREC scalef;
#endif
#ifdef FFT_INTEL
  FFT_DATA dummy;
#endif
#ifdef FFT_T3E
  FFT_DATA dummy[5];
  int isign,isys;
  double scalef;
#endif

  // query MPI info

  MPI_Comm_rank(comm,&me);
  MPI_Comm_size(comm,&nprocs);

  // allocate memory for plan data struct

  plan = (struct fft_plan_2d *) malloc(sizeof(struct fft_plan_2d));
  if (plan == NULL) return NULL;

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
    plan->pre_plan = NULL;
  }
  else {
    first_ilo = 0;
    first_ihi = nfast - 1;
    first_jlo = me*nslow/nprocs;
    first_jhi = (me+1)*nslow/nprocs - 1;
    plan->pre_plan =
      remap_2d_create_plan(comm,in_ilo,in_ihi,in_jlo,in_jhi,
			   first_ilo,first_ihi,first_jlo,first_jhi,
			   FFT_PRECISION,0,0,2);
    if (plan->pre_plan == NULL) return NULL;
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
    remap_2d_create_plan(comm,first_ilo,first_ihi,first_jlo,first_jhi,
			 second_ilo,second_ihi,second_jlo,second_jhi,
			 FFT_PRECISION,1,0,2);
  if (plan->mid_plan == NULL) return NULL;

  // 1d FFTs along slow axis

  plan->length2 = nslow;
  plan->total2 = (second_ihi-second_ilo+1) * nslow;
  
  // remap from 2nd FFT to final distribution
  // not needed if permute = 1 and second indices = out indices on all procs

  if (permute == 1 &&
      out_ilo == second_ilo && out_ihi == second_ihi &&
      out_jlo == second_jlo && out_jhi == second_jhi) flag = 0;
  else flag = 1;

  MPI_Allreduce(&flag,&remapflag,1,MPI_INT,MPI_MAX,comm);

  if (remapflag == 0) {
    plan->post_plan = NULL;
  }
  else {
    plan->post_plan =
      remap_2d_create_plan(comm,second_jlo,second_jhi,second_ilo,second_ihi,
			   out_jlo,out_jhi,out_ilo,out_ihi,
			   FFT_PRECISION,(permute+1)%2,0,2);
    if (plan->post_plan == NULL) return NULL;
  }

  // configure plan memory pointers and allocate work space
  // out_size = amount of memory given to FFT by user
  // first/second_size = amount of memory needed after pre, mid remaps
  // copy_size = amount needed internally for extra copy of data
  // scratch_size = amount needed internally for remap scratch space
  // for each remap:
  //   use out space for result if big enough, else require copy buffer
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
    plan->copy = (FFT_DATA *) malloc(copy_size*sizeof(FFT_DATA));
    if (plan->copy == NULL) return NULL;
  }
  else plan->copy = NULL;

  if (scratch_size) {
    plan->scratch = (FFT_DATA *) malloc(scratch_size*sizeof(FFT_DATA));
    if (plan->scratch == NULL) return NULL;
  }
  else plan->scratch = NULL;

  // system specific pre-computation of 1d FFT coeffs 
  // and scaling normalization

#if defined(FFT_SGI)

  plan->coeff1 = (FFT_DATA *) malloc((nfast+15)*sizeof(FFT_DATA));
  plan->coeff2 = (FFT_DATA *) malloc((nslow+15)*sizeof(FFT_DATA));

  if (plan->coeff1 == NULL || plan->coeff2 == NULL) return NULL;

  FFT_1D_INIT(nfast,plan->coeff1);
  FFT_1D_INIT(nslow,plan->coeff2);

  if (scaled == 0)
    plan->scaled = 0;
  else {
    plan->scaled = 1;
    plan->norm = 1.0/(nfast*nslow);
    plan->normnum = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1);
  }

#elif defined(FFT_SCSL)

  plan->coeff1 = (FFT_PREC *) malloc((2*nfast+30)*sizeof(FFT_PREC));
  plan->coeff2 = (FFT_PREC *) malloc((2*slow+30)*sizeof(FFT_PREC));

  if (plan->coeff1 == NULL || plan->coeff2 == NULL) return NULL;

  plan->work1 = (FFT_PREC *) malloc((2*nfast)*sizeof(FFT_PREC));
  plan->work2 = (FFT_PREC *) malloc((2*nslow)*sizeof(FFT_PREC));

  if (plan->work1 == NULL || plan->work2 == NULL) return NULL;

  isign = 0;
  scalef = 1.0;
  isys = 0;

  FFT_1D_INIT(isign,nfast,scalef,dummy_d,dummy_d,plan->coeff1,dummy_p,&isys);
  FFT_1D_INIT(isign,nslow,scalef,dummy_d,dummy_d,plan->coeff2,dummy_p,&isys);

  if (scaled == 0)
    plan->scaled = 0;
  else {
    plan->scaled = 1;
    plan->norm = 1.0/(nfast*nslow);
    plan->normnum = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1);
  }

#elif defined(FFT_ACML)

  plan->coeff1 = (FFT_DATA *) malloc((3*nfast+100)*sizeof(FFT_DATA));
  plan->coeff2 = (FFT_DATA *) malloc((3*nslow+100)*sizeof(FFT_DATA));

  if (plan->coeff1 == NULL || plan->coeff2 == NULL) return NULL;

  int isign = 100;
  int isys = 1;
  int info = 0;
  FFT_DATA *dummy = NULL;

  FFT_1D(&isign,&isys,&nfast,dummy,plan->coeff1,&info);
  FFT_1D(&isign,&isys,&nslow,dummy,plan->coeff2,&info);

  if (scaled == 0) {
    plan->scaled = 0;
    plan->norm = sqrt(nfast*nslow);
    plan->normnum = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1);
  } else {
    plan->scaled = 1;
    plan->norm = sqrt(nfast*nslow*nslow);
    plan->normnum = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1);
  }

#elif defined(FFT_INTEL)

  flag = 0;

  num = 0;
  factor_2d(nfast,&num,list);
  for (i = 0; i < num; i++)
    if (list[i] != 2 && list[i] != 3 && list[i] != 5) flag = 1;
  num = 0;
  factor_2d(nslow,&num,list);
  for (i = 0; i < num; i++)
    if (list[i] != 2 && list[i] != 3 && list[i] != 5) flag = 1;

  MPI_Allreduce(&flag,&fftflag,1,MPI_INT,MPI_MAX,comm);
  if (fftflag) {
    if (me == 0) printf("ERROR: FFTs are not power of 2,3,5\n");
    return NULL;
  }

  plan->coeff1 = (FFT_DATA *) malloc((3*nfast/2+1)*sizeof(FFT_DATA));
  plan->coeff2 = (FFT_DATA *) malloc((3*nslow/2+1)*sizeof(FFT_DATA));

  if (plan->coeff1 == NULL || plan->coeff2 == NULL) return NULL;

  flag = 0;
  FFT_1D_INIT(&dummy,&nfast,&flag,plan->coeff1);
  FFT_1D_INIT(&dummy,&nslow,&flag,plan->coeff2);

  if (scaled == 0) {
    plan->scaled = 1;
    plan->norm = nfast*nslow;
    plan->normnum = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1);
  }
  else
    plan->scaled = 0;

#elif defined(FFT_MKL)
  DftiCreateDescriptor( &(plan->handle_fast), FFT_MKL_PREC, DFTI_COMPLEX, 1, (MKL_LONG)nfast);
  DftiSetValue(plan->handle_fast, DFTI_NUMBER_OF_TRANSFORMS, (MKL_LONG)plan->total1/nfast);
  DftiSetValue(plan->handle_fast, DFTI_PLACEMENT,DFTI_INPLACE);
  DftiSetValue(plan->handle_fast, DFTI_INPUT_DISTANCE, (MKL_LONG)nfast);
  DftiSetValue(plan->handle_fast, DFTI_OUTPUT_DISTANCE, (MKL_LONG)nfast);
  DftiCommitDescriptor(plan->handle_fast);

  DftiCreateDescriptor( &(plan->handle_slow), FFT_MKL_PREC, DFTI_COMPLEX, 1, (MKL_LONG)nslow);
  DftiSetValue(plan->handle_slow, DFTI_NUMBER_OF_TRANSFORMS, (MKL_LONG)plan->total2/nslow);
  DftiSetValue(plan->handle_slow, DFTI_PLACEMENT,DFTI_INPLACE);
  DftiSetValue(plan->handle_slow, DFTI_INPUT_DISTANCE, (MKL_LONG)nslow);
  DftiSetValue(plan->handle_slow, DFTI_OUTPUT_DISTANCE, (MKL_LONG)nslow);
  DftiCommitDescriptor(plan->handle_slow);

  if (scaled == 0)
    plan->scaled = 0;
  else {
    plan->scaled = 1;
    plan->norm = 1.0/(nfast*nslow);
    plan->normnum = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1);
  }

#elif defined(FFT_DEC)

  if (scaled == 0) {
    plan->scaled = 1;
    plan->norm = nfast*nslow;
    plan->normnum = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1);
  }
  else
    plan->scaled = 0;

#elif defined(FFT_T3E)

  plan->coeff1 = (double *) malloc((12*nfast)*sizeof(double));
  plan->coeff2 = (double *) malloc((12*nslow)*sizeof(double));

  if (plan->coeff1 == NULL || plan->coeff2 == NULL) return NULL;

  plan->work1 = (double *) malloc((8*nfast)*sizeof(double));
  plan->work2 = (double *) malloc((8*nslow)*sizeof(double));

  if (plan->work1 == NULL || plan->work2 == NULL) return NULL;

  isign = 0;
  scalef = 1.0;
  isys = 0;

  FFT_1D_INIT(&isign,&nfast,&scalef,dummy,dummy,plan->coeff1,dummy,&isys);
  FFT_1D_INIT(&isign,&nslow,&scalef,dummy,dummy,plan->coeff2,dummy,&isys);

  if (scaled == 0)
    plan->scaled = 0;
  else {
    plan->scaled = 1;
    plan->norm = 1.0/(nfast*nslow);
    plan->normnum = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1);
  }

#elif defined(FFT_FFTW2)

  plan->plan_fast_forward =
    fftw_create_plan(nfast,FFTW_FORWARD,FFTW_ESTIMATE | FFTW_IN_PLACE);
  plan->plan_fast_backward =
    fftw_create_plan(nfast,FFTW_BACKWARD,FFTW_ESTIMATE | FFTW_IN_PLACE);

  if (nslow == nfast) {
    plan->plan_slow_forward = plan->plan_fast_forward;
    plan->plan_slow_backward = plan->plan_fast_backward;
  }
  else {
    plan->plan_slow_forward =
      fftw_create_plan(nslow,FFTW_FORWARD,FFTW_ESTIMATE | FFTW_IN_PLACE);
    plan->plan_slow_backward =
      fftw_create_plan(nslow,FFTW_BACKWARD,FFTW_ESTIMATE | FFTW_IN_PLACE);
  }

  if (scaled == 0)
    plan->scaled = 0;
  else {
    plan->scaled = 1;
    plan->norm = 1.0/(nfast*nslow);
    plan->normnum = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1);
  }

#elif defined(FFT_FFTW3)
  plan->plan_fast_forward =
    FFTW_API(plan_many_dft)(1, &nfast,plan->total1/plan->length1,
                            NULL,&nfast,1,plan->length1,
                            NULL,&nfast,1,plan->length1,
                            FFTW_FORWARD,FFTW_ESTIMATE);
  plan->plan_fast_backward =
    FFTW_API(plan_many_dft)(1, &nfast,plan->total1/plan->length1,
                            NULL,&nfast,1,plan->length1,
                            NULL,&nfast,1,plan->length1,
                            FFTW_BACKWARD,FFTW_ESTIMATE);
  plan->plan_slow_forward =
    FFTW_API(plan_many_dft)(1, &nslow,plan->total2/plan->length2,
                            NULL,&nslow,1,plan->length2,
                            NULL,&nslow,1,plan->length2,
                            FFTW_FORWARD,FFTW_ESTIMATE);
  plan->plan_slow_backward =
    FFTW_API(plan_many_dft)(1, &nslow,plan->total2/plan->length2,
                            NULL,&nslow,1,plan->length2,
                            NULL,&nslow,1,plan->length2,
                            FFTW_BACKWARD,FFTW_ESTIMATE);

  if (scaled == 0)
    plan->scaled = 0;
  else {
    plan->scaled = 1;
    plan->norm = 1.0/(nfast*nslow);
    plan->normnum = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1);
  }

#else
  plan->cfg_fast_forward = kiss_fft_alloc(nfast,0,NULL,NULL);
  plan->cfg_fast_backward = kiss_fft_alloc(nfast,1,NULL,NULL);

  if (nslow == nfast) {
    plan->cfg_slow_forward = plan->cfg_fast_forward;
    plan->cfg_slow_backward = plan->cfg_fast_backward;
  } else {
    plan->cfg_slow_forward = kiss_fft_alloc(nslow,0,NULL,NULL);
    plan->cfg_slow_backward = kiss_fft_alloc(nslow,1,NULL,NULL);
  }

  if (scaled == 0)
    plan->scaled = 0;
  else {
    plan->scaled = 1;
    plan->norm = 1.0/(nfast*nslow);
    plan->normnum = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1);
  }

#endif

  return plan;
}

/* ----------------------------------------------------------------------
   Destroy a 2d fft plan
------------------------------------------------------------------------- */

void fft_2d_destroy_plan(struct fft_plan_2d *plan)
{
  if (plan->pre_plan) remap_2d_destroy_plan(plan->pre_plan);
  if (plan->mid_plan) remap_2d_destroy_plan(plan->mid_plan);
  if (plan->post_plan) remap_2d_destroy_plan(plan->post_plan);

  if (plan->copy) free(plan->copy);
  if (plan->scratch) free(plan->scratch);

#if defined(FFT_SGI)
  free(plan->coeff1);
  free(plan->coeff2);
#elif defined(FFT_SCSL)
  free(plan->coeff1);
  free(plan->coeff2);
  free(plan->work1);
  free(plan->work2);
#elif defined(FFT_ACML)
  free(plan->coeff1);
  free(plan->coeff2);
#elif defined(FFT_INTEL)
  free(plan->coeff1);
  free(plan->coeff2);
#elif defined(FFT_MKL)
  DftiFreeDescriptor(&(plan->handle_fast));
  DftiFreeDescriptor(&(plan->handle_slow));
#elif defined(FFT_T3E)
  free(plan->coeff1);
  free(plan->coeff2);
  free(plan->work1);
  free(plan->work2);
#elif defined(FFT_FFTW2)
  if (plan->plan_slow_forward != plan->plan_fast_forward) {
    fftw_destroy_plan(plan->plan_slow_forward);
    fftw_destroy_plan(plan->plan_slow_backward);
  }
  fftw_destroy_plan(plan->plan_fast_forward);
  fftw_destroy_plan(plan->plan_fast_backward);
#elif defined(FFT_FFTW3)
  FFTW_API(destroy_plan)(plan->plan_slow_forward);
  FFTW_API(destroy_plan)(plan->plan_slow_backward);
  FFTW_API(destroy_plan)(plan->plan_fast_forward);
  FFTW_API(destroy_plan)(plan->plan_fast_backward);
#else
  if (plan->cfg_slow_forward != plan->cfg_fast_forward) {
    free(plan->cfg_slow_forward);
    free(plan->cfg_slow_backward);
  }
  free(plan->cfg_fast_forward);
  free(plan->cfg_fast_backward);
#endif

  free(plan);
}

/* ----------------------------------------------------------------------
   recursively divide n into small factors, return them in list
------------------------------------------------------------------------- */

void factor_2d(int n, int *num, int *list)
{
  if (n == 1) {
    return;
  }
  else if (n % 2 == 0) {
    *list = 2;
    (*num)++;
    factor_2d(n/2,num,list+1);
  }
  else if (n % 3 == 0) {
    *list = 3;
    (*num)++;
    factor_2d(n/3,num,list+1);
  }
  else if (n % 5 == 0) {
    *list = 5;
    (*num)++;
    factor_2d(n/5,num,list+1);
  }
  else if (n % 7 == 0) {
    *list = 7;
    (*num)++;
    factor_2d(n/7,num,list+1);
  }
  else if (n % 11 == 0) {
    *list = 11;
    (*num)++;
    factor_2d(n/11,num,list+1);
  }
  else if (n % 13 == 0) {
    *list = 13;
    (*num)++;
    factor_2d(n/13,num,list+1);
  }
  else {
    *list = n;
    (*num)++;
    return;
  }
}

/* ----------------------------------------------------------------------
   perform just the 1d FFTs needed by a 2d FFT, no data movement
   used for timing purposes

   Arguments:
   in           starting address of input data on this proc, all set to 0.0
   nsize        size of in
   flag         1 for forward FFT, -1 for inverse FFT
   plan         plan returned by previous call to fft_2d_create_plan
------------------------------------------------------------------------- */

void fft_2d_1d_only(FFT_DATA *data, int nsize, int flag,
                    struct fft_plan_2d *plan)
{
  int i,total,length,offset,num;
  FFT_SCALAR norm,*data_ptr;

  // system specific constants

#ifdef FFT_SCSL
  int isys = 0;
  FFT_PREC scalef = 1.0;
#endif
#ifdef FFT_DEC
  char c = 'C';
  char f = 'F';
  char b = 'B';
  int one = 1;
#endif
#ifdef FFT_T3E
  int isys = 0;
  double scalef = 1.0;
#endif

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
  if ((total1 > nsize) || (total2 > nsize)) return;
#endif
  if (total1 > nsize) total1 = (nsize/length1) * length1;
  if (total2 > nsize) total2 = (nsize/length2) * length2;

  // perform 1d FFTs in each of 3 dimensions
  // data is just an array of 0.0

#ifdef FFT_SGI
  for (offset = 0; offset < total1; offset += length1)
    FFT_1D(flag,length1,&data[offset],1,plan->coeff1);
  for (offset = 0; offset < total2; offset += length2)
    FFT_1D(flag,length2,&data[offset],1,plan->coeff2);

#elif defined(FFT_SCSL)
  for (offset = 0; offset < total1; offset += length1)
    FFT_1D(flag,length1,scalef,&data[offset],&data[offset],plan->coeff1,
           plan->work1,&isys);
  for (offset = 0; offset < total2; offset += length2)
    FFT_1D(flag,length2,scalef,&data[offset],&data[offset],plan->coeff2,
           plan->work2,&isys);

#elif defined(FFT_ACML)
  int info=0;
  num=total1/length1;
  FFT_1D(&flag,&num,&length1,data,plan->coeff1,&info);
  num=total2/length2;
  FFT_1D(&flag,&num,&length2,data,plan->coeff2,&info);

#elif defined(FFT_INTEL)
  for (offset = 0; offset < total1; offset += length1)
    FFT_1D(&data[offset],&length1,&flag,plan->coeff1);
  for (offset = 0; offset < total2; offset += length2)
    FFT_1D(&data[offset],&length2,&flag,plan->coeff2);

#elif defined(FFT_MKL)
  if (flag == -1) {
    DftiComputeForward(plan->handle_fast,data);
    DftiComputeForward(plan->handle_slow,data);
  } else {
    DftiComputeBackward(plan->handle_fast,data);
    DftiComputeBackward(plan->handle_slow,data);
  }

#elif defined(FFT_DEC)
  if (flag == -1) {
    for (offset = 0; offset < total1; offset += length1)
      FFT_1D(&c,&c,&f,&data[offset],&data[offset],&length1,&one);
    for (offset = 0; offset < total2; offset += length2)
      FFT_1D(&c,&c,&f,&data[offset],&data[offset],&length2,&one);
  } else {
    for (offset = 0; offset < total1; offset += length1)
      FFT_1D(&c,&c,&b,&data[offset],&data[offset],&length1,&one);
    for (offset = 0; offset < total2; offset += length2)
      FFT_1D(&c,&c,&b,&data[offset],&data[offset],&length2,&one);
  }

#elif defined(FFT_T3E)
  for (offset = 0; offset < total1; offset += length1)
    FFT_1D(&flag,&length1,&scalef,&data[offset],&data[offset],plan->coeff1,
           plan->work1,&isys);
  for (offset = 0; offset < total2; offset += length2)
    FFT_1D(&flag,&length2,&scalef,&data[offset],&data[offset],plan->coeff2,
           plan->work2,&isys);

#elif defined(FFT_FFTW2)
  if (flag == -1) {
    fftw(plan->plan_fast_forward,total1/length1,data,1,0,NULL,0,0);
    fftw(plan->plan_slow_forward,total2/length2,data,1,0,NULL,0,0);
  } else {
    fftw(plan->plan_fast_backward,total1/length1,data,1,0,NULL,0,0);
    fftw(plan->plan_slow_backward,total2/length2,data,1,0,NULL,0,0);
  }

#elif defined(FFT_FFTW3)
  FFTW_API(plan) theplan;
  if (flag == -1)
    theplan=plan->plan_fast_forward;
  else
    theplan=plan->plan_fast_backward;
  FFTW_API(execute_dft)(theplan,data,data);
  if (flag == -1)
    theplan=plan->plan_slow_forward;
  else
    theplan=plan->plan_slow_backward;
  FFTW_API(execute_dft)(theplan,data,data);

#else
  if (flag == -1) {
    for (offset = 0; offset < total1; offset += length1)
      kiss_fft(plan->cfg_fast_forward,&data[offset],&data[offset]);
    for (offset = 0; offset < total2; offset += length2)
      kiss_fft(plan->cfg_slow_forward,&data[offset],&data[offset]);
  } else {
    for (offset = 0; offset < total1; offset += length1)
      kiss_fft(plan->cfg_fast_backward,&data[offset],&data[offset]);
    for (offset = 0; offset < total2; offset += length2)
      kiss_fft(plan->cfg_slow_backward,&data[offset],&data[offset]);
  }
#endif

  // scaling if required
  // limit num to size of data

#ifndef FFT_T3E
  if (flag == 1 && plan->scaled) {
    norm = plan->norm;
    num = MIN(plan->normnum,nsize);
    data_ptr = (FFT_SCALAR *)data;
    for (i = 0; i < num; i++) {
#if defined(FFT_FFTW3)
      *(data_ptr++) *= norm;
      *(data_ptr++) *= norm;
#elif defined(FFT_MKL)
      data[i] *= norm;
#else
      data[i].re *= norm;
      data[i].im *= norm;
#endif
    }
  }
#endif

#ifdef FFT_T3E
  if (flag == 1 && plan->scaled) {
    norm = plan->norm;
    num = MIN(plan->normnum,nsize);
    for (i = 0; i < num; i++) data[i] *= (norm,norm);
  }
#endif
}

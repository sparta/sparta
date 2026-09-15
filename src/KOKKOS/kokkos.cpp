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

#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "ctype.h"
#include "kokkos.h"
#include "sparta.h"
#include "error.h"
#include "memory_kokkos.h"

using namespace SPARTA_NS;

// Kokkos may be initialized at most once per process and, once finalized, can
// never be initialized again. When SPARTA is embedded as a library (a host
// application reusing one process across many open/close cycles) Kokkos is
// initialized on first use and finalized by KokkosSPARTA::finalize(), not by
// the KokkosSPARTA destructor.  main() calls it via sparta_kokkos_finalize();
// an embedder calls that same library function when it is done with SPARTA.

int KokkosSPARTA::is_finalized = 0;

static int kokkos_initialized_nthreads = 0;

/* ---------------------------------------------------------------------- */

KokkosSPARTA::KokkosSPARTA(SPARTA *sparta, int narg, char **arg) : Pointers(sparta)
{
  kokkos_exists = 1;
  sparta->kokkos = this;

  delete memory;
  memory = new MemoryKokkos(sparta);
  memoryKK = (MemoryKokkos*) memory;

  int me = 0;
  MPI_Comm_rank(world,&me);
  if (me == 0) error->message(FLERR,"KOKKOS mode is enabled");

  // process any command-line args that invoke Kokkos settings

  ngpus = 0;
  int device = 0;
  nthreads = 1;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"d") == 0 || strcmp(arg[iarg],"device") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Invalid Kokkos command-line args");
      device = atoi(arg[iarg+1]);
      iarg += 2;

    } else if (strcmp(arg[iarg],"g") == 0 ||
               strcmp(arg[iarg],"gpus") == 0) {
#ifndef SPARTA_KOKKOS_GPU
      error->all(FLERR,"GPUs are requested but Kokkos has not been compiled with a GPU-enabled backend");
#endif
      if (iarg+2 > narg) error->all(FLERR,"Invalid Kokkos command-line args");
      ngpus = atoi(arg[iarg+1]);

      int skip_gpu = 9999;
      if (iarg+2 < narg && isdigit(arg[iarg+2][0])) {
        skip_gpu = atoi(arg[iarg+2]);
        iarg++;
      }
      iarg += 2;

      int set_flag = 0;
      char *str;
      if ((str = getenv("SLURM_LOCALID"))) {
        int local_rank = atoi(str);
        device = local_rank % ngpus;
        if (device >= skip_gpu) device++;
        set_flag = 1;
      }
      if ((str = getenv("FLUX_TASK_LOCAL_ID"))) {
        int local_rank = atoi(str);
        device = local_rank % ngpus;
        if (device >= skip_gpu) device++;
        set_flag = 1;
      }
      if ((str = getenv("MPT_LRANK"))) {
        int local_rank = atoi(str);
        device = local_rank % ngpus;
        if (device >= skip_gpu) device++;
        set_flag = 1;
      }
      if ((str = getenv("MV2_COMM_WORLD_LOCAL_RANK"))) {
        int local_rank = atoi(str);
        device = local_rank % ngpus;
        if (device >= skip_gpu) device++;
        set_flag = 1;
      }
      if ((str = getenv("OMPI_COMM_WORLD_LOCAL_RANK"))) {
        int local_rank = atoi(str);
        device = local_rank % ngpus;
        if (device >= skip_gpu) device++;
        set_flag = 1;
      }
      if ((str = getenv("PMI_LOCAL_RANK"))) {
        int local_rank = atoi(str);
        device = local_rank % ngpus;
        if (device >= skip_gpu) device++;
        set_flag = 1;
      }

      if (ngpus > 1 && !set_flag)
        error->all(FLERR,"Could not determine local MPI rank for multiple "
                           "GPUs with because MPI library not recognized");

    } else if (strcmp(arg[iarg],"t") == 0 ||
               strcmp(arg[iarg],"threads") == 0) {
      nthreads = atoi(arg[iarg+1]);
      iarg += 2;

    } else error->all(FLERR,"Invalid Kokkos command-line args");
  }

  // initialize Kokkos

  if (me == 0) {
    if (screen) fprintf(screen,"  requested %d GPU(s) per node\n",ngpus);
    if (logfile) fprintf(logfile,"  requested %d GPU(s) per node\n",ngpus);

    if (screen) fprintf(screen,"  requested %d thread(s) per MPI task\n",nthreads);
    if (logfile) fprintf(logfile,"  requested %d thread(s) per MPI task\n",nthreads);
  }

#ifdef SPARTA_KOKKOS_GPU
  if (ngpus <= 0)
    error->all(FLERR,"Kokkos has been compiled with a GPU-enabled backend but no GPUs are requested");
#endif

#ifndef KOKKOS_ENABLE_SERIAL
  if (nthreads == 1 && me == 0)
    error->warning(FLERR,"When using a single thread, the Kokkos Serial backend "
                         "(i.e. Makefile.kokkos_mpi_only) gives better performance "
                         "than the OpenMP backend");
#endif

  Kokkos::InitializationSettings args;
  args.set_num_threads(nthreads);
  args.set_device_id(device);

  // Initialize Kokkos only once per process (it can be initialized at most
  // once).  On any later re-open the requested thread count cannot be changed,
  // so keep the count Kokkos was actually initialized with -- otherwise the
  // atomics decision below would be made for the wrong number of threads.
  if (!Kokkos::is_initialized()) {
    if (is_finalized)
      error->all(FLERR,"Kokkos package already finalized, cannot re-initialize");
    Kokkos::initialize(args);
    kokkos_initialized_nthreads = nthreads;
  } else {
    if (nthreads != kokkos_initialized_nthreads && me == 0)
      error->warning(FLERR,"Kokkos is already initialized in this process; "
                     "ignoring the new thread count. Restart to change the "
                     "number of threads.");
    nthreads = kokkos_initialized_nthreads;
  }

  // default settings for package kokkos command

  prewrap = 1;
  auto_sync = 1;
  gpu_aware_flag = 1;

  if (ngpus > 0) {
    comm_serial = 0;

    // SPARTA_KOKKOS_REDUCE_ARCH (kokkos_type.h) is the single definition of
    //  which architectures UpdateKokkos::move() dispatches the
    //  ATOMIC_REDUCTION = -1 (parallel_reduce) kernel for; the counters are
    //  read back from the reduction result only when atomic_reduction is 0

#if SPARTA_KOKKOS_REDUCE_ARCH
    atomic_reduction = 0;
#else
    atomic_reduction = 1;
#endif
  } else {

    // on CPU the host migrate path beats the device pack/unpack kernels, so
    //   it stays the default.  Measured on 4 MPI ranks, Serial backend, 400k
    //   particles in a 20^3 grid over 400 steps with 2.65% of particles
    //   migrating per step, comparing the Comm section of the timing
    //   breakdown (median of 3):
    //
    //                      comm serial   comm threaded
    //     free molecular      0.228 s       0.246 s   (+7.6%)
    //     with VSS collide    0.211 s       0.259 s   (+22.5%)
    //
    //   There is no host/device transfer to avoid here, so the device path
    //   only adds the irregular-comm plan rebuild.  Users can still ask for
    //   it with "package kokkos comm threaded".

    comm_serial = 1;
    atomic_reduction = 0;
  }

  need_atomics = 1;
  if (nthreads == 1 && ngpus == 0)
    need_atomics = 0;

  react_retry_flag = 0;
  react_extra = 1.1;
}

/* ---------------------------------------------------------------------- */

KokkosSPARTA::~KokkosSPARTA()
{
  // Kokkos is finalized by KokkosSPARTA::finalize(), not here, so a library
  // embedder can destroy and re-create SPARTA in the same process without
  // tripping over Kokkos's initialize-at-most-once restriction.
}

/* ----------------------------------------------------------------------
   shut down Kokkos
   called from main() before it returns, and from sparta_kokkos_finalize()
   Kokkos has to be finalized while its own state is still intact, so this
     must not be deferred to a static destructor or an atexit handler
------------------------------------------------------------------------- */

void KokkosSPARTA::finalize()
{
  if (Kokkos::is_initialized() && !is_finalized)
    Kokkos::finalize();
  is_finalized = 1;
}

/* ----------------------------------------------------------------------
   invoked by package kokkos command
------------------------------------------------------------------------- */

void KokkosSPARTA::accelerator(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"comm") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package kokkos command");
      if (strcmp(arg[iarg+1],"serial") == 0) {
        comm_serial = 1;
      } else if (strcmp(arg[iarg+1],"classic") == 0) { // deprecated
        comm_serial = 1;
      } else if (strcmp(arg[iarg+1],"threaded") == 0) {
        comm_serial = 0;
      } else error->all(FLERR,"Illegal package kokkos command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"react/retry") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package kokkos command");
      if (strcmp(arg[iarg+1],"yes") == 0) {
        react_retry_flag = 1;
      } else if (strcmp(arg[iarg+1],"no") == 0) {
        react_retry_flag = 0;
      } else error->all(FLERR,"Illegal package kokkos command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"react/extra") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package kokkos command");
      react_extra = atof(arg[iarg+1]);
      iarg += 2;
    } else if ((strcmp(arg[iarg],"gpu/aware") == 0)
               || (strcmp(arg[iarg],"gpu/direct") == 0)) { // gpu/direct is deprecated
      if (iarg+2 > narg) error->all(FLERR,"Illegal package kokkos command");
      if (strcmp(arg[iarg+1],"yes") == 0) {
        gpu_aware_flag = 1;
      } else if (strcmp(arg[iarg+1],"no") == 0) {
        gpu_aware_flag = 0;
      } else error->all(FLERR,"Illegal package kokkos command");
      iarg += 2;
    } else error->all(FLERR,"Illegal package kokkos command");
  }
}

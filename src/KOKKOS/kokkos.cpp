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

#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "ctype.h"
#include "kokkos.h"
#include "sparta.h"
#include "error.h"
#include "memory_kokkos.h"

using namespace SPARTA_NS;

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
  numa = 1;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"d") == 0 || strcmp(arg[iarg],"device") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Invalid Kokkos command-line args");
      device = atoi(arg[iarg+1]);
      iarg += 2;

    } else if (strcmp(arg[iarg],"g") == 0 ||
               strcmp(arg[iarg],"gpus") == 0) {
#ifndef KOKKOS_HAVE_CUDA
      error->all(FLERR,"GPUs are requested but Kokkos has not been compiled for CUDA");
#endif
      if (iarg+2 > narg) error->all(FLERR,"Invalid Kokkos command-line args");
      ngpus = atoi(arg[iarg+1]);

      int skip_gpu = 9999;
      if (iarg+2 < narg && isdigit(arg[iarg+2][0])) {
        skip_gpu = atoi(arg[iarg+2]);
        iarg++;
      }
      iarg += 2;

      char *str;
      if ((str = getenv("SLURM_LOCALID"))) {
        int local_rank = atoi(str);
        device = local_rank % ngpus;
        if (device >= skip_gpu) device++;
      }
      if ((str = getenv("MV2_COMM_WORLD_LOCAL_RANK"))) {
        int local_rank = atoi(str);
        device = local_rank % ngpus;
        if (device >= skip_gpu) device++;
      }
      if ((str = getenv("OMPI_COMM_WORLD_LOCAL_RANK"))) {
        int local_rank = atoi(str);
        device = local_rank % ngpus;
        if (device >= skip_gpu) device++;
      }

    } else if (strcmp(arg[iarg],"t") == 0 ||
               strcmp(arg[iarg],"threads") == 0) {
      nthreads = atoi(arg[iarg+1]);
      iarg += 2;

    } else if (strcmp(arg[iarg],"n") == 0 ||
               strcmp(arg[iarg],"numa") == 0) {
      numa = atoi(arg[iarg+1]);
      iarg += 2;

    } else error->all(FLERR,"Invalid Kokkos command-line args");
  }

  // initialize Kokkos

  if (me == 0) {
    if (screen) fprintf(screen,"  using %d GPU(s) per MPI task\n",ngpus);
    if (logfile) fprintf(logfile,"  using %d GPU(s) per MPI task\n",ngpus);

    if (screen) fprintf(screen,"  using %d thread(s) per MPI task\n",nthreads);
    if (logfile) fprintf(logfile,"  using %d thread(s) per MPI task\n",nthreads);
  }

#ifdef KOKKOS_HAVE_CUDA
  if (ngpus <= 0)
    error->all(FLERR,"Kokkos has been compiled for CUDA but no GPUs are requested");
#endif

  Kokkos::InitArguments args;
  args.num_threads = nthreads;
  args.num_numa = numa;
  args.device_id = device;

  Kokkos::initialize(args);

  // default settings for package kokkos command

  comm_classic = 0;
  atomic_reduction = 0;
  prewrap = 1;
  auto_sync = 0;

  need_atomics = 1;
  if (nthreads == 1 && ngpus == 0)
    need_atomics = 0;

  collide_retry_flag = 0;
  collide_extra = 1.1;

  //if (need_atomics == 0) // prevent unnecessary parallel_reduce
  //  atomic_reduction = 1;
}

/* ---------------------------------------------------------------------- */

KokkosSPARTA::~KokkosSPARTA()
{
  // finalize Kokkos

  Kokkos::finalize();
}

/* ----------------------------------------------------------------------
   invoked by package kokkos command
------------------------------------------------------------------------- */

void KokkosSPARTA::accelerator(int narg, char **arg)
{
  // defaults

  comm_classic = 0;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"comm") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package kokkos command");
      if (strcmp(arg[iarg+1],"classic") == 0) {
        comm_classic = 1;
      } else if (strcmp(arg[iarg+1],"threaded") == 0) {
        comm_classic = 0;
      } else error->all(FLERR,"Illegal package kokkos command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"reduction") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package kokkos command");
      if (strcmp(arg[iarg+1],"atomic") == 0) {
        atomic_reduction = 1;
      } else if (strcmp(arg[iarg+1],"parallel/reduce") == 0) {
        atomic_reduction = 0;
      } else error->all(FLERR,"Illegal package kokkos command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"collide/retry") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package kokkos command");
      if (strcmp(arg[iarg+1],"yes") == 0) {
        collide_retry_flag = 1;
      } else if (strcmp(arg[iarg+1],"no") == 0) {
        collide_retry_flag = 0;
      } else error->all(FLERR,"Illegal package kokkos command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"collide/extra") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal package kokkos command");
        collide_extra = atof(arg[iarg+1]);
      iarg += 1;
    } else error->all(FLERR,"Illegal package kokkos command");
  }
}

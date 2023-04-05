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

#include <mpi.h>

#ifdef FFT_SINGLE
typedef float FFT_SCALAR;
#define MPI_FFT_SCALAR MPI_FLOAT
#else
typedef double FFT_SCALAR;
#define MPI_FFT_SCALAR MPI_DOUBLE
#endif

// details of how to do a 2d remap

struct remap_plan_2d {
  double *sendbuf;                  // buffer for MPI sends
  double *scratch;                  // scratch buffer for MPI recvs
  void (*pack)(FFT_SCALAR *, FFT_SCALAR *, struct pack_plan_2d *);
                                    // which pack function to use
  void (*unpack)(FFT_SCALAR *, FFT_SCALAR *, struct pack_plan_2d *);
                                    // which unpack function to use
  int *send_offset;                 // extraction loc for each send
  int *send_size;                   // size of each send message
  int *send_proc;                   // proc to send each message to
  struct pack_plan_2d *packplan;    // pack plan for each send message
  int *recv_offset;                 // insertion loc for each recv
  int *recv_size;                   // size of each recv message
  int *recv_proc;                   // proc to recv each message from
  int *recv_bufloc;                 // offset in scratch buf for each recv
  MPI_Request *request;             // MPI request for each posted recv
  struct pack_plan_2d *unpackplan;  // unpack plan for each recv message
  int nrecv;                        // # of recvs from other procs
  int nsend;                        // # of sends to other procs
  int self;                         // whether I send/recv with myself
  int memory;                       // user provides scratch space or not
  MPI_Comm comm;                    // group of procs performing remap
};

// collision between 2 regions

struct extent_2d {
  int ilo,ihi,isize;
  int jlo,jhi,jsize;
};

// function prototypes

void remap_2d(FFT_SCALAR *, FFT_SCALAR *, FFT_SCALAR *, struct remap_plan_2d *);
struct remap_plan_2d *remap_2d_create_plan(MPI_Comm,
                                           int, int, int, int,
                                           int, int, int, int,
                                           int, int, int, int);
void remap_2d_destroy_plan(struct remap_plan_2d *);
int remap_2d_collide(struct extent_2d *,
                     struct extent_2d *, struct extent_2d *);

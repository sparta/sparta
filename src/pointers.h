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

// Pointers class contains ptrs to master copy of
//   fundamental SPARTA class ptrs stored in sparta.h
// every SPARTA class inherits from Pointers to access sparta.h ptrs
// these variables are auto-initialized by Pointer class constructor
// *& variables are really pointers to the pointers in sparta.h
// & enables them to be accessed directly in any class, e.g. atom->x

#ifndef SPARTA_POINTERS_H
#define SPARTA_POINTERS_H

#include "spatype.h"
#include "mpi.h"
#include "stdio.h"
#include "sparta.h"

namespace SPARTA_NS {

// universal defines inside namespace

#define FLERR __FILE__,__LINE__

#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))

// roundup a char ptr to 8-byte boundary
// roundup an int to multiple of 8

#define ROUNDUP(A) (char *) (((uint64_t) (A) + 7) & ~7)
#define IROUNDUP(A) ((((int) (A) + 7) / 8) * 8)
#define BIROUNDUP(A) ((((bigint) (A) + 7) / 8) * 8)

class Pointers {
 public:
  Pointers(SPARTA *ptr) :
    sparta(ptr),
    memory(ptr->memory),
    error(ptr->error),
    universe(ptr->universe),
    input(ptr->input),
    particle(ptr->particle),
    update(ptr->update),
    comm(ptr->comm),
    domain(ptr->domain),
    modify(ptr->modify),
    grid(ptr->grid),
    surf(ptr->surf),
    collide(ptr->collide),
    react(ptr->react),
    output(ptr->output),
    timer(ptr->timer),
    memoryKK(ptr->memoryKK),
    world(ptr->world),
    infile(ptr->infile),
    screen(ptr->screen),
    logfile(ptr->logfile) {}

  virtual ~Pointers() {}

 protected:
  SPARTA *sparta;
  Memory *&memory;
  Error *&error;
  Universe *&universe;
  Input *&input;

  Particle *&particle;
  Update *&update;
  Comm *&comm;
  Domain *&domain;
  Modify *&modify;
  Grid *&grid;
  Surf *&surf;
  Collide *&collide;
  React *&react;
  Output *&output;
  Timer *&timer;

  MemoryKokkos *&memoryKK;

  MPI_Comm &world;
  FILE *&infile;
  FILE *&screen;
  FILE *&logfile;
};

}

#endif

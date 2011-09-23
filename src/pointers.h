/* ----------------------------------------------------------------------
   DSMC - Sandia parallel DSMC code
   www.sandia.gov/~sjplimp/dsmc.html
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2011) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level DSMC directory.
------------------------------------------------------------------------- */

// Pointers class contains ptrs to master copy of
//   fundamental DSMC class ptrs stored in dsmc.h
// every DSMC class inherits from Pointers to access dsmc.h ptrs
// these variables are auto-initialized by Pointer class constructor
// *& variables are really pointers to the pointers in dsmc.h
// & enables them to be accessed directly in any class, e.g. atom->x

#ifndef DSMC_POINTERS_H
#define DSMC_POINTERS_H

#include "dsmctype.h"
#include "mpi.h"
#include "dsmc.h"

namespace DSMC_NS {

// universal defines inside namespace

#define FLERR __FILE__,__LINE__

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

class Pointers {
 public:
  Pointers(DSMC *ptr) : 
    dsmc(ptr),
    memory(ptr->memory),
    error(ptr->error),
    universe(ptr->universe),
    input(ptr->input),
    particle(ptr->particle),
    update(ptr->update),
    comm(ptr->comm),
    domain(ptr->domain),
    grid(ptr->grid),
    surf(ptr->surf),
    output(ptr->output),
    timer(ptr->timer),
    world(ptr->world),
    infile(ptr->infile),
    screen(ptr->screen),
    logfile(ptr->logfile) {}
  virtual ~Pointers() {}

 protected:
  DSMC *dsmc;
  Memory *&memory;
  Error *&error;
  Universe *&universe;
  Input *&input;

  Particle *&particle;
  Update *&update;
  Comm *&comm;
  Domain *&domain;
  Grid *&grid;
  Surf *&surf;
  Output *&output;
  Timer *&timer;

  MPI_Comm &world;
  FILE *&infile;
  FILE *&screen;
  FILE *&logfile;
};

}

#endif

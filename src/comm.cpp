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

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "comm.h"
#include "irregular.h"
#include "particle.h"
#include "grid.h"
#include "update.h"
#include "memory.h"
#include "error.h"

using namespace DSMC_NS;

/* ---------------------------------------------------------------------- */

Comm::Comm(DSMC *dsmc) : Pointers(dsmc)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  ncomm = 0;

  irregular = new Irregular(dsmc);

  maxsend = 0;
  proclist = NULL;
  sbuf = rbuf = NULL;
}

/* ---------------------------------------------------------------------- */

Comm::~Comm()
{
  delete irregular;
  memory->destroy(proclist);
  memory->destroy(sbuf);
  memory->destroy(rbuf);
}

/* ----------------------------------------------------------------------
   migrate particles to new procs after move
------------------------------------------------------------------------- */

void Comm::migrate(int nmigrate, int *mlist)
{
  int i,j;

  Grid::OneCell *cells = grid->cells;
  Particle::OnePart *particles = particle->particles;
  int nbytes = sizeof(Particle::OnePart);

  // grow proclist and sbuf if necessary

  if (nmigrate > maxsend) {
    maxsend = nmigrate;
    memory->destroy(proclist);
    memory->create(proclist,maxsend,"comm:proclist");
    memory->destroy(sbuf);
    memory->create(sbuf,maxsend*nbytes,"comm:sbuf");
  }

  // fill proclist with procs to send to
  // pack sbuf with particles to migrate
  // if icell < 0, particle is deleted but not sent
  // nsend = particles that actually migrate

  int nsend = 0;
  int offset = 0;
  for (i = 0; i < nmigrate; i++) {
    j = mlist[i];
    if (particles[j].icell < 0) continue;
    proclist[nsend++] = cells[particles[j].icell].proc;
    memcpy(&sbuf[offset],&particles[j],nbytes);
    offset += nbytes;
  }

  // compress my list of particles

  particle->compress(nmigrate,mlist);

  // create irregular communication plan, perform comm, destroy plan
  // returned nrecv = size of buffer needed for incoming atoms

  int nrecv = irregular->create(nsend,proclist);

  // extend particle list if necessary

  particle->grow(nrecv);

  // perform irregular communication

  irregular->exchange(sbuf,nbytes,
		      (char *) &particle->particles[particle->nlocal]);
  particle->nlocal += nrecv;

  ncomm += nsend;
}

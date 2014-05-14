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

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "comm.h"
#include "irregular.h"
#include "particle.h"
#include "grid.h"
#include "update.h"
#include "output.h"
#include "dump.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{PKEEP,PINSERT,PDONE,PDISCARD,PENTRY,PEXIT};   // several files

/* ---------------------------------------------------------------------- */

Comm::Comm(SPARTA *sparta) : Pointers(sparta)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  ncomm = 0;
  commsortflag = 0;
  commpartstyle = 1;

  neighflag = 0;
  neighlist = NULL;

  irregular = new Irregular(sparta);
  irregular_grid = NULL;

  pproc = NULL;
  maxpproc = 0;
  gproc = gsize = NULL;
  maxgproc = 0;
  sbuf = rbuf = NULL;
  maxsendbuf = maxrecvbuf = 0;
}

/* ---------------------------------------------------------------------- */

Comm::~Comm()
{
  delete irregular;
  delete irregular_grid;
  memory->destroy(pproc);
  memory->destroy(gproc);
  memory->destroy(gsize);
  memory->destroy(sbuf);
  memory->destroy(rbuf);

  memory->destroy(neighlist);
}

/* ----------------------------------------------------------------------
   reset neighbor list used in particle comm and setup irregular for them
   invoked after grid decomposition changes
   no-op if commpartstyle not set or grid decomposition not clumped
     since different mode of irregular comm will be done
------------------------------------------------------------------------- */

void Comm::reset_neighbors()
{
  neighflag = 0;
  if (!commpartstyle || !grid->clumped) return;
  neighflag = 1;

  if (neighlist == NULL)
    memory->create(neighlist,nprocs,"comm:neighlist");
  for (int i = 0; i < nprocs; i++) neighlist[i] = 0;

  Grid::ChildCell *cells = grid->cells;
  int nglocal = grid->nlocal;
  int ntotal = nglocal + grid->nghost;
  
  for (int icell = nglocal; icell < ntotal; icell++)
    neighlist[cells[icell].proc] = 1;
  neighlist[me] = 0;
  
  nneigh = 0;
  for (int i = 0; i < nprocs; i++)
    if (neighlist[i]) neighlist[nneigh++] = i;
  
  irregular->create_procs(nneigh,neighlist,commsortflag);
}

/* ----------------------------------------------------------------------
   migrate particles to new procs after particle move
   return particle nlocal after compression, 
     so Update can iterate on particle move
------------------------------------------------------------------------- */

int Comm::migrate_particles(int nmigrate, int *plist)
{
  int i,j;

  Grid::ChildCell *cells = grid->cells;
  Particle::OnePart *particles = particle->particles;
  int nbytes = sizeof(Particle::OnePart);

  // grow pproc and sbuf if necessary

  if (nmigrate > maxpproc) {
    maxpproc = nmigrate;
    memory->destroy(pproc);
    memory->create(pproc,maxpproc,"comm:pproc");
  }
  if (nmigrate*nbytes > maxsendbuf) {
    maxsendbuf = nmigrate*nbytes;
    memory->destroy(sbuf);
    memory->create(sbuf,maxsendbuf,"comm:sbuf");
  }

  // fill proclist with procs to send to
  // pack sbuf with particles to migrate
  // if flag == PDISCARD, particle is deleted but not sent
  // change icell of migrated particle to owning cell on receiving proc
  // nsend = particles that actually migrate

  int nsend = 0;
  int offset = 0;
  for (i = 0; i < nmigrate; i++) {
    j = plist[i];
    if (particles[j].flag == PDISCARD) continue;
    pproc[nsend++] = cells[particles[j].icell].proc;
    particles[j].icell = cells[particles[j].icell].ilocal;
    memcpy(&sbuf[offset],&particles[j],nbytes);
    offset += nbytes;
  }

  // compress my list of particles

  particle->compress(nmigrate,plist);
  int ncompress = particle->nlocal;

  // create or augment irregular communication plan
  // nrecv = # of incoming particles
  
  int nrecv;
  if (neighflag) nrecv = irregular->augment_data_uniform(nsend,pproc);
  else nrecv = irregular->create_data_uniform(nsend,pproc,commsortflag);

  // extend particle list if necessary

  particle->grow(nrecv);

  // perform irregular communication
  // append received particles directly to particle list

  irregular->exchange_uniform(sbuf,nbytes,
                              (char *) &particle->particles[particle->nlocal]);
  particle->nlocal += nrecv;

  ncomm += nsend;
  return ncompress;
}

/* ----------------------------------------------------------------------
   migrate grid cells with their particles to new procs
   called from BalanceGrid and FixBalance
------------------------------------------------------------------------- */

void Comm::migrate_cells(int nmigrate)
{
  int i,n;

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  int nglocal = grid->nlocal;

  // grow proc and size lists if needed

  if (nmigrate > maxgproc) {
    maxgproc = nmigrate;
    memory->destroy(gproc);
    memory->destroy(gsize);
    memory->create(gproc,maxgproc,"comm:gproc");
    memory->create(gsize,maxgproc,"comm:gsize");
  }

  // fill proclist with procs to send to
  // compute byte count needed to pack cells

  int nsend = 0;
  bigint boffset = 0;
  for (int icell = 0; icell < nglocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    if (cells[icell].proc == me) continue;
    gproc[nsend] = cells[icell].proc;
    n = grid->pack_one(icell,NULL,1,1,0);
    gsize[nsend++] = n;
    boffset += n;
  }

  if (boffset > MAXSMALLINT) 
    error->one(FLERR,"Migrate cells send buffer exceeds 2 GB");
  int offset = boffset;

  // reallocate sbuf as needed

  if (offset > maxsendbuf) {
    memory->destroy(sbuf);
    maxsendbuf = offset;
    memory->create(sbuf,maxsendbuf,"comm:sbuf");
    memset(sbuf,0,maxsendbuf);
  }

  // pack cell info into sbuf
  // only called for unsplit and split cells I no longer own

  offset = 0;
  for (int icell = 0; icell < nglocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    if (cells[icell].proc == me) continue;
    offset += grid->pack_one(icell,&sbuf[offset],1,1,1);
  }

  // compress my list of owned grid cells to remove migrated cells
  // compress particle list to remove particles in migrating cells

  if (nmigrate) {
    grid->compress();
    particle->compress();
  }

  // create irregular communication plan with variable size datums
  // nrecv = # of incoming grid cells
  // recvsize = total byte size of incoming grid + particle info
  // DEBUG: append a sort=1 arg so that messages from other procs
  //        are received in repeatable order, thus grid cells stay in order

  if (!irregular_grid) irregular_grid = new Irregular(sparta);
  int recvsize;
  int nrecv = 
    irregular_grid->create_data_variable(nmigrate,gproc,gsize,
                                         recvsize,commsortflag);

  // reallocate rbuf as needed

  if (recvsize > maxrecvbuf) {
    memory->destroy(rbuf);
    maxrecvbuf = recvsize;
    memory->create(rbuf,maxrecvbuf,"comm:rbuf");
    memset(rbuf,0,maxrecvbuf);
  }

  // perform irregular communication

  irregular_grid->exchange_variable(sbuf,gsize,rbuf);

  // unpack received grid cells with their particles

  offset = 0;
  for (i = 0; i < nrecv; i++)
    offset += grid->unpack_one(&rbuf[offset],1,1);

  // signal any dump grid classes to resize their per-cell arrays

  for (i = 0; i < output->ndump; i++)
    output->dump[i]->reset_grid();
}

/* ----------------------------------------------------------------------
   communicate inbuf around full ring of processors with messtag
   nbytes = size of inbuf = n datums * nper bytes
   callback() is invoked to allow caller to process/update each proc's inbuf
   note that callback() is invoked on final iteration for original inbuf
   for non-NULL outbuf, final updated inbuf is copied to it
   outbuf = inbuf is OK
------------------------------------------------------------------------- */

void Comm::ring(int n, int nper, void *inbuf, int messtag,
                void (*callback)(int, char *), void *outbuf, int self)
{
  MPI_Request request;
  MPI_Status status;

  int nbytes = n*nper;
  int maxbytes;
  MPI_Allreduce(&nbytes,&maxbytes,1,MPI_INT,MPI_MAX,world);

  char *buf,*bufcopy;
  memory->create(buf,maxbytes,"comm:buf");
  memory->create(bufcopy,maxbytes,"comm:bufcopy");
  memcpy(buf,inbuf,nbytes);

  int next = me + 1;
  int prev = me - 1;
  if (next == nprocs) next = 0;
  if (prev < 0) prev = nprocs - 1;

  for (int loop = 0; loop < nprocs; loop++) {
    if (me != next) {
      MPI_Irecv(bufcopy,maxbytes,MPI_CHAR,prev,messtag,world,&request);
      MPI_Send(buf,nbytes,MPI_CHAR,next,messtag,world);
      MPI_Wait(&request,&status);
      MPI_Get_count(&status,MPI_CHAR,&nbytes);
      memcpy(buf,bufcopy,nbytes);
    }
    if (self || loop != nprocs-1) callback(nbytes/nper,buf);
  }

  if (outbuf) memcpy(outbuf,buf,nbytes);

  memory->destroy(buf);
  memory->destroy(bufcopy);
}

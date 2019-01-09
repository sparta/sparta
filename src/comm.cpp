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
#include "surf.h"
#include "update.h"
#include "adapt_grid.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{PKEEP,PINSERT,PDONE,PDISCARD,PENTRY,PEXIT,PSURF};   // several files

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

  iparticle = new Irregular(sparta);
  igrid = NULL;
  iuniform = NULL;

  pproc = NULL;
  maxpproc = 0;
  gproc = gsize = NULL;
  maxgproc = 0;
  sbuf = rbuf = NULL;
  maxsendbuf = maxrecvbuf = 0;

  copymode = 0;
}

/* ---------------------------------------------------------------------- */

Comm::~Comm()
{
  if (copymode) return;

  delete iparticle;
  delete igrid;
  delete iuniform;

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
  
  iparticle->create_procs(nneigh,neighlist,commsortflag);
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

  int ncustom = particle->ncustom;
  int nbytes_particle = sizeof(Particle::OnePart);
  int nbytes_custom = particle->sizeof_custom();
  int nbytes = nbytes_particle + nbytes_custom;

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
  // if no custom attributes, pack particles directly via memcpy()
  // else pack_custom() performs packing into sbuf

  int nsend = 0;
  int offset = 0;

  if (!ncustom) {
    for (i = 0; i < nmigrate; i++) {
      j = plist[i];
      if (particles[j].flag == PDISCARD) continue;
      pproc[nsend++] = cells[particles[j].icell].proc;
      particles[j].icell = cells[particles[j].icell].ilocal;
      memcpy(&sbuf[offset],&particles[j],nbytes_particle);
      offset += nbytes_particle;
    }
  } else {
    for (i = 0; i < nmigrate; i++) {
      j = plist[i];
      if (particles[j].flag == PDISCARD) continue;
      pproc[nsend++] = cells[particles[j].icell].proc;
      particles[j].icell = cells[particles[j].icell].ilocal;
      memcpy(&sbuf[offset],&particles[j],nbytes_particle);
      offset += nbytes_particle;
      particle->pack_custom(j,&sbuf[offset]);
      offset += nbytes_custom;
    }
  }

  // compress my list of particles

  particle->compress_migrate(nmigrate,plist);
  int ncompress = particle->nlocal;

  // create or augment irregular communication plan
  // nrecv = # of incoming particles
  
  int nrecv;
  if (neighflag)
    nrecv = iparticle->augment_data_uniform(nsend,pproc);
  else 
    nrecv = iparticle->create_data_uniform(nsend,pproc,commsortflag);

  // extend particle list if necessary

  particle->grow(nrecv);

  // perform irregular communication
  // if no custom attributes, append recv particles directly to particle list
  // else receive into rbuf, unpack particles one by one via unpack_custom()

  if (!ncustom)
    iparticle->
      exchange_uniform(sbuf,nbytes,
                       (char *) &particle->particles[particle->nlocal]);

  else {
    if (nrecv*nbytes > maxrecvbuf) {
      maxrecvbuf = nrecv*nbytes;
      memory->destroy(rbuf);
      memory->create(rbuf,maxrecvbuf,"comm:rbuf");
    }

    iparticle->exchange_uniform(sbuf,nbytes,rbuf);

    offset = 0;
    int nlocal = particle->nlocal;
    for (i = 0; i < nrecv; i++) {
      memcpy(&particle->particles[nlocal],&rbuf[offset],nbytes_particle);
      offset += nbytes_particle;
      particle->unpack_custom(&rbuf[offset],nlocal);
      offset += nbytes_custom;
      nlocal++;
    }
  }

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
  if (update->mem_limit_grid_flag)
    update->global_mem_limit = grid->nlocal*sizeof(Grid::ChildCell);
  if (update->global_mem_limit > 0)
    return migrate_cells_less_memory(nmigrate);

  int i,n;

  Grid::ChildCell *cells = grid->cells;
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

  // compress implicit surf list to remove surfs in migrating cells
  //   do before grid->compress() b/c uses grid hash
  // compress my list of owned grid cells to remove migrating cells
  // compress particle list to remove particles in migrating cells
  // unset particle sorted since compress_rebalance() does
  //   and may receive particles which will then be unsorted

  if (nmigrate) {
    if (surf->implicit) surf->compress_rebalance();
    grid->compress();
    particle->compress_rebalance();
  } else particle->sorted = 0;

  // create irregular communication plan with variable size datums
  // nrecv = # of incoming grid cells
  // recvsize = total byte size of incoming grid + particle info
  // DEBUG: append a sort=1 arg so that messages from other procs
  //        are received in repeatable order, thus grid cells stay in order

  if (!igrid) igrid = new Irregular(sparta);
  int recvsize;
  int nrecv = igrid->create_data_variable(nmigrate,gproc,gsize,
                                          recvsize,commsortflag);

  // reallocate rbuf as needed

  if (recvsize > maxrecvbuf) {
    memory->destroy(rbuf);
    maxrecvbuf = recvsize;
    memory->create(rbuf,maxrecvbuf,"comm:rbuf");
    memset(rbuf,0,maxrecvbuf);
  }

  // perform irregular communication

  igrid->exchange_variable(sbuf,gsize,rbuf);

  // unpack received grid cells with their particles

  offset = 0;
  for (i = 0; i < nrecv; i++)
    offset += grid->unpack_one(&rbuf[offset],1,1);
}

/* ----------------------------------------------------------------------
   migrate grid cells with their particles to new procs
   called from BalanceGrid and FixBalance
   uses multiple comm passes to reduce buffer size
------------------------------------------------------------------------- */

void Comm::migrate_cells_less_memory(int nmigrate)
{
  int i,n;

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

  int icell_start = 0;
  int icell_end = grid->nlocal;
  int not_done = 1;
  int nglocal = grid->nlocal;

  while (not_done) {
    Grid::ChildCell *cells = grid->cells;

    int nsend = 0;
    bigint boffset = 0;
    for (int icell = icell_start; icell < nglocal; icell++) {
      icell_end = icell+1;
      if (cells[icell].nsplit <= 0) continue;
      if (cells[icell].proc == me) continue;
      gproc[nsend] = cells[icell].proc;
      n = grid->pack_one(icell,NULL,1,1,0);
      if (n > 0 && boffset > 0 && boffset+n > update->global_mem_limit) {
        icell_end -= 1;
        break;
      }
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
    for (int icell = icell_start; icell < icell_end; icell++) {
      if (cells[icell].nsplit <= 0) continue;
      if (cells[icell].proc == me) continue;
      offset += grid->pack_one(icell,&sbuf[offset],1,1,1);
    }

    // compress particle list to remove particles in migrating cells

    if (nmigrate) particle->compress_rebalance_sorted();

    // create irregular communication plan with variable size datums
    // nrecv = # of incoming grid cells
    // recvsize = total byte size of incoming grid + particle info
    // DEBUG: append a sort=1 arg so that messages from other procs
    //        are received in repeatable order, thus grid cells stay in order

    if (!igrid) igrid = new Irregular(sparta);
    int recvsize;
    int nrecv = igrid->create_data_variable(nsend,gproc,gsize,
                                            recvsize,commsortflag);

    // reallocate rbuf as needed

    if (recvsize > maxrecvbuf) {
      memory->destroy(rbuf);
      maxrecvbuf = recvsize;
      memory->create(rbuf,maxrecvbuf,"comm:rbuf");
      memset(rbuf,0,maxrecvbuf);
    }

    // perform irregular communication

    igrid->exchange_variable(sbuf,gsize,rbuf);

    // unpack received grid cells with their particles
    // set unpack_one() sortflag arg to keep new particles sorted

    offset = 0;
    for (i = 0; i < nrecv; i++)
      offset += grid->unpack_one(&rbuf[offset],1,1,1);

    // deallocate large buffers to reduce memory footprint
    // also deallocate igrid for same reason

    if (sbuf) memory->destroy(sbuf);
    sbuf = NULL;
    maxsendbuf = 0;

    if (rbuf) memory->destroy(rbuf);
    rbuf = NULL;
    maxrecvbuf = 0;

    delete igrid;
    igrid = NULL;

    icell_start = icell_end;
    int not_done_local = icell_start < nglocal;
    MPI_Allreduce(&not_done_local,&not_done,1,MPI_INT,MPI_SUM,world); 
  }

  // compress my list of owned grid cells to remove migrated cells

  if (nmigrate) grid->compress();
}

/* ----------------------------------------------------------------------
   send grid cell info with their particles needed for possible grid adaptation
   return # of received cells and buf = ptr to received cell info
   called from AdaptGrid
------------------------------------------------------------------------- */

int Comm::send_cells_adapt(int nsend, int *procsend, char *inbuf, char **outbuf)
{
  int i,n;

  AdaptGrid::SendAdapt *sadapt = (AdaptGrid::SendAdapt *) inbuf;

  // grow size list if needed
  // don't use gproc, but needs to stay same size as gsize

  if (nsend > maxgproc) {
    maxgproc = nsend;
    memory->destroy(gproc);
    memory->destroy(gsize);
    memory->create(gproc,maxgproc,"comm:gproc");
    memory->create(gsize,maxgproc,"comm:gsize");
  }

  // compute byte count needed to pack cells

  bigint boffset = 0;
  for (i = 0; i < nsend; i++) {
    n = grid->pack_one_adapt((char *) &sadapt[i],NULL,0);
    gsize[i] = n;
    boffset += n;
  }

  if (boffset > MAXSMALLINT) 
    error->one(FLERR,"Adapt grid send buffer exceeds 2 GB");
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
  for (i = 0; i < nsend; i++)
    offset += grid->pack_one_adapt((char *) &sadapt[i],&sbuf[offset],1);

  // create irregular communication plan with variable size datums
  // nrecv = # of incoming grid cells
  // recvsize = total byte size of incoming grid + particle info
  // DEBUG: append a sort=1 arg so that messages from other procs
  //        are received in repeatable order, thus grid cells stay in order

  if (!igrid) igrid = new Irregular(sparta);
  int recvsize;
  int nrecv = 
    igrid->create_data_variable(nsend,procsend,gsize,
                                recvsize,commsortflag);

  // reallocate rbuf as needed

  if (recvsize > maxrecvbuf) {
    memory->destroy(rbuf);
    maxrecvbuf = recvsize;
    memory->create(rbuf,maxrecvbuf,"comm:rbuf");
    memset(rbuf,0,maxrecvbuf);
  }

  // perform irregular communication

  igrid->exchange_variable(sbuf,gsize,rbuf);

  // return rbuf and grid cell count

  *outbuf = rbuf;
  return nrecv;
}

/* ----------------------------------------------------------------------
   wrapper on irregular comm of datums on uniform size
   called from AdaptGrid
------------------------------------------------------------------------- */

int Comm::irregular_uniform(int nsend, int *procsend, 
                            char *inbuf, int nsize, char **outbuf)
{
  // create irregular communication plan with constant size datums
  // nrecv = # of incoming grid cells
  // DEBUG: append a sort=1 arg so that messages from other procs
  //        are received in repeatable order, thus grid cells stay in order

  if (!iuniform) iuniform = new Irregular(sparta);
  int nrecv = iuniform->create_data_uniform(nsend,procsend,commsortflag);

  // reallocate rbuf as needed

  if (nrecv*nsize > maxrecvbuf) {
    memory->destroy(rbuf);
    maxrecvbuf = nrecv*nsize;
    memory->create(rbuf,maxrecvbuf,"comm:rbuf");
    memset(rbuf,0,maxrecvbuf);
  }

  // perform irregular communication

  iuniform->exchange_uniform(inbuf,nsize,rbuf);

  // return rbuf and grid cell count

  *outbuf = rbuf;
  return nrecv;
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

/* ----------------------------------------------------------------------
   rendezvous communication operation
   three stages:
     first Irregular converts inbuf from caller decomp to rvous decomp
     callback operates on data in rendevous decomp
     last Irregular converts outbuf from rvous decomp back to caller decomp
   inputs:
     n = # of input datums
     proclist = proc that owns each input datum in rendezvous decomposition
     inbuf = list of input datums
     insize = size in bytes of each input datum
     callback = caller function to invoke in rendezvous decomposition
   outputs:
     nout = # of output datums (function return)
     outbuf = list of output datums
     outsize = size in bytes of each output datum
------------------------------------------------------------------------- */

int Comm::rendezvous(int n, int *proclist, char *inbuf, int insize,
                     int (*callback)(int, char *, int *&, char *&, void *),
                     char *&outbuf, int outsize, void *ptr)
{
  // comm inbuf from caller decomposition to rendezvous decomposition

  Irregular *irregular = new Irregular(sparta);

  int n_rvous = irregular->create_data_uniform(n,proclist);  // add sort
  char *inbuf_rvous = (char *) memory->smalloc((bigint) n_rvous*insize,
                                               "rendezvous:inbuf_rvous");
  irregular->exchange_uniform(inbuf,insize,inbuf_rvous);

  delete irregular;

  // peform rendezvous computation via callback()
  // callback() allocates/populates proclist_rvous and outbuf_rvous

  int *proclist_rvous;
  char *outbuf_rvous;

  int nout_rvous = 
    callback(n_rvous,inbuf_rvous,proclist_rvous,outbuf_rvous,ptr);

  memory->sfree(inbuf_rvous);

  // comm outbuf from rendezvous decomposition back to caller
  // caller will free outbuf

  irregular = new Irregular(sparta);
  
  int nout = irregular->create_data_uniform(nout_rvous,proclist_rvous);
  outbuf = (char *) memory->smalloc((bigint) nout*outsize,"rendezvous:outbuf");
  irregular->exchange_uniform(outbuf_rvous,outsize,outbuf);
  
  delete irregular;
  memory->destroy(proclist_rvous);
  memory->sfree(outbuf_rvous);

  // return number of datums

  return nout;
}

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

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "comm.h"
#include "irregular.h"
#include "particle.h"
#include "grid.h"
#include "surf.h"
#include "update.h"
#include "modify.h"
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

  if (me == 0) {
    if (screen) fprintf(screen,"Running on %d MPI task(s)\n",nprocs);
    if (logfile) fprintf(logfile,"Running on %d MPI task(s)\n",nprocs);
  }

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

  // for optimized particle moves, call compress_reactions rather than
  //  compress_migrate since mlist is not guaranteed to be in ascending
  //  order

  if (update->optmove_flag)
    particle->compress_reactions(nmigrate,plist);
  else
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
  if (update->have_mem_limit())
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
    n = grid->pack_one(icell,NULL,1,1,1,0);
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
    offset += grid->pack_one(icell,&sbuf[offset],1,1,1,1);
  }

  // compress my list of owned implicit surfs, resets csurfs in kept cells
  // compress my list of owned grid cells to remove migrating cells
  // compress my list of owned distributed/explicit surfs
  // compress particle list to remove particles in migrating cells
  // procs with no migrating cells must also unset particle sorted
  //   since compress_rebalance() unsets it

  if (nmigrate) {
    if (surf->implicit) surf->compress_implicit();
    grid->compress();
    if (surf->distributed && !surf->implicit) surf->compress_explicit();
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
    offset += grid->unpack_one(&rbuf[offset],1,1,1);
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
      n = grid->pack_one(icell,NULL,1,1,1,0);
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
      offset += grid->pack_one(icell,&sbuf[offset],1,1,1,1);
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
      offset += grid->unpack_one(&rbuf[offset],1,1,1,1);

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

  // compress my list of owned implicit surfs, resets csurfs in kept cells
  // compress my list of owned grid cells to remove migrated cells

  if (nmigrate) {
    if (surf->implicit) surf->compress_implicit();
    grid->compress();
    if (surf->distributed && !surf->implicit) surf->compress_explicit();
  }
}

/* ----------------------------------------------------------------------
   send grid cell info with their surfs/particles needed for grid adaptation
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
    igrid->create_data_variable(nsend,procsend,gsize,recvsize,commsortflag);

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
   wrapper on irregular comm of datums of uniform size
   receiving procs are neighbor procs of my owned grid cells, not including self
   so can use same comm calls as migrate_particles (if neighflag is set)
     via iparticle->augment_data_uniform() and exchange_uniform()
     otherwise same logic as irregular_uniform()
   called from FixAblate
------------------------------------------------------------------------- */

int Comm::irregular_uniform_neighs(int nsend, int *procsend,
                                   char *inbuf, int nsize, char **outbuf)
{
  // if neighflag, use iparticle
  // else one-time create of irregular comm plan with constant size datums
  // nrecv = # of incoming grid cells

  if (!neighflag && !iuniform) iuniform = new Irregular(sparta);

  int nrecv;
  if (neighflag)
    nrecv = iparticle->augment_data_uniform(nsend,procsend);
  else
    nrecv = iuniform->create_data_uniform(nsend,procsend,commsortflag);

  // reallocate rbuf as needed

  if (nrecv*nsize > maxrecvbuf) {
    memory->destroy(rbuf);
    maxrecvbuf = nrecv*nsize;
    memory->create(rbuf,maxrecvbuf,"comm:rbuf");
    memset(rbuf,0,maxrecvbuf);
  }

  // perform irregular communication

  if (neighflag)
    iparticle->exchange_uniform(inbuf,nsize,rbuf);
  else
    iuniform->exchange_uniform(inbuf,nsize,rbuf);

  // return rbuf and grid cell count

  *outbuf = rbuf;
  return nrecv;
}

/* ----------------------------------------------------------------------
   wrapper on irregular comm of datums of uniform size
   receiving procs can be anyone, including self
   use create_data_uniform() and exchange_uniform()
   called from AdaptGrid
------------------------------------------------------------------------- */

int Comm::irregular_uniform(int nsend, int *procsend,
                            char *inbuf, int nsize, char **outbuf)
{
  // one-time create of irregular comm plan with constant size datums
  // nrecv = # of incoming grid cells

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
                void (*callback)(int, char *, void *), void *outbuf, int self,
                void *ptr)
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
    if (self || loop != nprocs-1) callback(nbytes/nper,buf,ptr);
  }

  if (outbuf) memcpy(outbuf,buf,nbytes);

  memory->destroy(buf);
  memory->destroy(bufcopy);
}

/* ----------------------------------------------------------------------
   rendezvous communication operation
   three stages:
     first comm sends inbuf from caller decomp to rvous decomp
     callback operates on data in rendevous decomp
     second comm sends outbuf from rvous decomp back to caller decomp
   inputs:
     which = perform (0) irregular or (1) MPI_All2allv communication
     n = # of datums in inbuf
     inbuf = vector of input datums
     insize = byte size of each input datum
     inorder = 0 for inbuf in random proc order, 1 for datums ordered by proc
     procs: inorder 0 = proc to send each datum to, 1 = # of datums/proc,
     callback = caller function to invoke in rendezvous decomposition
                takes input datums, returns output datums
     outorder = same as inorder, but for datums returned by callback()
     ptr = pointer to caller class, passed to callback()
     statflag = 1 for stats output, else 0
   outputs:
     nout = # of output datums (function return)
     outbuf = vector of output datums
     outsize = byte size of each output datum
   callback inputs:
     nrvous = # of rvous decomp datums in inbuf_rvous
     inbuf_rvous = vector of rvous decomp input datums
     ptr = pointer to caller class
   callback outputs:
     nrvous_out = # of rvous decomp output datums (function return)
     flag = 0 for no second comm, 1 for outbuf_rvous = inbuf_rvous,
            2 for second comm with new outbuf_rvous
     procs_rvous = outorder 0 = proc to send each datum to, 1 = # of datums/proc
                   allocated
     outbuf_rvous = vector of rvous decomp output datums
   NOTE: could use MPI_INT or MPI_DOUBLE insead of MPI_CHAR
         to avoid checked-for overflow in MPI_Alltoallv?
------------------------------------------------------------------------- */

int Comm::
rendezvous(int which, int n, char *inbuf, int insize,
           int inorder, int *procs,
           int (*callback)(int, char *, int &, int *&, char *&, void *),
           int outorder, char *&outbuf, int outsize, void *ptr, int statflag)
{
  if (which == 0)
    return rendezvous_irregular(n,inbuf,insize,inorder,procs,callback,
                                outorder,outbuf,outsize,ptr,statflag);
  else
    return rendezvous_all2all(n,inbuf,insize,inorder,procs,callback,
                              outorder,outbuf,outsize,ptr,statflag);
}

/* ---------------------------------------------------------------------- */

int Comm::
rendezvous_irregular(int n, char *inbuf, int insize, int inorder, int *procs,
                     int (*callback)(int, char *, int &, int *&, char *&, void *),
                     int outorder, char *&outbuf,
                     int outsize, void *ptr, int statflag)
{
  // irregular comm of inbuf from caller decomp to rendezvous decomp

  Irregular *irregular = new Irregular(sparta);

  int nrvous;
  if (inorder) nrvous = irregular->create_data_uniform_grouped(n,procs);
  else nrvous = irregular->create_data_uniform(n,procs);

  char *inbuf_rvous = (char *) memory->smalloc((bigint) nrvous*insize,
                                               "rendezvous:inbuf");
  irregular->exchange_uniform(inbuf,insize,inbuf_rvous);

  bigint irregular1_bytes = 0;   // irregular->irregular_bytes;
  delete irregular;

  // done if callback is NULL, return inbuf_rvous

  if (!callback) {
    outbuf = inbuf_rvous;
    return nrvous;
  }

  // peform rendezvous computation via callback()
  // callback() allocates/populates proclist_rvous and outbuf_rvous

  int flag;
  int *procs_rvous;
  char *outbuf_rvous;
  int nrvous_out = callback(nrvous,inbuf_rvous,flag,
                            procs_rvous,outbuf_rvous,ptr);

  if (flag != 1) memory->sfree(inbuf_rvous);  // outbuf_rvous = inbuf_vous
  if (flag == 0) return 0;    // all nout_rvous are 0, no 2nd comm stage

  // irregular comm of outbuf from rendezvous decomp back to caller decomp
  // caller will free outbuf

  irregular = new Irregular(sparta);

  int nout;
  if (outorder)
    nout = irregular->create_data_uniform_grouped(nrvous_out,procs_rvous);
  else nout = irregular->create_data_uniform(nrvous_out,procs_rvous);

  outbuf = (char *) memory->smalloc((bigint) nout*outsize,
                                    "rendezvous:outbuf");
  irregular->exchange_uniform(outbuf_rvous,outsize,outbuf);

  bigint irregular2_bytes = 0;   // irregular->irregular_bytes;
  delete irregular;

  memory->destroy(procs_rvous);
  memory->sfree(outbuf_rvous);

  // return number of output datums

  if (!statflag) return nout;

  rendezvous_stats(n,insize,nout,outsize,nrvous,nrvous_out);

  /*
  rvous_bytes = 0;
  rvous_bytes += n*insize;                                // inbuf
  rvous_bytes += nout*outsize;                            // outbuf
  rvous_bytes += nrvous*insize;                           // inbuf_rvous
  rvous_bytes += nrvous_out*outsize;                      // outbuf_rvous
  rvous_bytes += nrvous_out*sizeof(int);                  // procs_rvous
  rvous_bytes += MAX(irregular1_bytes,irregular2_bytes);  // max of 2 comms
  */

  return nout;
}

/* ---------------------------------------------------------------------- */

int Comm::
rendezvous_all2all(int n, char *inbuf, int insize, int inorder, int *procs,
                   int (*callback)(int, char *, int &, int *&, char *&, void *),
                   int outorder, char *&outbuf, int outsize, void *ptr,
                   int statflag)
{
  int iproc;
  bigint all2all1_bytes,all2all2_bytes;
  int *sendcount,*sdispls,*recvcount,*rdispls;
  int *procs_a2a;
  bigint *offsets;
  char *inbuf_a2a,*outbuf_a2a;

  // create procs and inbuf for All2all if necesary

  if (!inorder) {
    memory->create(procs_a2a,nprocs,"rendezvous:procs");
    inbuf_a2a = (char *) memory->smalloc((bigint) n*insize,
                                         "rendezvous:inbuf");
    memory->create(offsets,nprocs,"rendezvous:offsets");

    for (int i = 0; i < nprocs; i++) procs_a2a[i] = 0;
    for (int i = 0; i < n; i++) procs_a2a[procs[i]]++;

    offsets[0] = 0;
    for (int i = 1; i < nprocs; i++)
      offsets[i] = offsets[i-1] + insize*procs_a2a[i-1];

    bigint offset = 0;
    for (int i = 0; i < n; i++) {
      iproc = procs[i];
      memcpy(&inbuf_a2a[offsets[iproc]],&inbuf[offset],insize);
      offsets[iproc] += insize;
      offset += insize;
    }

    all2all1_bytes = nprocs*sizeof(int) + nprocs*sizeof(bigint) + n*insize;

  } else {
    procs_a2a = procs;
    inbuf_a2a = inbuf;
    all2all1_bytes = 0;
  }

  // create args for MPI_Alltoallv() on input data

  memory->create(sendcount,nprocs,"rendezvous:sendcount");
  memcpy(sendcount,procs_a2a,nprocs*sizeof(int));

  memory->create(recvcount,nprocs,"rendezvous:recvcount");
  MPI_Alltoall(sendcount,1,MPI_INT,recvcount,1,MPI_INT,world);

  memory->create(sdispls,nprocs,"rendezvous:sdispls");
  memory->create(rdispls,nprocs,"rendezvous:rdispls");
  sdispls[0] = rdispls[0] = 0;
  for (int i = 1; i < nprocs; i++) {
    sdispls[i] = sdispls[i-1] + sendcount[i-1];
    rdispls[i] = rdispls[i-1] + recvcount[i-1];
  }
  int nrvous = rdispls[nprocs-1] + recvcount[nprocs-1];

  // test for overflow of input data due to imbalance or insize
  // means that individual sdispls or rdispls values overflow

  int overflow = 0;
  if ((bigint) n*insize > MAXSMALLINT) overflow = 1;
  if ((bigint) nrvous*insize > MAXSMALLINT) overflow = 1;
  int overflowall;
  MPI_Allreduce(&overflow,&overflowall,1,MPI_INT,MPI_MAX,world);
  if (overflowall) error->all(FLERR,"Overflow input size in rendezvous_a2a");

  for (int i = 0; i < nprocs; i++) {
    sendcount[i] *= insize;
    sdispls[i] *= insize;
    recvcount[i] *= insize;
    rdispls[i] *= insize;
  }

  // all2all comm of inbuf from caller decomp to rendezvous decomp

  char *inbuf_rvous = (char *) memory->smalloc((bigint) nrvous*insize,
                                               "rendezvous:inbuf");

  MPI_Alltoallv(inbuf_a2a,sendcount,sdispls,MPI_CHAR,
                inbuf_rvous,recvcount,rdispls,MPI_CHAR,world);

  if (!inorder) {
    memory->destroy(procs_a2a);
    memory->sfree(inbuf_a2a);
    memory->destroy(offsets);
  }

  // done if callback is NULL, return inbuf_rvous

  if (!callback) {
    memory->destroy(sendcount);
    memory->destroy(recvcount);
    memory->destroy(sdispls);
    memory->destroy(rdispls);
    outbuf = inbuf_rvous;
    return nrvous;
  }

  // peform rendezvous computation via callback()
  // callback() allocates/populates proclist_rvous and outbuf_rvous

  int flag;
  int *procs_rvous;
  char *outbuf_rvous;

  int nrvous_out = callback(nrvous,inbuf_rvous,flag,
                            procs_rvous,outbuf_rvous,ptr);

  if (flag != 1) memory->sfree(inbuf_rvous);  // outbuf_rvous = inbuf_vous
  if (flag == 0) {
    memory->destroy(sendcount);
    memory->destroy(recvcount);
    memory->destroy(sdispls);
    memory->destroy(rdispls);
    return 0;    // all nout_rvous are 0, no 2nd irregular
  }

  // create procs and outbuf for All2all if necesary

  if (!outorder) {
    memory->create(procs_a2a,nprocs,"rendezvous_a2a:procs");

    outbuf_a2a = (char *) memory->smalloc((bigint) nrvous_out*outsize,
                                          "rendezvous:outbuf");
    memory->create(offsets,nprocs,"rendezvous:offsets");

    for (int i = 0; i < nprocs; i++) procs_a2a[i] = 0;
    for (int i = 0; i < nrvous_out; i++) procs_a2a[procs_rvous[i]]++;

    offsets[0] = 0;
    for (int i = 1; i < nprocs; i++)
      offsets[i] = offsets[i-1] + outsize*procs_a2a[i-1];

    bigint offset = 0;
    for (int i = 0; i < nrvous_out; i++) {
      iproc = procs_rvous[i];
      memcpy(&outbuf_a2a[offsets[iproc]],&outbuf_rvous[offset],outsize);
      offsets[iproc] += outsize;
      offset += outsize;
    }

    all2all2_bytes = nprocs*sizeof(int) + nprocs*sizeof(bigint) +
      nrvous_out*outsize;

  } else {
    procs_a2a = procs_rvous;
    outbuf_a2a = outbuf_rvous;
    all2all2_bytes = 0;
  }

  // comm outbuf from rendezvous decomposition back to caller

  memcpy(sendcount,procs_a2a,nprocs*sizeof(int));

  MPI_Alltoall(sendcount,1,MPI_INT,recvcount,1,MPI_INT,world);

  sdispls[0] = rdispls[0] = 0;
  for (int i = 1; i < nprocs; i++) {
    sdispls[i] = sdispls[i-1] + sendcount[i-1];
    rdispls[i] = rdispls[i-1] + recvcount[i-1];
  }
  int nout = rdispls[nprocs-1] + recvcount[nprocs-1];

  // test for overflow of outbuf due to imbalance or outsize
  // means that individual sdispls or rdispls values overflow

  overflow = 0;
  if ((bigint) nrvous*outsize > MAXSMALLINT) overflow = 1;
  if ((bigint) nout*outsize > MAXSMALLINT) overflow = 1;
  MPI_Allreduce(&overflow,&overflowall,1,MPI_INT,MPI_MAX,world);
  if (overflowall) error->all(FLERR,"Overflow output in rendezvous_a2a");

  for (int i = 0; i < nprocs; i++) {
    sendcount[i] *= outsize;
    sdispls[i] *= outsize;
    recvcount[i] *= outsize;
    rdispls[i] *= outsize;
  }

  // all2all comm of outbuf from rendezvous decomp back to caller decomp
  // caller will free outbuf

  outbuf = (char *) memory->smalloc((bigint) nout*outsize,"rendezvous:outbuf");

  MPI_Alltoallv(outbuf_a2a,sendcount,sdispls,MPI_CHAR,
                outbuf,recvcount,rdispls,MPI_CHAR,world);

  memory->destroy(procs_rvous);
  memory->sfree(outbuf_rvous);

  if (!outorder) {
    memory->destroy(procs_a2a);
    memory->sfree(outbuf_a2a);
    memory->destroy(offsets);
  }

  // clean up

  memory->destroy(sendcount);
  memory->destroy(recvcount);
  memory->destroy(sdispls);
  memory->destroy(rdispls);

  // return number of datums

  if (!statflag) return nout;

  rendezvous_stats(n,insize,nout,outsize,nrvous,nrvous_out);

  /*
  rvous_bytes = 0;
  rvous_bytes += n*insize;                                // inbuf
  rvous_bytes += nout*outsize;                            // outbuf
  rvous_bytes += nrvous*insize;                           // inbuf_rvous
  rvous_bytes += nrvous_out*outsize;                      // outbuf_rvous
  rvous_bytes += nrvous_out*sizeof(int);                  // procs_rvous
  rvous_bytes += 4*nprocs*sizeof(int);                    // all2all vectors
  rvous_bytes += MAX(all2all1_bytes,all2all2_bytes);      // reorder ops
  */

  return nout;
}

/* ----------------------------------------------------------------------
   memory info for caller and rendezvous decompositions
------------------------------------------------------------------------- */

void Comm::rendezvous_stats(int n, int insize, int nout, int outsize,
                            int nrvous, int nrvous_out)
{
  bigint size_in_all,size_in_max,size_in_min;
  bigint size_out_all,size_out_max,size_out_min;
  bigint size_inrvous_all,size_inrvous_max,size_inrvous_min;
  bigint size_outrvous_all,size_outrvous_max,size_outrvous_min;

  bigint size = (bigint) n*insize;
  MPI_Allreduce(&size,&size_in_all,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  MPI_Allreduce(&size,&size_in_max,1,MPI_SPARTA_BIGINT,MPI_MAX,world);
  MPI_Allreduce(&size,&size_in_min,1,MPI_SPARTA_BIGINT,MPI_MIN,world);

  size = (bigint) nout*outsize;
  MPI_Allreduce(&size,&size_out_all,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  MPI_Allreduce(&size,&size_out_max,1,MPI_SPARTA_BIGINT,MPI_MAX,world);
  MPI_Allreduce(&size,&size_out_min,1,MPI_SPARTA_BIGINT,MPI_MIN,world);

  size = (bigint) nrvous*insize;
  MPI_Allreduce(&size,&size_inrvous_all,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  MPI_Allreduce(&size,&size_inrvous_max,1,MPI_SPARTA_BIGINT,MPI_MAX,world);
  MPI_Allreduce(&size,&size_inrvous_min,1,MPI_SPARTA_BIGINT,MPI_MIN,world);

  size = (bigint) nrvous_out*insize;
  MPI_Allreduce(&size,&size_outrvous_all,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  MPI_Allreduce(&size,&size_outrvous_max,1,MPI_SPARTA_BIGINT,MPI_MAX,world);
  MPI_Allreduce(&size,&size_outrvous_min,1,MPI_SPARTA_BIGINT,MPI_MIN,world);

  int mbytes = 1024*1024;

  if (me == 0) {
    if (screen) {
      fprintf(screen,"Rendezvous balance and memory info:\n");
      fprintf(screen,"  input datum count "
              "(tot,ave,max,min): " BIGINT_FORMAT " %g "
              BIGINT_FORMAT " " BIGINT_FORMAT "\n",
              size_in_all/insize,1.0*size_in_all/nprocs/insize,
              size_in_max/insize,size_in_min/insize);
      fprintf(screen,"  input data (MB) "
              "(tot,ave,max,min): %g %g %g %g\n",
              1.0*size_in_all/mbytes,1.0*size_in_all/nprocs/mbytes,
              1.0*size_in_max/mbytes,1.0*size_in_min/mbytes);
      fprintf(screen,"  output datum count "
              "(tot,ave,max,min): " BIGINT_FORMAT " %g "
              BIGINT_FORMAT " " BIGINT_FORMAT "\n",
              size_out_all/outsize,1.0*size_out_all/nprocs/outsize,
              size_out_max/outsize,size_out_min/outsize);
      fprintf(screen,"  output data (MB) "
              "(tot,ave,max,min): %g %g %g %g\n",
              1.0*size_out_all/mbytes,1.0*size_out_all/nprocs/mbytes,
              1.0*size_out_max/mbytes,1.0*size_out_min/mbytes);
      fprintf(screen,"  input rvous datum count "
              "(tot,ave,max,min): " BIGINT_FORMAT " %g "
              BIGINT_FORMAT " " BIGINT_FORMAT "\n",
              size_inrvous_all/insize,1.0*size_inrvous_all/nprocs/insize,
              size_inrvous_max/insize,size_inrvous_min/insize);
      fprintf(screen,"  input rvous data (MB) "
              "(tot,ave,max,min): %g %g %g %g\n",
              1.0*size_inrvous_all/mbytes,1.0*size_inrvous_all/nprocs/mbytes,
              1.0*size_inrvous_max/mbytes,1.0*size_inrvous_min/mbytes);
      fprintf(screen,"  output rvous datum count "
              "(tot,ave,max,min): " BIGINT_FORMAT " %g "
              BIGINT_FORMAT " " BIGINT_FORMAT "\n",
              size_outrvous_all/outsize,1.0*size_outrvous_all/nprocs/outsize,
              size_outrvous_max/outsize,size_outrvous_min/outsize);
      fprintf(screen,"  output rvous data (MB) "
              "(tot,ave,max,min): %g %g %g %g\n",
              1.0*size_outrvous_all/mbytes,1.0*size_outrvous_all/nprocs/mbytes,
              1.0*size_outrvous_max/mbytes,1.0*size_outrvous_min/mbytes);
    }
  }
}

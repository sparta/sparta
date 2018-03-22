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
#include "comm_kokkos.h"
#include "collide_vss_kokkos.h"
#include "irregular_kokkos.h"
#include "particle_kokkos.h"
#include "grid_kokkos.h"
#include "update.h"
//#include "adapt_grid.h"
#include "memory_kokkos.h"
#include "error.h"
#include "kokkos.h"
#include "sparta_masks.h"

using namespace SPARTA_NS;

enum{PKEEP,PINSERT,PDONE,PDISCARD,PENTRY,PEXIT,PSURF};   // several files

/* ---------------------------------------------------------------------- */

CommKokkos::CommKokkos(SPARTA *sparta) : Comm(sparta)
{
  delete iparticle;
  iparticle = new IrregularKokkos(sparta);

  k_nsend = DAT::tdual_int_scalar("comm:nsend");
  d_nsend = k_nsend.view<DeviceType>();
  h_nsend = k_nsend.h_view;
}

/* ---------------------------------------------------------------------- */

CommKokkos::~CommKokkos()
{
  if (copymode) return;

  if (!sparta->kokkos->comm_classic) {
    memoryKK->destroy_kokkos(k_pproc,pproc);
    pproc = NULL;
  }
}

/* ----------------------------------------------------------------------
   migrate particles to new procs after particle move
   return particle nlocal after compression, 
     so Update can iterate on particle move
------------------------------------------------------------------------- */

int CommKokkos::migrate_particles(int nmigrate, int *plist, DAT::t_int_1d d_plist_in)
{
  GridKokkos* grid_kk = (GridKokkos*) grid;
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;

  if (sparta->kokkos->comm_classic) {
    particle_kk->sync(Host,ALL_MASK);
    //grid_kk->sync(Host,ALL_MASK);
    int prev_auto_sync = sparta->kokkos->auto_sync;
    sparta->kokkos->auto_sync = 1;

    int ncompress = Comm::migrate_particles(nmigrate,plist);

    particle_kk->sync(Device,ALL_MASK);
    //grid_kk->sync(Device,ALL_MASK);
    sparta->kokkos->auto_sync = prev_auto_sync;

    return ncompress;
  }

  // int i,j;

  d_plist = d_plist_in;


  int ncustom = particle->ncustom;
  nbytes_particle = sizeof(Particle::OnePart);
  int nbytes_custom = particle->sizeof_custom();
  int nbytes = nbytes_particle + nbytes_custom;

  // Kokkos
  // memory access: cells, particles
  // local views: pproc, sbuf, rbuf
  // functions: particle->compress_migrate, grow,
  //  iparticle->augment_data_uniform, iparticle->create_data_uniform,
  //  iparticle->exchange_uniform
  // parallel_for: loop over nmigrate to pack buffer,
  //  loop over nrecv to unpack buffer
  // atomic variables: ?

  // grow pproc and sbuf if necessary

  if (nmigrate > maxpproc) {
    maxpproc = nmigrate;
    memoryKK->destroy_kokkos(k_pproc,pproc);
    memoryKK->create_kokkos(k_pproc,pproc,maxpproc,"comm:pproc");
    d_pproc = k_pproc.d_view;
  }
  //if (maxsendbuf == 0 || nmigrate*nbytes > maxsendbuf) { // this doesn't work, not sure why 
    maxsendbuf = nmigrate*nbytes;
    if (maxsendbuf > int(d_sbuf.extent(0)))
      d_sbuf = DAT::t_char_1d("comm:sbuf",maxsendbuf);
  //}

  // fill proclist with procs to send to
  // pack sbuf with particles to migrate
  // if flag == PDISCARD, particle is deleted but not sent
  // change icell of migrated particle to owning cell on receiving proc
  // nsend = particles that actually migrate
  // if no custom attributes, pack particles directly via memcpy()
  // else pack_custom() performs packing into sbuf

  int nsend = 0;
  //int offset = 0;

  h_nsend() = 0;
  k_nsend.modify<SPAHostType>();
  k_nsend.sync<DeviceType>();

  particle_kk->sync(Device,PARTICLE_MASK);
  grid_kk->sync(Device,CELL_MASK);

  d_cells = grid_kk->k_cells.d_view;
  d_particles = particle_kk->k_particles.d_view;

  if (ncustom)
    error->all(FLERR,"Custom per-particles attributes not yet supported with Kokkos");

  //if (!ncustom) {
    copymode = 1;
    if (sparta->kokkos->need_atomics)
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagCommMigrateParticles<1> >(0,nmigrate),*this);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagCommMigrateParticles<0> >(0,nmigrate),*this);
    DeviceType::fence();
    //pack_serial(0,nmigrate);
    copymode = 0;

  //} else {
  //  for (i = 0; i < nmigrate; i++) {
  //    j = plist[i];
  //    if (particles[j].flag == PDISCARD) continue;
  //    pproc[nsend++] = cells[particles[j].icell].proc;
  //    particles[j].icell = cells[particles[j].icell].ilocal;
  //    memcpy(&sbuf[offset],&particles[j],nbytes_particle);
  //    offset += nbytes_particle;
  //    particle->pack_custom(j,&sbuf[offset]);
  //    offset += nbytes_custom;
  //  }
  //}

  particle_kk->modify(Device,PARTICLE_MASK);
  d_particles = t_particle_1d(); // destroy reference to reduce memory use

  k_pproc.modify<DeviceType>();
  k_pproc.sync<SPAHostType>();

  k_nsend.modify<DeviceType>();
  k_nsend.sync<SPAHostType>();
  nsend = h_nsend();

  // compress my list of particles

  particle->compress_migrate(nmigrate,plist);
  int ncompress = particle->nlocal;

  // create or augment irregular communication plan
  // nrecv = # of incoming particles

  IrregularKokkos* iparticle_kk = (IrregularKokkos*) iparticle;

  int nrecv;
  if (neighflag)
    nrecv = iparticle_kk->augment_data_uniform(nsend,pproc);
  else 
    nrecv = iparticle_kk->create_data_uniform(nsend,pproc,commsortflag);

  // extend particle list if necessary

  particle->grow(nrecv);

  // perform irregular communication
  // if no custom attributes, append recv particles directly to particle list
  // else receive into rbuf, unpack particles one by one via unpack_custom()

  //if (!ncustom)
    particle_kk->sync(Device,PARTICLE_MASK);

    d_particles = particle_kk->k_particles.d_view;
    iparticle_kk->
      exchange_uniform(d_sbuf,nbytes,
                       (char *) (d_particles.data()+particle->nlocal));

    particle_kk->modify(Device,PARTICLE_MASK);
    d_particles = t_particle_1d(); // destroy reference to reduce memory use

  //else {
  //  if (nrecv*nbytes > maxrecvbuf) {
  //    maxrecvbuf = nrecv*nbytes;
  //    memory->destroy(rbuf);
  //    memory->create(rbuf,maxrecvbuf,"comm:rbuf");
  //  }
  //
  //  iparticle_kk->exchange_uniform(sbuf,nbytes,rbuf);
  //
  //  offset = 0;
  //  int nlocal = particle->nlocal;
  //  for (i = 0; i < nrecv; i++) {
  //    memcpy(&particle->particles[nlocal],&rbuf[offset],nbytes_particle);
  //    offset += nbytes_particle;
  //    particle->unpack_custom(&rbuf[offset],nlocal);
  //    offset += nbytes_custom;
  //    nlocal++;
  //  }
  //}

  particle->nlocal += nrecv;
  ncomm += nsend;
  return ncompress;
}

template<int NEED_ATOMICS>
KOKKOS_INLINE_FUNCTION
void CommKokkos::operator()(TagCommMigrateParticles<NEED_ATOMICS>, const int &i) const {
  const int j = d_plist[i];
  if (d_particles[j].flag == PDISCARD) return;
  int nsend;
  if (NEED_ATOMICS)
    nsend = Kokkos::atomic_fetch_add(&d_nsend(),1);
  else {
    nsend = d_nsend();
    d_nsend()++;
  }
  d_pproc[nsend] = d_cells[d_particles[j].icell].proc;
  d_particles[j].icell = d_cells[d_particles[j].icell].ilocal;
  const int offset = nsend*nbytes_particle;
  memcpy(&d_sbuf[offset],&d_particles[j],nbytes_particle);
}

inline
void CommKokkos::pack_serial(const int start, const int end) const {
  int nsend = 0;
  for (int i = start; i < end; i++) {
    const int j = d_plist[i];
    if (d_particles[j].flag == PDISCARD) return;
    d_pproc[nsend++] = d_cells[d_particles[j].icell].proc;
    d_particles[j].icell = d_cells[d_particles[j].icell].ilocal;
    const int offset = nsend*nbytes_particle;
    memcpy(&d_sbuf[offset],&d_particles[j],nbytes_particle);
  }
}


/* ----------------------------------------------------------------------
   migrate grid cells with their particles to new procs
   called from BalanceGrid and FixBalance
------------------------------------------------------------------------- */

void CommKokkos::migrate_cells(int nmigrate)
{
  CollideVSSKokkos* collide_kk = (CollideVSSKokkos*) collide;
  if (collide)
    collide_kk->sync(Host,ALL_MASK);

  Comm::migrate_cells(nmigrate);

  if (collide)
    collide_kk->modify(Host,ALL_MASK);
}

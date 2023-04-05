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

CommKokkos::CommKokkos(SPARTA *sparta) : Comm(sparta),
  particle_kk_copy(sparta)
{
  delete iparticle;
  iparticle = new IrregularKokkos(sparta);

  k_nsend = DAT::tdual_int_scalar("comm:nsend");
  d_nsend = k_nsend.d_view;
  h_nsend = k_nsend.h_view;
  d_nlocal = DAT::t_int_scalar("comm:nlocal");
}

/* ---------------------------------------------------------------------- */

CommKokkos::~CommKokkos()
{
  if (copymode) return;

  if (!sparta->kokkos->comm_serial) {
    pproc = NULL;
  }

  particle_kk_copy.uncopy();
}

/* ----------------------------------------------------------------------
   migrate particles to new procs after particle move
   return particle nlocal after compression,
     so Update can iterate on particle move
------------------------------------------------------------------------- */

int CommKokkos::migrate_particles(int nmigrate, int *plist, DAT::t_int_1d &d_plist_in)
{
  GridKokkos* grid_kk = (GridKokkos*) grid;
  ParticleKokkos* particle_kk = ((ParticleKokkos*)particle);
  particle_kk->update_class_variables();
  particle_kk_copy.copy(particle_kk);

  if (sparta->kokkos->comm_serial) {
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
  nbytes_total = nbytes_particle + nbytes_custom;

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
    d_pproc = DAT::t_int_1d(Kokkos::view_alloc("comm:pproc",Kokkos::WithoutInitializing),maxpproc);
    h_pproc = HAT::t_int_1d(Kokkos::view_alloc("comm:pproc_mirror",Kokkos::WithoutInitializing),maxpproc);
    pproc = h_pproc.data();
  }
  //if (maxsendbuf == 0 || nmigrate*nbytes_total > maxsendbuf) { // this doesn't work, not sure why
    int maxsendbuf = nmigrate*nbytes_total;
    if (maxsendbuf > int(d_sbuf.extent(0)))
      d_sbuf = DAT::t_char_1d(Kokkos::view_alloc("comm:sbuf",Kokkos::WithoutInitializing),maxsendbuf);
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
  k_nsend.modify_host();
  k_nsend.sync_device();

  particle_kk->sync(Device,PARTICLE_MASK);
  grid_kk->sync(Device,CELL_MASK);

  d_cells = grid_kk->k_cells.d_view;
  d_particles = particle_kk->k_particles.d_view;

  copymode = 1;
  if (!ncustom) {

    if (sparta->kokkos->need_atomics)
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagCommMigrateParticles<1,0> >(0,nmigrate),*this);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagCommMigrateParticles<0,0> >(0,nmigrate),*this);

  } else {

    if (sparta->kokkos->need_atomics)
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagCommMigrateParticles<1,1> >(0,nmigrate),*this);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagCommMigrateParticles<0,1> >(0,nmigrate),*this);

  }
  DeviceType().fence();
  copymode = 0;

  particle_kk->modify(Device,PARTICLE_MASK);
  d_particles = t_particle_1d(); // destroy reference to reduce memory use

  Kokkos::deep_copy(h_pproc,d_pproc);

  k_nsend.modify_device();
  k_nsend.sync_host();
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

  particle_kk->sync(Device,PARTICLE_MASK);
  d_particles = particle_kk->k_particles.d_view;

  if (sparta->kokkos->gpu_aware_flag && !ncustom) {
    iparticle_kk->
      exchange_uniform(d_sbuf,nbytes_total,
                       (char *) (d_particles.data()+particle->nlocal),d_rbuf);
  } else {

    // allocate exact buffer size to reduce GPU <--> CPU memory transfer

    int maxrecvbuf = nrecv*nbytes_total;
    d_rbuf = DAT::t_char_1d(Kokkos::view_alloc("comm:rbuf",Kokkos::WithoutInitializing),maxrecvbuf);

    Kokkos::deep_copy(d_nlocal,particle->nlocal);
    iparticle_kk->exchange_uniform(d_sbuf,nbytes_total,(char *)d_rbuf.data(),d_rbuf);

    copymode = 1;
    if (!ncustom) {

      if (sparta->kokkos->need_atomics)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagCommMigrateUnpackParticles<1,0> >(0,nrecv),*this);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagCommMigrateUnpackParticles<0,0> >(0,nrecv),*this);
      DeviceType().fence();
      copymode = 0;

    } else {

      if (sparta->kokkos->need_atomics)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagCommMigrateUnpackParticles<1,1> >(0,nrecv),*this);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagCommMigrateUnpackParticles<0,1> >(0,nrecv),*this);
      DeviceType().fence();
      copymode = 0;
    }

  }

  particle_kk->modify(Device,PARTICLE_MASK);
  d_particles = t_particle_1d(); // destroy reference to reduce memory use
  d_plist = decltype(d_plist)();

  particle->nlocal += nrecv;
  ncomm += nsend;
  return ncompress;
}

template<int NEED_ATOMICS, int HAVE_CUSTOM>
KOKKOS_INLINE_FUNCTION
void CommKokkos::operator()(TagCommMigrateParticles<NEED_ATOMICS, HAVE_CUSTOM>, const int &i) const {
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
  const int offset = nsend*nbytes_total;
  memcpy(&d_sbuf[offset],&d_particles[j],nbytes_particle);
  if (HAVE_CUSTOM)
    particle_kk_copy.obj.pack_custom_kokkos(j,(char*)(d_sbuf.data()+offset+nbytes_particle));
}

template<int NEED_ATOMICS, int HAVE_CUSTOM>
KOKKOS_INLINE_FUNCTION
void CommKokkos::operator()(TagCommMigrateUnpackParticles<NEED_ATOMICS,HAVE_CUSTOM>, const int &irecv) const {
  int i;
  if (NEED_ATOMICS)
    i = Kokkos::atomic_fetch_add(&d_nlocal(),1);
  else {
    i = d_nlocal();
    d_nlocal()++;
  }
  const int offset = irecv*nbytes_total;
  memcpy(&d_particles[i],&d_rbuf[offset],nbytes_particle);
  if (HAVE_CUSTOM)
    particle_kk_copy.obj.unpack_custom_kokkos((char*)(d_rbuf.data()+offset+nbytes_particle),i);
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
    collide_kk->modified(Host,ALL_MASK);
}

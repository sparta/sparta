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
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "ctype.h"
#include "particle_kokkos.h"
#include "grid_kokkos.h"
#include "update.h"
#include "comm.h"
#include "mixture.h"
#include "collide.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "memory_kokkos.h"
#include "error.h"
#include "kokkos.h"
#include "sparta_masks.h"

#include <Kokkos_Vector.hpp>

using namespace SPARTA_NS;

enum{NONE,DISCRETE,SMOOTH};            // several files
enum{INT,DOUBLE};                      // several files
enum{COPYPARTICLELIST,FIXEDMEMORY};

#define DELTA 16384
#define DELTACELLCOUNT 1

/* ---------------------------------------------------------------------- */

ParticleKokkos::ParticleKokkos(SPARTA *sparta) : Particle(sparta)
{
  d_resize = DAT::t_int_scalar("particle:resize");
  h_resize = HAT::t_int_scalar("particle:resize_mirror");

  k_reorder_pass = DAT::tdual_int_scalar("particle:reorder_pass");
  d_reorder_pass = k_reorder_pass.d_view;
  h_reorder_pass = k_reorder_pass.h_view;

  sorted_kk = 0;
  maxcellcount = 1;

  k_eivec = tdual_struct_tdual_int_1d_1d("particle:eivec",0);
  k_eiarray = tdual_struct_tdual_int_2d_1d("particle:eiarray",0);
  k_edvec = tdual_struct_tdual_float_1d_1d("particle:edvec",0);
  k_edarray = tdual_struct_tdual_float_2d_1d("particle:edarray",0);
}

/* ---------------------------------------------------------------------- */

ParticleKokkos::~ParticleKokkos()
{
  if (copy || copymode) return;

  particles = NULL;
  species = NULL;


  // deallocate views of views in serial to prevent race condition in profiling tools

  for (int i = 0; i < k_eivec.extent(0); i++)
    k_eivec.h_view(i).k_view = decltype(k_eivec.h_view(i).k_view)();

  for (int i = 0; i < k_eiarray.extent(0); i++)
    k_eiarray.h_view(i).k_view = decltype(k_eiarray.h_view(i).k_view)();

  for (int i = 0; i < k_edvec.extent(0); i++)
    k_edvec.h_view(i).k_view = decltype(k_edvec.h_view(i).k_view)();

  for (int i = 0; i < k_edarray.extent(0); i++)
    k_edarray.h_view(i).k_view = decltype(k_edarray.h_view(i).k_view)();

  ewhich = NULL;
  eicol = NULL;
  edcol = NULL;
}

#ifndef SPARTA_KOKKOS_EXACT
/* ----------------------------------------------------------------------
   compress particle list to remove particles with indices in dellist
   dellist indices can be in ANY order
------------------------------------------------------------------------- */

void ParticleKokkos::compress_migrate(int ndelete, int *dellist)
{
  // reallocate next list as needed

  if (maxsort < maxlocal) {
    maxsort = maxlocal;
    memory->destroy(next);
    memory->create(next,maxsort,"particle:next");
  }

  int i;

  nbytes = sizeof(OnePart);

  if (ndelete > d_lists.extent(1)) {
    d_lists = DAT::t_int_2d_lr(Kokkos::view_alloc("particle:lists",Kokkos::WithoutInitializing),2,ndelete);
    d_mlist = Kokkos::subview(d_lists,0,Kokkos::ALL);
    d_slist = Kokkos::subview(d_lists,1,Kokkos::ALL);

    h_lists = HAT::t_int_2d_lr(Kokkos::view_alloc("particle:lists_mirror",Kokkos::WithoutInitializing),2,ndelete);
    h_mlist = Kokkos::subview(h_lists,0,Kokkos::ALL);
    h_slist = Kokkos::subview(h_lists,1,Kokkos::ALL);
  }

  // use next as a scratch vector
  // next is only used for upper locs from nlocal-ndelete to nlocal
  // next[i] = 0 if deleted particle, 1 otherwise

  int upper = nlocal - ndelete;
  for (i = upper; i < nlocal; i++) next[i] = 1;

  for (int m = 0; m < ndelete; m++) {
    i = dellist[m];
    if (i >= upper)
      next[i] = 0;
  }

  int ncopy = 0;
  int ncount = 0;
  for (int j = upper; j < nlocal; j++) {
    if (!next[j]) continue;

    int i = dellist[ncount];
    while (i >= upper) {
      ncount++;
      i = dellist[ncount];
    }
    h_mlist[ncopy] = i;
    h_slist[ncopy] = j;
    ncopy++;
    ncount++;
  }

  nlocal = upper;

  Kokkos::deep_copy(d_lists,h_lists);

  this->sync(Device,PARTICLE_MASK|CUSTOM_MASK);
  d_particles = k_particles.d_view;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagParticleCompressReactions>(0,ncopy),*this);
  copymode = 0;

  this->modify(Device,PARTICLE_MASK|CUSTOM_MASK);
  d_particles = t_particle_1d(); // destroy reference to reduce memory use

  sorted = 0;
  sorted_kk = 0;
}
#endif

KOKKOS_INLINE_FUNCTION
void ParticleKokkos::operator()(TagParticleCompressReactions, const int &i) const {
  const int j = d_mlist[i];
  const int k = d_slist[i];
  d_particles[j] = d_particles[k];
  copy_custom_kokkos(j,k);
}

/* ----------------------------------------------------------------------
   sort particles into grid cells
   set cinfo.first = index of first particle in cell
   set cinfo.count = # of particles in cell
   next[] = index of next particle in same cell, -1 for no more
------------------------------------------------------------------------- */

void ParticleKokkos::sort_kokkos()
{
  sorted_kk = 1;
  int reorder_scheme = COPYPARTICLELIST;

  // FIXEDMEMORY reorder temporarily disabled due to bug on GPUs

  //if (update->have_mem_limit())
  //  reorder_scheme = FIXEDMEMORY;

  const int reorder_flag = (update->reorder_period &&
      (update->ntimestep % update->reorder_period == 0));

  ngrid = grid->nlocal;
  GridKokkos* grid_kk = (GridKokkos*)grid;
  d_cellcount = grid_kk->d_cellcount;
  d_plist = grid_kk->d_plist;

  if (ngrid > int(d_cellcount.extent(0))) {
    d_cellcount = decltype(d_cellcount)();
    MemKK::realloc_kokkos(grid_kk->d_cellcount,"particle:cellcount",ngrid);
    d_cellcount = grid_kk->d_cellcount;
  }

  Kokkos::deep_copy(d_cellcount,0);

  if (ngrid > int(d_plist.extent(0)) || maxcellcount > int(d_plist.extent(1))) {
    d_plist = decltype(d_plist)();
    MemKK::realloc_kokkos(grid_kk->d_plist,"particle:plist",ngrid,maxcellcount);
    d_plist = grid_kk->d_plist;
  }

  this->sync(Device,PARTICLE_MASK);
  d_particles = k_particles.d_view;

  if (reorder_flag && reorder_scheme == COPYPARTICLELIST) {
    if (d_particles.extent(0) > d_offsets_part.extent(0)) {
      MemKK::realloc_kokkos(d_offsets_part,"particle:offsets_part",d_particles.extent(0));
    }
  }

  // icell = global cell the particle is in

  // Cannot grow a Kokkos view in a parallel loop, so
  //  if the capacity of the list is exceeded, count the size
  //  needed, reallocate on the host, and then
  //  repeat the parallel loop again

  int resize = 1;
  while (resize) {
    resize = 0;

    copymode = 1;
    if (sparta->kokkos->need_atomics) {
      if (reorder_flag && reorder_scheme == COPYPARTICLELIST)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagParticleSort<1,1> >(0,nlocal),*this);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagParticleSort<1,0> >(0,nlocal),*this);
    } else {
      if (reorder_flag && reorder_scheme == COPYPARTICLELIST)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagParticleSort<0,1> >(0,nlocal),*this);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagParticleSort<0,0> >(0,nlocal),*this);
    }
    copymode = 0;

    Kokkos::deep_copy(h_resize,d_resize);
    resize = h_resize();

    if (resize) {
      Kokkos::deep_copy(d_cellcount,0);
      int old = maxcellcount;
      maxcellcount = MAX(maxcellcount+MAX(DELTACELLCOUNT,maxcellcount*0.1),resize);
      d_plist = decltype(d_plist)();
      MemKK::realloc_kokkos(grid_kk->d_plist,"particle:plist",ngrid,maxcellcount);
      d_plist = grid_kk->d_plist;

      Kokkos::deep_copy(d_resize,0);
    }
  }

  if (reorder_flag) {

    if (reorder_scheme == COPYPARTICLELIST) {
      if (d_particles.extent(0) > d_sorted.extent(0))
        MemKK::realloc_kokkos(d_sorted,"particle:sorted",d_particles.extent(0));

      if (d_particles.extent(0) > d_sorted_id.extent(0))
        MemKK::realloc_kokkos(d_sorted_id,"particle:sorted_id",d_particles.extent(0));
    } else if (reorder_scheme == FIXEDMEMORY && d_pswap1.size() == 0) {
      nParticlesWksp = MIN(nlocal,(double)update->global_mem_limit/sizeof(Particle::OnePart));
      d_pswap1 = t_particle_1d(Kokkos::view_alloc("particle:swap1",Kokkos::WithoutInitializing),nParticlesWksp);
      d_pswap2 = t_particle_1d(Kokkos::view_alloc("particle:swap2",Kokkos::WithoutInitializing),nParticlesWksp);
    }

    nbytes = sizeof(OnePart);

    if (reorder_scheme == COPYPARTICLELIST) {
      copymode = 1;
      Kokkos::parallel_scan(Kokkos::RangePolicy<DeviceType, TagParticleReorder_COPYPARTICLELIST1>(0,ngrid),*this);
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagParticleReorder_COPYPARTICLELIST2>(0,nlocal),*this);
      copymode = 0;
      //auto tmp = k_particles.d_view;
      //k_particles.d_view = d_sorted;
      //d_particles = k_particles.d_view;
      //d_sorted = tmp;
      Kokkos::deep_copy(d_particles,d_sorted);

      this->modify(Device,PARTICLE_MASK);
    }
    else if (reorder_scheme == FIXEDMEMORY) {
      // Copy particle destinations into the particle list cell locations
      //  (to avoid adding a "destination" integer in OnePart for the fixed memory reorder)
      // After the particle list has been reordered, reset the icell values to correctly reflect
      // the variable naming.
      copymode = 1;
      Kokkos::parallel_scan(Kokkos::RangePolicy<DeviceType, TagCopyParticleReorderDestinations>(0,ngrid),*this);
      copymode = 0;

      int npasses = (nlocal-1)/nParticlesWksp + 1;
      for (int ipass=0; ipass < npasses; ++ipass) {

        h_reorder_pass() = ipass;
        k_reorder_pass.modify_host();
        k_reorder_pass.sync_device();

        // identify next set of particles to reorder
        copymode = 1;
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixedMemoryReorderInit>(0,nParticlesWksp),*this);
        copymode = 0;

        // reorder this set of particles
        copymode = 1;
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixedMemoryReorder>(0,nParticlesWksp),*this);
        copymode = 0;
      }

      // reset the icell values in the particle list
      copymode = 1;
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagSetIcellFromPlist>(0,ngrid),*this);
      copymode = 0;
      this->modify(Device,PARTICLE_MASK);

      // destroy references to reduce memory use
      d_pswap1 = t_particle_1d();
      d_pswap2 = t_particle_1d();
    }
  }

  d_particles = t_particle_1d(); // destroy reference to reduce memory use
  d_plist = decltype(d_plist)();
}

KOKKOS_INLINE_FUNCTION
void ParticleKokkos::operator()(TagCopyParticleReorderDestinations, const int icell, int &m_fill, const bool &final) const
{
  // load icell values with reorder destination index for particle reordering
  //  (location in d_particles where particle is moving)
  for (int j = 0; j < d_cellcount[icell]; j++) {
    if (final) {
      const int iparticle = d_plist(icell,j);
      d_particles[iparticle].icell = m_fill;
      d_plist(icell,j) = m_fill; // the new plist after reordering
    }
    m_fill++;
  }
}

KOKKOS_INLINE_FUNCTION
void ParticleKokkos::operator()(TagFixedMemoryReorderInit, const int &i) const
{
  // note:  "i" is a thread id (cannot be greater than number of particles allocated in d_pswap* workspaces)

  // assign batch of threads to next batch of particles (batch size = nParticlesWksp)
  d_pswap1[i].icell = -999; // default to thread isn't moving a particle
  int nextParticleToCheckForReordering = d_reorder_pass() * nParticlesWksp;
  int n = nextParticleToCheckForReordering + i;
  if (n < nlocal) {
    if (d_particles[n].icell != n) {
      d_pswap1[i] = d_particles[n]; // copy the moving particle into the work space so the thread can move it later.
      d_particles[n].icell = -1;    // current location of moving particle is marked as vacant
    }
  }
}

KOKKOS_INLINE_FUNCTION
void ParticleKokkos::operator()(TagFixedMemoryReorder, const int &i) const
{
  // particle movement for this thread continues until a particle is moved to a vacant location (indicated by icell = -1)
  OnePart *movePtr;
  int newParticleLoc;
  bool iHaveAnotherParticle = false;
  if (d_pswap1[i].icell != -999)
    iHaveAnotherParticle = true;

  int count = 0;
  while (iHaveAnotherParticle){
    if (count % 2 == 0) { // even
      movePtr = &d_pswap1[i];
      newParticleLoc = movePtr->icell; // note:  this is a particle location, NOT a cell location
      d_pswap2[i] = d_particles[newParticleLoc];
      if (d_pswap2[i].icell == -1)
        iHaveAnotherParticle = false;
    }
    else { // odd
      movePtr = &d_pswap2[i];
      newParticleLoc = movePtr->icell; // note:  this is a particle location, NOT a cell location
      d_pswap1[i] = d_particles[newParticleLoc];
      if (d_pswap1[i].icell == -1)
        iHaveAnotherParticle = false;
    }
    d_particles[newParticleLoc] = *movePtr;
    count++;
  }
}

KOKKOS_INLINE_FUNCTION
void ParticleKokkos::operator()(TagSetIcellFromPlist, const int &icell) const
{
  for (int j = 0; j < d_cellcount[icell]; j++) {
    const int iparticle = d_plist(icell,j);
    d_particles[iparticle].icell = icell;
  }
}

template<int NEED_ATOMICS, int REORDER_FLAG>
KOKKOS_INLINE_FUNCTION
void ParticleKokkos::operator()(TagParticleSort<NEED_ATOMICS,REORDER_FLAG>, const int &i) const
{
  const int icell = d_particles[i].icell;
  int j;
  if (NEED_ATOMICS)
    j = Kokkos::atomic_fetch_add(&d_cellcount[icell],1);
  else {
    j = d_cellcount[icell];
    d_cellcount[icell]++;
  }

  if (j >= maxcellcount)
    d_resize() = MAX(d_resize(),j+1);
  else {
    d_plist(icell,j) = i;

    if (REORDER_FLAG)
      d_offsets_part[i] = j;
  }
}

KOKKOS_INLINE_FUNCTION
void ParticleKokkos::operator()(TagParticleReorder_COPYPARTICLELIST1, const int icell, int &m_fill, const bool &final) const
{
  if (final) {
    for (int j = 0; j < d_cellcount[icell]; j++) {
      const int iparticle = d_plist(icell,j);
      d_sorted_id[m_fill++] = iparticle;
    }
  } else
    m_fill += d_cellcount[icell];
}

KOKKOS_INLINE_FUNCTION
void ParticleKokkos::operator()(TagParticleReorder_COPYPARTICLELIST2, const int offset) const
{
  const int iparticle = d_sorted_id[offset];
  const Particle::OnePart &particle_i = d_particles[iparticle];
  d_sorted[offset] = particle_i;
  const int icell = particle_i.icell;
  const int j = d_offsets_part[iparticle];
  d_plist(icell,j) = offset;
}

/* ----------------------------------------------------------------------
   set the initial weight of each particle
   called by Update before particle move
   only called if particle weighting is enabled
   only grid-based weighting is currently implemented
------------------------------------------------------------------------- */

void ParticleKokkos::pre_weight()
{
  auto grid_kk = dynamic_cast<GridKokkos*>(grid);
  auto& k_cinfo = grid_kk->k_cinfo;
  grid_kk->sync(Device,CINFO_MASK);
  this->sync(Device,PARTICLE_MASK);
  auto d_cinfo = k_cinfo.d_view;
  auto d_particles = k_particles.d_view;

  Kokkos::parallel_for(nlocal, KOKKOS_LAMBDA(int i) {
    auto icell = d_particles[i].icell;
    d_particles[i].weight = d_cinfo[icell].weight;
  });
  this->modify(Device,PARTICLE_MASK);
  d_particles = t_particle_1d();
}

/* ----------------------------------------------------------------------
   clone/delete each particle based on ratio of its initial/final weights
   called by Update after particle move and migration
   only called if particle weighting is enabled
   only grid-based weighting is currently implemented
------------------------------------------------------------------------- */

struct PostWeightPair { int i; int id; };

void ParticleKokkos::post_weight()
{
  if (ncustom) {
    error->all(FLERR,"Custom per-particles attributes not yet supported with Kokkos");
  }
  constexpr int METHOD = 1;
  if (METHOD == 1) { // just call the host one
    this->sync(Host,PARTICLE_MASK);
    int prev_auto_sync = sparta->kokkos->auto_sync;
    sparta->kokkos->auto_sync = 1;
    Particle::post_weight();
    sparta->kokkos->auto_sync = prev_auto_sync;
    this->modify(Host,PARTICLE_MASK);
  } else if (METHOD == 2) { // Kokkos-parallel, gives same (correct) answer
    Kokkos::View<double*> d_ratios("post_weight:ratios", nlocal);
    auto h_ratios = Kokkos::create_mirror_view(d_ratios);
    auto grid_kk = dynamic_cast<GridKokkos*>(grid);
    auto& k_cinfo = grid_kk->k_cinfo;
    grid_kk->sync(Device,CINFO_MASK);
    this->sync(Device,PARTICLE_MASK);
    auto d_particles = k_particles.d_view;
    auto d_cinfo = k_cinfo.d_view;
    Kokkos::vector<PostWeightPair> map;
    map.set_overallocation(0.5);
    map.resize(nlocal);
    auto d_map = map.d_view;
    Kokkos::parallel_for(nlocal, KOKKOS_LAMBDA(int i) {
      auto icell = d_particles[i].icell;
      d_ratios[i] = d_particles[i].weight / d_cinfo[icell].weight;
      d_map[i].id = d_particles[i].id;
      d_map[i].i = i;
    });
    Kokkos::deep_copy(h_ratios, d_ratios);
    map.device_to_host();
    int nlocal_original = nlocal;
    int i = 0;
    while (i < nlocal_original) {
      auto ratio = h_ratios[map[i].i];
      if (ratio == 1.0) {
        i++;
      } else if (ratio < 1.0) {
        if (wrandom->uniform() > ratio) {
          map[i] = map.back();
          if (map.size() > nlocal_original) ++i;
          else --nlocal_original;
          map.pop_back();
        } else {
          ++i;
        }
      } else if (ratio > 1.0) {
        int nclone = int(ratio);
        double fraction = ratio - nclone;
        --nclone;
        if (wrandom->uniform() < fraction) ++nclone;
        for (int m = 0; m < nclone; ++m) {
          map.push_back(map[i]);
          map.back().id = MAXSMALLINT * wrandom->uniform();
        }
        ++i;
      }
    }
    map.host_to_device();
    nlocal = (int(map.size()));
    grow(0);
    Kokkos::View<OnePart*> d_newparticles("post_weight:newparticles", maxlocal);
    d_map = map.d_view;
    Kokkos::parallel_for(nlocal, KOKKOS_LAMBDA(int i) {
      d_newparticles[i] = d_particles[d_map[i].i];
      d_newparticles[i].id = d_map[i].id;
    });
    Kokkos::deep_copy(k_particles.d_view,d_newparticles);
    this->modify(Device,PARTICLE_MASK);
    d_particles = t_particle_1d();
  }
}

/* ---------------------------------------------------------------------- */

void ParticleKokkos::update_class_variables() {
  d_species = k_species.d_view;
  this->sync(Device,SPECIES_MASK);

  boltz = update->boltz;
  collide_rot = 0;
  vibstyle = NONE;
  if (collide) {
    vibstyle = collide->vibstyle;
    if (collide->rotstyle != NONE) collide_rot = 1;
  }
}

/* ----------------------------------------------------------------------
   insure particle list can hold nextra new particles
   if defined, also grow custom particle arrays and initialize with zeroes
------------------------------------------------------------------------- */

void ParticleKokkos::grow(int nextra)
{
  bigint target = (bigint) nlocal + nextra;
  if (target <= maxlocal) return;

  bigint newmax = maxlocal;
  while (newmax < target) newmax += MAX(DELTA, newmax*0.1);
  int oldmax = maxlocal;

  if (newmax > MAXSMALLINT)
    error->one(FLERR,"Per-processor particle count is too big");

  maxlocal = newmax;
  if (particles == NULL)
    k_particles = tdual_particle_1d("particle:particles",maxlocal);
  else {
    this->sync(Device,PARTICLE_MASK); // force resize on device
    k_particles.resize(maxlocal);
    this->modify(Device,PARTICLE_MASK); // needed for auto sync
  }
  d_particles = k_particles.d_view;
  particles = k_particles.h_view.data();

  if (ncustom == 0) return;

  for (int i = 0; i < ncustom; i++) {
    if (ename[i] == NULL) continue;
    grow_custom(i,oldmax,maxlocal);
  }
}

/* ----------------------------------------------------------------------
   insure species list can hold maxspecies species
   assumes that maxspecies has already been increased
------------------------------------------------------------------------- */

void ParticleKokkos::grow_species()
{
  if (sparta->kokkos->prewrap) {
    Particle::grow_species();
  } else {
    if (species == NULL)
      k_species = tdual_species_1d("particle:species",maxspecies);
    else
      k_species.resize(maxspecies);
    species = k_species.h_view.data();
  }
}

/* ---------------------------------------------------------------------- */

void ParticleKokkos::wrap_kokkos()
{
  // species

  if (species != k_species.h_view.data()) {
    memoryKK->wrap_kokkos(k_species,species,nspecies,"particle:species");
    k_species.modify_host();
    k_species.sync_device();
    memory->sfree(species);
    species = k_species.h_view.data();
  }

  // mixtures

  k_species2group = DAT::tdual_int_2d("particle:species2group",nmixture,nspecies);
  for (int i = 0; i < nmixture; i++)
    for (int j = 0; j < nspecies; j++)
      k_species2group.h_view(i,j) = mixture[i]->species2group[j];
  k_species2group.modify_host();
  k_species2group.sync_device();

  //if (mixtures != k_mixtures.h_view.data()) {
  //  memoryKK->wrap_kokkos(k_mixtures,mixture,nmixture,"particle:mixture");
  //  k_mixtures.modify_host();
  //  k_mixtures.sync_device();
  //  memory->sfree(mixtures);
  //  mixtures = k_mixtures.h_view.data();
  //}
}

/* ----------------------------------------------------------------------
   add a custom attribute with name
   assumes name does not already exist, except in case of restart
   type = 0/1 for int/double
   size = 0 for vector, size > 0 for array with size columns
   allocate the vector or array to current maxlocal via grow_custom()
   return index of its location;
------------------------------------------------------------------------- */

int ParticleKokkos::add_custom(char *name, int type, int size)
{
  ///modifies eivec,eiarray,edvec,edarray on either host or device, probably device since host isn't modified. May just want to use host
  ///modifies ewhich on host, sync to device here since it is never modified on the device
  //

  // force resize on host

  k_eivec.modify_host();
  k_eiarray.modify_host();
  k_edvec.modify_host();
  k_edarray.modify_host();

  k_ewhich.modify_host();
  k_eicol.modify_host();
  k_edcol.modify_host();

  int index;

  // if name already exists
  // just return index if a restart script and re-defining the name
  // else error

  index = find_custom(name);
  if (index >= 0) {
    if (custom_restart_flag == NULL || custom_restart_flag[index] == 1)
      error->all(FLERR,"Custom particle attribute name already exists");
    custom_restart_flag[index] = 1;
    return index;
  }

  // use first available NULL entry or allocate a new one

  for (index = 0; index < ncustom; index++)
    if (ename[index] == NULL) break;

  if (index == ncustom) {
    ncustom++;
    ename = (char **) memory->srealloc(ename,ncustom*sizeof(char *),
                                       "particle:ename");
    memory->grow(etype,ncustom,"particle:etype");
    memory->grow(esize,ncustom,"particle:esize");
    memoryKK->grow_kokkos(k_ewhich,ewhich,ncustom,"particle:ewhich");
  }

  int n = strlen(name) + 1;
  ename[index] = new char[n];
  strcpy(ename[index],name);
  etype[index] = type;
  esize[index] = size;

  if (type == INT) {
    if (size == 0) {
      ewhich[index] = ncustom_ivec++;
      eivec = (int **)
        memory->srealloc(eivec,ncustom_ivec*sizeof(int *),"particle:eivec");
      eivec[ncustom_ivec-1] = NULL;
      k_eivec.resize(ncustom_ivec);
      memory->grow(icustom_ivec,ncustom_ivec,"particle:icustom_ivec");
      icustom_ivec[ncustom_ivec-1] = index;
    } else {
      ewhich[index] = ncustom_iarray++;
      eiarray = (int ***)
        memory->srealloc(eiarray,ncustom_iarray*sizeof(int **),
                         "particle:eiarray");
      eiarray[ncustom_iarray-1] = NULL;
      k_eiarray.resize(ncustom_iarray);
      memory->grow(icustom_iarray,ncustom_iarray,"particle:icustom_iarray");
      icustom_iarray[ncustom_iarray-1] = index;
      memoryKK->grow_kokkos(k_eicol,eicol,ncustom_iarray,"particle:eicol");
      eicol[ncustom_iarray-1] = size;
    }
  } else if (type == DOUBLE) {
    if (size == 0) {
      ewhich[index] = ncustom_dvec++;
      edvec = (double **)
        memory->srealloc(edvec,ncustom_dvec*sizeof(double *),"particle:edvec");
      edvec[ncustom_dvec-1] = NULL;
      k_edvec.resize(ncustom_dvec);
      memory->grow(icustom_dvec,ncustom_dvec,"particle:icustom_dvec");
      icustom_dvec[ncustom_dvec-1] = index;
    } else {
      ewhich[index] = ncustom_darray++;
      edarray = (double ***)
        memory->srealloc(edarray,ncustom_darray*sizeof(double **),
                         "particle:edarray");
      edarray[ncustom_darray-1] = NULL;
      k_edarray.resize(ncustom_darray);
      memory->grow(icustom_darray,ncustom_darray,"particle:icustom_darray");
      icustom_darray[ncustom_darray-1] = index;
      memoryKK->grow_kokkos(k_edcol,edcol,ncustom_darray,"particle:edcol");
      edcol[ncustom_darray-1] = size;
    }
  }

  // ewhich,eicol,edcol never modified on the device, so sync here

  k_ewhich.modify_host();
  k_ewhich.sync_device();

  k_eicol.modify_host();
  k_eicol.sync_device();

  k_edcol.modify_host();
  k_edcol.sync_device();

  grow_custom(index,0,maxlocal);

  return index;
}

/* ----------------------------------------------------------------------
   grow the vector/array associated with custom attribute with index
   nold = old length, nnew = new length (typically maxlocal)
   set new values to 0 via memset()
------------------------------------------------------------------------- */

void ParticleKokkos::grow_custom(int index, int nold, int nnew)
{
  // modifies the inner part of eivec,eiarray,edvec,edarray on whatever, and the outer view on the host

  k_eivec.sync_host();
  k_eiarray.sync_host();
  k_edvec.sync_host();
  k_edarray.sync_host();

  if (etype[index] == INT) {
    if (esize[index] == 0) {
      int *ivector = eivec[ewhich[index]];
      auto k_ivector = k_eivec.h_view[ewhich[index]].k_view;
      k_ivector.modify_host(); // force resize on host
      memoryKK->grow_kokkos(k_ivector,ivector,nold+nnew,"particle:eivec");
      k_eivec.h_view[ewhich[index]].k_view = k_ivector;
      eivec[ewhich[index]] = ivector;
    } else {
      int **iarray = eiarray[ewhich[index]];
      auto k_iarray = k_eiarray.h_view[ewhich[index]].k_view;
      k_iarray.modify_host(); // force resize on host
      memoryKK->grow_kokkos(k_iarray,iarray,nold+nnew,esize[index],"particle:eiarray");
      k_eiarray.h_view[ewhich[index]].k_view = k_iarray;
      eiarray[ewhich[index]] = iarray;
    }

  } else {
    if (esize[index] == 0) {
      double *dvector = edvec[ewhich[index]];
      auto k_dvector = k_edvec.h_view[ewhich[index]].k_view;
      k_dvector.modify_host(); // force resize on host
      memoryKK->grow_kokkos(k_dvector,dvector,nold+nnew,"particle:edvec");
      k_edvec.h_view[ewhich[index]].k_view = k_dvector;
      edvec[ewhich[index]] = dvector;
    } else {
      double **darray = edarray[ewhich[index]];
      auto k_darray = k_edarray.h_view[ewhich[index]].k_view;
      k_darray.modify_host(); // force resize on host
      memoryKK->grow_kokkos(k_darray,darray,nold+nnew,esize[index],"particle:edarray");
      k_edarray.h_view[ewhich[index]].k_view = k_darray;
      edarray[ewhich[index]] = darray;
    }
  }

  k_eivec.modify_host();
  k_eiarray.modify_host();
  k_edvec.modify_host();
  k_edarray.modify_host();

  k_eivec.sync_device();
  k_eiarray.sync_device();
  k_edvec.sync_device();
  k_edarray.sync_device();
}

/* ----------------------------------------------------------------------
   remove a custom attribute at location index
   free memory for name and vector/array and set ptrs to NULL
   ncustom lists never shrink, but indices stored between
     the ncustom list and the dense vector/array lists must be reset
------------------------------------------------------------------------- */

void ParticleKokkos::remove_custom(int index)
{
  // modifies the outer host view, deletes the inner dual view
  //
  delete [] ename[index];
  ename[index] = NULL;

  if (etype[index] == INT) {
    if (esize[index] == 0) {
      ncustom_ivec--;
      for (int i = ewhich[index]; i < ncustom_ivec; i++) {
        icustom_ivec[i] = icustom_ivec[i+1];
        ewhich[icustom_ivec[i]] = i;
        eivec[i] = eivec[i+1];
        k_eivec.h_view[i] = k_eivec.h_view[i+1];
      }
    } else {
      ncustom_iarray--;
      for (int i = ewhich[index]; i < ncustom_iarray; i++) {
        icustom_iarray[i] = icustom_iarray[i+1];
        ewhich[icustom_iarray[i]] = i;
        eiarray[i] = eiarray[i+1];
        eicol[i] = eicol[i+1];
        k_eiarray.h_view[i] = k_eiarray.h_view[i+1];
      }
    }
  } else if (etype[index] == DOUBLE) {
    if (esize[index] == 0) {
      ncustom_dvec--;
      for (int i = ewhich[index]; i < ncustom_dvec; i++) {
        icustom_dvec[i] = icustom_dvec[i+1];
        ewhich[icustom_dvec[i]] = i;
        edvec[i] = edvec[i+1];
        k_edvec.h_view[i] = k_edvec.h_view[i+1];
      }
      k_edvec.modify_host();
    } else {
      ncustom_darray--;
      for (int i = ewhich[index]; i < ncustom_darray; i++) {
        icustom_darray[i] = icustom_darray[i+1];
        ewhich[icustom_darray[i]] = i;
        edarray[i] = edarray[i+1];
        edcol[i] = edcol[i+1];
        k_edarray.h_view[i] = k_edarray.h_view[i+1];
      }
      k_edarray.modify_host();
    }
  }

  // set ncustom = 0 if custom list is now entirely empty

  int empty = 1;
  for (int i = 0; i < ncustom; i++)
    if (ename[i]) empty = 0;
  if (empty) ncustom = 0;

  k_eivec.sync_device();
  k_eiarray.sync_device();
  k_edvec.sync_device();
  k_edarray.sync_device();
}

/* ----------------------------------------------------------------------
   copy info for one particle in custom attribute vectors/arrays
   into location I from location J
------------------------------------------------------------------------- */

void ParticleKokkos::copy_custom(int i, int j)
{
  this->sync(Host,CUSTOM_MASK);

  int m;

  // caller does not always check this
  // shouldn't be a problem, but valgrind can complain if memcpy to self
  // oddly memcpy(&particles[i],&particles[j],sizeof(OnePart)) seems OK

  if (i == j) return;

  // 4 flavors of vectors/arrays

  if (ncustom_ivec) {
    for (m = 0; m < ncustom_ivec; m++) eivec[m][i] = eivec[m][j];
  }
  if (ncustom_iarray) {
    for (m = 0; m < ncustom_iarray; m++)
      memcpy(eiarray[m][i],eiarray[m][j],eicol[m]*sizeof(int));
  }
  if (ncustom_dvec) {
    for (m = 0; m < ncustom_dvec; m++) edvec[m][i] = edvec[m][j];
  }
  if (ncustom_darray) {
    for (m = 0; m < ncustom_darray; m++)
      memcpy(edarray[m][i],edarray[m][j],edcol[m]*sizeof(double));
  }

  this->modify(Host,CUSTOM_MASK);
}

/* ----------------------------------------------------------------------
   copy info for one particle in custom attribute vectors/arrays
   into location I from location J
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ParticleKokkos::copy_custom_kokkos(int i, int j) const
{
  int m,ncol;

  // caller does not always check this
  // shouldn't be a problem, but valgrind can complain if memcpy to self
  // oddly memcpy(&particles[i],&particles[j],sizeof(OnePart)) seems OK

  if (i == j) return;

  // 4 flavors of vectors/arrays

  if (ncustom_ivec) {
    for (m = 0; m < ncustom_ivec; m++)
      k_eivec.d_view[m].k_view.d_view[i] = k_eivec.d_view[m].k_view.d_view[j];
  }
  if (ncustom_iarray) {
    for (m = 0; m < ncustom_iarray; m++)
      for (ncol = 0; ncol < k_eicol.d_view[m]; ncol++)
        k_eiarray.d_view[m].k_view.d_view(i,ncol) = k_eiarray.d_view[m].k_view.d_view(j,ncol);
  }
  if (ncustom_dvec) {
    for (m = 0; m < ncustom_dvec; m++)
      k_edvec.d_view[m].k_view.d_view[i] = k_edvec.d_view[m].k_view.d_view[j];
  }
  if (ncustom_darray) {
    for (m = 0; m < ncustom_darray; m++)
      for (ncol = 0; ncol < k_edcol.d_view[m]; ncol++)
        k_edarray.d_view[m].k_view.d_view(i,ncol) = k_edarray.d_view[m].k_view.d_view(j,ncol);
  }
}

/* ----------------------------------------------------------------------
   pack a custom attributes for a single particle N into buf
   this is done in order of 4 styles of vectors/arrays, not in ncustom order
------------------------------------------------------------------------- */

void ParticleKokkos::pack_custom(int n, char *buf)
{
  this->sync(Host,CUSTOM_MASK);
  Particle::pack_custom(n,buf);
}

/* ----------------------------------------------------------------------
   unpack custom attributes for a single particle N from buf
   this is done in order of 4 styles of vectors/arrays, not in ncustom order
------------------------------------------------------------------------- */

void ParticleKokkos::unpack_custom(char *buf, int n)
{
  Particle::unpack_custom(buf,n);
  this->modify(Host,CUSTOM_MASK);
}

/* ---------------------------------------------------------------------- */

void ParticleKokkos::sync(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if (sparta->kokkos->auto_sync)
      modify(Host,mask);
    if (mask & PARTICLE_MASK) k_particles.sync_device();
    if (mask & SPECIES_MASK) k_species.sync_device();
    if (mask & CUSTOM_MASK) {
      if (ncustom) {
        if (ncustom_ivec)
          for (int i = 0; i < ncustom_ivec; i++)
            k_eivec.h_view[i].k_view.sync_device();

        if (ncustom_iarray)
          for (int i = 0; i < ncustom_iarray; i++)
            k_eiarray.h_view[i].k_view.sync_device();

        if (ncustom_dvec)
          for (int i = 0; i < ncustom_dvec; i++)
            k_edvec.h_view[i].k_view.sync_device();

        if (ncustom_darray)
          for (int i = 0; i < ncustom_darray; i++)
            k_edarray.h_view[i].k_view.sync_device();
      }
    }
  } else {
    if (mask & PARTICLE_MASK) k_particles.sync_host();
    if (mask & SPECIES_MASK) k_species.sync_host();
    if (mask & CUSTOM_MASK) {
      if (ncustom_ivec)
        for (int i = 0; i < ncustom_ivec; i++)
          k_eivec.h_view[i].k_view.sync_host();

      if (ncustom_iarray)
        for (int i = 0; i < ncustom_iarray; i++)
          k_eiarray.h_view[i].k_view.sync_host();

      if (ncustom_dvec)
        for (int i = 0; i < ncustom_dvec; i++)
          k_edvec.h_view[i].k_view.sync_host();

      if (ncustom_darray)
        for (int i = 0; i < ncustom_darray; i++)
          k_edarray.h_view[i].k_view.sync_host();
    }
  }
}

/* ---------------------------------------------------------------------- */

void ParticleKokkos::modify(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if (mask & PARTICLE_MASK) k_particles.modify_device();
    if (mask & SPECIES_MASK) k_species.modify_device();
    if (mask & CUSTOM_MASK) {
      if (ncustom) {
        if (ncustom_ivec)
          for (int i = 0; i < ncustom_ivec; i++)
            k_eivec.h_view[i].k_view.modify_device();

        if (ncustom_iarray)
          for (int i = 0; i < ncustom_iarray; i++)
            k_eiarray.h_view[i].k_view.modify_device();

        if (ncustom_dvec)
          for (int i = 0; i < ncustom_dvec; i++)
            k_edvec.h_view[i].k_view.modify_device();

        if (ncustom_darray)
          for (int i = 0; i < ncustom_darray; i++)
            k_edarray.h_view[i].k_view.modify_device();
      }
    }
    if (sparta->kokkos->auto_sync)
      sync(Host,mask);
  } else {
    if (mask & PARTICLE_MASK) k_particles.modify_host();
    if (mask & SPECIES_MASK) k_species.modify_host();
    if (mask & CUSTOM_MASK) {
      if (ncustom) {
        if (ncustom_ivec)
          for (int i = 0; i < ncustom_ivec; i++)
            k_eivec.h_view[i].k_view.modify_host();

        if (ncustom_iarray)
          for (int i = 0; i < ncustom_iarray; i++)
            k_eiarray.h_view[i].k_view.modify_host();

        if (ncustom_dvec)
          for (int i = 0; i < ncustom_dvec; i++)
            k_edvec.h_view[i].k_view.modify_host();

        if (ncustom_darray)
          for (int i = 0; i < ncustom_darray; i++)
            k_edarray.h_view[i].k_view.modify_host();
      }
    }
  }
}

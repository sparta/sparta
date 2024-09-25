/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
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

//#include <Kokkos_Vector.hpp>

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
  if (!uncopy && (copy || copymode)) return;

  particles = NULL;
  species = NULL;

  deallocate_views_of_views(k_eivec.h_view);
  deallocate_views_of_views(k_eiarray.h_view);
  deallocate_views_of_views(k_edvec.h_view);
  deallocate_views_of_views(k_edarray.h_view);

  eivec = NULL;
  eiarray = NULL;
  edvec = NULL;
  edarray = NULL;

  ewhich = NULL;
  eicol = NULL;
  edcol = NULL;

  ncustom_ivec = ncustom_iarray = 0;
  ncustom_dvec = ncustom_darray = 0;
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
    d_cellcount = {};
    MemKK::realloc_kokkos(grid_kk->d_cellcount,"particle:cellcount",ngrid);
    d_cellcount = grid_kk->d_cellcount;
  }

  Kokkos::deep_copy(d_cellcount,0);

  if (ngrid > int(d_plist.extent(0)) || maxcellcount > int(d_plist.extent(1))) {
    d_plist = {};
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
      maxcellcount = MAX(maxcellcount+MAX(DELTACELLCOUNT,maxcellcount*0.1),resize);
      d_plist = {};
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
  d_plist = {};
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
  } /*else if (METHOD == 2) { // Kokkos-parallel, gives same (correct) answer
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
  }*/
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

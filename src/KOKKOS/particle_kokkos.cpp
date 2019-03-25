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
#include "random_park.h"
#include "memory_kokkos.h"
#include "error.h"
#include "kokkos.h"
#include "sparta_masks.h"

#include <Kokkos_Vector.hpp>

using namespace SPARTA_NS;

enum{PKEEP,PINSERT,PDONE,PDISCARD,PENTRY,PEXIT,PSURF};  // several files
enum{NONE,DISCRETE,SMOOTH};            // several files
enum{INT,DOUBLE};                      // several files
enum{COPYPARTICLELIST,FIXEDMEMORY};

#define DELTA 16384
#define DELTASPECIES 16
#define DELTAMIXTURE 8
#define DELTACELLCOUNT 10
#define MAXLINE 1024

// customize by adding an abbreviation string
// also add a check for the keyword in 2 places in add_species()

#define AIR "N O NO"

/* ---------------------------------------------------------------------- */

ParticleKokkos::ParticleKokkos(SPARTA *sparta) : Particle(sparta)
{
  k_fail_flag = DAT::tdual_int_scalar("particle:fail_flag");
  d_fail_flag = k_fail_flag.view<DeviceType>();
  h_fail_flag = k_fail_flag.h_view;

  k_reorder_pass = DAT::tdual_int_scalar("particle:reorder_pass");
  d_reorder_pass = k_reorder_pass.view<DeviceType>();
  h_reorder_pass = k_reorder_pass.h_view;

  sorted_kk = 0;
  maxcellcount = 10;
}

/* ---------------------------------------------------------------------- */

ParticleKokkos::~ParticleKokkos()
{
  if (copy || copymode) return;

  particles = NULL;
  species = NULL;
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
    d_lists = DAT::t_int_2d(Kokkos::view_alloc("particle:lists",Kokkos::WithoutInitializing),2,ndelete);
    d_mlist = Kokkos::subview(d_lists,0,Kokkos::ALL);
    d_slist = Kokkos::subview(d_lists,1,Kokkos::ALL);

    h_lists = HAT::t_int_2d(Kokkos::view_alloc("particle:lists_mirror",Kokkos::WithoutInitializing),2,ndelete);
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

  this->sync(Device,PARTICLE_MASK);
  d_particles = k_particles.d_view;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagParticleCompressReactions>(0,ncopy),*this);
  DeviceType::fence();
  copymode = 0;

  this->modify(Device,PARTICLE_MASK);
  d_particles = t_particle_1d(); // destroy reference to reduce memory use

  sorted = 0;
  sorted_kk = 0;
}
#endif

KOKKOS_INLINE_FUNCTION
void ParticleKokkos::operator()(TagParticleCompressReactions, const int &i) const {
  const int j = d_mlist[i];
  const int k = d_slist[i];
  //memcpy(&d_particles[j],&d_particles[k],nbytes);
  d_particles[j] = d_particles[k];
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
  if (update->mem_limit_grid_flag)
    update->global_mem_limit = grid->nlocal*sizeof(Grid::ChildCell);
  if (update->global_mem_limit > 0)
    reorder_scheme = FIXEDMEMORY;

  ngrid = grid->nlocal;
  GridKokkos* grid_kk = (GridKokkos*)grid;

  if (ngrid > int(d_cellcount.extent(0)))
    d_cellcount = typename AT::t_int_1d("particle:cellcount",ngrid);
  else {
    copymode = 1;
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagParticleZero_cellcount>(0,d_cellcount.extent(0)),*this);
    DeviceType::fence();
    copymode = 0;
  }

  if (ngrid > int(d_plist.extent(0)) || maxcellcount > int(d_plist.extent(1))) {
    grid_kk->d_plist = typename AT::t_int_2d(); // destroy reference to reduce memory use
    d_plist = typename AT::t_int_2d();
    d_plist = typename AT::t_int_2d(Kokkos::view_alloc("particle:plist",Kokkos::WithoutInitializing),ngrid,maxcellcount);
  }

  this->sync(Device,PARTICLE_MASK);
  d_particles = k_particles.d_view;

  // icell = global cell the particle is in

  // Cannot grow a Kokkos view in a parallel loop, so
  //  if the capacity of the list is exceeded, count the size
  //  needed, reallocate on the host, and then
  //  repeat the parallel loop again

  do {
    copymode = 1;
    if (sparta->kokkos->need_atomics)
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagParticleSort<1> >(0,nlocal),*this);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagParticleSort<0> >(0,nlocal),*this);
    DeviceType::fence();
    copymode = 0;

    k_fail_flag.modify<DeviceType>();
    k_fail_flag.sync<SPAHostType>();
    if (h_fail_flag()) {
      copymode = 1;
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagParticleZero_cellcount>(0,ngrid),*this);
      DeviceType::fence();
      copymode = 0;
      maxcellcount += DELTACELLCOUNT;
      grid_kk->d_plist = typename AT::t_int_2d(); // destroy reference to reduce memory use
      d_plist = typename AT::t_int_2d();
      d_plist = typename AT::t_int_2d(Kokkos::view_alloc("particle:plist",Kokkos::WithoutInitializing),ngrid,maxcellcount);

      h_fail_flag() = 0;
      k_fail_flag.modify<SPAHostType>();
      k_fail_flag.sync<DeviceType>();
    }
  } while (h_fail_flag());

  if (update->reorder_period &&
      (update->ntimestep % update->reorder_period == 0)) {

    if (reorder_scheme == COPYPARTICLELIST && d_particles.extent(0) > d_sorted.extent(0)) {
      d_sorted = t_particle_1d();
      d_sorted = t_particle_1d("particle:sorted",d_particles.extent(0));
    }
    else if (reorder_scheme == FIXEDMEMORY && d_pswap1.size() == 0){
      nParticlesWksp = (double)update->global_mem_limit/sizeof(Particle::OnePart);
      d_pswap1 = t_particle_1d(Kokkos::view_alloc("particle:swap1",Kokkos::WithoutInitializing),nParticlesWksp);
      d_pswap2 = t_particle_1d(Kokkos::view_alloc("particle:swap2",Kokkos::WithoutInitializing),nParticlesWksp);
    }

    nbytes = sizeof(OnePart);

    if (reorder_scheme == COPYPARTICLELIST) {
      copymode = 1;
      Kokkos::parallel_scan(Kokkos::RangePolicy<DeviceType, TagParticleReorder_COPYPARTICLELIST>(0,ngrid),*this);
      copymode = 0;
      Kokkos::deep_copy(k_particles.d_view,d_sorted);
      this->modify(Device,PARTICLE_MASK);
    }
    else if (reorder_scheme == FIXEDMEMORY) {
      // Copy particle destinations into the particle list cell locations
      //  (to avoid adding a "destination" integer in OnePart for the fixed memory reorder)
      // After the particle list has been reordered, reset the icell values to correctly reflect
      // the variable naming.
      copymode = 1;
      Kokkos::parallel_scan(Kokkos::RangePolicy<DeviceType, TagCopyParticleReorderDestinations>(0,ngrid),*this);
      DeviceType::fence();
      copymode = 0;

      int npasses = (nlocal-1)/nParticlesWksp + 1;
      for (int ipass=0; ipass < npasses; ++ipass) {

        h_reorder_pass() = ipass;
        k_reorder_pass.modify<SPAHostType>();
        k_reorder_pass.sync<DeviceType>();

        // identify next set of particles to reorder
        copymode = 1;
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixedMemoryReorderInit>(0,nParticlesWksp),*this);
        DeviceType::fence();
        copymode = 0;

        // reorder this set of particles
        copymode = 1;
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixedMemoryReorder>(0,nParticlesWksp),*this);
        DeviceType::fence();
        copymode = 0;
      }

      // reset the icell values in the particle list
      copymode = 1;
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagSetIcellFromPlist>(0,ngrid),*this);
      DeviceType::fence();
      copymode = 0;
      this->modify(Device,PARTICLE_MASK);

      // destroy references to reduce memory use
      d_pswap1 = t_particle_1d();
      d_pswap2 = t_particle_1d();
    }
  }

  grid_kk->d_cellcount = d_cellcount;
  grid_kk->d_plist = d_plist;

  d_particles = t_particle_1d(); // destroy reference to reduce memory use
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

template<int NEED_ATOMICS>
KOKKOS_INLINE_FUNCTION
void ParticleKokkos::operator()(TagParticleSort<NEED_ATOMICS>, const int &i) const
{
  const int icell = d_particles[i].icell;
  int j;
  if (NEED_ATOMICS)
    j = Kokkos::atomic_fetch_add(&d_cellcount[icell],1);
  else {
    j = d_cellcount[icell];
    d_cellcount[icell]++;
  }
  if (j+1 > maxcellcount)
    d_fail_flag() = 1;
  if (d_fail_flag()) return;
  d_plist(icell,j) = i;
}

KOKKOS_INLINE_FUNCTION
void ParticleKokkos::operator()(TagParticleReorder_COPYPARTICLELIST, const int icell, int &m_fill, const bool &final) const
{
  for (int j = 0; j < d_cellcount[icell]; j++) {
    if (final) {
      const int iparticle = d_plist(icell,j);
      //memcpy(&d_sorted[m_fill],&d_particles[iparticle],nbytes);
      d_sorted[m_fill] = d_particles[iparticle];
      d_plist(icell,j) = m_fill;
    }
    m_fill++;
  }
}

KOKKOS_INLINE_FUNCTION
void ParticleKokkos::operator()(TagParticleZero_cellcount, const int &i) const {
  d_cellcount[i] = 0.0;
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
    Particle::post_weight();
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
  if (ncustom)
    error->all(FLERR,"Custom per-particles attributes not yet supported with Kokkos");

  bigint target = (bigint) nlocal + nextra;
  if (target <= maxlocal) return;

  bigint newmax = maxlocal;
  while (newmax < target) newmax += MAX(DELTA, newmax*1.1);

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

  //for (int i = 0; i < ncustom; i++) {
  //  if (ename[i] == NULL) continue;
  //  grow_custom(i,oldmax,maxlocal);
  //}
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
    k_species.modify<SPAHostType>();
    k_species.sync<DeviceType>();
    memory->sfree(species);
    species = k_species.h_view.data();
  }

  // mixtures

  k_species2group = DAT::tdual_int_2d("particle:species2group",nmixture,nspecies);
  for (int i = 0; i < nmixture; i++)
    for (int j = 0; j < nspecies; j++)
      k_species2group.h_view(i,j) = mixture[i]->species2group[j];
  k_species2group.modify<SPAHostType>();
  k_species2group.sync<DeviceType>();

  //if (mixtures != k_mixtures.h_view.data()) {
  //  memoryKK->wrap_kokkos(k_mixtures,mixture,nmixture,"particle:mixture");
  //  k_mixtures.modify<SPAHostType>();
  //  k_mixtures.sync<DeviceType>();
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
    if (mask & PARTICLE_MASK) k_particles.sync<SPADeviceType>();
    if (mask & SPECIES_MASK) k_species.sync<SPADeviceType>();
  } else {
    if (mask & PARTICLE_MASK) k_particles.sync<SPAHostType>();
    if (mask & SPECIES_MASK) k_species.sync<SPAHostType>();
  }
}

/* ---------------------------------------------------------------------- */

void ParticleKokkos::modify(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if (mask & PARTICLE_MASK) k_particles.modify<SPADeviceType>();
    if (mask & SPECIES_MASK) k_species.modify<SPADeviceType>();
    if (sparta->kokkos->auto_sync)
      sync(Host,mask);
  } else {
    if (mask & PARTICLE_MASK) k_particles.modify<SPAHostType>();
    if (mask & SPECIES_MASK) k_species.modify<SPAHostType>();
  }
}

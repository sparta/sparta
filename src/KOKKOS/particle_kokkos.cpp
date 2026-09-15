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
#include <type_traits>

//#include <Kokkos_Vector.hpp>

using namespace SPARTA_NS;

enum{NONE,DISCRETE,SMOOTH};            // several files
enum{INT,DOUBLE};                      // several files
enum{COPYPARTICLELIST,FIXEDMEMORY};

#define DELTA 16384

// per-cell particle list sizing, see cellcount_target()

#define CELLCOUNT_MARGIN 4
#define CELLCOUNT_GROWTH 1.2

/* ----------------------------------------------------------------------
   spare slots to keep above the fullest cell's count
   scaled by sqrt of that count, not of the mean: where density is
   non-uniform the fullest cell holds many times the mean, and it is its own
   count that sets how far it can fluctuate
------------------------------------------------------------------------- */

static int cellcount_headroom(double count)
{
  if (count < 0.0) count = 0.0;
  return MAX(CELLCOUNT_MARGIN,static_cast<int> (3.0*sqrt(count)));
}

/* ----------------------------------------------------------------------
   first-sort seed for the largest per-cell count, before anything has been
   measured.  assumes a near-uniform gas, where the count is ~Poisson(mean)
   and its max over ngrid cells is near mean + sqrt(2*mean*ln(ngrid)).
   a flow with real density structure exceeds this; being wrong costs one
   retry in the first sort, after which the measured maximum governs
------------------------------------------------------------------------- */

static int cellcount_seed(double mean, int ngrid_in)
{
  double spread = sqrt(2.0*mean*log(1.0*MAX(ngrid_in,2)));
  return static_cast<int> (mean+spread) + 1;
}

/* ----------------------------------------------------------------------
   capacity for the per-cell particle list
   need = smallest capacity known to be required, 0 if nothing is known yet

   sizing to exactly what is required guarantees another realloc soon, as the
   largest per-cell count creeps up while later timesteps sample the tail.
   each realloc is ngrid x maxcellcount ints, GBs at production sizes, plus a
   repeat of the binning pass: ~90 ms on a unified-memory APU, and every other
   rank waits for it at the next collective in Update::move().

   LayoutLeft puts spare slots in columns past the used data, never touched,
   so headroom is free apart from memory.  LayoutRight makes maxcellcount the
   per-cell stride, so padding pushes consecutive cells' lists further apart
   and costs cache-line utilization in kernels that sweep cells in order:
   size to what is required and let the growth factor supply the headroom.
------------------------------------------------------------------------- */

static int cellcount_target(int need, int nlocal_in, int ngrid_in,
                           bool cell_contiguous)
{
  if (ngrid_in <= 0) return MAX(need,CELLCOUNT_MARGIN);

  // fall back on the mean only when nothing is known, i.e. the first sort

  int want = need;
  if (want <= 0) {
    double mean = 1.0*nlocal_in/ngrid_in;
    if (mean >= 1.0) want = cellcount_seed(mean,ngrid_in);
  }

  // below ~1 particle per cell the list is narrow, so a realloc is cheap and
  // padding it across ngrid cells is not worth the memory

  if (want < 1) return 1;

  if (cell_contiguous) want += cellcount_headroom(want);
  return want;
}

/* ---------------------------------------------------------------------- */

ParticleKokkos::ParticleKokkos(SPARTA *sparta) : Particle(sparta)
{
  // NOTE: the weight_rand_pool seed cannot be set here.  Every other Kokkos
  //   class seeds its pool in the constructor initializer list with
  //   12345 + comm->me, but those are all styles the input script creates,
  //   long after SPARTA::create() has finished.  ParticleKokkos is built by
  //   create() itself, at sparta.cpp:484, three lines BEFORE comm exists
  //   (:487), and comm is not NULL-initialized -- so reading comm->me there
  //   dereferences an uninitialized pointer.  Seed on first use instead.

#ifndef SPARTA_KOKKOS_EXACT
  weight_rand_pool_seeded = 0;
#endif


  d_resize = DAT::t_int_scalar("particle:resize");
  h_resize = HAT::t_int_scalar("particle:resize_mirror");

  k_reorder_pass = DAT::tdual_int_scalar("particle:reorder_pass");
  d_reorder_pass = k_reorder_pass.view_device();
  h_reorder_pass = k_reorder_pass.view_host();

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

  ewhich = NULL;
  eicol = NULL;
  edcol = NULL;

  for (int i = 0; i < ncustom_ivec; i++)
    memoryKK->destroy_kokkos(k_eivec.view_host()[i].k_view,eivec[i]);
  for (int i = 0; i < ncustom_iarray; i++)
    memoryKK->destroy_kokkos(k_eiarray.view_host()[i].k_view,eiarray[i]);
  for (int i = 0; i < ncustom_dvec; i++)
    memoryKK->destroy_kokkos(k_edvec.view_host()[i].k_view,edvec[i]);
  for (int i = 0; i < ncustom_darray; i++)
    memoryKK->destroy_kokkos(k_edarray.view_host()[i].k_view,edarray[i]);

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
  d_particles = k_particles.view_device();

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

  // reordering is a memory-locality optimization only, so it is safe to skip
  // custom per-particle data is permuted alongside the OnePart records below,
  //   through the same d_sorted_id map, which is what Particle::reorder()
  //   does on the host via copy_custom()

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

  // maxcellcount tracks the per-cell count that has to fit, not the capacity
  // that is allocated for it: CollideVSSKokkos sizes d_plist to
  // maxcellcount*react_extra, so folding the extent back in here would make
  // the next collide multiply its own padding again, growing d_plist by that
  // factor every timestep.  the binning kernel bounds against the extent
  // instead, so an allocation wider than maxcellcount is still used in full

  // pre-size before the binning pass, so the resize path below is only reached
  // when a cell needs more than this seed.  the resize path is itself
  // measurement-driven -- d_resize carries the largest count any cell actually
  // needed -- so a wrong seed costs one retry, not a wrong capacity

  // read the layout off d_plist itself, so changing its type in grid_kokkos.h
  // cannot silently leave the padding above sized for the other layout

  const bool cell_contiguous =
    std::is_same<decltype(d_plist)::array_layout,Kokkos::LayoutLeft>::value;

  maxcellcount =
    MAX(maxcellcount,cellcount_target(0,nlocal,ngrid,cell_contiguous));

  if (ngrid > int(d_plist.extent(0)) || maxcellcount > int(d_plist.extent(1))) {
    d_plist = {};
    MemKK::realloc_kokkos(grid_kk->d_plist,"particle:plist",ngrid,maxcellcount);
    d_plist = grid_kk->d_plist;
  }

  this->sync(Device,PARTICLE_MASK);
  d_particles = k_particles.view_device();

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

      // grow with headroom, not to exactly what this step needed

      maxcellcount =
        MAX(cellcount_target(resize,nlocal,ngrid,cell_contiguous),
            static_cast<int> (maxcellcount*CELLCOUNT_GROWTH));

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
      //auto tmp = k_particles.view_device();
      //k_particles.view_device() = d_sorted;
      //d_particles = k_particles.view_device();
      //d_sorted = tmp;
      Kokkos::deep_copy(d_particles,d_sorted);

      // permute the custom attributes through the same d_sorted_id map, or
      //   the reorder would silently decouple custom values from their
      //   particles.  gather through a temporary: not an in-place permutation

      if (ncustom) {
        this->sync(Device,CUSTOM_MASK);
        auto l_sorted_id = d_sorted_id;
        const int l_nlocal = nlocal;

        for (int m = 0; m < ncustom_ivec; m++) {
          auto d_src = k_eivec.view_host()[m].k_view.view_device();
          DAT::t_int_1d d_tmp(Kokkos::view_alloc("reorder:custom_ivec",
                                                 Kokkos::WithoutInitializing),l_nlocal);
          Kokkos::parallel_for(l_nlocal, KOKKOS_LAMBDA(int i) {
            d_tmp[i] = d_src[l_sorted_id[i]];
          });
          Kokkos::deep_copy(Kokkos::subview(d_src,Kokkos::make_pair(0,l_nlocal)),d_tmp);
        }

        for (int m = 0; m < ncustom_iarray; m++) {
          auto d_src = k_eiarray.view_host()[m].k_view.view_device();
          const int ncol = d_src.extent(1);
          DAT::t_int_2d d_tmp(Kokkos::view_alloc("reorder:custom_iarray",
                                                 Kokkos::WithoutInitializing),l_nlocal,ncol);
          Kokkos::parallel_for(l_nlocal, KOKKOS_LAMBDA(int i) {
            for (int k = 0; k < ncol; k++) d_tmp(i,k) = d_src(l_sorted_id[i],k);
          });
          Kokkos::deep_copy(Kokkos::subview(d_src,Kokkos::make_pair(0,l_nlocal),
                                            Kokkos::ALL()),d_tmp);
        }

        for (int m = 0; m < ncustom_dvec; m++) {
          auto d_src = k_edvec.view_host()[m].k_view.view_device();
          DAT::t_float_1d d_tmp(Kokkos::view_alloc("reorder:custom_dvec",
                                                   Kokkos::WithoutInitializing),l_nlocal);
          Kokkos::parallel_for(l_nlocal, KOKKOS_LAMBDA(int i) {
            d_tmp[i] = d_src[l_sorted_id[i]];
          });
          Kokkos::deep_copy(Kokkos::subview(d_src,Kokkos::make_pair(0,l_nlocal)),d_tmp);
        }

        for (int m = 0; m < ncustom_darray; m++) {
          auto d_src = k_edarray.view_host()[m].k_view.view_device();
          const int ncol = d_src.extent(1);
          DAT::t_float_2d d_tmp(Kokkos::view_alloc("reorder:custom_darray",
                                                   Kokkos::WithoutInitializing),l_nlocal,ncol);
          Kokkos::parallel_for(l_nlocal, KOKKOS_LAMBDA(int i) {
            for (int k = 0; k < ncol; k++) d_tmp(i,k) = d_src(l_sorted_id[i],k);
          });
          Kokkos::deep_copy(Kokkos::subview(d_src,Kokkos::make_pair(0,l_nlocal),
                                            Kokkos::ALL()),d_tmp);
        }

        this->modify(Device,CUSTOM_MASK);
      }

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

  if (j >= int(d_plist.extent(1)))
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
  auto d_cinfo = k_cinfo.view_device();
  auto d_particles = k_particles.view_device();

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

void ParticleKokkos::post_weight()
{
  // METHOD 1 is the host fallback.  it used to be taken whenever any custom
  //   per-particle attribute existed -- which is any run with fix ambipolar
  //   or fix vibmode -- costing a full particle+custom round trip on every
  //   timestep that grid weighting is active.  METHOD 2 now permutes the
  //   custom arrays with the same map it uses for the particles, so the
  //   fallback is only kept as a reference implementation

#ifndef SPARTA_KOKKOS_EXACT
  // the loop in METHOD 2 below is serial on the host because its
  //   delete-by-swap-from-the-end makes the RNG draw order matter, which
  //   SPARTA_KOKKOS_EXACT needs in order to reproduce Particle::post_weight()
  //   bit-for-bit.  the physics it implements is per-particle independent --
  //   survive with probability ratio, or replicate to 1+nclone copies -- so
  //   away from EXACT it is a prefix-sum scatter with no host round trip
  post_weight_device();
  return;
#endif

  int METHOD = 2;

  if (METHOD == 1) { // just call the host one
    this->sync(Host,PARTICLE_MASK|CUSTOM_MASK);

    auto grid_kk = (GridKokkos*) grid;
    grid_kk->sync(Host,CINFO_MASK);

    int prev_auto_sync = sparta->kokkos->auto_sync;
    sparta->kokkos->auto_sync = 1;
    Particle::post_weight();
    sparta->kokkos->auto_sync = prev_auto_sync;

    this->modify(Host,PARTICLE_MASK|CUSTOM_MASK);
  } else if (METHOD == 2) { // Kokkos-parallel, gives same (correct) answer

    auto grid_kk = dynamic_cast<GridKokkos*>(grid);
    auto& k_cinfo = grid_kk->k_cinfo;
    grid_kk->sync(Device,CINFO_MASK);
    this->sync(Device,PARTICLE_MASK);

    auto d_particles = k_particles.view_device();
    auto d_cinfo = k_cinfo.view_device();

    // k_map persists across calls and only ever grows: this runs every
    //   timestep that weighting is active, and reallocating an O(nlocal)
    //   DualView per step is pure overhead on a GPU

    if ((int) k_map.extent(0) < nlocal)
      MemKK::realloc_kokkos(k_map,"post_weight:map",(size_t)(nlocal*1.5));
    auto d_map = k_map.view_device();
    auto h_map = k_map.view_host();

    // count how many particles actually changed weight while filling the map
    // if none did, the host loop below would only walk the map without
    //   touching it and draw no random numbers, and the final gather would be
    //   an identity permutation -- so skip the whole round trip.  this is
    //   exactly equivalent, RNG stream included

    int nchanged = 0;
    Kokkos::parallel_reduce(nlocal, KOKKOS_LAMBDA(const int i, int &lsum) {
      const auto icell = d_particles[i].icell;
      const double ratio = d_particles[i].weight / d_cinfo[icell].weight;
      d_map[i].ratio = ratio;
      d_map[i].id = d_particles[i].id;
      d_map[i].i = i;
      if (ratio != 1.0) lsum++;
    },nchanged);

    if (!nchanged) {
      d_particles = t_particle_1d();
      return;
    }

    k_map.modify_device();
    k_map.sync_host();

    // nlocal_original-1 = index of last original particle

    int nlocal_original = nlocal;
    int i = 0;

    while (i < nlocal_original) {

      auto ratio = h_map[i].ratio;

      // next particle will be an original particle
      // skip it if no weight change

      if (ratio == 1.0) {
        i++;
        continue;
      }

      // ratio < 1.0 is candidate for deletion
      // if deleted and particle that takes its place is cloned (Nloc > Norig)
      //   then skip it via i++, else will examine it on next iteration

      if (ratio < 1.0) {
        if (wrandom->uniform() > ratio) {
          h_map[i] = h_map[nlocal-1];
          if (nlocal > nlocal_original) i++;
          else nlocal_original--;
          nlocal--;
        } else i++;
        continue;
      }

      // ratio > 1.0 is candidate for cloning
      // create Nclone new particles each with unique ID

      int nclone = static_cast<int>(ratio);
      double fraction = ratio - nclone;
      nclone--;
      if (wrandom->uniform() < fraction) nclone++;

      for (int m = 0; m < nclone; m++) {
        if (nlocal == MAXSMALLINT)
          error->one(FLERR,"Per-processor particle count is too big");
        if (k_map.extent(0) <= nlocal) {
          // 1.5x truncates back to the old size for the smallest views,
          //   so always leave room for at least the particle added below
          size_t newmax = k_map.extent(0)*1.5;
          if (newmax <= (size_t) nlocal) newmax = nlocal + 1;
          k_map.resize(newmax);
          // resize reallocates, so the previously bound host view now points
          //   at the freed buffer -- rebind before writing through it
          h_map = k_map.view_host();
        }

        h_map[nlocal] = h_map[i];
        h_map[nlocal].id = MAXSMALLINT*wrandom->uniform();
        nlocal++;
      }
      i++;
    }

    k_map.modify_host();
    k_map.sync_device();

    grow(0);

    // likewise persistent: this is a full maxlocal-sized particle array, and
    //   allocating and freeing it every timestep is the single largest
    //   avoidable cost in this routine on a GPU

    if ((int) d_newparticles.extent(0) < maxlocal)
      MemKK::realloc_kokkos(d_newparticles,"post_weight:newparticles",maxlocal);
    auto d_newparticles_l = d_newparticles;
    d_map = k_map.view_device();

    Kokkos::parallel_for(nlocal, KOKKOS_LAMBDA(int i) {
      d_newparticles_l[i] = d_particles[d_map[i].i];
      d_newparticles_l[i].id = d_map[i].id;
    });

    Kokkos::deep_copy(k_particles.view_device(),
                      Kokkos::subview(d_newparticles,Kokkos::make_pair(0,(int)k_particles.view_device().extent(0))));
    this->modify(Device,PARTICLE_MASK);

    // permute the custom attributes with the same map
    // a cloned particle carries its source's index in d_map, so it inherits
    //   that particle's custom values, matching Particle::post_weight()
    // gather through a temporary: the permutation is not in place

    if (ncustom) {
      this->sync(Device,CUSTOM_MASK);

      for (int m = 0; m < ncustom_ivec; m++) {
        auto d_src = k_eivec.view_host()[m].k_view.view_device();
        DAT::t_int_1d d_tmp(Kokkos::view_alloc("post_weight:custom_ivec",
                                               Kokkos::WithoutInitializing),nlocal);
        Kokkos::parallel_for(nlocal, KOKKOS_LAMBDA(int i) {
          d_tmp[i] = d_src[d_map[i].i];
        });
        Kokkos::deep_copy(Kokkos::subview(d_src,Kokkos::make_pair(0,nlocal)),d_tmp);
      }

      for (int m = 0; m < ncustom_iarray; m++) {
        auto d_src = k_eiarray.view_host()[m].k_view.view_device();
        const int ncol = d_src.extent(1);
        DAT::t_int_2d d_tmp(Kokkos::view_alloc("post_weight:custom_iarray",
                                               Kokkos::WithoutInitializing),nlocal,ncol);
        Kokkos::parallel_for(nlocal, KOKKOS_LAMBDA(int i) {
          for (int k = 0; k < ncol; k++) d_tmp(i,k) = d_src(d_map[i].i,k);
        });
        Kokkos::deep_copy(Kokkos::subview(d_src,Kokkos::make_pair(0,nlocal),
                                          Kokkos::ALL()),d_tmp);
      }

      for (int m = 0; m < ncustom_dvec; m++) {
        auto d_src = k_edvec.view_host()[m].k_view.view_device();
        DAT::t_float_1d d_tmp(Kokkos::view_alloc("post_weight:custom_dvec",
                                                 Kokkos::WithoutInitializing),nlocal);
        Kokkos::parallel_for(nlocal, KOKKOS_LAMBDA(int i) {
          d_tmp[i] = d_src[d_map[i].i];
        });
        Kokkos::deep_copy(Kokkos::subview(d_src,Kokkos::make_pair(0,nlocal)),d_tmp);
      }

      for (int m = 0; m < ncustom_darray; m++) {
        auto d_src = k_edarray.view_host()[m].k_view.view_device();
        const int ncol = d_src.extent(1);
        DAT::t_float_2d d_tmp(Kokkos::view_alloc("post_weight:custom_darray",
                                                 Kokkos::WithoutInitializing),nlocal,ncol);
        Kokkos::parallel_for(nlocal, KOKKOS_LAMBDA(int i) {
          for (int k = 0; k < ncol; k++) d_tmp(i,k) = d_src(d_map[i].i,k);
        });
        Kokkos::deep_copy(Kokkos::subview(d_src,Kokkos::make_pair(0,nlocal),
                                          Kokkos::ALL()),d_tmp);
      }

      this->modify(Device,CUSTOM_MASK);
    }

    d_particles = t_particle_1d();
  }
}

/* ---------------------------------------------------------------------- */

void ParticleKokkos::update_class_variables() {
  d_species = k_species.view_device();
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
    MemKK::realloc_kokkos(k_particles,"particle:particles",maxlocal);
  else {
    this->sync(Device,PARTICLE_MASK); // force resize on device
    Kokkos::resize(Kokkos::view_alloc(Kokkos::WithoutInitializing),
                   k_particles,maxlocal);
    this->modify(Device,PARTICLE_MASK); // needed for auto sync
  }
  d_particles = k_particles.view_device();
  particles = k_particles.view_host().data();

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
      MemKK::realloc_kokkos(k_species,"particle:species",maxspecies);
    else {
      this->sync(Device,SPECIES_MASK); // force resize on device
      Kokkos::resize(Kokkos::view_alloc(Kokkos::WithoutInitializing),
                     k_species,maxspecies);
      this->modify(Device,SPECIES_MASK); // needed for auto sync
    }
    species = k_species.view_host().data();
  }
}

/* ---------------------------------------------------------------------- */

void ParticleKokkos::wrap_kokkos()
{
  // species

  if (species != k_species.view_host().data()) {
    memoryKK->wrap_kokkos(k_species,species,nspecies,"particle:species");
    k_species.modify_host();
    k_species.sync_device();
    memory->sfree(species);
    species = k_species.view_host().data();
  }

  // mixtures

  k_species2group = DAT::tdual_int_2d("particle:species2group",nmixture,nspecies);
  for (int i = 0; i < nmixture; i++)
    for (int j = 0; j < nspecies; j++)
      k_species2group.view_host()(i,j) = mixture[i]->species2group[j];
  k_species2group.modify_host();
  k_species2group.sync_device();

  //if (mixtures != k_mixtures.view_host().data()) {
  //  memoryKK->wrap_kokkos(k_mixtures,mixture,nmixture,"particle:mixture");
  //  k_mixtures.modify_host();
  //  k_mixtures.sync_device();
  //  memory->sfree(mixtures);
  //  mixtures = k_mixtures.view_host().data();
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
            k_eivec.view_host()[i].k_view.sync_device();

        if (ncustom_iarray)
          for (int i = 0; i < ncustom_iarray; i++)
            k_eiarray.view_host()[i].k_view.sync_device();

        if (ncustom_dvec)
          for (int i = 0; i < ncustom_dvec; i++)
            k_edvec.view_host()[i].k_view.sync_device();

        if (ncustom_darray)
          for (int i = 0; i < ncustom_darray; i++)
            k_edarray.view_host()[i].k_view.sync_device();
      }
    }
  } else {
    if (mask & PARTICLE_MASK) k_particles.sync_host();
    if (mask & SPECIES_MASK) k_species.sync_host();
    if (mask & CUSTOM_MASK) {
      if (ncustom_ivec)
        for (int i = 0; i < ncustom_ivec; i++)
          k_eivec.view_host()[i].k_view.sync_host();

      if (ncustom_iarray)
        for (int i = 0; i < ncustom_iarray; i++)
          k_eiarray.view_host()[i].k_view.sync_host();

      if (ncustom_dvec)
        for (int i = 0; i < ncustom_dvec; i++)
          k_edvec.view_host()[i].k_view.sync_host();

      if (ncustom_darray)
        for (int i = 0; i < ncustom_darray; i++)
          k_edarray.view_host()[i].k_view.sync_host();
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
            k_eivec.view_host()[i].k_view.modify_device();

        if (ncustom_iarray)
          for (int i = 0; i < ncustom_iarray; i++)
            k_eiarray.view_host()[i].k_view.modify_device();

        if (ncustom_dvec)
          for (int i = 0; i < ncustom_dvec; i++)
            k_edvec.view_host()[i].k_view.modify_device();

        if (ncustom_darray)
          for (int i = 0; i < ncustom_darray; i++)
            k_edarray.view_host()[i].k_view.modify_device();
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
            k_eivec.view_host()[i].k_view.modify_host();

        if (ncustom_iarray)
          for (int i = 0; i < ncustom_iarray; i++)
            k_eiarray.view_host()[i].k_view.modify_host();

        if (ncustom_dvec)
          for (int i = 0; i < ncustom_dvec; i++)
            k_edvec.view_host()[i].k_view.modify_host();

        if (ncustom_darray)
          for (int i = 0; i < ncustom_darray; i++)
            k_edarray.view_host()[i].k_view.modify_host();
      }
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of Kokkos-managed data
   Particle::memory_usage() is deliberately not called: the host arrays it
     measures are the host mirrors of the DualViews below, so its formula
     would double count them.  next[] is the one plain host allocation it
     covers with no Kokkos counterpart, so it is carried over here
   the device half is added only when it is a distinct allocation; on a
     host-only backend the two views alias
------------------------------------------------------------------------- */

bigint ParticleKokkos::memory_usage()
{
  const bool device_distinct =
    !std::is_same<DeviceType::memory_space,Kokkos::HostSpace>::value;

  bigint bytes = (bigint) maxlocal * sizeof(int);   // next[]

  bytes += MemKK::memory_usage(k_particles.view_host());
  bytes += MemKK::memory_usage(k_species.view_host());
  bytes += MemKK::memory_usage(k_species2group.view_host());
  for (int i = 0; i < ncustom_ivec; i++)
    bytes += MemKK::memory_usage(k_eivec.view_host()[i].k_view.view_host());
  for (int i = 0; i < ncustom_iarray; i++)
    bytes += MemKK::memory_usage(k_eiarray.view_host()[i].k_view.view_host());
  for (int i = 0; i < ncustom_dvec; i++)
    bytes += MemKK::memory_usage(k_edvec.view_host()[i].k_view.view_host());
  for (int i = 0; i < ncustom_darray; i++)
    bytes += MemKK::memory_usage(k_edarray.view_host()[i].k_view.view_host());

  if (device_distinct) {
    bytes += MemKK::memory_usage(k_particles.view_device());
    bytes += MemKK::memory_usage(k_species.view_device());
    bytes += MemKK::memory_usage(k_species2group.view_device());
    for (int i = 0; i < ncustom_ivec; i++)
      bytes += MemKK::memory_usage(k_eivec.view_host()[i].k_view.view_device());
    for (int i = 0; i < ncustom_iarray; i++)
      bytes += MemKK::memory_usage(k_eiarray.view_host()[i].k_view.view_device());
    for (int i = 0; i < ncustom_dvec; i++)
      bytes += MemKK::memory_usage(k_edvec.view_host()[i].k_view.view_device());
    for (int i = 0; i < ncustom_darray; i++)
      bytes += MemKK::memory_usage(k_edarray.view_host()[i].k_view.view_device());
  }

  // device-only scratch for the sort/reorder path, with no host counterpart
  // in either backend

  bytes += MemKK::memory_usage(d_sorted);
  bytes += MemKK::memory_usage(d_sorted_id);
  bytes += MemKK::memory_usage(d_offsets_part);

  return bytes;
}

#ifndef SPARTA_KOKKOS_EXACT

/* ----------------------------------------------------------------------
   fully on-device post_weight()
   same physics as Particle::post_weight(): a particle whose weight ratio is
     below 1 survives with probability ratio, and one above 1 is replicated to
     1+nclone copies, nclone drawn from the fractional part.  the host version
     realizes that by walking a map and swapping deleted entries in from the
     end, which serializes it; here each particle decides independently, an
     exclusive scan turns the per-particle copy counts into output offsets,
     and one scatter writes the new list
   not bit-compatible with the host RNG stream, which is why EXACT builds keep
     the serial path
------------------------------------------------------------------------- */

void ParticleKokkos::post_weight_device()
{
  if (!nlocal) return;

  auto grid_kk = dynamic_cast<GridKokkos*>(grid);
  grid_kk->sync(Device,CINFO_MASK);
  this->sync(Device,PARTICLE_MASK|CUSTOM_MASK);

  auto d_particles_l = k_particles.view_device();
  auto d_cinfo = grid_kk->k_cinfo.view_device();

  // seed on first use: comm does not exist yet when this class is constructed
  //   (see the constructor).  the seed matches every other Kokkos style's

  if (!weight_rand_pool_seeded) {
    weight_rand_pool =
      Kokkos::Random_XorShift64_Pool<DeviceType>(12345 + comm->me);
    weight_rand_pool_seeded = 1;
  }

  auto l_pool = weight_rand_pool;
  const int nold = nlocal;

  // per-particle output count

  // plain Kokkos::View, not DAT::t_int_1d: offset_scan() takes
  //   Kokkos::View<int*,Device> and the DAT alias does not deduce against it
  //   (see the same note on FixEmitFaceKokkos::d_keep)

  Kokkos::View<int*, DeviceType> d_count("post_weight:count",nold);

  Kokkos::parallel_for(nold, KOKKOS_LAMBDA(const int i) {
    const int icell = d_particles_l[i].icell;
    const double ratio = d_particles_l[i].weight / d_cinfo[icell].weight;

    if (ratio == 1.0) { d_count[i] = 1; return; }

    rand_type rand_gen = l_pool.get_state();
    if (ratio < 1.0) {
      d_count[i] = (rand_gen.drand() > ratio) ? 0 : 1;
    } else {
      int nclone = static_cast<int>(ratio);
      const double fraction = ratio - nclone;
      nclone--;
      if (rand_gen.drand() < fraction) nclone++;
      d_count[i] = 1 + nclone;
    }
    l_pool.free_state(rand_gen);
  });

  // exclusive scan -> output offset of each particle's first copy.
  //   offset_scan() accumulates in 64-bit and aborts with a clear message if
  //   the total exceeds 2^31.  a hand-rolled scan carrying an int update
  //   would wrap to a negative nnew first, and the overflow guard that used
  //   to sit here could never fire: MAXSMALLINT is INT_MAX, so "nnew >
  //   MAXSMALLINT" is false for every int

  int nnew = 0;
  auto d_offset = offset_scan(d_count, nnew);

  // grow to the new count, then scatter

  const int nlocal_save = nlocal;
  nlocal = nnew;
  if (nnew > maxlocal) {
    nlocal = nlocal_save;
    grow(nnew - nlocal_save);
    nlocal = nnew;
  }

  if ((int) d_newparticles.extent(0) < maxlocal)
    MemKK::realloc_kokkos(d_newparticles,"post_weight:newparticles",maxlocal);

  auto d_new = d_newparticles;
  d_particles_l = k_particles.view_device();

  Kokkos::parallel_for(nold, KOKKOS_LAMBDA(const int i) {
    const int n = d_count[i];
    if (!n) return;
    const int base = d_offset[i];
    d_new[base] = d_particles_l[i];
    if (n > 1) {
      rand_type rand_gen = l_pool.get_state();
      for (int k = 1; k < n; k++) {
        d_new[base+k] = d_particles_l[i];
        d_new[base+k].id = MAXSMALLINT*rand_gen.drand();
      }
      l_pool.free_state(rand_gen);
    }
  });

  Kokkos::deep_copy(Kokkos::subview(k_particles.view_device(),Kokkos::make_pair(0,nnew)),
                    Kokkos::subview(d_new,Kokkos::make_pair(0,nnew)));
  this->modify(Device,PARTICLE_MASK);

  // permute the custom attributes through the same offsets

  if (ncustom) {
    for (int m = 0; m < ncustom_ivec; m++) {
      auto d_src = k_eivec.view_host()[m].k_view.view_device();
      DAT::t_int_1d d_tmp("post_weight:cust_iv",nnew);
      Kokkos::parallel_for(nold, KOKKOS_LAMBDA(const int i) {
        for (int k = 0; k < d_count[i]; k++) d_tmp[d_offset[i]+k] = d_src[i];
      });
      Kokkos::deep_copy(Kokkos::subview(d_src,Kokkos::make_pair(0,nnew)),d_tmp);
    }
    for (int m = 0; m < ncustom_iarray; m++) {
      auto d_src = k_eiarray.view_host()[m].k_view.view_device();
      const int ncol = d_src.extent(1);
      DAT::t_int_2d d_tmp("post_weight:cust_ia",nnew,ncol);
      Kokkos::parallel_for(nold, KOKKOS_LAMBDA(const int i) {
        for (int k = 0; k < d_count[i]; k++)
          for (int c = 0; c < ncol; c++) d_tmp(d_offset[i]+k,c) = d_src(i,c);
      });
      Kokkos::deep_copy(Kokkos::subview(d_src,Kokkos::make_pair(0,nnew),Kokkos::ALL()),d_tmp);
    }
    for (int m = 0; m < ncustom_dvec; m++) {
      auto d_src = k_edvec.view_host()[m].k_view.view_device();
      DAT::t_float_1d d_tmp("post_weight:cust_dv",nnew);
      Kokkos::parallel_for(nold, KOKKOS_LAMBDA(const int i) {
        for (int k = 0; k < d_count[i]; k++) d_tmp[d_offset[i]+k] = d_src[i];
      });
      Kokkos::deep_copy(Kokkos::subview(d_src,Kokkos::make_pair(0,nnew)),d_tmp);
    }
    for (int m = 0; m < ncustom_darray; m++) {
      auto d_src = k_edarray.view_host()[m].k_view.view_device();
      const int ncol = d_src.extent(1);
      DAT::t_float_2d d_tmp("post_weight:cust_da",nnew,ncol);
      Kokkos::parallel_for(nold, KOKKOS_LAMBDA(const int i) {
        for (int k = 0; k < d_count[i]; k++)
          for (int c = 0; c < ncol; c++) d_tmp(d_offset[i]+k,c) = d_src(i,c);
      });
      Kokkos::deep_copy(Kokkos::subview(d_src,Kokkos::make_pair(0,nnew),Kokkos::ALL()),d_tmp);
    }
    this->modify(Device,CUSTOM_MASK);
  }

  sorted_kk = 0;
}

#endif

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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "collide_vss_kokkos.h"
#include "grid.h"
#include "update.h"
#include "particle_kokkos.h"
#include "mixture.h"
#include "collide.h"
#include "react.h"
#include "comm.h"
#include "random_knuth.h"
#include "random_mars.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "error.h"
#include "kokkos.h"
#include "sparta_masks.h"
#include "modify.h"
#include "fix.h"
#include "fix_ambipolar.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{NONE,DISCRETE,SMOOTH};            // several files
enum{CONSTANT,VARIABLE};

#define DELTAGRID 1000            // must be bigger than split cells per cell
#define DELTADELETE 1024
#define DELTAELECTRON 128
#define DELTACELLCOUNT 2

#define MAXLINE 1024
#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

CollideVSSKokkos::CollideVSSKokkos(SPARTA *sparta, int narg, char **arg) :
  CollideVSS(sparta, narg, arg),
  rand_pool(12345 + comm->me
#ifdef SPARTA_KOKKOS_EXACT
            , sparta
#endif
            ),
  grid_kk_copy(sparta),
  react_kk_copy(sparta)
{
  kokkos_flag = 1;

  // use 1D view for scalars to reduce GPU memory operations

  d_scalars = t_int_11("collide:scalars");
  h_scalars = t_host_int_11("collide:scalars_mirror");

  d_nattempt_one = Kokkos::subview(d_scalars,0);
  d_ncollide_one = Kokkos::subview(d_scalars,1);
  d_nreact_one   = Kokkos::subview(d_scalars,2);
  d_error_flag   = Kokkos::subview(d_scalars,3);
  d_retry        = Kokkos::subview(d_scalars,4);
  d_maxdelete    = Kokkos::subview(d_scalars,5);
  d_maxcellcount = Kokkos::subview(d_scalars,6);
  d_part_grow    = Kokkos::subview(d_scalars,7);
  d_ndelete      = Kokkos::subview(d_scalars,8);
  d_nlocal       = Kokkos::subview(d_scalars,9);
  d_maxelectron  = Kokkos::subview(d_scalars,10);

  h_nattempt_one = Kokkos::subview(h_scalars,0);
  h_ncollide_one = Kokkos::subview(h_scalars,1);
  h_nreact_one   = Kokkos::subview(h_scalars,2);
  h_error_flag   = Kokkos::subview(h_scalars,3);
  h_retry        = Kokkos::subview(h_scalars,4);
  h_maxdelete    = Kokkos::subview(h_scalars,5);
  h_maxcellcount = Kokkos::subview(h_scalars,6);
  h_part_grow    = Kokkos::subview(h_scalars,7);
  h_ndelete      = Kokkos::subview(h_scalars,8);
  h_nlocal       = Kokkos::subview(h_scalars,9);
  h_maxelectron  = Kokkos::subview(h_scalars,10);

  random_backup = NULL;
  react_defined = 0;

  maxdelete = DELTADELETE;
}

/* ---------------------------------------------------------------------- */

CollideVSSKokkos::~CollideVSSKokkos()
{
  if (copymode) return;

  grid_kk_copy.uncopy();
  react_kk_copy.uncopy();

  memoryKK->destroy_kokkos(k_dellist,dellist);

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.destroy();
  if (random_backup)
    delete random_backup;
#endif
}

/* ---------------------------------------------------------------------- */

void CollideVSSKokkos::init()
{
  // error check

  // initially read-in per-species params must match current species list

  if (nparams != particle->nspecies)
    error->all(FLERR,"VSS parameters do not match current species");

  if (ambiflag && nearcp)
    error->all(FLERR,"Ambipolar collision model does not yet support "
               "near-neighbor collisions");

  // require mixture to contain all species

  int imix = particle->find_mixture(mixID);
  if (imix < 0) error->all(FLERR,"Collision mixture does not exist");
  mixture = particle->mixture[imix];

  if (mixture->nspecies != particle->nspecies)
    error->all(FLERR,"Collision mixture does not contain all species");

  // if rotstyle or vibstyle = DISCRETE,
  // check that extra rotation/vibration info is defined
  // for species that require it

  if (vibstyle == DISCRETE) {
    index_vibmode = particle->find_custom((char *) "vibmode");

    Particle::Species *species = particle->species;
    int nspecies = particle->nspecies;

    int flag = 0;
    for (int isp = 0; isp < nspecies; isp++) {
      if (species[isp].vibdof <= 2) continue;
      if (index_vibmode < 0)
        error->all(FLERR,
                   "Fix vibmode must be used with discrete vibrational modes");
      if (species[isp].nvibmode != species[isp].vibdof / 2) flag++;
    }
    if (flag) {
      char str[128];
      sprintf(str,"%d species do not define correct vibrational "
              "modes for discrete model",flag);
      error->all(FLERR,str);
    }
  }

  // reallocate one-cell data structs for one or many groups

  oldgroups = ngroups;
  ngroups = mixture->ngroup;

  if (ngroups != oldgroups) {
    if (oldgroups == 1) {
      memory->destroy(plist);
      npmax = 0;
      plist = NULL;
    }
    if (oldgroups > 1) {
      delete [] ngroup;
      delete [] maxgroup;
      for (int i = 0; i < oldgroups; i++) memory->destroy(glist[i]);
      delete [] glist;
      memory->destroy(gpair);
      ngroup = NULL;
      maxgroup = NULL;
      glist = NULL;
      gpair = NULL;
    }

    if (ngroups == 1) {
      npmax = DELTAPART;
      memory->create(plist,npmax,"collide:plist");
    }
    if (ngroups > 1) {
      ngroup = new int[ngroups];
      maxgroup = new int[ngroups];
      glist = new int*[ngroups];
      for (int i = 0; i < ngroups; i++) {
        maxgroup[i] = DELTAPART;
        memory->create(glist[i],DELTAPART,"collide:glist");
      }
      memory->create(gpair,ngroups*ngroups,3,"collide:gpair");
    }
  }

  // allocate vremax,remain if group count changed
  // will always be allocated on first run since oldgroups = 0
  // set vremax_intitial via values calculated by collide style

  if (ngroups != oldgroups) {
    memory->destroy(vremax_initial);
    nglocal = grid->nlocal;
    nglocalmax = nglocal;
    memory->create(vremax_initial,ngroups,ngroups,"collide:vremax_initial");

    k_vremax_initial = DAT::tdual_float_2d("collide:vremax_initial",ngroups,ngroups);
    k_vremax = DAT::tdual_float_3d("collide:vremax",nglocalmax,ngroups,ngroups);
    d_vremax = k_vremax.d_view;
    k_remain = DAT::tdual_float_3d("collide:remain",nglocalmax,ngroups,ngroups);
    d_remain = k_remain.d_view;

    for (int igroup = 0; igroup < ngroups; igroup++) {
      for (int jgroup = 0; jgroup < ngroups; jgroup++) {
        vremax_initial[igroup][jgroup] = vremax_init(igroup,jgroup);
        k_vremax_initial.h_view(igroup,jgroup) = vremax_initial[igroup][jgroup];
      }
    }

    k_vremax_initial.modify_host();
    k_vremax_initial.sync_device();
    d_vremax_initial = k_vremax_initial.d_view;
  }

  // if recombination reactions exist, set flags per species pair

  recombflag = 0;
  if (react) {
    react_defined = 1;
    recombflag = react->recombflag;
    recomb_boost_inverse = react->recomb_boost_inverse;
  }

  if (recombflag) {
    int nspecies = particle->nspecies;
    //memory->destroy(recomb_ijflag);
    //memory->create(recomb_ijflag,nspecies,nspecies,"collide:recomb_ijflag");
    d_recomb_ijflag = DAT::t_float_2d("collide:recomb_ijflag",nspecies,nspecies);
    auto h_recomb_ijflag = Kokkos::create_mirror_view(d_recomb_ijflag);
    for (int i = 0; i < nspecies; i++)
      for (int j = 0; j < nspecies; j++)
        h_recomb_ijflag(i,j) = react->recomb_exist(i,j);
    Kokkos::deep_copy(d_recomb_ijflag,h_recomb_ijflag);
  }

  // find ambipolar fix
  // set ambipolar vector/array indices
  // if reactions defined, check that they are valid ambipolar reactions

  if (ambiflag) {
    index_ionambi = particle->find_custom((char *) "ionambi");
    index_velambi = particle->find_custom((char *) "velambi");
    if (index_ionambi < 0 || index_velambi < 0)
      error->all(FLERR,"Collision ambipolar without fix ambipolar");
    if (react) react->ambi_check();

    int ifix;
    for (ifix = 0; ifix < modify->nfix; ifix++)
      if (strcmp(modify->fix[ifix]->style,"ambipolar") == 0) break;
    FixAmbipolar *afix = (FixAmbipolar *) modify->fix[ifix];
    ambispecies = afix->especies;
  }

  // if ambipolar and multiple groups in mixture, ambispecies must be its own group

  if (ambiflag && mixture->ngroup > 1) {
    int *species2group = mixture->species2group;
    int egroup = species2group[ambispecies];
    if (mixture->groupsize[egroup] != 1)
      error->all(FLERR,"Multigroup ambipolar collisions require "
                 "electrons be their own group");
  }

  // vre_next = next timestep to zero vremax & remain, based on vre_every

  if (vre_every) vre_next = (update->ntimestep/vre_every)*vre_every + vre_every;
  else vre_next = update->laststep + 1;

  // if requested reset vremax & remain
  // must be after per-species vremax_initial is setup

  if (vre_first || vre_start) {
    reset_vremax();
    vre_first = 0;
  }

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.init(random);
#endif

  // VSS specific

  k_params = tdual_params_2d("collide_vss:params",nparams,nparams);
  k_prefactor = DAT::tdual_float_2d("collide_vss:prefactor",nparams,nparams);

  for (int i = 0; i < nparams; i++) {
    for (int j = 0; j < nparams; j++){
      k_params.h_view(i,j) = params[i][j];
      k_prefactor.h_view(i,j) = prefactor[i][j];
    }
  }

  k_params.modify_host();
  k_params.sync_device();
  d_params = k_params.d_view;

  k_prefactor.modify_host();
  k_prefactor.sync_device();
  d_prefactor = k_prefactor.d_view;

  // initialize running stats before each run

  ncollide_running = nattempt_running = nreact_running = 0;
}

/* ----------------------------------------------------------------------
   reset vremax to initial species-based values
   reset remain to 0.0
------------------------------------------------------------------------- */

void CollideVSSKokkos::reset_vremax()
{
  grid_kk_copy.copy((GridKokkos*)grid);

  this->sync(Device,ALL_MASK);

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagCollideResetVremax>(0,nglocal),*this);
  copymode = 0;

  this->modified(Device,ALL_MASK);
}

KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::operator()(TagCollideResetVremax, const int &icell) const {
  for (int igroup = 0; igroup < ngroups; igroup++)
    for (int jgroup = 0; jgroup < ngroups; jgroup++) {
      d_vremax(icell,igroup,jgroup) = d_vremax_initial(igroup,jgroup);
      if (remainflag) d_remain(icell,igroup,jgroup) = 0.0;
    }
}

/* ----------------------------------------------------------------------
  NTC algorithm
------------------------------------------------------------------------- */

void CollideVSSKokkos::collisions()
{
  // if requested, reset vrwmax & remain

  if (update->ntimestep == vre_next) {
    reset_vremax();
    vre_next += vre_every;
  }

  // counters

  ncollide_one = nattempt_one = nreact_one = 0;
  h_ndelete() = 0;

  if (sparta->kokkos->atomic_reduction) {
    h_nattempt_one() = 0;
    h_ncollide_one() = 0;
    h_nreact_one() = 0;
  }

  dt = update->dt;
  fnum = update->fnum;
  boltz = update->boltz;

  // perform collisions:
  // variant for single group or multiple groups (not yet supported)
  // variant for nearcp flag or not
  // variant for ambipolar approximation or not

  if (ngroups != 1)
    error->all(FLERR,"Group collisions not yet supported with Kokkos");

  COLLIDE_REDUCE reduce;

  if (!ambiflag) {
    if (nearcp == 0)
      collisions_one<0>(reduce);
    else
      collisions_one<1>(reduce);
  } else {
    collisions_one_ambipolar(reduce);
  }

  // remove any particles deleted in chemistry reactions
  // if particles deleted/created by chemistry, particles are no longer sorted

  if (ndelete) {
    k_dellist.modify_device();
    k_dellist.sync_host();
    ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
#ifndef SPARTA_KOKKOS_EXACT
    particle_kk->compress_migrate(ndelete,dellist);
#else
    particle->compress_reactions(ndelete,dellist);
#endif
  }
  if (react) {
    particle->sorted = 0;
    ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
    particle_kk->sorted_kk = 0;
  }

  // accumulate running totals

  if (sparta->kokkos->atomic_reduction) {
    nattempt_one = h_nattempt_one();
    ncollide_one = h_ncollide_one();
    nreact_one = h_nreact_one();
  } else {
    nattempt_one += reduce.nattempt_one;
    ncollide_one += reduce.ncollide_one;
    nreact_one += reduce.nreact_one;
  }

  nattempt_running += nattempt_one;
  ncollide_running += ncollide_one;
  nreact_running += nreact_one;
}

/* ----------------------------------------------------------------------
   NTC algorithm for a single group
------------------------------------------------------------------------- */

template < int NEARCP > void CollideVSSKokkos::collisions_one(COLLIDE_REDUCE &reduce)
{
  // loop over cells I own

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,PARTICLE_MASK|SPECIES_MASK);
  if (vibstyle == DISCRETE) particle_kk->sync(Device,CUSTOM_MASK);
  d_particles = particle_kk->k_particles.d_view;
  d_species = particle_kk->k_species.d_view;
  d_ewhich = particle_kk->k_ewhich.d_view;
  k_eiarray = particle_kk->k_eiarray;

  GridKokkos* grid_kk = (GridKokkos*) grid;
  grid_kk->sync(Device,CINFO_MASK);
  d_plist = grid_kk->d_plist;

  grid_kk_copy.copy(grid_kk);

  if (react) {
    ReactTCEKokkos* react_kk = (ReactTCEKokkos*) react;
    if (!react_kk)
      error->all(FLERR,"Must use TCE reactions with Kokkos");
    react_kk_copy.copy(react_kk);
  }

  copymode = 1;

  if (NEARCP) {
    if (int(d_nn_last_partner.extent(0)) < nglocal || int(d_nn_last_partner.extent(1)) < d_plist.extent(1))
      MemKK::realloc_kokkos(d_nn_last_partner,"collide:nn_last_partner",nglocal,d_plist.extent(1));
    //Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagCollideZeroNN>(0,nglocal),*this);
  }

  /* ATOMIC_REDUCTION: 1 = use atomics
                       0 = don't need atomics
                      -1 = use parallel_reduce
  */

  // Reactions may create or delete more particles than existing views can hold.
  //  Cannot grow a Kokkos view in a parallel loop, so
  //  if the capacity of the view is exceeded, break out of parallel loop,
  //  reallocate on the host, and then repeat the parallel loop again.
  //  Unfortunately this leads to really messy code.

  h_retry() = 1;

  if (react) {
    double extra_factor = 1.0;
    if (sparta->kokkos->react_retry_flag)
      extra_factor = sparta->kokkos->react_extra;

    auto maxdelete_extra = maxdelete*extra_factor;
    if (d_dellist.extent(0) < maxdelete_extra) {
      memoryKK->destroy_kokkos(k_dellist,dellist);
      memoryKK->create_kokkos(k_dellist,dellist,maxdelete_extra,"collide:dellist");
      d_dellist = k_dellist.d_view;
    }

    maxcellcount = particle_kk->get_maxcellcount();
    auto maxcellcount_extra = maxcellcount*extra_factor;
    if (d_plist.extent(1) < maxcellcount_extra) {
      d_plist = decltype(d_plist)();
      Kokkos::resize(grid_kk->d_plist,nglocal,maxcellcount_extra);
      d_plist = grid_kk->d_plist;
      if (NEARCP)
        MemKK::realloc_kokkos(d_nn_last_partner,"collide:nn_last_partner",nglocal,maxcellcount_extra);
    }

    auto nlocal_extra = particle->nlocal*extra_factor;
    if (d_particles.extent(0) < nlocal_extra) {
      particle->grow(nlocal_extra - particle->nlocal);
      d_particles = particle_kk->k_particles.d_view;
      k_eiarray = particle_kk->k_eiarray;
    }
  }

  while (h_retry()) {

    if (react && sparta->kokkos->react_retry_flag)
      backup();

    h_retry() = 0;
    h_maxdelete() = maxdelete;
    h_maxcellcount() = maxcellcount;
    h_part_grow() = 0;
    h_ndelete() = 0;
    h_nlocal() = particle->nlocal;

    Kokkos::deep_copy(d_scalars,h_scalars);

    if (sparta->kokkos->atomic_reduction) {
      if (sparta->kokkos->need_atomics)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagCollideCollisionsOne<NEARCP,1> >(0,nglocal),*this);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagCollideCollisionsOne<NEARCP,0> >(0,nglocal),*this);
    } else
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagCollideCollisionsOne<NEARCP,-1> >(0,nglocal),*this,reduce);

    Kokkos::deep_copy(h_scalars,d_scalars);

    if (h_retry()) {
      //printf("Retrying, reason %i %i %i !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n",h_maxdelete() > d_dellist.extent(0),h_maxcellcount() > d_plist.extent(1),h_part_grow());
      if (!sparta->kokkos->react_retry_flag) {
        error->one(FLERR,"Ran out of space in Kokkos collisions, increase collide/extra"
                         " or use collide/retry");
      } else
        restore();

      reduce = COLLIDE_REDUCE();

      maxdelete = h_maxdelete();
      if (d_dellist.extent(0) < maxdelete) {
        memoryKK->destroy_kokkos(k_dellist,dellist);
        memoryKK->grow_kokkos(k_dellist,dellist,maxdelete,"collide:dellist");
        d_dellist = k_dellist.d_view;
      }

      maxcellcount = h_maxcellcount();
      particle_kk->set_maxcellcount(maxcellcount);
      if (d_plist.extent(1) < maxcellcount) {
        d_plist = decltype(d_plist)();
        Kokkos::resize(grid_kk->d_plist,nglocal,maxcellcount);
        d_plist = grid_kk->d_plist;
      }

      auto nlocal_new = h_nlocal();
      if (d_particles.extent(0) < nlocal_new) {
        particle->grow(nlocal_new - particle->nlocal);
        d_particles = particle_kk->k_particles.d_view;
        k_eiarray = particle_kk->k_eiarray;
      }
    }
  }

  ndelete = h_ndelete();

  particle->nlocal = h_nlocal();

  copymode = 0;

  if (h_error_flag())
    error->one(FLERR,"Collision cell volume is zero");

  particle_kk->modify(Device,PARTICLE_MASK);

  d_particles = t_particle_1d(); // destroy reference to reduce memory use
  d_nn_last_partner = decltype(d_nn_last_partner)();
  d_plist = decltype(d_nn_last_partner)();
}

KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::operator()(TagCollideZeroNN, const int &icell) const {
  const int np = grid_kk_copy.obj.d_cellcount[icell];
  for (int i = 0; i < np; i++)
    d_nn_last_partner(icell,i) = 0;
}

template < int NEARCP, int ATOMIC_REDUCTION >
KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::operator()(TagCollideCollisionsOne< NEARCP, ATOMIC_REDUCTION >, const int &icell) const {
  COLLIDE_REDUCE reduce;
  this->template operator()< NEARCP, ATOMIC_REDUCTION >(TagCollideCollisionsOne< NEARCP, ATOMIC_REDUCTION >(), icell, reduce);
}

template < int NEARCP, int ATOMIC_REDUCTION >
KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::operator()(TagCollideCollisionsOne< NEARCP, ATOMIC_REDUCTION >, const int &icell, COLLIDE_REDUCE &reduce) const {
  if (d_retry()) return;

  int np = grid_kk_copy.obj.d_cellcount[icell];
  if (np <= 1) return;

  if (NEARCP) {
    for (int i = 0; i < np; i++)
      d_nn_last_partner(icell,i) = 0;
  }

  const double volume = grid_kk_copy.obj.k_cinfo.d_view[icell].volume / grid_kk_copy.obj.k_cinfo.d_view[icell].weight;
  if (volume == 0.0) d_error_flag() = 1;

  struct State precoln;       // state before collision
  struct State postcoln;      // state after collision

  rand_type rand_gen = rand_pool.get_state();

  // attempt = exact collision attempt count for a pair of groups
  // nattempt = rounded attempt with RN

  const double attempt = attempt_collision_kokkos(icell,np,volume,rand_gen);
  const int nattempt = static_cast<int> (attempt);
  if (!nattempt){
    rand_pool.free_state(rand_gen);
    return;
  }
  if (ATOMIC_REDUCTION == 1)
    Kokkos::atomic_add(&d_nattempt_one(),nattempt);
  else if (ATOMIC_REDUCTION == 0)
    d_nattempt_one() += nattempt;
  else
    reduce.nattempt_one += nattempt;

  // perform collisions
  // select random pair of particles, cannot be same
  // test if collision actually occurs

  for (int m = 0; m < nattempt; m++) {
    const int i = np * rand_gen.drand();
    int j;
    if (NEARCP) j = find_nn(rand_gen,i,np,icell);
    else {
      j = np * rand_gen.drand();
      while (i == j) j = np * rand_gen.drand();
    }

    Particle::OnePart* ipart = &d_particles[d_plist(icell,i)];
    Particle::OnePart* jpart = &d_particles[d_plist(icell,j)];
    Particle::OnePart* kpart;

    // test if collision actually occurs, then perform it
    // ijspecies = species before collision chemistry
    // continue to next collision if no reaction

    if (!test_collision_kokkos(icell,0,0,ipart,jpart,precoln,rand_gen)) continue;

    if (NEARCP) {
      d_nn_last_partner(icell,i) = j+1;
      d_nn_last_partner(icell,j) = i+1;
    }

    // if recombination reaction is possible for this IJ pair
    // pick a 3rd particle to participate and set cell number density
    // unless boost factor turns it off, or there is no 3rd particle

    Particle::OnePart* recomb_part3 = NULL;
    int recomb_species = -1;
    double recomb_density = 0.0;
    if (recombflag && d_recomb_ijflag(ipart->ispecies,jpart->ispecies)) {
      if (rand_gen.drand() > recomb_boost_inverse)
        //react->recomb_species = -1;
        recomb_species = -1;
      else if (np <= 2)
        //react->recomb_species = -1;
        recomb_species = -1;
      else {
        int k = np * rand_gen.drand();
        while (k == i || k == j) k = np * rand_gen.drand();
        // NOT thread safe
        //react->recomb_part3 = &particles[plist[k]];
        //react->recomb_species = react->recomb_part3->ispecies;
        //react->recomb_density = np * update->fnum / volume;
        recomb_part3 = &d_particles[d_plist(icell,k)];
        recomb_species = recomb_part3->ispecies;
        recomb_density = np * fnum / volume;
      }
    }

    // perform collision and possible reaction

    int index_kpart;

    setup_collision_kokkos(ipart,jpart,precoln,postcoln);
    const int reactflag = perform_collision_kokkos(ipart,jpart,kpart,precoln,postcoln,rand_gen,
                                                   recomb_part3,recomb_species,recomb_density,index_kpart);

    if (ATOMIC_REDUCTION == 1)
      Kokkos::atomic_increment(&d_ncollide_one());
    else if (ATOMIC_REDUCTION == 0)
      d_ncollide_one()++;
    else
      reduce.ncollide_one++;

    if (reactflag) {
      if (ATOMIC_REDUCTION == 1)
        Kokkos::atomic_increment(&d_nreact_one());
      else if (ATOMIC_REDUCTION == 0)
        d_nreact_one()++;
      else
        reduce.nreact_one++;
    } else {
      rand_pool.free_state(rand_gen);
      continue;
    }

    // if jpart destroyed, delete from plist
    // also add particle to deletion list
    // exit attempt loop if only single particle left

    if (!jpart) {
      int ndelete = Kokkos::atomic_fetch_add(&d_ndelete(),1);
      if (ndelete < d_dellist.extent(0)) {
        d_dellist(ndelete) = d_plist(icell,j);
      } else {
        d_retry() = 1;
        d_maxdelete() += DELTADELETE;
        rand_pool.free_state(rand_gen);
        return;
      }
      np--;
      d_plist(icell,j) = d_plist(icell,np);
      if (NEARCP) d_nn_last_partner(icell,j) = d_nn_last_partner(icell,np);
      if (np < 2) break;
    }

    // if kpart created, add to plist
    // kpart was just added to particle list, so index = nlocal-1
    // particle data structs may have been realloced by kpart

    if (kpart) {
      if (np < d_plist.extent(1)) {
        if (NEARCP) d_nn_last_partner(icell,np) = 0;
        d_plist(icell,np++) = index_kpart;
      } else {
        d_retry() = 1;
        d_maxcellcount() += DELTACELLCOUNT;
        rand_pool.free_state(rand_gen);
        return;
      }

    }
  }
  rand_pool.free_state(rand_gen);
}

/* ----------------------------------------------------------------------
   NTC algorithm for a single group with ambipolar approximation
------------------------------------------------------------------------- */

void CollideVSSKokkos::collisions_one_ambipolar(COLLIDE_REDUCE &reduce)
{
  // ambipolar vectors

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,PARTICLE_MASK|SPECIES_MASK|CUSTOM_MASK);
  d_particles = particle_kk->k_particles.d_view;
  d_species = particle_kk->k_species.d_view;
  d_ewhich = particle_kk->k_ewhich.d_view;
  auto h_ewhich = particle_kk->k_ewhich.h_view;
  k_eivec = particle_kk->k_eivec;
  k_eiarray = particle_kk->k_eiarray;
  k_edarray = particle_kk->k_edarray;
  d_ionambi = k_eivec.h_view[h_ewhich[index_ionambi]].k_view.d_view;
  d_velambi = k_edarray.h_view[h_ewhich[index_velambi]].k_view.d_view;

  GridKokkos* grid_kk = (GridKokkos*) grid;
  grid_kk->sync(Device,CINFO_MASK);
  d_plist = grid_kk->d_plist;

  grid_kk_copy.copy(grid_kk);

  if (react) {
    ReactTCEKokkos* react_kk = (ReactTCEKokkos*) react;
    if (!react_kk)
      error->all(FLERR,"Must use TCE reactions with Kokkos");
    react_kk_copy.copy(react_kk);
  }

  copymode = 1;

  /* ATOMIC_REDUCTION: 1 = use atomics
                       0 = don't need atomics
                      -1 = use parallel_reduce
  */

  // Reactions may create or delete more particles than existing views can hold.
  //  Cannot grow a Kokkos view in a parallel loop, so
  //  if the capacity of the view is exceeded, break out of parallel loop,
  //  reallocate on the host, and then repeat the parallel loop again.
  //  Unfortunately this leads to really messy code.

  h_retry() = 1;

  if (react) {
    double extra_factor = 1.0;
    if (sparta->kokkos->react_retry_flag)
      extra_factor = sparta->kokkos->react_extra;
 
    maxcellcount = particle_kk->get_maxcellcount();

    auto maxelectron_extra = maxcellcount*extra_factor;
    if (d_elist.extent(0) < nglocal || d_elist.extent(1) < maxelectron_extra) {
      d_elist = t_particle_2d(); // reduce memory use by deallocating first
      d_elist = t_particle_2d(Kokkos::view_alloc("collide:elist",Kokkos::WithoutInitializing),nglocal,maxelectron_extra);
    }

    auto maxdelete_extra = maxdelete*extra_factor;
    if (d_dellist.extent(0) < maxdelete_extra) {
      memoryKK->destroy_kokkos(k_dellist,dellist);
      memoryKK->grow_kokkos(k_dellist,dellist,maxdelete_extra,"collide:dellist");
      d_dellist = k_dellist.d_view;
    }

    auto maxcellcount_extra = maxcellcount*extra_factor;
    if (d_plist.extent(1) < maxcellcount_extra) {
      d_plist = decltype(d_plist)();
      Kokkos::resize(grid_kk->d_plist,nglocal,maxcellcount_extra);
      d_plist = grid_kk->d_plist;
    }

    auto nlocal_extra = particle->nlocal*extra_factor;
    if (d_particles.extent(0) < nlocal_extra) {
      particle->grow(nlocal_extra - particle->nlocal);
      d_particles = particle_kk->k_particles.d_view;
      auto h_ewhich = particle_kk->k_ewhich.h_view;
      k_eivec = particle_kk->k_eivec;
      k_eiarray = particle_kk->k_eiarray;
      k_edarray = particle_kk->k_edarray;
      d_ionambi = k_eivec.h_view[h_ewhich[index_ionambi]].k_view.d_view;
      d_velambi = k_edarray.h_view[h_ewhich[index_velambi]].k_view.d_view;
    }
  }

  while (h_retry()) {

    if (react && sparta->kokkos->react_retry_flag)
      backup();

    h_retry() = 0;
    h_maxelectron() = maxelectron;
    h_maxdelete() = maxdelete;
    h_maxcellcount() = maxcellcount;
    h_part_grow() = 0;
    h_ndelete() = 0;
    h_nlocal() = particle->nlocal;

    Kokkos::deep_copy(d_scalars,h_scalars);

    if (sparta->kokkos->atomic_reduction) {
      if (sparta->kokkos->need_atomics)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagCollideCollisionsOneAmbipolar<1> >(0,nglocal),*this);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagCollideCollisionsOneAmbipolar<0> >(0,nglocal),*this);
    } else
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagCollideCollisionsOneAmbipolar<-1> >(0,nglocal),*this,reduce);

    Kokkos::deep_copy(h_scalars,d_scalars);

    if (h_retry()) {
      //printf("Retrying, reason %i %i %i %i !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n",h_maxelectron() > d_elist.extent(1),h_maxdelete() > d_dellist.extent(0),h_maxcellcount() > d_plist.extent(1),h_part_grow());
      //printf("%i %i %i %i %i %i %i !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n",h_maxelectron(),d_elist.extent(1),h_maxdelete(),d_dellist.extent(0),h_maxcellcount(),d_plist.extent(1),h_part_grow());

      if (!sparta->kokkos->react_retry_flag) {
        error->one(FLERR,"Ran out of space in Kokkos collisions, increase collide/extra"
                         " or use collide/retry");
      } else
        restore();

      reduce = COLLIDE_REDUCE();

      maxelectron = h_maxelectron();
      if (d_elist.extent(1) < maxelectron) {
        d_elist = t_particle_2d(); // reduce memory use by deallocating first
        d_elist = t_particle_2d(Kokkos::view_alloc("collide:elist",Kokkos::WithoutInitializing),nglocal,maxelectron);
      }

      maxdelete = h_maxdelete();
      if (d_dellist.extent(0) < maxdelete) {
        memoryKK->destroy_kokkos(k_dellist,dellist);
        memoryKK->grow_kokkos(k_dellist,dellist,maxdelete,"collide:dellist");
        d_dellist = k_dellist.d_view;
      }

      maxcellcount = h_maxcellcount();
      particle_kk->set_maxcellcount(maxcellcount);
      if (d_plist.extent(1) < maxcellcount) {
        d_plist = decltype(d_plist)();
        Kokkos::resize(grid_kk->d_plist,nglocal,maxcellcount);
        d_plist = grid_kk->d_plist;
      }

      auto nlocal_new = h_nlocal();
      if (d_particles.extent(0) < nlocal_new) {
        particle->grow(nlocal_new - particle->nlocal);
        d_particles = particle_kk->k_particles.d_view;
        auto h_ewhich = particle_kk->k_ewhich.h_view;
        k_eivec = particle_kk->k_eivec;
        k_eiarray = particle_kk->k_eiarray;
        k_edarray = particle_kk->k_edarray;
        d_ionambi = k_eivec.h_view[h_ewhich[index_ionambi]].k_view.d_view;
        d_velambi = k_edarray.h_view[h_ewhich[index_velambi]].k_view.d_view;
      }
    }
  }

  ndelete = h_ndelete();

  particle->nlocal = h_nlocal();

  copymode = 0;

  if (h_error_flag() == 1)
    error->one(FLERR,"Collision cell volume is zero");
  else if (h_error_flag() == 2)
    error->one(FLERR,"Collisions in cell did not conserve electron count");

  particle_kk->modify(Device,PARTICLE_MASK|CUSTOM_MASK);

  d_particles = t_particle_1d(); // destroy reference to reduce memory use
  d_plist = decltype(d_nn_last_partner)();
}

template < int ATOMIC_REDUCTION >
KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::operator()(TagCollideCollisionsOneAmbipolar< ATOMIC_REDUCTION >, const int &icell) const {
  COLLIDE_REDUCE reduce;
  this->template operator()< ATOMIC_REDUCTION >(TagCollideCollisionsOneAmbipolar< ATOMIC_REDUCTION >(), icell, reduce);
}

template < int ATOMIC_REDUCTION >
KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::operator()(TagCollideCollisionsOneAmbipolar< ATOMIC_REDUCTION >, const int &icell, COLLIDE_REDUCE &reduce) const {
  if (d_retry()) return;

  int np = grid_kk_copy.obj.d_cellcount[icell];
  if (np <= 1) return;

  const double volume = grid_kk_copy.obj.k_cinfo.d_view[icell].volume / grid_kk_copy.obj.k_cinfo.d_view[icell].weight;
  if (volume == 0.0) d_error_flag() = 1;

  struct State precoln;       // state before collision
  struct State postcoln;      // state after collision

  int i,j;
  Particle::OnePart *ipart,*jpart,*kpart,*p,*ep;

  rand_type rand_gen = rand_pool.get_state();

  // setup elist of ionized electrons for this cell
  // create them in separate array since will never become real particles
  // create electrons for ambipolar ions

  int nelectron = 0;
  for (i = 0; i < np; i++) {
    if (d_ionambi[d_plist(icell,i)]) {
      p = &d_particles[d_plist(icell,i)];
      ep = &d_elist(icell,nelectron);
      // memcpy(ep,p,nbytes);
      *ep = *p;
      //memcpy(ep->v,d_velambi[d_plist(icell,i)],3*sizeof(double));
      ep->v[0] = d_velambi(d_plist(icell,i),0);
      ep->v[1] = d_velambi(d_plist(icell,i),1);
      ep->v[2] = d_velambi(d_plist(icell,i),2);
      ep->ispecies = ambispecies;
      nelectron++;
    }
  }

  // attempt = exact collision attempt count for all particles in cell
  // nptotal = includes neutrals, ions, electrons
  // nattempt = rounded attempt with RN

  int nptotal = np + nelectron;
  const double attempt = attempt_collision_kokkos(icell,nptotal,volume,rand_gen);
  const int nattempt = static_cast<int> (attempt);
  if (!nattempt) {
    rand_pool.free_state(rand_gen);
    return;
  }
  if (ATOMIC_REDUCTION == 1)
    Kokkos::atomic_fetch_add(&d_nattempt_one(),nattempt);
  else if (ATOMIC_REDUCTION == 0)
    d_nattempt_one() += nattempt;
  else
    reduce.nattempt_one += nattempt;

  // perform collisions
  // select random pair of particles, cannot be same
  // test if collision actually occurs
  // if chemistry occurs, exit attempt loop if group count goes to 0

  for (int iattempt = 0; iattempt < nattempt; iattempt++) {
    i = nptotal * rand_gen.drand();
    j = nptotal * rand_gen.drand();
    while (i == j) j = nptotal * rand_gen.drand();

    // ipart,jpart = heavy particles or electrons

    if (i < np) ipart = &d_particles[d_plist(icell,i)];
    else ipart = &d_elist(icell,i-np);
    if (j < np) jpart = &d_particles[d_plist(icell,j)];
    else jpart = &d_elist(icell,j-np);

    // check for e/e pair
    // count as collision, but do not perform it

    if (ipart->ispecies == ambispecies && jpart->ispecies == ambispecies) {
      if (ATOMIC_REDUCTION == 1)
        Kokkos::atomic_fetch_add(&d_ncollide_one(),1);
      else if (ATOMIC_REDUCTION == 0)
        d_ncollide_one()++;
      else
        reduce.ncollide_one++;

      continue;
    }

    // if particle I is electron
    // swap with J, since electron must be 2nd in any ambipolar reaction
    // just need to swap i/j, ipart/jpart
    // don't have to worry if an ambipolar ion is I or J

    if (ipart->ispecies == ambispecies) {
      int tmp = i;
      i = j;
      j = tmp;
      p = ipart;
      ipart = jpart;
      jpart = p;
    }

    // test if collision actually occurs

    if (!test_collision_kokkos(icell,0,0,ipart,jpart,precoln,rand_gen)) continue;

    // if recombination reaction is possible for this IJ pair
    // pick a 3rd particle to participate and set cell number density
    // unless boost factor turns it off, or there is no 3rd particle
    // 3rd particle cannot be an electron, so select from Np

    Particle::OnePart* recomb_part3 = NULL;
    int recomb_species = -1;
    double recomb_density = 0.0;
    if (recombflag && d_recomb_ijflag(ipart->ispecies,jpart->ispecies)) {
      if (rand_gen.drand() > recomb_boost_inverse)
        //react->recomb_species = -1;
        recomb_species = -1;
      else if (np <= 2)
        //react->recomb_species = -1;
        recomb_species = -1;
      else {
        int k = np * rand_gen.drand();
        while (k == i || k == j) k = np * rand_gen.drand();
        // NOT thread safe
        //react->recomb_part3 = &particles[plist[k]];
        //react->recomb_species = react->recomb_part3->ispecies;
        //react->recomb_density = np * update->fnum / volume;
        recomb_part3 = &d_particles[d_plist(icell,k)];
        recomb_species = recomb_part3->ispecies;
        recomb_density = np * fnum / volume;
      }
    }

    // perform collision
    // ijspecies = species before collision chemistry
    // continue to next collision if no reaction

    int index_kpart = 0;

    const int jspecies = jpart->ispecies;
    setup_collision_kokkos(ipart,jpart,precoln,postcoln);
    const int reactflag = perform_collision_kokkos(ipart,jpart,kpart,precoln,postcoln,rand_gen,
                                                   recomb_part3,recomb_species,recomb_density,index_kpart);

    if (ATOMIC_REDUCTION == 1)
      Kokkos::atomic_fetch_add(&d_ncollide_one(),1);
    else if (ATOMIC_REDUCTION == 0)
      d_ncollide_one()++;
    else
      reduce.ncollide_one++;

    if (reactflag) {
      if (ATOMIC_REDUCTION == 1)
        Kokkos::atomic_fetch_add(&d_nreact_one(),1);
      else if (ATOMIC_REDUCTION == 0)
        d_nreact_one()++;
      else
        reduce.nreact_one++;
    } else {
      rand_pool.free_state(rand_gen);
      continue;
    }

    // reset ambipolar ion flags due to collision
    // must do now before particle count reset below can break out of loop
    // first reset ionambi if kpart was added since ambi_reset() uses it

    if (jspecies == ambispecies)
      ambi_reset_kokkos(d_plist(icell,i),-1,jspecies,index_kpart,ipart,jpart,kpart,d_ionambi);
    else
      ambi_reset_kokkos(d_plist(icell,i),d_plist(icell,j),jspecies,index_kpart,ipart,jpart,kpart,d_ionambi);

    // if kpart created:
    // particles and custom data structs may have been realloced by kpart
    // add kpart to plist or elist
    // kpart was just added to particle list, so index = nlocal-1
    // must come before jpart code below since it modifies nlocal

    if (kpart) {
      if (kpart->ispecies != ambispecies) {
        if (np < d_plist.extent(1)) {
          d_plist(icell,np++) = index_kpart;
        } else {
          d_retry() = 1;
          d_maxcellcount() += DELTACELLCOUNT;
          rand_pool.free_state(rand_gen);
          return;
        }

      } else {

        if (nelectron < d_elist.extent(1)) {
          ep = &d_elist(icell,nelectron);
          //memcpy(ep,kpart,nbytes);
          *ep = *kpart;
          ep->ispecies = ambispecies;
          nelectron++;
#ifdef SPARTA_KOKKOS_EXACT
          d_nlocal()--;
#else
          int ndelete = Kokkos::atomic_fetch_add(&d_ndelete(),1);
          if (ndelete < d_dellist.extent(0)) {
            d_dellist(ndelete) = index_kpart;
          } else {
            d_retry() = 1;
            d_maxdelete() += DELTADELETE;
            rand_pool.free_state(rand_gen);
            return;
          }
#endif
        } else {
          d_retry() = 1;
          d_maxelectron() += DELTACELLCOUNT;
          rand_pool.free_state(rand_gen);
          return;
        }
      }
    }

    // if jpart exists, was originally not an electron, now is an electron:
    //   ionization reaction converted 2 neutrals to one ion
    //   add to elist, remove from plist, flag J for deletion
    // if jpart exists, was originally an electron, now is not an electron:
    //   exchange reaction converted ion + electron to two neutrals
    //   add neutral J to master particle list, remove from elist, add to plist
    // if jpart destroyed, was an electron:
    //   recombination reaction converted ion + electron to one neutral
    //   remove electron from elist
    // else if jpart destroyed:
    //   non-ambipolar recombination reaction
    //   remove from plist, flag J for deletion

    if (jpart) {
      if (jspecies != ambispecies && jpart->ispecies == ambispecies) {
        if (nelectron < d_elist.extent(1)) {
          ep = &d_elist(icell,nelectron);
          //memcpy(ep,jpart,nbytes);
          *ep = *jpart;
          ep->ispecies = ambispecies;
          nelectron++;
          jpart = NULL;
        } else {
          d_retry() = 1;
          d_maxelectron() += DELTACELLCOUNT;
          rand_pool.free_state(rand_gen);
          return;
        }

      } else if (jspecies == ambispecies && jpart->ispecies != ambispecies) {
        int index = Kokkos::atomic_fetch_add(&d_nlocal(),1);
        int reallocflag = ParticleKokkos::add_particle_kokkos(d_particles,index,0,jspecies,icell,jpart->x,jpart->v,0.0,0.0);
        if (reallocflag) {
          d_retry() = 1;
          d_part_grow() = 1;
          rand_pool.free_state(rand_gen);
          return;
        }

        //memcpy(&particles[index],jpart,nbytes);
        d_particles[index] = *jpart;
        d_particles[index].id = MAXSMALLINT*rand_gen.drand();
        d_ionambi[index] = 0;

        //if (nelectron-1 != j-np) memcpy(&d_elist(icell,j-np),&d_elist(icell,nelectron-1),nbytes);
        if (nelectron-1 != j-np) d_elist(icell,j-np) = d_elist(icell,nelectron-1);
        nelectron--;

        if (np < d_plist.extent(1)) {
          d_plist(icell,np++) = index;
        } else {
          d_retry() = 1;
          d_maxcellcount() += DELTACELLCOUNT;
          rand_pool.free_state(rand_gen);
          return;
        }

      }
    }

    if (!jpart && jspecies == ambispecies) {
      //if (nelectron-1 != j-np) memcpy(&d_elist(icell,j-np),&d_elist(icell,nelectron-1),nbytes);
      if (nelectron-1 != j-np) d_elist(icell,j-np) = d_elist(icell,nelectron-1);
      nelectron--;

    } else if (!jpart) {
      int ndelete = Kokkos::atomic_fetch_add(&d_ndelete(),1);
      if (ndelete < d_dellist.extent(0)) {
        d_dellist(ndelete) = d_plist(icell,j);
      } else {
        d_retry() = 1;
        d_maxdelete() += DELTADELETE;
        rand_pool.free_state(rand_gen);
        return;
      }
      d_plist(icell,j) = d_plist(icell,np-1);
      np--;
    }

    // update particle counts
    // quit if no longer enough particles for another collision

    nptotal = np + nelectron;
    if (nptotal < 2) break;
  }

  // done with collisions/chemistry for one grid cell
  // recombine ambipolar ions with their matching electrons
  //   by copying electron velocity into velambi
  // which ion is combined with which electron does not matter
  // error if ion count does not match electron count

  int melectron = 0;
  for (int n = 0; n < np; n++) {
    const int i = d_plist(icell,n);
    if (d_ionambi[i]) {
      if (melectron < nelectron) {
        ep = &d_elist(icell,melectron);
        //memcpy(d_velambi[i],ep->v,3*sizeof(double));
        d_velambi(i,0) = ep->v[0];
        d_velambi(i,1) = ep->v[1];
        d_velambi(i,2) = ep->v[2];
      }
      melectron++;
    }
  }
  if (melectron != nelectron)
    d_error_flag() = 2;

  rand_pool.free_state(rand_gen);
}


/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
double CollideVSSKokkos::attempt_collision_kokkos(int icell, int np, double volume, rand_type &rand_gen) const
{
 double nattempt;

 if (remainflag) {
   nattempt = 0.5 * np * (np-1) *
     d_vremax(icell,0,0) * dt * fnum / volume + d_remain(icell,0,0);
   d_remain(icell,0,0) = nattempt - static_cast<int> (nattempt);
 } else {
   nattempt = 0.5 * np * (np-1) *
     d_vremax(icell,0,0) * dt * fnum / volume + rand_gen.drand();
 }

 // DEBUG
 //nattempt = 10;

  return nattempt;
}

/* ----------------------------------------------------------------------
   determine if collision actually occurs
   1 = yes, 0 = no
   update vremax either way
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
int CollideVSSKokkos::test_collision_kokkos(int icell, int igroup, int jgroup,
                                     Particle::OnePart *ip, Particle::OnePart *jp,
                                     struct State &precoln, rand_type &rand_gen) const
{
  double *vi = ip->v;
  double *vj = jp->v;
  int ispecies = ip->ispecies;
  int jspecies = jp->ispecies;
  double du  = vi[0] - vj[0];
  double dv  = vi[1] - vj[1];
  double dw  = vi[2] - vj[2];
  double vr2 = du*du + dv*dv + dw*dw;
  double vro  = pow(vr2,1.0-d_params(ispecies,jspecies).omega);

  // although the vremax is calcualted for the group,
  // the individual collisions calculated species dependent vre

  double vre = vro*d_prefactor(ispecies,jspecies);
  d_vremax(icell,igroup,jgroup) = MAX(vre,d_vremax(icell,igroup,jgroup));
  if (vre/d_vremax(icell,igroup,jgroup) < rand_gen.drand()) return 0;
  precoln.vr2 = vr2;
  return 1;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::setup_collision_kokkos(Particle::OnePart *ip, Particle::OnePart *jp,
                                       struct State &precoln, struct State &postcoln) const
{
  int isp = ip->ispecies;
  int jsp = jp->ispecies;

  precoln.vr = sqrt(precoln.vr2);

  precoln.ave_rotdof = 0.5 * (d_species[isp].rotdof + d_species[jsp].rotdof);
  precoln.ave_vibdof = 0.5 * (d_species[isp].vibdof + d_species[jsp].vibdof);
  precoln.ave_dof = (precoln.ave_rotdof  + precoln.ave_vibdof)/2.;

  precoln.imass = d_species[isp].mass;
  precoln.jmass = d_species[jsp].mass;

  precoln.etrans = 0.5 * d_params(isp,jsp).mr * precoln.vr2;
  precoln.erot = ip->erot + jp->erot;
  precoln.evib = ip->evib + jp->evib;

  precoln.eint   = precoln.erot + precoln.evib;
  precoln.etotal = precoln.etrans + precoln.eint;

  // COM velocity calculated using reactant masses

  double divisor = 1.0 / (d_species[isp].mass + d_species[jsp].mass);
  double *vi = ip->v;
  double *vj = jp->v;
  precoln.ucmf = ((d_species[isp].mass*vi[0])+(d_species[jsp].mass*vj[0]))*divisor;
  precoln.vcmf = ((d_species[isp].mass*vi[1])+(d_species[jsp].mass*vj[1]))*divisor;
  precoln.wcmf = ((d_species[isp].mass*vi[2])+(d_species[jsp].mass*vj[2]))*divisor;

  postcoln.etrans = precoln.etrans;
  postcoln.erot = 0.0;
  postcoln.evib = 0.0;
  postcoln.eint = 0.0;
  postcoln.etotal = precoln.etotal;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
int CollideVSSKokkos::perform_collision_kokkos(Particle::OnePart *&ip,
                                  Particle::OnePart *&jp,
                                  Particle::OnePart *&kp,
                                  struct State &precoln, struct State &postcoln, rand_type &rand_gen,
                                  Particle::OnePart *&p3, int &recomb_species, double &recomb_density,
                                  int &index_kpart) const
{
  int reactflag,kspecies;
  double x[3],v[3];

  // if gas-phase chemistry defined, attempt and perform reaction
  // if a 3rd particle is created, its kspecies >= 0 is returned
  // if 2nd particle is removed, its jspecies is set to -1

  if (react_defined)
    reactflag = react_kk_copy.obj.attempt_kk(ip,jp,
                                             precoln.etrans,precoln.erot,
                                             precoln.evib,postcoln.etotal,kspecies,
                                             recomb_species,recomb_density,d_species);
  else reactflag = 0;

  // repartition energy and perform velocity scattering for I,J,K particles
  // reaction may have changed species of I,J particles
  // J,K particles may have been removed or created by reaction

  kp = NULL;

  if (reactflag) {

    // add 3rd K particle if reaction created it
    // index of new K particle = nlocal-1
    // if add_particle() performs a realloc:
    //   make copy of x,v, then repoint ip,jp to new particles data struct

    if (kspecies >= 0) {
      int id = MAXSMALLINT*rand_gen.drand();

      memcpy(x,ip->x,3*sizeof(double));
      memcpy(v,ip->v,3*sizeof(double));
      index_kpart = Kokkos::atomic_fetch_add(&d_nlocal(),1);
      int reallocflag =
        ParticleKokkos::add_particle_kokkos(d_particles,index_kpart,id,kspecies,ip->icell,x,v,0.0,0.0);
      if (reallocflag) {
        d_retry() = 1;
        d_part_grow() = 1;
        return 0;
      }

      kp = &d_particles[index_kpart];
      EEXCHANGE_ReactingEDisposal(ip,jp,kp,precoln,postcoln,rand_gen);
      SCATTER_ThreeBodyScattering(ip,jp,kp,precoln,postcoln,rand_gen);

    // remove 2nd J particle if recombination reaction removed it
    // p3 is 3rd particle participating in energy exchange

    } else if (jp->ispecies < 0) {
      double *vi = ip->v;
      double *vj = jp->v;

      const double divisor = 1.0 / (precoln.imass + precoln.jmass);
      const double ucmf = ((precoln.imass*vi[0]) + (precoln.jmass*vj[0])) * divisor;
      const double vcmf = ((precoln.imass*vi[1]) + (precoln.jmass*vj[1])) * divisor;
      const double wcmf = ((precoln.imass*vi[2]) + (precoln.jmass*vj[2])) * divisor;

      vi[0] = ucmf;
      vi[1] = vcmf;
      vi[2] = wcmf;

      jp = NULL;

      // account for 3rd body energy via another call to setup_collision()
      // set precoln.vr2 = relative velocity between ip and 3rd body p3

      const double *vp3 = p3->v;
      const double du  = vi[0] - vp3[0];
      const double dv  = vi[1] - vp3[1];
      const double dw  = vi[2] - vp3[2];
      const double vr2 = du*du + dv*dv + dw*dw;
      precoln.vr2 = vr2;

      // save postcoln.etotal from previous setup_collision()
      // add 3rd body internal energy to it
      // ip internal energy is already included in postcoln.etotal

      double partial_energy =  postcoln.etotal + p3->erot + p3->evib;
      ip->erot = 0.0;
      ip->evib = 0.0;
      p3->erot = 0.0;
      p3->evib = 0.0;

      // 2nd call to setup_collision() sets new postcoln.etotal
      // then add saved partial_energy to it

      setup_collision_kokkos(ip,p3,precoln,postcoln);
      postcoln.etotal += partial_energy;

      if (precoln.ave_dof > 0.0) EEXCHANGE_ReactingEDisposal(ip,p3,jp,precoln,postcoln,rand_gen);
      SCATTER_TwoBodyScattering(ip,p3,precoln,postcoln,rand_gen);

    } else {
      EEXCHANGE_ReactingEDisposal(ip,jp,kp,precoln,postcoln,rand_gen);
      SCATTER_TwoBodyScattering(ip,jp,precoln,postcoln,rand_gen);
    }

  } else {
    kp = NULL;
    if (precoln.ave_dof > 0.0) EEXCHANGE_NonReactingEDisposal(ip,jp,precoln,postcoln,rand_gen);
    SCATTER_TwoBodyScattering(ip,jp,precoln,postcoln,rand_gen);
  }

  return reactflag;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::SCATTER_TwoBodyScattering(Particle::OnePart *ip,
                                                 Particle::OnePart *jp,
                                                 struct State &precoln, struct State &postcoln,
                                                 rand_type &rand_gen) const
{
  double ua,vb,wc;
  double vrc[3];

  double *vi = ip->v;
  double *vj = jp->v;
  int isp = ip->ispecies;
  int jsp = jp->ispecies;
  double mass_i = d_species[isp].mass;
  double mass_j = d_species[jsp].mass;

  double alpha_r = 1.0 / d_params(isp,jsp).alpha;

  double eps = rand_gen.drand() * 2*MY_PI;
  if (fabs(alpha_r - 1.0) < 0.001) {
    double vr = sqrt(2.0 * postcoln.etrans / d_params(isp,jsp).mr);
    double cosX = 2.0*rand_gen.drand() - 1.0;
    double sinX = sqrt(1.0 - cosX*cosX);
    ua = vr*cosX;
    vb = vr*sinX*cos(eps);
    wc = vr*sinX*sin(eps);
  } else {
    double scale = sqrt((2.0 * postcoln.etrans) / (d_params(isp,jsp).mr * precoln.vr2));
    double cosX = 2.0*pow(rand_gen.drand(),alpha_r) - 1.0;
    double sinX = sqrt(1.0 - cosX*cosX);
    vrc[0] = vi[0]-vj[0];
    vrc[1] = vi[1]-vj[1];
    vrc[2] = vi[2]-vj[2];
    double d = sqrt(vrc[1]*vrc[1]+vrc[2]*vrc[2]);
    if (d > 1.0e-6) {
      ua = scale * ( cosX*vrc[0] + sinX*d*sin(eps) );
      vb = scale * ( cosX*vrc[1] + sinX*(precoln.vr*vrc[2]*cos(eps) -
                                         vrc[0]*vrc[1]*sin(eps))/d );
      wc = scale * ( cosX*vrc[2] - sinX*(precoln.vr*vrc[1]*cos(eps) +
                                         vrc[0]*vrc[2]*sin(eps))/d );
    } else {
      ua = scale * ( cosX*vrc[0] );
      vb = scale * ( sinX*vrc[0]*cos(eps) );
      wc = scale * ( sinX*vrc[0]*sin(eps) );
    }
  }

  // new velocities for the products

  double divisor = 1.0 / (mass_i + mass_j);
  vi[0] = precoln.ucmf + (mass_j*divisor)*ua;
  vi[1] = precoln.vcmf + (mass_j*divisor)*vb;
  vi[2] = precoln.wcmf + (mass_j*divisor)*wc;
  vj[0] = precoln.ucmf - (mass_i*divisor)*ua;
  vj[1] = precoln.vcmf - (mass_i*divisor)*vb;
  vj[2] = precoln.wcmf - (mass_i*divisor)*wc;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::EEXCHANGE_NonReactingEDisposal(Particle::OnePart *ip,
                                                      Particle::OnePart *jp,
                                                      struct State &precoln, struct State &postcoln,
                                                      rand_type &rand_gen) const
{
  double State_prob,Fraction_Rot,Fraction_Vib,E_Dispose;
  int i,rotdof,vibdof,max_level,ivib;

  Particle::OnePart *p;

  double AdjustFactor = 0.99999999;
  postcoln.erot = 0.0;
  postcoln.evib = 0.0;
  double pevib = 0.0;

  // handle each kind of energy disposal for non-reacting reactants

  if (precoln.ave_dof == 0) {
    ip->erot = 0.0;
    jp->erot = 0.0;
    ip->evib = 0.0;
    jp->evib = 0.0;

  } else {
    E_Dispose = precoln.etrans;

    for (i = 0; i < 2; i++) {
      if (i == 0) p = ip;
      else p = jp;

      int sp = p->ispecies;
      rotdof = d_species[sp].rotdof;
      double rotn_phi = d_species[sp].rotrel;

      if (rotdof) {
        if (relaxflag == VARIABLE) rotn_phi = rotrel(sp,E_Dispose+p->erot);
        if (rotn_phi >= rand_gen.drand()) {
          if (rotstyle == NONE) {
            p->erot = 0.0 ;

          } else if (rotstyle != NONE && rotdof == 2) {
            E_Dispose += p->erot;
            Fraction_Rot =
              1- pow(rand_gen.drand(),
                     (1/(2.5-d_params(ip->ispecies,jp->ispecies).omega)));
            p->erot = Fraction_Rot * E_Dispose;
            E_Dispose -= p->erot;
          } else {
            E_Dispose += p->erot;
            p->erot = E_Dispose *
              sample_bl(rand_gen,0.5*d_species[sp].rotdof-1.0,
                        1.5-d_params(ip->ispecies,jp->ispecies).omega);
            E_Dispose -= p->erot;
          }
        }
      }
      postcoln.erot += p->erot;

      vibdof = d_species[sp].vibdof;
      double vibn_phi = d_species[sp].vibrel[0];

      if (vibdof) {
        if (relaxflag == VARIABLE) vibn_phi = vibrel(sp,E_Dispose+p->evib);
        if (vibn_phi >= rand_gen.drand()) {
          if (vibstyle == NONE) {
            p->evib = 0.0;

          } else if (vibdof == 2) {
            if (vibstyle == SMOOTH) {
              E_Dispose += p->evib;
              Fraction_Vib =
                1.0 - pow(rand_gen.drand(),(1.0/(2.5-d_params(ip->ispecies,jp->ispecies).omega)));
              p->evib= Fraction_Vib * E_Dispose;
              E_Dispose -= p->evib;

            } else if (vibstyle == DISCRETE) {
              E_Dispose += p->evib;
              max_level = static_cast<int>
                (E_Dispose / (boltz * d_species[sp].vibtemp[0]));
              do {
                ivib = static_cast<int>
                  (rand_gen.drand()*(max_level+AdjustFactor));
                p->evib = ivib * boltz * d_species[sp].vibtemp[0];
                State_prob = pow((1.0 - p->evib / E_Dispose),
                                 (1.5 - d_params(ip->ispecies,jp->ispecies).omega));
              } while (State_prob < rand_gen.drand());
              E_Dispose -= p->evib;
            }
          } else if (vibdof > 2) {
            if (vibstyle == SMOOTH) {
              E_Dispose += p->evib;
              p->evib = E_Dispose *
                sample_bl(rand_gen,0.5*d_species[sp].vibdof-1.0,
                          1.5-d_params(ip->ispecies,jp->ispecies).omega);
              E_Dispose -= p->evib;

            } else if (vibstyle == DISCRETE) {
              p->evib = 0.0;

              int nmode = d_species[sp].nvibmode;
              auto &d_vibmode = k_eiarray.d_view[d_ewhich[index_vibmode]].k_view.d_view;
              int pindex = p - d_particles.data();

              for (int imode = 0; imode < nmode; imode++) {
                ivib = d_vibmode(pindex,imode);
                E_Dispose += ivib * boltz *
                  d_species[sp].vibtemp[imode];
                max_level = static_cast<int>
                  (E_Dispose / (boltz * d_species[sp].vibtemp[imode]));

                do {
                  ivib = static_cast<int>
                    (rand_gen.drand()*(max_level+AdjustFactor));
                  pevib = ivib * boltz * d_species[sp].vibtemp[imode];
                  State_prob = pow((1.0 - pevib / E_Dispose),
                                   (1.5 - d_params(ip->ispecies,jp->ispecies).omega));
                } while (State_prob < rand_gen.drand());

                d_vibmode(pindex,imode) = ivib;
                p->evib += pevib;
                E_Dispose -= pevib;
              }
            }
          } // end of vibstyle/vibdof if
        }
        postcoln.evib += p->evib;
      } // end of vibdof if
    }
  }

  // compute portion of energy left over for scattering

  postcoln.eint = postcoln.erot + postcoln.evib;
  postcoln.etrans = E_Dispose;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::SCATTER_ThreeBodyScattering(Particle::OnePart *ip,
                                                   Particle::OnePart *jp,
                                                   Particle::OnePart *kp,
                                                   struct State &precoln, struct State &postcoln,
                                                   rand_type &rand_gen) const
{
  double vrc[3],ua,vb,wc;

  int isp = ip->ispecies;
  int jsp = jp->ispecies;
  int ksp = kp->ispecies;
  double mass_i = d_species[isp].mass;
  double mass_j = d_species[jsp].mass;
  double mass_k = d_species[ksp].mass;
  double mass_ij = mass_i + mass_j;
  double *vi = ip->v;
  double *vj = jp->v;
  double *vk = kp->v;

  double alpha_r = 1.0 / d_params(isp,jsp).alpha;
  double mr = mass_ij * mass_k / (mass_ij + mass_k);
  postcoln.eint = ip->erot + jp->erot + ip->evib + jp->evib
                + kp->erot + kp->evib;

  double cosX = 2.0*pow(rand_gen.drand(), alpha_r) - 1.0;
  double sinX = sqrt(1.0 - cosX*cosX);
  double eps = rand_gen.drand() * 2*MY_PI;

  if (fabs(alpha_r - 1.0) < 0.001) {
    double vr = sqrt(2*postcoln.etrans/mr);
    ua = vr*cosX;
    vb = vr*sinX*cos(eps);
    wc = vr*sinX*sin(eps);
  } else {
    double scale = sqrt((2.0*postcoln.etrans) / (mr*precoln.vr2));
    vrc[0] = vi[0]-vj[0];
    vrc[1] = vi[1]-vj[1];
    vrc[2] = vi[2]-vj[2];
    double d = sqrt(vrc[1]*vrc[1]+vrc[2]*vrc[2]);
    if (d > 1.E-6 ) {
      ua = scale * (cosX*vrc[0] + sinX*d*sin(eps));
      vb = scale * (cosX*vrc[1] + sinX*(precoln.vr*vrc[2]*cos(eps) -
                                        vrc[0]*vrc[1]*sin(eps))/d);
      wc = scale * (cosX*vrc[2] - sinX*(precoln.vr*vrc[1]*cos(eps) +
                                        vrc[0]*vrc[2]*sin(eps))/d);
    } else {
      ua = scale * cosX*vrc[0];
      vb = scale * sinX*vrc[0]*cos(eps);
      wc = scale * sinX*vrc[0]*sin(eps);
    }
  }

  // new velocities for the products

  double divisor = 1.0 / (mass_ij + mass_k);
  vi[0] = precoln.ucmf + (mass_k*divisor)*ua;
  vi[1] = precoln.vcmf + (mass_k*divisor)*vb;
  vi[2] = precoln.wcmf + (mass_k*divisor)*wc;
  vk[0] = precoln.ucmf - (mass_ij*divisor)*ua;
  vk[1] = precoln.vcmf - (mass_ij*divisor)*vb;
  vk[2] = precoln.wcmf - (mass_ij*divisor)*wc;
  vj[0] = vi[0];
  vj[1] = vi[1];
  vj[2] = vi[2];
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::EEXCHANGE_ReactingEDisposal(Particle::OnePart *ip,
                                                   Particle::OnePart *jp,
                                                   Particle::OnePart *kp,
                                                   struct State &precoln, struct State &postcoln,
                                                   rand_type &rand_gen) const
{
  double State_prob,Fraction_Rot,Fraction_Vib;
  int i,numspecies,rotdof,vibdof,max_level,ivib;
  double aveomega,pevib;

  Particle::OnePart *p;
  double AdjustFactor = 0.99999999;

  if (!kp) {
    ip->erot = 0.0;
    jp->erot = 0.0;
    ip->evib = 0.0;
    jp->evib = 0.0;
    numspecies = 2;
    aveomega = d_params(ip->ispecies,jp->ispecies).omega;
  } else {
    ip->erot = 0.0;
    jp->erot = 0.0;
    kp->erot = 0.0;
    ip->evib = 0.0;
    jp->evib = 0.0;
    kp->evib = 0.0;
    numspecies = 3;
    aveomega = (d_params(ip->ispecies,ip->ispecies).omega + d_params(jp->ispecies,jp->ispecies).omega +
                d_params(kp->ispecies,kp->ispecies).omega)/3.0;
  }

  // handle each kind of energy disposal for non-reacting reactants
  // clean up memory for the products

  double E_Dispose = postcoln.etotal;

  for (i = 0; i < numspecies; i++) {
    if (i == 0) p = ip;
    else if (i == 1) p = jp;
    else p = kp;

    int sp = p->ispecies;
    rotdof = d_species[sp].rotdof;

    if (rotdof) {
      if (rotstyle == NONE) {
        p->erot = 0.0 ;
      } else if (rotdof == 2) {
        Fraction_Rot =
          1- pow(rand_gen.drand(),(1/(2.5-aveomega)));
        p->erot = Fraction_Rot * E_Dispose;
        E_Dispose -= p->erot;

      } else if (rotdof > 2) {
        p->erot = E_Dispose *
          sample_bl(rand_gen,0.5*d_species[sp].rotdof-1.0,
                    1.5-aveomega);
        E_Dispose -= p->erot;
      }
    }

    vibdof = d_species[sp].vibdof;

    if (vibdof) {
      if (vibstyle == NONE) {
        p->evib = 0.0;
      } else if (vibdof == 2 && vibstyle == DISCRETE) {
        max_level = static_cast<int>
          (E_Dispose / (boltz * d_species[sp].vibtemp[0]));
        do {
          ivib = static_cast<int>
            (rand_gen.drand()*(max_level+AdjustFactor));
          p->evib = (double)
            (ivib * boltz * d_species[sp].vibtemp[0]);
          State_prob = pow((1.0 - p->evib / E_Dispose),
                           (1.5 - aveomega));
        } while (State_prob < rand_gen.drand());
        E_Dispose -= p->evib;

      } else if (vibdof == 2 && vibstyle == SMOOTH) {
        Fraction_Vib =
          1.0 - pow(rand_gen.drand(),(1.0 / (2.5-aveomega)));
        p->evib = Fraction_Vib * E_Dispose;
        E_Dispose -= p->evib;

      } else if (vibdof > 2 && vibstyle == SMOOTH) {
        p->evib = E_Dispose *
          sample_bl(rand_gen,0.5*d_species[sp].vibdof-1.0,
                    1.5-aveomega);
        E_Dispose -= p->evib;
      } else if (vibdof > 2 && vibstyle == DISCRETE) {
        p->evib = 0.0;

        int nmode = d_species[sp].nvibmode;
        auto &d_vibmode = k_eiarray.d_view[d_ewhich[index_vibmode]].k_view.d_view;
        int pindex = p - d_particles.data();

        for (int imode = 0; imode < nmode; imode++) {
          ivib = d_vibmode(pindex,imode);
          E_Dispose += ivib * boltz *
          d_species[sp].vibtemp[imode];
          max_level = static_cast<int>
          (E_Dispose / (boltz * d_species[sp].vibtemp[imode]));
          do {
            ivib = static_cast<int>
            (rand_gen.drand()*(max_level+AdjustFactor));
            pevib = ivib * boltz * d_species[sp].vibtemp[imode];
            State_prob = pow((1.0 - pevib / E_Dispose),
                             (1.5 - aveomega));
          } while (State_prob < rand_gen.drand());

          d_vibmode(pindex,imode) = ivib;
          p->evib += pevib;
          E_Dispose -= pevib;
        }
      }
    }
  }

  // compute post-collision internal energies

  postcoln.erot = ip->erot + jp->erot;
  postcoln.evib = ip->evib + jp->evib;

  if (kp) {
    postcoln.erot += kp->erot;
    postcoln.evib += kp->evib;
  }

  // compute portion of energy left over for scattering

  postcoln.eint = postcoln.erot + postcoln.evib;
  postcoln.etrans = E_Dispose;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
double CollideVSSKokkos::sample_bl(rand_type &rand_gen, double Exp_1, double Exp_2) const
{
  double Exp_s = Exp_1 + Exp_2;
  double x,y;
  do {
    x = rand_gen.drand();
    y = pow(x*Exp_s/Exp_1, Exp_1)*pow((1.0-x)*Exp_s/Exp_2, Exp_2);
  } while (y < rand_gen.drand());
  return x;
}

/* ----------------------------------------------------------------------
   compute a variable rotational relaxation parameter
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
double CollideVSSKokkos::rotrel(int isp, double Ec) const
{
  // Because we are only relaxing one of the particles in each call, we only
  //  include its DoF, consistent with Bird 2013 (3.32)

  double Tr = Ec /(boltz * (2.5-d_params(isp,isp).omega + d_species[isp].rotdof/2.0));
  double rotphi = (1.0+d_params(isp,isp).rotc2/sqrt(Tr) + d_params(isp,isp).rotc3/Tr)
                / d_params(isp,isp).rotc1;
  return rotphi;
}

/* ----------------------------------------------------------------------
   compute a variable vibrational relaxation parameter
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
double CollideVSSKokkos::vibrel(int isp, double Ec) const
{
  double Tr = Ec /(boltz * (3.5-d_params(isp,isp).omega));
  double omega = d_params(isp,isp).omega;
  double vibphi = 1.0 / (d_params(isp,isp).vibc1/pow(Tr,omega) *
                         exp(d_params(isp,isp).vibc2/pow(Tr,1.0/3.0)));
  return vibphi;
}

/* ----------------------------------------------------------------------
   for particle I, find collision partner J via near neighbor algorithm
   always returns a J neighbor, even if not that near
   near neighbor algorithm:
     check up to nearlimit particles, starting with random particle
     as soon as find one within distance moved by particle I, return it
     else return the closest one found
     also exclude an I,J pair if both most recently collided with each other
   this version is for single group collisions
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
int CollideVSSKokkos::find_nn(rand_type &rand_gen, int i, int np, int icell) const
{
  int jneigh;
  double dx,dy,dz,rsq;
  double *xj;

  // if np = 2, just return J = non-I particle
  // np is never < 2

  if (np == 2) return (i+1) % 2;

  Particle::OnePart *ipart,*jpart;

  // thresh = distance particle I moves in this timestep

  ipart = &d_particles[d_plist(icell,i)];
  double *vi = ipart->v;
  double *xi = ipart->x;
  double threshsq =  dt*dt * (vi[0]*vi[0]+vi[1]*vi[1]+vi[2]*vi[2]);
  double minrsq = BIG;

  // nlimit = max # of J candidates to consider

  int nlimit = MIN(nearlimit,np-1);
  int count = 0;

  // pick a random starting J
  // jneigh = collision partner when exit loop
  //   set to initial J as default in case no Nlimit J meets criteria

  int j = np * rand_gen.drand();
  while (i == j) j = np * rand_gen.drand();
  jneigh = j;

  while (count < nlimit) {
    count++;

    // skip this J if I,J last collided with each other

    if (d_nn_last_partner(icell,i) == j+1 && d_nn_last_partner(icell,j) == i+1) {
      j++;
      if (j == np) j = 0;
      continue;
    }

    // rsq = squared distance between particles I and J
    // if rsq = 0.0, skip this J
    //   could be I = J, or a cloned J at same position as I
    // if rsq <= threshsq, this J is collision partner
    // if rsq = smallest yet seen, this J is tentative collision partner

    jpart = &d_particles[d_plist(icell,j)];
    xj = jpart->x;
    dx = xi[0] - xj[0];
    dy = xi[1] - xj[1];
    dz = xi[2] - xj[2];
    rsq = dx*dx + dy*dy + dz*dz;

    if (rsq > 0.0) {
      if (rsq <= threshsq) {
        jneigh = j;
        break;
      }
      if (rsq < minrsq) {
        minrsq = rsq;
        jneigh = j;
      }
    }
    j++;
    if (j == np) j = 0;
  }

  return jneigh;
}

/* ----------------------------------------------------------------------
   reset ionambi flags if ambipolar reaction occurred
   this operates independent of cell particle counts and plist/elist data structs
     caller will adjust those after this method returns
   i/j = indices of I,J reactants
   isp/jsp = pre-reaction species of I,J
     both will not be electrons, if one is electron it will be jsp
   reactants i,j and isp/jsp will always be in order listed below
   products ip,jp,kp will always be in order listed below
   logic must be valid for all ambipolar AND non-ambipolar reactions
   check for 3 versions of 2 -> 3: dissociation or ionization
     all have J product = electron
     D: AB + e -> A + e + B
        if I reactant = neutral and K product not electron:
        set K product = neutral
     D: AB+ + e -> A+ + e + B
        if I reactant = ion:
        set K product = neutral
     I: A + e -> A+ + e + e
        if I reactant = neutral and K product = electron:
        set I product = ion
     all other 2 -> 3 cases, set K product = neutral
   check for 4 versions of 2 -> 2: ionization or exchange
     I: A + B -> AB+ + e
        if J product = electron:
        set I product to ion
     E: AB+ + e -> A + B
        if I reactant = ion and J reactant = elecrton
        set I/J products to neutral
     E: AB+ + C -> A + BC+
        if I reactant = ion:
        set I/J products to neutral/ion
     E: C + AB+ -> A + BC+
        if J reactant = ion:
        nothing to change for products
     all other 2 -> 2 cases, no changes
   check for one version of 2 -> 1: recombination
     R: A+ + e -> A
        if ej = elec, set I product to neutral
     all other 2 -> 1 cases, no changes
   WARNING:
     do not index by I,J if could be e, since may be negative I,J index
     do not access ionambi if could be e, since e may be in elist
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::ambi_reset_kokkos(int i, int j, int jsp, int index_kpart,
                                      Particle::OnePart *ip, Particle::OnePart *jp,
                                      Particle::OnePart *kp, const DAT::t_int_1d &d_ionambi) const
{
  int e = ambispecies;

  // 2 reactants become 3 products
  // in all ambi reactions with an electron reactant, it is J

  if (kp) {
    int k = index_kpart;
    d_ionambi[k] = 0;
    if (jsp != e) return;

    if (d_ionambi[i]) {                // nothing to change
    } else if (kp->ispecies == e) {
      d_ionambi[i] = 1;                // 1st reactant is now 1st product ion
    }

  // 2 reactants become 2 products
  // ambi reaction if J product is electron or either reactant is ion

  } else if (jp) {
    if (jp->ispecies == e) {
      d_ionambi[i] = 1;         // 1st reactant is now 1st product ion
    } else if (d_ionambi[i] && jsp == e) {
      d_ionambi[i] = 0;         // 1st reactant is now 1st product neutral
    } else if (d_ionambi[i]) {
      d_ionambi[i] = 0;         // 1st reactant is now 1st product neutral
      d_ionambi[j] = 1;         // 2nd reactant is now 2nd product ion
    }

  // 2 reactants become 1 product
  // ambi reaction if J reactant is electron

  } else if (!jp) {
    if (jsp == e) d_ionambi[i] = 0;   // 1st reactant is now 1st product neutral
  }
}

/* ----------------------------------------------------------------------
   pack icell values for per-cell arrays into buf
   if icell is a split cell, also pack all sub cell values
   return byte count of amount packed
   if memflag, only return count, do not fill buf
   NOTE: why packing/unpacking parent cell if a split cell?
------------------------------------------------------------------------- */

int CollideVSSKokkos::pack_grid_one(int icell, char *buf_char, int memflag)
{
  double* buf = (double*) buf_char;

  Grid::ChildCell *cells = grid->cells;

  int n = 0;
  if (memflag) {
    for (int igroup = 0; igroup < ngroups; igroup++) {
      for (int jgroup = 0; jgroup < ngroups; jgroup++) {
        buf[n++] = k_vremax.h_view(icell,igroup,jgroup);
        if (remainflag)
          buf[n++] = k_remain.h_view(icell,igroup,jgroup);
      }
    }
  } else {
    n += ngroups*ngroups;
    if (remainflag)
      n += ngroups*ngroups;
  }

  if (cells[icell].nsplit > 1) {
    int isplit = cells[icell].isplit;
    int nsplit = cells[icell].nsplit;
    for (int i = 0; i < nsplit; i++) {
      int m = grid->sinfo[isplit].csubs[i];
      if (memflag) {
        for (int igroup = 0; igroup < ngroups; igroup++) {
          for (int jgroup = 0; jgroup < ngroups; jgroup++) {
            buf[n++] = k_vremax.h_view(m,igroup,jgroup);
            if (remainflag)
              buf[n++] = k_remain.h_view(m,igroup,jgroup);
          }
        }
      } else {
        n += ngroups*ngroups;
        if (remainflag)
          n += ngroups*ngroups;
      }
    }
  }

  return n*sizeof(double);
}

/* ----------------------------------------------------------------------
   unpack icell values for per-cell arrays from buf
   if icell is a split cell, also unpack all sub cell values
   return byte count of amount unpacked
------------------------------------------------------------------------- */

int CollideVSSKokkos::unpack_grid_one(int icell, char *buf_char)
{
  double* buf = (double*) buf_char;

  Grid::ChildCell *cells = grid->cells;
  Grid::SplitInfo *sinfo = grid->sinfo;

  grow_percell(1);
  int n = 0;
  for (int igroup = 0; igroup < ngroups; igroup++) {
    for (int jgroup = 0; jgroup < ngroups; jgroup++) {
      k_vremax.h_view(icell,igroup,jgroup) = buf[n++];
      if (remainflag)
        k_remain.h_view(icell,igroup,jgroup) = buf[n++];
    }
  }
  nglocal++;

  if (cells[icell].nsplit > 1) {
    int isplit = cells[icell].isplit;
    int nsplit = cells[icell].nsplit;
    grow_percell(nsplit);
    for (int i = 0; i < nsplit; i++) {
      int m = sinfo[isplit].csubs[i];
      for (int igroup = 0; igroup < ngroups; igroup++) {
        for (int jgroup = 0; jgroup < ngroups; jgroup++) {
          k_vremax.h_view(m,igroup,jgroup) = buf[n++];
          if (remainflag)
            k_remain.h_view(m,igroup,jgroup) = buf[n++];
        }
      }
    }
    nglocal += nsplit;
  }

  return n*sizeof(double);
}

/* ----------------------------------------------------------------------
   copy per-cell collision info from Icell to Jcell
   called whenever a grid cell is removed from this processor's list
   caller checks that Icell != Jcell
------------------------------------------------------------------------- */

void CollideVSSKokkos::copy_grid_one(int icell, int jcell)
{
  this->sync(Host,ALL_MASK);
  for (int igroup = 0; igroup < ngroups; igroup++) {
    for (int jgroup = 0; jgroup < ngroups; jgroup++) {
      k_vremax.h_view(jcell,igroup,jgroup) = k_vremax.h_view(icell,igroup,jgroup);
      if (remainflag)
        k_remain.h_view(jcell,igroup,jgroup) = k_remain.h_view(icell,igroup,jgroup);
    }
  }
  this->modified(Host,ALL_MASK);
}

/* ----------------------------------------------------------------------
   reset final grid cell count after grid cell removals
------------------------------------------------------------------------- */

void CollideVSSKokkos::reset_grid_count(int nlocal)
{
  nglocal = nlocal;
}

/* ----------------------------------------------------------------------
   add a grid cell
   called when a grid cell is added to this processor's list
   initialize values to 0.0
------------------------------------------------------------------------- */

void CollideVSSKokkos::add_grid_one()
{
  grow_percell(1);

  this->sync(Host,ALL_MASK);
  for (int igroup = 0; igroup < ngroups; igroup++)
    for (int jgroup = 0; jgroup < ngroups; jgroup++) {
      k_vremax.h_view(nglocal,igroup,jgroup) = vremax_initial[igroup][jgroup];
      if (remainflag) k_remain.h_view(nglocal,igroup,jgroup) = 0.0;
    }
  this->modified(Host,ALL_MASK);

  nglocal++;
}

/* ----------------------------------------------------------------------
   reinitialize per-cell arrays due to grid cell adaptation
   count of owned grid cells has changed
   called from adapt_grid
------------------------------------------------------------------------- */

void CollideVSSKokkos::adapt_grid()
{
  int nglocal_old = nglocal;
  nglocal = grid->nlocal;

  // reallocate vremax and remain
  // initialize only new added locations
  // this leaves vremax/remain for non-adapted cells the same

  this->sync(Host,ALL_MASK);
  this->modified(Host,ALL_MASK); // force resize on host

  nglocalmax = nglocal;
  k_vremax.resize(nglocalmax,ngroups,ngroups);
  d_vremax = k_vremax.d_view;
  if (remainflag) {
    k_remain.resize(nglocalmax,ngroups,ngroups);
    d_remain = k_remain.d_view;
  }
  this->sync(Host,ALL_MASK);
  for (int icell = nglocal_old; icell < nglocal; icell++)
    for (int igroup = 0; igroup < ngroups; igroup++)
      for (int jgroup = 0; jgroup < ngroups; jgroup++) {
        k_vremax.h_view(icell,igroup,jgroup) = vremax_initial[igroup][jgroup];
        if (remainflag) k_remain.h_view(icell,igroup,jgroup) = 0.0;
      }

  this->modified(Host,ALL_MASK);
}

/* ----------------------------------------------------------------------
   insure per-cell arrays are allocated long enough for N new cells
------------------------------------------------------------------------- */

void CollideVSSKokkos::grow_percell(int n)
{
  if (nglocal+n < nglocalmax || !ngroups) return;
  while (nglocal+n >= nglocalmax) nglocalmax += DELTAGRID;

  this->sync(Device,ALL_MASK); // force resize on device

  k_vremax.resize(nglocalmax,ngroups,ngroups);
  d_vremax = k_vremax.d_view;
  if (remainflag) {
    k_remain.resize(nglocalmax,ngroups,ngroups);
    d_remain = k_remain.d_view;
  }

  this->modified(Device,ALL_MASK); // needed for auto sync
}

/* ---------------------------------------------------------------------- */

void CollideVSSKokkos::sync(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if (sparta->kokkos->auto_sync)
      modified(Host,mask);
    if (mask & VREMAX_MASK) k_vremax.sync_device();
    if (remainflag)
      if (mask & REMAIN_MASK) k_remain.sync_device();
  } else {
    if (mask & VREMAX_MASK) k_vremax.sync_host();
    if (remainflag)
      if (mask & REMAIN_MASK) k_remain.sync_host();
  }
}

/* ---------------------------------------------------------------------- */

void CollideVSSKokkos::modified(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if (mask & VREMAX_MASK) k_vremax.modify_device();
    if (remainflag)
      if (mask & REMAIN_MASK) k_remain.modify_device();
    if (sparta->kokkos->auto_sync)
      sync(Host,mask);
  } else {
    if (mask & VREMAX_MASK) k_vremax.modify_host();
    if (remainflag)
      if (mask & REMAIN_MASK) k_remain.modify_host();
  }
}

/* ---------------------------------------------------------------------- */

void CollideVSSKokkos::backup()
{
  d_particles_backup = decltype(d_particles)(Kokkos::view_alloc("collide:particles_backup",Kokkos::WithoutInitializing),d_particles.extent(0));
  d_plist_backup = decltype(d_plist)(Kokkos::view_alloc("collide:plist_backup",Kokkos::WithoutInitializing),d_plist.extent(0),d_plist.extent(1));
  d_vremax_backup = decltype(d_vremax)(Kokkos::view_alloc("collide:vremax_backup",Kokkos::WithoutInitializing),d_vremax.extent(0),d_vremax.extent(1),d_vremax.extent(2));
  d_remain_backup = decltype(d_remain)(Kokkos::view_alloc("collide:remain_backup",Kokkos::WithoutInitializing),d_remain.extent(0),d_remain.extent(1),d_remain.extent(2));

  if (ambiflag) {
    d_ionambi_backup = decltype(d_ionambi)(Kokkos::view_alloc("collide:ionambi_backup",Kokkos::WithoutInitializing),d_ionambi.extent(0));
    d_velambi_backup = decltype(d_velambi)(Kokkos::view_alloc("collide:velambi_backup",Kokkos::WithoutInitializing),d_velambi.extent(0),d_velambi.extent(1));
  }

  Kokkos::deep_copy(d_particles_backup,d_particles);
  Kokkos::deep_copy(d_plist_backup,d_plist);
  Kokkos::deep_copy(d_vremax_backup,d_vremax);
  Kokkos::deep_copy(d_remain_backup,d_remain);

  if (ambiflag) {
    Kokkos::deep_copy(d_ionambi_backup,d_ionambi);
    Kokkos::deep_copy(d_velambi_backup,d_velambi);
  }

  if (react) {
    ReactBirdKokkos* react_kk = (ReactBirdKokkos*) react;
    react_kk->backup();
  }

#ifdef SPARTA_KOKKOS_EXACT
  if (!random_backup)
    random_backup = new RanKnuth(12345 + comm->me);
  memcpy(random_backup,random,sizeof(RanKnuth));
#endif

}

/* ---------------------------------------------------------------------- */

void CollideVSSKokkos::restore()
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  Kokkos::deep_copy(particle_kk->k_particles.d_view,d_particles_backup);
  d_particles = particle_kk->k_particles.d_view;

  GridKokkos* grid_kk = (GridKokkos*) grid;
  Kokkos::deep_copy(grid_kk->d_plist,d_plist_backup);
  d_plist = grid_kk->d_plist;

  Kokkos::deep_copy(d_vremax,d_vremax_backup);
  Kokkos::deep_copy(d_remain,d_remain_backup);

  if (ambiflag) {
    auto h_ewhich = particle_kk->k_ewhich.h_view;

    Kokkos::deep_copy(particle_kk->k_eivec.h_view[h_ewhich[index_ionambi]].k_view.d_view,d_ionambi_backup);
    Kokkos::deep_copy(particle_kk->k_edarray.h_view[h_ewhich[index_velambi]].k_view.d_view,d_velambi_backup);

    k_eivec = particle_kk->k_eivec;
    k_edarray = particle_kk->k_edarray;
    d_ionambi = k_eivec.h_view[h_ewhich[index_ionambi]].k_view.d_view;
    d_velambi = k_edarray.h_view[h_ewhich[index_velambi]].k_view.d_view;
  }

  if (react) {
    ReactBirdKokkos* react_kk = (ReactBirdKokkos*) react;
    react_kk->restore();
  }

#ifdef SPARTA_KOKKOS_EXACT
  memcpy(random,random_backup,sizeof(RanKnuth));
#endif

  //  reset counters

  if (sparta->kokkos->atomic_reduction) {
    h_nattempt_one() = 0;
    h_ncollide_one() = 0;
    h_nreact_one() = 0;
  }

  // deallocate references to reduce memory use

  d_particles_backup = decltype(d_particles_backup)();
  d_plist_backup = decltype(d_plist_backup)();
  d_vremax_backup = decltype(d_vremax_backup)();
  d_remain_backup = decltype(d_remain_backup)();

  if (ambiflag) {
    d_ionambi_backup = decltype(d_ionambi_backup)();
    d_velambi_backup = decltype(d_velambi_backup)();
  }
}

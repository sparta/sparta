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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "collide_vss_kokkos.h"
#include "grid.h"
#include "domain.h"
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
#include "fix_ambipolar_kokkos.h"

using namespace SPARTA_NS;
using namespace MathConst;

#define VAL_1(X) X
#define VAL_2(X) VAL_1(X), VAL_1(X)
#define VAL_4(X) VAL_2(X), VAL_2(X)
// blit one active gas tally compute into its per-type device buffer
// same operation and rationale as KKCopy::copy() (kokkos_copy.h:71): the
//   object is only read on device, through KOKKOS_INLINE_FUNCTION members,
//   so its vtable pointer is never used and the View handles it carries stay
//   alive in the original that update->glist_active holds

#ifndef SPARTA_KOKKOS_FIXED_LISTS
namespace {

  template<class T>
  void gas_buf_resize(DAT::tdual_char_1d &k, DAT::t_char_1d &d, int n)
  {
    const size_t need = (size_t) MAX(n,1) * sizeof(T);
    if (k.view_device().extent(0) < need) {
      k = DAT::tdual_char_1d("collide:gas_tally_models",need);
      d = k.view_device();
    }
  }

  template<class T>
  void gas_buf_blit(DAT::tdual_char_1d &k, int slot, T *obj)
  {
    char *dst = k.view_host().data() + (size_t) slot*sizeof(T);
    memcpy((void*) dst, (const void*) obj, sizeof(T));
    ((T *) dst)->copy = 1;
  }

  void gas_buf_sync(DAT::tdual_char_1d &k, DAT::t_char_1d &d)
  {
    if (k.view_device().extent(0) == 0) return;
    k.modify_host();
    k.sync_device();
    d = k.view_device();
  }
}
#endif

enum{NONE,DISCRETE,SMOOTH};            // several files
enum{CONSTANT,VARIABLE};

#define DELTAGRID 1000            // must be bigger than split cells per cell
#define DELTADELETE 1024
#define DELTAELECTRON 128
#define DELTACELLCOUNT 2

#define EPSZERO 1.0e-14
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
  react_kk_copy(sparta),
  react_qk_kk_copy(sparta),
  react_tceqk_kk_copy(sparta)
#ifdef SPARTA_KOKKOS_FIXED_LISTS
  , glist_collision_copy{VAL_4(KKCopy<ComputeGasCollisionGridKokkos>(sparta))}
  , glist_coll_tally_copy{VAL_4(KKCopy<ComputeGasCollisionTallyKokkos>(sparta))}
  , glist_react_tally_copy{VAL_4(KKCopy<ComputeGasReactionTallyKokkos>(sparta))}
  , glist_reaction_copy{VAL_4(KKCopy<ComputeGasReactionGridKokkos>(sparta))}
  , tmp_compute_gas_collision_kk(sparta)
  , tmp_compute_gas_reaction_kk(sparta)
  , tmp_compute_gas_coll_tally_kk(sparta)
  , tmp_compute_gas_react_tally_kk(sparta)
#endif
{
  kokkos_flag = 1;
  react_style = 0;
  nglist_collision = nglist_reaction = 0;
  nglist_coll_tally = nglist_react_tally = 0;
  egroup = -1;

  // use 1D view for scalars to reduce GPU memory operations

  // int view = flags and view-size counters, must stay int
  // bigint view = per-step statistics counters, can exceed 2^31
  //   in one step at large per-proc particle counts

  d_scalars = t_int_8("collide:scalars");
  h_scalars = t_host_int_8("collide:scalars_mirror");

  d_scalars_big = t_bigint_3("collide:scalars_big");
  h_scalars_big = t_host_bigint_3("collide:scalars_big_mirror");

  d_error_flag   = Kokkos::subview(d_scalars,0);
  d_retry        = Kokkos::subview(d_scalars,1);
  d_maxdelete    = Kokkos::subview(d_scalars,2);
  d_maxcellcount = Kokkos::subview(d_scalars,3);
  d_part_grow    = Kokkos::subview(d_scalars,4);
  d_ndelete      = Kokkos::subview(d_scalars,5);
  d_nlocal       = Kokkos::subview(d_scalars,6);
  d_maxelectron  = Kokkos::subview(d_scalars,7);
  d_tally_overflow = Kokkos::subview(d_scalars,8);

  d_nattempt_one = Kokkos::subview(d_scalars_big,0);
  d_ncollide_one = Kokkos::subview(d_scalars_big,1);
  d_nreact_one   = Kokkos::subview(d_scalars_big,2);

  h_error_flag   = Kokkos::subview(h_scalars,0);
  h_retry        = Kokkos::subview(h_scalars,1);
  h_maxdelete    = Kokkos::subview(h_scalars,2);
  h_maxcellcount = Kokkos::subview(h_scalars,3);
  h_part_grow    = Kokkos::subview(h_scalars,4);
  h_ndelete      = Kokkos::subview(h_scalars,5);
  h_nlocal       = Kokkos::subview(h_scalars,6);
  h_maxelectron  = Kokkos::subview(h_scalars,7);
  h_tally_overflow = Kokkos::subview(h_scalars,8);

  h_nattempt_one = Kokkos::subview(h_scalars_big,0);
  h_ncollide_one = Kokkos::subview(h_scalars_big,1);
  h_nreact_one   = Kokkos::subview(h_scalars_big,2);

  random_backup = NULL;
  react_defined = 0;

  maxdelete = DELTADELETE;
}

/* ---------------------------------------------------------------------- */

CollideVSSKokkos::~CollideVSSKokkos()
{
  if (copymode) return;

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

  // CollideVSSKokkos::init() does not call the host base, so it carries its
  //   own copy of the host's restriction checks.  These three mirror
  //   Collide::init() (collide.cpp:166-176) condition-for-condition and
  //   message-for-message: the ambipolar model has no near-neighbor or subcell
  //   implementation on the CPU either, so they are host restrictions being
  //   reproduced, not Kokkos limitations.

  if (ambiflag && nearcp)
    error->all(FLERR,"Ambipolar collision model does not yet support "
               "near-neighbor collisions");

  if (ambiflag && subcellflag)
    error->all(FLERR,"Ambipolar collision model does not yet support "
               "subcell collisions");

  if (nearcp && subcellflag)
    error->all(FLERR,"Cannot use both nearcp and subcell collision partners");

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
      snprintf(str,sizeof(str),"%d species do not define correct vibrational "
              "modes for discrete model",flag);
      error->all(FLERR,str);
    }
  }

  if (elecstyle == DISCRETE) {
    index_elecstate = particle->find_custom((char *) "elecstate");
    index_eelec = particle->find_custom((char *) "eelec");

    if (index_elecstate < 0 || index_eelec < 0) {
        error->all(FLERR,
                   "Fix elecmode must be used with discrete electronic modes");
    }

    // rebuild the flattened electronic-data views if they are stale,
    // e.g. species read from a restart file without a new species command;
    // without this the views are zero-length and indexed out of bounds

    ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
    if ((int)particle_kk->d_nelecstates.extent(0) != particle->nspecies)
      particle_kk->update_elec_views();
  }

  // reallocate one-cell data structs for one or many groups

  oldgroups = ngroups;
  ngroups = mixture->ngroup;

  // must follow ngroups assignment above, as in Collide::init()
  // mirrors collide.cpp:274-276; the host has no multigroup subcell algorithm
  //   either (Collide::subcell_alloc() refuses to allocate at collide.cpp:988),
  //   so there is no host result for a Kokkos version to reproduce

  if (subcellflag && ngroups > 1)
    error->all(FLERR,"Cannot yet use subcell collisions with "
               "multiple collision groups");

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
    MemKK::realloc_kokkos(k_vremax,"collide:vremax",nglocalmax,ngroups,ngroups);
    d_vremax = k_vremax.view_device();
    MemKK::realloc_kokkos(k_remain,"collide:remain",nglocalmax,ngroups,ngroups);
    d_remain = k_remain.view_device();

    for (int igroup = 0; igroup < ngroups; igroup++) {
      for (int jgroup = 0; jgroup < ngroups; jgroup++) {
        vremax_initial[igroup][jgroup] = vremax_init(igroup,jgroup);
        k_vremax_initial.view_host()(igroup,jgroup) = vremax_initial[igroup][jgroup];
      }
    }

    k_vremax_initial.modify_host();
    k_vremax_initial.sync_device();
    d_vremax_initial = k_vremax_initial.view_device();
  }

  // device copy of species-to-group mapping for group collisions

  if (ngroups > 1) {
    int nspecies = particle->nspecies;
    int *species2group = mixture->species2group;
    d_species2group = DAT::t_int_1d("collide:species2group",nspecies);
    auto h_species2group = Kokkos::create_mirror_view(d_species2group);
    for (int i = 0; i < nspecies; i++) h_species2group(i) = species2group[i];
    Kokkos::deep_copy(d_species2group,h_species2group);
  }

  // if recombination reactions exist, set flags per species pair

  recombflag = 0;
  if (react) {
    react_defined = 1;

    // the collision kernels blind-cast react to a Kokkos react type (e.g.
    //   ReactTCEKokkos) and byte-copy it into the functor; a host-only react
    //   style would be reinterpreted as a device object -> UB.  The collide
    //   style string is "vss" (suffix-created), so react's own init() cannot
    //   catch this; require a Kokkos-enabled (ReactBirdKokkos-derived) react.

    if (!dynamic_cast<ReactBirdKokkos*>(react))
      error->all(FLERR,"Must use a Kokkos-enabled reaction style with collide vss/kk");

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
    FixAmbipolarKokkos *afix_kk = (FixAmbipolarKokkos *) afix;
    d_ions = afix_kk->d_ions;
  }

  // if ambipolar and multiple groups in mixture, ambispecies must be its own group

  if (ambiflag && mixture->ngroup > 1) {
    int *species2group = mixture->species2group;
    egroup = species2group[ambispecies];
    if (mixture->groupsize[egroup] != 1)
      error->all(FLERR,"Multigroup ambipolar collisions require "
                 "electrons be their own group");
  }

  // warn if ambipolar and a single group (e.g. collide ... all)
  // the light electrons inflate the single-group vremax, so many more
  //   collision attempts are made than with a per-species grouping
  // grouping electrons separately (e.g. collide ... species) is far faster

  if (ambiflag && mixture->ngroup == 1)
    error->warning(FLERR,"Single-group ambipolar collisions are inefficient; "
                   "grouping electrons separately (e.g. collide ... species) "
                   "is recommended");

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
      k_params.view_host()(i,j) = params[i][j];
      k_prefactor.view_host()(i,j) = prefactor[i][j];
    }
  }

  k_params.modify_host();
  k_params.sync_device();
  d_params = k_params.view_device();
  d_params_const = k_params.view_device();

  k_prefactor.modify_host();
  k_prefactor.sync_device();
  d_prefactor = k_prefactor.view_device();

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

  k_vremax.clear_sync_state();
  if (remainflag) k_remain.clear_sync_state();

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

  if (elecstyle == DISCRETE &&
      (grid->maxlocal > (int)d_cumulative_probabilities.extent(0) ||
       particle->maxelecstate > (int)d_cumulative_probabilities.extent(1)))
    MemKK::realloc_kokkos(d_cumulative_probabilities,"collide:cumulative_probabilities",grid->maxlocal,particle->maxelecstate);

  // copy Update count of gas/gas collision computes active on this timestep

  ngas_tally = update->ngas_tally;

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
  // variant for ambipolar approximation or not
  // variant for nearcp flag or not
  // variant for ngas_tally active or not
  // variant for single group or multiple groups

  // partition active gas/gas tally computes by type into typed KKCopy lists
  // each must be a Kokkos gas tally compute; call pre_gas_tally() on it
  // covers the per-grid gas/collision/grid and gas/reaction/grid, and the
  //   per-event gas/collision/tally and gas/reaction/tally

  if (ngas_tally) setup_gas_tally();

  COLLIDE_REDUCE reduce;

  if (ngroups == 1) {
    if (!ambiflag) {
      if (subcellflag) {
        if (!ngas_tally) {
          if (domain->dimension == 2) collisions_one_subcell<2,0>(reduce);
          else collisions_one_subcell<3,0>(reduce);
        } else {
          if (domain->dimension == 2) collisions_one_subcell<2,1>(reduce);
          else collisions_one_subcell<3,1>(reduce);
        }
      } else if (!nearcp) {
        if (!ngas_tally) {
          collisions_one<0,0>(reduce);
        } else if (ngas_tally) {
          collisions_one<0,1>(reduce);
        }
      } else if (nearcp) {
        if (!ngas_tally) {
          collisions_one<1,0>(reduce);
        } else if (ngas_tally) {
          collisions_one<1,1>(reduce);
        }
      }
    } else if (ambiflag) {
      if (!ngas_tally) {
        collisions_one_ambipolar<0>(reduce);
      } else if (ngas_tally) {
        collisions_one_ambipolar<1>(reduce);
      }
    }

  // multiple groups
  // both the plain and the ambipolar group paths support reactions; the
  //   plain one also supports near-neighbor selection

  } else {
    // unreachable: init() above already aborts this combination, matching the
    //   host.  Kept as a belt-and-braces assert in case the dispatch is ever
    //   reached by another path

    if (subcellflag)
      error->all(FLERR,"Cannot yet use subcell collisions with "
                 "multiple collision groups");
    if (!ambiflag) {
      if (!nearcp) {
        if (!ngas_tally) collisions_group<0,0>(reduce);
        else collisions_group<0,1>(reduce);
      } else {
        if (!ngas_tally) collisions_group<1,0>(reduce);
        else collisions_group<1,1>(reduce);
      }
    } else if (ambiflag) {
      if (!ngas_tally) {
        collisions_group_ambipolar<0>(reduce);
      } else if (ngas_tally) {
        collisions_group_ambipolar<1>(reduce);
      }
    }
  }

  // finalize active gas/gas tally computes: contribute and sync to host

  if (ngas_tally) finish_gas_tally();

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
   partition the active gas/gas tally computes (update->glist_active) into
     typed KKCopy lists and call pre_gas_tally() on each
   only the per-grid Kokkos computes are supported; the per-event
     gas/collision/tally and gas/reaction/tally computes are not
------------------------------------------------------------------------- */

void CollideVSSKokkos::setup_gas_tally()
{
  nglist_collision = nglist_reaction = 0;
  nglist_coll_tally = nglist_react_tally = 0;

  // dispatch by dynamic_cast, not by style string, so a compute the user
  //   typed with the explicit "/kk" suffix is still recognized
  // count first: the buffers have to be sized before anything is blitted in

  for (int i = 0; i < ngas_tally; i++) {
    Compute *c = update->glist_active[i];
    if (dynamic_cast<ComputeGasCollisionGridKokkos*>(c)) nglist_collision++;
    else if (dynamic_cast<ComputeGasReactionGridKokkos*>(c)) nglist_reaction++;
    else if (dynamic_cast<ComputeGasCollisionTallyKokkos*>(c)) nglist_coll_tally++;
    else if (dynamic_cast<ComputeGasReactionTallyKokkos*>(c)) nglist_react_tally++;
    else
      error->all(FLERR,"Kokkos does not (yet) support this gas tally compute; "
                       "use a Kokkos-enabled gas tally compute (-sf kk)");
  }

#ifdef SPARTA_KOKKOS_FIXED_LISTS
  if (nglist_collision > KOKKOS_MAX_GLIST || nglist_reaction > KOKKOS_MAX_GLIST ||
      nglist_coll_tally > KOKKOS_MAX_GLIST || nglist_react_tally > KOKKOS_MAX_GLIST)
    error->all(FLERR,"Kokkos supports at most KOKKOS_MAX_GLIST instances of each gas tally compute");
#else
  gas_buf_resize<ComputeGasCollisionGridKokkos>(k_glist_collision,d_glist_collision,nglist_collision);
  gas_buf_resize<ComputeGasReactionGridKokkos>(k_glist_reaction,d_glist_reaction,nglist_reaction);
  gas_buf_resize<ComputeGasCollisionTallyKokkos>(k_glist_coll_tally,d_glist_coll_tally,nglist_coll_tally);
  gas_buf_resize<ComputeGasReactionTallyKokkos>(k_glist_react_tally,d_glist_react_tally,nglist_react_tally);
#endif

  int ncg = 0, nrg = 0, nct = 0, nrt = 0;

  for (int i = 0; i < ngas_tally; i++) {
    Compute *c = update->glist_active[i];
    if (ComputeGasCollisionGridKokkos *ckk =
          dynamic_cast<ComputeGasCollisionGridKokkos*>(c)) {
      ckk->pre_gas_tally();
#ifdef SPARTA_KOKKOS_FIXED_LISTS
      glist_collision_copy[ncg++].copy(ckk);
#else
      gas_buf_blit(k_glist_collision,ncg++,ckk);
#endif
    } else if (ComputeGasReactionGridKokkos *ckk =
                 dynamic_cast<ComputeGasReactionGridKokkos*>(c)) {
      ckk->pre_gas_tally();
#ifdef SPARTA_KOKKOS_FIXED_LISTS
      glist_reaction_copy[nrg++].copy(ckk);
#else
      gas_buf_blit(k_glist_reaction,nrg++,ckk);
#endif
    } else if (ComputeGasCollisionTallyKokkos *ckk =
                 dynamic_cast<ComputeGasCollisionTallyKokkos*>(c)) {
      ckk->pre_gas_tally();
      ckk->d_overflow = d_tally_overflow;
#ifdef SPARTA_KOKKOS_FIXED_LISTS
      glist_coll_tally_copy[nct++].copy(ckk);
#else
      gas_buf_blit(k_glist_coll_tally,nct++,ckk);
#endif
    } else if (ComputeGasReactionTallyKokkos *ckk =
                 dynamic_cast<ComputeGasReactionTallyKokkos*>(c)) {
      ckk->pre_gas_tally();
      ckk->d_overflow = d_tally_overflow;
#ifdef SPARTA_KOKKOS_FIXED_LISTS
      glist_react_tally_copy[nrt++].copy(ckk);
#else
      gas_buf_blit(k_glist_react_tally,nrt++,ckk);
#endif
    }
  }

#ifdef SPARTA_KOKKOS_FIXED_LISTS
  for (int i = ncg; i < KOKKOS_MAX_GLIST; i++) glist_collision_copy[i].copy(&tmp_compute_gas_collision_kk);
  for (int i = nrg; i < KOKKOS_MAX_GLIST; i++) glist_reaction_copy[i].copy(&tmp_compute_gas_reaction_kk);
  for (int i = nct; i < KOKKOS_MAX_GLIST; i++) glist_coll_tally_copy[i].copy(&tmp_compute_gas_coll_tally_kk);
  for (int i = nrt; i < KOKKOS_MAX_GLIST; i++) glist_react_tally_copy[i].copy(&tmp_compute_gas_react_tally_kk);
#else
  gas_buf_sync(k_glist_collision,d_glist_collision);
  gas_buf_sync(k_glist_reaction,d_glist_reaction);
  gas_buf_sync(k_glist_coll_tally,d_glist_coll_tally);
  gas_buf_sync(k_glist_react_tally,d_glist_react_tally);
#endif
}

/* ----------------------------------------------------------------------
   finalize the active gas/gas tally computes
   call post_gas_tally() on the real compute objects (not the copies)
------------------------------------------------------------------------- */

void CollideVSSKokkos::finish_gas_tally()
{
  for (int i = 0; i < ngas_tally; i++) {
    Compute *c = update->glist_active[i];
    if (ComputeGasCollisionGridKokkos *ckk = dynamic_cast<ComputeGasCollisionGridKokkos*>(c))
      ckk->post_gas_tally();
    else if (ComputeGasReactionGridKokkos *ckk = dynamic_cast<ComputeGasReactionGridKokkos*>(c))
      ckk->post_gas_tally();
    else if (ComputeGasCollisionTallyKokkos *ckk = dynamic_cast<ComputeGasCollisionTallyKokkos*>(c))
      ckk->post_gas_tally();
    else if (ComputeGasReactionTallyKokkos *ckk = dynamic_cast<ComputeGasReactionTallyKokkos*>(c))
      ckk->post_gas_tally();
  }
}

/* ----------------------------------------------------------------------
   re-zero the active gas tally per-grid arrays
   called on the react/retry rollback path so tally events from the aborted
     collision pass are not double-counted when the kernel re-runs
------------------------------------------------------------------------- */

void CollideVSSKokkos::clear_gas_tally()
{
  for (int i = 0; i < ngas_tally; i++) {
    Compute *c = update->glist_active[i];
    if (ComputeGasCollisionGridKokkos *ckk = dynamic_cast<ComputeGasCollisionGridKokkos*>(c))
      ckk->clear();
    else if (ComputeGasReactionGridKokkos *ckk = dynamic_cast<ComputeGasReactionGridKokkos*>(c))
      ckk->clear();
  }
}

/* ----------------------------------------------------------------------
   NTC algorithm for a single group
------------------------------------------------------------------------- */

template < int NEARCP, int GASTALLY > void CollideVSSKokkos::collisions_one(COLLIDE_REDUCE &reduce)
{
  // loop over cells I own

  this->sync(Device,ALL_MASK);

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,PARTICLE_MASK|SPECIES_MASK);
  if (vibstyle == DISCRETE || elecstyle == DISCRETE)
    particle_kk->sync(Device,CUSTOM_MASK);
  d_particles = particle_kk->k_particles.view_device();
  d_species = particle_kk->k_species.view_device();
  d_nelecstates = particle_kk->d_nelecstates;
  d_elecstates = particle_kk->d_elecstates;
  d_elec_default_rels = particle_kk->d_elec_default_rels;
  d_elec_species_rels = particle_kk->d_elec_species_rels;
  d_enforce_spin_conservation = particle_kk->d_enforce_spin_conservation;
  d_ewhich = particle_kk->k_ewhich.view_device();
  k_eivec = particle_kk->k_eivec;
  k_eiarray = particle_kk->k_eiarray;
  k_edvec = particle_kk->k_edvec;

  GridKokkos* grid_kk = (GridKokkos*) grid;
  grid_kk->sync(Device,CINFO_MASK);
  d_plist = grid_kk->d_plist;

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

    // form the product in double and check it before it becomes an int,
    //   dellist is indexed by an int

    if (maxdelete*extra_factor > MAXSMALLINT)
      error->one(FLERR,"Per-processor delete count is too big");
    int maxdelete_extra = maxdelete*extra_factor;
    if (d_dellist.extent(0) < maxdelete_extra) {
      memoryKK->destroy_kokkos(k_dellist,dellist);
      memoryKK->create_kokkos(k_dellist,dellist,maxdelete_extra,"collide:dellist");
      d_dellist = k_dellist.view_device();
    }

    maxcellcount = particle_kk->get_maxcellcount();
    int maxcellcount_extra = maxcellcount*extra_factor;
    if (d_plist.extent(1) < maxcellcount_extra) {
      d_plist = {};
      Kokkos::resize(grid_kk->d_plist,nglocal,maxcellcount_extra);
      d_plist = grid_kk->d_plist;
      if (NEARCP)
        MemKK::realloc_kokkos(d_nn_last_partner,"collide:nn_last_partner",nglocal,maxcellcount_extra);
    }

    bigint nlocal_extra = static_cast<bigint> (particle->nlocal*extra_factor);
    if (nlocal_extra > MAXSMALLINT)
      error->one(FLERR,"Per-processor particle count is too big");
    if ((bigint) d_particles.extent(0) < nlocal_extra) {
      particle->grow(nlocal_extra - particle->nlocal);
      d_particles = particle_kk->k_particles.view_device();
      k_eivec = particle_kk->k_eivec;
      k_eiarray = particle_kk->k_eiarray;
      k_edvec = particle_kk->k_edvec;
    }
  }

  // a per-event gas tally compute can force a retry of its own, and a retry
  //   re-runs the collision pass over the same particles.  that is only sound
  //   if the particle list can be rolled back first, so the backup is not
  //   gated on react/retry when one of those computes is active

  const int tally_backup = (nglist_coll_tally || nglist_react_tally);
  const int do_backup =
    (react && sparta->kokkos->react_retry_flag) || tally_backup;

  if (tally_backup) rewind_gas_tally_computes(1);

  while (h_retry()) {

    if (do_backup) backup();

    // discard the rows an aborted attempt appended, including an attempt
    //   repeated for a reaction overflow rather than a tally overflow

    if (tally_backup) rewind_gas_tally_computes(0);

    h_retry() = 0;
    h_maxdelete() = maxdelete;
    h_maxcellcount() = maxcellcount;
    h_part_grow() = 0;
    h_ndelete() = 0;
    h_nlocal() = particle->nlocal;

    // h_tally_overflow is not zeroed anywhere else on this path: the reaction
    //   retry branch below reads maxdelete/maxcellcount/nlocal back out of
    //   h_scalars, so unlike UpdateKokkos it cannot bulk-zero the array.  A
    //   pass that raised both flags would otherwise push a stale 1 back to the
    //   device and trigger a spurious grow plus a wasted sweep next attempt

    h_tally_overflow() = 0;

    Kokkos::deep_copy(d_scalars,h_scalars);
    Kokkos::deep_copy(d_scalars_big,h_scalars_big);

    grid_kk_copy.copy(grid_kk);
    if (react) {
      ReactQKKokkos* react_qk = dynamic_cast<ReactQKKokkos*>(react);
      ReactTCEQKKokkos* react_tceqk = dynamic_cast<ReactTCEQKKokkos*>(react);
      if (react_tceqk) {
        react_style = 2;
        react_tceqk_kk_copy.copy(react_tceqk);
      } else if (react_qk) {
        react_style = 1;
        react_qk_kk_copy.copy(react_qk);
      } else {
        react_style = 0;
        react_kk_copy.copy((ReactTCEKokkos*) react);
      }
    }

    // zero the custom attributes of the slots a reaction can fill
    // must precede the kernel, not follow it: EEXCHANGE_ReactingEDisposal()
    //   sets the vibrational mode levels of the third product it just created
    // repeated on each retry, since a rolled back attempt leaves values
    //   behind in those slots

    if (react) particle_kk->zero_custom_kokkos();

    if (sparta->kokkos->atomic_reduction) {
      if (sparta->kokkos->need_atomics)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagCollideCollisionsOne<NEARCP,GASTALLY,1> >(0,nglocal),*this);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagCollideCollisionsOne<NEARCP,GASTALLY,0> >(0,nglocal),*this);
    } else
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagCollideCollisionsOne<NEARCP,GASTALLY,-1> >(0,nglocal),*this,reduce);

    Kokkos::deep_copy(h_scalars,d_scalars);
    Kokkos::deep_copy(h_scalars_big,d_scalars_big);

    // a per-event gas tally compute ran out of room: grow it and repeat the
    //   pass.  unlike a reaction overflow this needs no react/retry opt-in,
    //   and clear_gas_tally() below already discards the aborted pass

    if (h_tally_overflow() && !h_retry()) {
      grow_gas_tally_computes();
      if (do_backup) restore();
      if (ngas_tally) clear_gas_tally();
      Kokkos::deep_copy(h_scalars,0);
      Kokkos::deep_copy(h_scalars_big,0);
      reduce = COLLIDE_REDUCE();
      h_retry() = 1;
      continue;
    }

    if (h_retry()) {
      //printf("Retrying, reason %i %i %i !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n",h_maxdelete() > d_dellist.extent(0),h_maxcellcount() > d_plist.extent(1),h_part_grow());
      if (!do_backup) {
        error->one(FLERR,"Ran out of space in Kokkos collisions, increase react/extra"
                         " or use react/retry");
      } else
        restore();

      // undo gas tally events from the aborted pass before the kernel re-runs
      if (ngas_tally) clear_gas_tally();

      reduce = COLLIDE_REDUCE();

      maxdelete = h_maxdelete();
      if (d_dellist.extent(0) < maxdelete) {
        memoryKK->destroy_kokkos(k_dellist,dellist);
        memoryKK->grow_kokkos(k_dellist,dellist,maxdelete,"collide:dellist");
        d_dellist = k_dellist.view_device();
      }

      maxcellcount = h_maxcellcount();
      particle_kk->set_maxcellcount(maxcellcount);
      if (d_plist.extent(1) < maxcellcount) {
        d_plist = {};
        Kokkos::resize(grid_kk->d_plist,nglocal,maxcellcount);
        d_plist = grid_kk->d_plist;
      }

      auto nlocal_new = h_nlocal();
      if (d_particles.extent(0) < nlocal_new) {
        particle->grow(nlocal_new - particle->nlocal);
        d_particles = particle_kk->k_particles.view_device();
        k_eivec = particle_kk->k_eivec;
        k_eiarray = particle_kk->k_eiarray;
        k_edvec = particle_kk->k_edvec;
      }
    }
  }

  ndelete = h_ndelete();

  particle->nlocal = h_nlocal();

  copymode = 0;

  if (h_error_flag())
    error->one(FLERR,"Collision cell volume is zero");

  this->modified(Device,ALL_MASK);
  particle_kk->modify(Device,PARTICLE_MASK);
  if (vibstyle == DISCRETE || elecstyle == DISCRETE)
    particle_kk->modify(Device,CUSTOM_MASK);

  d_particles = t_particle_1d(); // destroy reference to reduce memory use
  d_nn_last_partner = {};
  d_plist = {};
}

KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::operator()(TagCollideZeroNN, const int &icell) const {
  const int np = grid_kk_copy.obj.d_cellcount[icell];
  for (int i = 0; i < np; i++)
    d_nn_last_partner(icell,i) = 0;
}

template < int NEARCP, int GASTALLY, int ATOMIC_REDUCTION >
KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::operator()(TagCollideCollisionsOne< NEARCP, GASTALLY, ATOMIC_REDUCTION >, const int &icell) const {
  COLLIDE_REDUCE reduce;
  this->template operator()< NEARCP, GASTALLY, ATOMIC_REDUCTION >(TagCollideCollisionsOne< NEARCP, GASTALLY, ATOMIC_REDUCTION >(), icell, reduce);
}

template < int NEARCP, int GASTALLY, int ATOMIC_REDUCTION >
KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::operator()(TagCollideCollisionsOne< NEARCP, GASTALLY, ATOMIC_REDUCTION >, const int &icell, COLLIDE_REDUCE &reduce) const {
  if (d_retry()) return;

  int np = grid_kk_copy.obj.d_cellcount[icell];
  if (np <= 1) return;

  if (NEARCP) {
    for (int i = 0; i < np; i++)
      d_nn_last_partner(icell,i) = 0;
  }

  const double volume = grid_kk_copy.obj.k_cinfo.view_device()[icell].volume / grid_kk_copy.obj.k_cinfo.view_device()[icell].weight;
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
    // if GASTALLY: tally prep with iorig/jorig, then trigger tally

    Particle::OnePart iorig,jorig;

    if (GASTALLY) {
      iorig = *ipart;
      jorig = *jpart;
    }

    int index_kpart;

    setup_collision_kokkos(ipart,jpart,precoln,postcoln);
    const int reactflag = perform_collision_kokkos(icell,ipart,jpart,kpart,precoln,postcoln,rand_gen,
                                                   recomb_part3,recomb_species,recomb_density,index_kpart);

    if (ATOMIC_REDUCTION == 1)
      Kokkos::atomic_inc(&d_ncollide_one());
    else if (ATOMIC_REDUCTION == 0)
      d_ncollide_one()++;
    else
      reduce.ncollide_one++;

    if (GASTALLY) {
      for (int m = 0; m < nglist_collision; m++)
        CVK_GLIST_COLLISION(m).template gas_tally_kk<ATOMIC_REDUCTION>(icell,reactflag,&iorig,&jorig,ipart,jpart,kpart);
      for (int m = 0; m < nglist_reaction; m++)
        CVK_GLIST_REACTION(m).template gas_tally_kk<ATOMIC_REDUCTION>(icell,reactflag,&iorig,&jorig,ipart,jpart,kpart);
      for (int m = 0; m < nglist_coll_tally; m++)
        CVK_GLIST_COLL_TALLY(m).template gas_tally_kk<ATOMIC_REDUCTION>(icell,reactflag,&iorig,&jorig,ipart,jpart,kpart);
      for (int m = 0; m < nglist_react_tally; m++)
        CVK_GLIST_REACT_TALLY(m).template gas_tally_kk<ATOMIC_REDUCTION>(icell,reactflag,&iorig,&jorig,ipart,jpart,kpart);
    }

    if (reactflag) {
      if (ATOMIC_REDUCTION == 1)
        Kokkos::atomic_inc(&d_nreact_one());
      else if (ATOMIC_REDUCTION == 0)
        d_nreact_one()++;
      else
        reduce.nreact_one++;
    } else {
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
   (re)allocate the per-cell transient subcell scratch views
   all are sized (nglocal, maxcellcount), same as d_nn_last_partner
------------------------------------------------------------------------- */

void CollideVSSKokkos::grow_subcell_views(int n1, int n2)
{
  if (int(d_nn_last_partner.extent(0)) < n1 ||
      int(d_nn_last_partner.extent(1)) < n2)
    MemKK::realloc_kokkos(d_nn_last_partner,"collide:nn_last_partner",n1,n2);

  if (int(d_subcell_id.extent(0)) < n1 || int(d_subcell_id.extent(1)) < n2) {
    MemKK::realloc_kokkos(d_subcell_id,"collide:subcell_id",n1,n2);
    MemKK::realloc_kokkos(d_subcell_count,"collide:subcell_count",n1,n2);
    MemKK::realloc_kokkos(d_subcell_first,"collide:subcell_first",n1,n2);
    MemKK::realloc_kokkos(d_subcell_next,"collide:subcell_next",n1,n2);
    MemKK::realloc_kokkos(d_subcell_ring,"collide:subcell_ring",n1,n2);
  }
}

/* ----------------------------------------------------------------------
   NTC algorithm for a single group with the transient subcell method
   Kokkos port of Collide::collisions_one_subcell()
   per-cell subcell binning is thread-private via the (icell,*) view row
------------------------------------------------------------------------- */

template < int DIM, int GASTALLY > void CollideVSSKokkos::collisions_one_subcell(COLLIDE_REDUCE &reduce)
{
  // loop over cells I own

  this->sync(Device,ALL_MASK);

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,PARTICLE_MASK|SPECIES_MASK);
  if (vibstyle == DISCRETE || elecstyle == DISCRETE)
    particle_kk->sync(Device,CUSTOM_MASK);
  d_particles = particle_kk->k_particles.view_device();
  d_species = particle_kk->k_species.view_device();
  d_nelecstates = particle_kk->d_nelecstates;
  d_elecstates = particle_kk->d_elecstates;
  d_elec_default_rels = particle_kk->d_elec_default_rels;
  d_elec_species_rels = particle_kk->d_elec_species_rels;
  d_enforce_spin_conservation = particle_kk->d_enforce_spin_conservation;
  d_ewhich = particle_kk->k_ewhich.view_device();
  k_eivec = particle_kk->k_eivec;
  k_eiarray = particle_kk->k_eiarray;
  k_edvec = particle_kk->k_edvec;

  GridKokkos* grid_kk = (GridKokkos*) grid;
  grid_kk->sync(Device,CINFO_MASK|CELL_MASK);
  d_plist = grid_kk->d_plist;

  copymode = 1;

  grow_subcell_views(nglocal,d_plist.extent(1));

  /* ATOMIC_REDUCTION: 1 = use atomics
                       0 = don't need atomics
                      -1 = use parallel_reduce
  */

  // Reactions may create or delete more particles than existing views can hold.
  //  Cannot grow a Kokkos view in a parallel loop, so
  //  if the capacity of the view is exceeded, break out of parallel loop,
  //  reallocate on the host, and then repeat the parallel loop again.

  h_retry() = 1;

  if (react) {
    double extra_factor = 1.0;
    if (sparta->kokkos->react_retry_flag)
      extra_factor = sparta->kokkos->react_extra;

    // form the product in double and check it before it becomes an int,
    //   dellist is indexed by an int

    if (maxdelete*extra_factor > MAXSMALLINT)
      error->one(FLERR,"Per-processor delete count is too big");
    int maxdelete_extra = maxdelete*extra_factor;
    if (d_dellist.extent(0) < maxdelete_extra) {
      memoryKK->destroy_kokkos(k_dellist,dellist);
      memoryKK->create_kokkos(k_dellist,dellist,maxdelete_extra,"collide:dellist");
      d_dellist = k_dellist.view_device();
    }

    maxcellcount = particle_kk->get_maxcellcount();
    int maxcellcount_extra = maxcellcount*extra_factor;
    if (d_plist.extent(1) < maxcellcount_extra) {
      d_plist = {};
      Kokkos::resize(grid_kk->d_plist,nglocal,maxcellcount_extra);
      d_plist = grid_kk->d_plist;
      grow_subcell_views(nglocal,maxcellcount_extra);
    }

    bigint nlocal_extra = static_cast<bigint> (particle->nlocal*extra_factor);
    if (nlocal_extra > MAXSMALLINT)
      error->one(FLERR,"Per-processor particle count is too big");
    if ((bigint) d_particles.extent(0) < nlocal_extra) {
      particle->grow(nlocal_extra - particle->nlocal);
      d_particles = particle_kk->k_particles.view_device();
      k_eiarray = particle_kk->k_eiarray;
    }
  }

  // a per-event gas tally compute can force a retry of its own, and a retry
  //   re-runs the collision pass over the same particles.  that is only sound
  //   if the particle list can be rolled back first, so the backup is not
  //   gated on react/retry when one of those computes is active

  const int tally_backup = (nglist_coll_tally || nglist_react_tally);
  const int do_backup =
    (react && sparta->kokkos->react_retry_flag) || tally_backup;

  if (tally_backup) rewind_gas_tally_computes(1);

  while (h_retry()) {

    if (do_backup) backup();

    // discard the rows an aborted attempt appended, including an attempt
    //   repeated for a reaction overflow rather than a tally overflow

    if (tally_backup) rewind_gas_tally_computes(0);

    h_retry() = 0;
    h_maxdelete() = maxdelete;
    h_maxcellcount() = maxcellcount;
    h_part_grow() = 0;
    h_ndelete() = 0;
    h_nlocal() = particle->nlocal;

    // h_tally_overflow is not zeroed anywhere else on this path: the reaction
    //   retry branch below reads maxdelete/maxcellcount/nlocal back out of
    //   h_scalars, so unlike UpdateKokkos it cannot bulk-zero the array.  A
    //   pass that raised both flags would otherwise push a stale 1 back to the
    //   device and trigger a spurious grow plus a wasted sweep next attempt

    h_tally_overflow() = 0;

    Kokkos::deep_copy(d_scalars,h_scalars);
    Kokkos::deep_copy(d_scalars_big,h_scalars_big);

    grid_kk_copy.copy(grid_kk);
    if (react) {
      ReactQKKokkos* react_qk = dynamic_cast<ReactQKKokkos*>(react);
      ReactTCEQKKokkos* react_tceqk = dynamic_cast<ReactTCEQKKokkos*>(react);
      if (react_tceqk) {
        react_style = 2;
        react_tceqk_kk_copy.copy(react_tceqk);
      } else if (react_qk) {
        react_style = 1;
        react_qk_kk_copy.copy(react_qk);
      } else {
        react_style = 0;
        react_kk_copy.copy((ReactTCEKokkos*) react);
      }
    }

    // zero the custom attributes of the slots a reaction can fill
    // must precede the kernel, not follow it: EEXCHANGE_ReactingEDisposal()
    //   sets the vibrational mode levels of the third product it just created
    // repeated on each retry, since a rolled back attempt leaves values
    //   behind in those slots

    if (react) particle_kk->zero_custom_kokkos();

    if (sparta->kokkos->atomic_reduction) {
      if (sparta->kokkos->need_atomics)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagCollideCollisionsOneSubcell<DIM,GASTALLY,1> >(0,nglocal),*this);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagCollideCollisionsOneSubcell<DIM,GASTALLY,0> >(0,nglocal),*this);
    } else
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagCollideCollisionsOneSubcell<DIM,GASTALLY,-1> >(0,nglocal),*this,reduce);

    Kokkos::deep_copy(h_scalars,d_scalars);
    Kokkos::deep_copy(h_scalars_big,d_scalars_big);

    // a per-event gas tally compute ran out of room: grow it and repeat the
    //   pass.  unlike a reaction overflow this needs no react/retry opt-in,
    //   and clear_gas_tally() below already discards the aborted pass

    if (h_tally_overflow() && !h_retry()) {
      grow_gas_tally_computes();
      if (do_backup) restore();
      if (ngas_tally) clear_gas_tally();
      Kokkos::deep_copy(h_scalars,0);
      Kokkos::deep_copy(h_scalars_big,0);
      reduce = COLLIDE_REDUCE();
      h_retry() = 1;
      continue;
    }

    if (h_retry()) {
      if (!do_backup) {
        error->one(FLERR,"Ran out of space in Kokkos collisions, increase react/extra"
                         " or use react/retry");
      } else
        restore();

      // undo gas tally events from the aborted pass before the kernel re-runs

      if (ngas_tally) clear_gas_tally();

      reduce = COLLIDE_REDUCE();

      maxdelete = h_maxdelete();
      if (d_dellist.extent(0) < maxdelete) {
        memoryKK->destroy_kokkos(k_dellist,dellist);
        memoryKK->grow_kokkos(k_dellist,dellist,maxdelete,"collide:dellist");
        d_dellist = k_dellist.view_device();
      }

      maxcellcount = h_maxcellcount();
      particle_kk->set_maxcellcount(maxcellcount);
      if (d_plist.extent(1) < maxcellcount) {
        d_plist = {};
        Kokkos::resize(grid_kk->d_plist,nglocal,maxcellcount);
        d_plist = grid_kk->d_plist;
        grow_subcell_views(nglocal,maxcellcount);
      }

      auto nlocal_new = h_nlocal();
      if (d_particles.extent(0) < nlocal_new) {
        particle->grow(nlocal_new - particle->nlocal);
        d_particles = particle_kk->k_particles.view_device();
        k_eiarray = particle_kk->k_eiarray;
      }
    }
  }

  ndelete = h_ndelete();

  particle->nlocal = h_nlocal();

  copymode = 0;

  if (h_error_flag())
    error->one(FLERR,"Collision cell volume is zero");

  this->modified(Device,ALL_MASK);
  particle_kk->modify(Device,PARTICLE_MASK);
  if (vibstyle == DISCRETE || elecstyle == DISCRETE)
    particle_kk->modify(Device,CUSTOM_MASK);

  d_particles = t_particle_1d(); // destroy reference to reduce memory use
  d_nn_last_partner = {};
  d_subcell_id = {};
  d_subcell_count = {};
  d_subcell_first = {};
  d_subcell_next = {};
  d_subcell_ring = {};
  d_plist = {};
}

template < int DIM, int GASTALLY, int ATOMIC_REDUCTION >
KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::operator()(TagCollideCollisionsOneSubcell< DIM, GASTALLY, ATOMIC_REDUCTION >, const int &icell) const {
  COLLIDE_REDUCE reduce;
  this->template operator()< DIM, GASTALLY, ATOMIC_REDUCTION >(TagCollideCollisionsOneSubcell< DIM, GASTALLY, ATOMIC_REDUCTION >(), icell, reduce);
}

template < int DIM, int GASTALLY, int ATOMIC_REDUCTION >
KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::operator()(TagCollideCollisionsOneSubcell< DIM, GASTALLY, ATOMIC_REDUCTION >, const int &icell, COLLIDE_REDUCE &reduce) const {
  if (d_retry()) return;

  int np = grid_kk_copy.obj.d_cellcount[icell];
  if (np <= 1) return;

  // zero previous-partner records for this cell
  // used to avoid an immediate 2nd collision of the same pair

  for (int ii = 0; ii < np; ii++)
    d_nn_last_partner(icell,ii) = 0;

  const double volume = grid_kk_copy.obj.k_cinfo.view_device()[icell].volume / grid_kk_copy.obj.k_cinfo.view_device()[icell].weight;
  if (volume == 0.0) d_error_flag() = 1;

  struct State precoln;       // state before collision
  struct State postcoln;      // state after collision

  rand_type rand_gen = rand_pool.get_state();

  // attempt = exact collision attempt count for this cell
  // nattempt = rounded attempt with RN

  const double attempt = attempt_collision_kokkos(icell,np,volume,rand_gen);
  const int nattempt = static_cast<int> (attempt);
  if (!nattempt) {
    rand_pool.free_state(rand_gen);
    return;
  }
  if (ATOMIC_REDUCTION == 1)
    Kokkos::atomic_add(&d_nattempt_one(),nattempt);
  else if (ATOMIC_REDUCTION == 0)
    d_nattempt_one() += nattempt;
  else
    reduce.nattempt_one += nattempt;

  // subcell grid: nsub subcells per dim so # subcells <= np
  //   small tolerance insures exact roots are not rounded down

  int nsub;
  if (DIM == 2) nsub = static_cast<int> (sqrt((double) np) + 1.0e-9);
  else nsub = static_cast<int> (cbrt((double) np) + 1.0e-9);
  const int nsubsq = nsub*nsub;

  auto cell = grid_kk_copy.obj.k_cells.view_device()[icell];
  double lo[3],ood[3];
  lo[0] = cell.lo[0];
  lo[1] = cell.lo[1];
  lo[2] = cell.lo[2];
  ood[0] = nsub / (cell.hi[0] - lo[0]);
  ood[1] = nsub / (cell.hi[1] - lo[1]);
  if (DIM == 3) ood[2] = nsub / (cell.hi[2] - lo[2]);
  else ood[2] = 0.0;

  rebin_subcell<DIM>(icell,np,nsub,lo,ood);

  // perform collisions
  // select random first particle, partner from same or nearby subcell
  // test if collision actually occurs

  for (int m = 0; m < nattempt; m++) {
    const int i = np * rand_gen.drand();
    const int j = find_nn_subcell<DIM>(rand_gen,i,np,icell,nsub,nsubsq);

    Particle::OnePart* ipart = &d_particles[d_plist(icell,i)];
    Particle::OnePart* jpart = &d_particles[d_plist(icell,j)];
    Particle::OnePart* kpart;

    // test if collision actually occurs, then perform it
    // continue to next collision if no reaction

    if (!test_collision_kokkos(icell,0,0,ipart,jpart,precoln,rand_gen)) continue;

    d_nn_last_partner(icell,i) = j+1;
    d_nn_last_partner(icell,j) = i+1;

    // if recombination reaction is possible for this IJ pair
    // pick a 3rd particle to participate and set cell number density
    // unless boost factor turns it off, or there is no 3rd particle

    Particle::OnePart* recomb_part3 = NULL;
    int recomb_species = -1;
    double recomb_density = 0.0;
    if (recombflag && d_recomb_ijflag(ipart->ispecies,jpart->ispecies)) {
      if (rand_gen.drand() > recomb_boost_inverse)
        recomb_species = -1;
      else if (np <= 2)
        recomb_species = -1;
      else {
        int k = np * rand_gen.drand();
        while (k == i || k == j) k = np * rand_gen.drand();
        recomb_part3 = &d_particles[d_plist(icell,k)];
        recomb_species = recomb_part3->ispecies;
        recomb_density = np * fnum / volume;
      }
    }

    // perform collision and possible reaction

    Particle::OnePart iorig,jorig;

    if (GASTALLY) {
      iorig = *ipart;
      jorig = *jpart;
    }

    int index_kpart;

    setup_collision_kokkos(ipart,jpart,precoln,postcoln);
    const int reactflag = perform_collision_kokkos(icell,ipart,jpart,kpart,precoln,postcoln,rand_gen,
                                                   recomb_part3,recomb_species,recomb_density,index_kpart);

    if (ATOMIC_REDUCTION == 1)
      Kokkos::atomic_inc(&d_ncollide_one());
    else if (ATOMIC_REDUCTION == 0)
      d_ncollide_one()++;
    else
      reduce.ncollide_one++;

    if (GASTALLY) {
      for (int m = 0; m < nglist_collision; m++)
        CVK_GLIST_COLLISION(m).template gas_tally_kk<ATOMIC_REDUCTION>(icell,reactflag,&iorig,&jorig,ipart,jpart,kpart);
      for (int m = 0; m < nglist_reaction; m++)
        CVK_GLIST_REACTION(m).template gas_tally_kk<ATOMIC_REDUCTION>(icell,reactflag,&iorig,&jorig,ipart,jpart,kpart);
      for (int m = 0; m < nglist_coll_tally; m++)
        CVK_GLIST_COLL_TALLY(m).template gas_tally_kk<ATOMIC_REDUCTION>(icell,reactflag,&iorig,&jorig,ipart,jpart,kpart);
      for (int m = 0; m < nglist_react_tally; m++)
        CVK_GLIST_REACT_TALLY(m).template gas_tally_kk<ATOMIC_REDUCTION>(icell,reactflag,&iorig,&jorig,ipart,jpart,kpart);
    }

    if (reactflag) {
      if (ATOMIC_REDUCTION == 1)
        Kokkos::atomic_inc(&d_nreact_one());
      else if (ATOMIC_REDUCTION == 0)
        d_nreact_one()++;
      else
        reduce.nreact_one++;
    } else
      continue;

    // if jpart destroyed, delete from plist, add to deletion list
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
      d_nn_last_partner(icell,j) = d_nn_last_partner(icell,np);
      unbin_one_subcell(icell,j,np);
      if (np < 2) break;
    }

    // if kpart created, add to plist
    // kpart was just added to particle list, index = index_kpart

    if (kpart) {
      if (np < d_plist.extent(1)) {
        d_nn_last_partner(icell,np) = 0;
        d_plist(icell,np++) = index_kpart;
      } else {
        d_retry() = 1;
        d_maxcellcount() += DELTACELLCOUNT;
        rand_pool.free_state(rand_gen);
        return;
      }
    }

    // if plist was changed by a reaction, rebin particles into subcells
    // so subcell vectors stay consistent with plist
    // keep same subcell grid even though np changed by one

    // a deleted particle was already unbound above by unbin_one_subcell()
    // a created particle is appended at the end of plist and moves no
    //   other particle, so it can be binned by itself in O(1)

    if (kpart) bin_one_subcell<DIM>(icell,np-1,nsub,lo,ood);
  }

  rand_pool.free_state(rand_gen);
}

/* ----------------------------------------------------------------------
   bin np particles of cell icell into transient subcell linked lists
   Kokkos port of Collide::subcell_rebin()
------------------------------------------------------------------------- */

template < int DIM >
KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::rebin_subcell(int icell, int np, int nsub,
                                     const double *lo, const double *ood) const
{
  int nsubcell = nsub*nsub;
  if (DIM == 3) nsubcell *= nsub;

  for (int isc = 0; isc < nsubcell; isc++) {
    d_subcell_count(icell,isc) = 0;
    d_subcell_first(icell,isc) = -1;
  }

  for (int n = 0; n < np; n++) {
    double *x = d_particles[d_plist(icell,n)].x;
    int ix = static_cast<int> ((x[0]-lo[0])*ood[0]);
    ix = MIN(MAX(ix,0),nsub-1);
    int iy = static_cast<int> ((x[1]-lo[1])*ood[1]);
    iy = MIN(MAX(iy,0),nsub-1);
    int iz;
    if (DIM == 3) {
      iz = static_cast<int> ((x[2]-lo[2])*ood[2]);
      iz = MIN(MAX(iz,0),nsub-1);
    } else iz = 0;

    int isc = (iz*nsub + iy)*nsub + ix;
    d_subcell_id(icell,n) = isc;
    d_subcell_next(icell,n) = d_subcell_first(icell,isc);
    d_subcell_first(icell,isc) = n;
    d_subcell_count(icell,isc)++;
  }
}

/* ----------------------------------------------------------------------
   bin the single particle at index n of cell icell's plist
   Kokkos port of Collide::subcell_bin_one()
   cell icell is owned by one thread, and every view row touched here is
     indexed by icell, so no other thread can be in these chains
------------------------------------------------------------------------- */

template < int DIM >
KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::bin_one_subcell(int icell, int n, int nsub,
                                       const double *lo, const double *ood) const
{
  const double *x = d_particles[d_plist(icell,n)].x;

  int ix = static_cast<int> ((x[0]-lo[0])*ood[0]);
  ix = MIN(MAX(ix,0),nsub-1);
  int iy = static_cast<int> ((x[1]-lo[1])*ood[1]);
  iy = MIN(MAX(iy,0),nsub-1);
  int iz;
  if (DIM == 3) {
    iz = static_cast<int> ((x[2]-lo[2])*ood[2]);
    iz = MIN(MAX(iz,0),nsub-1);
  } else iz = 0;

  int isc = (iz*nsub + iy)*nsub + ix;
  d_subcell_id(icell,n) = isc;
  d_subcell_next(icell,n) = d_subcell_first(icell,isc);
  d_subcell_first(icell,isc) = n;
  d_subcell_count(icell,isc)++;
}

/* ----------------------------------------------------------------------
   remove plist index j of cell icell from the subcell chains
   Kokkos port of Collide::subcell_unbin_one()
   caller has already done np-- and d_plist(icell,j) = d_plist(icell,np)
   chains are kept ordered by decreasing plist index, as rebin_subcell()
     leaves them, so this reproduces a full rebin exactly
   only row icell is touched, and that row belongs to this thread alone
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::unbin_one_subcell(int icell, int j, int np) const
{
  // unlink j from its own chain

  int isc = d_subcell_id(icell,j);
  int prev = -1;
  int k = d_subcell_first(icell,isc);
  while (k != j) {
    prev = k;
    k = d_subcell_next(icell,k);
  }
  if (prev < 0) d_subcell_first(icell,isc) = d_subcell_next(icell,j);
  else d_subcell_next(icell,prev) = d_subcell_next(icell,j);
  d_subcell_count(icell,isc)--;

  // if j was the last particle there is nothing to relabel

  if (np == j) return;

  // unlink old index np, relink it under its new index j in sorted order

  int jsc = d_subcell_id(icell,np);
  prev = -1;
  k = d_subcell_first(icell,jsc);
  while (k != np) {
    prev = k;
    k = d_subcell_next(icell,k);
  }
  if (prev < 0) d_subcell_first(icell,jsc) = d_subcell_next(icell,np);
  else d_subcell_next(icell,prev) = d_subcell_next(icell,np);

  prev = -1;
  k = d_subcell_first(icell,jsc);
  while (k >= 0 && k > j) {
    prev = k;
    k = d_subcell_next(icell,k);
  }
  d_subcell_next(icell,j) = k;
  if (prev < 0) d_subcell_first(icell,jsc) = j;
  else d_subcell_next(icell,prev) = j;
  d_subcell_id(icell,j) = jsc;
}

/* ----------------------------------------------------------------------
   for particle I, find collision partner J via the transient subcell method
   partner is random from same subcell, else expanding shells of subcells
   excludes an I,J pair that most recently collided with each other
   Kokkos port of the partner-selection logic in Collide::collisions_one_subcell()
------------------------------------------------------------------------- */

template < int DIM >
KOKKOS_INLINE_FUNCTION
int CollideVSSKokkos::find_nn_subcell(rand_type &rand_gen, int i, int np, int icell,
                                      int nsub, int nsubsq) const
{
  int isc = d_subcell_id(icell,i);
  int jexcl = -1;
  int j = -1;
  int jcand;

  // if another particle is in same subcell, select partner randomly from it
  // if that partner most recently collided with I, pick a different one,
  //   else fall thru to shell search for next-nearest partner

  int scount = d_subcell_count(icell,isc);
  if (scount >= 2) {
    do {
      jcand = static_cast<int> (scount*rand_gen.drand());
      j = d_subcell_first(icell,isc);
      while (jcand--) j = d_subcell_next(icell,j);
    } while (j == i);

    if (d_nn_last_partner(icell,i) == j+1 && d_nn_last_partner(icell,j) == i+1) {
      jexcl = j;
      if (scount > 2) {
        do {
          jcand = static_cast<int> (scount*rand_gen.drand());
          j = d_subcell_first(icell,isc);
          while (jcand--) j = d_subcell_next(icell,j);
        } while (j == i || j == jexcl);
      } else j = -1;
    }
  }

  // search shells of neighbor subcells with increasing radius
  //   until one or more candidate partners found
  // select partner randomly from all particles in the shell
  // shell list of subcells is clipped to bounds of subcell grid

  if (j < 0) {
    int ibox = isc % nsub;
    int jbox = (isc / nsub) % nsub;
    int kbox = isc / nsubsq;     // 0 for DIM = 2

    for (int radius = 1; radius < nsub; radius++) {
      int nring = 0;
      int ilo = MAX(ibox-radius,0);
      int ihi = MIN(ibox+radius,nsub-1);
      int jlo = MAX(jbox-radius+1,0);
      int jhi = MIN(jbox+radius-1,nsub-1);

      if (DIM == 2) {
        if (jbox-radius >= 0)
          for (int i2 = ilo; i2 <= ihi; i2++)
            d_subcell_ring(icell,nring++) = (jbox-radius)*nsub + i2;
        if (jbox+radius < nsub)
          for (int i2 = ilo; i2 <= ihi; i2++)
            d_subcell_ring(icell,nring++) = (jbox+radius)*nsub + i2;
        if (ibox-radius >= 0)
          for (int j2 = jlo; j2 <= jhi; j2++)
            d_subcell_ring(icell,nring++) = j2*nsub + (ibox-radius);
        if (ibox+radius < nsub)
          for (int j2 = jlo; j2 <= jhi; j2++)
            d_subcell_ring(icell,nring++) = j2*nsub + (ibox+radius);
      } else {
        int jflo = MAX(jbox-radius,0);
        int jfhi = MIN(jbox+radius,nsub-1);
        int klo = MAX(kbox-radius+1,0);
        int khi = MIN(kbox+radius-1,nsub-1);

        if (kbox-radius >= 0)
          for (int j2 = jflo; j2 <= jfhi; j2++)
            for (int i2 = ilo; i2 <= ihi; i2++)
              d_subcell_ring(icell,nring++) = (kbox-radius)*nsubsq + j2*nsub + i2;
        if (kbox+radius < nsub)
          for (int j2 = jflo; j2 <= jfhi; j2++)
            for (int i2 = ilo; i2 <= ihi; i2++)
              d_subcell_ring(icell,nring++) = (kbox+radius)*nsubsq + j2*nsub + i2;
        if (jbox-radius >= 0)
          for (int k2 = klo; k2 <= khi; k2++)
            for (int i2 = ilo; i2 <= ihi; i2++)
              d_subcell_ring(icell,nring++) = k2*nsubsq + (jbox-radius)*nsub + i2;
        if (jbox+radius < nsub)
          for (int k2 = klo; k2 <= khi; k2++)
            for (int i2 = ilo; i2 <= ihi; i2++)
              d_subcell_ring(icell,nring++) = k2*nsubsq + (jbox+radius)*nsub + i2;
        if (ibox-radius >= 0)
          for (int k2 = klo; k2 <= khi; k2++)
            for (int j2 = jlo; j2 <= jhi; j2++)
              d_subcell_ring(icell,nring++) = k2*nsubsq + j2*nsub + (ibox-radius);
        if (ibox+radius < nsub)
          for (int k2 = klo; k2 <= khi; k2++)
            for (int j2 = jlo; j2 <= jhi; j2++)
              d_subcell_ring(icell,nring++) = k2*nsubsq + j2*nsub + (ibox+radius);
      }

      // ncand = # of candidate partners in shell subcells
      // if none, expand search to next shell

      int ncand = 0;
      for (int mm = 0; mm < nring; mm++)
        ncand += d_subcell_count(icell,d_subcell_ring(icell,mm));
      if (!ncand) continue;

      // select random particle from all candidates in shell

      jcand = static_cast<int> (ncand*rand_gen.drand());
      int jsc = d_subcell_ring(icell,0);
      for (int mm = 0; mm < nring; mm++) {
        jsc = d_subcell_ring(icell,mm);
        if (jcand < d_subcell_count(icell,jsc)) break;
        jcand -= d_subcell_count(icell,jsc);
      }
      j = d_subcell_first(icell,jsc);
      while (jcand--) j = d_subcell_next(icell,j);

      // if partner most recently collided with I:
      // pick a different one from shell if it has others,
      //   else expand search to next shell for next-nearest partner

      if (d_nn_last_partner(icell,i) == j+1 && d_nn_last_partner(icell,j) == i+1) {
        jexcl = j;
        if (ncand > 1) {
          do {
            jcand = static_cast<int> (ncand*rand_gen.drand());
            for (int mm = 0; mm < nring; mm++) {
              jsc = d_subcell_ring(icell,mm);
              if (jcand < d_subcell_count(icell,jsc)) break;
              jcand -= d_subcell_count(icell,jsc);
            }
            j = d_subcell_first(icell,jsc);
            while (jcand--) j = d_subcell_next(icell,j);
          } while (j == jexcl);
        } else {
          j = -1;
          continue;
        }
      }
      break;
    }

    // only remaining partner is the one just collided with: accept it

    if (j < 0) j = jexcl;
  }

  return j;
}

/* ----------------------------------------------------------------------
   resize the per-group lists to match the current d_plist capacity
   a reaction can move every particle of a cell into one group, so each
     group region must be able to hold the whole cell
------------------------------------------------------------------------- */

void CollideVSSKokkos::grow_group_lists()
{
  MemKK::realloc_kokkos(d_glist,"collide:glist",nglocal,ngroups,d_plist.extent(1));
  MemKK::realloc_kokkos(d_p2g,"collide:p2g",nglocal,d_plist.extent(1),2);
  if (nearcp) {
    MemKK::realloc_kokkos(d_nn_igroup,"collide:nn_igroup",nglocal,d_plist.extent(1));
    MemKK::realloc_kokkos(d_nn_jgroup,"collide:nn_jgroup",nglocal,d_plist.extent(1));
  }
}

/* ----------------------------------------------------------------------
   NTC algorithm for multiple groups
   supports reactions and near-neighbor selection; group membership changes
     inside the kernel as reactions rebin, create and destroy particles
------------------------------------------------------------------------- */

template < int NEARCP, int GASTALLY >
void CollideVSSKokkos::collisions_group(COLLIDE_REDUCE &reduce)
{
  // loop over cells I own

  this->sync(Device,ALL_MASK);

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,PARTICLE_MASK|SPECIES_MASK);
  if (vibstyle == DISCRETE || elecstyle == DISCRETE)
    particle_kk->sync(Device,CUSTOM_MASK);
  d_particles = particle_kk->k_particles.view_device();
  d_species = particle_kk->k_species.view_device();
  d_nelecstates = particle_kk->d_nelecstates;
  d_elecstates = particle_kk->d_elecstates;
  d_elec_default_rels = particle_kk->d_elec_default_rels;
  d_elec_species_rels = particle_kk->d_elec_species_rels;
  d_enforce_spin_conservation = particle_kk->d_enforce_spin_conservation;
  d_ewhich = particle_kk->k_ewhich.view_device();
  k_eivec = particle_kk->k_eivec;
  k_eiarray = particle_kk->k_eiarray;
  k_edvec = particle_kk->k_edvec;

  GridKokkos* grid_kk = (GridKokkos*) grid;
  grid_kk->sync(Device,CINFO_MASK);
  d_plist = grid_kk->d_plist;

  // allocate per-cell group scratch arrays
  // d_glist holds plist indices laid out group-contiguous per cell
  // d_nattempt_pair holds the pre-computed attempt count per group pair

  // one region per group, each able to hold the whole cell: a reaction can
  //   move every particle of a cell into the same group

  if (int(d_glist.extent(0)) < nglocal ||
      int(d_glist.extent(1)) < ngroups ||
      int(d_glist.extent(2)) < int(d_plist.extent(1))) {
    MemKK::realloc_kokkos(d_glist,"collide:glist",nglocal,ngroups,d_plist.extent(1));
    MemKK::realloc_kokkos(d_p2g,"collide:p2g",nglocal,d_plist.extent(1),2);
  }
  if (nearcp &&
      (int(d_nn_igroup.extent(0)) < nglocal ||
       int(d_nn_igroup.extent(1)) < int(d_plist.extent(1)))) {
    MemKK::realloc_kokkos(d_nn_igroup,"collide:nn_igroup",nglocal,d_plist.extent(1));
    MemKK::realloc_kokkos(d_nn_jgroup,"collide:nn_jgroup",nglocal,d_plist.extent(1));
  }
  if (int(d_nattempt_pair.extent(0)) < nglocal ||
      int(d_nattempt_pair.extent(1)) < ngroups)
    MemKK::realloc_kokkos(d_nattempt_pair,"collide:nattempt_pair",nglocal,ngroups,ngroups);

  // d_gcount holds the per-group particle counts the kernel used to keep in a
  //   per-thread stack array with a compile-time group cap.  One row per cell,
  //   so the work item's icell is the row index: no token, no contention.
  //   Checked separately from d_glist because it does not scale with the
  //   d_plist capacity, so a reaction retry that grows d_plist leaves it alone.

  if (int(d_gcount.extent(0)) < nglocal ||
      int(d_gcount.extent(1)) < ngroups)
    MemKK::realloc_kokkos(d_gcount,"collide:gcount",nglocal,ngroups);

  copymode = 1;

  // reactions can create or delete particles, so this needs the same
  //   grow-and-repeat loop collisions_one() uses: a Kokkos view cannot be
  //   grown inside a parallel loop, so the kernel raises d_retry and returns,
  //   the host reallocates, and the pass runs again

  h_retry() = 1;

  if (react) {
    double extra_factor = 1.0;
    if (sparta->kokkos->react_retry_flag)
      extra_factor = sparta->kokkos->react_extra;

    if (maxdelete*extra_factor > MAXSMALLINT)
      error->one(FLERR,"Per-processor delete count is too big");
    int maxdelete_extra = maxdelete*extra_factor;
    if (d_dellist.extent(0) < maxdelete_extra) {
      memoryKK->destroy_kokkos(k_dellist,dellist);
      memoryKK->create_kokkos(k_dellist,dellist,maxdelete_extra,"collide:dellist");
      d_dellist = k_dellist.view_device();
    }

    maxcellcount = particle_kk->get_maxcellcount();
    int maxcellcount_extra = maxcellcount*extra_factor;
    if (d_plist.extent(1) < maxcellcount_extra) {
      d_plist = {};
      Kokkos::resize(grid_kk->d_plist,nglocal,maxcellcount_extra);
      d_plist = grid_kk->d_plist;
      grow_group_lists();
    }

    bigint nlocal_extra = static_cast<bigint> (particle->nlocal*extra_factor);
    if (nlocal_extra > MAXSMALLINT)
      error->one(FLERR,"Per-processor particle count is too big");
    if ((bigint) d_particles.extent(0) < nlocal_extra) {
      particle->grow(nlocal_extra - particle->nlocal);
      d_particles = particle_kk->k_particles.view_device();
      k_eiarray = particle_kk->k_eiarray;
    }
  }

  const int tally_backup = (nglist_coll_tally || nglist_react_tally);
  const int do_backup =
    (react && sparta->kokkos->react_retry_flag) || tally_backup;

  if (tally_backup) rewind_gas_tally_computes(1);

  while (h_retry()) {

    if (do_backup) backup();
    if (tally_backup) rewind_gas_tally_computes(0);

    h_retry() = 0;
    h_maxdelete() = maxdelete;
    h_maxcellcount() = maxcellcount;
    h_part_grow() = 0;
    h_ndelete() = 0;
    h_nlocal() = particle->nlocal;

    // h_tally_overflow is not zeroed anywhere else on this path: the reaction
    //   retry branch below reads maxdelete/maxcellcount/nlocal back out of
    //   h_scalars, so unlike UpdateKokkos it cannot bulk-zero the array.  A
    //   pass that raised both flags would otherwise push a stale 1 back to the
    //   device and trigger a spurious grow plus a wasted sweep next attempt

    h_tally_overflow() = 0;
    h_error_flag() = 0;

    Kokkos::deep_copy(d_scalars,h_scalars);
    Kokkos::deep_copy(d_scalars_big,h_scalars_big);

    grid_kk_copy.copy(grid_kk);
    if (react) {
      ReactQKKokkos* react_qk = dynamic_cast<ReactQKKokkos*>(react);
      ReactTCEQKKokkos* react_tceqk = dynamic_cast<ReactTCEQKKokkos*>(react);
      if (react_tceqk) {
        react_style = 2;
        react_tceqk_kk_copy.copy(react_tceqk);
      } else if (react_qk) {
        react_style = 1;
        react_qk_kk_copy.copy(react_qk);
      } else {
        react_style = 0;
        react_kk_copy.copy((ReactTCEKokkos*) react);
      }
    }

    if (react) particle_kk->zero_custom_kokkos();

    if (sparta->kokkos->atomic_reduction) {
      if (sparta->kokkos->need_atomics)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagCollideCollisionsGroup<NEARCP,GASTALLY,1> >(0,nglocal),*this);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagCollideCollisionsGroup<NEARCP,GASTALLY,0> >(0,nglocal),*this);
    } else
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagCollideCollisionsGroup<NEARCP,GASTALLY,-1> >(0,nglocal),*this,reduce);

    Kokkos::deep_copy(h_scalars,d_scalars);
    Kokkos::deep_copy(h_scalars_big,d_scalars_big);

    if (h_tally_overflow() && !h_retry()) {
      grow_gas_tally_computes();
      if (do_backup) restore();
      if (ngas_tally) clear_gas_tally();
      Kokkos::deep_copy(h_scalars,0);
      Kokkos::deep_copy(h_scalars_big,0);
      reduce = COLLIDE_REDUCE();
      h_retry() = 1;
      continue;
    }

    if (h_retry()) {
      if (!do_backup) {
        error->one(FLERR,"Ran out of space in Kokkos collisions, increase react/extra"
                         " or use react/retry");
      } else
        restore();

      if (ngas_tally) clear_gas_tally();

      reduce = COLLIDE_REDUCE();

      maxdelete = h_maxdelete();
      if (d_dellist.extent(0) < maxdelete) {
        memoryKK->destroy_kokkos(k_dellist,dellist);
        memoryKK->grow_kokkos(k_dellist,dellist,maxdelete,"collide:dellist");
        d_dellist = k_dellist.view_device();
      }

      maxcellcount = h_maxcellcount();
      particle_kk->set_maxcellcount(maxcellcount);
      if (d_plist.extent(1) < maxcellcount) {
        d_plist = {};
        Kokkos::resize(grid_kk->d_plist,nglocal,maxcellcount);
        d_plist = grid_kk->d_plist;
        grow_group_lists();
      }

      auto nlocal_new = h_nlocal();
      if (d_particles.extent(0) < nlocal_new) {
        particle->grow(nlocal_new - particle->nlocal);
        d_particles = particle_kk->k_particles.view_device();
        k_eiarray = particle_kk->k_eiarray;
      }
    }
  }

  ndelete = h_ndelete();

  // publish the particles the reactions created: the kernel appended them to
  //   the device list and counted them in d_nlocal, but until nlocal is
  //   carried back the host cannot see them

  particle->nlocal = h_nlocal();

  copymode = 0;

  if (h_error_flag())
    error->one(FLERR,"Collision cell volume is zero");

  this->modified(Device,ALL_MASK);
  particle_kk->modify(Device,PARTICLE_MASK);
  if (vibstyle == DISCRETE || elecstyle == DISCRETE)
    particle_kk->modify(Device,CUSTOM_MASK);

  d_particles = t_particle_1d(); // destroy reference to reduce memory use

  // d_nn_igroup/d_nn_jgroup are owned by this class, not borrowed from grid
  //   or particle, so they are not released here.  Freeing them made the
  //   guard above reallocate two nglocal x maxcellcount int arrays on every
  //   timestep of every nearcp multigroup run; they are sized by that guard
  //   and reused instead.  (d_particles and d_plist are references into
  //   other classes' allocations, so those are still dropped.)

  d_plist = {};
}

template < int NEARCP, int GASTALLY, int ATOMIC_REDUCTION >
KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::operator()(TagCollideCollisionsGroup< NEARCP, GASTALLY, ATOMIC_REDUCTION >, const int &icell) const {
  COLLIDE_REDUCE reduce;
  this->template operator()< NEARCP, GASTALLY, ATOMIC_REDUCTION >(TagCollideCollisionsGroup< NEARCP, GASTALLY, ATOMIC_REDUCTION >(), icell, reduce);
}

template < int NEARCP, int GASTALLY, int ATOMIC_REDUCTION >
KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::operator()(TagCollideCollisionsGroup< NEARCP, GASTALLY, ATOMIC_REDUCTION >, const int &icell, COLLIDE_REDUCE &reduce) const {
  // once any cell has raised d_retry the whole pass is going to be rolled
  //   back and re-run, so the remaining cells should not do a full reacting
  //   pass whose results are only going to be discarded.  the other four
  //   collision kernels all start this way

  if (d_retry()) return;

  int np = grid_kk_copy.obj.d_cellcount[icell];
  if (np <= 1) return;

  const double volume = grid_kk_copy.obj.k_cinfo.view_device()[icell].volume / grid_kk_copy.obj.k_cinfo.view_device()[icell].weight;
  if (volume == 0.0) d_error_flag() = 1;

  // build per-group particle lists for this cell
  // d_gcount(icell,g) = # of particles in group g
  // d_glist(icell,g,k) = plist index of the kth particle of group g
  // built with addgroup_kk in plist order, as the non-Kokkos version does

  for (int g = 0; g < ngroups; g++) d_gcount(icell,g) = 0;
  for (int n = 0; n < np; n++) {
    const int isp = d_particles[d_plist(icell,n)].ispecies;
    addgroup_kk(icell,d_species2group[isp],n);
  }

  struct State precoln;       // state before collision
  struct State postcoln;      // state after collision

  rand_type rand_gen = rand_pool.get_state();

  // pre-compute # of attempts for each pair of groups
  // double loop over N^2 / 2 pairs of groups
  // draw RN for every pair to match non-Kokkos collision ordering

  for (int ig = 0; ig < ngroups; ig++)
    for (int jg = ig; jg < ngroups; jg++) {
      const double attempt =
        attempt_collision_kokkos(icell,ig,jg,d_gcount(icell,ig),
                                 d_gcount(icell,jg),volume,rand_gen);
      const int nattempt = static_cast<int> (attempt);
      d_nattempt_pair(icell,ig,jg) = nattempt;
      if (nattempt) {
        if (ATOMIC_REDUCTION == 1)
          Kokkos::atomic_add(&d_nattempt_one(),nattempt);
        else if (ATOMIC_REDUCTION == 0)
          d_nattempt_one() += nattempt;
        else
          reduce.nattempt_one += nattempt;
      }
    }

  // perform collisions for each pair of groups
  // select random particle in each group, cannot be same if igroup == jgroup
  // test if collision actually occurs, then perform it

  for (int ig = 0; ig < ngroups; ig++)
    for (int jg = ig; jg < ngroups; jg++) {
      const int nattempt = d_nattempt_pair(icell,ig,jg);
      if (!nattempt) continue;
      if (d_gcount(icell,ig) == 0 || d_gcount(icell,jg) == 0) continue;
      if (ig == jg && d_gcount(icell,ig) == 1) continue;

      // near-neighbor bookkeeping is per group pair and starts cleared,
      //   as Collide::collisions_group() does via set_nn_group()

      if (NEARCP) {
        const int nclear_i = d_gcount(icell,ig);
        for (int k = 0; k < nclear_i; k++) d_nn_igroup(icell,k) = 0;
        if (ig != jg) {
          const int nclear_j = d_gcount(icell,jg);
          for (int k = 0; k < nclear_j; k++) d_nn_jgroup(icell,k) = 0;
        }
      }

      for (int iattempt = 0; iattempt < nattempt; iattempt++) {
        const int ni = d_gcount(icell,ig);
        const int nj = d_gcount(icell,jg);

        int i = ni * rand_gen.drand();
        int j;
        if (NEARCP) j = find_nn_group(rand_gen,icell,i,ig,jg,ni,nj);
        else {
          j = nj * rand_gen.drand();
          if (ig == jg)
            while (i == j) j = nj * rand_gen.drand();
        }

        const int ii = d_glist(icell,ig,i);
        const int jj = d_glist(icell,jg,j);

        Particle::OnePart* ipart = &d_particles[d_plist(icell,ii)];
        Particle::OnePart* jpart = &d_particles[d_plist(icell,jj)];

        // test if collision actually occurs

        if (!test_collision_kokkos(icell,ig,jg,ipart,jpart,precoln,rand_gen)) continue;

        if (NEARCP) {
          d_nn_igroup(icell,i) = j+1;
          if (ig == jg) d_nn_igroup(icell,j) = i+1;
          else d_nn_jgroup(icell,j) = i+1;
        }

        // if recombination is possible for this IJ pair, pick a 3rd particle
        //   and set the cell number density, unless the boost factor turns it
        //   off or there is no 3rd particle

        Particle::OnePart* recomb_part3 = NULL;
        int recomb_species = -1;
        double recomb_density = 0.0;
        if (recombflag && d_recomb_ijflag(ipart->ispecies,jpart->ispecies)) {
          if (rand_gen.drand() > recomb_boost_inverse)
            recomb_species = -1;
          else if (np <= 2)
            recomb_species = -1;
          else {
            int k = np * rand_gen.drand();
            while (k == ii || k == jj) k = np * rand_gen.drand();
            recomb_part3 = &d_particles[d_plist(icell,k)];
            recomb_species = recomb_part3->ispecies;
            recomb_density = np * fnum / volume;
          }
        }

        Particle::OnePart iorig,jorig;
        if (GASTALLY) {
          iorig = *ipart;
          jorig = *jpart;
        }

        Particle::OnePart* kpart = NULL;
        int index_kpart = 0;

        setup_collision_kokkos(ipart,jpart,precoln,postcoln);
        const int reactflag =
          perform_collision_kokkos(icell,ipart,jpart,kpart,precoln,postcoln,rand_gen,
                                   recomb_part3,recomb_species,recomb_density,index_kpart);

        if (ATOMIC_REDUCTION == 1)
          Kokkos::atomic_inc(&d_ncollide_one());
        else if (ATOMIC_REDUCTION == 0)
          d_ncollide_one()++;
        else
          reduce.ncollide_one++;

        if (GASTALLY) {
          for (int m = 0; m < nglist_collision; m++)
            CVK_GLIST_COLLISION(m).template gas_tally_kk<ATOMIC_REDUCTION>(icell,reactflag,&iorig,&jorig,ipart,jpart,kpart);
          for (int m = 0; m < nglist_reaction; m++)
            CVK_GLIST_REACTION(m).template gas_tally_kk<ATOMIC_REDUCTION>(icell,reactflag,&iorig,&jorig,ipart,jpart,kpart);
          for (int m = 0; m < nglist_coll_tally; m++)
            CVK_GLIST_COLL_TALLY(m).template gas_tally_kk<ATOMIC_REDUCTION>(icell,reactflag,&iorig,&jorig,ipart,jpart,kpart);
          for (int m = 0; m < nglist_react_tally; m++)
            CVK_GLIST_REACT_TALLY(m).template gas_tally_kk<ATOMIC_REDUCTION>(icell,reactflag,&iorig,&jorig,ipart,jpart,kpart);
        }

        if (reactflag) {
          if (ATOMIC_REDUCTION == 1)
            Kokkos::atomic_inc(&d_nreact_one());
          else if (ATOMIC_REDUCTION == 0)
            d_nreact_one()++;
          else
            reduce.nreact_one++;
        } else continue;

        // ipart may now belong to a different group

        int newgroup = d_species2group[ipart->ispecies];
        if (newgroup != ig) {
          addgroup_kk(icell,newgroup,ii);
          delgroup_kk(icell,ig,i);
          // needed if jg == ig and delgroup moved the J particle
          if (jg == ig && j == d_gcount(icell,ig)) j = i;
        }

        // jpart may now belong to a different group, or have been destroyed

        if (jpart) {
          newgroup = d_species2group[jpart->ispecies];
          if (newgroup != jg) {
            addgroup_kk(icell,newgroup,jj);
            delgroup_kk(icell,jg,j);
          }

        } else {
          const int ndelete = Kokkos::atomic_fetch_add(&d_ndelete(),1);
          if (ndelete < d_dellist.extent(0)) {
            d_dellist(ndelete) = d_plist(icell,jj);
          } else {
            d_retry() = 1;
            d_maxdelete() += DELTADELETE;
            rand_pool.free_state(rand_gen);
            return;
          }

          delgroup_kk(icell,jg,j);

          // swap-remove jj from plist and repair the moved entry's group entry
          //   through the reverse map, as Collide does with p2g

          np--;
          d_plist(icell,jj) = d_plist(icell,np);
          if (jj < np) {
            const int mg = d_p2g(icell,np,0);
            const int mk = d_p2g(icell,np,1);
            d_glist(icell,mg,mk) = jj;
            d_p2g(icell,jj,0) = mg;
            d_p2g(icell,jj,1) = mk;
          }

          if (NEARCP) {
            if (ig == jg) d_nn_igroup(icell,j) = d_nn_igroup(icell,d_gcount(icell,jg));
            else d_nn_jgroup(icell,j) = d_nn_jgroup(icell,d_gcount(icell,jg));
          }
        }

        // if kpart was created, append it to plist and to its group

        if (kpart) {
          newgroup = d_species2group[kpart->ispecies];

          if (np < d_plist.extent(1)) {
            // the host clears the new particle's slot in BOTH nn arrays of
            //   the current pair (collide.cpp:1379-1390); when ig == jg the
            //   two alias, so one write covers it

            if (NEARCP) {
              if (newgroup == ig || newgroup == jg) {
                const int n = d_gcount(icell,newgroup);
                d_nn_igroup(icell,n) = 0;
                if (ig != jg) d_nn_jgroup(icell,n) = 0;
              }
            }
            d_plist(icell,np) = index_kpart;
            addgroup_kk(icell,newgroup,np);
            np++;
          } else {
            d_retry() = 1;
            d_maxcellcount() += DELTACELLCOUNT;
            rand_pool.free_state(rand_gen);
            return;
          }
        }

        // stop attempting if either group has become too small

        const int nig = d_gcount(icell,ig);
        if (nig <= 1) {
          if (nig == 0) break;
          if (ig == jg) break;
        }
        const int njg = d_gcount(icell,jg);
        if (njg <= 1) {
          if (njg == 0) break;
          if (ig == jg) break;
        }
      }
    }

  rand_pool.free_state(rand_gen);
}

/* ----------------------------------------------------------------------
   NTC algorithm for multiple groups with ambipolar approximation
   supports reactions: group membership, the electron list and the cell
     particle list all change inside the kernel as reactions rebin, create
     and destroy particles
   ports Collide::collisions_group_ambipolar() (collide.cpp:1727-2135); the
     order in which the lists are mutated is load bearing, because rebinning
     changes which index a later random draw lands on, so every add/del is
     placed exactly where the host places it
------------------------------------------------------------------------- */

template < int GASTALLY >
void CollideVSSKokkos::collisions_group_ambipolar(COLLIDE_REDUCE &reduce)
{
  // ambipolar vectors

  this->sync(Device,ALL_MASK);

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,PARTICLE_MASK|SPECIES_MASK|CUSTOM_MASK);
  d_particles = particle_kk->k_particles.view_device();
  d_species = particle_kk->k_species.view_device();
  d_nelecstates = particle_kk->d_nelecstates;
  d_elecstates = particle_kk->d_elecstates;
  d_elec_default_rels = particle_kk->d_elec_default_rels;
  d_elec_species_rels = particle_kk->d_elec_species_rels;
  d_enforce_spin_conservation = particle_kk->d_enforce_spin_conservation;
  d_ewhich = particle_kk->k_ewhich.view_device();
  auto h_ewhich = particle_kk->k_ewhich.view_host();
  k_eivec = particle_kk->k_eivec;
  k_eiarray = particle_kk->k_eiarray;
  k_edvec = particle_kk->k_edvec;
  k_edarray = particle_kk->k_edarray;
  d_ionambi = k_eivec.view_host()[h_ewhich[index_ionambi]].k_view.view_device();
  d_velambi = k_edarray.view_host()[h_ewhich[index_velambi]].k_view.view_device();

  GridKokkos* grid_kk = (GridKokkos*) grid;
  grid_kk->sync(Device,CINFO_MASK);
  d_plist = grid_kk->d_plist;

  // allocate per-cell group scratch arrays (see collisions_group)
  // one region per group, each able to hold the whole cell: a reaction can
  //   move every particle of a cell into the same group

  if (int(d_glist.extent(0)) < nglocal ||
      int(d_glist.extent(1)) < ngroups ||
      int(d_glist.extent(2)) < int(d_plist.extent(1)))
    grow_group_lists();
  if (int(d_nattempt_pair.extent(0)) < nglocal ||
      int(d_nattempt_pair.extent(1)) < ngroups)
    MemKK::realloc_kokkos(d_nattempt_pair,"collide:nattempt_pair",nglocal,ngroups,ngroups);

  // per-cell group counters, formerly per-thread stack arrays with a
  //   compile-time group cap.  shared with collisions_group()

  if (int(d_gcount.extent(0)) < nglocal ||
      int(d_gcount.extent(1)) < ngroups)
    MemKK::realloc_kokkos(d_gcount,"collide:gcount",nglocal,ngroups);

  copymode = 1;

  // reactions can create or delete particles and electrons, so this needs the
  //   same grow-and-repeat loop the other reacting paths use: a Kokkos view
  //   cannot be grown inside a parallel loop, so the kernel raises d_retry and
  //   returns, the host reallocates, and the pass runs again

  h_retry() = 1;

  // the elist of split-off ambipolar electrons must be allocated whether or
  //   not reactions are defined: ambipolar collisions create a temporary
  //   electron for every ambipolar ion on every timestep.  only the extra
  //   sizing for reaction-created particles and deletions is react-specific

  double extra_factor = 1.0;
  if (react && sparta->kokkos->react_retry_flag)
    extra_factor = sparta->kokkos->react_extra;

  maxcellcount = particle_kk->get_maxcellcount();

  int maxelectron_extra = maxcellcount*extra_factor;
  if (int(d_elist.extent(0)) < nglocal || int(d_elist.extent(1)) < maxelectron_extra) {
    d_elist = t_particle_2d(); // reduce memory use by deallocating first
    d_elist = t_particle_2d(Kokkos::view_alloc("collide:elist",Kokkos::WithoutInitializing),nglocal,maxelectron_extra);
  }

  if (react) {
    // form the product in double and check it before it becomes an int,
    //   dellist is indexed by an int

    if (maxdelete*extra_factor > MAXSMALLINT)
      error->one(FLERR,"Per-processor delete count is too big");
    int maxdelete_extra = maxdelete*extra_factor;
    if (d_dellist.extent(0) < maxdelete_extra) {
      memoryKK->destroy_kokkos(k_dellist,dellist);
      memoryKK->grow_kokkos(k_dellist,dellist,maxdelete_extra,"collide:dellist");
      d_dellist = k_dellist.view_device();
    }

    int maxcellcount_extra = maxcellcount*extra_factor;
    if (d_plist.extent(1) < maxcellcount_extra) {
      d_plist = {};
      Kokkos::resize(grid_kk->d_plist,nglocal,maxcellcount_extra);
      d_plist = grid_kk->d_plist;
      grow_group_lists();
    }

    bigint nlocal_extra = static_cast<bigint> (particle->nlocal*extra_factor);
    if (nlocal_extra > MAXSMALLINT)
      error->one(FLERR,"Per-processor particle count is too big");
    if ((bigint) d_particles.extent(0) < nlocal_extra) {
      particle->grow(nlocal_extra - particle->nlocal);
      particle_kk->sync(Device,PARTICLE_MASK|SPECIES_MASK|CUSTOM_MASK);
      d_particles = particle_kk->k_particles.view_device();
      auto h_ewhich = particle_kk->k_ewhich.view_host();
      k_eivec = particle_kk->k_eivec;
      k_eiarray = particle_kk->k_eiarray;
      k_edarray = particle_kk->k_edarray;
      d_ionambi = k_eivec.view_host()[h_ewhich[index_ionambi]].k_view.view_device();
      d_velambi = k_edarray.view_host()[h_ewhich[index_velambi]].k_view.view_device();
    }
  }

  // a per-event gas tally compute can force a retry of its own, and a retry
  //   re-runs the collision pass over the same particles.  that is only sound
  //   if the particle list can be rolled back first, so the backup is not
  //   gated on react/retry when one of those computes is active

  const int tally_backup = (nglist_coll_tally || nglist_react_tally);
  const int do_backup =
    (react && sparta->kokkos->react_retry_flag) || tally_backup;

  if (tally_backup) rewind_gas_tally_computes(1);

  while (h_retry()) {

    if (do_backup) backup();

    // discard the rows an aborted attempt appended, including an attempt
    //   repeated for a reaction overflow rather than a tally overflow

    if (tally_backup) rewind_gas_tally_computes(0);

    h_retry() = 0;
    // seed from the allocation, not from the Collide member: maxelectron
    //   starts at 0 and the kernel only overflows against d_elist.extent(1),
    //   so seeding it with 0 makes each retry bump the request by
    //   DELTACELLCOUNT from zero and re-run the whole pass until it finally
    //   exceeds the extent -- dozens of wasted sweeps before the realloc

    maxelectron = d_elist.extent(1);
    h_maxelectron() = maxelectron;
    h_maxdelete() = maxdelete;
    h_maxcellcount() = maxcellcount;
    h_part_grow() = 0;
    h_ndelete() = 0;
    h_nlocal() = particle->nlocal;

    // h_tally_overflow is not zeroed anywhere else on this path: the reaction
    //   retry branch below reads maxdelete/maxcellcount/nlocal back out of
    //   h_scalars, so unlike UpdateKokkos it cannot bulk-zero the array.  A
    //   pass that raised both flags would otherwise push a stale 1 back to the
    //   device and trigger a spurious grow plus a wasted sweep next attempt

    h_tally_overflow() = 0;
    h_error_flag() = 0;

    Kokkos::deep_copy(d_scalars,h_scalars);
    Kokkos::deep_copy(d_scalars_big,h_scalars_big);

    grid_kk_copy.copy(grid_kk);
    if (react) {
      ReactQKKokkos* react_qk = dynamic_cast<ReactQKKokkos*>(react);
      ReactTCEQKKokkos* react_tceqk = dynamic_cast<ReactTCEQKKokkos*>(react);
      if (react_tceqk) {
        react_style = 2;
        react_tceqk_kk_copy.copy(react_tceqk);
      } else if (react_qk) {
        react_style = 1;
        react_qk_kk_copy.copy(react_qk);
      } else {
        react_style = 0;
        react_kk_copy.copy((ReactTCEKokkos*) react);
      }
    }

    // zero the custom attributes of the slots a reaction can fill
    // must precede the kernel, not follow it: ambi_reset_kokkos() sets the
    //   ion flag of the third product the reaction just created, and
    //   EEXCHANGE_ReactingEDisposal() sets its vibrational mode levels

    if (react) particle_kk->zero_custom_kokkos();

    if (sparta->kokkos->atomic_reduction) {
      if (sparta->kokkos->need_atomics)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagCollideCollisionsGroupAmbipolar<GASTALLY,1> >(0,nglocal),*this);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagCollideCollisionsGroupAmbipolar<GASTALLY,0> >(0,nglocal),*this);
    } else
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagCollideCollisionsGroupAmbipolar<GASTALLY,-1> >(0,nglocal),*this,reduce);

    Kokkos::deep_copy(h_scalars,d_scalars);
    Kokkos::deep_copy(h_scalars_big,d_scalars_big);

    // a per-event gas tally compute ran out of room: grow it and repeat the
    //   pass.  unlike a reaction overflow this needs no react/retry opt-in,
    //   and clear_gas_tally() below already discards the aborted pass

    if (h_tally_overflow() && !h_retry()) {
      grow_gas_tally_computes();
      if (do_backup) restore();
      if (ngas_tally) clear_gas_tally();
      Kokkos::deep_copy(h_scalars,0);
      Kokkos::deep_copy(h_scalars_big,0);
      reduce = COLLIDE_REDUCE();
      h_retry() = 1;
      continue;
    }

    if (h_retry()) {
      if (!do_backup) {
        error->one(FLERR,"Ran out of space in Kokkos collisions, increase react/extra"
                         " or use react/retry");
      } else
        restore();

      // undo gas tally events from the aborted pass before the kernel re-runs

      if (ngas_tally) clear_gas_tally();

      reduce = COLLIDE_REDUCE();

      maxelectron = h_maxelectron();
      if (int(d_elist.extent(1)) < maxelectron) {
        d_elist = t_particle_2d(); // reduce memory use by deallocating first
        d_elist = t_particle_2d(Kokkos::view_alloc("collide:elist",Kokkos::WithoutInitializing),nglocal,maxelectron);
      }

      maxdelete = h_maxdelete();
      if (d_dellist.extent(0) < maxdelete) {
        memoryKK->destroy_kokkos(k_dellist,dellist);
        memoryKK->grow_kokkos(k_dellist,dellist,maxdelete,"collide:dellist");
        d_dellist = k_dellist.view_device();
      }

      maxcellcount = h_maxcellcount();
      particle_kk->set_maxcellcount(maxcellcount);
      if (d_plist.extent(1) < maxcellcount) {
        d_plist = {};
        Kokkos::resize(grid_kk->d_plist,nglocal,maxcellcount);
        d_plist = grid_kk->d_plist;
        grow_group_lists();
      }

      auto nlocal_new = h_nlocal();
      if (d_particles.extent(0) < nlocal_new) {
        particle->grow(nlocal_new - particle->nlocal);
        particle_kk->sync(Device,PARTICLE_MASK|SPECIES_MASK|CUSTOM_MASK);
        d_particles = particle_kk->k_particles.view_device();
        auto h_ewhich = particle_kk->k_ewhich.view_host();
        k_eivec = particle_kk->k_eivec;
        k_eiarray = particle_kk->k_eiarray;
        k_edarray = particle_kk->k_edarray;
        d_ionambi = k_eivec.view_host()[h_ewhich[index_ionambi]].k_view.view_device();
        d_velambi = k_edarray.view_host()[h_ewhich[index_velambi]].k_view.view_device();
      }
    }
  }

  ndelete = h_ndelete();

  // publish the particles the reactions created: the kernel appended them to
  //   the device list and counted them in d_nlocal, but until nlocal is
  //   carried back the host cannot see them

  particle->nlocal = h_nlocal();

  copymode = 0;

  if (h_error_flag() == 1)
    error->one(FLERR,"Collision cell volume is zero");
  else if (h_error_flag() == 2)
    error->one(FLERR,"Collisions in cell did not conserve electron count");

  this->modified(Device,ALL_MASK);
  particle_kk->modify(Device,PARTICLE_MASK|CUSTOM_MASK);

  d_particles = t_particle_1d(); // destroy reference to reduce memory use
  d_plist = {};
}

template < int GASTALLY, int ATOMIC_REDUCTION >
KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::operator()(TagCollideCollisionsGroupAmbipolar< GASTALLY, ATOMIC_REDUCTION >, const int &icell) const {
  COLLIDE_REDUCE reduce;
  this->template operator()< GASTALLY, ATOMIC_REDUCTION >(TagCollideCollisionsGroupAmbipolar< GASTALLY, ATOMIC_REDUCTION >(), icell, reduce);
}

template < int GASTALLY, int ATOMIC_REDUCTION >
KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::operator()(TagCollideCollisionsGroupAmbipolar< GASTALLY, ATOMIC_REDUCTION >, const int &icell, COLLIDE_REDUCE &reduce) const {
  if (d_retry()) return;

  int np = grid_kk_copy.obj.d_cellcount[icell];
  if (np <= 1) return;

  const double volume = grid_kk_copy.obj.k_cinfo.view_device()[icell].volume / grid_kk_copy.obj.k_cinfo.view_device()[icell].weight;
  if (volume == 0.0) d_error_flag() = 1;

  // build the per-group particle lists for this cell and the electron list,
  //   in one pass over plist, exactly as collide.cpp:1792-1824 does
  // d_glist(icell,g,k) = plist index of the kth particle of group g
  // d_p2g(icell,n,*)   = reverse map, group and slot within it, for plist n
  // d_gcount(icell,g)  = particle count in group g
  // d_elist(icell,e)   = the eth ionized electron, split off from its ion
  //
  // the electron group egroup is the one exception: its d_gcount is
  //   maintained but its d_glist row is not written, because the host's
  //   glist[egroup][k] is always k (electrons are appended in order and
  //   removed by swapping the last one down), and both the host and this
  //   kernel index elist by the drawn index directly rather than through
  //   glist.  writing it would also need a bound the elist capacity does not
  //   share, since maxelectron can exceed the per-group row width

  for (int g = 0; g < ngroups; g++) d_gcount(icell,g) = 0;

  int nelectron = 0;
  for (int n = 0; n < np; n++) {
    const int ip = d_plist(icell,n);
    const int isp = d_particles[ip].ispecies;
    addgroup_kk(icell,d_species2group[isp],n);

    if (d_ionambi[ip]) {
      Particle::OnePart* ep = &d_elist(icell,nelectron);
      *ep = d_particles[ip];
      ep->v[0] = d_velambi(ip,0);
      ep->v[1] = d_velambi(ip,1);
      ep->v[2] = d_velambi(ip,2);
      ep->ispecies = ambispecies;
      nelectron++;
      d_gcount(icell,egroup)++;
    }
  }

  struct State precoln;       // state before collision
  struct State postcoln;      // state after collision

  rand_type rand_gen = rand_pool.get_state();

  // pre-compute # of attempts for each pair of groups
  // skip electron/electron pairs (no e/e collisions in the ambipolar model)
  // draw RN for every other pair to match non-Kokkos collision ordering

  for (int ig = 0; ig < ngroups; ig++)
    for (int jg = ig; jg < ngroups; jg++) {
      if (ig == egroup && jg == egroup) {
        d_nattempt_pair(icell,ig,jg) = 0;
        continue;
      }
      const double attempt =
        attempt_collision_kokkos(icell,ig,jg,d_gcount(icell,ig),
                                 d_gcount(icell,jg),volume,rand_gen);
      const int nattempt = static_cast<int> (attempt);
      d_nattempt_pair(icell,ig,jg) = nattempt;
      if (nattempt) {
        if (ATOMIC_REDUCTION == 1)
          Kokkos::atomic_add(&d_nattempt_one(),nattempt);
        else if (ATOMIC_REDUCTION == 0)
          d_nattempt_one() += nattempt;
        else
          reduce.nattempt_one += nattempt;
      }
    }

  // perform collisions for each pair of groups
  // electron group is always the J side, so ipart is never an electron
  //   (matches the non-Kokkos gpair igroup/jgroup flip)

  for (int ig = 0; ig < ngroups; ig++)
    for (int jg = ig; jg < ngroups; jg++) {
      if (ig == egroup && jg == egroup) continue;
      const int nattempt = d_nattempt_pair(icell,ig,jg);
      if (!nattempt) continue;

      int aig,ajg;
      if (ig == egroup) { aig = jg; ajg = ig; }
      else { aig = ig; ajg = jg; }

      // group counts are re-read from d_gcount every time, because a
      //   reaction in an earlier pair may have emptied a group

      if (d_gcount(icell,aig) == 0 || d_gcount(icell,ajg) == 0) continue;
      if (aig == ajg && d_gcount(icell,aig) == 1) continue;

      for (int iattempt = 0; iattempt < nattempt; iattempt++) {
        const int ni = d_gcount(icell,aig);
        const int nj = d_gcount(icell,ajg);

        int i = ni * rand_gen.drand();
        int j = nj * rand_gen.drand();
        if (aig == ajg)
          while (i == j) j = nj * rand_gen.drand();

        // ii/jj are plist indices, captured before any regrouping moves them
        // for the electron side there is no plist entry: elist is indexed by
        //   the drawn index itself

        const int ii = d_glist(icell,aig,i);
        const int jj = (ajg == egroup) ? -1 : d_glist(icell,ajg,j);

        Particle::OnePart* ipart = &d_particles[d_plist(icell,ii)];
        Particle::OnePart* jpart;
        if (ajg == egroup) jpart = &d_elist(icell,j);
        else jpart = &d_particles[d_plist(icell,jj)];

        // test if collision actually occurs

        if (!test_collision_kokkos(icell,aig,ajg,ipart,jpart,precoln,rand_gen)) continue;

        // if recombination reaction is possible for this IJ pair
        // pick a 3rd particle to participate and set cell number density
        // unless boost factor turns it off, or there is no 3rd particle
        // 3rd particle is never an electron since plist has no electrons
        // if ajg == egroup, no need to check k for match to jj

        Particle::OnePart* recomb_part3 = NULL;
        int recomb_species = -1;
        double recomb_density = 0.0;
        if (recombflag && d_recomb_ijflag(ipart->ispecies,jpart->ispecies)) {
          if (rand_gen.drand() > recomb_boost_inverse)
            recomb_species = -1;
          else if (np <= 2)
            recomb_species = -1;
          else {
            int k = np * rand_gen.drand();
            while (k == ii || k == jj) k = np * rand_gen.drand();
            recomb_part3 = &d_particles[d_plist(icell,k)];
            recomb_species = recomb_part3->ispecies;
            recomb_density = np * fnum / volume;
          }
        }

        // perform collision
        // if GASTALLY: save iorig/jorig, then trigger the tally

        Particle::OnePart iorig,jorig;
        if (GASTALLY) {
          iorig = *ipart;
          jorig = *jpart;
        }

        Particle::OnePart* kpart = NULL;
        int index_kpart = 0;

        const int jspecies = jpart->ispecies;
        setup_collision_kokkos(ipart,jpart,precoln,postcoln);
        const int reactflag =
          perform_collision_kokkos(icell,ipart,jpart,kpart,precoln,postcoln,rand_gen,
                                   recomb_part3,recomb_species,recomb_density,index_kpart);

        if (ATOMIC_REDUCTION == 1)
          Kokkos::atomic_inc(&d_ncollide_one());
        else if (ATOMIC_REDUCTION == 0)
          d_ncollide_one()++;
        else
          reduce.ncollide_one++;

        if (GASTALLY) {
          for (int m = 0; m < nglist_collision; m++)
            CVK_GLIST_COLLISION(m).template gas_tally_kk<ATOMIC_REDUCTION>(icell,reactflag,&iorig,&jorig,ipart,jpart,kpart);
          for (int m = 0; m < nglist_reaction; m++)
            CVK_GLIST_REACTION(m).template gas_tally_kk<ATOMIC_REDUCTION>(icell,reactflag,&iorig,&jorig,ipart,jpart,kpart);
          for (int m = 0; m < nglist_coll_tally; m++)
            CVK_GLIST_COLL_TALLY(m).template gas_tally_kk<ATOMIC_REDUCTION>(icell,reactflag,&iorig,&jorig,ipart,jpart,kpart);
          for (int m = 0; m < nglist_react_tally; m++)
            CVK_GLIST_REACT_TALLY(m).template gas_tally_kk<ATOMIC_REDUCTION>(icell,reactflag,&iorig,&jorig,ipart,jpart,kpart);
        }

        if (reactflag) {
          if (ATOMIC_REDUCTION == 1)
            Kokkos::atomic_inc(&d_nreact_one());
          else if (ATOMIC_REDUCTION == 0)
            d_nreact_one()++;
          else
            reduce.nreact_one++;
        } else continue;

        // reset ambipolar ion flags due to reaction
        // must do now before the group reset below can break out of the loop

        if (ajg == egroup)
          ambi_reset_kokkos(d_plist(icell,ii),-1,jspecies,index_kpart,
                            ipart,jpart,kpart,d_ionambi);
        else
          ambi_reset_kokkos(d_plist(icell,ii),d_plist(icell,jj),jspecies,index_kpart,
                            ipart,jpart,kpart,d_ionambi);

        // ipart may now belong to a different group
        // ipart is never an electron, so aig is never egroup here

        int newgroup = d_species2group[ipart->ispecies];
        if (newgroup != aig) {
          addgroup_kk(icell,newgroup,ii);
          delgroup_kk(icell,aig,i);
          // needed if ajg == aig and delgroup moved the J particle
          if (ajg == aig && j == d_gcount(icell,aig)) j = i;
        }

        // if kpart was created, add it to plist or elist and to its group
        // must come before the jpart code below, since that also appends

        if (kpart) {
          newgroup = d_species2group[kpart->ispecies];

          if (newgroup != egroup) {
            if (np < int(d_plist.extent(1))) {
              d_plist(icell,np) = index_kpart;
              addgroup_kk(icell,newgroup,np);
              np++;
            } else {
              d_retry() = 1;
              d_maxcellcount() += DELTACELLCOUNT;
              rand_pool.free_state(rand_gen);
              return;
            }

          } else {
            if (nelectron < int(d_elist.extent(1))) {
              Particle::OnePart* ep = &d_elist(icell,nelectron);
              *ep = *kpart;
              ep->ispecies = ambispecies;
              nelectron++;
              d_gcount(icell,egroup)++;
#ifdef SPARTA_KOKKOS_EXACT
              d_nlocal()--;
#else
              const int ndelete = Kokkos::atomic_fetch_add(&d_ndelete(),1);
              if (ndelete < int(d_dellist.extent(0))) {
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

        // jpart may now be in a different group, have become an electron,
        //   have stopped being one, or have been destroyed.  the four cases
        //   are Collide::collisions_group_ambipolar()'s, in its order

        if (jpart) {
          newgroup = d_species2group[jpart->ispecies];

          if (newgroup == ajg) {
            // nothing to do

          } else if (ajg != egroup && newgroup != egroup) {
            addgroup_kk(icell,newgroup,jj);
            delgroup_kk(icell,ajg,j);

          } else if (ajg != egroup && jpart->ispecies == ambispecies) {

            // ionization: two neutrals became an ion plus an electron.
            //   the electron goes to elist; jpart is nulled so the block
            //   below removes its now-stale plist and group entries

            if (nelectron < int(d_elist.extent(1))) {
              Particle::OnePart* ep = &d_elist(icell,nelectron);
              *ep = *jpart;
              ep->ispecies = ambispecies;
              nelectron++;
              d_gcount(icell,egroup)++;
              jpart = NULL;
            } else {
              d_retry() = 1;
              d_maxelectron() += DELTACELLCOUNT;
              rand_pool.free_state(rand_gen);
              return;
            }

          } else if (ajg == egroup && jpart->ispecies != ambispecies) {

            // exchange: an ion plus an electron became two neutrals, so the
            //   electron becomes a real particle

            const int index = Kokkos::atomic_fetch_add(&d_nlocal(),1);
            const int reallocflag =
              ParticleKokkos::add_particle_kokkos(d_particles,index,0,jspecies,icell,
                                                  jpart->x,jpart->v,0.0,0.0);
            if (reallocflag) {
              d_retry() = 1;
              d_part_grow() = 1;
              rand_pool.free_state(rand_gen);
              return;
            }

            d_particles[index] = *jpart;
            d_particles[index].id = MAXSMALLINT*rand_gen.drand();
            d_ionambi[index] = 0;

            if (nelectron-1 != j) d_elist(icell,j) = d_elist(icell,nelectron-1);
            nelectron--;
            d_gcount(icell,egroup)--;

            if (np < int(d_plist.extent(1))) {
              d_plist(icell,np) = index;
              addgroup_kk(icell,newgroup,np);
              np++;
            } else {
              d_retry() = 1;
              d_maxcellcount() += DELTACELLCOUNT;
              rand_pool.free_state(rand_gen);
              return;
            }
          }
        }

        if (!jpart && jspecies == ambispecies) {

          // recombination consumed the electron: swap the last one down,
          //   which keeps the host's glist[egroup][k] == k invariant

          if (nelectron-1 != j) d_elist(icell,j) = d_elist(icell,nelectron-1);
          nelectron--;
          d_gcount(icell,egroup)--;

        } else if (!jpart) {

          // jpart was a real particle and is gone: delete it, drop it from
          //   its group, and swap-remove it from plist, repairing the moved
          //   entry's group slot through the reverse map as Collide does

          const int ndelete = Kokkos::atomic_fetch_add(&d_ndelete(),1);
          if (ndelete < int(d_dellist.extent(0))) {
            d_dellist(ndelete) = d_plist(icell,jj);
          } else {
            d_retry() = 1;
            d_maxdelete() += DELTADELETE;
            rand_pool.free_state(rand_gen);
            return;
          }

          delgroup_kk(icell,ajg,j);

          np--;
          d_plist(icell,jj) = d_plist(icell,np);
          if (jj < np) {
            const int mg = d_p2g(icell,np,0);
            const int mk = d_p2g(icell,np,1);
            d_glist(icell,mg,mk) = jj;
            d_p2g(icell,jj,0) = mg;
            d_p2g(icell,jj,1) = mk;
          }
        }

        // stop attempting if either group has become too small

        const int nig = d_gcount(icell,aig);
        if (nig <= 1) {
          if (nig == 0) break;
          if (aig == ajg) break;
        }
        const int njg = d_gcount(icell,ajg);
        if (njg <= 1) {
          if (njg == 0) break;
          if (aig == ajg) break;
        }
      }
    }

  // recombine ambipolar ions with their matching electrons
  //   by copying the (possibly scattered) electron velocity back into velambi
  // which ion is paired with which electron does not matter

  int melectron = 0;
  for (int n = 0; n < np; n++) {
    const int i = d_plist(icell,n);
    if (d_ionambi[i]) {
      if (melectron < nelectron) {
        Particle::OnePart* ep = &d_elist(icell,melectron);
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

/* ----------------------------------------------------------------------
   NTC algorithm for a single group with ambipolar approximation
------------------------------------------------------------------------- */

template < int GASTALLY >
void CollideVSSKokkos::collisions_one_ambipolar(COLLIDE_REDUCE &reduce)
{
  // ambipolar vectors

  this->sync(Device,ALL_MASK);

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,PARTICLE_MASK|SPECIES_MASK|CUSTOM_MASK);
  d_particles = particle_kk->k_particles.view_device();
  d_species = particle_kk->k_species.view_device();
  d_nelecstates = particle_kk->d_nelecstates;
  d_elecstates = particle_kk->d_elecstates;
  d_elec_default_rels = particle_kk->d_elec_default_rels;
  d_elec_species_rels = particle_kk->d_elec_species_rels;
  d_enforce_spin_conservation = particle_kk->d_enforce_spin_conservation;
  d_ewhich = particle_kk->k_ewhich.view_device();
  auto h_ewhich = particle_kk->k_ewhich.view_host();
  k_eivec = particle_kk->k_eivec;
  k_eiarray = particle_kk->k_eiarray;
  k_edvec = particle_kk->k_edvec;
  k_edarray = particle_kk->k_edarray;
  d_ionambi = k_eivec.view_host()[h_ewhich[index_ionambi]].k_view.view_device();
  d_velambi = k_edarray.view_host()[h_ewhich[index_velambi]].k_view.view_device();

  GridKokkos* grid_kk = (GridKokkos*) grid;
  grid_kk->sync(Device,CINFO_MASK);
  d_plist = grid_kk->d_plist;

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

  // the elist of split-off ambipolar electrons must be allocated whether
  // or not reactions are defined: ambipolar collisions create a temporary
  // electron for every ambipolar ion on every timestep.  Only the extra
  // sizing for reaction-created particles/deletions is react-specific.

  double extra_factor = 1.0;
  if (react && sparta->kokkos->react_retry_flag)
    extra_factor = sparta->kokkos->react_extra;

  maxcellcount = particle_kk->get_maxcellcount();

  int maxelectron_extra = maxcellcount*extra_factor;
  if (d_elist.extent(0) < nglocal || d_elist.extent(1) < maxelectron_extra) {
    d_elist = t_particle_2d(); // reduce memory use by deallocating first
    d_elist = t_particle_2d(Kokkos::view_alloc("collide:elist",Kokkos::WithoutInitializing),nglocal,maxelectron_extra);
  }

  if (react) {
    // form the product in double and check it before it becomes an int,
    //   dellist is indexed by an int

    if (maxdelete*extra_factor > MAXSMALLINT)
      error->one(FLERR,"Per-processor delete count is too big");
    int maxdelete_extra = maxdelete*extra_factor;
    if (d_dellist.extent(0) < maxdelete_extra) {
      memoryKK->destroy_kokkos(k_dellist,dellist);
      memoryKK->grow_kokkos(k_dellist,dellist,maxdelete_extra,"collide:dellist");
      d_dellist = k_dellist.view_device();
    }

    int maxcellcount_extra = maxcellcount*extra_factor;
    if (d_plist.extent(1) < maxcellcount_extra) {
      d_plist = {};
      Kokkos::resize(grid_kk->d_plist,nglocal,maxcellcount_extra);
      d_plist = grid_kk->d_plist;
    }

    bigint nlocal_extra = static_cast<bigint> (particle->nlocal*extra_factor);
    if (nlocal_extra > MAXSMALLINT)
      error->one(FLERR,"Per-processor particle count is too big");
    if ((bigint) d_particles.extent(0) < nlocal_extra) {
      particle->grow(nlocal_extra - particle->nlocal);
      particle_kk->sync(Device,PARTICLE_MASK|SPECIES_MASK|CUSTOM_MASK);
      d_particles = particle_kk->k_particles.view_device();
      auto h_ewhich = particle_kk->k_ewhich.view_host();
      k_eivec = particle_kk->k_eivec;
      k_eiarray = particle_kk->k_eiarray;
      k_edvec = particle_kk->k_edvec;
      k_edarray = particle_kk->k_edarray;
      d_ionambi = k_eivec.view_host()[h_ewhich[index_ionambi]].k_view.view_device();
      d_velambi = k_edarray.view_host()[h_ewhich[index_velambi]].k_view.view_device();
    }
  }

  // a per-event gas tally compute can force a retry of its own, and a retry
  //   re-runs the collision pass over the same particles.  that is only sound
  //   if the particle list can be rolled back first, so the backup is not
  //   gated on react/retry when one of those computes is active

  const int tally_backup = (nglist_coll_tally || nglist_react_tally);
  const int do_backup =
    (react && sparta->kokkos->react_retry_flag) || tally_backup;

  if (tally_backup) rewind_gas_tally_computes(1);

  while (h_retry()) {

    if (do_backup) backup();

    // discard the rows an aborted attempt appended, including an attempt
    //   repeated for a reaction overflow rather than a tally overflow

    if (tally_backup) rewind_gas_tally_computes(0);

    h_retry() = 0;
    // seed from the allocation, not from the Collide member: maxelectron
    //   starts at 0 and the kernel only overflows against d_elist.extent(1),
    //   so seeding it with 0 makes each retry bump the request by
    //   DELTACELLCOUNT from zero and re-run the whole pass until it finally
    //   exceeds the extent -- dozens of wasted sweeps before the realloc

    maxelectron = d_elist.extent(1);
    h_maxelectron() = maxelectron;
    h_maxdelete() = maxdelete;
    h_maxcellcount() = maxcellcount;
    h_part_grow() = 0;
    h_ndelete() = 0;
    h_nlocal() = particle->nlocal;

    // h_tally_overflow is not zeroed anywhere else on this path: the reaction
    //   retry branch below reads maxdelete/maxcellcount/nlocal back out of
    //   h_scalars, so unlike UpdateKokkos it cannot bulk-zero the array.  A
    //   pass that raised both flags would otherwise push a stale 1 back to the
    //   device and trigger a spurious grow plus a wasted sweep next attempt

    h_tally_overflow() = 0;

    Kokkos::deep_copy(d_scalars,h_scalars);
    Kokkos::deep_copy(d_scalars_big,h_scalars_big);

    grid_kk_copy.copy(grid_kk);
    if (react) {
      ReactQKKokkos* react_qk = dynamic_cast<ReactQKKokkos*>(react);
      ReactTCEQKKokkos* react_tceqk = dynamic_cast<ReactTCEQKKokkos*>(react);
      if (react_tceqk) {
        react_style = 2;
        react_tceqk_kk_copy.copy(react_tceqk);
      } else if (react_qk) {
        react_style = 1;
        react_qk_kk_copy.copy(react_qk);
      } else {
        react_style = 0;
        react_kk_copy.copy((ReactTCEKokkos*) react);
      }
    }

    // zero the custom attributes of the slots a reaction can fill
    // must precede the kernel, not follow it: ambi_reset_kokkos() sets the
    //   ion flag of the third product the reaction just created, and
    //   EEXCHANGE_ReactingEDisposal() sets its vibrational mode levels
    // repeated on each retry, since a rolled back attempt leaves values
    //   behind in those slots

    if (react) particle_kk->zero_custom_kokkos();

    if (sparta->kokkos->atomic_reduction) {
      if (sparta->kokkos->need_atomics)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagCollideCollisionsOneAmbipolar<GASTALLY,1> >(0,nglocal),*this);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagCollideCollisionsOneAmbipolar<GASTALLY,0> >(0,nglocal),*this);
    } else
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagCollideCollisionsOneAmbipolar<GASTALLY,-1> >(0,nglocal),*this,reduce);

    Kokkos::deep_copy(h_scalars,d_scalars);
    Kokkos::deep_copy(h_scalars_big,d_scalars_big);

    // a per-event gas tally compute ran out of room: grow it and repeat the
    //   pass.  unlike a reaction overflow this needs no react/retry opt-in,
    //   and clear_gas_tally() below already discards the aborted pass

    if (h_tally_overflow() && !h_retry()) {
      grow_gas_tally_computes();
      if (do_backup) restore();
      if (ngas_tally) clear_gas_tally();
      Kokkos::deep_copy(h_scalars,0);
      Kokkos::deep_copy(h_scalars_big,0);
      reduce = COLLIDE_REDUCE();
      h_retry() = 1;
      continue;
    }

    if (h_retry()) {
      //printf("Retrying, reason %i %i %i %i !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n",h_maxelectron() > d_elist.extent(1),h_maxdelete() > d_dellist.extent(0),h_maxcellcount() > d_plist.extent(1),h_part_grow());
      //printf("%i %i %i %i %i %i %i !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n",h_maxelectron(),d_elist.extent(1),h_maxdelete(),d_dellist.extent(0),h_maxcellcount(),d_plist.extent(1),h_part_grow());

      if (!do_backup) {
        error->one(FLERR,"Ran out of space in Kokkos collisions, increase react/extra"
                         " or use react/retry");
      } else
        restore();

      // undo gas tally events from the aborted pass before the kernel re-runs
      if (ngas_tally) clear_gas_tally();

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
        d_dellist = k_dellist.view_device();
      }

      maxcellcount = h_maxcellcount();
      particle_kk->set_maxcellcount(maxcellcount);
      if (d_plist.extent(1) < maxcellcount) {
        d_plist = {};
        Kokkos::resize(grid_kk->d_plist,nglocal,maxcellcount);
        d_plist = grid_kk->d_plist;
      }

      auto nlocal_new = h_nlocal();
      if (d_particles.extent(0) < nlocal_new) {
        particle->grow(nlocal_new - particle->nlocal);
        particle_kk->sync(Device,PARTICLE_MASK|SPECIES_MASK|CUSTOM_MASK);
        d_particles = particle_kk->k_particles.view_device();
        auto h_ewhich = particle_kk->k_ewhich.view_host();
        k_eivec = particle_kk->k_eivec;
        k_eiarray = particle_kk->k_eiarray;
        k_edvec = particle_kk->k_edvec;
        k_edarray = particle_kk->k_edarray;
        d_ionambi = k_eivec.view_host()[h_ewhich[index_ionambi]].k_view.view_device();
        d_velambi = k_edarray.view_host()[h_ewhich[index_velambi]].k_view.view_device();
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

  this->modified(Device,ALL_MASK);
  particle_kk->modify(Device,PARTICLE_MASK|CUSTOM_MASK);

  d_particles = t_particle_1d(); // destroy reference to reduce memory use
  d_plist = {};
}

template < int GASTALLY, int ATOMIC_REDUCTION >
KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::operator()(TagCollideCollisionsOneAmbipolar< GASTALLY, ATOMIC_REDUCTION >, const int &icell) const {
  COLLIDE_REDUCE reduce;
  this->template operator()< GASTALLY, ATOMIC_REDUCTION >(TagCollideCollisionsOneAmbipolar< GASTALLY, ATOMIC_REDUCTION >(), icell, reduce);
}

template < int GASTALLY, int ATOMIC_REDUCTION >
KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::operator()(TagCollideCollisionsOneAmbipolar< GASTALLY, ATOMIC_REDUCTION >, const int &icell, COLLIDE_REDUCE &reduce) const {
  if (d_retry()) return;

  int np = grid_kk_copy.obj.d_cellcount[icell];
  if (np <= 1) return;

  const double volume = grid_kk_copy.obj.k_cinfo.view_device()[icell].volume / grid_kk_copy.obj.k_cinfo.view_device()[icell].weight;
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
    // no collision is performed, so it must not be counted as one; see the
    //   same test in Collide::collisions_one_ambipolar()

    if (ipart->ispecies == ambispecies && jpart->ispecies == ambispecies)
      continue;

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
    // if GASTALLY: tally prep with iorig/jorig, then trigger tally

    Particle::OnePart iorig,jorig;

    if (GASTALLY) {
      iorig = *ipart;
      jorig = *jpart;
    }

    int index_kpart = 0;

    const int jspecies = jpart->ispecies;
    setup_collision_kokkos(ipart,jpart,precoln,postcoln);
    const int reactflag = perform_collision_kokkos(icell,ipart,jpart,kpart,precoln,postcoln,rand_gen,
                                                   recomb_part3,recomb_species,recomb_density,index_kpart);

    if (ATOMIC_REDUCTION == 1)
      Kokkos::atomic_fetch_add(&d_ncollide_one(),1);
    else if (ATOMIC_REDUCTION == 0)
      d_ncollide_one()++;
    else
      reduce.ncollide_one++;

    if (GASTALLY) {
      for (int m = 0; m < nglist_collision; m++)
        CVK_GLIST_COLLISION(m).template gas_tally_kk<ATOMIC_REDUCTION>(icell,reactflag,&iorig,&jorig,ipart,jpart,kpart);
      for (int m = 0; m < nglist_reaction; m++)
        CVK_GLIST_REACTION(m).template gas_tally_kk<ATOMIC_REDUCTION>(icell,reactflag,&iorig,&jorig,ipart,jpart,kpart);
      for (int m = 0; m < nglist_coll_tally; m++)
        CVK_GLIST_COLL_TALLY(m).template gas_tally_kk<ATOMIC_REDUCTION>(icell,reactflag,&iorig,&jorig,ipart,jpart,kpart);
      for (int m = 0; m < nglist_react_tally; m++)
        CVK_GLIST_REACT_TALLY(m).template gas_tally_kk<ATOMIC_REDUCTION>(icell,reactflag,&iorig,&jorig,ipart,jpart,kpart);
    }

    if (reactflag) {
      if (ATOMIC_REDUCTION == 1)
        Kokkos::atomic_fetch_add(&d_nreact_one(),1);
      else if (ATOMIC_REDUCTION == 0)
        d_nreact_one()++;
      else
        reduce.nreact_one++;
    } else {
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

 // MCF scheme: attempt count is a Poisson variate whose mean is the
 //   majorant collision frequency x timestep, remain is not used

 if (mcflag)
   return poisson_kokkos(0.5 * np * (np-1) *
                         d_vremax(icell,0,0) * dt * fnum / volume, rand_gen);

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
   attempt count for a pair of groups
   ni,nj = particle counts in igroup,jgroup
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
double CollideVSSKokkos::attempt_collision_kokkos(int icell, int igroup, int jgroup,
                                                  int ni, int nj, double volume,
                                                  rand_type &rand_gen) const
{
  double nattempt;

  // return 2x the value for igroup != jgroup, since no J,I pairing

  // compute npairs in double, else the igroup != jgroup int*int product
  //   can overflow a 32-bit int for large per-cell group counts

  double npairs;
  if (igroup == jgroup) npairs = 0.5 * ni * (ni-1);
  else npairs = (double) ni * nj;

  nattempt = npairs * d_vremax(icell,igroup,jgroup) * dt * fnum / volume;

  // MCF scheme: attempt count is a Poisson variate whose mean is the
  //   majorant collision frequency x timestep, remain is not used

  if (mcflag) return poisson_kokkos(nattempt,rand_gen);

  if (remainflag) {
    nattempt += d_remain(icell,igroup,jgroup);
    d_remain(icell,igroup,jgroup) = nattempt - static_cast<int> (nattempt);
  } else nattempt += rand_gen.drand();

  return nattempt;
}

/* ----------------------------------------------------------------------
   Poisson RN with specified mean, on device
   returned as a double with an exact integer value
   Knuth multiplication method for small mean,
   else normal approximation with continuity correction
   mirrors RanKnuth::poisson() used by the non-Kokkos path
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
double CollideVSSKokkos::poisson_kokkos(double mean, rand_type &rand_gen) const
{
  if (mean <= 0.0) return 0.0;

  if (mean < 30.0) {
    double L = exp(-mean);
    double p = 1.0;
    int k = 0;
    do {
      k++;
      p *= rand_gen.drand();
    } while (p > L);
    return (double) (k-1);
  }

  double value = floor(mean + sqrt(mean)*rand_gen.normal() + 0.5);
  if (value < 0.0) return 0.0;
  return value;
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

  // prevent division by zero

  if (vr2 < EPSZERO && d_params(ispecies,jspecies).omega >= 1.0)
    return 0;

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
  precoln.ave_dof = (precoln.ave_rotdof + precoln.ave_vibdof)/2.;

  precoln.imass = d_species[isp].mass;
  precoln.jmass = d_species[jsp].mass;

  precoln.etrans = 0.5 * d_params(isp,jsp).mr * precoln.vr2;
  precoln.erot = ip->erot + jp->erot;
  precoln.evib = ip->evib + jp->evib;
  precoln.eelec = 0.0;
  if (elecstyle == DISCRETE) {
    auto &d_eelecs = k_edvec.view_device()[d_ewhich[index_eelec]].k_view.view_device();
    if (d_nelecstates[isp] > 0)
      precoln.eelec += d_eelecs[ip - d_particles.data()];
    if (d_nelecstates[jsp] > 0)
      precoln.eelec += d_eelecs[jp - d_particles.data()];
  }

  precoln.eint   = precoln.erot + precoln.evib + precoln.eelec;
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
  postcoln.eelec = 0.0;
  postcoln.eint = 0.0;
  postcoln.etotal = precoln.etotal;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
int CollideVSSKokkos::perform_collision_kokkos(int icell,
                                  Particle::OnePart *&ip,
                                  Particle::OnePart *&jp,
                                  Particle::OnePart *&kp,
                                  struct State &precoln, struct State &postcoln, rand_type &rand_gen,
                                  Particle::OnePart *&p3, int &recomb_species, double &recomb_density,
                                  int &index_kpart) const
{
  int reaction,kspecies;
  double x[3],v[3];

  // if gas-phase chemistry defined, attempt and perform reaction
  // if a 3rd particle is created, its kspecies >= 0 is returned
  // if 2nd particle is removed, its jspecies is set to -1
  // reaction = 0 if no reaction occurs
  // reaction = 1 to N for which reaction occurs
  // reaction is returned to caller

  // effective electronic DoF for the TCE reaction model
  // read from per-state input (dof), 0.0 if species has no electronic states

  double idof = 0.0, jdof = 0.0;
  if (react_defined && elecstyle == DISCRETE) {
    auto& d_estate = k_eivec.view_device()[d_ewhich[index_elecstate]].k_view.view_device();
    if (d_nelecstates[ip->ispecies] > 0)
      idof = d_elecstates(ip->ispecies,d_estate[ip - d_particles.data()]).dof;
    if (d_nelecstates[jp->ispecies] > 0)
      jdof = d_elecstates(jp->ispecies,d_estate[jp - d_particles.data()]).dof;
  }

  if (react_defined) {
    if (react_style == 1)
      reaction = react_qk_kk_copy.obj.attempt_kk(ip,jp,
                                             precoln.etrans,precoln.erot,
                                             precoln.evib,precoln.eelec,postcoln.etotal,kspecies,
                                             recomb_species,recomb_density,d_species);
    else if (react_style == 2)
      reaction = react_tceqk_kk_copy.obj.attempt_kk(ip,jp,
                                             precoln.etrans,precoln.erot,
                                             precoln.evib,precoln.eelec,postcoln.etotal,kspecies,
                                             recomb_species,recomb_density,d_species);
    else
      reaction = react_kk_copy.obj.attempt_kk(ip,jp,
                                             precoln.etrans,precoln.erot,
                                             precoln.evib,precoln.eelec,postcoln.etotal,kspecies,
                                             recomb_species,recomb_density,d_species,idof,jdof);
  } else reaction = 0;

  // just collision, no reaction

  if (!reaction) {
    // ave_dof counts only rot/vib DOF, so also call the energy disposal
    // when either species has electronic states (e.g. two atoms),
    // else their electronic modes would never relax
    if (precoln.ave_dof > 0.0 || elec_exchange(ip,jp))
      EEXCHANGE_NonReactingEDisposal(icell,ip,jp,precoln,postcoln,rand_gen);
    SCATTER_TwoBodyScattering(ip,jp,precoln,postcoln,rand_gen);
    return reaction;
  }

  // reaction took place
  // repartition energy and perform velocity scattering for I,J,K particles
  // reaction may have changed species of I,J particles
  // J,K particles may have been removed or created by reaction

  kp = NULL;

  // add 3rd K particle if reaction created it
  // index of new K particle = nlocal-1
  // if add_particle() performs a realloc:
  //   make copy of x,v

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
    EEXCHANGE_ReactingEDisposal(icell,ip,jp,kp,precoln,postcoln,rand_gen);
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
    if (elecstyle == DISCRETE) {
      auto &d_eelecs = k_edvec.view_device()[d_ewhich[index_eelec]].k_view.view_device();
      if (d_nelecstates[p3->ispecies] > 0)
        partial_energy += d_eelecs[p3 - d_particles.data()];
    }

    ip->erot = 0.0;
    ip->evib = 0.0;
    p3->erot = 0.0;
    p3->evib = 0.0;

    // zero electronic energy like erot/evib above: it is already part of
    // partial_energy, so setup_collision() must not count it a second time

    if (elecstyle == DISCRETE) {
      zero_elec(ip);
      zero_elec(p3);
    }

    // 2nd call to setup_collision() sets new postcoln.etotal
    // then add saved partial_energy to it

    setup_collision_kokkos(ip,p3,precoln,postcoln);
    postcoln.etotal += partial_energy;

    if (precoln.ave_dof > 0.0 || elec_exchange(ip,p3))
      EEXCHANGE_ReactingEDisposal(icell,ip,p3,jp,precoln,postcoln,rand_gen);
    SCATTER_TwoBodyScattering(ip,p3,precoln,postcoln,rand_gen);

  } else {
    EEXCHANGE_ReactingEDisposal(icell,ip,jp,kp,precoln,postcoln,rand_gen);
    SCATTER_TwoBodyScattering(ip,jp,precoln,postcoln,rand_gen);
  }

  return reaction;
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
void CollideVSSKokkos::EEXCHANGE_NonReactingEDisposal(int icell,
                                                      Particle::OnePart *ip,
                                                      Particle::OnePart *jp,
                                                      struct State &precoln, struct State &postcoln,
                                                      rand_type &rand_gen) const
{
  double State_prob,Fraction_Rot,Fraction_Vib,E_Dispose;
  int i,rotdof,vibdof,max_level,ivib;

  Particle::OnePart *p, *p2;

  double AdjustFactor = 0.99999999;
  postcoln.erot = 0.0;
  postcoln.evib = 0.0;
  postcoln.eelec = 0.0;
  double pevib = 0.0;

  // handle each kind of energy disposal for non-reacting reactants
  // enter the disposal loop even if ave_dof (rot/vib) is zero when either
  // species has electronic states, so those can still relax;
  // the rot/vib blocks below are skipped naturally via rotdof/vibdof = 0

  if (precoln.ave_dof == 0 && !elec_exchange(ip,jp)) {
    ip->erot = 0.0;
    jp->erot = 0.0;
    ip->evib = 0.0;
    jp->evib = 0.0;

  } else {
    E_Dispose = precoln.etrans;

    // This is pairwise Borgnakke-Larsen relaxation: each internal mode that
    // relaxes adds back only its OWN energy (E_Dispose += p->erot; sample;
    // E_Dispose -= p->erot), so it exchanges energy with the translational
    // pool alone.  No shared multi-mode pool is formed, each exchange
    // independently satisfies detailed balance, and the exponent is the plain
    // translational one.  Do NOT add the reacting path's remaining_dof
    // (Dirichlet stick-breaking) correction here -- that correction exists
    // only because EEXCHANGE_ReactingEDisposal splits the full collision
    // energy among all modes at once from a single depleting pool.

    for (i = 0; i < 2; i++) {
      if (i == 0) {
        p = ip;
        p2 = jp;
      } else {
        p = jp;
        p2 = ip;
      }

      int sp = p->ispecies;
      rotdof = d_species[sp].rotdof;
      double rotn_phi = d_species[sp].rotrel;

      if (rotdof) {
        if (relaxflag == VARIABLE) rotn_phi = rotrel(sp,E_Dispose+p->erot);
        if (rotn_phi >= rand_gen.drand()) {
          if (rotstyle == NONE) {
            p->erot = 0.0;
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
              const auto &d_vibmode = k_eiarray.view_device()[d_ewhich[index_vibmode]].k_view.view_device();
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

      if (elecstyle == DISCRETE && d_nelecstates[sp] > 0) {
        auto &d_estates = k_eivec.view_device()[d_ewhich[index_elecstate]].k_view.view_device();
        double elec_phi = get_elec_phi(p->ispecies, p2->ispecies, d_estates[p - d_particles.data()], E_Dispose);
        if (elec_phi >= rand_gen.drand()) {
          relax_electronic_mode(icell, p, p2, E_Dispose,
                                d_params(p->ispecies,p2->ispecies).omega,
                                rand_gen, false);
        }
      }
    }
  }

  // compute post-collision internal energies

  postcoln.erot = ip->erot + jp->erot;
  postcoln.evib = ip->evib + jp->evib;
  postcoln.eelec = 0.0;
  if (elecstyle == DISCRETE) {
    auto &d_eelecs = k_edvec.view_device()[d_ewhich[index_eelec]].k_view.view_device();
    if (d_nelecstates[ip->ispecies] > 0)
      postcoln.eelec += d_eelecs[ip - d_particles.data()];
    if (d_nelecstates[jp->ispecies] > 0)
      postcoln.eelec += d_eelecs[jp - d_particles.data()];
  }

  // compute portion of energy left over for scattering

  postcoln.eint = postcoln.erot + postcoln.evib + postcoln.eelec;
  postcoln.etrans = E_Dispose;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::relax_electronic_mode(int icell,
                                             Particle::OnePart *p,
                                             Particle::OnePart *jp,
                                             double& E_Dispose,
                                             double omega,
                                             rand_type &rand_gen,
                                             bool reacting) const
{
  auto &d_eelecs = k_edvec.view_device()[d_ewhich[index_eelec]].k_view.view_device();
  E_Dispose += d_eelecs[p - d_particles.data()];

  int ielec = select_elec_state(
    icell, p, jp, E_Dispose, omega,
    d_enforce_spin_conservation(p->ispecies,jp->ispecies) && !reacting,
    rand_gen, reacting);
  double eelec = d_elecstates(p->ispecies,ielec).temp*boltz;

  d_eelecs[p - d_particles.data()] = eelec;

  auto &d_estates = k_eivec.view_device()[d_ewhich[index_elecstate]].k_view.view_device();
  d_estates[p - d_particles.data()] = ielec;
  E_Dispose -= d_eelecs[p - d_particles.data()];
}

/* ----------------------------------------------------------------------
   reset the electronic state/energy of particle p to the ground state
   skip ambipolar electrons: they live in a separate scratch array (elist),
   not in d_particles, so they have no custom storage to reset
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::zero_elec(Particle::OnePart *p) const
{
  if (ambiflag && p->ispecies == ambispecies) return;
  auto &d_estates = k_eivec.view_device()[d_ewhich[index_elecstate]].k_view.view_device();
  auto &d_eelecs = k_edvec.view_device()[d_ewhich[index_eelec]].k_view.view_device();
  d_eelecs[p - d_particles.data()] = 0.0;
  d_estates[p - d_particles.data()] = 0;
}

/* ----------------------------------------------------------------------
   return 1 if electronic energy exchange is possible between two particles,
   i.e. discrete electronic modes are enabled and either species has
   electronic states defined
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
int CollideVSSKokkos::elec_exchange(Particle::OnePart *ip, Particle::OnePart *jp) const
{
  if (elecstyle != DISCRETE) return 0;
  if (d_nelecstates[ip->ispecies] > 0 || d_nelecstates[jp->ispecies] > 0)
    return 1;
  return 0;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
double CollideVSSKokkos::get_elec_phi(int ispec1, int ispec2, int ielec, double) const
{
  double species_rel = d_elec_species_rels(ispec1,ispec2,ielec);

  if (species_rel >= 0.0)
    return species_rel;
  else
    return d_elec_default_rels(ispec1,ielec);
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
int CollideVSSKokkos::select_elec_state(int icell,Particle::OnePart *p,
                                        Particle::OnePart *jp,
                                        double E_Dispose, double omega,
                                        bool enforce_spin_conservation,
                                        rand_type &rand_gen,
                                        bool reacting) const
{
  double State_prob;
  int max_level;
  auto &d_estates = k_eivec.view_device()[d_ewhich[index_elecstate]].k_view.view_device();
  // Find the maximum electronic level it can be in, given the current E_dispose
  max_level = 0;
  while (max_level < d_nelecstates[p->ispecies] && E_Dispose > d_elecstates(p->ispecies,max_level).temp*boltz) {
    ++max_level;
  }
  --max_level;

  // not enough energy to reach even the ground state: stay in ground state
  // (guards against max_level == -1 indexing d_state_probability[-1])
  if (max_level < 0) return 0;

  auto &d_state_probability = d_cumulative_probabilities;

  // Calculate number of total states, including degeneracies
  //
  // IMPORTANT: phi appears TWICE in a transition (gate on phi(current) in
  // the caller, weight by phi(candidate) below); both factors are required
  // for detailed balance with state-dependent relaxation numbers -- see the
  // matching comment in CollideVSS::select_elec_state
  for (int state = 0; state <= max_level; ++state) {
    if (state != 0) {
      d_state_probability(icell,state) = d_state_probability(icell,state-1);
    } else {
      d_state_probability(icell,state) = 0.0;
    }
    if (!enforce_spin_conservation ||
           d_elecstates(p->ispecies,state).spin == d_elecstates(p->ispecies,d_estates[p - d_particles.data()]).spin) {
      // Note we can use E_Dispose here since the current implementation requires the collision numbers
      // to be collision invariant (and therefore depend on E_Dispose, the trans + elec energy) but this
      // algorithm allows that to be relaxed. If other models are needed, the correct data would need passed
      // in here.
      if (reacting) {
        // We are distributing energy after a reaction, so relaxation probability is 100%
        d_state_probability(icell,state) += d_elecstates(p->ispecies,state).degen;
      } else {
        d_state_probability(icell,state) += d_elecstates(p->ispecies,state).degen*get_elec_phi(p->ispecies, jp->ispecies, state, E_Dispose);
      }
    }
  }
  // if no selectable state has any weight (e.g. every allowed state has a
  // zero relaxation probability), leave the particle in its current state

  if (d_state_probability(icell,max_level) <= 0.0)
    return d_estates[p - d_particles.data()];

  // Select a state from the distribution
  int ielec,ilast;
  double eelec = 0.0;
  do {
    double rand_state = rand_gen.drand()*d_state_probability(icell,max_level);
    ielec = 0;
    ilast = -1;
    // bound by max_level: roundoff can leave rand_state >= 0 after the last
    // included state, which would index d_elecstates past max_level/nelecstate
    while (rand_state >= 0 && ielec <= max_level) {
      if (!enforce_spin_conservation ||
             d_elecstates(p->ispecies,ielec).spin == d_elecstates(p->ispecies,d_estates[p - d_particles.data()]).spin) {
        if (reacting) {
          // We are distributing energy after a reaction, so relaxation probability is 100%
          rand_state -= d_elecstates(p->ispecies,ielec).degen;
        } else {
          rand_state -= d_elecstates(p->ispecies,ielec).degen*get_elec_phi(p->ispecies, jp->ispecies, ielec, E_Dispose);
        }
        ilast = ielec;
      }
      ++ielec;
    }
    // floating-point round-off can leave rand_state non-negative after all
    // weights are subtracted, so clamp to the last spin-allowed state
    ielec = ilast;
    eelec = d_elecstates(p->ispecies,ielec).temp*boltz;
    State_prob = pow((1.0 - eelec / E_Dispose),
                     (1.5 - omega));
  } while (State_prob < rand_gen.drand());
  return ielec;
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
  postcoln.eint = ip->erot + jp->erot + kp->erot
                + ip->evib + jp->evib + kp->evib;
  if (elecstyle == DISCRETE) {
    auto &d_eelecs = k_edvec.view_device()[d_ewhich[index_eelec]].k_view.view_device();
    if (d_nelecstates[ip->ispecies] > 0)
      postcoln.eint += d_eelecs[ip - d_particles.data()];
    if (d_nelecstates[jp->ispecies] > 0)
      postcoln.eint += d_eelecs[jp - d_particles.data()];
    if (d_nelecstates[kp->ispecies] > 0)
      postcoln.eint += d_eelecs[kp - d_particles.data()];
  }

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
void CollideVSSKokkos::EEXCHANGE_ReactingEDisposal(int icell,
                                                   Particle::OnePart *ip,
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

  // zero electronic state/energy of all products, not just those whose
  // species has electronic data: the reactant electronic energy is already
  // part of postcoln.etotal, so leaving a stale eelec on a product whose
  // (new) species has no electronic data would duplicate that energy

  if (!kp) {
    ip->erot = 0.0;
    jp->erot = 0.0;
    ip->evib = 0.0;
    jp->evib = 0.0;
    if (elecstyle == DISCRETE) {
      zero_elec(ip);
      zero_elec(jp);
    }
    numspecies = 2;
    aveomega = d_params(ip->ispecies,jp->ispecies).omega;
  } else {
    ip->erot = 0.0;
    jp->erot = 0.0;
    kp->erot = 0.0;
    ip->evib = 0.0;
    jp->evib = 0.0;
    kp->evib = 0.0;
    if (elecstyle == DISCRETE) {
      zero_elec(ip);
      zero_elec(jp);
      zero_elec(kp);
    }
    numspecies = 3;
    aveomega = (d_params(ip->ispecies,ip->ispecies).omega + d_params(jp->ispecies,jp->ispecies).omega +
                d_params(kp->ispecies,kp->ispecies).omega)/3.0;
  }

  // Phase 1: total effective internal DOF competing for the shared energy pool,
  // used to correct the Larsen-Borgnakke exponent for sequential sampling
  // (Dirichlet stick-breaking).  A discrete vibrational mode holds less energy
  // than a classical 2-DOF oscillator, so it is counted by its instantaneous
  // effective DOF zeta_m = eff_vib_dof(theta_m,Tcoll), evaluated at the
  // collision temperature Tcoll of the whole pool, found self-consistently from
  //   E = (2.5-aveomega + sum_classical_dof/2)*kB*Tcoll
  //       + sum_m kB*theta_m/(exp(theta_m/Tcoll) - 1).
  // Counting discrete modes as a static 2 DOF instead overstates the competing
  // pool and starves rotation of energy.

  double E_Dispose = postcoln.etotal;
  Particle::OnePart *plist[3] = {ip,jp,kp};

  double shape_classical = 2.5 - aveomega;   // translational shape (2.5-omega)
  double remaining_dof = 0.0;                // effective internal DOF left to draw
  int ndiscrete = 0;

  for (i = 0; i < numspecies; i++) {
    int sp = plist[i]->ispecies;
    if ((d_species[sp].rotdof > 0) && (rotstyle != NONE)) {
      shape_classical += 0.5 * d_species[sp].rotdof;
      remaining_dof += d_species[sp].rotdof;
    }
    if ((d_species[sp].vibdof > 0) && (vibstyle != NONE)) {
      if (vibstyle == DISCRETE) ndiscrete += d_species[sp].nvibmode;
      else {
        shape_classical += 0.5 * d_species[sp].vibdof;
        remaining_dof += d_species[sp].vibdof;
      }
    }
  }

  // collision temperature of the pool (classical unless discrete modes present)

  double tcoll = (shape_classical > 0.0) ? E_Dispose/(boltz*shape_classical) : 0.0;

  if (ndiscrete && E_Dispose > 0.0) {

    // flatten the discrete-mode frequencies once, skipping any theta <= 0
    // (a zero-frequency mode carries no energy and would make x/(exp(x)-1) NaN)

    double theta[3*Particle::MAXVIBMODE];
    int nflat = 0;
    for (i = 0; i < numspecies; i++) {
      int sp = plist[i]->ispecies;
      if ((d_species[sp].vibdof > 0) && (vibstyle == DISCRETE))
        for (int m = 0; m < d_species[sp].nvibmode; m++)
          if (d_species[sp].vibtemp[m] > 0.0) theta[nflat++] = d_species[sp].vibtemp[m];
    }

    // solve for the pool collision temperature, then add each discrete mode's
    // effective DOF at that temperature to the competing pool

    tcoll = vib_pool_temp(shape_classical,nflat,theta,E_Dispose);
    for (int m = 0; m < nflat; m++)
      remaining_dof += eff_vib_dof(theta[m],tcoll);
  }

  // Phase 2: Handle energy disposal for products with remaining_dof correction
  // to account for sequential sampling from shared pool (Dirichlet stick-breaking)

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
        double b_rot = (1.5 - aveomega) + 0.5 * (remaining_dof - rotdof);
        Fraction_Rot =
          1.0 - pow(rand_gen.drand(),(1.0/(1.0 + b_rot)));
        p->erot = Fraction_Rot * E_Dispose;
        E_Dispose -= p->erot;
        remaining_dof -= rotdof;

      } else if (rotdof > 2) {
        double b_rot = (1.5 - aveomega) + 0.5 * (remaining_dof - rotdof);
        p->erot = E_Dispose *
          sample_bl(rand_gen,0.5*d_species[sp].rotdof-1.0, b_rot);
        E_Dispose -= p->erot;
        remaining_dof -= rotdof;
      }
    }

    vibdof = d_species[sp].vibdof;

    if (vibdof) {
      if (vibstyle == NONE) {
        p->evib = 0.0;
      } else if (vibdof == 2 && vibstyle == DISCRETE) {
        double zeta = eff_vib_dof(d_species[sp].vibtemp[0],tcoll);
        double b_vib = (1.5 - aveomega) + 0.5 * (remaining_dof - zeta);
        max_level = static_cast<int>
          (E_Dispose / (boltz * d_species[sp].vibtemp[0]));
        do {
          ivib = static_cast<int>
            (rand_gen.drand()*(max_level+AdjustFactor));
          p->evib = (double)
            (ivib * boltz * d_species[sp].vibtemp[0]);
          State_prob = pow((1.0 - p->evib / E_Dispose), b_vib);
        } while (State_prob < rand_gen.drand());
        E_Dispose -= p->evib;
        remaining_dof -= zeta;

      } else if (vibdof == 2 && vibstyle == SMOOTH) {
        double b_vib = (1.5 - aveomega) + 0.5 * (remaining_dof - vibdof);
        Fraction_Vib =
          1.0 - pow(rand_gen.drand(),(1.0 / (1.0 + b_vib)));
        p->evib = Fraction_Vib * E_Dispose;
        E_Dispose -= p->evib;
        remaining_dof -= vibdof;

      } else if (vibdof > 2 && vibstyle == SMOOTH) {
        double b_vib = (1.5 - aveomega) + 0.5 * (remaining_dof - vibdof);
        p->evib = E_Dispose *
          sample_bl(rand_gen,0.5*d_species[sp].vibdof-1.0, b_vib);
        E_Dispose -= p->evib;
        remaining_dof -= vibdof;

      } else if (vibdof > 2 && vibstyle == DISCRETE) {
        p->evib = 0.0;

        int nmode = d_species[sp].nvibmode;
        const auto &d_vibmode = k_eiarray.view_device()[d_ewhich[index_vibmode]].k_view.view_device();
        int pindex = p - d_particles.data();

        for (int imode = 0; imode < nmode; imode++) {
          double zeta = eff_vib_dof(d_species[sp].vibtemp[imode],tcoll);
          max_level = static_cast<int>
          (E_Dispose / (boltz * d_species[sp].vibtemp[imode]));
          double b_vib = (1.5 - aveomega) + 0.5 * (remaining_dof - zeta);
          do {
            ivib = static_cast<int>
            (rand_gen.drand()*(max_level+AdjustFactor));
            pevib = ivib * boltz * d_species[sp].vibtemp[imode];
            State_prob = pow((1.0 - pevib / E_Dispose), b_vib);
          } while (State_prob < rand_gen.drand());

          d_vibmode(pindex,imode) = ivib;
          p->evib += pevib;
          E_Dispose -= pevib;
          remaining_dof -= zeta;
        }
      }
    }
    // use aveomega for the LB exponent, consistent with the rot/vib
    // redistribution above (the partner argument is only meaningful for
    // the relaxation-number lookup, which the reacting path skips)

    if (elecstyle == DISCRETE && d_nelecstates[sp] > 0)
      relax_electronic_mode(icell, p, p, E_Dispose, aveomega, rand_gen, true);
  }

  // compute post-collision internal energies

  postcoln.erot = ip->erot + jp->erot;
  postcoln.evib = ip->evib + jp->evib;
  postcoln.eelec = 0.0;
  if (elecstyle == DISCRETE) {
    auto &d_eelecs = k_edvec.view_device()[d_ewhich[index_eelec]].k_view.view_device();
    if (d_nelecstates[ip->ispecies] > 0)
      postcoln.eelec += d_eelecs[ip - d_particles.data()];
    if (d_nelecstates[jp->ispecies] > 0)
      postcoln.eelec += d_eelecs[jp - d_particles.data()];
  }

  if (kp) {
    postcoln.erot += kp->erot;
    postcoln.evib += kp->evib;
    if (elecstyle == DISCRETE) {
      auto &d_eelecs = k_edvec.view_device()[d_ewhich[index_eelec]].k_view.view_device();
      if (d_nelecstates[kp->ispecies] > 0)
        postcoln.eelec += d_eelecs[kp - d_particles.data()];
    }
  }

  // compute portion of energy left over for scattering

  postcoln.eint = postcoln.erot + postcoln.evib + postcoln.eelec;
  postcoln.etrans = E_Dispose;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
double CollideVSSKokkos::eff_vib_dof(double theta, double tcoll) const
{
  if (theta <= 0.0 || tcoll <= 0.0) return 0.0;
  double x = theta / tcoll;
  return 2.0 * x / (exp(x) - 1.0);
}

/* ----------------------------------------------------------------------
   collision temperature Tcoll of an energy pool E shared by shape_classical
   translational+classical-internal shape and nmode discrete SHO modes of
   characteristic temperatures theta[]:
     E = kB*( shape_classical*Tcoll + sum_m theta_m/(exp(theta_m/Tcoll)-1) )
   [0, E/(kB*shape_classical)] brackets the single root; solved with a
   safeguarded Newton iteration (bisection fallback) that converges
   quadratically in the typical case and cannot overshoot to a nonphysical
   temperature.  Requires shape_classical > 0 and E > 0 (caller guaranteed).
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
double CollideVSSKokkos::vib_pool_temp(double shape_classical, int nmode,
                                       double *theta, double E) const
{
  double Thi = E / (boltz * shape_classical);
  double Tlo = 0.0;
  double T = Thi;

  for (int iter = 0; iter < 30; iter++) {
    double f = boltz * shape_classical * T - E;
    double df = boltz * shape_classical;
    for (int m = 0; m < nmode; m++) {
      double x = theta[m] / T;
      if (x > 200.0) continue;             // frozen mode: exp overflow, ~0 term
      double ex = exp(x);
      double den = ex - 1.0;
      f  += boltz * theta[m] / den;
      df += boltz * theta[m]*theta[m] * ex / (T*T * den*den);
    }
    if (f > 0.0) Thi = T; else Tlo = T;    // keep [Tlo,Thi] bracketing the root
    double Tnew = T - f/df;                // Newton step
    if (!(Tnew > Tlo && Tnew < Thi))       // ... but stay inside the bracket
      Tnew = 0.5 * (Tlo + Thi);
    double delta = fabs(Tnew - T);
    T = Tnew;
    if (delta < 1.0e-4 * T) break;
  }
  return T;
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
   near neighbor search for a group pair
   mirrors Collide::find_nn_group() (collide.cpp:2634).  ni/nj are the two
     group counts; when ig == jg the host passes the same nn array for both,
     which is d_nn_igroup here
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
int CollideVSSKokkos::find_nn_group(rand_type &rand_gen, int icell, int i,
                                    int ig, int jg, int ni, int nj) const
{
  int jneigh;
  double dx,dy,dz,rsq;
  double *xj;

  const int same = (ig == jg);

  // if same group and nj = 2, just return J = non-I particle

  if (same && nj == 2) return (i+1) % 2;

  Particle::OnePart *ipart,*jpart;

  // thresh = distance particle I moves in this timestep

  ipart = &d_particles[d_plist(icell,d_glist(icell,ig,i))];
  double *vi = ipart->v;
  double *xi = ipart->x;
  double threshsq = dt*dt * (vi[0]*vi[0]+vi[1]*vi[1]+vi[2]*vi[2]);
  double minrsq = BIG;

  // nlimit = max # of J candidates to consider

  int nlimit = MIN(nearlimit,nj-1);
  int count = 0;

  // pick a random starting J
  // jneigh = collision partner when exit loop
  //   set to initial J as default in case no Nlimit J meets criteria

  int j = nj * rand_gen.drand();
  if (same)
    while (i == j) j = nj * rand_gen.drand();
  jneigh = j;

  while (count < nlimit) {
    count++;

    // skip this J if I,J last collided with each other

    const int nnj = same ? d_nn_igroup(icell,j) : d_nn_jgroup(icell,j);
    if (d_nn_igroup(icell,i) == j+1 && nnj == i+1) {
      j++;
      if (j == nj) j = 0;
      continue;
    }

    // rsq = squared distance between particles I and J
    // if rsq = 0.0, skip this J
    // if rsq <= threshsq, this J is collision partner
    // if rsq = smallest yet seen, this J is tentative collision partner

    jpart = &d_particles[d_plist(icell,d_glist(icell,jg,j))];
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
    if (j == nj) j = 0;
  }

  return jneigh;
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

    // no electron reactant: I/J order is not canonical if an ion is the
    // third body (e.g. AB + C+ -> A + C+ + B), so sync each product's
    // ion flag to its post-reaction species
    // also correct for all-neutral dissociation, where flags stay 0

    if (jsp != e) {
      d_ionambi[i] = d_ions[ip->ispecies];
      d_ionambi[j] = d_ions[jp->ispecies];
      d_ionambi[k] = d_ions[kp->ispecies];
      return;
    }

    d_ionambi[k] = 0;
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
    if (jsp == e) d_ionambi[i] = 0;   // R: A+ + e -> A, 1st product neutral
    else d_ionambi[i] = d_ions[ip->ispecies];  // sync product to its species
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

  this->sync(Host,ALL_MASK);

  int n = 0;
  if (memflag) {
    for (int igroup = 0; igroup < ngroups; igroup++) {
      for (int jgroup = 0; jgroup < ngroups; jgroup++) {
        buf[n++] = k_vremax.view_host()(icell,igroup,jgroup);
        if (remainflag)
          buf[n++] = k_remain.view_host()(icell,igroup,jgroup);
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
            buf[n++] = k_vremax.view_host()(m,igroup,jgroup);
            if (remainflag)
              buf[n++] = k_remain.view_host()(m,igroup,jgroup);
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

  this->sync(Host,ALL_MASK);

  int n = 0;
  for (int igroup = 0; igroup < ngroups; igroup++) {
    for (int jgroup = 0; jgroup < ngroups; jgroup++) {
      k_vremax.view_host()(icell,igroup,jgroup) = buf[n++];
      if (remainflag)
        k_remain.view_host()(icell,igroup,jgroup) = buf[n++];
    }
  }
  nglocal++;

  this->modified(Host,ALL_MASK);

  if (cells[icell].nsplit > 1) {
    int isplit = cells[icell].isplit;
    int nsplit = cells[icell].nsplit;
    grow_percell(nsplit);

    this->sync(Host,ALL_MASK);

    for (int i = 0; i < nsplit; i++) {
      int m = sinfo[isplit].csubs[i];
      for (int igroup = 0; igroup < ngroups; igroup++) {
        for (int jgroup = 0; jgroup < ngroups; jgroup++) {
          k_vremax.view_host()(m,igroup,jgroup) = buf[n++];
          if (remainflag)
            k_remain.view_host()(m,igroup,jgroup) = buf[n++];
        }
      }
    }
    nglocal += nsplit;

    this->modified(Host,ALL_MASK);
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
      k_vremax.view_host()(jcell,igroup,jgroup) = k_vremax.view_host()(icell,igroup,jgroup);
      if (remainflag)
        k_remain.view_host()(jcell,igroup,jgroup) = k_remain.view_host()(icell,igroup,jgroup);
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
      k_vremax.view_host()(nglocal,igroup,jgroup) = vremax_initial[igroup][jgroup];
      if (remainflag) k_remain.view_host()(nglocal,igroup,jgroup) = 0.0;
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
  k_vremax.resize(Kokkos::view_alloc(Kokkos::WithoutInitializing),
                  nglocalmax,ngroups,ngroups);
  d_vremax = k_vremax.view_device();
  if (remainflag) {
    k_remain.resize(Kokkos::view_alloc(Kokkos::WithoutInitializing),
                    nglocalmax,ngroups,ngroups);
    d_remain = k_remain.view_device();
  }
  this->sync(Host,ALL_MASK);
  for (int icell = nglocal_old; icell < nglocal; icell++)
    for (int igroup = 0; igroup < ngroups; igroup++)
      for (int jgroup = 0; jgroup < ngroups; jgroup++) {
        k_vremax.view_host()(icell,igroup,jgroup) = vremax_initial[igroup][jgroup];
        if (remainflag) k_remain.view_host()(icell,igroup,jgroup) = 0.0;
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

  k_vremax.resize(Kokkos::view_alloc(Kokkos::WithoutInitializing),
                  nglocalmax,ngroups,ngroups);
  d_vremax = k_vremax.view_device();
  if (remainflag) {
    k_remain.resize(Kokkos::view_alloc(Kokkos::WithoutInitializing),
                    nglocalmax,ngroups,ngroups);
    d_remain = k_remain.view_device();
  }

  this->modified(Device,ALL_MASK); // needed for auto sync
}

/* ---------------------------------------------------------------------- */

void CollideVSSKokkos::sync(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if (sparta->kokkos->auto_sync) {
      // Automatic syncing exists because non-Kokkos code may have written the
      // plain vremax/remain pointers, so the host is declared modified and
      // copied down.  Declaring it while the device still holds a claim is
      // both a lie -- the host copy is the older one -- and fatal: Kokkos
      // aborts a DualView claimed on both sides at once.  collisions() claims
      // the device every step and nothing calls sync(Host) in the run loop, so
      // the claim stands for the whole run; the first sync(Device) made with
      // auto_sync on then aborts.  fix adapt is exactly that caller -- it sets
      // kokkos_flag = 0, which makes ModifyKokkos turn auto_sync on around
      // end_of_step(), and CollideVSSKokkos::grow_percell() syncs to the
      // device from inside AdaptGrid::perform_refine().  Refresh the host
      // first: a no-op when the device is clean, and the copy the host is owed
      // when it is not.
      sync(Host,mask);
      modified(Host,mask);
    }
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

  // custom per-particle arrays mutated during collisions must also be
  // backed up, else a react/retry pass re-runs on top of modified values

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  auto h_ewhich = particle_kk->k_ewhich.view_host();

  if (vibstyle == DISCRETE && index_vibmode >= 0) {
    auto d_vibmode = particle_kk->k_eiarray.view_host()[h_ewhich[index_vibmode]].k_view.view_device();
    d_vibmode_backup = decltype(d_vibmode_backup)(Kokkos::view_alloc("collide:vibmode_backup",Kokkos::WithoutInitializing),d_vibmode.extent(0),d_vibmode.extent(1));
    Kokkos::deep_copy(d_vibmode_backup,d_vibmode);
  }

  if (elecstyle == DISCRETE) {
    auto d_eelec = particle_kk->k_edvec.view_host()[h_ewhich[index_eelec]].k_view.view_device();
    auto d_elecstate = particle_kk->k_eivec.view_host()[h_ewhich[index_elecstate]].k_view.view_device();
    d_eelec_backup = decltype(d_eelec_backup)(Kokkos::view_alloc("collide:eelec_backup",Kokkos::WithoutInitializing),d_eelec.extent(0));
    d_elecstate_backup = decltype(d_elecstate_backup)(Kokkos::view_alloc("collide:elecstate_backup",Kokkos::WithoutInitializing),d_elecstate.extent(0));
    Kokkos::deep_copy(d_eelec_backup,d_eelec);
    Kokkos::deep_copy(d_elecstate_backup,d_elecstate);
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
  Kokkos::deep_copy(particle_kk->k_particles.view_device(),d_particles_backup);
  d_particles = particle_kk->k_particles.view_device();

  GridKokkos* grid_kk = (GridKokkos*) grid;
  Kokkos::deep_copy(grid_kk->d_plist,d_plist_backup);
  d_plist = grid_kk->d_plist;

  Kokkos::deep_copy(d_vremax,d_vremax_backup);
  Kokkos::deep_copy(d_remain,d_remain_backup);

  if (ambiflag) {
    auto h_ewhich = particle_kk->k_ewhich.view_host();

    Kokkos::deep_copy(particle_kk->k_eivec.view_host()[h_ewhich[index_ionambi]].k_view.view_device(),d_ionambi_backup);
    Kokkos::deep_copy(particle_kk->k_edarray.view_host()[h_ewhich[index_velambi]].k_view.view_device(),d_velambi_backup);

    k_eivec = particle_kk->k_eivec;
    k_edarray = particle_kk->k_edarray;
    d_ionambi = k_eivec.view_host()[h_ewhich[index_ionambi]].k_view.view_device();
    d_velambi = k_edarray.view_host()[h_ewhich[index_velambi]].k_view.view_device();
  }

  if (vibstyle == DISCRETE && index_vibmode >= 0) {
    auto h_ewhich = particle_kk->k_ewhich.view_host();
    Kokkos::deep_copy(particle_kk->k_eiarray.view_host()[h_ewhich[index_vibmode]].k_view.view_device(),d_vibmode_backup);
    k_eiarray = particle_kk->k_eiarray;
  }

  if (elecstyle == DISCRETE) {
    auto h_ewhich = particle_kk->k_ewhich.view_host();
    Kokkos::deep_copy(particle_kk->k_edvec.view_host()[h_ewhich[index_eelec]].k_view.view_device(),d_eelec_backup);
    Kokkos::deep_copy(particle_kk->k_eivec.view_host()[h_ewhich[index_elecstate]].k_view.view_device(),d_elecstate_backup);
    k_eivec = particle_kk->k_eivec;
    k_edvec = particle_kk->k_edvec;
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

  d_particles_backup = {};
  d_plist_backup = {};
  d_vremax_backup = {};
  d_remain_backup = {};

  if (ambiflag) {
    d_ionambi_backup = {};
    d_velambi_backup = {};
  }

  d_vibmode_backup = {};
  d_eelec_backup = {};
  d_elecstate_backup = {};
}

/* ----------------------------------------------------------------------
   grow every per-event gas tally compute past what the failed attempt
     needed; see the same helper in UpdateKokkos
------------------------------------------------------------------------- */

void CollideVSSKokkos::grow_gas_tally_computes()
{
  int ncoll = 0, nreact = 0;

  for (int i = 0; i < ngas_tally; i++) {
    Compute *c = update->glist_active[i];
    if (ComputeGasCollisionTallyKokkos *ckk = dynamic_cast<ComputeGasCollisionTallyKokkos*>(c)) {
      ckk->grow_after_overflow();

      // growing reallocated the compute's row buffer, so the copy the kernel
      //   reads still points at the old, too-small one.  Without re-blitting
      //   it the repeated attempt overflows on the same row and the retry
      //   loop never terminates

#ifdef SPARTA_KOKKOS_FIXED_LISTS
      glist_coll_tally_copy[ncoll++].copy(ckk);
#else
      gas_buf_blit(k_glist_coll_tally,ncoll++,ckk);
#endif
    } else if (ComputeGasReactionTallyKokkos *ckk = dynamic_cast<ComputeGasReactionTallyKokkos*>(c)) {
      ckk->grow_after_overflow();
#ifdef SPARTA_KOKKOS_FIXED_LISTS
      glist_react_tally_copy[nreact++].copy(ckk);
#else
      gas_buf_blit(k_glist_react_tally,nreact++,ckk);
#endif
    }
  }

#ifndef SPARTA_KOKKOS_FIXED_LISTS
  gas_buf_sync(k_glist_coll_tally,d_glist_coll_tally);
  gas_buf_sync(k_glist_react_tally,d_glist_react_tally);
#endif
}

/* ----------------------------------------------------------------------
   mark (mark=1) or rewind to (mark=0) the append position of every per-event
     gas tally compute
   an aborted attempt of the retry loop has to take back the rows it appended
     before the pass re-runs; clear_gas_tally() only resets the per-grid
     computes, so it does not cover these
------------------------------------------------------------------------- */

void CollideVSSKokkos::rewind_gas_tally_computes(int mark)
{
  for (int i = 0; i < ngas_tally; i++) {
    Compute *c = update->glist_active[i];
    if (ComputeGasCollisionTallyKokkos *ckk =
          dynamic_cast<ComputeGasCollisionTallyKokkos*>(c))
      { if (mark) ckk->mark_ntally(); else ckk->rewind_ntally(); }
    else if (ComputeGasReactionTallyKokkos *ckk =
               dynamic_cast<ComputeGasReactionTallyKokkos*>(c))
      { if (mark) ckk->mark_ntally(); else ckk->rewind_ntally(); }
  }
}

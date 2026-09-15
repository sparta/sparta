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

#include "spatype.h"
#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "update_kokkos.h"
#include "math_const.h"
#include "particle_kokkos.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "domain.h"
#include "comm_kokkos.h"
#include "collide.h"
#include "collide_vss_kokkos.h"
#include "grid_kokkos.h"
#include "surf_kokkos.h"
#include "surf_collide.h"
#include "surf_react.h"
#include "output.h"
#include "geometry_kokkos.h"
#include "random_mars.h"
#include "timer.h"
#include "math_extra.h"
#include "memory_kokkos.h"
#include "error.h"
#include <unistd.h>
#include "kokkos.h"
#include "sparta_masks.h"
#include "surf_collide_specular_kokkos.h"
#include "kokkos_base.h"

using namespace SPARTA_NS;

enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};         // same as Domain
enum{PERIODIC,OUTFLOW,REFLECT,SURFACE,AXISYM};  // same as Domain
//enum{OUTSIDE,INSIDE,ONSURF2OUT,ONSURF2IN};      // several files
enum{PKEEP,PINSERT,PDONE,PDISCARD,PENTRY,PEXIT,PSURF};   // several files
enum{NCHILD,NPARENT,NUNKNOWN,NPBCHILD,NPBPARENT,NPBUNKNOWN,NBOUND};  // Grid
enum{TALLYAUTO,TALLYREDUCE,TALLYLOCAL};         // same as Surf
enum{PERAUTO,PERCELL,PERSURF};                  // several files
enum{NOFIELD,CFIELD,PFIELD,GFIELD};             // several files
enum{BCSTD,BCWRAP,BCMIRROR,BCEXIT};             // Update::bcopt values

#define MAXSTUCK 20
#define EPSPARAM 1.0e-7

// either set ID or PROC/INDEX, set other to -1

//#define MOVE_DEBUG 1              // un-comment to debug one particle
#define MOVE_DEBUG_ID 308143534  // particle ID
#define MOVE_DEBUG_PROC -1        // owning proc
#define MOVE_DEBUG_INDEX -1   // particle index on owning proc
#define MOVE_DEBUG_STEP 4107    // timestep

#define VAL_1(X) X
#define VAL_2(X) VAL_1(X), VAL_1(X)

/* ----------------------------------------------------------------------
   blit one active tally compute into its per-type device buffer
   same operation and same rationale as KKCopy::copy() (kokkos_copy.h:71):
     the object is only read on device, through KOKKOS_INLINE_FUNCTION
     members, so its vtable pointer is never used and the View handles it
     carries stay alive in the original the compute list holds
------------------------------------------------------------------------- */

#ifndef SPARTA_KOKKOS_FIXED_LISTS
namespace {

  template<class T>
  void tally_buf_resize(DAT::tdual_char_1d &k, DAT::t_char_1d &d, int n)
  {
    const size_t need = (size_t) MAX(n,1) * sizeof(T);
    if (k.view_device().extent(0) < need) {
      k = DAT::tdual_char_1d("update:tally_models",need);
      d = k.view_device();
    }
  }

  template<class T>
  void tally_buf_blit(DAT::tdual_char_1d &k, int slot, T *obj)
  {
    char *dst = k.view_host().data() + (size_t) slot*sizeof(T);
    memcpy((void*) dst, (const void*) obj, sizeof(T));
    ((T *) dst)->copy = 1;
  }

  void tally_buf_sync(DAT::tdual_char_1d &k, DAT::t_char_1d &d)
  {
    if (k.view_device().extent(0) == 0) return;
    k.modify_host();
    k.sync_device();
    d = k.view_device();
  }
}
#endif

/* ---------------------------------------------------------------------- */


/* ---------------------------------------------------------------------- */

UpdateKokkos::UpdateKokkos(SPARTA *sparta) : Update(sparta),
  grid_kk_copy(sparta),
  domain_kk_copy(sparta)
#ifdef SPARTA_KOKKOS_FIXED_LISTS
  // Virtual functions are not yet supported on the GPU, which leads to pain:
  , slist_active_copy{VAL_2(KKCopy<ComputeSurfKokkos>(sparta))}
  , slist_active_isurf_copy{VAL_2(KKCopy<ComputeISurfGridKokkos>(sparta))}
  , slist_active_coll_tally_copy{VAL_2(KKCopy<ComputeSurfCollisionTallyKokkos>(sparta))}
  , slist_active_react_tally_copy{VAL_2(KKCopy<ComputeSurfReactionTallyKokkos>(sparta))}
  , slist_active_react_isurf_copy{VAL_2(KKCopy<ComputeReactISurfGridKokkos>(sparta))}
  , slist_active_react_surf_copy{VAL_2(KKCopy<ComputeReactSurfKokkos>(sparta))}
  , blist_active_copy{VAL_2(KKCopy<ComputeBoundaryKokkos>(sparta))}
  , blist_active_react_copy{VAL_2(KKCopy<ComputeReactBoundaryKokkos>(sparta))}
  , tmp_compute_boundary_kk(sparta)
  , tmp_compute_react_boundary_kk(sparta)
  , tmp_compute_surf_kk(sparta)
  , tmp_compute_isurf_grid_kk(sparta)
  , tmp_compute_react_isurf_grid_kk(sparta)
  , tmp_compute_react_surf_kk(sparta)
#endif
{
  nslist_surf = nslist_isurf = nslist_react_isurf = nslist_react_surf = 0;
  nslist_coll_tally = nslist_react_tally = 0;
  nsc_index_cached = -1;

  // the Kokkos views of Particle/Grid/Surf are populated from the host data
  //   once, by setup() when prewrap is set, which then clears prewrap
  // a "clear" command destroys and recreates all of these classes but not
  //   KokkosSPARTA, so prewrap has to be re-armed here or the second problem
  //   runs with empty device views
  // this is a no-op on the first construction, KokkosSPARTA sets prewrap = 1

  sparta->kokkos->prewrap = 1;

  // use 1D views for scalars to reduce GPU memory operations
  // int view = flags and view-index counters, must stay int
  // bigint view = per-step statistics counters, can exceed 2^31
  //   in one step at large per-proc particle counts

  d_scalars = t_int_7("update:scalars");
  h_scalars = t_host_int_7("update:scalars_mirror");

  d_scalars_big = t_bigint_7("update:scalars_big");
  h_scalars_big = t_host_bigint_7("update:scalars_big_mirror");

  d_nmigrate      = Kokkos::subview(d_scalars,0);
  d_entryexit     = Kokkos::subview(d_scalars,1);
  d_nstuck        = Kokkos::subview(d_scalars,2);
  d_naxibad       = Kokkos::subview(d_scalars,3);
  d_error_flag    = Kokkos::subview(d_scalars,4);
  d_retry         = Kokkos::subview(d_scalars,5);
  d_nlocal        = Kokkos::subview(d_scalars,6);
  d_tally_overflow = Kokkos::subview(d_scalars,7);

  d_ncomm_one     = Kokkos::subview(d_scalars_big,0);
  d_nexit_one     = Kokkos::subview(d_scalars_big,1);
  d_nboundary_one = Kokkos::subview(d_scalars_big,2);
  d_ntouch_one    = Kokkos::subview(d_scalars_big,3);
  d_nscheck_one   = Kokkos::subview(d_scalars_big,4);
  d_nscollide_one = Kokkos::subview(d_scalars_big,5);
  d_nreact_one    = Kokkos::subview(d_scalars_big,6);

  h_nmigrate      = Kokkos::subview(h_scalars,0);
  h_entryexit     = Kokkos::subview(h_scalars,1);
  h_nstuck        = Kokkos::subview(h_scalars,2);
  h_naxibad       = Kokkos::subview(h_scalars,3);
  h_error_flag    = Kokkos::subview(h_scalars,4);
  h_retry         = Kokkos::subview(h_scalars,5);
  h_nlocal        = Kokkos::subview(h_scalars,6);
  h_tally_overflow = Kokkos::subview(h_scalars,7);

  h_ncomm_one     = Kokkos::subview(h_scalars_big,0);
  h_nexit_one     = Kokkos::subview(h_scalars_big,1);
  h_nboundary_one = Kokkos::subview(h_scalars_big,2);
  h_ntouch_one    = Kokkos::subview(h_scalars_big,3);
  h_nscheck_one   = Kokkos::subview(h_scalars_big,4);
  h_nscollide_one = Kokkos::subview(h_scalars_big,5);
  h_nreact_one    = Kokkos::subview(h_scalars_big,6);

  nboundary_tally = 0;

  for (int f = 0; f < 6; f++) bcopt[f] = BCSTD;

  d_bcmirror = DAT::t_bigint_1d("update:bcmirror",6);
  h_bcmirror = HAT::t_bigint_1d("update:bcmirror_mirror",6);
}

/* ---------------------------------------------------------------------- */

UpdateKokkos::~UpdateKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_mlist,mlist);
  mlist = NULL;
}

/* ---------------------------------------------------------------------- */

void UpdateKokkos::init()
{
  // the surf_collide style list is fixed within a run but can change between
  //   them, so force setup_surf_collide_models() to rebuild its index maps

  nsc_index_cached = -1;

  // init the UpdateKokkos class if performing a run, else just return
  // only set first_update if a run is being performed

  if (runflag == 0) return;
  first_update = 1;

  if (optmove_flag) {
    if (!grid->uniform)
      error->all(FLERR,"Cannot use optimized move with non-uniform grid");
    else if (surf->exist)
      error->all(FLERR,"Cannot use optimized move when surfaces are defined");
    else {
      for (int ifix = 0; ifix < modify->nfix; ifix++) {
        if (strstr(modify->fix[ifix]->style,"adapt") != NULL)
          error->all(FLERR,"Cannot use optimized move with fix adapt");
      }
    }

    // the dense cell index is built by rehash(), which skips it unless
    //   optmove is on, so build it here in case optmove was turned on after
    //   the last rehash.  during a run rehash() keeps it in step

    grid->update_halo_index();
  }

  optmove_surf_init();

  // choose the appropriate move method

  // REACT=1 is also needed without explicit surfs when box-face/boundary
  //   reactions are defined (e.g. surf_react adsorb in face mode)

  if (domain->dimension == 3) {
    if (surf->exist) {
      if (surf->nsr) moveptr = &UpdateKokkos::move<3,1,1,0>;
      else moveptr = &UpdateKokkos::move<3,1,0,0>;
    } else {
      if (surf->nsr) moveptr = &UpdateKokkos::move<3,0,1,0>;
      else if (optmove_flag) moveptr = &UpdateKokkos::move<3,0,0,1>;
      else moveptr = &UpdateKokkos::move<3,0,0,0>;
    }
  } else if (domain->axisymmetric) {
    if (surf->exist) {
      if (surf->nsr) moveptr = &UpdateKokkos::move<1,1,1,0>;
      else moveptr = &UpdateKokkos::move<1,1,0,0>;
    } else {
      if (surf->nsr) moveptr = &UpdateKokkos::move<1,0,1,0>;
      else if (optmove_flag) moveptr = &UpdateKokkos::move<1,0,0,1>;
      else moveptr = &UpdateKokkos::move<1,0,0,0>;
    }
  } else if (domain->dimension == 2) {
    if (surf->exist) {
      if (surf->nsr) moveptr = &UpdateKokkos::move<2,1,1,0>;
      else moveptr = &UpdateKokkos::move<2,1,0,0>;
    } else {
      if (surf->nsr) moveptr = &UpdateKokkos::move<2,0,1,0>;
      else if (optmove_flag) moveptr = &UpdateKokkos::move<2,0,0,1>;
      else moveptr = &UpdateKokkos::move<2,0,0,0>;
    }
  }

  // checks on external field options

  if (fstyle == CFIELD) {
    if (domain->dimension == 2 && field[2] != 0.0)
      error->all(FLERR,"External field in z not allowed for 2d");
    if (domain->axisymmetric && field[1] != 0.0)
      error->all(FLERR,
                 "External field in y not allowed for axisymmetric model");
  } else if (fstyle == PFIELD) {
    ifieldfix = modify->find_fix(fieldID);
    if (ifieldfix < 0) error->all(FLERR,"External field fix ID not found");
    if (!modify->fix[ifieldfix]->per_particle_field)
      error->all(FLERR,"External field fix does not compute necessary field");
  } else if (fstyle == GFIELD) {
    ifieldfix = modify->find_fix(fieldID);
    if (ifieldfix < 0) error->all(FLERR,"External field fix ID not found");
    if (!modify->fix[ifieldfix]->per_grid_field)
      error->all(FLERR,"External field fix does not compute necessary field");
  }

  if (optmove_flag) {
    xlo = domain->boxlo[0];
    ylo = domain->boxlo[1];
    zlo = domain->boxlo[2];
    xhi = domain->boxhi[0];
    yhi = domain->boxhi[1];
    zhi = domain->boxhi[2];
    Lx = xhi-xlo;
    Ly = yhi-ylo;
    Lz = zhi-zlo;
    ncx = grid->unx;
    ncy = grid->uny;
    ncz = grid->unz;
    dx = Lx/ncx;
    dy = Ly/ncy;
    dz = Lz/ncz;
  }

  if (fstyle == PFIELD) {
    field_active[0] = modify->fix[ifieldfix]->field_active[0];
    field_active[1] = modify->fix[ifieldfix]->field_active[1];
    field_active[2] = modify->fix[ifieldfix]->field_active[2];
    KKBaseFieldFix = dynamic_cast<KokkosBase*>(modify->fix[ifieldfix]);
    if (!KKBaseFieldFix)
      error->all(FLERR,"External field fix is not Kokkos-enabled");
  } else if (fstyle == GFIELD) {
    field_active[0] = modify->fix[ifieldfix]->field_active[0];
    field_active[1] = modify->fix[ifieldfix]->field_active[1];
    field_active[2] = modify->fix[ifieldfix]->field_active[2];
    KKBaseFieldFix = dynamic_cast<KokkosBase*>(modify->fix[ifieldfix]);
    if (!KKBaseFieldFix)
      error->all(FLERR,"External field fix is not Kokkos-enabled");
  }
}

/* ---------------------------------------------------------------------- */

void UpdateKokkos::setup()
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  GridKokkos* grid_kk = (GridKokkos*) grid;
  SurfKokkos* surf_kk = (SurfKokkos*) surf;

  particle_kk->sync(Device,ALL_MASK);
  particle_kk->sorted_kk = 0;

  if (sparta->kokkos->prewrap) {

    // particle

    particle_kk->wrap_kokkos();

    // grid

    grid_kk->wrap_kokkos();
    grid_kk->update_hash();

    // surf

    if (surf->exist)
      surf_kk->wrap_kokkos();

    sparta->kokkos->prewrap = 0;
  } else {
    grid_kk->modify(Host,ALL_MASK);
    grid_kk->update_hash();

    if (surf->exist) {
      surf_kk->modify(Host,ALL_MASK);
      grid_kk->wrap_kokkos_graphs();
    }
  }
  grid_index_refresh();

  // device grid/surf graphs are now current; clear any pending change flag so
  // the run loop does not do a redundant resync on the first step
  grid->changed = 0;

  Update::setup(); // must come after prewrap since computes are called by setup()

  // For MPI debugging
  //
  //  volatile int i = 0;
  //  char hostname[256];
  //  gethostname(hostname, sizeof(hostname));
  //  printf("PID %d on %s ready for attach, i = %i\n", getpid(), hostname, i);
  //  fflush(stdout);
  //  sleep(30);
  //  printf("Continuing...\n");
}

/* ----------------------------------------------------------------------
   take the cell lookups used by the optimized move from the grid
   both are rebuilt from scratch by GridKokkos::update_hash(), which fix
     balance calls mid-run, so the copies held here have to be retaken
     whenever the cell views are, not once at setup
------------------------------------------------------------------------- */

void UpdateKokkos::grid_index_refresh()
{
  GridKokkos* grid_kk = (GridKokkos*) grid;

  hash_kk = grid_kk->hash_kk;

  d_halo_index = grid_kk->d_halo_index;
  halo_ilo = grid_kk->halo_ilo; halo_nx = grid_kk->halo_nx;
  halo_jlo = grid_kk->halo_jlo; halo_ny = grid_kk->halo_ny;
  halo_klo = grid_kk->halo_klo; halo_nz = grid_kk->halo_nz;
}

/* ---------------------------------------------------------------------- */

void UpdateKokkos::run(int nsteps)
{
  sparta->kokkos->auto_sync = 0;

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;

  int n_start_of_step = modify->n_start_of_step;
  int n_end_of_step = modify->n_end_of_step;

  // external per grid cell field
  // only evaluate once at beginning of run b/c time-independent
  // fix calculates field acting at center point of all grid cells

  if (fstyle == GFIELD && fieldfreq == 0) {
    modify->fix[ifieldfix]->compute_field();
    d_fieldfix_array_grid = KKBaseFieldFix->d_array_grid;
  }

  // cellweightflag = 1 if grid-based particle weighting is ON

  int cellweightflag = 0;
  if (grid->cellweightflag) cellweightflag = 1;

  // loop over timesteps

  for (int i = 0; i < nsteps; i++) {

    if (timer->check_timeout(i)) {
      update->nsteps = i;
      break;
    }

    ntimestep++;

    if (collide_react) collide_react_reset();
    if (tallyflag) tally_set(ntimestep);
    if (dynamic) dynamic_update();

    timer->stamp();

    // start of step fixes

    if (n_start_of_step) {
      modify->start_of_step();
      timer->stamp(TIME_MODIFY);
    }

    // establish surf-tally compute copies for the move kernel here, after
    //   start-of-step fixes have run.  A fix such as fix emit/surf performs
    //   its own surf-tally session during start_of_step that reallocates the
    //   shared dup_array_surf_tally scatter views, so they must be recreated
    //   now to give the move kernel a live, freshly zeroed scatter view.

    if (tallyflag) setup_surf_tally_copies();

    // move particles

    if (cellweightflag) particle->pre_weight();
    (this->*moveptr)();
    timer->stamp(TIME_MOVE);

    // communicate particles

    if (nmigrate) {
      k_mlist_small = Kokkos::subview(k_mlist,std::make_pair(0,nmigrate));
      k_mlist_small.sync_host();
    }
    auto mlist_small = k_mlist_small.view_host().data();

    ((CommKokkos*)comm)->migrate_particles(nmigrate,mlist_small,k_mlist_small.view_device());
    if (cellweightflag) particle->post_weight();
    timer->stamp(TIME_COMM);

    const int reorder_flag = (update->reorder_period &&
        (update->ntimestep % update->reorder_period == 0));

    if (collide || reorder_flag) {
      particle_kk->sort_kokkos();
      timer->stamp(TIME_SORT);
    }

    if (collide) {
      collide->collisions();
      timer->stamp(TIME_COLLIDE);
    }

    if (collide_react) collide_react_update();

    // diagnostic fixes

    if (n_end_of_step) {
      modify->end_of_step();
      timer->stamp(TIME_MODIFY);
    }

    // if an end-of-step fix changed the grid/surf topology (e.g. fix ablate
    // regenerated implicit surfaces), the host grid is now authoritative but
    // the device per-cell surf graphs (d_csurfs/d_csplits/d_csubs) are stale.
    // Resync them to the device before the next move, mirroring setup().
    // Safe here: grid_kk_copy from this step's move is no longer in use and is
    // refreshed at the start of the next move.

    // safety net: a grid change from anywhere other than the end-of-step
    //   batch (which resyncs per fix in ModifyKokkos) is handled here before
    //   the next move reads the device grid

    if (grid->changed) ((GridKokkos*) grid)->resync_after_host_change();

    // all output

    if (ntimestep == output->next) {
      particle_kk->sync(Host,ALL_MASK);
      output->write(ntimestep);
      timer->stamp(TIME_OUTPUT);
    }
  }

  modify->post_run();

  sparta->kokkos->auto_sync = 1;
  particle_kk->sync(Host,ALL_MASK);
}

/* ----------------------------------------------------------------------
   advect particles thru grid
   DIM = 2/3 for 2d/3d, 1 for 2d axisymmetric
   SURF = 0/1 for no surfs or surfs
   use multiple iterations of move/comm if necessary
------------------------------------------------------------------------- */

template < int DIM, int SURF, int REACT, int OPT > void UpdateKokkos::move()
{
  int pstart,pstop,entryexit,any_entryexit;
  int continue_loop_flag = 0;

  // extend migration list if necessary

  int maxlocal = particle->maxlocal;

  if (particle->nlocal > maxmigrate) {
    maxmigrate = maxlocal;
    memoryKK->destroy_kokkos(k_mlist,mlist);
    memoryKK->create_kokkos(k_mlist,mlist,maxmigrate,"particle:mlist");
  }

  // counters

  niterate = 0;
  ntouch_one = ncomm_one = 0;
  nboundary_one = nexit_one = 0;
  nscheck_one = nscollide_one = 0;
  surf->nreact_one = 0;

  if (!sparta->kokkos->need_atomics || sparta->kokkos->atomic_reduction) {
    h_ntouch_one() = 0;
    h_nexit_one() = 0;
    h_nboundary_one() = 0;
    h_ncomm_one() = 0;
    h_nscheck_one() = 0;
    h_nscollide_one() = 0;
    h_nreact_one() = 0;
  }

  h_error_flag() = 0;

  // move/migrate iterations

  dt = update->dt;

  // which global boundary faces the fast path may handle itself, see the OPT
  //   block in the move kernel below
  // a compute boundary tallies every kind of crossing, and only the standard
  //   move calls the tally, so give the optimization up on any step where one
  //   is active.  nboundary_tally is set per step by tally_set()
  // a surface face is left to the standard move, since it runs a collision
  //   model.  an outflow face only deletes the particle, so the fast path
  //   does that itself

  if (OPT) {
    for (int f = 0; f < 6; f++) {
      if (nboundary_tally) bcopt[f] = BCSTD;
      else if (domain->bflag[f] == PERIODIC) bcopt[f] = BCWRAP;
      else if (domain->bflag[f] == REFLECT) bcopt[f] = BCMIRROR;
      else if (domain->bflag[f] == OUTFLOW) bcopt[f] = BCEXIT;
      else if (bcmirror_surf[f]) bcopt[f] = BCMIRROR;
      else bcopt[f] = BCSTD;
    }
    // axisymmetric: a mirror at the outer radial face is not a mirror in the
    //   (x,r) plane.  reflecting off that cylinder turns the particle in 3d,
    //   and the radial path after the turn is not the continuation of the one
    //   before it, which is what mirroring r about the face would assume --
    //   the error is percent-level, not round-off.  leave it to the standard
    //   move, which walks to the face and reflects there
    // outflow at that face is still exact: r(t)^2 is a parabola in t, so r is
    //   unimodal, and a particle that starts inside can cross the face only
    //   once, upward.  ending outside therefore means it left
    // the axis itself needs nothing: axi_remap() returns r >= 0 = boxlo[1],
    //   so the fast path never sees a particle below it

    if (domain->axisymmetric && bcopt[YHI] == BCMIRROR) bcopt[YHI] = BCSTD;


    if (bcmirror_any) Kokkos::deep_copy(d_bcmirror,0);
  }

  ParticleKokkos* particle_kk = ((ParticleKokkos*)particle);

  // external per particle field
  // fix calculates field acting on all owned particles

  if (fstyle == PFIELD) {
    modify->fix[ifieldfix]->compute_field();
    d_fieldfix_array_particle = KKBaseFieldFix->d_array_particle;
  }

  // external per grid cell field
  // evaluate once every fieldfreq steps b/c time-dependent
  // fix calculates field acting at center point of all grid cells

  if (fstyle == GFIELD && fieldfreq && ((ntimestep-1) % fieldfreq == 0)) {
    modify->fix[ifieldfix]->compute_field();
    d_fieldfix_array_grid = KKBaseFieldFix->d_array_grid;
  }

  // one or more loops over particles
  // first iteration = all my particles
  // subsequent iterations = received particles

  while (1) {

    if (!continue_loop_flag)
      niterate++;

    d_particles = particle_kk->k_particles.view_device();

    GridKokkos* grid_kk = ((GridKokkos*)grid);
    d_cells = grid_kk->k_cells.view_device();
    d_sinfo = grid_kk->k_sinfo.view_device();
    d_pcells = grid_kk->k_pcells.view_device();

    // GridKokkos::update_hash() builds a brand new UnorderedMap rather than
    //   updating in place, so a copy taken at setup() goes stale as soon as
    //   anything rehashes mid-run (fix adapt/balance/move surf, or the
    //   grid/surf resync after fix ablate regenerates implicit surfaces).
    //   The optimized-move kernel looks up cell IDs in it, so refresh it here
    //   with the other grid handles or it maps IDs to pre-change cell indices

    hash_kk = grid_kk->hash_kk;

    d_csurfs = grid_kk->d_csurfs;
    d_csplits = grid_kk->d_csplits;
    d_csubs = grid_kk->d_csubs;

    // the cell lookups are refreshed here alongside the cell views, not left
    //   at what setup() copied: fix balance rebuilds them mid-run, into fresh
    //   views, so a copy taken once at setup goes stale after the first
    //   rebalance while d_cells above does not

    grid_index_refresh();

    if (surf->exist) {
      SurfKokkos* surf_kk = ((SurfKokkos*)surf);
      surf_kk->sync(Device,ALL_MASK);
      d_lines = surf_kk->k_lines.view_device();
      d_tris = surf_kk->k_tris.view_device();
    }

    if (surf->nsr) {
      double extra_factor = 1.0;
      if (!sparta->kokkos->react_retry_flag)
        extra_factor = sparta->kokkos->react_extra;

      // compute in bigint and guard: the double->int conversion of
      //   nlocal*extra_factor is UB once it exceeds 2^31

      bigint nlocal_extra = static_cast<bigint> (particle->nlocal*extra_factor);
      if (nlocal_extra > MAXSMALLINT)
        error->one(FLERR,"Per-processor particle count is too big");
      if ((bigint) d_particles.extent(0) < nlocal_extra) {
        particle->grow(nlocal_extra - particle->nlocal); // this!
        d_particles = particle_kk->k_particles.view_device();
      }
    }

    particle_kk->sync(Device,PARTICLE_MASK);
    grid_kk->sync(Device,CELL_MASK|PCELL_MASK|SINFO_MASK|PLEVEL_MASK);

    // may be able to move this outside of the while loop
    grid_kk_copy.copy(grid_kk);
    domain_kk_copy.copy((DomainKokkos*)domain);

    setup_surf_collide_models();

    Kokkos::deep_copy(h_scalars,0);
    Kokkos::deep_copy(h_scalars_big,0);

    if (!continue_loop_flag) {
      nmigrate = 0;
      entryexit = 0;
    }

    if (niterate == 1 && !continue_loop_flag) {
      pstart = 0;
      pstop = particle->nlocal;
    }

    UPDATE_REDUCE reduce;

    // Reactions may create or delete more particles than existing views can hold.
    //  Cannot grow a Kokkos view in a parallel loop, so
    //  if the capacity of the view is exceeded, break out of parallel loop,
    //  reallocate on the host, and then repeat the parallel loop again.
    //  Unfortunately this leads to really messy code.

    h_retry() = 1;

    // a per-event tally compute can force a retry of its own, and a retry
    //   re-runs the move over the same particles.  that is only sound if the
    //   particle list can be rolled back first, so the backup is not gated on
    //   react/retry when one of those computes is active: without it the
    //   second attempt would move already-moved particles

    const int tally_backup = (nslist_coll_tally || nslist_react_tally);
    const int do_backup =
      (surf->nsr && sparta->kokkos->react_retry_flag) || tally_backup;

    // rows already tallied by earlier migration iterations of this step stay;
    //   an attempt of this iteration takes back only its own

    if (tally_backup) rewind_tally_computes(1);

    while (h_retry()) {

      if (do_backup) backup();

      // discard the rows an aborted attempt appended, including an attempt
      //   repeated for a reaction overflow rather than a tally overflow

      if (tally_backup) rewind_tally_computes(0);

      h_retry() = 0;
      h_nlocal() = particle->nlocal;
      if (continue_loop_flag) h_nmigrate() = nmigrate;

      Kokkos::deep_copy(d_scalars,h_scalars);
      Kokkos::deep_copy(d_scalars_big,h_scalars_big);

      // zero the custom attributes of the slots a surf reaction can fill
      // must precede the kernel, not follow it: SurfCollide calls
      //   update_custom_kokkos() for a particle the reaction just created
      // repeated on each retry, since a rolled back attempt leaves values
      //   behind in those slots

      if (surf->nsr) particle_kk->zero_custom_kokkos();

      copymode = 1;

    /* ATOMIC_REDUCTION: 1 = use atomics
                         0 = don't need atomics
                        -1 = use parallel_reduce
    */

#if defined SPARTA_KOKKOS_GPU
  #if SPARTA_KOKKOS_REDUCE_ARCH
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagUpdateMove<DIM,SURF,REACT,OPT,-1> >(pstart,pstop),*this,reduce);
  #else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagUpdateMove<DIM,SURF,REACT,OPT,1> >(pstart,pstop),*this);
  #endif
#elif defined KOKKOS_ENABLE_SERIAL
      if constexpr(std::is_same<DeviceType,Kokkos::Serial>::value)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagUpdateMove<DIM,SURF,REACT,OPT,0> >(pstart,pstop),*this);
      else {
        if (!sparta->kokkos->need_atomics)
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagUpdateMove<DIM,SURF,REACT,OPT,0> >(pstart,pstop),*this);
        else
          Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagUpdateMove<DIM,SURF,REACT,OPT,-1> >(pstart,pstop),*this,reduce);
      }
#else
      if (!sparta->kokkos->need_atomics)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagUpdateMove<DIM,SURF,REACT,OPT,0> >(pstart,pstop),*this);
      else
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagUpdateMove<DIM,SURF,REACT,OPT,-1> >(pstart,pstop),*this,reduce);
#endif

      copymode = 0;

      Kokkos::deep_copy(h_scalars,d_scalars);
      Kokkos::deep_copy(h_scalars_big,d_scalars_big);

      // a per-event surf tally compute ran out of room.  the row count is
      //   only knowable by running the move, so grow every such compute to
      //   what this attempt actually needed and repeat, exactly as a
      //   reaction overflow does.  unlike a reaction overflow this needs no
      //   react/extra opt-in: nothing about the particle state forced it,
      //   and truncating a tally would silently corrupt dump tally output

      if (h_tally_overflow() && !h_retry()) {
        grow_tally_computes();
        if (do_backup) restore();
        Kokkos::deep_copy(h_scalars,0);
        Kokkos::deep_copy(h_scalars_big,0);
        reduce = UPDATE_REDUCE();
        h_retry() = 1;
        continue;
      }

      if (h_retry()) {
        int nlocal_new = h_nlocal();

        if (!do_backup) {
          error->one(FLERR,"Ran out of space for Kokkos reactions, increase react/extra"
                           " or use react/retry");
        } else
          restore();

        //  reset counters

        Kokkos::deep_copy(h_scalars,0);
        Kokkos::deep_copy(h_scalars_big,0);
        reduce = UPDATE_REDUCE();
        h_retry() = 1;

        if (d_particles.extent(0) < nlocal_new) {
          particle->grow(nlocal_new - particle->nlocal);
          d_particles = particle_kk->k_particles.view_device();
        }
      }
    }

    particle_kk->modify(Device,PARTICLE_MASK);
    d_particles = t_particle_1d(); // destroy reference to reduce memory use

    k_mlist.modify_device();

    // END of pstart/pstop loop advecting all particles

    nmigrate = h_nmigrate();

    particle->nlocal = h_nlocal();

    int error_flag;

    if (!sparta->kokkos->need_atomics || sparta->kokkos->atomic_reduction) {
      ntouch_one += h_ntouch_one();
      nexit_one += h_nexit_one();
      nboundary_one += h_nboundary_one();
      ncomm_one += h_ncomm_one();
      nscheck_one += h_nscheck_one();
      nscollide_one += h_nscollide_one();
      surf->nreact_one += h_nreact_one();
      nstuck += h_nstuck();
      naxibad += h_naxibad();
    } else {
      ntouch_one       += reduce.ntouch_one   ;
      nexit_one        += reduce.nexit_one    ;
      nboundary_one    += reduce.nboundary_one;
      ncomm_one        += reduce.ncomm_one    ;
      nscheck_one      += reduce.nscheck_one  ;
      nscollide_one    += reduce.nscollide_one;
      surf->nreact_one += reduce.nreact_one   ;
      nstuck           += reduce.nstuck       ;
      naxibad          += reduce.naxibad      ;
    }

    entryexit += h_entryexit();

    error_flag = h_error_flag();

    if (error_flag) {
      char str[128];
      snprintf(str,sizeof(str),
              "Particle being sent to self proc "
              "on step " BIGINT_FORMAT,
              update->ntimestep);
      error->one(FLERR,str);
    }

    for (int n = 0; n < surf->nsc; n++) sc_phase(surf->sc[n],SC_POST);

    // move newly created particles from surface reactions

    continue_loop_flag = 0;

    if (surf->nsr && pstop < particle->nlocal) {
      pstart = pstop;
      pstop = particle->nlocal;
      continue_loop_flag = 1;
      continue;
    }

    // if gridcut >= 0.0, check if another iteration of move is required
    // only the case if some particle flag = PENTRY/PEXIT
    //   in which case perform particle migration
    // if not, move is done and final particle comm will occur in run()
    // if iterating, reset pstart/pstop and extend migration list if necessary

    if (grid->cutoff < 0.0) break;

    timer->stamp(TIME_MOVE);
    MPI_Allreduce(&entryexit,&any_entryexit,1,MPI_INT,MPI_MAX,world);
    timer->stamp(TIME_SYNC);

    if (any_entryexit) {
      if (nmigrate) {
        k_mlist_small = Kokkos::subview(k_mlist,std::make_pair(0,nmigrate));
        k_mlist_small.sync_host();
      }
      auto mlist_small = k_mlist_small.view_host().data();
      timer->stamp(TIME_MOVE);
      pstart = ((CommKokkos*)comm)->migrate_particles(nmigrate,mlist_small,k_mlist_small.view_device());
      timer->stamp(TIME_COMM);
      pstop = particle->nlocal;
      if (pstop-pstart > maxmigrate) {
        maxmigrate = pstop-pstart;
        memoryKK->destroy_kokkos(k_mlist,mlist);
        memoryKK->create_kokkos(k_mlist,mlist,maxmigrate,"particle:mlist");
      }
    } else break;

    // END of single move/migrate iteration
  }

  // END of all move/migrate iterations

  // the retry-loop particle backup is reused across this step's migration
  //   iterations; release it now so peak memory matches the old behaviour

  free_particle_backup();

  particle->sorted = 0;
  particle_kk->sorted_kk = 0;

  // hand any {s} face mirrors the fast path did back to their collide models

  if (OPT && bcmirror_any) {
    Kokkos::deep_copy(h_bcmirror,d_bcmirror);
    for (int f = 0; f < 6; f++) bcmirror_one[f] = h_bcmirror[f];
    optmove_surf_tally();
  }

  // accumulate running totals

  niterate_running += niterate;
  nmove_running += particle->nlocal;
  ntouch_running += ntouch_one;
  ncomm_running += ncomm_one;
  nboundary_running += nboundary_one;
  nexit_running += nexit_one;
  nscheck_running += nscheck_one;
  nscollide_running += nscollide_one;
  surf->nreact_running += surf->nreact_one;

  // dispatch by dynamic_cast, not by style string, and in the same order as
  //   setup_surf_tally_copies(): the styles are also registered under explicit
  //   "/kk" names (e.g. react/isurf/grid/kk), so a style compare would miss a
  //   compute the user typed with the suffix and fall through to a wrong cast

  if (nsurf_tally) {
    for (int m = 0; m < nsurf_tally; m++) {
      if (ComputeISurfGridKokkos* compute_isurf_kk =
            dynamic_cast<ComputeISurfGridKokkos*>(slist_active[m])) {
        compute_isurf_kk->post_surf_tally();
      } else if (ComputeReactISurfGridKokkos* compute_react_isurf_kk =
                   dynamic_cast<ComputeReactISurfGridKokkos*>(slist_active[m])) {
        compute_react_isurf_kk->post_surf_tally();
      } else if (ComputeReactSurfKokkos* compute_react_surf_kk =
                   dynamic_cast<ComputeReactSurfKokkos*>(slist_active[m])) {
        compute_react_surf_kk->post_surf_tally();
      } else if (ComputeSurfKokkos* compute_surf_kk =
                   dynamic_cast<ComputeSurfKokkos*>(slist_active[m])) {
        compute_surf_kk->post_surf_tally();
      } else if (ComputeSurfCollisionTallyKokkos* compute_ct_kk =
                   dynamic_cast<ComputeSurfCollisionTallyKokkos*>(slist_active[m])) {
        compute_ct_kk->post_surf_tally();
      } else if (ComputeSurfReactionTallyKokkos* compute_rt_kk =
                   dynamic_cast<ComputeSurfReactionTallyKokkos*>(slist_active[m])) {
        compute_rt_kk->post_surf_tally();
      } else {
        error->all(FLERR,"Kokkos does not (yet) support this surf tally compute; "
                         "use a Kokkos-enabled surf tally compute (-sf kk)");
      }
    }
  }

  // dispatch by dynamic_cast for the same reason as the surf tally list above,
  //   and because compute boundary and compute react/boundary both set
  //   boundary_tally_flag but are unrelated classes: a static cast would call
  //   one's methods on the other

  if (nboundary_tally) {
    for (int m = 0; m < nboundary_tally; m++) {
      if (ComputeBoundaryKokkos* c =
            dynamic_cast<ComputeBoundaryKokkos*>(blist_active[m]))
        c->post_boundary_tally();
      else if (ComputeReactBoundaryKokkos* c =
                 dynamic_cast<ComputeReactBoundaryKokkos*>(blist_active[m]))
        c->post_boundary_tally();
      else
        error->all(FLERR,"Kokkos does not (yet) support this boundary tally compute; "
                         "use a Kokkos-enabled boundary tally compute (-sf kk)");
    }
  }
}

/* ----------------------------------------------------------------------
   first step of the optimized move: apply the global boundary condition to an
     end-of-step position
   xnew = position at the end of a straight-line move, not modified
   xp = xnew after any boundary condition, valid only if this returns 1
   flip = bit 0/1/2 per dimension whose velocity component a mirror negated,
     plus bit 3/4/5 when it was that dimension's upper face, so the caller can
     name the face it reflected off.  the caller applies it only once it
     decides to keep the particle: one that falls through must reach the
     standard move with its velocity untouched, which is also why xnew is left
     alone (and why xnew+L-L will not do, since that is not xnew in floating
     point)
   bcopt[face] says what this face may do -- BCSTD nothing, BCWRAP translate,
     BCMIRROR mirror, BCEXIT delete -- so a face the fast path cannot handle
     leaves the position outside the box and the bound tests below reject it
   one translation or mirror per dimension covers any particle that did not
     cross a whole domain in a single step; anything still outside is rejected
   return 1 if xp is inside the global box, 0 to use the standard move,
     -1 if the particle left through an outflow face and is to be deleted
------------------------------------------------------------------------- */

template < int DIM >
KOKKOS_INLINE_FUNCTION
int UpdateKokkos::optmove_bc(const double *xnew, double *xp, int &flip) const
{
  xp[0] = xnew[0];
  xp[1] = xnew[1];
  xp[2] = xnew[2];
  flip = 0;

  // exitbit records, per dimension, that the face the particle went out of is
  //   an outflow face.  the position is left alone for those, so the bound
  //   tests below still see the dimension as outside and can tell the two
  //   reasons for that apart

  int exitbit = 0;

  if (xp[0] < xlo) {
    if (bcopt[XLO] == BCWRAP) xp[0] += Lx;
    else if (bcopt[XLO] == BCMIRROR) { xp[0] = xlo + (xlo-xp[0]); flip |= 1; }
    else if (bcopt[XLO] == BCEXIT) exitbit |= 1;
  } else if (xp[0] >= xhi) {
    if (bcopt[XHI] == BCWRAP) xp[0] -= Lx;
    else if (bcopt[XHI] == BCMIRROR) { xp[0] = xhi - (xp[0]-xhi); flip |= 1|8; }
    else if (bcopt[XHI] == BCEXIT) exitbit |= 1;
  }

  if (xp[1] < ylo) {
    if (bcopt[YLO] == BCWRAP) xp[1] += Ly;
    else if (bcopt[YLO] == BCMIRROR) { xp[1] = ylo + (ylo-xp[1]); flip |= 2; }
    else if (bcopt[YLO] == BCEXIT) exitbit |= 2;
  } else if (xp[1] >= yhi) {
    if (bcopt[YHI] == BCWRAP) xp[1] -= Ly;
    else if (bcopt[YHI] == BCMIRROR) { xp[1] = yhi - (xp[1]-yhi); flip |= 2|16; }
    else if (bcopt[YHI] == BCEXIT) exitbit |= 2;
  }

  if (DIM == 3) {
    if (xp[2] < zlo) {
      if (bcopt[ZLO] == BCWRAP) xp[2] += Lz;
      else if (bcopt[ZLO] == BCMIRROR) { xp[2] = zlo + (zlo-xp[2]); flip |= 4; }
      else if (bcopt[ZLO] == BCEXIT) exitbit |= 4;
    } else if (xp[2] >= zhi) {
      if (bcopt[ZHI] == BCWRAP) xp[2] -= Lz;
      else if (bcopt[ZHI] == BCMIRROR) { xp[2] = zhi - (xp[2]-zhi); flip |= 4|32; }
      else if (bcopt[ZHI] == BCEXIT) exitbit |= 4;
    }
  }

  // cell bounds are half open, [lo,hi), so a particle exactly on an upper face
  // belongs to no cell and has to be rejected too.  with > instead of >= it
  // would reach the lookup with an index of ncx/ncy/ncz and alias onto an
  // unrelated cell
  //
  // a dimension still outside because of an outflow face means the particle
  //   left the domain, but only if every other dimension is accounted for: a
  //   face the fast path does not handle runs a surface collision model that
  //   can turn the particle around before it ever reaches the outflow face, so
  //   one of those anywhere sends the particle to the standard move instead.
  //   a wrap or a mirror cannot, since neither changes the motion in the
  //   dimension that exits

  int exited = 0;

  if (xp[0] < xlo || xp[0] >= xhi) {
    if (!(exitbit & 1)) return 0;
    exited = 1;
  }
  if (xp[1] < ylo || xp[1] >= yhi) {
    if (!(exitbit & 2)) return 0;
    exited = 1;
  }
  if (DIM == 3)
    if (xp[2] < zlo || xp[2] >= zhi) {
      if (!(exitbit & 4)) return 0;
      exited = 1;
    }

  // a particle that mirrored off one face and left through another cannot be
  //   deleted here.  whether the standard move counts that boundary collision
  //   depends on which face the particle reached first, and this does not
  //   determine that: reaching the mirror first reflects it and then it still
  //   leaves, since a mirror in one dimension does not change the motion in
  //   the one it exits through, but reaching the outflow face first means the
  //   reflection never happened.  the fate is the same either way, the tally
  //   is not, so hand it to the standard move
  // a wrap alongside an exit is fine and stays here, since a periodic crossing
  //   tallies nothing

  if (exited) return flip ? 0 : -1;

  return 1;
}

/* ----------------------------------------------------------------------
   second step of the optimized move: map a position inside the global box to
     the local index of the cell holding it
   caller must have established that xp is inside the box, via optmove_bc()
   preferred path is one indexed load into d_halo_index, keyed on the cell's
     position within this proc's halo arc.  the conditional adds fold a
     periodically wrapped ghost layer back into the arc
   return the local cell index, or -1 if this proc does not hold that cell, in
     which case the caller uses the standard move
------------------------------------------------------------------------- */

template < int DIM >
KOKKOS_INLINE_FUNCTION
int UpdateKokkos::optmove_cell(const double *xp) const
{
  const int ip = static_cast<int> ((xp[0] - xlo)/dx);
  const int jp = static_cast<int> ((xp[1] - ylo)/dy);
  int kp = 0;
  if (DIM == 3) kp = static_cast<int> ((xp[2] - zlo)/dz);

  // optmove_bc() has already put xp inside the box, but dx is a rounded
  //   quotient, so a position within an ulp of an upper face can still divide
  //   to ncx.  that index is not a miss to be caught by the lookup: the cell
  //   ID it forms is the valid ID of the first cell of the next row, on the
  //   far side of the box, and the hash would return it.  reject it here, for
  //   both lookups

  if (ip >= ncx || jp >= ncy || kp >= ncz) return -1;

  if (d_halo_index.extent(0)) {
    int il = ip - halo_ilo; if (il < 0) il += ncx;
    int jl = jp - halo_jlo; if (jl < 0) jl += ncy;
    int kl = kp - halo_klo; if (kl < 0) kl += ncz;
    if (il < halo_nx && jl < halo_ny && kl < halo_nz)
      return d_halo_index[((size_t) kl*halo_ny + jl)*halo_nx + il];
    return -1;
  }

  // no dense index for this decomposition: hash on the global cell ID
  // must accumulate in cellint, since ncx/ncy/ncz are int and an int
  //   expression overflows once the global cell count passes 2^31, silently
  //   disabling this fast path for every particle above that point in the grid

  const cellint cellIdx = ((cellint) kp*ncy + jp)*ncx + ip + 1;
  auto index = hash_kk.find(static_cast<GridKokkos::key_type>(cellIdx));
  if (hash_kk.valid_at(index)) return static_cast<int> (hash_kk.value_at(index));

  return -1;
}

/* ---------------------------------------------------------------------- */

template<int DIM, int SURF, int REACT, int OPT, int ATOMIC_REDUCTION>
KOKKOS_INLINE_FUNCTION
void UpdateKokkos::operator()(TagUpdateMove<DIM,SURF,REACT,OPT,ATOMIC_REDUCTION>, const int &i) const {
  UPDATE_REDUCE reduce;
  this->template operator()<DIM,SURF,REACT,OPT,ATOMIC_REDUCTION>(TagUpdateMove<DIM,SURF,REACT,OPT,ATOMIC_REDUCTION>(), i, reduce);
}

/*-----------------------------------------------------------------------------*/

template<int DIM, int SURF, int REACT, int OPT, int ATOMIC_REDUCTION>
KOKKOS_INLINE_FUNCTION
void UpdateKokkos::operator()(TagUpdateMove<DIM,SURF,REACT,OPT,ATOMIC_REDUCTION>, const int &i, UPDATE_REDUCE &reduce) const {
  if (d_error_flag() || d_retry()) return;

  // int m;
  bool hitflag;
  int icell,icell_original,outface,bflag,nflag,pflag,itmp;
  int side,minsurf,nsurf,cflag,isurf,exclude,stuck_iterate;
  double dtremain,frac,newfrac,param,minparam,rnew,dtsurf,tc,tmp;
  double xnew[3],xhold[3],xc[3],vc[3],minxc[3],minvc[3];
  double *x,*v;
  Surf::Tri *tri;
  Surf::Line *line;
  int reaction;

  Particle::OnePart &particle_i = d_particles[i];
  pflag = particle_i.flag;

  Particle::OnePart iorig;
  Particle::OnePart *ipart,*jpart;
  jpart = NULL;

  // received from another proc and move is done
  // if first iteration, PDONE is from a previous step,
  //   set pflag to PKEEP so move the particle on this step
  // else do nothing

  if (pflag == PDONE) {
    pflag = particle_i.flag = PKEEP;
    if (niterate > 1) return;
  }

  x = particle_i.x;
  v = particle_i.v;
  exclude = -1;

  // for 2d and axisymmetry only
  // xnew,xc passed to geometry routines which use or set z component

  if (DIM < 3) xnew[2] = xc[2] = 0.0;

  // apply moveperturb() to PKEEP and PINSERT since are computing xnew
  // not to PENTRY,PEXIT since are just re-computing xnew of sender
  // set xnew[2] to linear move for axisymmetry, will be remapped later
  // let pflag = PEXIT persist to check during axisymmetric cell crossing

  if (DIM < 3) xnew[2] = 0.0;
  if (pflag == PKEEP) {
    dtremain = dt;
    xnew[0] = x[0] + dtremain*v[0];
    xnew[1] = x[1] + dtremain*v[1];
    if (DIM != 2) xnew[2] = x[2] + dtremain*v[2];
    if (fstyle == CFIELD) {
      // DIM == 1 is the axisymmetric model, which the host treats as 2d here:
      //   Update::init() selects field2d on domain->dimension == 2, which is
      //   true for axisymmetric.  Do not narrow this back to DIM == 2
      if (DIM == 3) field3d(dtremain,xnew,v);
      else field2d(dtremain,xnew,v);
    } else if (fstyle == PFIELD) field_per_particle(i,particle_i.icell,dtremain,xnew,v);
    else if (fstyle == GFIELD) field_per_grid(i,particle_i.icell,dtremain,xnew,v);
  } else if (pflag == PINSERT) {
    dtremain = particle_i.dtremain;
    xnew[0] = x[0] + dtremain*v[0];
    xnew[1] = x[1] + dtremain*v[1];
    if (DIM != 2) xnew[2] = x[2] + dtremain*v[2];
    if (fstyle == CFIELD) {
      // DIM == 1 is the axisymmetric model, which the host treats as 2d here:
      //   Update::init() selects field2d on domain->dimension == 2, which is
      //   true for axisymmetric.  Do not narrow this back to DIM == 2
      if (DIM == 3) field3d(dtremain,xnew,v);
      else field2d(dtremain,xnew,v);
    } else if (fstyle == PFIELD) field_per_particle(i,particle_i.icell,dtremain,xnew,v);
    else if (fstyle == GFIELD) field_per_grid(i,particle_i.icell,dtremain,xnew,v);
  } else if (pflag == PENTRY) {
    icell = particle_i.icell;
    if (d_cells[icell].nsplit > 1) {
      if (DIM == 3 && SURF) icell = split3d(icell,x);
      if (DIM < 3 && SURF) icell = split2d(icell,x);
      particle_i.icell = icell;
    }
    dtremain = particle_i.dtremain;
    xnew[0] = x[0] + dtremain*v[0];
    xnew[1] = x[1] + dtremain*v[1];
    if (DIM != 2) xnew[2] = x[2] + dtremain*v[2];
  } else if (pflag == PEXIT) {
    dtremain = particle_i.dtremain;
    xnew[0] = x[0] + dtremain*v[0];
    xnew[1] = x[1] + dtremain*v[1];
    if (DIM != 2) xnew[2] = x[2] + dtremain*v[2];
  } else if (pflag >= PSURF) {
    dtremain = particle_i.dtremain;
    xnew[0] = x[0] + dtremain*v[0];
    xnew[1] = x[1] + dtremain*v[1];
    if (DIM != 2) xnew[2] = x[2] + dtremain*v[2];
    if (pflag > PSURF) exclude = pflag - PSURF - 1;
  }

  // optimized move
  // resolve the particle's end-of-step cell in one step instead of walking the
  //   grid cell by cell.  optmove_bc() applies whatever global boundary
  //   condition is needed and optmove_cell() maps the resulting position to a
  //   local cell index; if either declines, the particle falls through to the
  //   standard move below with its state untouched

  if (OPT) {
    double xp[3];
    int flip;

    // axisymmetry: fold the linear end-of-step position back into the (x,r)
    //   plane before anything else looks at it.  one remap of the whole step
    //   is enough, and gives the same answer as the standard move's remap at
    //   every cell crossing: each remap is a rotation about the x axis applied
    //   to position and velocity together, so the trajectory is unchanged and
    //   only the frame moves, and r is invariant under it.  the intermediate
    //   remaps are there so the cell-by-cell walk can follow the curve in
    //   (x,r), which the fast path does not need to do
    // remap a copy: axi_remap() rotates the velocity, and a particle that
    //   falls through has to reach the standard move with v untouched

    const double *xin = xnew;
    double xaxi[3],vaxi[3];

    if (DIM == 1) {
      xaxi[0] = xnew[0]; xaxi[1] = xnew[1]; xaxi[2] = xnew[2];
      vaxi[0] = v[0];    vaxi[1] = v[1];    vaxi[2] = v[2];
      axi_remap(xaxi,vaxi);
      xin = xaxi;
    }

    const int bc = optmove_bc<DIM>(xin,xp,flip);

    // left through an outflow face: the standard move would walk it to the
    //   face and delete it there, which is the same particle gone and the
    //   same counter, so do it here
    // a discarded particle still goes on the migrate list, since that is what
    //   deletes it -- migration drops a PDISCARD rather than sending it.  it
    //   is not counted in ncomm_one, which counts particles sent

    if (bc < 0) {
      particle_i.flag = PDISCARD;

      int indx;
      if (ATOMIC_REDUCTION == 0) {
        indx = d_nmigrate();
        d_nmigrate()++;
      } else {
        indx = Kokkos::atomic_fetch_add(&d_nmigrate(),1);
      }
      k_mlist.view_device()[indx] = i;

      if (ATOMIC_REDUCTION == 1)
        Kokkos::atomic_inc(&d_nexit_one());
      else if (ATOMIC_REDUCTION == 0)
        d_nexit_one()++;
      else
        reduce.nexit_one++;
      return;
    }

    if (bc) {
      const int icell = optmove_cell<DIM>(xp);

      if (icell >= 0) {

        // reset particle cell and coordinates

        particle_i.icell = icell;
        particle_i.flag = PKEEP;
        x[0] = xp[0];
        x[1] = xp[1];
        x[2] = xp[2];

        // axisymmetry: the particle is committed, so the rotated velocity from
        //   the remap becomes its velocity.  x[2] is 0, which xp already
        //   carries through from the remap

        if (DIM == 1) {
          v[0] = vaxi[0];
          v[1] = vaxi[1];
          v[2] = vaxi[2];
        }

        // specular reflection off a global boundary: now that the particle is
        //   committed, negate the velocity components the mirrors flipped and
        //   count the boundary collisions the standard move would have counted
        // a particle that reaches two faces in one step reflects off both, and
        //   mirroring one dimension does not change the motion in the other,
        //   so the two are independent and each is tallied

        if (flip) {
          int nb = 0;
          if (flip & 1) { v[0] = -v[0]; nb++; }
          if (flip & 2) { v[1] = -v[1]; nb++; }
          if (flip & 4) { v[2] = -v[2]; nb++; }
          if (ATOMIC_REDUCTION == 1)
            Kokkos::atomic_add(&d_nboundary_one(),nb);
          else if (ATOMIC_REDUCTION == 0)
            d_nboundary_one() += nb;
          else
            reduce.nboundary_one += nb;

          // an {s} face's surf collide model keeps its own count of the
          //   collisions it handled, so count them per face here and give them
          //   to it after the kernel.  a {r} face has no model, and this is
          //   skipped entirely unless some face is an {s} one, since every
          //   reflecting particle would otherwise contend on the same counter

          if (bcmirror_any) {
            if (flip & 1)
              Kokkos::atomic_inc(&d_bcmirror[(flip & 8) ? XHI : XLO]);
            if (flip & 2)
              Kokkos::atomic_inc(&d_bcmirror[(flip & 16) ? YHI : YLO]);
            if (flip & 4)
              Kokkos::atomic_inc(&d_bcmirror[(flip & 32) ? ZHI : ZLO]);
          }
        }

        if (d_cells[icell].proc != me) {
          int indx;
          if (ATOMIC_REDUCTION == 0) {
            indx = d_nmigrate();
            d_nmigrate()++;
          } else {
            indx = Kokkos::atomic_fetch_add(&d_nmigrate(),1);
          }
          k_mlist.view_device()[indx] = i;

          particle_i.flag = PDONE;

          if (ATOMIC_REDUCTION == 1)
            Kokkos::atomic_inc(&d_ncomm_one());
          else if (ATOMIC_REDUCTION == 0)
            d_ncomm_one()++;
          else
            reduce.ncomm_one++;
        }

        return;
      }
    }
  }

  particle_i.flag = PKEEP;
  icell = particle_i.icell;
  double* lo = d_cells[icell].lo;
  double* hi = d_cells[icell].hi;
  cellint* neigh = d_cells[icell].neigh;
  int nmask = d_cells[icell].nmask;
  stuck_iterate = 0;
  if (ATOMIC_REDUCTION == 1)
    Kokkos::atomic_inc(&d_ntouch_one());
  else if (ATOMIC_REDUCTION == 0)
    d_ntouch_one()++;
  else
    reduce.ntouch_one++;

  // advect one particle from cell to cell and thru surf collides til done

  while (1) {

#ifdef MOVE_DEBUG
    if (DIM == 3) {
      if (ntimestep == MOVE_DEBUG_STEP &&
          (MOVE_DEBUG_ID == d_particles[i].id ||
           (me == MOVE_DEBUG_PROC && i == MOVE_DEBUG_INDEX)))
        printf("PARTICLE %d %ld: %d %d: %d: x %g %g %g: xnew %g %g %g: %d "
               CELLINT_FORMAT ": lo %g %g %g: hi %g %g %g: DTR %g\n",
               me,ntimestep,i,d_particles[i].id,
               d_cells[icell].nsurf,
               x[0],x[1],x[2],xnew[0],xnew[1],xnew[2],
               icell,d_cells[icell].id,
               lo[0],lo[1],lo[2],hi[0],hi[1],hi[2],dtremain);
    }
    if (DIM == 2) {
      if (ntimestep == MOVE_DEBUG_STEP &&
          (MOVE_DEBUG_ID == d_particles[i].id ||
           (me == MOVE_DEBUG_PROC && i == MOVE_DEBUG_INDEX)))
        printf("PARTICLE %d %ld: %d %d: %d: x %g %g: xnew %g %g: %d "
               CELLINT_FORMAT ": lo %g %g: hi %g %g: DTR: %g\n",
               me,ntimestep,i,d_particles[i].id,
               d_cells[icell].nsurf,
               x[0],x[1],xnew[0],xnew[1],
               icell,d_cells[icell].id,
               lo[0],lo[1],hi[0],hi[1],dtremain);
    }
    if (DIM == 1) {
      if (ntimestep == MOVE_DEBUG_STEP &&
          (MOVE_DEBUG_ID == d_particles[i].id ||
           (me == MOVE_DEBUG_PROC && i == MOVE_DEBUG_INDEX)))
        printf("PARTICLE %d %ld: %d %d: %d: x %g %g: xnew %g %g: %d "
               CELLINT_FORMAT ": lo %g %g: hi %g %g: DTR: %g\n",
               me,ntimestep,i,d_particles[i].id,
               d_cells[icell].nsurf,
               x[0],x[1],xnew[0],sqrt(xnew[1]*xnew[1]+xnew[2]*xnew[2]),
               icell,d_cells[icell].id,
               lo[0],lo[1],hi[0],hi[1],dtremain);
    }
#endif

    // check if particle crosses any cell face
    // frac = fraction of move completed before hitting cell face
    // this section should be as efficient as possible,
    //   since most particles won't do anything else
    // axisymmetric cell face crossings:
    //   use linear xnew to check vertical faces
    //   must always check move against curved lower y face of cell
    //   use remapped rnew to check horizontal lines
    //   for y faces, if pflag = PEXIT, particle was just received
    //     from another proc and is exiting this cell from face:
    //       axi_horizontal_line() will not detect correct crossing,
    //       so set frac and outface directly to move into adjacent cell,
    //       then unset pflag so not checked again for this particle

    outface = INTERIOR;
    frac = 1.0;

    if (xnew[0] < lo[0]) {
      if (xnew[0] != x[0]) frac = (lo[0]-x[0]) / (xnew[0]-x[0]);
      else frac = 0.0;
      if (frac < 0.0) frac = 0.0;
      else if (frac > 1.0) frac = 1.0;
      outface = XLO;
    } else if (xnew[0] >= hi[0]) {
      if (xnew[0] != x[0]) frac = (hi[0]-x[0]) / (xnew[0]-x[0]);
      else frac = 0.0;
      if (frac < 0.0) frac = 0.0;
      else if (frac > 1.0) frac = 1.0;
      outface = XHI;
    }

    if (DIM != 1) {
      if (xnew[1] < lo[1]) {
        if (xnew[1] != x[1]) newfrac = (lo[1]-x[1]) / (xnew[1]-x[1]);
        else newfrac = 0.0;
        if (newfrac < 0.0) newfrac = 0.0;
        else if (newfrac > 1.0) newfrac = 1.0;
        if (newfrac < frac) {
          frac = newfrac;
          outface = YLO;
        }
      } else if (xnew[1] >= hi[1]) {
        if (xnew[1] != x[1]) newfrac = (hi[1]-x[1]) / (xnew[1]-x[1]);
        else newfrac = 0.0;
        if (newfrac < 0.0) newfrac = 0.0;
        else if (newfrac > 1.0) newfrac = 1.0;
        if (newfrac < frac) {
          frac = newfrac;
          outface = YHI;
        }
      }
    }

    if (DIM == 1) {
      if (x[1] == lo[1] && (pflag == PEXIT || v[1] < 0.0)) {
        frac = 0.0;
        outface = YLO;
      } else if (GeometryKokkos::
                 axi_horizontal_line(dtremain,x,v,lo[1],itmp,tc,tmp)) {
        newfrac = tc/dtremain;
        if (newfrac < frac) {
          frac = newfrac;
          outface = YLO;
        }
      }

      if (x[1] == hi[1] && (pflag == PEXIT || v[1] > 0.0)) {
        frac = 0.0;
        outface = YHI;
      } else {
        rnew = sqrt(xnew[1]*xnew[1] + xnew[2]*xnew[2]);
        if (rnew >= hi[1]) {
          if (GeometryKokkos::
              axi_horizontal_line(dtremain,x,v,hi[1],itmp,tc,tmp)) {
            newfrac = tc/dtremain;
            if (newfrac < frac) {
              frac = newfrac;
              outface = YHI;
            }
          }
        }
      }

      pflag = 0;
    }

    if (DIM == 3) {
      if (xnew[2] < lo[2]) {
        if (xnew[2] != x[2]) newfrac = (lo[2]-x[2]) / (xnew[2]-x[2]);
        else newfrac = 0.0;
        if (newfrac < 0.0) newfrac = 0.0;
        else if (newfrac > 1.0) newfrac = 1.0;
        if (newfrac < frac) {
          frac = newfrac;
          outface = ZLO;
        }
      } else if (xnew[2] >= hi[2]) {
        if (xnew[2] != x[2]) newfrac = (hi[2]-x[2]) / (xnew[2]-x[2]);
        else newfrac = 0.0;
        if (newfrac < 0.0) newfrac = 0.0;
        else if (newfrac > 1.0) newfrac = 1.0;
        if (newfrac < frac) {
          frac = newfrac;
          outface = ZHI;
        }
      }
    }

#ifdef MOVE_DEBUG
    if (ntimestep == MOVE_DEBUG_STEP &&
        (MOVE_DEBUG_ID == d_particles[i].id ||
         (me == MOVE_DEBUG_PROC && i == MOVE_DEBUG_INDEX))) {
      if (outface != INTERIOR)
        printf("  OUTFACE %d out: %d %d, frac %g\n",
               outface,grid_kk_copy.obj.neigh_decode(nmask,outface),
               neigh[outface],frac);
      else
        printf("  INTERIOR %d %d\n",outface,INTERIOR);
    }
#endif

    // START of code specific to surfaces

    if (SURF) {

      // skip surf checks if particle flagged as EXITing this cell
      // then unset pflag so not checked again for this particle

      nsurf = d_cells[icell].nsurf;
      if (pflag == PEXIT) {
        nsurf = 0;
        pflag = 0;
      }

      if (ATOMIC_REDUCTION == 1)
        Kokkos::atomic_add(&d_nscheck_one(),nsurf);
      else if (ATOMIC_REDUCTION == 0)
        d_nscheck_one() += nsurf;
      else
        reduce.nscheck_one += nsurf;

      if (nsurf) {

        // particle crosses cell face, reset xnew exactly on face of cell
        // so surface check occurs only for particle path within grid cell
        // xhold = saved xnew so can restore below if no surf collision

        if (outface != INTERIOR) {
          xhold[0] = xnew[0];
          xhold[1] = xnew[1];
          if (DIM != 2) xhold[2] = xnew[2];

          xnew[0] = x[0] + frac*(xnew[0]-x[0]);
          xnew[1] = x[1] + frac*(xnew[1]-x[1]);
          if (DIM != 2) xnew[2] = x[2] + frac*(xnew[2]-x[2]);

          if (outface == XLO) xnew[0] = lo[0];
          else if (outface == XHI) xnew[0] = hi[0];
          else if (outface == YLO) xnew[1] = lo[1];
          else if (outface == YHI) xnew[1] = hi[1];
          else if (outface == ZLO) xnew[2] = lo[2];
          else if (outface == ZHI) xnew[2] = hi[2];
        }

        // for axisymmetric, dtsurf = time that particle stays in cell
        // used as arg to axi_line_intersect()

        if (DIM == 1) {
          if (outface == INTERIOR) dtsurf = dtremain;
          else dtsurf = dtremain * frac;
        }

        // check for collisions with triangles or lines in cell
        // find 1st surface hit via minparam
        // skip collisions with previous surf, but not for axisymmetric
        // not considered collision if 2 params are tied and one INSIDE surf
        // if collision occurs, perform collision with surface model
        // reset x,v,xnew,dtremain and continue single particle trajectory

        cflag = 0;
        minparam = 2.0;
        auto csurfs_begin = d_csurfs.row_map(icell);

        for (int m = 0; m < nsurf; m++) {
          isurf = d_csurfs.entries(csurfs_begin + m);

          if (DIM > 1) {
            if (isurf == exclude) continue;
          }
          if (DIM == 3) {
            tri = &d_tris[isurf];
            hitflag = GeometryKokkos::
              line_tri_intersect(x,xnew,
                                 tri->p1,tri->p2,
                                 tri->p3,tri->norm,xc,param,side);
          }
          if (DIM == 2) {
            line = &d_lines[isurf];
            hitflag = GeometryKokkos::
              line_line_intersect(x,xnew,
                                  line->p1,line->p2,
                                  line->norm,xc,param,side);
          }
          if (DIM == 1) {
            line = &d_lines[isurf];
            hitflag = GeometryKokkos::
              axi_line_intersect(dtsurf,x,v,outface,lo,hi,
                                 line->p1,line->p2,
                                 line->norm,exclude == isurf,
                                 xc,vc,param,side);
          }

#ifdef MOVE_DEBUG
          if (DIM == 3) {
            if (hitflag && ntimestep == MOVE_DEBUG_STEP &&
                (MOVE_DEBUG_ID == d_particles[i].id ||
                 (me == MOVE_DEBUG_PROC && i == MOVE_DEBUG_INDEX)))
              printf("SURF COLLIDE: %d %d %d %d: "
                     "P1 %g %g %g: P2 %g %g %g: "
                     "T1 %g %g %g: T2 %g %g %g: T3 %g %g %g: "
                     "TN %g %g %g: XC %g %g %g: "
                     "Param %g: Side %d\n",
                     MOVE_DEBUG_INDEX,icell,nsurf,isurf,
                     x[0],x[1],x[2],xnew[0],xnew[1],xnew[2],
                     tri->p1[0],tri->p1[1],tri->p1[2],
                     tri->p2[0],tri->p2[1],tri->p2[2],
                     tri->p3[0],tri->p3[1],tri->p3[2],
                     tri->norm[0],tri->norm[1],tri->norm[2],
                     xc[0],xc[1],xc[2],param,side);
          }
          if (DIM == 2) {
            if (hitflag && ntimestep == MOVE_DEBUG_STEP &&
                (MOVE_DEBUG_ID == d_particles[i].id ||
                 (me == MOVE_DEBUG_PROC && i == MOVE_DEBUG_INDEX)))
              printf("SURF COLLIDE: %d %d %d %d: P1 %g %g: P2 %g %g: "
                     "L1 %g %g: L2 %g %g: LN %g %g: XC %g %g: "
                     "Param %g: Side %d\n",
                     MOVE_DEBUG_INDEX,icell,nsurf,isurf,
                     x[0],x[1],xnew[0],xnew[1],
                     line->p1[0],line->p1[1],line->p2[0],line->p2[1],
                     line->norm[0],line->norm[1],
                     xc[0],xc[1],param,side);
          }
          if (DIM == 1) {
            if (hitflag && ntimestep == MOVE_DEBUG_STEP &&
                (MOVE_DEBUG_ID == d_particles[i].id ||
                 (me == MOVE_DEBUG_PROC && i == MOVE_DEBUG_INDEX)))
              printf("SURF COLLIDE %d %ld: %d %d %d %d: P1 %g %g: P2 %g %g: "
                     "L1 %g %g: L2 %g %g: LN %g %g: XC %g %g: "
                     "VC %g %g %g: Param %g: Side %d\n",
                     hitflag,ntimestep,MOVE_DEBUG_INDEX,icell,nsurf,isurf,
                     x[0],x[1],
                     xnew[0],sqrt(xnew[1]*xnew[1]+xnew[2]*xnew[2]),
                     line->p1[0],line->p1[1],line->p2[0],line->p2[1],
                     line->norm[0],line->norm[1],
                     xc[0],xc[1],vc[0],vc[1],vc[2],param,side);
            double edge1[3],edge2[3],xfinal[3],cross[3];
            MathExtraKokkos::sub3(line->p2,line->p1,edge1);
            MathExtraKokkos::sub3(x,line->p1,edge2);
            MathExtraKokkos::cross3(edge2,edge1,cross);
            if (hitflag && ntimestep == MOVE_DEBUG_STEP &&
                MOVE_DEBUG_ID == d_particles[i].id)
              printf("CROSSSTART %g %g %g\n",cross[0],cross[1],cross[2]);
            xfinal[0] = xnew[0];
            xfinal[1] = sqrt(xnew[1]*xnew[1]+xnew[2]*xnew[2]);
            xfinal[2] = 0.0;
            MathExtraKokkos::sub3(xfinal,line->p1,edge2);
            MathExtraKokkos::cross3(edge2,edge1,cross);
            if (hitflag && ntimestep == MOVE_DEBUG_STEP &&
                MOVE_DEBUG_ID == d_particles[i].id)
              printf("CROSSFINAL %g %g %g\n",cross[0],cross[1],cross[2]);
          }
#endif

          if (hitflag && param < minparam && side == OUTSIDE) {

            // NOTE: these were the old checks
            //       think it is now sufficient to test for particle
            //       in an INSIDE cell in fix grid/check

          //if (hitflag && side != ONSURF2OUT && param <= minparam)

            // this if test is to avoid case where particle
            // previously hit 1 of 2 (or more) touching angled surfs at
            // common edge/corner, on this iteration first surf
            // is excluded, but others may be hit on inside:
            // param will be epsilon and exclude must be set
            // skip the hits of other touching surfs

            //if (side == INSIDE && param < EPSPARAM && exclude >= 0)
            // continue;

            // this if test is to avoid case where particle
            // hits 2 touching angled surfs at common edge/corner
            // from far away:
            // param is same, but hits one on outside, one on inside
            // only keep surf hit on outside

            //if (param == minparam && side == INSIDE) continue;

            cflag = 1;
            minparam = param;
            // minside = side;
            minsurf = isurf;
            minxc[0] = xc[0];
            minxc[1] = xc[1];
            if (DIM == 3) minxc[2] = xc[2];
            if (DIM == 1) {
              minvc[1] = vc[1];
              minvc[2] = vc[2];
            }
          }

        } // END of for loop over surfs

        // tri/line = surf that particle hit first

        if (cflag) {
          if (DIM == 3) tri = &d_tris[minsurf];
          if (DIM != 3) line = &d_lines[minsurf];

          // set x to collision point
          // if axisymmetric, set v to remapped velocity at collision pt

          x[0] = minxc[0];
          x[1] = minxc[1];
          if (DIM == 3) x[2] = minxc[2];
          if (DIM == 1) {
            v[1] = minvc[1];
            v[2] = minvc[2];
          }

          // perform surface collision using surface collision model
          // surface chemistry may destroy particle or create new one
          // must update particle's icell to current icell so that
          //   if jpart is created, it will be added to correct cell
          // if jpart, add new particle to this iteration via pstop++
          // tally surface collision stats if requested using iorig

          ipart = &particle_i;
          ipart->icell = icell;
          dtremain *= 1.0 - minparam*frac;

          if (nsurf_tally)
            iorig = particle_i;
          const int n = DIM == 3 ? tri->isc : line->isc;

          if (DIM == 3) {
            jpart = surf_collide_dispatch<REACT,ATOMIC_REDUCTION>
              (n,ipart,dtremain,minsurf,tri->norm,tri->isr,reaction,d_retry,d_nlocal);
          }

          if (DIM != 3) {
            jpart = surf_collide_dispatch<REACT,ATOMIC_REDUCTION>
              (n,ipart,dtremain,minsurf,line->norm,line->isr,reaction,d_retry,d_nlocal);
          }

          if (jpart) {
            x = particle_i.x;
            v = particle_i.v;
            jpart->flag = PSURF + 1 + minsurf;
            jpart->dtremain = dtremain;
            jpart->weight = particle_i.weight;
          }

          if (nsurf_tally) {
            for (int m = 0; m < nslist_surf; m++)
              UK_SLIST_SURF(m).
                    surf_tally_kk<ATOMIC_REDUCTION>(dtremain,minsurf,icell,reaction,&iorig,ipart,jpart);
            for (int m = 0; m < nslist_isurf; m++)
              UK_SLIST_ISURF(m).
                    surf_tally_kk<ATOMIC_REDUCTION>(dtremain,minsurf,icell,reaction,&iorig,ipart,jpart);
            for (int m = 0; m < nslist_coll_tally; m++)
              UK_SLIST_COLL_TALLY(m).
                    surf_tally_kk(dtremain,minsurf,icell,reaction,&iorig,ipart,jpart);
            for (int m = 0; m < nslist_react_tally; m++)
              UK_SLIST_REACT_TALLY(m).
                    surf_tally_kk(dtremain,minsurf,icell,reaction,&iorig,ipart,jpart);
            for (int m = 0; m < nslist_react_isurf; m++)
              UK_SLIST_REACT_ISURF(m).
                    surf_tally_kk<ATOMIC_REDUCTION>(dtremain,minsurf,icell,reaction,&iorig,ipart,jpart);
            for (int m = 0; m < nslist_react_surf; m++)
              UK_SLIST_REACT_SURF(m).
                    surf_tally_kk<ATOMIC_REDUCTION>(dtremain,minsurf,icell,reaction,&iorig,ipart,jpart);
          }

          // stuck_iterate = consecutive iterations particle is immobile

          if (minparam <= 1.0e-14) stuck_iterate++;
          else stuck_iterate = 0;

          // reset post-bounce xnew

          xnew[0] = x[0] + dtremain*v[0];
          xnew[1] = x[1] + dtremain*v[1];
          if (DIM != 2) xnew[2] = x[2] + dtremain*v[2];

          exclude = minsurf;
          if (ATOMIC_REDUCTION == 1)
            Kokkos::atomic_inc(&d_nscollide_one());
          else if (ATOMIC_REDUCTION == 0)
            d_nscollide_one()++;
          else
            reduce.nscollide_one++;

#ifdef MOVE_DEBUG
          if (DIM == 3) {
            if (ntimestep == MOVE_DEBUG_STEP &&
                (MOVE_DEBUG_ID == d_particles[i].id ||
                 (me == MOVE_DEBUG_PROC && i == MOVE_DEBUG_INDEX)))
              printf("POST COLLISION %d: %g %g %g: %g %g %g: %g %g %g\n",
                     MOVE_DEBUG_INDEX,
                     x[0],x[1],x[2],xnew[0],xnew[1],xnew[2],
                     minparam,frac,dtremain);
          }
          if (DIM == 2) {
            if (ntimestep == MOVE_DEBUG_STEP &&
                (MOVE_DEBUG_ID == d_particles[i].id ||
                 (me == MOVE_DEBUG_PROC && i == MOVE_DEBUG_INDEX)))
              printf("POST COLLISION %d: %g %g: %g %g: %g %g %g\n",
                     MOVE_DEBUG_INDEX,
                     x[0],x[1],xnew[0],xnew[1],
                     minparam,frac,dtremain);
          }
          if (DIM == 1) {
            if (ntimestep == MOVE_DEBUG_STEP &&
                (MOVE_DEBUG_ID == d_particles[i].id ||
                 (me == MOVE_DEBUG_PROC && i == MOVE_DEBUG_INDEX)))
              printf("POST COLLISION %d: %g %g: %g %g: vel %g %g %g: %g %g %g\n",
                     MOVE_DEBUG_INDEX,
                     x[0],x[1],
                     xnew[0],sqrt(xnew[1]*xnew[1]+xnew[2]*xnew[2]),
                     v[0],v[1],v[2],
                     minparam,frac,dtremain);
          }
#endif

          // if ipart = NULL, particle discarded due to surface chem
          // else if particle not stuck, continue advection while loop
          // if stuck, mark for DISCARD, and drop out of SURF code

          if (ipart == NULL) particle_i.flag = PDISCARD;
          else if (stuck_iterate < MAXSTUCK) continue;
          else {
            particle_i.flag = PDISCARD;
            if (ATOMIC_REDUCTION == 1)
              Kokkos::atomic_inc(&d_nstuck());
            else if (ATOMIC_REDUCTION == 0)
              d_nstuck()++;
            else
              reduce.nstuck++;
          }

        } // END of cflag if section that performed collision

        // no collision, so restore saved xnew if changed it above

        if (outface != INTERIOR) {
          xnew[0] = xhold[0];
          xnew[1] = xhold[1];
          if (DIM != 2) xnew[2] = xhold[2];
        }

      } // END of if test for any surfs in this cell
    } // END of code specific to surfaces

    // break from advection loop if discarding particle

    if (particle_i.flag == PDISCARD) break;

    // no cell crossing
    // set final particle position to xnew, then break from advection loop
    // for axisymmetry, must first remap linear xnew and v
    // for axisymmetry, check if final particle position is within cell
    //   can be rare epsilon round-off cases where particle ends up outside
    //     of final cell curved surf when move logic thinks it is inside
    //   example is when Geom::axi_horizontal_line() says no crossing of cell edge
    //     but axi_remap() puts particle outside the cell
    //   in this case, just DISCARD particle and tally it to naxibad
    // if migrating to another proc,
    //   flag as PDONE so new proc won't move it more on this step

    if (outface == INTERIOR) {
      if (DIM == 1) axi_remap(xnew,v);
      x[0] = xnew[0];
      x[1] = xnew[1];
      if (DIM == 3) x[2] = xnew[2];
      if (DIM == 1) {
        if (x[1] < lo[1] || x[1] > hi[1]) {
          particle_i.flag = PDISCARD;
          if (ATOMIC_REDUCTION == 1)
            Kokkos::atomic_inc(&d_naxibad());
          else if (ATOMIC_REDUCTION == 0)
            d_naxibad()++;
          else
            reduce.naxibad++;
          break;
        }
      }
      if (d_cells[icell].proc != me) particle_i.flag = PDONE;
      break;
    }

    // particle crosses cell face
    // decrement dtremain in case particle is passed to another proc
    // for axisymmetry, must then remap linear x and v
    // reset particle x to be exactly on cell face
    // for axisymmetry, must reset xnew for next iteration since v changed

    dtremain *= 1.0-frac;
    exclude = -1;

    x[0] += frac * (xnew[0]-x[0]);
    x[1] += frac * (xnew[1]-x[1]);
    if (DIM != 2) x[2] += frac * (xnew[2]-x[2]);
    if (DIM == 1) axi_remap(x,v);

    if (outface == XLO) x[0] = lo[0];
    else if (outface == XHI) x[0] = hi[0];
    else if (outface == YLO) x[1] = lo[1];
    else if (outface == YHI) x[1] = hi[1];
    else if (outface == ZLO) x[2] = lo[2];
    else if (outface == ZHI) x[2] = hi[2];

    if (DIM == 1) {
      xnew[0] = x[0] + dtremain*v[0];
      xnew[1] = x[1] + dtremain*v[1];
      xnew[2] = x[2] + dtremain*v[2];
    }

    // nflag = type of neighbor cell: child, parent, unknown, boundary
    // if parent, use id_find_child to identify child cell
    //   result can be -1 for unknown cell, occurs when:
    //   (a) particle hits face of ghost child cell
    //   (b) the ghost cell extends beyond ghost halo
    //   (c) cell on other side of face is a parent
    //   (d) its child, which the particle is in, is entirely beyond my halo
    // if new cell is child and surfs exist, check if a split cell

    nflag = grid_kk_copy.obj.neigh_decode(nmask,outface);
    icell_original = icell;

    if (nflag == NCHILD) {
      icell = neigh[outface];
      if (DIM == 3 && SURF) {
        if (d_cells[icell].nsplit > 1 && d_cells[icell].nsurf >= 0)
          icell = split3d(icell,x);
      }
      if (DIM < 3 && SURF) {
        if (d_cells[icell].nsplit > 1 && d_cells[icell].nsurf >= 0)
          icell = split2d(icell,x);
      }
    } else if (nflag == NPARENT) {
      auto pcell = &d_pcells[neigh[outface]];
      icell = grid_kk_copy.obj.id_find_child(pcell->id,d_cells[icell].level,
                                             pcell->lo,pcell->hi,x);
      if (icell >= 0) {
        if (DIM == 3 && SURF) {
          if (d_cells[icell].nsplit > 1 && d_cells[icell].nsurf >= 0)
            icell = split3d(icell,x);
        }
        if (DIM < 3 && SURF) {
          if (d_cells[icell].nsplit > 1 && d_cells[icell].nsurf >= 0)
            icell = split2d(icell,x);
        }
      }
    } else if (nflag == NUNKNOWN) icell = -1;

    // neighbor cell is global boundary
    // tally boundary stats if requested using iorig
    // collide() updates x,v,xnew as needed due to boundary interaction
    //   may also update dtremain (piston BC)
    // for axisymmetric, must recalculate xnew since v may have changed
    // surface chemistry may destroy particle or create new one
    // if jpart, add new particle to this iteration via pstop++
    // OUTFLOW: exit with particle flag = PDISCARD
    // PERIODIC: new cell via same logic as above for child/parent/unknown
    // OTHER = reflected particle stays in same grid cell

    else {
      ipart = &particle_i;

      Particle::OnePart iorig;
      if (nboundary_tally)
        memcpy(&iorig,&particle_i,sizeof(Particle::OnePart));

      // from Domain:

      Particle::OnePart* ipart = &particle_i;
      lo = d_cells[icell].lo;
      hi = d_cells[icell].hi;
      if (domain_kk_copy.obj.bflag[outface] == SURFACE) {
        // treat global boundary as a surface
        // particle velocity is changed by surface collision model
        // dtremain may be changed by collision model
        // reset all components of xnew, in case dtremain changed
        // if axisymmetric, caller will reset again, including xnew[2]

        const int n = domain_kk_copy.obj.surf_collide[outface];

        jpart = surf_collide_dispatch<REACT,ATOMIC_REDUCTION>
          (n,ipart,dtremain,-(outface+1),domain_kk_copy.obj.norm[outface],
           domain_kk_copy.obj.surf_react[outface],reaction,d_retry,d_nlocal);

        if (ipart) {
          double *x = ipart->x;
          double *v = ipart->v;
          xnew[0] = x[0] + dtremain*v[0];
          xnew[1] = x[1] + dtremain*v[1];
          if (domain_kk_copy.obj.dimension == 3) xnew[2] = x[2] + dtremain*v[2];
        }
        bflag = SURFACE;
      } else {
        bflag = domain_kk_copy.obj.collide_kokkos(ipart,outface,lo,hi,xnew/*,dtremain*/,reaction);
      }

      if (jpart) {
        x = particle_i.x;
        v = particle_i.v;
      }

      if (nboundary_tally) {
        for (int m = 0; m < nblist_boundary; m++)
          UK_BLIST(m).
            boundary_tally_kk<ATOMIC_REDUCTION>(dtremain,outface,bflag,reaction,&iorig,ipart,jpart,domain_kk_copy.obj.norm[outface]);
        for (int m = 0; m < nblist_react; m++)
          UK_BLIST_REACT(m).
            boundary_tally_kk<ATOMIC_REDUCTION>(dtremain,outface,bflag,reaction,&iorig,ipart,jpart,domain_kk_copy.obj.norm[outface]);
      }

      if (DIM == 1) {
        xnew[0] = x[0] + dtremain*v[0];
        xnew[1] = x[1] + dtremain*v[1];
        xnew[2] = x[2] + dtremain*v[2];
      }

      if (bflag == OUTFLOW) {
        particle_i.flag = PDISCARD;
        if (ATOMIC_REDUCTION == 1)
          Kokkos::atomic_inc(&d_nexit_one());
        else if (ATOMIC_REDUCTION == 0)
          d_nexit_one()++;
        else
          reduce.nexit_one++;
        break;
      } else if (bflag == PERIODIC) {
        if (nflag == NPBCHILD) {
          icell = neigh[outface];
          if (DIM == 3 && SURF) {
            if (d_cells[icell].nsplit > 1 && d_cells[icell].nsurf >= 0)
              icell = split3d(icell,x);
          }
          if (DIM < 3 && SURF) {
            if (d_cells[icell].nsplit > 1 && d_cells[icell].nsurf >= 0)
              icell = split2d(icell,x);
          }
        } else if (nflag == NPBPARENT) {
          auto pcell = &d_pcells[neigh[outface]];
          icell = grid_kk_copy.obj.id_find_child(pcell->id,d_cells[icell].level,
                                                 pcell->lo,pcell->hi,x);
          if (icell >= 0) {
            if (DIM == 3 && SURF) {
              if (d_cells[icell].nsplit > 1 && d_cells[icell].nsurf >= 0)
                icell = split3d(icell,x);
            }
            if (DIM < 3 && SURF) {
              if (d_cells[icell].nsplit > 1 && d_cells[icell].nsurf >= 0)
                icell = split2d(icell,x);
            }
          } else domain_kk_copy.obj.uncollide_kokkos(outface,x);
        } else if (nflag == NPBUNKNOWN) {
          icell = -1;
          domain_kk_copy.obj.uncollide_kokkos(outface,x);
        }

      } else if (bflag == SURFACE) {
        if (ipart == NULL) {
          particle_i.flag = PDISCARD;
          break;
        } else if (jpart) {
          jpart->flag = PSURF;
          jpart->dtremain = dtremain;
          jpart->weight = particle_i.weight;
        }

        if (ATOMIC_REDUCTION == 1) {
          Kokkos::atomic_inc(&d_nboundary_one());
          Kokkos::atomic_dec(&d_ntouch_one());    // decrement here since will increment below
        } else if (ATOMIC_REDUCTION == 0) {
          d_nboundary_one()++;
          d_ntouch_one()--;    // decrement here since will increment below
        } else {
          reduce.nboundary_one++;
          reduce.ntouch_one--;    // decrement here since will increment below
        }

      } else {
        if (ATOMIC_REDUCTION == 1) {
          Kokkos::atomic_inc(&d_nboundary_one());
          Kokkos::atomic_dec(&d_ntouch_one());    // decrement here since will increment below
        } else if (ATOMIC_REDUCTION == 0) {
          d_nboundary_one()++;
          d_ntouch_one()--;    // decrement here since will increment below
        } else {
          reduce.nboundary_one++;
          reduce.ntouch_one--;    // decrement here since will increment below
        }
      }
    }

    // neighbor cell is unknown
    // reset icell to original icell which must be a ghost cell
    // exit with particle flag = PEXIT, so receiver can identify neighbor

    if (icell < 0) {
      icell = icell_original;
      particle_i.flag = PEXIT;
      particle_i.dtremain = dtremain;
      d_entryexit() = 1;
      break;
    }

    // if nsurf < 0, new cell is EMPTY ghost
    // exit with particle flag = PENTRY, so receiver can continue move

    if (d_cells[icell].nsurf < 0) {
      particle_i.flag = PENTRY;
      particle_i.dtremain = dtremain;
      d_entryexit() = 1;
      break;
    }

    // move particle into new grid cell for next stage of move

    lo = d_cells[icell].lo;
    hi = d_cells[icell].hi;
    neigh = d_cells[icell].neigh;
    nmask = d_cells[icell].nmask;
    if (ATOMIC_REDUCTION == 1)
      Kokkos::atomic_inc(&d_ntouch_one());
    else if (ATOMIC_REDUCTION == 0)
      d_ntouch_one()++;
    else
      reduce.ntouch_one++;
  }

  // END of while loop over advection of single particle

#ifdef MOVE_DEBUG
  if (ntimestep == MOVE_DEBUG_STEP &&
      (MOVE_DEBUG_ID == d_particles[i].id ||
       (me == MOVE_DEBUG_PROC && i == MOVE_DEBUG_INDEX)))
    printf("MOVE DONE %d %d %d: %g %g %g: DTR %g\n",
           MOVE_DEBUG_INDEX,d_particles[i].flag,icell,
           x[0],x[1],x[2],dtremain);
#endif

  // move is complete, or as much as can be done on this proc
  // update particle's grid cell
  // if particle flag set, add particle to migrate list
  // if discarding, migration will delete particle

  particle_i.icell = icell;

  if (particle_i.flag != PKEEP) {
    int index;
    if (ATOMIC_REDUCTION == 0) {
      index = d_nmigrate();
      d_nmigrate()++;
    } else {
      index = Kokkos::atomic_fetch_add(&d_nmigrate(),1);
    }
    k_mlist.view_device()[index] = i;
    if (particle_i.flag != PDISCARD) {
      if (d_cells[icell].proc == me && !d_error_flag()) {
        d_error_flag() = 1;
        return;
      }
      if (ATOMIC_REDUCTION == 1)
        Kokkos::atomic_inc(&d_ncomm_one());
      else if (ATOMIC_REDUCTION == 0)
        d_ncomm_one()++;
      else
        reduce.ncomm_one++;
    }
  }
} // end of Kokkos parallel_reduce

/* ----------------------------------------------------------------------
   particle is entering split parent icell at x
   determine which split child cell it is in
   return index of sub-cell in ChildCell
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
int UpdateKokkos::split3d(int icell, double *x) const
{
  int m,cflag,isurf,hitflag,side,minsurfindex;
  double param,minparam;
  double xc[3];
  Surf::Tri *tri;

  // check for collisions with lines in cell
  // find 1st surface hit via minparam
  // only consider tris that are mapped via csplits to a split cell
  //   unmapped tris only touch cell surf at xnew
  //   another mapped tri should include same xnew
  // NOTE: these next 2 lines do not seem correct compared to code
  // not considered a collision if particles starts on surf, moving out
  // not considered a collision if 2 params are tied and one is INSIDE surf

  int nsurf = d_cells[icell].nsurf;
  int isplit = d_cells[icell].isplit;
  double *xnew = d_sinfo[isplit].xsplit;

  cflag = 0;
  minparam = 2.0;

  auto csplits_begin = d_csplits.row_map(isplit);
  auto csurfs_begin = d_csurfs.row_map(icell);
  for (m = 0; m < nsurf; m++) {
    if (d_csplits.entries(csplits_begin + m) < 0) continue;
    isurf = d_csurfs.entries(csurfs_begin + m);
    tri = &d_tris[isurf];
    hitflag = GeometryKokkos::
      line_tri_intersect(x,xnew,
                         tri->p1,tri->p2,tri->p3,
                         tri->norm,xc,param,side);

    if (hitflag && side != INSIDE && param < minparam) {
      cflag = 1;
      minparam = param;
      minsurfindex = m;
    }
  }

  auto csubs_begin = d_csubs.row_map(isplit);
  if (!cflag) return d_csubs.entries(csubs_begin + d_sinfo[isplit].xsub);
  int index = d_csplits.entries(csplits_begin + minsurfindex);
  return d_csubs.entries(csubs_begin + index);
}

/* ----------------------------------------------------------------------
   particle is entering split ICELL at X
   determine which split sub-cell it is in
   return index of sub-cell in ChildCell
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
int UpdateKokkos::split2d(int icell, double *x) const
{
  int m,cflag,isurf,hitflag,side,minsurfindex;
  double param,minparam;
  double xc[3];
  Surf::Line *line;

  // check for collisions with lines in cell
  // find 1st surface hit via minparam
  // only consider lines that are mapped via csplits to a split cell
  //   unmapped lines only touch cell surf at xnew
  //   another mapped line should include same xnew
  // NOTE: these next 2 lines do not seem correct compared to code
  // not considered a collision if particle starts on surf, moving out
  // not considered a collision if 2 params are tied and one is INSIDE surf

  int nsurf = d_cells[icell].nsurf;
  int isplit = d_cells[icell].isplit;
  double *xnew = d_sinfo[isplit].xsplit;

  cflag = 0;
  minparam = 2.0;
  auto csplits_begin = d_csplits.row_map(isplit);
  auto csurfs_begin = d_csurfs.row_map(icell);
  for (m = 0; m < nsurf; m++) {
    if (d_csplits.entries(csplits_begin + m) < 0) continue;
    isurf = d_csurfs.entries(csurfs_begin + m);
    line = &d_lines[isurf];
    hitflag = GeometryKokkos::
      line_line_intersect(x,xnew,
                          line->p1,line->p2,line->norm,
                          xc,param,side);

    if (hitflag && side != INSIDE && param < minparam) {
      cflag = 1;
      minparam = param;
      minsurfindex = m;
    }
  }

  auto csubs_begin = d_csubs.row_map(isplit);
  if (!cflag) return d_csubs.entries(csubs_begin + d_sinfo[isplit].xsub);
  int index = d_csplits.entries(csplits_begin + minsurfindex);
  return d_csubs.entries(csubs_begin + index);
}

/* ----------------------------------------------------------------------
   set bounce tally flags for current timestep
   nsurf_tally = # of computes needing bounce info on this step
   clear accumulators in computes that will be invoked this step
------------------------------------------------------------------------- */

void UpdateKokkos::tally_set(bigint ntimestep)
{
  Update::tally_set(ntimestep);

  int i;

  // dispatch by dynamic_cast, as setup_surf_tally_copies() does: compute
  //   boundary and compute react/boundary both set boundary_tally_flag but
  //   are unrelated class hierarchies, so a static cast would call one's
  //   methods on the other.  The cast also fails for a plain compute boundary
  //   under "-k on" without "-sf kk", which is likewise not the Kokkos class

  // count first: the buffers have to be sized before anything is blitted in

  nblist_boundary = nblist_react = 0;
  for (i = 0; i < nboundary_tally; i++) {
    if (dynamic_cast<ComputeBoundaryKokkos*>(blist_active[i])) nblist_boundary++;
    else if (dynamic_cast<ComputeReactBoundaryKokkos*>(blist_active[i])) nblist_react++;
    else
      error->all(FLERR,"Kokkos does not (yet) support this boundary tally compute; "
                       "use a Kokkos-enabled boundary tally compute (-sf kk)");
  }

#ifdef SPARTA_KOKKOS_FIXED_LISTS
  if (nblist_boundary > KOKKOS_MAX_BLIST || nblist_react > KOKKOS_MAX_BLIST)
    error->all(FLERR,"Kokkos currently only supports two instances of compute boundary");
#else
  tally_buf_resize<ComputeBoundaryKokkos>(k_blist,d_blist,nblist_boundary);
  tally_buf_resize<ComputeReactBoundaryKokkos>(k_blist_react,d_blist_react,nblist_react);
#endif

  nblist_boundary = nblist_react = 0;

  for (i = 0; i < nboundary_tally; i++) {
    if (ComputeBoundaryKokkos* c =
          dynamic_cast<ComputeBoundaryKokkos*>(blist_active[i])) {
      c->pre_boundary_tally();
#ifdef SPARTA_KOKKOS_FIXED_LISTS
      blist_active_copy[nblist_boundary].copy(c);
#else
      tally_buf_blit(k_blist,nblist_boundary,c);
#endif
      nblist_boundary++;
    } else if (ComputeReactBoundaryKokkos* c =
                 dynamic_cast<ComputeReactBoundaryKokkos*>(blist_active[i])) {
      c->pre_boundary_tally();
#ifdef SPARTA_KOKKOS_FIXED_LISTS
      blist_active_react_copy[nblist_react].copy(c);
#else
      tally_buf_blit(k_blist_react,nblist_react,c);
#endif
      nblist_react++;
    } else
      error->all(FLERR,"Kokkos does not (yet) support this boundary tally compute; "
                       "use a Kokkos-enabled boundary tally compute (-sf kk)");
  }

#ifdef SPARTA_KOKKOS_FIXED_LISTS
  for (i = nblist_boundary; i < KOKKOS_MAX_BLIST; i++)
    blist_active_copy[i].copy(&tmp_compute_boundary_kk);
  for (i = nblist_react; i < KOKKOS_MAX_BLIST; i++)
    blist_active_react_copy[i].copy(&tmp_compute_react_boundary_kk);
#else
  tally_buf_sync(k_blist,d_blist);
  tally_buf_sync(k_blist_react,d_blist_react);
#endif

  // surf-tally compute scatter views (slist_active_copy et al.) are
  //   (re)established in setup_surf_tally_copies(), which run() calls after
  //   start-of-step fixes have executed.  Doing it here would be unsafe: a
  //   start-of-step fix such as fix emit/surf runs its own surf-tally session
  //   that reallocates each compute's dup_array_surf_tally scatter view,
  //   freeing the one captured here and leaving the move kernel's functor
  //   copy pointing at freed memory.  See UpdateKokkos::run().
}

/* ----------------------------------------------------------------------
   set up per-compute surf-tally copies used on-device by the move kernel
   must be called after start-of-step fixes run (see tally_set/run)
------------------------------------------------------------------------- */

void UpdateKokkos::setup_surf_tally_copies()
{
  // partition the active surf tally computes by type, one runtime-sized
  //   device buffer each; all of them tally on-device via surf_tally_kk(),
  //   invoked from the move kernel's surface collision loop
  // dispatch by dynamic_cast, not by style string: the styles are also
  //   registered under explicit "/kk" names (e.g. isurf/grid/kk), so a
  //   style-string compare would reject a compute the user typed with the
  //   suffix.  The Kokkos tally computes are unrelated class hierarchies,
  //   so the casts are mutually exclusive and order-independent.

  nslist_surf = nslist_isurf = nslist_react_isurf = nslist_react_surf = 0;
  nslist_coll_tally = nslist_react_tally = 0;

  // count first: the buffers have to be sized before anything is blitted in

  for (int i = 0; i < nsurf_tally; i++) {
    if (dynamic_cast<ComputeISurfGridKokkos*>(slist_active[i])) nslist_isurf++;
    else if (dynamic_cast<ComputeReactISurfGridKokkos*>(slist_active[i])) nslist_react_isurf++;
    else if (dynamic_cast<ComputeReactSurfKokkos*>(slist_active[i])) nslist_react_surf++;
    else if (dynamic_cast<ComputeSurfKokkos*>(slist_active[i])) nslist_surf++;
    else if (dynamic_cast<ComputeSurfCollisionTallyKokkos*>(slist_active[i])) nslist_coll_tally++;
    else if (dynamic_cast<ComputeSurfReactionTallyKokkos*>(slist_active[i])) nslist_react_tally++;
    else
      error->all(FLERR,"Kokkos does not (yet) support this surf tally compute; "
                       "use a Kokkos-enabled surf tally compute (-sf kk)");
  }

#ifdef SPARTA_KOKKOS_FIXED_LISTS
  if (nslist_isurf > KOKKOS_MAX_SLIST || nslist_react_isurf > KOKKOS_MAX_SLIST ||
      nslist_react_surf > KOKKOS_MAX_SLIST || nslist_surf > KOKKOS_MAX_SLIST ||
      nslist_coll_tally > KOKKOS_MAX_SLIST || nslist_react_tally > KOKKOS_MAX_SLIST)
    error->all(FLERR,"Kokkos currently only supports two instances of each surf tally compute");
#else
  tally_buf_resize<ComputeISurfGridKokkos>(k_slist_isurf,d_slist_isurf,nslist_isurf);
  tally_buf_resize<ComputeReactISurfGridKokkos>(k_slist_react_isurf,d_slist_react_isurf,nslist_react_isurf);
  tally_buf_resize<ComputeReactSurfKokkos>(k_slist_react_surf,d_slist_react_surf,nslist_react_surf);
  tally_buf_resize<ComputeSurfKokkos>(k_slist_surf,d_slist_surf,nslist_surf);
  tally_buf_resize<ComputeSurfCollisionTallyKokkos>(k_slist_coll_tally,d_slist_coll_tally,nslist_coll_tally);
  tally_buf_resize<ComputeSurfReactionTallyKokkos>(k_slist_react_tally,d_slist_react_tally,nslist_react_tally);
#endif

  // then run each compute's pre_surf_tally() in list order, as before, and
  //   blit it into its type's buffer

  int nisurf = 0, nrisurf = 0, nrsurf = 0, nsurf = 0, nct = 0, nrt = 0;

  for (int i = 0; i < nsurf_tally; i++) {
    if (ComputeISurfGridKokkos* c =
          dynamic_cast<ComputeISurfGridKokkos*>(slist_active[i])) {
      c->pre_surf_tally();
#ifdef SPARTA_KOKKOS_FIXED_LISTS
      slist_active_isurf_copy[nisurf++].copy(c);
#else
      tally_buf_blit(k_slist_isurf,nisurf++,c);
#endif
    } else if (ComputeReactISurfGridKokkos* c =
                 dynamic_cast<ComputeReactISurfGridKokkos*>(slist_active[i])) {
      c->pre_surf_tally();
#ifdef SPARTA_KOKKOS_FIXED_LISTS
      slist_active_react_isurf_copy[nrisurf++].copy(c);
#else
      tally_buf_blit(k_slist_react_isurf,nrisurf++,c);
#endif
    } else if (ComputeReactSurfKokkos* c =
                 dynamic_cast<ComputeReactSurfKokkos*>(slist_active[i])) {
      c->pre_surf_tally();
#ifdef SPARTA_KOKKOS_FIXED_LISTS
      slist_active_react_surf_copy[nrsurf++].copy(c);
#else
      tally_buf_blit(k_slist_react_surf,nrsurf++,c);
#endif
    } else if (ComputeSurfKokkos* c =
                 dynamic_cast<ComputeSurfKokkos*>(slist_active[i])) {
      c->pre_surf_tally();
#ifdef SPARTA_KOKKOS_FIXED_LISTS
      slist_active_copy[nsurf++].copy(c);
#else
      tally_buf_blit(k_slist_surf,nsurf++,c);
#endif
    } else if (ComputeSurfCollisionTallyKokkos* c =
                 dynamic_cast<ComputeSurfCollisionTallyKokkos*>(slist_active[i])) {
      c->pre_surf_tally();
      c->d_overflow = d_tally_overflow;
#ifdef SPARTA_KOKKOS_FIXED_LISTS
      slist_active_coll_tally_copy[nct++].copy(c);
#else
      tally_buf_blit(k_slist_coll_tally,nct++,c);
#endif
    } else if (ComputeSurfReactionTallyKokkos* c =
                 dynamic_cast<ComputeSurfReactionTallyKokkos*>(slist_active[i])) {
      c->pre_surf_tally();
      c->d_overflow = d_tally_overflow;
#ifdef SPARTA_KOKKOS_FIXED_LISTS
      slist_active_react_tally_copy[nrt++].copy(c);
#else
      tally_buf_blit(k_slist_react_tally,nrt++,c);
#endif
    }
  }

#ifdef SPARTA_KOKKOS_FIXED_LISTS
  for (int i = nsurf; i < KOKKOS_MAX_SLIST; i++) slist_active_copy[i].copy(&tmp_compute_surf_kk);
  for (int i = nisurf; i < KOKKOS_MAX_SLIST; i++) slist_active_isurf_copy[i].copy(&tmp_compute_isurf_grid_kk);
  for (int i = nrisurf; i < KOKKOS_MAX_SLIST; i++) slist_active_react_isurf_copy[i].copy(&tmp_compute_react_isurf_grid_kk);
  for (int i = nrsurf; i < KOKKOS_MAX_SLIST; i++) slist_active_react_surf_copy[i].copy(&tmp_compute_react_surf_kk);
#else
  tally_buf_sync(k_slist_isurf,d_slist_isurf);
  tally_buf_sync(k_slist_react_isurf,d_slist_react_isurf);
  tally_buf_sync(k_slist_react_surf,d_slist_react_surf);
  tally_buf_sync(k_slist_surf,d_slist_surf);
  tally_buf_sync(k_slist_coll_tally,d_slist_coll_tally);
  tally_buf_sync(k_slist_react_tally,d_slist_react_tally);
#endif

  // gas/gas tally computes are validated and set up by CollideVSSKokkos,
  //   which invokes their on-device gas_tally_kk() from the collision kernel
}

/* ---------------------------------------------------------------------- */

void UpdateKokkos::backup()
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  d_particles = particle_kk->k_particles.view_device();

  // reuse the buffer across the migration iterations of a step.  backup() is
  //   called once per iteration of the retry loop, and the loop runs once per
  //   migration iteration, so reallocating here churned a full particle-sized
  //   allocation several times per timestep whenever a per-event surf tally
  //   compute was active.  restore() no longer frees it; free_particle_backup()
  //   does, once the step's migration is done, so peak memory is unchanged.
  //   The extents must stay equal because restore() deep_copies between them.

  if (d_particles_backup.extent(0) != d_particles.extent(0))
    d_particles_backup = decltype(d_particles)(Kokkos::view_alloc("update:particles_backup",Kokkos::WithoutInitializing),d_particles.extent(0));

  Kokkos::deep_copy(d_particles_backup,d_particles);

  for (int n = 0; n < surf->nsc; n++) sc_phase(surf->sc[n],SC_BACKUP);
  upload_surf_collide_models();
}

/* ---------------------------------------------------------------------- */

void UpdateKokkos::restore()
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  Kokkos::deep_copy(particle_kk->k_particles.view_device(),d_particles_backup);
  d_particles = particle_kk->k_particles.view_device();

  for (int n = 0; n < surf->nsc; n++) sc_phase(surf->sc[n],SC_RESTORE);
  upload_surf_collide_models();

  // the buffer stays allocated for the next attempt of this step;
  //   free_particle_backup() releases it once the step is done
}

/* ----------------------------------------------------------------------
   release the particle backup buffer at the end of a step's migration
   keeps peak memory the same as when restore() freed it, without
     reallocating on every migration iteration
------------------------------------------------------------------------- */

void UpdateKokkos::free_particle_backup()
{
  d_particles_backup = {};
}

/* ----------------------------------------------------------------------
   grow every per-event surf tally compute past what the failed attempt
     needed, then let the caller repeat the move
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   surf_collide model plumbing
   the nine styles share no Kokkos base class -- pre_collide(), post_collide(),
     backup() and restore() are declared on the concrete classes, not on
     SurfCollide -- so every host-side pass over the models has to name all
     nine types.  It used to be spelled out four times as a strcmp ladder;
     the list lives here once instead
------------------------------------------------------------------------- */

#define SC_FOREACH(F)                                \
  F(SC_SPECULAR,SurfCollideSpecularKokkos)           \
  F(SC_DIFFUSE,SurfCollideDiffuseKokkos)             \
  F(SC_VANISH,SurfCollideVanishKokkos)               \
  F(SC_PISTON,SurfCollidePistonKokkos)               \
  F(SC_TRANSPARENT,SurfCollideTransparentKokkos)     \
  F(SC_ADIABATIC,SurfCollideAdiabaticKokkos)         \
  F(SC_IMPULSIVE,SurfCollideImpulsiveKokkos)         \
  F(SC_TD,SurfCollideTDKokkos)                       \
  F(SC_CLL,SurfCollideCLLKokkos)

namespace {

  template<class T> void sc_run(SurfCollide *base, int phase)
  {
    T *m = (T *) base;
    if (phase == SC_PRE) m->pre_collide();
    else if (phase == SC_POST) m->post_collide();
    else if (phase == SC_BACKUP) m->backup();
    else m->restore();
  }

  template<class T> void sc_blit(char *dst, SurfCollide *base)
  {
    memcpy((void*) dst, (const void*) ((T *) base), sizeof(T));

    // the image in the buffer is read on device and never destructed, so
    //   mark it non-owning exactly as KKCopy::copy() does

    ((T *) dst)->copy = 1;
  }
}

/* ---------------------------------------------------------------------- */

int UpdateKokkos::surf_collide_style_tag(SurfCollide *sc)
{
  if (strcmp(sc->style,"specular") == 0) return SC_SPECULAR;
  if (strcmp(sc->style,"diffuse") == 0) return SC_DIFFUSE;
  if (strcmp(sc->style,"vanish") == 0) return SC_VANISH;
  if (strcmp(sc->style,"piston") == 0) return SC_PISTON;
  if (strcmp(sc->style,"transparent") == 0) return SC_TRANSPARENT;
  if (strcmp(sc->style,"adiabatic") == 0) return SC_ADIABATIC;
  if (strcmp(sc->style,"impulsive") == 0) return SC_IMPULSIVE;
  if (strcmp(sc->style,"td") == 0) return SC_TD;
  if (strcmp(sc->style,"cll") == 0) return SC_CLL;
  return -1;
}

/* ---------------------------------------------------------------------- */

void UpdateKokkos::sc_phase(SurfCollide *sc, int phase)
{
  switch (surf_collide_style_tag(sc)) {
#define SC_RUN(TAG,TYPE) case TAG: sc_run<TYPE>(sc,phase); break;
    SC_FOREACH(SC_RUN)
#undef SC_RUN
  }
}

/* ----------------------------------------------------------------------
   size of one blitted model of each style
------------------------------------------------------------------------- */

size_t UpdateKokkos::sc_sizeof(int tag)
{
  switch (tag) {
#define SC_SIZE(TAG,TYPE) case TAG: return sizeof(TYPE);
    SC_FOREACH(SC_SIZE)
#undef SC_SIZE
  }
  return 0;
}

/* ----------------------------------------------------------------------
   count the models, run their pre_collide(), and blit them to the device
   called once per move(), before the retry loop
------------------------------------------------------------------------- */

void UpdateKokkos::setup_surf_collide_models()
{
  if (surf->nsc == 0) {
    for (int t = 0; t < SC_NSTYLE; t++) nsc_style[t] = 0;
    nsc_index_cached = -1;
    return;
  }

  // index maps: which style each surf_collide is, and its slot within it.
  //   this function is called from inside move()'s migration loop, but the
  //   maps depend only on the surf_collide style list, which cannot change
  //   during a run -- so build them once and keep nsc_style[] and the host
  //   mirrors, rather than re-deriving and re-uploading both every iteration.
  //   init() clears nsc_index_cached so a new run rebuilds

  if (nsc_index_cached != surf->nsc) {
    for (int t = 0; t < SC_NSTYLE; t++) nsc_style[t] = 0;

    if ((int) d_sc_type.extent(0) < surf->nsc) {
      d_sc_type = DAT::t_int_1d("update:sc_type",surf->nsc);
      d_sc_map = DAT::t_int_1d("update:sc_map",surf->nsc);
      h_sc_type = Kokkos::create_mirror_view(d_sc_type);
      h_sc_map = Kokkos::create_mirror_view(d_sc_map);
    }

    for (int n = 0; n < surf->nsc; n++) {
      if (!surf->sc[n]->kokkosable)
        error->all(FLERR,"Must use Kokkos-enabled surface collide method with Kokkos");
      const int tag = surf_collide_style_tag(surf->sc[n]);
      if (tag < 0) error->all(FLERR,"Unknown Kokkos surface collide method");
      h_sc_type(n) = tag;
      h_sc_map(n) = nsc_style[tag]++;
    }

    Kokkos::deep_copy(d_sc_type,h_sc_type);
    Kokkos::deep_copy(d_sc_map,h_sc_map);

    nsc_index_cached = surf->nsc;
  }

  // one buffer per style, grown to hold every instance of it

  for (int t = 0; t < SC_NSTYLE; t++) {
    if (!nsc_style[t]) continue;
    const size_t need = (size_t) nsc_style[t] * sc_sizeof(t);
    if (k_sc[t].view_device().extent(0) < need) {
      k_sc[t] = DAT::tdual_char_1d("update:sc_models",need);
      d_sc[t] = k_sc[t].view_device();
    }
  }

  for (int n = 0; n < surf->nsc; n++) sc_phase(surf->sc[n],SC_PRE);

  upload_surf_collide_models();
}

/* ----------------------------------------------------------------------
   re-blit the models and push them to the device
   pre_collide(), backup() and restore() all rewrite members of the live
     model -- d_particles above all, which a grow reallocates -- so the
     device image is stale until this runs again
------------------------------------------------------------------------- */

void UpdateKokkos::upload_surf_collide_models()
{
  if (surf->nsc == 0) return;

  int slot[SC_NSTYLE];
  for (int t = 0; t < SC_NSTYLE; t++) slot[t] = 0;

  for (int n = 0; n < surf->nsc; n++) {
    const int tag = surf_collide_style_tag(surf->sc[n]);
    char *dst = k_sc[tag].view_host().data() + (size_t) slot[tag]*sc_sizeof(tag);
    switch (tag) {
#define SC_BLIT(TAG,TYPE) case TAG: sc_blit<TYPE>(dst,surf->sc[n]); break;
      SC_FOREACH(SC_BLIT)
#undef SC_BLIT
    }
    slot[tag]++;
  }

  for (int t = 0; t < SC_NSTYLE; t++) {
    if (!nsc_style[t]) continue;
    k_sc[t].modify_host();
    k_sc[t].sync_device();
    d_sc[t] = k_sc[t].view_device();
  }
}

/* ---------------------------------------------------------------------- */

void UpdateKokkos::grow_tally_computes()
{
  int ncoll = 0, nreact = 0;

  for (int m = 0; m < nsurf_tally; m++) {
    if (ComputeSurfCollisionTallyKokkos* c =
          dynamic_cast<ComputeSurfCollisionTallyKokkos*>(slist_active[m])) {
      c->grow_after_overflow();

      // growing reallocated the compute's row buffer, so the copy the kernel
      //   reads still points at the old, too-small one.  Without re-blitting
      //   it the repeated attempt overflows on the same row and the retry
      //   loop never terminates

#ifdef SPARTA_KOKKOS_FIXED_LISTS
      slist_active_coll_tally_copy[ncoll++].copy(c);
#else
      tally_buf_blit(k_slist_coll_tally,ncoll++,c);
#endif
    } else if (ComputeSurfReactionTallyKokkos* c =
               dynamic_cast<ComputeSurfReactionTallyKokkos*>(slist_active[m])) {
      c->grow_after_overflow();
#ifdef SPARTA_KOKKOS_FIXED_LISTS
      slist_active_react_tally_copy[nreact++].copy(c);
#else
      tally_buf_blit(k_slist_react_tally,nreact++,c);
#endif
    }
  }

#ifndef SPARTA_KOKKOS_FIXED_LISTS
  tally_buf_sync(k_slist_coll_tally,d_slist_coll_tally);
  tally_buf_sync(k_slist_react_tally,d_slist_react_tally);
#endif
}

/* ----------------------------------------------------------------------
   mark (mark=1) or rewind to (mark=0) the append position of every per-event
     tally compute
   the move kernel runs once per migration iteration and the tally accumulates
     across all of them, so a retried attempt must not zero the counter -- it
     rewinds to where the current iteration started, discarding only the rows
     the aborted attempt appended
------------------------------------------------------------------------- */

void UpdateKokkos::rewind_tally_computes(int mark)
{
  for (int m = 0; m < nsurf_tally; m++) {
    if (ComputeSurfCollisionTallyKokkos* c =
          dynamic_cast<ComputeSurfCollisionTallyKokkos*>(slist_active[m]))
      { if (mark) c->mark_ntally(); else c->rewind_ntally(); }
    else if (ComputeSurfReactionTallyKokkos* c =
               dynamic_cast<ComputeSurfReactionTallyKokkos*>(slist_active[m]))
      { if (mark) c->mark_ntally(); else c->rewind_ntally(); }
  }
}

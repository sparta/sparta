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

/* ---------------------------------------------------------------------- */

UpdateKokkos::UpdateKokkos(SPARTA *sparta) : Update(sparta),
  grid_kk_copy(sparta),
  domain_kk_copy(sparta),
  // Virtual functions are not yet supported on the GPU, which leads to pain:
  sc_kk_specular_copy{VAL_2(KKCopy<SurfCollideSpecularKokkos>(sparta))},
  sc_kk_diffuse_copy{VAL_2(KKCopy<SurfCollideDiffuseKokkos>(sparta))},
  sc_kk_vanish_copy{VAL_2(KKCopy<SurfCollideVanishKokkos>(sparta))},
  sc_kk_piston_copy{VAL_2(KKCopy<SurfCollidePistonKokkos>(sparta))},
  sc_kk_transparent_copy{VAL_2(KKCopy<SurfCollideTransparentKokkos>(sparta))},
  blist_active_copy{VAL_2(KKCopy<ComputeBoundaryKokkos>(sparta))},
  slist_active_copy{VAL_2(KKCopy<ComputeSurfKokkos>(sparta))}
{

  // use 1D view for scalars to reduce GPU memory operations

  d_scalars = t_int_14("collide:scalars");
  h_scalars = t_host_int_14("collide:scalars_mirror");

  d_ncomm_one     = Kokkos::subview(d_scalars,0);
  d_nexit_one     = Kokkos::subview(d_scalars,1);
  d_nboundary_one = Kokkos::subview(d_scalars,2);
  d_nmigrate      = Kokkos::subview(d_scalars,3);
  d_entryexit     = Kokkos::subview(d_scalars,4);
  d_ntouch_one    = Kokkos::subview(d_scalars,5);
  d_nscheck_one   = Kokkos::subview(d_scalars,6);
  d_nscollide_one = Kokkos::subview(d_scalars,7);
  d_nreact_one    = Kokkos::subview(d_scalars,8);
  d_nstuck        = Kokkos::subview(d_scalars,9);
  d_naxibad       = Kokkos::subview(d_scalars,10);
  d_error_flag    = Kokkos::subview(d_scalars,11);
  d_retry         = Kokkos::subview(d_scalars,12);
  d_nlocal        = Kokkos::subview(d_scalars,13);

  h_ncomm_one     = Kokkos::subview(h_scalars,0);
  h_nexit_one     = Kokkos::subview(h_scalars,1);
  h_nboundary_one = Kokkos::subview(h_scalars,2);
  h_nmigrate      = Kokkos::subview(h_scalars,3);
  h_entryexit     = Kokkos::subview(h_scalars,4);
  h_ntouch_one    = Kokkos::subview(h_scalars,5);
  h_nscheck_one   = Kokkos::subview(h_scalars,6);
  h_nscollide_one = Kokkos::subview(h_scalars,7);
  h_nreact_one    = Kokkos::subview(h_scalars,8);
  h_nstuck        = Kokkos::subview(h_scalars,9);
  h_naxibad       = Kokkos::subview(h_scalars,10);
  h_error_flag    = Kokkos::subview(h_scalars,11);
  h_retry         = Kokkos::subview(h_scalars,12);
  h_nlocal        = Kokkos::subview(h_scalars,13);

  nboundary_tally = 0;
}

/* ---------------------------------------------------------------------- */

UpdateKokkos::~UpdateKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_mlist,mlist);
  mlist = NULL;

  grid_kk_copy.uncopy();
  domain_kk_copy.uncopy();

  for (int i=0; i<KOKKOS_MAX_SURF_COLL_PER_TYPE; i++) {
    sc_kk_specular_copy[i].uncopy();
    sc_kk_diffuse_copy[i].uncopy();
    sc_kk_vanish_copy[i].uncopy();
    sc_kk_piston_copy[i].uncopy();
    sc_kk_transparent_copy[i].uncopy();
  }

  for (int i=0; i<KOKKOS_MAX_BLIST; i++) {
    blist_active_copy[i].uncopy();
  }

  for (int i=0; i<KOKKOS_MAX_SLIST; i++) {
    slist_active_copy[i].uncopy();
  }
}

/* ---------------------------------------------------------------------- */

void UpdateKokkos::init()
{
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
  }

  // choose the appropriate move method

  if (domain->dimension == 3) {
    if (surf->exist) {
      if (surf->nsr) moveptr = &UpdateKokkos::move<3,1,1,0>;
      else moveptr = &UpdateKokkos::move<3,1,0,0>;
    } else {
      if (optmove_flag) moveptr = &UpdateKokkos::move<3,0,0,1>;
      else moveptr = &UpdateKokkos::move<3,0,0,0>;
    }
  } else if (domain->axisymmetric) {
    if (surf->exist) {
      if (surf->nsr) moveptr = &UpdateKokkos::move<1,1,1,0>;
      else moveptr = &UpdateKokkos::move<1,1,0,0>;
    } else {
      if (optmove_flag) moveptr = &UpdateKokkos::move<1,0,0,1>;
      else moveptr = &UpdateKokkos::move<1,0,0,0>;
    }
  } else if (domain->dimension == 2) {
    if (surf->exist) {
      if (surf->nsr) moveptr = &UpdateKokkos::move<2,1,1,0>;
      else moveptr = &UpdateKokkos::move<2,1,0,0>;
    } else {
      if (optmove_flag) moveptr = &UpdateKokkos::move<2,0,0,1>;
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
  hash_kk = grid_kk->hash_kk;

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

    ntimestep++;

    if (collide_react) collide_react_reset();
    if (bounce_tally) bounce_set(ntimestep);

    timer->stamp();

    // dynamic parameter updates

    if (dynamic) dynamic_update();

    // start of step fixes

    if (n_start_of_step) {
      modify->start_of_step();
      timer->stamp(TIME_MODIFY);
    }

    // move particles

    if (cellweightflag) particle->pre_weight();
    (this->*moveptr)();
    timer->stamp(TIME_MOVE);

    // communicate particles

    if (nmigrate) {
      k_mlist_small = Kokkos::subview(k_mlist,std::make_pair(0,nmigrate));
      k_mlist_small.sync_host();
    }
    auto mlist_small = k_mlist_small.h_view.data();

    ((CommKokkos*)comm)->migrate_particles(nmigrate,mlist_small,k_mlist_small.d_view);
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

    // all output

    if (ntimestep == output->next) {
      particle_kk->sync(Host,ALL_MASK);
      output->write(ntimestep);
      timer->stamp(TIME_OUTPUT);
    }
  }
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

    d_particles = particle_kk->k_particles.d_view;

    GridKokkos* grid_kk = ((GridKokkos*)grid);
    d_cells = grid_kk->k_cells.d_view;
    d_sinfo = grid_kk->k_sinfo.d_view;
    d_pcells = grid_kk->k_pcells.d_view;

    d_csurfs = grid_kk->d_csurfs;
    d_csplits = grid_kk->d_csplits;
    d_csubs = grid_kk->d_csubs;

    if (surf->exist) {
      SurfKokkos* surf_kk = ((SurfKokkos*)surf);
      surf_kk->sync(Device,ALL_MASK);
      d_lines = surf_kk->k_lines.d_view;
      d_tris = surf_kk->k_tris.d_view;
    }

    if (surf->nsr) {
      double extra_factor = 1.0;
      if (!sparta->kokkos->react_retry_flag)
        extra_factor = sparta->kokkos->react_extra;

      int nlocal_extra = particle->nlocal*extra_factor;
      if (d_particles.extent(0) < nlocal_extra) {
        particle->grow(nlocal_extra - particle->nlocal); // this!
        d_particles = particle_kk->k_particles.d_view;
      }
    }

    particle_kk->sync(Device,PARTICLE_MASK);
    grid_kk->sync(Device,CELL_MASK|PCELL_MASK|SINFO_MASK|PLEVEL_MASK);

    // may be able to move this outside of the while loop
    grid_kk_copy.copy(grid_kk);
    domain_kk_copy.copy((DomainKokkos*)domain);

    if (surf->nsc > KOKKOS_MAX_TOT_SURF_COLL)
      error->all(FLERR,"Kokkos currently supports two instances of each surface collide method");

    if (surf->nsc > 0) {
      int nspec,ndiff,nvan,npist,ntrans;
      nspec = ndiff = nvan = npist = ntrans = 0;
      for (int n = 0; n < surf->nsc; n++) {
        if (!surf->sc[n]->kokkosable)
          error->all(FLERR,"Must use Kokkos-enabled surface collide method with Kokkos");
        if (strcmp(surf->sc[n]->style,"specular") == 0) {
          sc_kk_specular_copy[nspec].copy((SurfCollideSpecularKokkos*)(surf->sc[n]));
          sc_kk_specular_copy[nspec].obj.pre_collide();
          sc_type_list[n] = 0;
          sc_map[n] = nspec;
          nspec++;
        } else if (strcmp(surf->sc[n]->style,"diffuse") == 0) {
          sc_kk_diffuse_copy[ndiff].copy((SurfCollideDiffuseKokkos*)(surf->sc[n]));
          sc_kk_diffuse_copy[ndiff].obj.pre_collide();
          sc_type_list[n] = 1;
          sc_map[n] = ndiff;
          ndiff++;
        } else if (strcmp(surf->sc[n]->style,"vanish") == 0) {
          sc_kk_vanish_copy[nvan].copy((SurfCollideVanishKokkos*)(surf->sc[n]));
          sc_kk_vanish_copy[nvan].obj.pre_collide();
          sc_type_list[n] = 2;
          sc_map[n] = nvan;
          nvan++;
        } else if (strcmp(surf->sc[n]->style,"piston") == 0) {
          sc_kk_piston_copy[npist].copy((SurfCollidePistonKokkos*)(surf->sc[n]));
          sc_kk_piston_copy[npist].obj.pre_collide();
          sc_type_list[n] = 3;
          sc_map[n] = npist;
          npist++;
        } else if (strcmp(surf->sc[n]->style,"transparent") == 0) {
          sc_kk_transparent_copy[ntrans].copy((SurfCollideTransparentKokkos*)(surf->sc[n]));
          sc_kk_transparent_copy[ntrans].obj.pre_collide();
          sc_type_list[n] = 4;
          sc_map[n] = ntrans;
          ntrans++;
        } else {
          error->all(FLERR,"Unknown Kokkos surface collide method");
        }
      }
      if (nspec > KOKKOS_MAX_SURF_COLL_PER_TYPE || ndiff > KOKKOS_MAX_SURF_COLL_PER_TYPE ||
          nvan > KOKKOS_MAX_SURF_COLL_PER_TYPE || npist > KOKKOS_MAX_SURF_COLL_PER_TYPE ||
          ntrans > KOKKOS_MAX_SURF_COLL_PER_TYPE)
        error->all(FLERR,"Kokkos currently supports two instances of each surface collide method");
    }

    Kokkos::deep_copy(h_scalars,0);

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

    while (h_retry()) {

      if (surf->nsr && sparta->kokkos->react_retry_flag)
        backup();

      h_retry() = 0;
      h_nlocal() = particle->nlocal;
      if (continue_loop_flag) h_nmigrate() = nmigrate;

      Kokkos::deep_copy(d_scalars,h_scalars);

      copymode = 1;

    /* ATOMIC_REDUCTION: 1 = use atomics
                         0 = don't need atomics
                        -1 = use parallel_reduce
    */

#if defined SPARTA_KOKKOS_GPU
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagUpdateMove<DIM,SURF,REACT,OPT,1> >(pstart,pstop),*this);
#elif defined KOKKOS_ENABLE_SERIAL
      if constexpr(std::is_same<DeviceType,Kokkos::Serial>::value)
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagUpdateMove<DIM,SURF,REACT,OPT,-1> >(pstart,pstop),*this,reduce);
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

      if (h_retry()) {
        int nlocal_new = h_nlocal();

        if (!sparta->kokkos->react_retry_flag) {
          error->one(FLERR,"Ran out of space for Kokkos reactions, increase react/extra"
                           " or use react/retry");
        } else
          restore();

        //  reset counters

        Kokkos::deep_copy(h_scalars,0);
        reduce = UPDATE_REDUCE();
        h_retry() = 1;

        if (d_particles.extent(0) < nlocal_new) {
          particle->grow(nlocal_new - particle->nlocal);
          d_particles = particle_kk->k_particles.d_view;
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
      sprintf(str,
              "Particle being sent to self proc "
              "on step " BIGINT_FORMAT,
              update->ntimestep);
      error->one(FLERR,str);
    }

    if (surf->nsc > 0) {
      int nspec,ndiff,nvan,npist,ntrans;
      nspec = ndiff = nvan = npist = ntrans = 0;
      for (int n = 0; n < surf->nsc; n++) {
        if (strcmp(surf->sc[n]->style,"specular") == 0) {
          sc_kk_specular_copy[nspec].obj.post_collide();
          nspec++;
        } else if (strcmp(surf->sc[n]->style,"diffuse") == 0) {
          sc_kk_diffuse_copy[ndiff].obj.post_collide();
          ndiff++;
        } else if (strcmp(surf->sc[n]->style,"vanish") == 0) {
          sc_kk_vanish_copy[nvan].obj.post_collide();
          nvan++;
        } else if (strcmp(surf->sc[n]->style,"piston") == 0) {
          sc_kk_piston_copy[npist].obj.post_collide();
          npist++;
        } else if (strcmp(surf->sc[n]->style,"transparent") == 0) {
          sc_kk_transparent_copy[ntrans].obj.post_collide();
          ntrans++;
        }
      }
    }

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
    timer->stamp();

    if (any_entryexit) {
      if (nmigrate) {
        k_mlist_small = Kokkos::subview(k_mlist,std::make_pair(0,nmigrate));
        k_mlist_small.sync_host();
      }
      auto mlist_small = k_mlist_small.h_view.data();
      timer->stamp(TIME_MOVE);
      pstart = ((CommKokkos*)comm)->migrate_particles(nmigrate,mlist_small,k_mlist_small.d_view);
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

  particle->sorted = 0;
  particle_kk->sorted_kk = 0;

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

  if (nsurf_tally) {
    for (int m = 0; m < nsurf_tally; m++) {
      ComputeSurfKokkos* compute_surf_kk = (ComputeSurfKokkos*)(slist_active[m]);
      compute_surf_kk->post_surf_tally();
    }
  }

  if (nboundary_tally) {
    for (int m = 0; m < nboundary_tally; m++) {
      ComputeBoundaryKokkos* compute_boundary_kk = (ComputeBoundaryKokkos*)(blist_active[m]);
      compute_boundary_kk->post_boundary_tally();
    }
  }
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
      if (DIM == 3) field3d(dtremain,xnew,v);
      else if (DIM == 2) field2d(dtremain,xnew,v);
    } else if (fstyle == PFIELD) field_per_particle(i,particle_i.icell,dtremain,xnew,v);
    else if (fstyle == GFIELD) field_per_grid(i,particle_i.icell,dtremain,xnew,v);
  } else if (pflag == PINSERT) {
    dtremain = particle_i.dtremain;
    xnew[0] = x[0] + dtremain*v[0];
    xnew[1] = x[1] + dtremain*v[1];
    if (DIM != 2) xnew[2] = x[2] + dtremain*v[2];
    if (fstyle == CFIELD) {
      if (DIM == 3) field3d(dtremain,xnew,v);
      else if (DIM == 2) field2d(dtremain,xnew,v);
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

  if (OPT) {
    int optmove = 1;

    if (xnew[0] < xlo || xnew[0] > xhi)
      optmove = 0;

    if (xnew[1] < ylo || xnew[1] > yhi)
      optmove = 0;

    if (DIM == 3) {
      if (xnew[2] < zlo || xnew[2] > zhi)
        optmove = 0;
    }

    if (optmove) {

      const int ip = static_cast<int>((xnew[0] - xlo)/dx);
      const int jp = static_cast<int>((xnew[1] - ylo)/dy);
      int kp = 0;
      if (DIM == 3) kp = static_cast<int>((xnew[2] - zlo)/dz);

      int cellIdx = (kp*ncy + jp)*ncx + ip + 1;
      auto index = hash_kk.find(static_cast<GridKokkos::key_type>(cellIdx));

      // particle moving outside ghost halo will be flagged for standard move

      if (hash_kk.valid_at(index)) {

        int icell = static_cast<int>(hash_kk.value_at(index));

        // reset particle cell and coordinates

        particle_i.icell = icell;
        particle_i.flag = PKEEP;
        x[0] = xnew[0];
        x[1] = xnew[1];
        x[2] = xnew[2];

        if (d_cells[icell].proc != me) {
          int indx;
          if (ATOMIC_REDUCTION == 0) {
            indx = d_nmigrate();
            d_nmigrate()++;
          } else {
            indx = Kokkos::atomic_fetch_add(&d_nmigrate(),1);
          }
          k_mlist.d_view[indx] = i;

          particle_i.flag = PDONE;

	  if (ATOMIC_REDUCTION == 1)
	    Kokkos::atomic_increment(&d_ncomm_one());
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
    Kokkos::atomic_increment(&d_ntouch_one());
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
      frac = (lo[0]-x[0]) / (xnew[0]-x[0]);
      outface = XLO;
    } else if (xnew[0] >= hi[0]) {
      frac = (hi[0]-x[0]) / (xnew[0]-x[0]);
      outface = XHI;
    }

    if (DIM != 1) {
      if (xnew[1] < lo[1]) {
        newfrac = (lo[1]-x[1]) / (xnew[1]-x[1]);
        if (newfrac < frac) {
          frac = newfrac;
          outface = YLO;
        }
      } else if (xnew[1] >= hi[1]) {
        newfrac = (hi[1]-x[1]) / (xnew[1]-x[1]);
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
        newfrac = (lo[2]-x[2]) / (xnew[2]-x[2]);
        if (newfrac < frac) {
          frac = newfrac;
          outface = ZLO;
        }
      } else if (xnew[2] >= hi[2]) {
        newfrac = (hi[2]-x[2]) / (xnew[2]-x[2]);
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
          // tally surface statistics if requested using iorig

          ipart = &particle_i;
          ipart->icell = icell;
          dtremain *= 1.0 - minparam*frac;

          if (nsurf_tally)
            iorig = particle_i;
          int n = DIM == 3 ? tri->isc : line->isc;
          int sc_type = sc_type_list[n];
          int m = sc_map[n];

          if (DIM == 3) {
            if (sc_type == 0) {
              jpart = sc_kk_specular_copy[m].obj.
                collide_kokkos<REACT,ATOMIC_REDUCTION>(ipart,dtremain,minsurf,tri->norm,tri->isr,reaction,d_retry,d_nlocal);
            } else if (sc_type == 1) {
              jpart = sc_kk_diffuse_copy[m].obj.
                collide_kokkos<REACT,ATOMIC_REDUCTION>(ipart,dtremain,minsurf,tri->norm,tri->isr,reaction,d_retry,d_nlocal);
            } else if (sc_type == 2) {
              jpart = sc_kk_vanish_copy[m].obj.
                collide_kokkos<REACT,ATOMIC_REDUCTION>(ipart,dtremain,minsurf,tri->norm,tri->isr,reaction,d_retry,d_nlocal);
            } else if (sc_type == 3) {
              jpart = sc_kk_piston_copy[m].obj.
                collide_kokkos<REACT,ATOMIC_REDUCTION>(ipart,dtremain,minsurf,tri->norm,tri->isr,reaction,d_retry,d_nlocal);
            } else if (sc_type == 4) {
              jpart = sc_kk_transparent_copy[m].obj.
                collide_kokkos<REACT,ATOMIC_REDUCTION>(ipart,dtremain,minsurf,tri->norm,tri->isr,reaction,d_retry,d_nlocal);
            }
          }

          if (DIM != 3) {
            if (sc_type == 0) {
              jpart = sc_kk_specular_copy[m].obj.
                collide_kokkos<REACT,ATOMIC_REDUCTION>(ipart,dtremain,minsurf,line->norm,line->isr,reaction,d_retry,d_nlocal);
            } else if (sc_type == 1) {
              jpart = sc_kk_diffuse_copy[m].obj.
                collide_kokkos<REACT,ATOMIC_REDUCTION>(ipart,dtremain,minsurf,line->norm,line->isr,reaction,d_retry,d_nlocal);
            } else if (sc_type == 2) {
              jpart = sc_kk_vanish_copy[m].obj.
                collide_kokkos<REACT,ATOMIC_REDUCTION>(ipart,dtremain,minsurf,line->norm,line->isr,reaction,d_retry,d_nlocal);
            } else if (sc_type == 3) {
              jpart = sc_kk_piston_copy[m].obj.
                collide_kokkos<REACT,ATOMIC_REDUCTION>(ipart,dtremain,minsurf,line->norm,line->isr,reaction,d_retry,d_nlocal);
            } else if (sc_type == 4) {
              jpart = sc_kk_transparent_copy[m].obj.
                collide_kokkos<REACT,ATOMIC_REDUCTION>(ipart,dtremain,minsurf,line->norm,line->isr,reaction,d_retry,d_nlocal);
            }
          }

          if (jpart) {
            x = particle_i.x;
            v = particle_i.v;
            jpart->flag = PSURF + 1 + minsurf;
            jpart->dtremain = dtremain;
            jpart->weight = particle_i.weight;
          }

          if (nsurf_tally)
            for (m = 0; m < nsurf_tally; m++)
              slist_active_copy[m].obj.
                    surf_tally_kk<ATOMIC_REDUCTION>(minsurf,icell,reaction,&iorig,ipart,jpart);

          // stuck_iterate = consecutive iterations particle is immobile

          if (minparam == 0.0) stuck_iterate++;
          else stuck_iterate = 0;

          // reset post-bounce xnew

          xnew[0] = x[0] + dtremain*v[0];
          xnew[1] = x[1] + dtremain*v[1];
          if (DIM != 2) xnew[2] = x[2] + dtremain*v[2];

          exclude = minsurf;
          if (ATOMIC_REDUCTION == 1)
            Kokkos::atomic_increment(&d_nscollide_one());
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
              Kokkos::atomic_increment(&d_nstuck());
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
            Kokkos::atomic_increment(&d_naxibad());
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

        int n = domain_kk_copy.obj.surf_collide[outface];
        int sc_type = sc_type_list[n];
        int m = sc_map[n];

        if (sc_type == 0)
          jpart = sc_kk_specular_copy[m].obj.
            collide_kokkos<REACT,ATOMIC_REDUCTION>(ipart,dtremain,-(outface+1),domain_kk_copy.obj.norm[outface],domain_kk_copy.obj.surf_react[outface],reaction,d_retry,d_nlocal);
        else if (sc_type == 1)
          jpart = sc_kk_diffuse_copy[m].obj.
            collide_kokkos<REACT,ATOMIC_REDUCTION>(ipart,dtremain,-(outface+1),domain_kk_copy.obj.norm[outface],domain_kk_copy.obj.surf_react[outface],reaction,d_retry,d_nlocal);
        else if (sc_type == 2)
          jpart = sc_kk_vanish_copy[m].obj.
            collide_kokkos<REACT,ATOMIC_REDUCTION>(ipart,dtremain,-(outface+1),domain_kk_copy.obj.norm[outface],domain_kk_copy.obj.surf_react[outface],reaction,d_retry,d_nlocal);
        else if (sc_type == 3)
          jpart = sc_kk_piston_copy[m].obj.
            collide_kokkos<REACT,ATOMIC_REDUCTION>(ipart,dtremain,-(outface+1),domain_kk_copy.obj.norm[outface],domain_kk_copy.obj.surf_react[outface],reaction,d_retry,d_nlocal);
        else if (sc_type == 4)
          jpart = sc_kk_transparent_copy[m].obj.
            collide_kokkos<REACT,ATOMIC_REDUCTION>(ipart,dtremain,-(outface+1),domain_kk_copy.obj.norm[outface],domain_kk_copy.obj.surf_react[outface],reaction,d_retry,d_nlocal);

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

      if (nboundary_tally)
        for (int m = 0; m < nboundary_tally; m++)
          blist_active_copy[m].obj.
            boundary_tally_kk<ATOMIC_REDUCTION>(outface,bflag,reaction,&iorig,ipart,jpart,domain_kk_copy.obj.norm[outface]);

      if (DIM == 1) {
        xnew[0] = x[0] + dtremain*v[0];
        xnew[1] = x[1] + dtremain*v[1];
        xnew[2] = x[2] + dtremain*v[2];
      }

      if (bflag == OUTFLOW) {
        particle_i.flag = PDISCARD;
        if (ATOMIC_REDUCTION == 1)
          Kokkos::atomic_increment(&d_nexit_one());
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

        Kokkos::atomic_increment(&d_nboundary_one());
        Kokkos::atomic_decrement(&d_ntouch_one());    // decrement here since will increment below
      } else {
        if (ATOMIC_REDUCTION == 1) {
          Kokkos::atomic_increment(&d_nboundary_one());
          Kokkos::atomic_decrement(&d_ntouch_one());    // decrement here since will increment below
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
      Kokkos::atomic_increment(&d_ntouch_one());
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
    k_mlist.d_view[index] = i;
    if (particle_i.flag != PDISCARD) {
      if (d_cells[icell].proc == me && !d_error_flag()) {
        d_error_flag() = 1;
        return;
      }
      if (ATOMIC_REDUCTION == 1)
        Kokkos::atomic_increment(&d_ncomm_one());
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

void UpdateKokkos::bounce_set(bigint ntimestep)
{
  Update::bounce_set(ntimestep);

  int i;

  if (nboundary_tally > KOKKOS_MAX_BLIST)
    error->all(FLERR,"Kokkos currently only supports two instances of compute boundary");

  if (nboundary_tally) {
    for (i = 0; i < nboundary_tally; i++) {
      ComputeBoundaryKokkos* compute_boundary_kk = (ComputeBoundaryKokkos*)(blist_active[i]);
      compute_boundary_kk->pre_boundary_tally();
      blist_active_copy[i].copy(compute_boundary_kk);
    }
  }

  if (nsurf_tally > KOKKOS_MAX_SLIST)
    error->all(FLERR,"Kokkos currently only supports two instances of compute surface");

  if (nsurf_tally) {
    for (i = 0; i < nsurf_tally; i++) {
      if (strcmp(slist_active[i]->style,"isurf/grid") == 0)
        error->all(FLERR,"Kokkos doesn't yet support compute isurf/grid");
      ComputeSurfKokkos* compute_surf_kk = (ComputeSurfKokkos*)(slist_active[i]);
      compute_surf_kk->pre_surf_tally();
      slist_active_copy[i].copy(compute_surf_kk);
    }
  }
}

/* ---------------------------------------------------------------------- */

void UpdateKokkos::backup()
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  d_particles = particle_kk->k_particles.d_view;
  d_particles_backup = decltype(d_particles)(Kokkos::view_alloc("update:particles_backup",Kokkos::WithoutInitializing),d_particles.extent(0));

  Kokkos::deep_copy(d_particles_backup,d_particles);

  if (surf->nsc > 0) {
    int nspec,ndiff,npist;
    nspec = ndiff = npist = 0;
    for (int n = 0; n < surf->nsc; n++) {
      if (strcmp(surf->sc[n]->style,"specular") == 0) {
        sc_kk_specular_copy[nspec].obj.backup();
        nspec++;
      } else if (strcmp(surf->sc[n]->style,"diffuse") == 0) {
        sc_kk_diffuse_copy[ndiff].obj.backup();
        ndiff++;
      } else if (strcmp(surf->sc[n]->style,"piston") == 0) {
        sc_kk_piston_copy[npist].obj.backup();
        npist++;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void UpdateKokkos::restore()
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  Kokkos::deep_copy(particle_kk->k_particles.d_view,d_particles_backup);
  d_particles = particle_kk->k_particles.d_view;

  if (surf->nsc > 0) {
    int nspec,ndiff,npist;
    nspec = ndiff = npist = 0;
    for (int n = 0; n < surf->nsc; n++) {
      if (strcmp(surf->sc[n]->style,"specular") == 0) {
        sc_kk_specular_copy[nspec].obj.restore();
        nspec++;
      } else if (strcmp(surf->sc[n]->style,"diffuse") == 0) {
        sc_kk_diffuse_copy[ndiff].obj.restore();
        ndiff++;
      } else if (strcmp(surf->sc[n]->style,"piston") == 0) {
        sc_kk_piston_copy[npist].obj.restore();
        npist++;
      }
    }
  }

  // deallocate references to reduce memory use

  d_particles_backup = decltype(d_particles_backup)();
}

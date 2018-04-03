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

using namespace SPARTA_NS;

enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};         // same as Domain
enum{PERIODIC,OUTFLOW,REFLECT,SURFACE,AXISYM};  // same as Domain
//enum{OUTSIDE,INSIDE,ONSURF2OUT,ONSURF2IN};      // several files
enum{PKEEP,PINSERT,PDONE,PDISCARD,PENTRY,PEXIT,PSURF};   // several files
enum{NCHILD,NPARENT,NUNKNOWN,NPBCHILD,NPBPARENT,NPBUNKNOWN,NBOUND};  // Grid
enum{TALLYAUTO,TALLYREDUCE,TALLYLOCAL};         // same as Surf

#define MAXSTUCK 20
#define EPSPARAM 1.0e-7

// either set ID or PROC/INDEX, set other to -1

//#define MOVE_DEBUG 1              // un-comment to debug one particle
#define MOVE_DEBUG_ID 537898343   // particle ID
#define MOVE_DEBUG_PROC 0        // owning proc
#define MOVE_DEBUG_INDEX 16617       // particle index on owning proc
#define MOVE_DEBUG_STEP 16    // timestep

#define VAL_1(X) X
#define VAL_2(X) VAL_1(X), VAL_1(X)

/* ---------------------------------------------------------------------- */

UpdateKokkos::UpdateKokkos(SPARTA *sparta) : Update(sparta),
  grid_kk_copy(sparta),
  domain_kk_copy(sparta),
  //// Virtual functions are not yet supported on the GPU, which leads to pain:
  sc_kk_specular_copy{VAL_2(KKCopy<SurfCollideSpecularKokkos>(sparta))},
  sc_kk_diffuse_copy{VAL_2(KKCopy<SurfCollideDiffuseKokkos>(sparta))},
  sc_kk_vanish_copy{VAL_2(KKCopy<SurfCollideVanishKokkos>(sparta))},
  blist_active_copy{VAL_2(KKCopy<ComputeBoundaryKokkos>(sparta))}
{
  k_ncomm_one = DAT::tdual_int_scalar("UpdateKokkos:ncomm_one");
  d_ncomm_one = k_ncomm_one.view<DeviceType>();
  h_ncomm_one = k_ncomm_one.h_view;

  k_nexit_one = DAT::tdual_int_scalar("UpdateKokkos:nexit_one");
  d_nexit_one = k_nexit_one.view<DeviceType>();
  h_nexit_one = k_nexit_one.h_view;

  k_nboundary_one = DAT::tdual_int_scalar("UpdateKokkos:nboundary_one");
  d_nboundary_one = k_nboundary_one.view<DeviceType>();
  h_nboundary_one = k_nboundary_one.h_view;

  k_nmigrate = DAT::tdual_int_scalar("UpdateKokkos:nmigrate");
  d_nmigrate = k_nmigrate.view<DeviceType>();
  h_nmigrate = k_nmigrate.h_view;

  k_entryexit = DAT::tdual_int_scalar("UpdateKokkos:entryexit");
  d_entryexit = k_entryexit.view<DeviceType>();
  h_entryexit = k_entryexit.h_view;

  k_ntouch_one = DAT::tdual_int_scalar("UpdateKokkos:ntouch_one");
  d_ntouch_one = k_ntouch_one.view<DeviceType>();
  h_ntouch_one = k_ntouch_one.h_view;

  k_nscheck_one = DAT::tdual_int_scalar("UpdateKokkos:nscheck_one");
  d_nscheck_one = k_nscheck_one.view<DeviceType>();
  h_nscheck_one = k_nscheck_one.h_view;

  k_nscollide_one = DAT::tdual_int_scalar("UpdateKokkos:nscollide_one");
  d_nscollide_one = k_nscollide_one.view<DeviceType>();
  h_nscollide_one = k_nscollide_one.h_view;

  k_nreact_one = DAT::tdual_int_scalar("UpdateKokkos:nreact_one");
  d_nreact_one = k_nreact_one.view<DeviceType>();
  h_nreact_one = k_nreact_one.h_view;

  k_nstuck = DAT::tdual_int_scalar("UpdateKokkos:nstuck");
  d_nstuck = k_nstuck.view<DeviceType>();
  h_nstuck = k_nstuck.h_view;

  k_error_flag = DAT::tdual_int_scalar("UpdateKokkos:error_flag");
  d_error_flag = k_error_flag.view<DeviceType>();
  h_error_flag = k_error_flag.h_view;

  boundary_tally = 0;
}

/* ---------------------------------------------------------------------- */

UpdateKokkos::~UpdateKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_mlist,mlist);
  mlist = NULL;

  grid_kk_copy.uncopy();
  domain_kk_copy.uncopy();

  if (surf->nsc > 0) {
    for (int i=0; i<surf->nsc; i++) {
      sc_kk_specular_copy[i].uncopy();
      sc_kk_diffuse_copy[i].uncopy();
      sc_kk_vanish_copy[i].uncopy();
    }
  }

  if (nboundary_tally > 0) {
    for (int i=0; i<nboundary_tally; i++) {
      blist_active_copy[i].uncopy();
    }
  }

}

/* ---------------------------------------------------------------------- */

void UpdateKokkos::init()
{
  // init the UpdateKokkos class if performing a run, else just return
  // only set first_update if a run is being performed

  if (runflag == 0) return;
  first_update = 1;

  // choose the appropriate move method

  moveptr = NULL;
  if (domain->dimension == 3) {
    if (surf->exist) moveptr = &UpdateKokkos::move<3,1>;
    else moveptr = &UpdateKokkos::move<3,0>;
  } else if (domain->axisymmetric) {
    if (surf->exist) moveptr = &UpdateKokkos::move<1,1>;
    else moveptr = &UpdateKokkos::move<1,0>;
  } else if (domain->dimension == 2) {
    if (surf->exist) moveptr = &UpdateKokkos::move<2,1>;
    else moveptr = &UpdateKokkos::move<2,0>;
  }

  // check gravity vector

  if (domain->dimension == 2 && gravity[2] != 0.0)
    error->all(FLERR,"Gravity in z not allowed for 2d");
  if (domain->axisymmetric && gravity[1] != 0.0) 
    error->all(FLERR,"Gravity in y not allowed for axi-symmetric model");

  // moveperturb method is set if particle motion is perturbed

  gravity_3d_flag = gravity_2d_flag = 0;
  moveperturb = NULL;
  if (gravity[0] != 0.0 || gravity[1] != 0.0 || gravity[2] != 0.0) {
    if (domain->dimension == 3) {
      moveperturb = &UpdateKokkos::gravity3d;
      gravity_3d_flag = 1;
    }
    if (domain->dimension == 2){
      moveperturb = &UpdateKokkos::gravity2d;
      gravity_2d_flag = 1;
    }
  }

  if (moveperturb) perturbflag = 1;
  else perturbflag = 0;

}

/* ---------------------------------------------------------------------- */

void UpdateKokkos::setup()
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  GridKokkos* grid_kk = (GridKokkos*) grid;
  SurfKokkos* surf_kk = (SurfKokkos*) surf;

  if (sparta->kokkos->prewrap) {

    // particle

    particle_kk->wrap_kokkos();

    // grid

    grid_kk->wrap_kokkos();

    // surf

    if (surf->exist)
      surf_kk->wrap_kokkos();

    sparta->kokkos->prewrap = 0;
  } else {
    particle_kk->modify(Host,ALL_MASK);
    grid_kk->modify(Host,ALL_MASK);
    if (surf->exist) {
      surf_kk->modify(Host,ALL_MASK);
      grid_kk->update_hash();
      grid_kk->wrap_kokkos_graphs();
    }
  }

  Update::setup(); // must come after prewrap since computes are called by setup()

  //// For MPI debugging
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
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;

  int n_start_of_step = modify->n_start_of_step;
  int n_end_of_step = modify->n_end_of_step;
  //int dynamic = 0;

  // cellweightflag = 1 if grid-based particle weighting is ON

  int cellweightflag = 0;
  if (grid->cellweightflag) cellweightflag = 1;

  // loop over timesteps

  for (int i = 0; i < nsteps; i++) {

    ntimestep++;
    if (collide_react) collide_react_update();
    if (bounce_tally) bounce_set(ntimestep);

    timer->stamp();

    // start of step fixes

    if (n_start_of_step) {
      modify->start_of_step();
      timer->stamp(TIME_MODIFY);
    }


    //if (dynamic) domain->dynamic();

    // move particles

    if (cellweightflag) particle->pre_weight();
    (this->*moveptr)();
    timer->stamp(TIME_MOVE);

    // communicate particles

    DAT::t_int_1d d_mlist_small = Kokkos::subview(k_mlist.d_view,std::make_pair(0,nmigrate));
    HAT::t_int_1d h_mlist_small = HAT::t_int_1d(Kokkos::view_alloc("mlist_mirror",Kokkos::WithoutInitializing),nmigrate);
    Kokkos::deep_copy(h_mlist_small, d_mlist_small);
    auto mlist_small = h_mlist_small.data();

    ((CommKokkos*)comm)->migrate_particles(nmigrate,mlist_small,d_mlist_small);
    if (cellweightflag) particle->post_weight();
    timer->stamp(TIME_COMM);

    if (collide) {
      particle_kk->sort_kokkos();
      timer->stamp(TIME_SORT);

      //CollideVSSKokkos* collide_kk = (CollideVSSKokkos*) collide;
      //collide_kk->wrap_kokkos_sort();

      collide->collisions();
      timer->stamp(TIME_COLLIDE);
    }

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
  particle_kk->sync(Host,ALL_MASK);

}

//make randomread versions of d_particles?

/* ----------------------------------------------------------------------
   advect particles thru grid
   DIM = 2/3 for 2d/3d, 1 for 2d axisymmetric
   SURF = 0/1 for no surfs or surfs
   use multiple iterations of move/comm if necessary
------------------------------------------------------------------------- */

template < int DIM, int SURF > void UpdateKokkos::move()
{
  //bool hitflag;
  //int m,icell,icell_original,nmask,outface,bflag,nflag,pflag,itmp;
  //int side,minside,minsurf,nsurf,cflag,isurf,exclude,stuck_iterate;
  //int pstart,pstop,entryexit,any_entryexit;
  //cellint *neigh;
  //double dtremain,frac,newfrac,param,minparam,rnew,dtsurf,tc,tmp;
  //double xnew[3],xhold[3],xc[3],vc[3],minxc[3],minvc[3];
  //double *x,*v,*lo,*hi;
  //Particle::OnePart iorig;
  //Particle::OnePart *particles;
  //Particle::OnePart *ipart,*jpart;

  int pstart,pstop,entryexit,any_entryexit;

  // extend migration list if necessary

  int nlocal = particle->nlocal;
  int maxlocal = particle->maxlocal;

  if (nlocal > maxmigrate) {
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

  if (sparta->kokkos->atomic_reduction) {
    h_ntouch_one() = 0;
    k_ntouch_one.modify<SPAHostType>();
    k_ntouch_one.sync<DeviceType>();

    h_nexit_one() = 0;
    k_nexit_one.modify<SPAHostType>();
    k_nexit_one.sync<DeviceType>();

    h_nboundary_one() = 0;
    k_nboundary_one.modify<SPAHostType>();
    k_nboundary_one.sync<DeviceType>();

    h_ncomm_one() = 0;
    k_ncomm_one.modify<SPAHostType>();
    k_ncomm_one.sync<DeviceType>();

    h_nscheck_one() = 0;
    k_nscheck_one.modify<SPAHostType>();
    k_nscheck_one.sync<DeviceType>();

    h_nscollide_one() = 0;
    k_nscollide_one.modify<SPAHostType>();
    k_nscollide_one.sync<DeviceType>();

    h_nreact_one() = 0;
    k_nreact_one.modify<SPAHostType>();
    k_nreact_one.sync<DeviceType>();
  }

  h_error_flag() = 0;
  k_error_flag.modify<SPAHostType>();
  k_error_flag.sync<DeviceType>();

  // move/migrate iterations

  dt = update->dt;
  int notfirst = 0;

  ParticleKokkos* particle_kk = ((ParticleKokkos*)particle);

  while (1) {

    // loop over particles
    // first iteration = all my particles
    // subsequent iterations = received particles

    niterate++;

    d_particles = particle_kk->k_particles.d_view;
    
    GridKokkos* grid_kk = ((GridKokkos*)grid);
    d_cells = grid_kk->k_cells.d_view;
    d_sinfo = grid_kk->k_sinfo.d_view;

    d_csurfs = grid_kk->d_csurfs;
    d_csplits = grid_kk->d_csplits;
    d_csubs = grid_kk->d_csubs;

    if (surf->exist) {
      SurfKokkos* surf_kk = ((SurfKokkos*)surf);
      surf_kk->sync(Device,ALL_MASK);
      d_pts = surf_kk->k_pts.d_view;
      d_lines = surf_kk->k_lines.d_view;
      d_tris = surf_kk->k_tris.d_view;
    }

    particle_kk->sync(Device,PARTICLE_MASK);
    grid_kk->sync(Device,CELL_MASK|SINFO_MASK);

    // may be able to move this outside of the while loop
    grid_kk_copy.copy(grid_kk);
    domain_kk_copy.copy((DomainKokkos*)domain);

    if (surf->nsc > KOKKOS_TOT_SURF_COLL)
      error->all(FLERR,"Kokkos currently supports two instances of each surface collide method");

    if (surf->nsc > 0) {
      int nspec,ndiff,nvan;
      nspec = ndiff = nvan = 0;
      for (int n = 0; n < surf->nsc; n++) {
        if (strcmp(surf->sc[n]->style,"specular") == 0) {
          sc_kk_specular_copy[nspec].copy((SurfCollideSpecularKokkos*)(surf->sc[n]));
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
          sc_type_list[n] = 2;
          sc_map[nvan] = nvan;
          nvan++;
        } else {
          error->all(FLERR,"Unknown Kokkos surface collide method");
        }
      }
      if (nspec > KOKKOS_MAX_SURF_COLL_PER_TYPE || ndiff > KOKKOS_MAX_SURF_COLL_PER_TYPE || nvan > KOKKOS_MAX_SURF_COLL_PER_TYPE)
        error->all(FLERR,"Kokkos currently supports two instances of each surface collide method");
    }

    if (nsurf_tally)
      error->all(FLERR,"Cannot (yet) use surface tallying compute with Kokkos");

    if (surf->nsr)
      error->all(FLERR,"Cannot (yet) use surface reactions with Kokkos");

    h_nmigrate() = 0;
    k_nmigrate.modify<SPAHostType>();
    k_nmigrate.sync<DeviceType>();

    h_entryexit() = 0;
    k_entryexit.modify<SPAHostType>();
    k_entryexit.sync<DeviceType>();

    nmigrate = 0;
    entryexit = 0;

    if (notfirst == 0) {
      notfirst = 1;
      pstart = 0;
      pstop = nlocal;
    }

    UPDATE_REDUCE reduce;

    /* ATOMIC_REDUCTION: 1 = use atomics
                         0 = don't need atomics
                        -1 = use parallel_reduce
    */

    k_mlist.sync<SPADeviceType>();
    copymode = 1;
    if (sparta->kokkos->atomic_reduction) {
      if (sparta->kokkos->need_atomics) {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagUpdateMove<DIM,SURF,1> >(pstart,pstop),*this);
      } else {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagUpdateMove<DIM,SURF,0> >(pstart,pstop),*this);
      }
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagUpdateMove<DIM,SURF,-1> >(pstart,pstop),*this,reduce);
    }
    copymode = 0;

    particle_kk->modify(Device,PARTICLE_MASK);
    d_particles = t_particle_1d(); // destroy reference to reduce memory use

    // END of pstart/pstop loop advecting all particles

    k_nmigrate.modify<DeviceType>();
    k_nmigrate.sync<SPAHostType>();
    nmigrate = h_nmigrate();

    DAT::t_int_1d d_mlist_small = Kokkos::subview(k_mlist.d_view,std::make_pair(0,nmigrate));
    HAT::t_int_1d h_mlist_small = HAT::t_int_1d(Kokkos::view_alloc("mlist_mirror",Kokkos::WithoutInitializing),nmigrate);
    Kokkos::deep_copy(h_mlist_small, d_mlist_small);
    auto mlist_small = h_mlist_small.data();

    int error_flag;

    if (sparta->kokkos->atomic_reduction) {
      k_ntouch_one.modify<DeviceType>();
      k_ntouch_one.sync<SPAHostType>();
      ntouch_one = h_ntouch_one();

      k_nexit_one.modify<DeviceType>();
      k_nexit_one.sync<SPAHostType>();
      nexit_one = h_nexit_one();

      k_nboundary_one.modify<DeviceType>();
      k_nboundary_one.sync<SPAHostType>();
      nboundary_one = h_nboundary_one();

      k_ncomm_one.modify<DeviceType>();
      k_ncomm_one.sync<SPAHostType>();
      ncomm_one = h_ncomm_one();

      k_nscheck_one.modify<DeviceType>();
      k_nscheck_one.sync<SPAHostType>();
      nscheck_one = h_nscheck_one();

      k_nscollide_one.modify<DeviceType>();
      k_nscollide_one.sync<SPAHostType>();
      nscollide_one = h_nscollide_one();

      k_nreact_one.modify<DeviceType>();
      k_nreact_one.sync<SPAHostType>();
      surf->nreact_one = h_nreact_one();

      k_nstuck.modify<DeviceType>();
      k_nstuck.sync<SPAHostType>();
      nstuck = h_nstuck();
    } else {
      ntouch_one       += reduce.ntouch_one   ;
      nexit_one        += reduce.nexit_one    ;
      nboundary_one    += reduce.nboundary_one;
      ncomm_one        += reduce.ncomm_one    ;
      nscheck_one      += reduce.nscheck_one  ;
      nscollide_one    += reduce.nscollide_one;
      surf->nreact_one += reduce.nreact_one   ;
      nstuck           += reduce.nstuck       ;
    }

    k_entryexit.modify<DeviceType>();
    k_entryexit.sync<SPAHostType>();
    entryexit = h_entryexit();

    k_error_flag.modify<DeviceType>();
    k_error_flag.sync<SPAHostType>();
    error_flag = h_error_flag();

    if (error_flag) {
      char str[128];
      sprintf(str,
              "Particle being sent to self proc "
              "on step " BIGINT_FORMAT,
              update->ntimestep);
      error->one(FLERR,str);
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
      timer->stamp(TIME_MOVE);
      pstart = ((CommKokkos*)comm)->migrate_particles(nmigrate,mlist_small,d_mlist_small);
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
  nmove_running += nlocal;
  ntouch_running += ntouch_one;
  ncomm_running += ncomm_one;
  nboundary_running += nboundary_one;
  nexit_running += nexit_one;
  nscheck_running += nscheck_one;
  nscollide_running += nscollide_one;
  surf->nreact_running += surf->nreact_one;
}

/* ---------------------------------------------------------------------- */

template<int DIM, int SURF, int ATOMIC_REDUCTION>
KOKKOS_INLINE_FUNCTION
void UpdateKokkos::operator()(TagUpdateMove<DIM,SURF,ATOMIC_REDUCTION>, const int &i) const {
  UPDATE_REDUCE reduce;
  this->template operator()<DIM,SURF,ATOMIC_REDUCTION>(TagUpdateMove<DIM,SURF,ATOMIC_REDUCTION>(), i, reduce);
}

/*-----------------------------------------------------------------------------*/

template<int DIM, int SURF, int ATOMIC_REDUCTION>
KOKKOS_INLINE_FUNCTION
void UpdateKokkos::operator()(TagUpdateMove<DIM,SURF,ATOMIC_REDUCTION>, const int &i, UPDATE_REDUCE &reduce) const {
  if (d_error_flag()) return;

  // int m;
  bool hitflag;
  int icell,icell_original,outface,bflag,nflag,pflag,itmp;
  int side,minsurf,nsurf,cflag,isurf,exclude,stuck_iterate;
  double dtremain,frac,newfrac,param,minparam,rnew,dtsurf,tc,tmp;
  double xnew[3],xhold[3],xc[3],vc[3],minxc[3],minvc[3];
  double *x,*v;
  Surf::Tri *tri;
  Surf::Line *line;

  Particle::OnePart &particle_i = d_particles[i];
  pflag = particle_i.flag;

  // Particle::OnePart iorig;
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
    if (gravity_3d_flag) gravity3d(dtremain,xnew,v);
    else if (gravity_2d_flag) gravity2d(dtremain,xnew,v);
  } else if (pflag == PINSERT) {
    dtremain = particle_i.dtremain;
    xnew[0] = x[0] + dtremain*v[0];
    xnew[1] = x[1] + dtremain*v[1];
    if (DIM != 2) xnew[2] = x[2] + dtremain*v[2];
    if (gravity_3d_flag) gravity3d(dtremain,xnew,v);
    else if (gravity_2d_flag) gravity2d(dtremain,xnew,v);
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

  particle_i.flag = PKEEP;
  icell = particle_i.icell;
  double* lo = d_cells[icell].lo;
  double* hi = d_cells[icell].hi;
  cellint* neigh = d_cells[icell].neigh;
  int nmask = d_cells[icell].nmask;
  stuck_iterate = 0;
  if (ATOMIC_REDUCTION == 1)
    Kokkos::atomic_fetch_add(&d_ntouch_one(),1);
  else if (ATOMIC_REDUCTION == 0)
    d_ntouch_one()++;
  else
    reduce.ntouch_one++;

  // advect one particle from cell to cell and thru surf collides til done

  while (1) {

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
      if (pflag == PEXIT && x[1] == lo[1]) {
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

      if (pflag == PEXIT && x[1] == hi[1]) {
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

    /** Kokkos functions **/
    // slist_active[m]->surf_tally

    /** Kokkos views **/
    // slist_active


//    ////Need to error out for now if surface reactions create (or destroy?) particles////
//
//    if (nsurf_tally)
//      for (m = 0; m < nsurf_tally; m++)
//        slist_active[m]->surf_tally(minsurf,&iorig,ipart,jpart);////
//

    if (SURF) {

      nsurf = d_cells[icell].nsurf;
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
                                 d_pts[tri->p1].x,d_pts[tri->p2].x,
                                 d_pts[tri->p3].x,tri->norm,xc,param,side);
          }
          if (DIM == 2) {
            line = &d_lines[isurf];
            hitflag = GeometryKokkos::
              line_line_intersect(x,xnew,
                                  d_pts[line->p1].x,d_pts[line->p2].x,
                                  line->norm,xc,param,side);
          }
          if (DIM == 1) {
            line = &d_lines[isurf];
            hitflag = GeometryKokkos::
              axi_line_intersect(dtsurf,x,v,outface,lo,hi,
                                 d_pts[line->p1].x,d_pts[line->p2].x,
                                 line->norm,exclude == isurf,
                                 xc,vc,param,side);
          }

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


        if (ATOMIC_REDUCTION == 1)
          Kokkos::atomic_fetch_add(&d_nscheck_one(),nsurf);
        else if (ATOMIC_REDUCTION == 0)
          d_nscheck_one() += nsurf;
        else
          reduce.nscheck_one += nsurf;

        if (cflag) {
          // NOTE: this check is no longer needed?
          /**if (minside == INSIDE) {
            char str[128];
            sprintf(str,
                    "Particle %d on proc %d hit inside of "
                    "surf %d on step " BIGINT_FORMAT,
                    i,me,minsurf,update->ntimestep);
            error->one(FLERR,str);
          }**/

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

          //if (nsurf_tally) 
          //  memcpy(&iorig,&particle_i,sizeof(Particle::OnePart));
          int n = DIM == 3 ? tri->isc : line->isc;
          int sc_type = sc_type_list[n];
          int m = sc_map[n];

          if (DIM == 3)
            if (sc_type == 0)
              jpart = sc_kk_specular_copy[m].obj.
                collide_kokkos(ipart,tri->norm,dtremain,tri->isr);/////
            else if (sc_type == 1)
              jpart = sc_kk_diffuse_copy[m].obj.
                collide_kokkos(ipart,tri->norm,dtremain,tri->isr);/////
            else if (sc_type == 2)
              jpart = sc_kk_vanish_copy[m].obj.
                collide_kokkos(ipart,tri->norm,dtremain,tri->isr);/////
          if (DIM != 3)
            if (sc_type == 0)
              jpart = sc_kk_specular_copy[m].obj.
                collide_kokkos(ipart,line->norm,dtremain,line->isr);////
            else if (sc_type == 1)
              jpart = sc_kk_diffuse_copy[m].obj.
                collide_kokkos(ipart,line->norm,dtremain,line->isr);////
            else if (sc_type == 2)
              jpart = sc_kk_vanish_copy[m].obj.
                collide_kokkos(ipart,line->norm,dtremain,line->isr);////

          ////Need to error out for now if surface reactions create (or destroy?) particles////
          //if (jpart) {
          //  particles = particle->particles;
          //  x = particle_i.x;
          //  v = particle_i.v;
          //  jpart->flag = PSURF + 1 + minsurf;
          //  jpart->dtremain = dtremain;
          //  jpart->weight = particle_i.weight;
          //  pstop++;
          //}

          //if (nsurf_tally)
          //  for (m = 0; m < nsurf_tally; m++)
          //    slist_active[m]->surf_tally(minsurf,&iorig,ipart,jpart);////
          
          // nstuck = consective iterations particle is immobile

          if (minparam == 0.0) stuck_iterate++;
          else stuck_iterate = 0;

          // reset post-bounce xnew

          xnew[0] = x[0] + dtremain*v[0];
          xnew[1] = x[1] + dtremain*v[1];
          if (DIM != 2) xnew[2] = x[2] + dtremain*v[2];

          exclude = minsurf;
          if (ATOMIC_REDUCTION == 1)
            Kokkos::atomic_fetch_add(&d_nscollide_one(),1);
          else if (ATOMIC_REDUCTION == 0)
            d_nscollide_one()++;
          else
            reduce.nscollide_one++;

          // if ipart = NULL, particle discarded due to surface chem
          // else if particle not stuck, continue advection while loop
          // if stuck, mark for DISCARD, and drop out of SURF code

          if (ipart == NULL) particle_i.flag = PDISCARD;
          else if (stuck_iterate < MAXSTUCK) continue;
          else {
            particle_i.flag = PDISCARD;
            if (ATOMIC_REDUCTION == 1)
              Kokkos::atomic_fetch_add(&d_nstuck(),1);
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
    // if migrating to another proc,
    //   flag as PDONE so new proc won't move it more on this step
    
    if (outface == INTERIOR) {
      if (DIM == 1) axi_remap(xnew,v);
      x[0] = xnew[0];
      x[1] = xnew[1];
      if (DIM == 3) x[2] = xnew[2];
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
    //   result of id_find_child could be unknown:
    //     particle is hitting face of a ghost child cell which extends
    //     beyond my ghost halo, cell on other side of face is a parent,
    //     it's child which the particle is in is entirely beyond my halo
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
      icell = grid_kk_copy.obj.id_find_child(neigh[outface],x);
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
    // other = reflected particle stays in same grid cell

    else {
      ipart = &particle_i;

      Particle::OnePart iorig;
      if (nboundary_tally) 
        memcpy(&iorig,&particle_i,sizeof(Particle::OnePart));

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
            collide_kokkos(ipart,domain_kk_copy.obj.norm[outface],dtremain,domain_kk_copy.obj.surf_react[outface]);/////
        else if (sc_type == 1)
          jpart = sc_kk_diffuse_copy[m].obj.
            collide_kokkos(ipart,domain_kk_copy.obj.norm[outface],dtremain,domain_kk_copy.obj.surf_react[outface]);/////
        else if (sc_type == 2)
          jpart = sc_kk_vanish_copy[m].obj.
            collide_kokkos(ipart,domain_kk_copy.obj.norm[outface],dtremain,domain_kk_copy.obj.surf_react[outface]);/////

        if (ipart) {
          double *x = ipart->x;
          double *v = ipart->v;
          xnew[0] = x[0] + dtremain*v[0];
          xnew[1] = x[1] + dtremain*v[1];
          if (domain_kk_copy.obj.dimension == 3) xnew[2] = x[2] + dtremain*v[2];
        }
        bflag = SURFACE;
      } else {
        bflag = domain_kk_copy.obj.collide_kokkos(ipart,outface,lo,hi,xnew/*,dtremain*/);
      }

      //if (jpart) {
      //  particles = particle->particles;
      //  x = particle_i.x;
      //  v = particle_i.v;
      //}

      if (nboundary_tally)
        for (int m = 0; m < nboundary_tally; m++)
          blist_active_copy[m].obj.
            boundary_tally_kk(outface,bflag,&iorig,ipart,jpart,domain_kk_copy.obj.norm[outface]);

      if (DIM == 1) {
        xnew[0] = x[0] + dtremain*v[0];
        xnew[1] = x[1] + dtremain*v[1];
        xnew[2] = x[2] + dtremain*v[2];
      }

      if (bflag == OUTFLOW) {
        particle_i.flag = PDISCARD;
        if (ATOMIC_REDUCTION == 1)
          Kokkos::atomic_fetch_add(&d_nexit_one(),1);
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
          icell = grid_kk_copy.obj.id_find_child(neigh[outface],x);
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
          //jpart->flag = PSURF;
          //jpart->dtremain = dtremain;
          //jpart->weight = particle_i.weight;
          //pstop++;
        }
        Kokkos::atomic_fetch_add(&d_nboundary_one(),1);
        Kokkos::atomic_fetch_add(&d_ntouch_one(),-1);    // decrement here since will increment below
      } else {
        if (ATOMIC_REDUCTION == 1) {
          Kokkos::atomic_fetch_add(&d_nboundary_one(),1);
          Kokkos::atomic_fetch_add(&d_ntouch_one(),-1);    // decrement here since will increment below
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
      Kokkos::atomic_fetch_add(&d_ntouch_one(),1);
    else if (ATOMIC_REDUCTION == 0)
      d_ntouch_one()++;
    else
      reduce.ntouch_one++;
  }

  // END of while loop over advection of single particle

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
        Kokkos::atomic_fetch_add(&d_ncomm_one(),1);
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
  // not considered a collision if particles starts on surf, moving out
  // not considered a collision if 2 params are tied and one is INSIDE surf

  int nsurf = d_cells[icell].nsurf;
  int isplit = d_cells[icell].isplit;
  double *xnew = d_sinfo[isplit].xsplit;

  cflag = 0;
  minparam = 2.0;
  auto csplits_begin = d_csplits.row_map(isplit);
  auto csurfs_begin = d_csurfs.row_map(isplit);
  for (m = 0; m < nsurf; m++) {
    if (d_csplits.entries(csplits_begin + m) < 0) continue;
    isurf = d_csurfs.entries(csurfs_begin + m);
    tri = &d_tris[isurf];
    hitflag = GeometryKokkos::
      line_tri_intersect(x,xnew,
                         d_pts[tri->p1].x,d_pts[tri->p2].x,d_pts[tri->p3].x,
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
                          d_pts[line->p1].x,d_pts[line->p2].x,line->norm,
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

  //if (nslist_compute) {
  //  for (i = 0; i < nslist_compute; i++)
  //    if (slist_compute[i]->matchstep(ntimestep)) {
  //      slist_active[nsurf_tally++] = slist_compute[i];
  //    }
  //}

  if (nboundary_tally > KOKKOS_MAX_BLIST)
    error->all(FLERR,"Kokkos currently only supports two instances of compute boundary");

  if (nboundary_tally) {
    for (i = 0; i < nboundary_tally; i++) {
      ComputeBoundaryKokkos* compute_boundary_kk = (ComputeBoundaryKokkos*)(blist_active[i]);
      compute_boundary_kk->pre_boundary_tally();
      blist_active_copy[i].copy(compute_boundary_kk);
    }
  }
}

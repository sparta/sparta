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

#ifndef SPARTA_UPDATE_KOKKOS_H
#define SPARTA_UPDATE_KOKKOS_H

#include "update.h"
#include "kokkos_type.h"
#include "particle.h"
#include "grid_kokkos.h"
#include "domain_kokkos.h"
#include "kokkos_copy.h"
#include "surf_collide_diffuse_kokkos.h"
#include "surf_collide_specular_kokkos.h"
#include "surf_collide_vanish_kokkos.h"
#include "surf_collide_piston_kokkos.h"
#include "surf_collide_transparent_kokkos.h"
#include "surf_collide_adiabatic_kokkos.h"
#include "surf_collide_impulsive_kokkos.h"
#include "surf_collide_td_kokkos.h"
#include "surf_collide_cll_kokkos.h"
#include "compute_boundary_kokkos.h"
#include "compute_react_boundary_kokkos.h"
#include "compute_surf_kokkos.h"
#include "compute_surf_collision_tally_kokkos.h"
#include "compute_surf_reaction_tally_kokkos.h"
#include "compute_isurf_grid_kokkos.h"
#include "compute_react_isurf_grid_kokkos.h"
#include "compute_react_surf_kokkos.h"

namespace SPARTA_NS {

// surf_collide style tags, used to dispatch on device where the host's
//   virtual SurfCollide::collide() is not available

enum{SC_SPECULAR,SC_DIFFUSE,SC_VANISH,SC_PISTON,SC_TRANSPARENT,
     SC_ADIABATIC,SC_IMPULSIVE,SC_TD,SC_CLL,SC_NSTYLE};

// host-side phases every surf_collide model is run through

enum{SC_PRE,SC_POST,SC_BACKUP,SC_RESTORE};

struct s_UPDATE_REDUCE {
  // per-step counters are bigint since they can exceed 2^31
  //   in one step at large per-proc particle counts
  bigint ntouch_one,nexit_one,nboundary_one,ncomm_one,
         nscheck_one,nscollide_one,nreact_one;
  int entryexit,nstuck,naxibad,error_flag;
  KOKKOS_INLINE_FUNCTION
  s_UPDATE_REDUCE() {
    ntouch_one    = 0;
    nexit_one     = 0;
    nboundary_one = 0;
    ncomm_one     = 0;
    nscheck_one   = 0;
    nscollide_one = 0;
    nreact_one    = 0;
    nstuck        = 0;
    naxibad       = 0;
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(const s_UPDATE_REDUCE &rhs) {
    ntouch_one    += rhs.ntouch_one   ;
    nexit_one     += rhs.nexit_one    ;
    nboundary_one += rhs.nboundary_one;
    ncomm_one     += rhs.ncomm_one    ;
    nscheck_one   += rhs.nscheck_one  ;
    nscollide_one += rhs.nscollide_one;
    nreact_one    += rhs.nreact_one   ;
    nstuck        += rhs.nstuck       ;
    naxibad       += rhs.naxibad      ;
  }
};
typedef struct s_UPDATE_REDUCE UPDATE_REDUCE;

template<int DIM, int SURF, int REACT, int OPT, int ATOMIC_REDUCTION>
struct TagUpdateMove{};

class UpdateKokkos : public Update {
 public:
  typedef UPDATE_REDUCE value_type;

  DAT::tdual_int_1d k_mlist;
  DAT::tdual_int_1d k_mlist_small;
  //DAT::t_int_1d d_mlist_small;
  //HAT::t_int_scalar h_mlist_small;
  //int* mlist_small;

  UpdateKokkos(class SPARTA *);
  ~UpdateKokkos();
  void init();
  void setup();
  void run(int);

  template<int DIM, int SURF, int REACT, int OPT, int ATOMIC_REDUCTION>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagUpdateMove<DIM,SURF,REACT,OPT,ATOMIC_REDUCTION>, const int&) const;

  template<int DIM, int SURF, int REACT, int OPT, int ATOMIC_REDUCTION>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagUpdateMove<DIM,SURF,REACT,OPT,ATOMIC_REDUCTION>, const int&, UPDATE_REDUCE&) const;

 private:

  double dt;
  int field_active[3];

  // data for optimized particle moves

  double dx,dy,dz,Lx,Ly,Lz;
  double xlo,ylo,zlo,xhi,yhi,zhi;
  int ncx,ncy,ncz;

  // what the fast path may do for itself on each global boundary face, indexed
  // XLO..ZHI as domain->bflag is.  0 = hand the particle to the standard move,
  // 1 = periodic, translate by the box length, 2 = specular, mirror about the
  // face and negate the normal velocity component.
  // set in move() rather than init() because it depends on nboundary_tally,
  // which is per-step

  int bcopt[6];

  GridKokkos::hash_type hash_kk;

  // dense cell lookup, see GridKokkos::update_halo_index().  extent 0 when
  // unavailable, in which case the fast path uses hash_kk

  DAT::t_int_1d d_halo_index;
  int halo_ilo,halo_jlo,halo_klo;
  int halo_nx,halo_ny,halo_nz;

  // retake hash_kk and d_halo_index from the grid, see its definition

  void grid_index_refresh();

  t_cell_1d d_cells;
  t_sinfo_1d d_sinfo;
  t_pcell_1d d_pcells;

  Kokkos::Crs<int, DeviceType, void, crs_size_type> d_csurfs;
  Kokkos::Crs<int, DeviceType, void, crs_size_type> d_csplits;
  Kokkos::Crs<int, DeviceType, void, crs_size_type> d_csubs;

  t_line_1d d_lines;
  t_tri_1d d_tris;

  t_particle_1d d_particles;

  DAT::t_float_2d_lr d_fieldfix_array_particle;
  DAT::t_float_2d_lr d_fieldfix_array_grid;

  class KokkosBase* KKBaseFieldFix;

  KKCopy<GridKokkos> grid_kk_copy;
  KKCopy<DomainKokkos> domain_kk_copy;

  // surf_collide models used to sit in fixed-size KKCopy arrays here, nine
  //   styles x two instances.  This class is itself the functor handed by
  //   value to every move kernel, and each model nests its own surf_react
  //   copies, so those arrays made sizeof(UpdateKokkos) 224 KB -- copied to
  //   the device on every launch -- while capping a run at two instances of
  //   each style, which an ordinary model with three wall temperatures hits.
  // The models now live in device memory instead, one buffer per style sized
  //   at run time, and the functor carries only the buffers and two index
  //   maps.  The bytes are blitted in rather than constructed there, which is
  //   what KKCopy::copy() already does (see kokkos_copy.h) and is sound for
  //   the same reason: on device the models are only read, through
  //   KOKKOS_INLINE_FUNCTION members, so the vtable pointer is never used and
  //   the View handles they carry are kept alive by the originals in surf->sc.

  DAT::t_int_1d d_sc_type;              // surf_collide index -> style tag
  DAT::t_int_1d d_sc_map;               // surf_collide index -> slot in style
  DAT::tdual_char_1d k_sc[SC_NSTYLE];   // blitted models, one buffer per style
  DAT::t_char_1d d_sc[SC_NSTYLE];

  int nsc_style[SC_NSTYLE];             // # of instances of each style

  // the surf_collide index maps depend only on surf->sc[n]->style, which is
  //   fixed for a run (surf_collide is a between-runs command), but
  //   setup_surf_collide_models() runs once per migration iteration.  Build
  //   them once per run and keep the host mirrors, instead of allocating a
  //   fresh mirror and re-uploading both maps every iteration.
  //   init() invalidates, so a new run picks up any style change

  int nsc_index_cached;                 // surf->nsc the maps were built for
  DAT::t_int_1d::host_mirror_type h_sc_type;
  DAT::t_int_1d::host_mirror_type h_sc_map;

  static int surf_collide_style_tag(class SurfCollide *);
  static size_t sc_sizeof(int);
  void sc_phase(class SurfCollide *, int);
  void setup_surf_collide_models();     // count, blit and upload, once per move
  void upload_surf_collide_models();    // re-blit after backup()/restore()

  // dispatch a surface collision to the model at surf_collide index n
  //   the nine-way switch is what the host's virtual call becomes on device

  template<int REACT, int ATOMIC_REDUCTION>
  KOKKOS_INLINE_FUNCTION
  Particle::OnePart* surf_collide_dispatch(const int n, Particle::OnePart *&ip,
                                           double &dtremain, const int isurf,
                                           const double *norm, const int isr,
                                           int &reaction,
                                           const DAT::t_int_scalar &d_retry,
                                           const DAT::t_int_scalar &d_nlocal) const
  {
    const int m = d_sc_map[n];

#define SC_CASE(TAG,TYPE)                                               \
    case TAG:                                                           \
      return ((const TYPE *) d_sc[TAG].data())[m].                      \
        template collide_kokkos<REACT,ATOMIC_REDUCTION>                 \
          (ip,dtremain,isurf,norm,isr,reaction,d_retry,d_nlocal);

    switch (d_sc_type[n]) {
      SC_CASE(SC_SPECULAR,SurfCollideSpecularKokkos)
      SC_CASE(SC_DIFFUSE,SurfCollideDiffuseKokkos)
      SC_CASE(SC_VANISH,SurfCollideVanishKokkos)
      SC_CASE(SC_PISTON,SurfCollidePistonKokkos)
      SC_CASE(SC_TRANSPARENT,SurfCollideTransparentKokkos)
      SC_CASE(SC_ADIABATIC,SurfCollideAdiabaticKokkos)
      SC_CASE(SC_IMPULSIVE,SurfCollideImpulsiveKokkos)
      SC_CASE(SC_TD,SurfCollideTDKokkos)
      SC_CASE(SC_CLL,SurfCollideCLLKokkos)
    }

#undef SC_CASE

    return NULL;
  }

  // the active tally computes, partitioned by type.  Two representations,
  //   selected by SPARTA_KOKKOS_FIXED_LISTS (see kokkos_type.h):
  //   - default: one runtime-sized device buffer per type, objects blitted in
  //     exactly as KKCopy::copy() does (kokkos_copy.h:71).  No instance cap.
  //   - SPARTA_KOKKOS_FIXED_LISTS: the original fixed-size KKCopy arrays, held
  //     by value in this functor, capped at two instances of each type.
  //   The dispatch sites in the move kernel are written once, against the
  //   UK_* accessors below, so only these declarations and the setup routine
  //   differ between the two.

#ifdef SPARTA_KOKKOS_FIXED_LISTS
  KKCopy<ComputeSurfKokkos> slist_active_copy[KOKKOS_MAX_SLIST];
  KKCopy<ComputeISurfGridKokkos> slist_active_isurf_copy[KOKKOS_MAX_SLIST];
  KKCopy<ComputeSurfCollisionTallyKokkos> slist_active_coll_tally_copy[KOKKOS_MAX_SLIST];
  KKCopy<ComputeSurfReactionTallyKokkos> slist_active_react_tally_copy[KOKKOS_MAX_SLIST];
  KKCopy<ComputeReactISurfGridKokkos> slist_active_react_isurf_copy[KOKKOS_MAX_SLIST];
  KKCopy<ComputeReactSurfKokkos> slist_active_react_surf_copy[KOKKOS_MAX_SLIST];
  KKCopy<ComputeBoundaryKokkos> blist_active_copy[KOKKOS_MAX_BLIST];
  KKCopy<ComputeReactBoundaryKokkos> blist_active_react_copy[KOKKOS_MAX_BLIST];
  ComputeReactBoundaryKokkos tmp_compute_react_boundary_kk;

  // unused fixed slots must not alias a compute that may be reallocated or
  //   deleted while they still reference count it

  ComputeBoundaryKokkos tmp_compute_boundary_kk;
  ComputeSurfKokkos tmp_compute_surf_kk;
  ComputeISurfGridKokkos tmp_compute_isurf_grid_kk;
  ComputeReactISurfGridKokkos tmp_compute_react_isurf_grid_kk;
  ComputeReactSurfKokkos tmp_compute_react_surf_kk;

#define UK_SLIST_SURF(m)        slist_active_copy[m].obj
#define UK_SLIST_ISURF(m)       slist_active_isurf_copy[m].obj
#define UK_SLIST_COLL_TALLY(m)  slist_active_coll_tally_copy[m].obj
#define UK_SLIST_REACT_TALLY(m) slist_active_react_tally_copy[m].obj
#define UK_SLIST_REACT_ISURF(m) slist_active_react_isurf_copy[m].obj
#define UK_SLIST_REACT_SURF(m)  slist_active_react_surf_copy[m].obj
#define UK_BLIST(m)             blist_active_copy[m].obj
#define UK_BLIST_REACT(m)       blist_active_react_copy[m].obj

#else
  DAT::tdual_char_1d k_slist_surf, k_slist_isurf, k_slist_coll_tally,
                     k_slist_react_tally, k_slist_react_isurf,
                     k_slist_react_surf, k_blist, k_blist_react;
  DAT::t_char_1d d_slist_surf, d_slist_isurf, d_slist_coll_tally,
                 d_slist_react_tally, d_slist_react_isurf,
                 d_slist_react_surf, d_blist, d_blist_react;

#define UK_SLIST_SURF(m)        ((const ComputeSurfKokkos *) d_slist_surf.data())[m]
#define UK_SLIST_ISURF(m)       ((const ComputeISurfGridKokkos *) d_slist_isurf.data())[m]
#define UK_SLIST_COLL_TALLY(m)  ((const ComputeSurfCollisionTallyKokkos *) d_slist_coll_tally.data())[m]
#define UK_SLIST_REACT_TALLY(m) ((const ComputeSurfReactionTallyKokkos *) d_slist_react_tally.data())[m]
#define UK_SLIST_REACT_ISURF(m) ((const ComputeReactISurfGridKokkos *) d_slist_react_isurf.data())[m]
#define UK_SLIST_REACT_SURF(m)  ((const ComputeReactSurfKokkos *) d_slist_react_surf.data())[m]
#define UK_BLIST(m)             ((const ComputeBoundaryKokkos *) d_blist.data())[m]
#define UK_BLIST_REACT(m)       ((const ComputeReactBoundaryKokkos *) d_blist_react.data())[m]
#endif

  // partition of slist_active (set in tally_set):
  //   nslist_surf        = # of compute surf style tallies (slist_active_copy)
  //   nslist_isurf       = # of compute isurf/grid tallies (slist_active_isurf_copy)
  //   nslist_react_isurf = # of compute react/isurf/grid tallies
  // nslist_surf + nslist_isurf + nslist_react_isurf == nsurf_tally

  int nslist_surf,nslist_isurf,nslist_react_isurf,nslist_react_surf;
  int nslist_coll_tally,nslist_react_tally;
  int nblist_boundary,nblist_react;

  // grow every per-event tally compute after an overflowed attempt
  void grow_tally_computes();
  void rewind_tally_computes(int);



  // int scalars = flags and view-index counters, must stay int
  // bigint scalars = per-step statistics counters, can exceed 2^31
  //   in one step at large per-proc particle counts

  typedef Kokkos::DualView<int[8], DeviceType::array_layout, DeviceType> tdual_int_7;
  typedef tdual_int_7::t_dev t_int_7;
  typedef tdual_int_7::t_host t_host_int_7;
  t_int_7 d_scalars;
  t_host_int_7 h_scalars;

  typedef Kokkos::DualView<bigint[7], DeviceType::array_layout, DeviceType> tdual_bigint_7;
  typedef tdual_bigint_7::t_dev t_bigint_7;
  typedef tdual_bigint_7::t_host t_host_bigint_7;
  t_bigint_7 d_scalars_big;
  t_host_bigint_7 h_scalars_big;

  DAT::t_bigint_scalar d_ntouch_one;
  HAT::t_bigint_scalar h_ntouch_one;

  DAT::t_bigint_scalar d_nexit_one;
  HAT::t_bigint_scalar h_nexit_one;

  DAT::t_bigint_scalar d_nboundary_one;
  HAT::t_bigint_scalar h_nboundary_one;

  DAT::t_int_scalar d_nmigrate;
  HAT::t_int_scalar h_nmigrate;

  DAT::t_int_scalar d_entryexit;
  HAT::t_int_scalar h_entryexit;

  DAT::t_bigint_scalar d_ncomm_one;
  HAT::t_bigint_scalar h_ncomm_one;

  DAT::t_bigint_scalar d_nscheck_one;
  HAT::t_bigint_scalar h_nscheck_one;

  DAT::t_bigint_scalar d_nscollide_one;
  HAT::t_bigint_scalar h_nscollide_one;

  DAT::t_bigint_scalar d_nreact_one;
  HAT::t_bigint_scalar h_nreact_one;

  DAT::t_int_scalar d_nstuck;
  HAT::t_int_scalar h_nstuck;

  DAT::t_int_scalar d_naxibad;
  HAT::t_int_scalar h_naxibad;

  DAT::t_int_scalar d_error_flag;
  HAT::t_int_scalar h_error_flag;

  DAT::t_int_scalar d_retry;
  DAT::t_int_scalar d_tally_overflow;
  HAT::t_int_scalar h_retry;
  HAT::t_int_scalar h_tally_overflow;

  DAT::t_int_scalar d_nlocal;
  HAT::t_int_scalar h_nlocal;

  // per-face count of the {s} face mirrors the fast path did, see
  // Update::optmove_surf_tally().  kept out of d_scalars because it is only
  // touched when a face qualifies, which most runs have none of

  DAT::t_bigint_1d d_bcmirror;
  HAT::t_bigint_1d h_bcmirror;

  void backup();
  void free_particle_backup();
  void restore();
  t_particle_1d d_particles_backup;

  void tally_set(bigint);
  void setup_surf_tally_copies();

  // remap x and v components into axisymmetric plane
  // input x at end of linear move (x = xold + dt*v)
  // change x[1] = sqrt(x[1]^2 + x[2]^2), x[2] = 0.0
  // change vy,vz by rotation into axisymmetric plane

  KOKKOS_INLINE_FUNCTION
  void axi_remap(double *x, double *v) const {
    double ynew = x[1];
    double znew = x[2];
    x[1] = sqrt(ynew*ynew + znew*znew);
    x[2] = 0.0;
    if (x[1] > 0.0) {
      double rn = ynew / x[1];
      double wn = znew / x[1];
      double vy = v[1];
      double vz = v[2];
      v[1] = vy*rn + vz*wn;
      v[2] = -vy*wn + vz*rn;
    }
  };

  typedef void (UpdateKokkos::*FnPtr)();
  FnPtr moveptr;             // ptr to move method
  template <int, int, int, int> void move();

  //typedef void (UpdateKokkos::*FnPtr2)(int, int, double, double *, double *) const;
  //FnPtr2 moveperturb;        // ptr to moveperturb method
  //
  //// variants of moveperturb method
  //// adjust end-of-move x,v due to perturbation on straight-line advection

  KOKKOS_INLINE_FUNCTION
  int split3d(int, double*) const;

  KOKKOS_INLINE_FUNCTION
  int split2d(int, double*) const;

  // the two steps of the optimized move, see their definitions and the OPT
  // block of the move kernel

  template < int DIM > KOKKOS_INLINE_FUNCTION
  int optmove_bc(const double*, double*, int&) const;

  template < int DIM > KOKKOS_INLINE_FUNCTION
  int optmove_cell(const double*) const;

  // variants of moveperturb method
  // adjust end-of-move x,v due to perturbation on straight-line advection

  KOKKOS_INLINE_FUNCTION
  void field2d(double dt, double *x, double *v) const {
    const double dtsq = 0.5*dt*dt;
    x[0] += dtsq*field[0];
    x[1] += dtsq*field[1];
    v[0] += dt*field[0];
    v[1] += dt*field[1];
  };

  KOKKOS_INLINE_FUNCTION
  void field3d(double dt, double *x, double *v) const {
    const double dtsq = 0.5*dt*dt;
    x[0] += dtsq*field[0];
    x[1] += dtsq*field[1];
    x[2] += dtsq*field[2];
    v[0] += dt*field[0];
    v[1] += dt*field[1];
    v[2] += dt*field[2];
  };

  /* ----------------------------------------------------------------------
     calculate motion perturbation for a single particle I
       due to external per particle field
     array in fix[ifieldfix] stores per particle perturbations for x and v
  ------------------------------------------------------------------------- */

  KOKKOS_INLINE_FUNCTION
  void field_per_particle(int i, int icell, double dt, double *x, double *v) const
  {
    const double dtsq = 0.5*dt*dt;
    auto &d_array = d_fieldfix_array_particle;

    int icol = 0;
    if (field_active[0]) {
      x[0] += dtsq*d_array(i,icol);
      v[0] += dt*d_array(i,icol);
      icol++;
    }
    if (field_active[1]) {
      x[1] += dtsq*d_array(i,icol);
      v[1] += dt*d_array(i,icol);
      icol++;
    }
    if (field_active[2]) {
      x[2] += dtsq*d_array(i,icol);
      v[2] += dt*d_array(i,icol);
      icol++;
    }
  };

  /* ----------------------------------------------------------------------
     calculate motion perturbation for a single particle I in grid cell Icell
       due to external per grid cell field
     array in fix[ifieldfix] stores per grid cell perturbations for x and v
  ------------------------------------------------------------------------- */

  KOKKOS_INLINE_FUNCTION
  void field_per_grid(int i, int icell, double dt, double *x, double *v) const
  {
    const double dtsq = 0.5*dt*dt;
    auto &d_array = d_fieldfix_array_grid;

    int icol = 0;
    if (field_active[0]) {
      x[0] += dtsq*d_array(icell,icol);
      v[0] += dt*d_array(icell,icol);
      icol++;
    }
    if (field_active[1]) {
      x[1] += dtsq*d_array(icell,icol);
      v[1] += dt*d_array(icell,icol);
      icol++;
    }
    if (field_active[2]) {
      x[2] += dtsq*d_array(icell,icol);
      v[2] += dt*d_array(icell,icol);
      icol++;
    }
  };
};

}

#endif

/* ERROR/WARNING messages:

*/

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
#include "compute_boundary_kokkos.h"
#include "compute_surf_kokkos.h"

namespace SPARTA_NS {

#define KOKKOS_MAX_SURF_COLL_PER_TYPE 2
#define KOKKOS_MAX_TOT_SURF_COLL 10
#define KOKKOS_MAX_BLIST 2
#define KOKKOS_MAX_SLIST 2

struct s_UPDATE_REDUCE {
  int ntouch_one,nexit_one,nboundary_one,
      entryexit,ncomm_one,
      nscheck_one,nscollide_one,nreact_one,nstuck,
      naxibad,error_flag;
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
  GridKokkos::hash_type hash_kk;

  t_cell_1d d_cells;
  t_sinfo_1d d_sinfo;
  t_pcell_1d d_pcells;

  Kokkos::Crs<int, DeviceType, void, int> d_csurfs;
  Kokkos::Crs<int, DeviceType, void, int> d_csplits;
  Kokkos::Crs<int, DeviceType, void, int> d_csubs;

  t_line_1d d_lines;
  t_tri_1d d_tris;

  t_particle_1d d_particles;

  DAT::t_float_2d_lr d_fieldfix_array_particle;
  DAT::t_float_2d_lr d_fieldfix_array_grid;

  class KokkosBase* KKBaseFieldFix;

  KKCopy<GridKokkos> grid_kk_copy;
  KKCopy<DomainKokkos> domain_kk_copy;

  int sc_type_list[KOKKOS_MAX_TOT_SURF_COLL];
  int sc_map[KOKKOS_MAX_TOT_SURF_COLL];
  KKCopy<SurfCollideSpecularKokkos> sc_kk_specular_copy[KOKKOS_MAX_SURF_COLL_PER_TYPE];
  KKCopy<SurfCollideDiffuseKokkos> sc_kk_diffuse_copy[KOKKOS_MAX_SURF_COLL_PER_TYPE];
  KKCopy<SurfCollideVanishKokkos> sc_kk_vanish_copy[KOKKOS_MAX_SURF_COLL_PER_TYPE];
  KKCopy<SurfCollidePistonKokkos> sc_kk_piston_copy[KOKKOS_MAX_SURF_COLL_PER_TYPE];
  KKCopy<SurfCollideTransparentKokkos> sc_kk_transparent_copy[KOKKOS_MAX_SURF_COLL_PER_TYPE];
  KKCopy<ComputeBoundaryKokkos> blist_active_copy[KOKKOS_MAX_BLIST];
  KKCopy<ComputeSurfKokkos> slist_active_copy[KOKKOS_MAX_SLIST];


  typedef Kokkos::DualView<int[14], DeviceType::array_layout, DeviceType> tdual_int_14;
  typedef tdual_int_14::t_dev t_int_14;
  typedef tdual_int_14::t_host t_host_int_14;
  t_int_14 d_scalars;
  t_host_int_14 h_scalars;

  DAT::t_int_scalar d_ntouch_one;
  HAT::t_int_scalar h_ntouch_one;

  DAT::t_int_scalar d_nexit_one;
  HAT::t_int_scalar h_nexit_one;

  DAT::t_int_scalar d_nboundary_one;
  HAT::t_int_scalar h_nboundary_one;

  DAT::t_int_scalar d_nmigrate;
  HAT::t_int_scalar h_nmigrate;

  DAT::t_int_scalar d_entryexit;
  HAT::t_int_scalar h_entryexit;

  DAT::t_int_scalar d_ncomm_one;
  HAT::t_int_scalar h_ncomm_one;

  DAT::t_int_scalar d_nscheck_one;
  HAT::t_int_scalar h_nscheck_one;

  DAT::t_int_scalar d_nscollide_one;
  HAT::t_int_scalar h_nscollide_one;

  DAT::t_int_scalar d_nreact_one;
  HAT::t_int_scalar h_nreact_one;

  DAT::t_int_scalar d_nstuck;
  HAT::t_int_scalar h_nstuck;

  DAT::t_int_scalar d_naxibad;
  HAT::t_int_scalar h_naxibad;

  DAT::t_int_scalar d_error_flag;
  HAT::t_int_scalar h_error_flag;

  DAT::t_int_scalar d_retry;
  HAT::t_int_scalar h_retry;

  DAT::t_int_scalar d_nlocal;
  HAT::t_int_scalar h_nlocal;

  void backup();
  void restore();
  t_particle_1d d_particles_backup;

  void bounce_set(bigint);

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
    double rn = ynew / x[1];
    double wn = znew / x[1];
    double vy = v[1];
    double vz = v[2];
    v[1] = vy*rn + vz*wn;
    v[2] = -vy*wn + vz*rn;
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

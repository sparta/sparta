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
#include "compute_boundary_kokkos.h"
#include "compute_surf_kokkos.h"

namespace SPARTA_NS {

#define KOKKOS_SURF_COLL_TYPES 4
#define KOKKOS_MAX_SURF_COLL_PER_TYPE 2
#define KOKKOS_TOT_SURF_COLL 6
#define KOKKOS_MAX_BLIST 2
#define KOKKOS_MAX_SLIST 2

struct s_UPDATE_REDUCE {
  int ntouch_one,nexit_one,nboundary_one,
      entryexit,ncomm_one,
      nscheck_one,nscollide_one,nreact_one,nstuck,error_flag;
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
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(const volatile s_UPDATE_REDUCE &rhs) volatile {
    ntouch_one    += rhs.ntouch_one   ;
    nexit_one     += rhs.nexit_one    ;
    nboundary_one += rhs.nboundary_one;
    ncomm_one     += rhs.ncomm_one    ;
    nscheck_one   += rhs.nscheck_one  ;
    nscollide_one += rhs.nscollide_one;
    nreact_one    += rhs.nreact_one   ;
    nstuck        += rhs.nstuck       ;
  }
};
typedef struct s_UPDATE_REDUCE UPDATE_REDUCE;

template<int DIM, int SURF, int ATOMIC_REDUCTION>
struct TagUpdateMove{};

class UpdateKokkos : public Update {
 public:
  typedef ArrayTypes<DeviceType> AT;
  typedef UPDATE_REDUCE value_type;

  DAT::tdual_int_1d k_mlist;
  //DAT::t_int_1d d_mlist_small;
  //HAT::t_int_scalar h_mlist_small;
  //int* mlist_small;

  UpdateKokkos(class SPARTA *);
  ~UpdateKokkos();
  void init();
  void setup();
  void run(int);

  template<int DIM, int SURF, int ATOMIC_REDUCTION>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagUpdateMove<DIM,SURF,ATOMIC_REDUCTION>, const int&) const;

  template<int DIM, int SURF, int ATOMIC_REDUCTION>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagUpdateMove<DIM,SURF,ATOMIC_REDUCTION>, const int&, UPDATE_REDUCE&) const;

 private:

  double dt;

  t_cell_1d d_cells;
  t_sinfo_1d d_sinfo;

  Kokkos::Crs<int, SPADeviceType, void, int> d_csurfs;
  Kokkos::Crs<int, SPADeviceType, void, int> d_csplits;
  Kokkos::Crs<int, SPADeviceType, void, int> d_csubs;

  t_line_1d d_lines;
  t_tri_1d d_tris;

  t_particle_1d d_particles;

  KKCopy<GridKokkos> grid_kk_copy;
  KKCopy<DomainKokkos> domain_kk_copy;

  int sc_type_list[KOKKOS_TOT_SURF_COLL];
  int sc_map[KOKKOS_TOT_SURF_COLL];
  KKCopy<SurfCollideSpecularKokkos> sc_kk_specular_copy[KOKKOS_MAX_SURF_COLL_PER_TYPE];
  KKCopy<SurfCollideDiffuseKokkos> sc_kk_diffuse_copy[KOKKOS_MAX_SURF_COLL_PER_TYPE];
  KKCopy<SurfCollideVanishKokkos> sc_kk_vanish_copy[KOKKOS_MAX_SURF_COLL_PER_TYPE];
  KKCopy<SurfCollidePistonKokkos> sc_kk_piston_copy[KOKKOS_MAX_SURF_COLL_PER_TYPE];
  KKCopy<ComputeBoundaryKokkos> blist_active_copy[KOKKOS_MAX_BLIST];
  KKCopy<ComputeSurfKokkos> slist_active_copy[KOKKOS_MAX_SLIST];


  typedef Kokkos::DualView<int[11], SPADeviceType::array_layout, SPADeviceType> tdual_int_11;
  typedef tdual_int_11::t_dev t_int_11;
  typedef tdual_int_11::t_host t_host_int_11;
  t_int_11 d_scalars;
  t_host_int_11 h_scalars;

  typename AT::t_int_scalar d_ntouch_one;
  HAT::t_int_scalar h_ntouch_one;

  typename AT::t_int_scalar d_nexit_one;
  HAT::t_int_scalar h_nexit_one;

  typename AT::t_int_scalar d_nboundary_one;
  HAT::t_int_scalar h_nboundary_one;

  typename AT::t_int_scalar d_nmigrate;
  HAT::t_int_scalar h_nmigrate;

  typename AT::t_int_scalar d_entryexit;
  HAT::t_int_scalar h_entryexit;

  typename AT::t_int_scalar d_ncomm_one;
  HAT::t_int_scalar h_ncomm_one;

  typename AT::t_int_scalar d_nscheck_one;
  HAT::t_int_scalar h_nscheck_one;

  typename AT::t_int_scalar d_nscollide_one;
  HAT::t_int_scalar h_nscollide_one;

  typename AT::t_int_scalar d_nreact_one;
  HAT::t_int_scalar h_nreact_one;

  typename AT::t_int_scalar d_nstuck;
  HAT::t_int_scalar h_nstuck;

  typename AT::t_int_scalar d_error_flag;
  HAT::t_int_scalar h_error_flag;

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
  template < int, int > void move();

  //
  //int perturbflag;
  typedef void (UpdateKokkos::*FnPtr2)(double, double *, double *) const;
  FnPtr2 moveperturb;        // ptr to moveperturb method
  //
  //// variants of moveperturb method
  //// adjust end-of-move x,v due to perturbation on straight-line advection

  KOKKOS_INLINE_FUNCTION
  int split3d(int, double*) const;

  KOKKOS_INLINE_FUNCTION
  int split2d(int, double*) const;

  int gravity_3d_flag,gravity_2d_flag;

  // variants of moveperturb method
  // adjust end-of-move x,v due to perturbation on straight-line advection

  KOKKOS_INLINE_FUNCTION
  void gravity2d(double dt, double *x, double *v) const {
    double dtsq = 0.5*dt*dt;
    x[0] += dtsq*gravity[0];
    x[1] += dtsq*gravity[1];
    v[0] += dt*gravity[0];
    v[1] += dt*gravity[1];
  };

  KOKKOS_INLINE_FUNCTION
  void gravity3d(double dt, double *x, double *v) const {
    double dtsq = 0.5*dt*dt;
    x[0] += dtsq*gravity[0];
    x[1] += dtsq*gravity[1];
    x[2] += dtsq*gravity[2];
    v[0] += dt*gravity[0];
    v[1] += dt*gravity[1];
    v[2] += dt*gravity[2];
  };
};

}

#endif

/* ERROR/WARNING messages:

*/

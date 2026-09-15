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

#ifdef REGION_CLASS

RegionStyle(sphere/kk,RegSphereKokkos)

#else

#ifndef SPARTA_REGION_SPHERE_KOKKOS_H
#define SPARTA_REGION_SPHERE_KOKKOS_H

#include "region_sphere.h"

#include "kokkos_base.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

struct TagRegSphereMatchAll{};

class RegSphereKokkos : public RegSphere, public KokkosBase {

 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  RegSphereKokkos(class SPARTA *, int, char **);
  RegSphereKokkos(class SPARTA *sparta);
  ~RegSphereKokkos() override;

  void match_all_kokkos(DAT::tdual_int_1d) override;

  // flatten to a single-token postfix stream; see region_prim_kokkos.h

  int flatten_region_kokkos(tdual_region_token_1d &k_tokens) override
  {
    region_token_grow(k_tokens,1);
    RegionTokenKK &t = k_tokens.view_host()[0];
    t.type = RKK_TOK_PRIM;
    RegionPrimKK &p = t.prim;
    p.style = RKK_SPHERE;
    p.interior = interior;
    p.axis = 0;
    p.a = p.b = p.c = p.d = p.e = p.f = 0.0;
    p.n0 = p.n1 = p.n2 = 0.0;
    p.a = xc; p.b = yc; p.c = zc; p.d = radius;
    k_tokens.modify_host();
    k_tokens.sync_device();
    return 1;
  }


  KOKKOS_INLINE_FUNCTION
  void operator()(TagRegSphereMatchAll, const int&) const;

  KOKKOS_INLINE_FUNCTION
  int match_kokkos(double x, double y, double z) const
  {
    return !(k_inside(x,y,z) ^ interior);
  }

 private:
  int groupbit;
  typename AT::t_int_1d d_match;
  t_particle_1d d_particles;

  KOKKOS_INLINE_FUNCTION
  int k_inside(double x, double y, double z) const
  {
    double delx = x - xc;
    double dely = y - yc;
    double delz = z - zc;
    double r = sqrt(delx*delx + dely*dely + delz*delz);

    if (r <= radius) return 1;
    return 0;
  }
};

}

#endif
#endif


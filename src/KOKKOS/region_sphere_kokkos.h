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


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

RegionStyle(cylinder/kk,RegCylinderKokkos)

#else

#ifndef SPARTA_REGION_CYLINDER_KOKKOS_H
#define SPARTA_REGION_CYLINDER_KOKKOS_H

#include "region_cylinder.h"

#include "kokkos_base.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

struct TagRegCylinderMatchAll{};

class RegCylinderKokkos : public RegCylinder, public KokkosBase {

 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  RegCylinderKokkos(class SPARTA *, int, char **);
  RegCylinderKokkos(class SPARTA *sparta);
  ~RegCylinderKokkos() override;

  void match_all_kokkos(DAT::tdual_int_1d) override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagRegCylinderMatchAll, const int&) const;

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
    double del1,del2,dist;
    int inside;

    if (axis == 'x') {
      del1 = y - c1;
      del2 = z - c2;
      dist = sqrt(del1*del1 + del2*del2);
      if (dist <= radius && x >= lo && x <= hi) inside = 1;
      else inside = 0;
    } else if (axis == 'y') {
      del1 = x - c1;
      del2 = z - c2;
      dist = sqrt(del1*del1 + del2*del2);
      if (dist <= radius && y >= lo && y <= hi) inside = 1;
      else inside = 0;
    } else {
      del1 = x - c1;
      del2 = y - c2;
      dist = sqrt(del1*del1 + del2*del2);
      if (dist <= radius && z >= lo && z <= hi) inside = 1;
      else inside = 0;
    }

    return inside;
  }
};

}

#endif
#endif


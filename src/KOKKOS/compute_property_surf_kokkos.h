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

#ifdef COMPUTE_CLASS

ComputeStyle(property/surf/kk,ComputePropertySurfKokkos)

#else

#ifndef SPARTA_COMPUTE_PROPERTY_SURF_KOKKOS_H
#define SPARTA_COMPUTE_PROPERTY_SURF_KOKKOS_H

#include "compute_property_surf.h"
#include "kokkos_base.h"
#include "kokkos_type.h"
#include "math_extra_kokkos.h"

namespace SPARTA_NS {

class ComputePropertySurfKokkos : public ComputePropertySurf, public KokkosBase {
 public:
  enum{ID,V1X,V1Y,V1Z,V2X,V2Y,V2Z,V3X,V3Y,V3Z,XC,YC,ZC,AREA,NORMX,NORMY,NORMZ};

  ComputePropertySurfKokkos(class SPARTA *, int, char **);
  ~ComputePropertySurfKokkos();
  void init();
  void compute_per_surf();
  void compute_per_surf_kokkos();

  KOKKOS_INLINE_FUNCTION
  double pack_one(int m, int field) const
  {
    const double THIRD = 1.0/3.0;
    if (dim == 2) {
      const auto &L = d_lines[m];
      switch (field) {
      case ID:    return (double) L.id;
      case V1X:   return L.p1[0];
      case V1Y:   return L.p1[1];
      case V2X:   return L.p2[0];
      case V2Y:   return L.p2[1];
      case XC:    return 0.5*(L.p1[0]+L.p2[0]);
      case YC:    return 0.5*(L.p1[1]+L.p2[1]);
      case AREA:  { double p12[3]; MathExtraKokkos::sub3(L.p2,L.p1,p12);
                    return MathExtraKokkos::len3(p12); }
      case NORMX: return L.norm[0];
      case NORMY: return L.norm[1];
      }
    } else {
      const auto &T = d_tris[m];
      switch (field) {
      case ID:    return (double) T.id;
      case V1X:   return T.p1[0];
      case V1Y:   return T.p1[1];
      case V1Z:   return T.p1[2];
      case V2X:   return T.p2[0];
      case V2Y:   return T.p2[1];
      case V2Z:   return T.p2[2];
      case V3X:   return T.p3[0];
      case V3Y:   return T.p3[1];
      case V3Z:   return T.p3[2];
      case XC:    return THIRD*(T.p1[0]+T.p2[0]+T.p3[0]);
      case YC:    return THIRD*(T.p1[1]+T.p2[1]+T.p3[1]);
      case ZC:    return THIRD*(T.p1[2]+T.p2[2]+T.p3[2]);
      case AREA:  { double p12[3],p13[3],cross[3];
                    MathExtraKokkos::sub3(T.p2,T.p1,p12);
                    MathExtraKokkos::sub3(T.p3,T.p1,p13);
                    MathExtraKokkos::cross3(p12,p13,cross);
                    return 0.5*MathExtraKokkos::len3(cross); }
      case NORMX: return T.norm[0];
      case NORMY: return T.norm[1];
      case NORMZ: return T.norm[2];
      }
    }
    return 0.0;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int &i) const
  {
    int m = d_cglobal[i];
    if (nvalues == 1) d_vector_surf[i] = pack_one(m,d_index[0]);
    else
      for (int n = 0; n < nvalues; n++)
        d_array_surf(i,n) = pack_one(m,d_index[n]);
  }

  DAT::tdual_float_1d k_vector_surf;
  DAT::tdual_float_2d_lr k_array_surf;

 private:
  int dim;
  DAT::t_int_1d d_index;
  DAT::t_int_1d d_cglobal;
  DAT::t_float_1d d_vector_surf;
  DAT::t_float_2d_lr d_array_surf;
  t_line_1d d_lines;
  t_tri_1d d_tris;
};

}

#endif
#endif

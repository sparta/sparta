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

ComputeStyle(surf/kk,ComputeSurfKokkos)

#else

#ifndef SPARTA_COMPUTE_SURF_KOKKOS_H
#define SPARTA_COMPUTE_SURF_KOKKOS_H

#include "compute_surf.h"
#include "kokkos_type.h"
#include "math_extra_kokkos.h"
#include "kokkos_copy.h"
#include "surf_react_global_kokkos.h"
#include "surf_react_prob_kokkos.h"

namespace SPARTA_NS {

struct TagComputeSurf_clear{};

class ComputeSurfKokkos : public ComputeSurf {
 public:
  ComputeSurfKokkos(class SPARTA *, int, char **);
  ComputeSurfKokkos(class SPARTA *);
  ~ComputeSurfKokkos();
  void init();
  void init_normflux();
  void clear();
  int tallyinfo(surfint *&);
  void update_hash();
  void pre_surf_tally();
  void post_surf_tally();

enum{NUM,NUMWT,NFLUX,NFLUXIN,MFLUX,MFLUXIN,FX,FY,FZ,TX,TY,TZ,
  PRESS,XPRESS,YPRESS,ZPRESS,XSHEAR,YSHEAR,ZSHEAR,KE,EROT,EVIB,ECHEM,ETOT};

/* ----------------------------------------------------------------------
   tally values for a single particle in icell
     colliding with surface element isurf, performing reaction (1 to N)
   iorig = particle ip before collision
   ip,jp = particles after collision
   ip = NULL means no particles after collision
   jp = NULL means one particle after collision
   jp != NULL means two particles after collision
------------------------------------------------------------------------- */

template <int ATOMIC_REDUCTION>
KOKKOS_INLINE_FUNCTION
void surf_tally_kk(double /*dtremain*/, int isurf, int icell, int reaction,
                   Particle::OnePart *iorig,
                   Particle::OnePart *ip, Particle::OnePart *jp) const
{
  // skip if no original particle and a reaction is taking place
  //   called by SurfReactAdsorb for on-surf reaction
  // FixEmitSurf also calls with no original particle but no reaction

  if (!iorig && reaction) return;

  // skip if isurf not in surface group

  if (dim == 2) {
    if (!(d_lines(isurf).mask & groupbit)) return;
  } else {
    if (!(d_tris(isurf).mask & groupbit)) return;
  }

  // skip if colliding/emitting species not in mixture group

  int origspecies = -1;
  int igroup;
  if (iorig) {
    origspecies = iorig->ispecies;
    igroup = d_s2g(imix,origspecies);
    if (igroup < 0) return;
  } else {
    igroup = d_s2g(imix,ip->ispecies);
    if (igroup < 0) return;
  }

  // itally = tally index of isurf
  // grow tally list if needed

  int itally,transparent,isr;

  surfint surfID;
  if (dim == 2) {
    surfID = d_lines[isurf].id;
    transparent = d_lines[isurf].transparent;
    isr = d_lines[isurf].isr;
  } else {
    surfID = d_tris[isurf].id;
    transparent = d_tris[isurf].transparent;
    isr = d_tris[isurf].isr;
  }

  // thread-safe, tally array will be compressed later

  itally = isurf;
  d_tally2surf(itally) = surfID;
  d_surf2tally(isurf) = isurf;

  double fluxscale = d_normflux(isurf);

  // tally all values associated with group into array
  // set fflag after force computation is done once
  // set tqflag after torque computation is done once
  // set nflag and tflag after normal and tangent computation is done once
  // particle weight used for all keywords except NUM
  // forcescale factor applied for keywords FX,FY,FZ
  // fluxscale factor applied for all keywords except NUM,FX,FY,FZ
  // if surf is transparent, all flux tallying is for incident particle only

  double vsqpre,ivsqpost,jvsqpost;
  double ierot,jerot,ievib,jevib,iother,jother,otherpre,etot;
  double pdelta[3],pnorm[3],ptang[3],pdelta_force[3],rdelta[3],torque[3];
  double *xcollide;

  double *norm;
  if (dim == 2) norm = d_lines(isurf).norm;
  else norm = d_tris(isurf).norm;

  double weight = 1.0;
  double origmass = 0.0;
  double imass,jmass;
  if (weightflag && iorig) weight = iorig->weight;
  else if (weightflag) weight = ip->weight;
  if (origspecies >= 0) origmass = d_species[origspecies].mass * weight;
  if (ip) imass = d_species(ip->ispecies).mass * weight;
  if (jp) jmass = d_species(jp->ispecies).mass * weight;

  double *vorig = NULL;
  double oerot,oevib;
  if (iorig) {
    vorig = iorig->v;
    oerot = iorig->erot;
    oevib = iorig->evib;
  } else {
    oerot = 0.0;
    oevib = 0.0;
  }

  int k = igroup*nvalue;

  int fflag = 0;
  int tqflag = 0;
  int nflag = 0;
  int tflag = 0;

  auto v_array_surf_tally = ScatterViewHelper<typename NeedDup<ATOMIC_REDUCTION,DeviceType>::value,decltype(dup_array_surf_tally),decltype(ndup_array_surf_tally)>::get(dup_array_surf_tally,ndup_array_surf_tally);
  auto a_array_surf_tally = v_array_surf_tally.template access<typename AtomicDup<ATOMIC_REDUCTION,DeviceType>::value>();

  for (int m = 0; m < nvalue; m++) {

    switch (d_which(m)) {

    // counts and fluxes

    case NUM:
      a_array_surf_tally(itally,k++) += 1.0;
      break;
    case NUMWT:
      a_array_surf_tally(itally,k++) += weight;
      break;
    case NFLUX:
      if (iorig) a_array_surf_tally(itally,k) += weight * fluxscale;
      if (!transparent) {
        if (ip) a_array_surf_tally(itally,k) -= weight * fluxscale;
        if (jp) a_array_surf_tally(itally,k) -= weight * fluxscale;
      }
      k++;
      break;
    case NFLUXIN:
      a_array_surf_tally(itally,k) += weight * fluxscale;
      k++;
      break;
    case MFLUX:
      if (iorig) a_array_surf_tally(itally,k) += origmass * fluxscale;
      if (!transparent) {
        if (ip) a_array_surf_tally(itally,k) -= imass * fluxscale;
        if (jp) a_array_surf_tally(itally,k) -= jmass * fluxscale;
      }
      k++;
      break;
    case MFLUXIN:
      a_array_surf_tally(itally,k) += origmass * fluxscale;
      k++;
      break;

    // forces and torques

    case FX:
      if (!fflag) {
        fflag = 1;
        pdelta_force[0] = pdelta_force[1] = pdelta_force[2] = 0.0;
        if (iorig) MathExtraKokkos::axpy3(-origmass,vorig,pdelta_force);
        if (ip) MathExtraKokkos::axpy3(imass,ip->v,pdelta_force);
        if (jp) MathExtraKokkos::axpy3(jmass,jp->v,pdelta_force);
      }
      a_array_surf_tally(itally,k++) -= pdelta_force[0] * nfactor_inverse;
      break;
    case FY:
      if (!fflag) {
        fflag = 1;
        pdelta_force[0] = pdelta_force[1] = pdelta_force[2] = 0.0;
        if (iorig) MathExtraKokkos::axpy3(-origmass,vorig,pdelta_force);
        if (ip) MathExtraKokkos::axpy3(imass,ip->v,pdelta_force);
        if (jp) MathExtraKokkos::axpy3(jmass,jp->v,pdelta_force);
      }
      a_array_surf_tally(itally,k++) -= pdelta_force[1] * nfactor_inverse;
      break;
    case FZ:
      if (!fflag) {
        fflag = 1;
        pdelta_force[0] = pdelta_force[1] = pdelta_force[2] = 0.0;
        if (iorig) MathExtraKokkos::axpy3(-origmass,vorig,pdelta_force);
        if (ip) MathExtraKokkos::axpy3(imass,ip->v,pdelta_force);
        if (jp) MathExtraKokkos::axpy3(jmass,jp->v,pdelta_force);
      }
      a_array_surf_tally(itally,k++) -= pdelta_force[2] * nfactor_inverse;
      break;

    // for torque, xcollide should be same for any non-NULL particle

    case TX:
      if (!fflag) {
        fflag = 1;
        MathExtraKokkos::scale3(-origmass,vorig,pdelta_force);
        if (ip) MathExtraKokkos::axpy3(imass,ip->v,pdelta_force);
        if (jp) MathExtraKokkos::axpy3(jmass,jp->v,pdelta_force);
      }
      if (!tqflag) {
        tqflag = 1;
        if (ip) xcollide = ip->x;
        else xcollide = iorig->x;
        MathExtraKokkos::sub3(xcollide,com,rdelta);
        MathExtraKokkos::cross3(rdelta,pdelta_force,torque);
      }
      a_array_surf_tally(itally,k++) -= torque[0] * nfactor_inverse;
      break;
    case TY:
      if (!fflag) {
        fflag = 1;
        MathExtraKokkos::scale3(-origmass,vorig,pdelta_force);
        if (ip) MathExtraKokkos::axpy3(imass,ip->v,pdelta_force);
        if (jp) MathExtraKokkos::axpy3(jmass,jp->v,pdelta_force);
      }
      if (!tqflag) {
        tqflag = 1;
        if (ip) xcollide = ip->x;
        else xcollide = iorig->x;
        MathExtraKokkos::sub3(xcollide,com,rdelta);
        MathExtraKokkos::cross3(rdelta,pdelta_force,torque);
      }
      a_array_surf_tally(itally,k++) -= torque[1] * nfactor_inverse;
      break;
    case TZ:
      if (!fflag) {
        fflag = 1;
        MathExtraKokkos::scale3(-origmass,vorig,pdelta_force);
        if (ip) MathExtraKokkos::axpy3(imass,ip->v,pdelta_force);
        if (jp) MathExtraKokkos::axpy3(jmass,jp->v,pdelta_force);
      }
      if (!tqflag) {
        tqflag = 1;
        if (ip) xcollide = ip->x;
        else xcollide = iorig->x;
        MathExtraKokkos::sub3(xcollide,com,rdelta);
        MathExtraKokkos::cross3(rdelta,pdelta_force,torque);
      }
      a_array_surf_tally(itally,k++) -= torque[2] * nfactor_inverse;
      break;

    // pressures

    case PRESS:
      if (!nflag && !tflag) {
        pdelta[0] = pdelta[1] = pdelta[2] = 0.0;
        if (iorig) MathExtraKokkos::axpy3(-origmass,vorig,pdelta);
        if (ip) MathExtraKokkos::axpy3(imass,ip->v,pdelta);
        if (jp) MathExtraKokkos::axpy3(jmass,jp->v,pdelta);
      }
      a_array_surf_tally(itally,k++) += MathExtraKokkos::dot3(pdelta,norm) * fluxscale;
      break;

    case XPRESS:
      if (!nflag) {
        nflag = 1;
        pdelta[0] = pdelta[1] = pdelta[2] = 0.0;
        if (iorig) MathExtraKokkos::axpy3(-origmass,vorig,pdelta);
        if (ip) MathExtraKokkos::axpy3(imass,ip->v,pdelta);
        if (jp) MathExtraKokkos::axpy3(jmass,jp->v,pdelta);
        MathExtraKokkos::scale3(MathExtraKokkos::dot3(pdelta,norm),norm,pnorm);
      }
      a_array_surf_tally(itally,k++) -= pnorm[0] * fluxscale;
      break;
    case YPRESS:
      if (!nflag) {
        nflag = 1;
        pdelta[0] = pdelta[1] = pdelta[2] = 0.0;
        if (iorig) MathExtraKokkos::axpy3(-origmass,vorig,pdelta);
        if (ip) MathExtraKokkos::axpy3(imass,ip->v,pdelta);
        if (jp) MathExtraKokkos::axpy3(jmass,jp->v,pdelta);
        MathExtraKokkos::scale3(MathExtraKokkos::dot3(pdelta,norm),norm,pnorm);
      }
      a_array_surf_tally(itally,k++) -= pnorm[1] * fluxscale;
      break;
    case ZPRESS:
      if (!nflag) {
        nflag = 1;
        pdelta[0] = pdelta[1] = pdelta[2] = 0.0;
        if (iorig) MathExtraKokkos::axpy3(-origmass,vorig,pdelta);
        if (ip) MathExtraKokkos::axpy3(imass,ip->v,pdelta);
        if (jp) MathExtraKokkos::axpy3(jmass,jp->v,pdelta);
        MathExtraKokkos::scale3(MathExtraKokkos::dot3(pdelta,norm),norm,pnorm);
      }
      a_array_surf_tally(itally,k++) -= pnorm[2] * fluxscale;
      break;

    case XSHEAR:
      if (!tflag) {
        tflag = 1;
        pdelta[0] = pdelta[1] = pdelta[2] = 0.0;
        if (iorig) MathExtraKokkos::axpy3(-origmass,vorig,pdelta);
        if (ip) MathExtraKokkos::axpy3(imass,ip->v,pdelta);
        if (jp) MathExtraKokkos::axpy3(jmass,jp->v,pdelta);
        MathExtraKokkos::scale3(MathExtraKokkos::dot3(pdelta,norm),norm,pnorm);
        MathExtraKokkos::sub3(pdelta,pnorm,ptang);
      }
      a_array_surf_tally(itally,k++) -= ptang[0] * fluxscale;
      break;
    case YSHEAR:
      if (!tflag) {
        tflag = 1;
        pdelta[0] = pdelta[1] = pdelta[2] = 0.0;
        if (iorig) MathExtraKokkos::axpy3(-origmass,vorig,pdelta);
        if (ip) MathExtraKokkos::axpy3(imass,ip->v,pdelta);
        if (jp) MathExtraKokkos::axpy3(jmass,jp->v,pdelta);
        MathExtraKokkos::scale3(MathExtraKokkos::dot3(pdelta,norm),norm,pnorm);
        MathExtraKokkos::sub3(pdelta,pnorm,ptang);
      }
      a_array_surf_tally(itally,k++) -= ptang[1] * fluxscale;
      break;
    case ZSHEAR:
      if (!tflag) {
        tflag = 1;
        pdelta[0] = pdelta[1] = pdelta[2] = 0.0;
        if (iorig) MathExtraKokkos::axpy3(-origmass,vorig,pdelta);
        if (ip) MathExtraKokkos::axpy3(imass,ip->v,pdelta);
        if (jp) MathExtraKokkos::axpy3(jmass,jp->v,pdelta);
        MathExtraKokkos::scale3(MathExtraKokkos::dot3(pdelta,norm),norm,pnorm);
        MathExtraKokkos::sub3(pdelta,pnorm,ptang);
      }
      a_array_surf_tally(itally,k++) -= ptang[2] * fluxscale;
      break;

    // energies

    case KE:
      if (iorig) vsqpre = origmass * MathExtraKokkos::lensq3(vorig);
      else vsqpre = 0.0;
      if (ip) ivsqpost = imass * MathExtraKokkos::lensq3(ip->v);
      else ivsqpost = 0.0;
      if (jp) jvsqpost = jmass * MathExtraKokkos::lensq3(jp->v);
      else jvsqpost = 0.0;
      if (transparent)
        a_array_surf_tally(itally,k++) += 0.5*mvv2e * vsqpre * fluxscale;
      else
        a_array_surf_tally(itally,k++) -= 0.5*mvv2e * (ivsqpost + jvsqpost - vsqpre) * fluxscale;
      break;
    case EROT:
      if (ip) ierot = ip->erot;
      else ierot = 0.0;
      if (jp) jerot = jp->erot;
      else jerot = 0.0;
      if (transparent)
        a_array_surf_tally(itally,k++) += weight * oerot * fluxscale;
      else
        a_array_surf_tally(itally,k++) -= weight * (ierot + jerot - oerot) * fluxscale;
      break;
    case EVIB:
      if (ip) ievib = ip->evib;
      else ievib = 0.0;
      if (jp) jevib = jp->evib;
      else jevib = 0.0;
      if (transparent)
        a_array_surf_tally(itally,k++) += weight * oevib * fluxscale;
      else
        a_array_surf_tally(itally,k++) -= weight * (ievib + jevib - oevib) * fluxscale;
      break;
    case ECHEM:
      if (reaction && !transparent) {
        int sr_type = sr_type_list[isr];
        int m = sr_map[isr];
        double r_coeff = 0.0;
        if (sr_type == 1)
          r_coeff = sr_kk_prob_copy[m].obj.d_coeffs(reaction-1,1);
        a_array_surf_tally(itally,k++) += weight * r_coeff * fluxscale;
      }
      break;
    case ETOT:
      if (iorig) vsqpre = origmass * MathExtraKokkos::lensq3(vorig);
      else vsqpre = 0.0;
      otherpre = oerot + oevib;
      if (ip) {
        ivsqpost = imass * MathExtraKokkos::lensq3(ip->v);
        iother = ip->erot + ip->evib;
      } else ivsqpost = iother = 0.0;
      if (jp) {
        jvsqpost = jmass * MathExtraKokkos::lensq3(jp->v);
        jother = jp->erot + jp->evib;
      } else jvsqpost = jother = 0.0;
      if (transparent)
        etot = -0.5*mvv2e*vsqpre - weight*otherpre;
      else {
        etot = 0.5*mvv2e*(ivsqpost + jvsqpost - vsqpre) +
          weight * (iother + jother - otherpre);
        if (reaction) {
          int sr_type = sr_type_list[isr];
          int m = sr_map[isr];
          double r_coeff = 0.0;
          if (sr_type == 1)
            r_coeff = sr_kk_prob_copy[m].obj.d_coeffs(reaction-1,1);
          etot -= weight * r_coeff;
        }
      }
      a_array_surf_tally(itally,k++) -= etot * fluxscale;
      break;
    }
  }
}

 private:
  int mvv2e;

  DAT::t_int_1d d_which;

  DAT::tdual_float_2d_lr k_array_surf_tally;
  DAT::t_float_2d_lr d_array_surf_tally;  // tally values for local surfs

  int need_dup;
  Kokkos::Experimental::ScatterView<F_FLOAT**, typename DAT::t_float_2d_lr::array_layout,DeviceType,typename Kokkos::Experimental::ScatterSum,typename Kokkos::Experimental::ScatterDuplicated> dup_array_surf_tally;
  Kokkos::Experimental::ScatterView<F_FLOAT**, typename DAT::t_float_2d_lr::array_layout,DeviceType,typename Kokkos::Experimental::ScatterSum,typename Kokkos::Experimental::ScatterNonDuplicated> ndup_array_surf_tally;

  DAT::t_surfint_1d d_tally2surf;           // tally2surf[I] = surf ID of Ith tally
  DAT::tdual_surfint_1d k_tally2surf;
  DAT::t_int_1d d_surf2tally;         // using Kokkos::UnorderedMap::insert uses too many registers on GPUs

  DAT::t_float_1d d_normflux;         // normalization factor for each surf element

  t_species_1d d_species;
  DAT::t_int_2d d_s2g;

  t_line_1d d_lines;
  t_tri_1d d_tris;

  int sr_type_list[KOKKOS_MAX_TOT_SURF_REACT];
  int sr_map[KOKKOS_MAX_TOT_SURF_REACT];
  KKCopy<SurfReactGlobalKokkos> sr_kk_global_copy[KOKKOS_MAX_SURF_REACT_PER_TYPE];
  KKCopy<SurfReactProbKokkos> sr_kk_prob_copy[KOKKOS_MAX_SURF_REACT_PER_TYPE];

  void grow_tally();
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/

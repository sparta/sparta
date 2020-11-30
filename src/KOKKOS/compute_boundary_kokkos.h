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

#ifdef COMPUTE_CLASS

ComputeStyle(boundary/kk,ComputeBoundaryKokkos)

#else

#ifndef SPARTA_BOUNDARY_SURF_KOKKOS_H
#define SPARTA_BOUNDARY_SURF_KOKKOS_H

#include "compute_boundary.h"
#include "kokkos_base.h"
#include "kokkos_type.h"
#include "math_extra_kokkos.h"

namespace SPARTA_NS {

class ComputeBoundaryKokkos : public ComputeBoundary, public KokkosBase {
 public:
  ComputeBoundaryKokkos(class SPARTA *, int, char **);
  ComputeBoundaryKokkos(class SPARTA* sparta);
  ~ComputeBoundaryKokkos();
  void init();
  void compute_array();
  void clear();
  void pre_boundary_tally();
  void post_boundary_tally();

  enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};         // same as Domain
  enum{PERIODIC,OUTFLOW,REFLECT,SURFACE,AXISYM};  // same as Domain
  enum{NUM,NUMWT,MFLUX,PRESS,XSHEAR,YSHEAR,ZSHEAR,KE,EROT,EVIB,ETOT};

/* ----------------------------------------------------------------------
   tally values for a single particle colliding with boundary iface/istyle
   iorig = particle ip before collision
   ip,jp = particles after collision
   ip = NULL means no particles after collision
   jp = NULL means one particle after collision
   jp != NULL means two particles after collision
------------------------------------------------------------------------- */

template <int ATOMIC_REDUCTION>
KOKKOS_INLINE_FUNCTION
void boundary_tally_kk(int iface, int istyle, int reaction,
                       Particle::OnePart *iorig,
                       Particle::OnePart *ip,
                       Particle::OnePart *jp,
                       const double* norm) const
{
  // skip if species not in mixture group

  int origspecies = iorig->ispecies;
  int igroup = d_s2g(imix,origspecies);
  if (igroup < 0) return;

  // tally all values associated with group into array
  // set nflag and tflag if normal and tangent computation already done once
  // particle weight used for all keywords except NUM
  // styles PERIODIC and OUTFLOW do not have post-bounce velocity

  auto v_myarray = ScatterViewHelper<typename NeedDup<ATOMIC_REDUCTION,DeviceType>::value,decltype(dup_myarray),decltype(ndup_myarray)>::get(dup_myarray,ndup_myarray);
  auto a_myarray = v_myarray.template access<typename AtomicDup<ATOMIC_REDUCTION,DeviceType>::value>();

  double vsqpre,ivsqpost,jvsqpost;
  double ierot,jerot,ievib,jevib,iother,jother,otherpre;
  double vnorm[3],vtang[3],pdelta[3],pnorm[3],ptang[3];

  double origmass,imass,jmass,pre;
  double weight = weightflag ? iorig->weight : 1.0;
  origmass = d_species[origspecies].mass * weight;
  if (ip) imass = d_species[ip->ispecies].mass * weight;
  if (jp) jmass = d_species[jp->ispecies].mass * weight;

  double *vorig = iorig->v;

  int k = igroup*nvalue;
  int nflag = 0;
  int tflag = 0;

  if (istyle == PERIODIC) {
    for (int m = 0; m < nvalue; m++) {
      if (d_which[m] == NUM) a_myarray(iface,k++) += 1.0;
      else if (d_which[m] == NUMWT) a_myarray(iface,k++) += weight;
      else k++;
    }

  } else if (istyle == OUTFLOW) {
    for (int m = 0; m < nvalue; m++) {
      switch (d_which[m]) {
      case NUM:
        a_myarray(iface,k++) += 1.0;
        break;
      case NUMWT:
        a_myarray(iface,k++) += weight;
        break;
      case MFLUX:
        a_myarray(iface,k++) += origmass;
        break;
      case PRESS:
        if (!nflag) {
          nflag = 1;
          pre = MathExtraKokkos::dot3(vorig,norm);
        }
        a_myarray(iface,k++) -= origmass * pre;
        break;
      case XSHEAR:
        if (!tflag) {
          tflag = 1;
          MathExtraKokkos::scale3(MathExtraKokkos::dot3(vorig,norm),norm,vnorm);
          MathExtraKokkos::sub3(vorig,vnorm,vtang);
        }
        a_myarray(iface,k++) += origmass * vtang[0];
        break;
      case YSHEAR:
        if (!tflag) {
          tflag = 1;
          MathExtraKokkos::scale3(MathExtraKokkos::dot3(vorig,norm),norm,vnorm);
          MathExtraKokkos::sub3(vorig,vnorm,vtang);
        }
        a_myarray(iface,k++) += origmass * vtang[1];
        break;
      case ZSHEAR:
        if (!tflag) {
          tflag = 1;
          MathExtraKokkos::scale3(MathExtraKokkos::dot3(vorig,norm),norm,vnorm);
          MathExtraKokkos::sub3(vorig,vnorm,vtang);
        }
        a_myarray(iface,k++) += origmass * vtang[2];
        break;
      case KE:
        vsqpre = MathExtraKokkos::lensq3(vorig);
        a_myarray(iface,k++) += 0.5 * mvv2e * origmass * vsqpre;
        break;
      case EROT:
        a_myarray(iface,k++) += weight * iorig->erot;
        break;
      case EVIB:
        a_myarray(iface,k++) += weight * iorig->evib;
        break;
      case ETOT:
        vsqpre = MathExtraKokkos::lensq3(vorig);
        a_myarray(iface,k++) += 0.5*mvv2e*origmass*vsqpre +
          weight*(iorig->erot+iorig->evib);
        break;
      }
    }

  } else {
    for (int m = 0; m < nvalue; m++) {
      switch (d_which[m]) {
      case NUM:
        a_myarray(iface,k++) += 1.0;
        break;
      case NUMWT:
        a_myarray(iface,k++) += weight;
        break;
      case MFLUX:
        a_myarray(iface,k++) += origmass;
        if (ip) a_myarray(iface,k++) -= imass;
        if (jp) a_myarray(iface,k++) -= jmass;
        break;
      case PRESS:
        MathExtraKokkos::scale3(-origmass,vorig,pdelta);
        if (ip) MathExtraKokkos::axpy3(imass,ip->v,pdelta);
        if (jp) MathExtraKokkos::axpy3(jmass,jp->v,pdelta);
        a_myarray(iface,k++) += MathExtraKokkos::dot3(pdelta,norm);
        break;
      case XSHEAR:
        if (!tflag) {
          tflag = 1;
          MathExtraKokkos::scale3(-origmass,vorig,pdelta);
          if (ip) MathExtraKokkos::axpy3(imass,ip->v,pdelta);
          if (jp) MathExtraKokkos::axpy3(jmass,jp->v,pdelta);
          MathExtraKokkos::scale3(MathExtraKokkos::dot3(pdelta,norm),norm,pnorm);
          MathExtraKokkos::sub3(pdelta,pnorm,ptang);
        }
        a_myarray(iface,k++) -= ptang[0];
        break;
      case YSHEAR:
        if (!tflag) {
          tflag = 1;
          MathExtraKokkos::scale3(-origmass,vorig,pdelta);
          if (ip) MathExtraKokkos::axpy3(imass,ip->v,pdelta);
          if (jp) MathExtraKokkos::axpy3(jmass,jp->v,pdelta);
          MathExtraKokkos::scale3(MathExtraKokkos::dot3(pdelta,norm),norm,pnorm);
          MathExtraKokkos::sub3(pdelta,pnorm,ptang);
        }
        a_myarray(iface,k++) -= ptang[1];
        break;
      case ZSHEAR:
        if (!tflag) {
          tflag = 1;
          MathExtraKokkos::scale3(-origmass,vorig,pdelta);
          if (ip) MathExtraKokkos::axpy3(imass,ip->v,pdelta);
          if (jp) MathExtraKokkos::axpy3(jmass,jp->v,pdelta);
          MathExtraKokkos::scale3(MathExtraKokkos::dot3(pdelta,norm),norm,pnorm);
          MathExtraKokkos::sub3(pdelta,pnorm,ptang);
        }
        a_myarray(iface,k++) -= ptang[2];
        break;
      case KE:
        vsqpre = origmass * MathExtraKokkos::lensq3(vorig);
        if (ip) ivsqpost = imass * MathExtraKokkos::lensq3(ip->v);
        else ivsqpost = 0.0;
        if (jp) jvsqpost = jmass * MathExtraKokkos::lensq3(jp->v);
        else jvsqpost = 0.0;
        a_myarray(iface,k++) -= 0.5*mvv2e * (ivsqpost + jvsqpost - vsqpre);
        break;
      case EROT:
        if (ip) ierot = ip->erot;
        else ierot = 0.0;
        if (jp) jerot = jp->erot;
        else jerot = 0.0;
        a_myarray(iface,k++) -= weight * (ierot + jerot - iorig->erot);
        break;
      case EVIB:
        if (ip) ievib = ip->evib;
        else ievib = 0.0;
        if (jp) jevib = jp->evib;
        else jevib = 0.0;
        a_myarray(iface,k++) -= weight * (ievib + jevib - iorig->evib);
        break;
      case ETOT:
        vsqpre = origmass * MathExtraKokkos::lensq3(vorig);
        otherpre = iorig->erot + iorig->evib;
        if (ip) {
          ivsqpost = imass * MathExtraKokkos::lensq3(ip->v);
          iother = ip->erot + ip->evib;
        } else ivsqpost = iother = 0.0;
        if (jp) {
          jvsqpost = jmass * MathExtraKokkos::lensq3(jp->v);
          jother = jp->erot + jp->evib;
        } else jvsqpost = jother = 0.0;
        a_myarray(iface,k++) -= 0.5*mvv2e*(ivsqpost + jvsqpost - vsqpre) +
          weight * (iother + jother - otherpre);
        break;
      }
    }
  }
}


 private:
  int mvv2e;
  DAT::t_int_1d d_which;

  DAT::tdual_float_2d_lr k_myarray; // local accumulator array
  DAT::t_float_2d_lr d_myarray;

  int need_dup;
  Kokkos::Experimental::ScatterView<F_FLOAT**, typename DAT::t_float_2d_lr::array_layout,DeviceType,typename Kokkos::Experimental::ScatterSum,typename Kokkos::Experimental::ScatterDuplicated> dup_myarray;
  Kokkos::Experimental::ScatterView<F_FLOAT**, typename DAT::t_float_2d_lr::array_layout,DeviceType,typename Kokkos::Experimental::ScatterSum,typename Kokkos::Experimental::ScatterNonDuplicated> ndup_myarray;

  t_species_1d d_species;
  DAT::t_int_2d d_s2g;
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/

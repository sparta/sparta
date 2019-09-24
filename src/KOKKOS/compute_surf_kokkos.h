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

ComputeStyle(surf/kk,ComputeSurfKokkos)

#else

#ifndef SPARTA_COMPUTE_SURF_KOKKOS_H
#define SPARTA_COMPUTE_SURF_KOKKOS_H

#include "compute_surf.h"
#include "kokkos_type.h"
#include "math_extra_kokkos.h"
#include <Kokkos_UnorderedMap.hpp>

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

enum{NUM,NUMWT,MFLUX,FX,FY,FZ,PRESS,XPRESS,YPRESS,ZPRESS,
     XSHEAR,YSHEAR,ZSHEAR,KE,EROT,EVIB,ETOT};

/* ----------------------------------------------------------------------
   tally values for a single particle colliding with surface element isurf
   iorig = particle ip before collision
   ip,jp = particles after collision
   ip = NULL means no particles after collision
   jp = NULL means one particle after collision
   jp != NULL means two particles after collision
------------------------------------------------------------------------- */

template <int ATOMIC_REDUCTION>
KOKKOS_INLINE_FUNCTION
void surf_tally_kk(int isurf, Particle::OnePart *iorig, 
                   Particle::OnePart *ip, Particle::OnePart *jp) const
{
  // skip if isurf not in surface group

  if (dim == 2) {
    if (!(d_lines(isurf).mask & groupbit)) return;
  } else {
    if (!(d_tris(isurf).mask & groupbit)) return;
  }

  // skip if species not in mixture group

  int origspecies = iorig->ispecies;
  int igroup = d_s2g(imix,origspecies);
  if (igroup < 0) return;

  // itally = tally index of isurf
  // grow tally list if needed

  int itally;

  surfint surfID;
  if (dim == 2) surfID = d_lines[isurf].id;
  else surfID = d_tris[isurf].id;

  int h_index = d_surf2tally(surfID);
  if (h_index != -1)
    itally = h_index;
  else {
    itally = d_ntally();

    d_surf2tally(surfID) = itally;

    d_tally2surf[itally] = surfID;
    if (ATOMIC_REDUCTION != 0)
      Kokkos::atomic_fetch_add(&d_ntally(),1);
    else
      d_ntally()++;
  }

  double fluxscale = d_normflux(isurf);

  // tally all values associated with group into array
  // set nflag and tflag after normal and tangent computation is done once
  // particle weight used for all keywords except NUM
  // forcescale factor applied for keywords FX,FY,FZ
  // fluxscale factor applied for all keywords except NUM,FX,FY,FZ

  double vsqpre,ivsqpost,jvsqpost;
  double ierot,jerot,ievib,jevib,iother,jother,otherpre,etot;
  double pdelta[3],pnorm[3],ptang[3],pdelta_force[3];

  double *norm;
  if (dim == 2) norm = d_lines(isurf).norm;
  else norm = d_tris(isurf).norm;

  double origmass,imass,jmass;

  double weight = 1.0;
  if (weightflag) weight = iorig->weight;
  origmass = d_species(origspecies).mass * weight;
  if (ip) imass = d_species(ip->ispecies).mass * weight;
  if (jp) jmass = d_species(jp->ispecies).mass * weight;

  double *vorig = iorig->v;

  int k = igroup*nvalue;
  int fflag = 0;
  int nflag = 0;
  int tflag = 0;

  for (int m = 0; m < nvalue; m++) {
    switch (d_which(m)) {
    case NUM:
      d_array_surf_tally(itally,k++) += 1.0;
      break;
    case NUMWT:
      d_array_surf_tally(itally,k++) += weight;
      break;
    case MFLUX:
      d_array_surf_tally(itally,k++) += origmass;
      if (ip) d_array_surf_tally(itally,k++) -= imass;
      if (jp) d_array_surf_tally(itally,k++) -= jmass;
      break;
    case FX:
      if (!fflag) {
        fflag = 1;
        MathExtraKokkos::scale3(-origmass,vorig,pdelta_force);
        if (ip) MathExtraKokkos::axpy3(imass,ip->v,pdelta_force);
        if (jp) MathExtraKokkos::axpy3(jmass,jp->v,pdelta_force);
      }
      d_array_surf_tally(itally,k++) -= pdelta_force[0] * nfactor_inverse;
      break;
    case FY:
      if (!fflag) {
        fflag = 1;
        MathExtraKokkos::scale3(-origmass,vorig,pdelta_force);
        if (ip) MathExtraKokkos::axpy3(imass,ip->v,pdelta_force);
        if (jp) MathExtraKokkos::axpy3(jmass,jp->v,pdelta_force);
      }
      d_array_surf_tally(itally,k++) -= pdelta_force[1] * nfactor_inverse;
      break;
    case FZ:
      if (!fflag) {
        fflag = 1;
        MathExtraKokkos::scale3(-origmass,vorig,pdelta_force);
        if (ip) MathExtraKokkos::axpy3(imass,ip->v,pdelta_force);
        if (jp) MathExtraKokkos::axpy3(jmass,jp->v,pdelta_force);
      }
      d_array_surf_tally(itally,k++) -= pdelta_force[2] * nfactor_inverse;
      break;
    case PRESS:
      MathExtraKokkos::scale3(-origmass,vorig,pdelta);
      if (ip) MathExtraKokkos::axpy3(imass,ip->v,pdelta);
      if (jp) MathExtraKokkos::axpy3(jmass,jp->v,pdelta);
      d_array_surf_tally(itally,k++) += MathExtraKokkos::dot3(pdelta,norm) * fluxscale;
      break;
    case XPRESS:
      if (!nflag) {
        nflag = 1;
        MathExtraKokkos::scale3(-origmass,vorig,pdelta);
        if (ip) MathExtraKokkos::axpy3(imass,ip->v,pdelta);
        if (jp) MathExtraKokkos::axpy3(jmass,jp->v,pdelta);
        MathExtraKokkos::scale3(MathExtraKokkos::dot3(pdelta,norm),norm,pnorm);
      }
      d_array_surf_tally(itally,k++) -= pnorm[0] * fluxscale;
      break;
    case YPRESS:
      if (!nflag) {
        nflag = 1;
        MathExtraKokkos::scale3(-origmass,vorig,pdelta);
        if (ip) MathExtraKokkos::axpy3(imass,ip->v,pdelta);
        if (jp) MathExtraKokkos::axpy3(jmass,jp->v,pdelta);
        MathExtraKokkos::scale3(MathExtraKokkos::dot3(pdelta,norm),norm,pnorm);
      }
      d_array_surf_tally(itally,k++) -= pnorm[1] * fluxscale;
      break;
    case ZPRESS:
      if (!nflag) {
        nflag = 1;
        MathExtraKokkos::scale3(-origmass,vorig,pdelta);
        if (ip) MathExtraKokkos::axpy3(imass,ip->v,pdelta);
        if (jp) MathExtraKokkos::axpy3(jmass,jp->v,pdelta);
        MathExtraKokkos::scale3(MathExtraKokkos::dot3(pdelta,norm),norm,pnorm);
      }
      d_array_surf_tally(itally,k++) -= pnorm[2] * fluxscale;
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
      d_array_surf_tally(itally,k++) -= ptang[0] * fluxscale;
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
      d_array_surf_tally(itally,k++) -= ptang[1] * fluxscale;
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
      d_array_surf_tally(itally,k++) -= ptang[2] * fluxscale;
      break;
    case KE:
      vsqpre = origmass * MathExtraKokkos::lensq3(vorig);
      if (ip) ivsqpost = imass * MathExtraKokkos::lensq3(ip->v);
      else ivsqpost = 0.0;
      if (jp) jvsqpost = jmass * MathExtraKokkos::lensq3(jp->v);
      else jvsqpost = 0.0;
      d_array_surf_tally(itally,k++) -= 0.5*mvv2e * (ivsqpost + jvsqpost - vsqpre) * fluxscale;
      break;
    case EROT:
      if (ip) ierot = ip->erot;
      else ierot = 0.0;
      if (jp) jerot = jp->erot;
      else jerot = 0.0;
      d_array_surf_tally(itally,k++) -= weight * (ierot + jerot - iorig->erot) * fluxscale;
      break;
    case EVIB:
      if (ip) ievib = ip->evib;
      else ievib = 0.0;
      if (jp) jevib = jp->evib;
      else jevib = 0.0;
      d_array_surf_tally(itally,k++) -= weight * (ievib + jevib - iorig->evib) * fluxscale;
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
      etot = 0.5*mvv2e*(ivsqpost + jvsqpost - vsqpre) + 
        weight * (iother + jother - otherpre);
      d_array_surf_tally(itally,k++) -= etot * fluxscale;
      break;
    }
  }
}

 private:
  int mvv2e;

  DAT::t_int_scalar d_ntally;
  DAT::t_int_1d d_which;

  DAT::t_float_2d d_array_surf_tally;  // tally values for local surfs
  DAT::tdual_float_2d k_array_surf_tally;
  DAT::t_surfint_1d d_tally2surf;           // tally2surf[I] = surf ID of Ith tally
  DAT::tdual_surfint_1d k_tally2surf;
  DAT::t_int_1d d_surf2tally;         // using Kokkos::UnorderedMap::insert uses too many registers on GPUs

  DAT::t_float_1d d_normflux;         // normalization factor for each surf element

  t_species_1d d_species;
  DAT::t_int_2d d_s2g;

  t_line_1d d_lines;
  t_tri_1d d_tris;
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/

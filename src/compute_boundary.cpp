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

#include "string.h"
#include "compute_boundary.h"
#include "particle.h"
#include "mixture.h"
#include "grid.h"
#include "surf.h"
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"
#include "math_const.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};         // same as Domain
enum{PERIODIC,OUTFLOW,REFLECT,SURFACE,AXISYM};  // same as Domain
enum{NUM,NUMWT,NFLUX,MFLUX,PRESS,XSHEAR,YSHEAR,ZSHEAR,KE,EROT,EVIB,ETOT};

/* ---------------------------------------------------------------------- */

ComputeBoundary::ComputeBoundary(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute boundary command");

  imix = particle->find_mixture(arg[2]);
  if (imix < 0) error->all(FLERR,"Compute boundary mixture ID does not exist");

  nvalue = narg - 3;
  which = new int[nvalue];

  nvalue = 0;
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"n") == 0) which[nvalue++] = NUM;
    else if (strcmp(arg[iarg],"nwt") == 0) which[nvalue++] = NUMWT;
    else if (strcmp(arg[iarg],"nflux") == 0) which[nvalue++] = NFLUX;
    else if (strcmp(arg[iarg],"mflux") == 0) which[nvalue++] = MFLUX;
    else if (strcmp(arg[iarg],"press") == 0) which[nvalue++] = PRESS;
    else if (strcmp(arg[iarg],"shx") == 0) which[nvalue++] = XSHEAR;
    else if (strcmp(arg[iarg],"shy") == 0) which[nvalue++] = YSHEAR;
    else if (strcmp(arg[iarg],"shz") == 0) which[nvalue++] = ZSHEAR;
    else if (strcmp(arg[iarg],"ke") == 0) which[nvalue++] = KE;
    else if (strcmp(arg[iarg],"erot") == 0) which[nvalue++] = EROT;
    else if (strcmp(arg[iarg],"evib") == 0) which[nvalue++] = EVIB;
    else if (strcmp(arg[iarg],"etot") == 0) which[nvalue++] = ETOT;
    else error->all(FLERR,"Illegal compute boundary command");
    iarg++;
  }

  boundary_tally_flag = 1;
  timeflag = 1;
  array_flag = 1;
  ngroup = particle->mixture[imix]->ngroup;
  ntotal = ngroup*nvalue;
  nrow = 2 * domain->dimension;
  size_array_rows = nrow;
  size_array_cols = ntotal;

  memory->create(array,size_array_rows,size_array_cols,"boundary:array");
  memory->create(myarray,size_array_rows,size_array_cols,"boundary:array");
}

/* ---------------------------------------------------------------------- */

ComputeBoundary::~ComputeBoundary()
{
  if (copy || copymode) return;

  delete [] which;
  memory->destroy(array);
  memory->destroy(myarray);
}

/* ---------------------------------------------------------------------- */

void ComputeBoundary::init()
{
  if (ngroup != particle->mixture[imix]->ngroup)
    error->all(FLERR,
               "Number of groups in compute boundary mixture has changed");

  // set normflux based on box face area and timestep size

  double nfactor = update->dt/update->fnum;
  if (domain->dimension == 2) {
    normflux[XLO] = normflux[XHI] = domain->yprd * nfactor;

    if (!domain->axisymmetric)
      normflux[YLO] = normflux[YHI] = domain->xprd * nfactor;
    else {
      // normflux[YLO] is actually 0 in axisymmetric case
      //  but is used in tally normalization even though numerator will be 0
      normflux[YLO] = normflux[YHI] = 2 * MY_PI * domain->xprd * domain->yprd * nfactor;
    }
  } else if (domain->dimension == 3) {
    normflux[XLO] = normflux[XHI] = domain->yprd*domain->zprd * nfactor;
    normflux[YLO] = normflux[YHI] = domain->xprd*domain->zprd * nfactor;
    normflux[ZLO] = normflux[ZHI] = domain->xprd*domain->yprd * nfactor;
  }

  // set weightflag if cell weighting is enabled
  // else weight = 1.0 for all particles

  weight = 1.0;
  if (grid->cellweightflag) weightflag = 1;
  else weightflag = 0;

  // initialize tally array in case accessed before a tally timestep

  clear();
}

/* ---------------------------------------------------------------------- */

void ComputeBoundary::compute_array()
{
  invoked_array = update->ntimestep;

  // sum tallies across processors

  MPI_Allreduce(&myarray[0][0],&array[0][0],nrow*ntotal,
                MPI_DOUBLE,MPI_SUM,world);

  // normalize tallies

  int m;
  for (int j = 0; j < ntotal; j++) {
    m = j % nvalue;
    if (which[m] != NUM && which[m] != NUMWT)
      for (int i = 0; i < size_array_rows; i++)
        array[i][j] /= normflux[i];
  }
}

/* ---------------------------------------------------------------------- */

void ComputeBoundary::clear()
{
  // reset tally values to zero
  // called by Update at beginning of timesteps boundary tallying is done

  for (int i = 0; i < size_array_rows; i++)
    for (int j = 0; j < ntotal; j++)
      myarray[i][j] = 0.0;
}

/* ----------------------------------------------------------------------
   tally values for a single particle colliding with boundary iface/istyle,
     performing reaction (1 to N)
   iorig = particle ip before collision
   ip,jp = particles after collision
   ip = NULL means no particles after collision
   jp = NULL means one particle after collision
   jp != NULL means two particles after collision
------------------------------------------------------------------------- */

void ComputeBoundary::boundary_tally(int iface, int istyle, int reaction,
                                     Particle::OnePart *iorig,
                                     Particle::OnePart *ip,
                                     Particle::OnePart *jp)
{
  // skip if species not in mixture group

  int origspecies = iorig->ispecies;
  int igroup = particle->mixture[imix]->species2group[origspecies];
  if (igroup < 0) return;

  // tally all values associated with group into array
  // set nflag and tflag if normal and tangent computation already done once
  // particle weight used for all keywords except NUM
  // styles PERIODIC and OUTFLOW do not have post-bounce velocity

  double vsqpre,ivsqpost,jvsqpost;
  double ierot,jerot,ievib,jevib,iother,jother,otherpre;
  double vnorm[3],vtang[3],pdelta[3],pnorm[3],ptang[3];

  double *norm = domain->norm[iface];

  double origmass,imass,jmass,pre;
  if (weightflag) weight = iorig->weight;
  origmass = particle->species[origspecies].mass * weight;
  if (ip) imass = particle->species[ip->ispecies].mass * weight;
  if (jp) jmass = particle->species[jp->ispecies].mass * weight;

  double *vorig = iorig->v;
  double mvv2e = update->mvv2e;

  double *vec = myarray[iface];
  int k = igroup*nvalue;
  int nflag = 0;
  int tflag = 0;

  if (istyle == PERIODIC) {
    for (int m = 0; m < nvalue; m++) {
      if (which[m] == NUM) myarray[iface][k++] += 1.0;
      else if (which[m] == NUMWT) myarray[iface][k++] += weight;
      else k++;
    }

  } else if (istyle == OUTFLOW) {
    for (int m = 0; m < nvalue; m++) {
      switch (which[m]) {
      case NUM:
        vec[k++] += 1.0;
        break;
      case NUMWT:
        vec[k++] += weight;
        break;
      case NFLUX:
        vec[k++] += weight;
        break;
      case MFLUX:
        vec[k++] += origmass;
        break;
      case PRESS:
        if (!nflag) {
          nflag = 1;
          pre = MathExtra::dot3(vorig,norm);
        }
        vec[k++] -= origmass * pre;
        break;
      case XSHEAR:
        if (!tflag) {
          tflag = 1;
          MathExtra::scale3(MathExtra::dot3(vorig,norm),norm,vnorm);
          MathExtra::sub3(vorig,vnorm,vtang);
        }
        vec[k++] += origmass * vtang[0];
        break;
      case YSHEAR:
        if (!tflag) {
          tflag = 1;
          MathExtra::scale3(MathExtra::dot3(vorig,norm),norm,vnorm);
          MathExtra::sub3(vorig,vnorm,vtang);
        }
        vec[k++] += origmass * vtang[1];
        break;
      case ZSHEAR:
        if (!tflag) {
          tflag = 1;
          MathExtra::scale3(MathExtra::dot3(vorig,norm),norm,vnorm);
          MathExtra::sub3(vorig,vnorm,vtang);
        }
        vec[k++] += origmass * vtang[2];
        break;
      case KE:
        vsqpre = MathExtra::lensq3(vorig);
        vec[k++] += 0.5 * mvv2e * origmass * vsqpre;
        break;
      case EROT:
        vec[k++] += weight * iorig->erot;
        break;
      case EVIB:
        vec[k++] += weight * iorig->evib;
        break;
      case ETOT:
        vsqpre = MathExtra::lensq3(vorig);
        vec[k++] += 0.5*mvv2e*origmass*vsqpre +
          weight*(iorig->erot+iorig->evib);
        break;
      }
    }

  } else {
    for (int m = 0; m < nvalue; m++) {
      switch (which[m]) {
      case NUM:
        vec[k++] += 1.0;
        break;
      case NUMWT:
        vec[k++] += weight;
        break;
      case NFLUX:
        vec[k] += weight;
        if (ip) vec[k] -= weight;
        if (jp) vec[k] -= weight;
        k++;
        break;
      case MFLUX:
        vec[k] += origmass;
        if (ip) vec[k] -= imass;
        if (jp) vec[k] -= jmass;
        k++;
        break;
      case PRESS:
        MathExtra::scale3(-origmass,vorig,pdelta);
        if (ip) MathExtra::axpy3(imass,ip->v,pdelta);
        if (jp) MathExtra::axpy3(jmass,jp->v,pdelta);
        vec[k++] += MathExtra::dot3(pdelta,norm);
        break;
      case XSHEAR:
        if (!tflag) {
          tflag = 1;
          MathExtra::scale3(-origmass,vorig,pdelta);
          if (ip) MathExtra::axpy3(imass,ip->v,pdelta);
          if (jp) MathExtra::axpy3(jmass,jp->v,pdelta);
          MathExtra::scale3(MathExtra::dot3(pdelta,norm),norm,pnorm);
          MathExtra::sub3(pdelta,pnorm,ptang);
        }
        vec[k++] -= ptang[0];
        break;
      case YSHEAR:
        if (!tflag) {
          tflag = 1;
          MathExtra::scale3(-origmass,vorig,pdelta);
          if (ip) MathExtra::axpy3(imass,ip->v,pdelta);
          if (jp) MathExtra::axpy3(jmass,jp->v,pdelta);
          MathExtra::scale3(MathExtra::dot3(pdelta,norm),norm,pnorm);
          MathExtra::sub3(pdelta,pnorm,ptang);
        }
        vec[k++] -= ptang[1];
        break;
      case ZSHEAR:
        if (!tflag) {
          tflag = 1;
          MathExtra::scale3(-origmass,vorig,pdelta);
          if (ip) MathExtra::axpy3(imass,ip->v,pdelta);
          if (jp) MathExtra::axpy3(jmass,jp->v,pdelta);
          MathExtra::scale3(MathExtra::dot3(pdelta,norm),norm,pnorm);
          MathExtra::sub3(pdelta,pnorm,ptang);
        }
        vec[k++] -= ptang[2];
        break;
      case KE:
        vsqpre = origmass * MathExtra::lensq3(vorig);
        if (ip) ivsqpost = imass * MathExtra::lensq3(ip->v);
        else ivsqpost = 0.0;
        if (jp) jvsqpost = jmass * MathExtra::lensq3(jp->v);
        else jvsqpost = 0.0;
        vec[k++] -= 0.5*mvv2e * (ivsqpost + jvsqpost - vsqpre);
        break;
      case EROT:
        if (ip) ierot = ip->erot;
        else ierot = 0.0;
        if (jp) jerot = jp->erot;
        else jerot = 0.0;
        vec[k++] -= weight * (ierot + jerot - iorig->erot);
        break;
      case EVIB:
        if (ip) ievib = ip->evib;
        else ievib = 0.0;
        if (jp) jevib = jp->evib;
        else jevib = 0.0;
        vec[k++] -= weight * (ievib + jevib - iorig->evib);
        break;
      case ETOT:
        vsqpre = origmass * MathExtra::lensq3(vorig);
        otherpre = iorig->erot + iorig->evib;
        if (ip) {
          ivsqpost = imass * MathExtra::lensq3(ip->v);
          iother = ip->erot + ip->evib;
        } else ivsqpost = iother = 0.0;
        if (jp) {
          jvsqpost = jmass * MathExtra::lensq3(jp->v);
          jother = jp->erot + jp->evib;
        } else jvsqpost = jother = 0.0;
        vec[k++] -= 0.5*mvv2e*(ivsqpost + jvsqpost - vsqpre) +
          weight * (iother + jother - otherpre);
        break;
      }
    }
  }
}

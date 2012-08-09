/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2012) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "string.h"
#include "compute_boundary.h"
#include "particle.h"
#include "mixture.h"
#include "surf.h"
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{NUM,PRESS,XSHEAR,YSHEAR,ZSHEAR,KE};
enum{PERIODIC,OUTFLOW,REFLECT,SURFACE};     // same as Domain
enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};     // same as Domain

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
    else if (strcmp(arg[iarg],"press") == 0) which[nvalue++] = PRESS;
    else if (strcmp(arg[iarg],"shx") == 0) which[nvalue++] = XSHEAR;
    else if (strcmp(arg[iarg],"shy") == 0) which[nvalue++] = YSHEAR;
    else if (strcmp(arg[iarg],"shz") == 0) which[nvalue++] = ZSHEAR;
    else if (strcmp(arg[iarg],"ke") == 0) which[nvalue++] = KE;
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

  double dt = update->dt;
  if (domain->dimension == 2) {
    normflux[XLO] = normflux[XHI] = domain->yprd * dt;
    normflux[YLO] = normflux[YHI] = domain->xprd * dt;
  } else if (domain->dimension == 3) {
    normflux[XLO] = normflux[XHI] = domain->yprd*domain->zprd * dt;
    normflux[YLO] = normflux[YHI] = domain->xprd*domain->xprd * dt;
    normflux[ZLO] = normflux[ZHI] = domain->xprd*domain->yprd * dt;
  }

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
    if (which[m] != NUM) 
      for (int i = 0; i < size_array_rows; i++) {
        array[i][j] /= normflux[i];
      }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeBoundary::clear()
{
  for (int i = 0; i < size_array_rows; i++)
    for (int j = 0; j < ntotal; j++)
      myarray[i][j] = 0.0;
}

/* ---------------------------------------------------------------------- */

void ComputeBoundary::boundary_tally(int iface, int istyle,
                                     double *vold, Particle::OnePart *p)
{
  // skip species not in mixture group

  int ispecies = p->ispecies;
  int igroup = particle->mixture[imix]->species2group[ispecies];
  if (igroup < 0) return;

  // tally all values associated with group into array
  // styles PERIODIC and OUTFLOW do not have post-bounce velocity
  // nflag and tflag set if normal and tangent computation already done once

  double pre,post,vsqpre,vsqpost;
  double vnorm[3],vdelta[3],vtang[3];
  double *vec;

  double *norm = domain->norm[iface];
  double mass = particle->species[ispecies].mass;
  double mvv2e = update->mvv2e;

  vec = myarray[iface];
  int k = igroup*nvalue;
  int nflag = 0;
  int tflag = 0;

  if (istyle == PERIODIC) {
    for (int m = 0; m < nvalue; m++) {
      if (which[m] == NUM) myarray[iface][k++] += 1.0;
      else k++;
    }

  } else if (istyle == OUTFLOW) {
    for (int m = 0; m < nvalue; m++) {
      switch (which[m]) {
      case NUM:
        vec[k++] += 1.0;
        break;
      case PRESS:
        if (!nflag) {
          nflag = 1;
          pre = MathExtra::dot3(vold,norm);
        }
        vec[k++] -= mass * pre;
        break;
      case XSHEAR:
        if (!tflag) {
          tflag = 1;
          MathExtra::scale3(MathExtra::dot3(vold,norm),norm,vnorm);
          MathExtra::sub3(vold,vnorm,vtang);
        }
        vec[k++] += mass * vtang[0];
        break;
      case YSHEAR:
        if (!tflag) {
          tflag = 1;
          MathExtra::scale3(MathExtra::dot3(vold,norm),norm,vnorm);
          MathExtra::sub3(vold,vnorm,vtang);
        }
        vec[k++] += mass * vtang[1];
        break;
      case ZSHEAR:
        if (!tflag) {
          tflag = 1;
          MathExtra::scale3(MathExtra::dot3(vold,norm),norm,vnorm);
          MathExtra::sub3(vold,vnorm,vtang);
        }
        vec[k++] += mass * vtang[2];
        break;
      case KE:
        vsqpre = MathExtra::lensq3(vold);
        vec[k++] += 0.5 * mvv2e * mass * vsqpre;
        break;
      }
    }

  } else {
    for (int m = 0; m < nvalue; m++) {
      switch (which[m]) {
      case NUM:
        vec[k++] += 1.0;
        break;
      case PRESS:
        if (!nflag) {
          nflag = 1;
          pre = MathExtra::dot3(vold,norm);
          post = MathExtra::dot3(p->v,norm);
        }
        vec[k++] += mass * (post-pre);
        break;
      case XSHEAR:
        if (!tflag) {
          tflag = 1;
          MathExtra::sub3(p->v,vold,vdelta);
          MathExtra::scale3(MathExtra::dot3(vdelta,norm),norm,vnorm);
          MathExtra::sub3(vdelta,vnorm,vtang);
        }
        vec[k++] -= mass * vtang[0];
        break;
      case YSHEAR:
        if (!tflag) {
          tflag = 1;
          MathExtra::sub3(p->v,vold,vdelta);
          MathExtra::scale3(MathExtra::dot3(vdelta,norm),norm,vnorm);
          MathExtra::sub3(vdelta,vnorm,vtang);
        }
        vec[k++] -= mass * vtang[1];
        break;
      case ZSHEAR:
        if (!tflag) {
          tflag = 1;
          MathExtra::sub3(p->v,vold,vdelta);
          MathExtra::scale3(MathExtra::dot3(vdelta,norm),norm,vnorm);
          MathExtra::sub3(vdelta,vnorm,vtang);
        }
        vec[k++] -= mass * vtang[2];
        break;
      case KE:
        vsqpre = MathExtra::lensq3(vold);
        vsqpost = MathExtra::lensq3(p->v);
        vec[k++] -= 0.5*mvv2e*mass * (vsqpost-vsqpre);
        break;
      }
    }
  }
}

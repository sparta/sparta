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
#include "compute_surf.h"
#include "particle.h"
#include "mixture.h"
#include "surf.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{NUM,PRESS,XPRESS,YPRESS,ZPRESS,XSHEAR,YSHEAR,ZSHEAR,KE};
enum{NONE,COUNT,MASSWT,TEMPWT};

#define DELTA 10

/* ---------------------------------------------------------------------- */

ComputeSurf::ComputeSurf(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute surf command");

  imix = particle->find_mixture(arg[2]);
  if (imix < 0) error->all(FLERR,"Compute surf mixture ID does not exist");

  nvalue = narg - 3;
  which = new int[nvalue];

  nvalue = 0;
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"n") == 0) which[nvalue++] = NUM;
    else if (strcmp(arg[iarg],"press") == 0) which[nvalue++] = PRESS;
    else if (strcmp(arg[iarg],"px") == 0) which[nvalue++] = XPRESS;
    else if (strcmp(arg[iarg],"py") == 0) which[nvalue++] = YPRESS;
    else if (strcmp(arg[iarg],"pz") == 0) which[nvalue++] = ZPRESS;
    else if (strcmp(arg[iarg],"shx") == 0) which[nvalue++] = XSHEAR;
    else if (strcmp(arg[iarg],"shy") == 0) which[nvalue++] = YSHEAR;
    else if (strcmp(arg[iarg],"shz") == 0) which[nvalue++] = ZSHEAR;
    else if (strcmp(arg[iarg],"ke") == 0) which[nvalue++] = KE;
    else error->all(FLERR,"Illegal compute surf command");
    iarg++;
  }

  surf_tally_flag = 1;
  timeflag = 1;
  per_surf_flag = 1;
  ngroup = particle->mixture[imix]->ngroup;
  ntotal = ngroup*nvalue;
  size_per_surf_cols = ntotal;

  nlocal = maxlocal = 0;
  glob2loc = loc2glob = NULL;
  array_surf = NULL;

  memory->create(value_norm_style,ngroup,nvalue,"surf:value_norm_style");
  norm_count = new double*[ngroup];
  norm_mass = new double*[ngroup];
  norm_temp = new double*[ngroup];
}

/* ---------------------------------------------------------------------- */

ComputeSurf::~ComputeSurf()
{
  delete [] which;
  memory->destroy(glob2loc);
  memory->destroy(loc2glob);
  memory->destroy(array_surf);

  memory->destroy(value_norm_style);
  for (int i = 0; i < ngroup; i++) {
    memory->destroy(norm_count[i]);
    memory->destroy(norm_mass[i]);
    memory->destroy(norm_temp[i]);
  }
  delete [] norm_count;
  delete [] norm_mass;
  delete [] norm_temp;
}

/* ---------------------------------------------------------------------- */

void ComputeSurf::init()
{
  if (ngroup != particle->mixture[imix]->ngroup)
    error->all(FLERR,"Number of groups in compute surf mixture has changed");

  // one-time allocation of accumulators and norms
  // cannot allocate norms until now since depends on group sizes

  if (array_surf == NULL) {
    memory->create(array_surf,surf->nlocal,ntotal,"surf:array_surf");

    for (int i = 0; i < ngroup; i++)
      for (int j = 0; j < nvalue; j++) {
        if (which[j] == NUM) value_norm_style[i][j] = NONE;
        else if (which[j] == KE) value_norm_style[i][j] = COUNT; 
        else if (particle->mixture[imix]->groupsize[i] == 1)
          value_norm_style[i][j] = COUNT; 
        else value_norm_style[i][j] = MASSWT;
      }

    for (int i = 0; i < ngroup; i++) {
      norm_count[i] = norm_mass[i] = norm_temp[i] = NULL;
      for (int j = 0; j < nvalue; j++) {
        if (value_norm_style[i][j] == NONE) continue;
        if (value_norm_style[i][j] == COUNT && norm_count[i] == NULL)
          memory->create(norm_count[i],surf->nlocal,"surf:norm_count");
        if (value_norm_style[i][j] == MASSWT && norm_mass[i] == NULL)
          memory->create(norm_mass[i],surf->nlocal,"surf:norm_mass");
        if (value_norm_style[i][j] == TEMPWT && norm_temp[i] == NULL)
          memory->create(norm_temp[i],surf->nlocal,"surf:norm_temp");
      }
    }
  }

  // local copies

  dimension = domain->dimension;
  lines = surf->lines;
  tris = surf->tris;

  // allocate and initialize glob2loc indices

  if (dimension == 2) nsurf = surf->nline;
  else nsurf = surf->ntri;

  memory->destroy(glob2loc);
  memory->create(glob2loc,nsurf,"surf:glob2loc");
  for (int i = 0; i < nsurf; i++) glob2loc[i] = -1;

  // initialize tally array in case accessed without a tally timestep

  clear();
}

/* ---------------------------------------------------------------------- */

void ComputeSurf::compute_per_surf()
{
  invoked_per_surf = update->ntimestep;

  // maybe should normalize by area and timestep here?

}

/* ---------------------------------------------------------------------- */

void ComputeSurf::clear()
{
  int i,j;
  double *norm;

  // clear out local list and tally values
  // reset all set glob2loc values to -1

  for (i = 0; i < nlocal; i++) {
    for (j = 0; j < ntotal; j++) array_surf[i][j] = 0.0;
    glob2loc[loc2glob[i]] = -1;
  }

  nlocal = 0;

  /*
  for (j = 0; j < ngroup; j++) {
    if (norm = norm_count[j])
      for (i = 0; i < nslocal; i++) norm[i] = 0.0;
    if (norm = norm_mass[j])
      for (i = 0; i < nslocal; i++) norm[i] = 0.0;
    if (norm = norm_temp[j])
      for (i = 0; i < nslocal; i++) norm[i] = 0.0;
  }
  */
}

/* ---------------------------------------------------------------------- */

void ComputeSurf::surf_tally(int isurf, double *vold, Particle::OnePart *p)
{
  // skip species not in mixture group

  int ispecies = p->ispecies;
  int igroup = particle->mixture[imix]->species2group[ispecies];
  if (igroup < 0) return;

  // ilocal = local index of global isurf
  // if 1st particle hitting isurf, add isurf to local list
  // grow local list if needed

  int ilocal = glob2loc[isurf];
  if (ilocal < 0) {
    if (nlocal == maxlocal) grow();
    ilocal = nlocal++;
    loc2glob[ilocal] = isurf;
  }

  // tally all values associated with group into array
  // nflag and tflag set if normal and tangent computation already done once

  double vnorm[3],vdelta[3],vtang[3];
  double pre,post,vsqpre,vsqpost;

  double *norm;
  if (dimension == 2) norm = lines[isurf].norm;
  else norm = tris[isurf].norm;

  double mass = particle->species[ispecies].mass;
  double mvv2e = update->mvv2e;

  int k = igroup*nvalue;
  int nflag = 0;
  int tflag = 0;

  for (int m = 0; m < nvalue; m++) {
    switch (which[m]) {
    case NUM:
      array_surf[ilocal][k++] += 1.0;
      break;
    case PRESS:
      pre = MathExtra::dot3(vold,norm);
      post = MathExtra::dot3(p->v,norm);
      array_surf[ilocal][k++] += mass * (post-pre);
      break;
    case XPRESS:
      if (!nflag) {
        nflag = 1;
        MathExtra::sub3(p->v,vold,vdelta);
        MathExtra::scale3(MathExtra::dot3(vdelta,norm),norm,vnorm);
      }
      array_surf[ilocal][k++] -= mass * vnorm[0];
      break;
    case YPRESS:
      if (!nflag) {
        nflag = 1;
        MathExtra::sub3(p->v,vold,vdelta);
        MathExtra::scale3(MathExtra::dot3(vdelta,norm),norm,vnorm);
      }
      array_surf[ilocal][k++] -= mass * vnorm[1];
      break;
    case ZPRESS:
      if (!nflag) {
        nflag = 1;
        MathExtra::sub3(p->v,vold,vdelta);
        MathExtra::scale3(MathExtra::dot3(vdelta,norm),norm,vnorm);
      }
      array_surf[ilocal][k++] -= mass * vnorm[2];
      break;
    case XSHEAR:
      if (!tflag) {
        tflag = 1;
        MathExtra::sub3(p->v,vold,vdelta);
        MathExtra::scale3(MathExtra::dot3(vdelta,norm),norm,vnorm);
        MathExtra::sub3(vdelta,vnorm,vtang);
      }
      array_surf[ilocal][k++] -= mass * vtang[0];
      break;
    case YSHEAR:
      if (!tflag) {
        tflag = 1;
        MathExtra::sub3(p->v,vold,vdelta);
        MathExtra::scale3(MathExtra::dot3(vdelta,norm),norm,vnorm);
        MathExtra::sub3(vdelta,vnorm,vtang);
      }
      array_surf[ilocal][k++] -= mass * vtang[1];
      break;
    case ZSHEAR:
      if (!tflag) {
        tflag = 1;
        MathExtra::sub3(p->v,vold,vdelta);
        MathExtra::scale3(MathExtra::dot3(vdelta,norm),norm,vnorm);
        MathExtra::sub3(vdelta,vnorm,vtang);
      }
      array_surf[ilocal][k++] -= mass * vtang[2];
      break;
    case KE:
      vsqpre = MathExtra::lensq3(vold);
      vsqpost = MathExtra::lensq3(p->v);
      array_surf[ilocal][k++] -= 0.5*mvv2e*mass * (vsqpost-vsqpre);
      break;
    }
  }
}

/* ----------------------------------------------------------------------
   return ptr to norm vector used by column N
------------------------------------------------------------------------- */

double *ComputeSurf::normptr(int n)
{
  /*
  int igroup = n / nvalue;
  int ivalue = n % nvalue;
  if (value_norm_style[igroup][ivalue] == COUNT) return norm_count[igroup];
  if (value_norm_style[igroup][ivalue] == MASSWT) return norm_mass[igroup];
  if (value_norm_style[igroup][ivalue] == TEMPWT) return norm_temp[igroup];
  */
  return NULL;
}

/* ---------------------------------------------------------------------- */

void ComputeSurf::grow()
{
  maxlocal += DELTA;
  memory->grow(loc2glob,maxlocal,"surf:loc2glob");
  memory->grow(array_surf,maxlocal,ntotal,"surf:array_surf");
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

bigint ComputeSurf::memory_usage()
{
  bigint bytes = 0;
  bytes += ntotal*maxlocal * sizeof(double);
  bytes += maxlocal * sizeof(int);
  bytes += nsurf * sizeof(int);
  /*
  for (int i = 0; i < ngroup; i++) {
  if (norm_count[i]) bytes += surf->nlocal * sizeof(double);
    if (norm_mass[i]) bytes += surf->nlocal * sizeof(double);
   if (norm_temp[i]) bytes += surf->nlocal * sizeof(double);
  }
  */
  return bytes;
}

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

#include "string.h"
#include "compute_surf.h"
#include "particle.h"
#include "mixture.h"
#include "surf.h"
#include "grid.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{NUM,NUMWT,MFLUX,FX,FY,FZ,PRESS,XPRESS,YPRESS,ZPRESS,
     XSHEAR,YSHEAR,ZSHEAR,KE,EROT,EVIB,ETOT};

#define DELTA 4096

/* ---------------------------------------------------------------------- */

ComputeSurf::ComputeSurf(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal compute surf command");

  int igroup = surf->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Compute surf group ID does not exist");
  groupbit = surf->bitmask[igroup];

  imix = particle->find_mixture(arg[3]);
  if (imix < 0) error->all(FLERR,"Compute surf mixture ID does not exist");
  ngroup = particle->mixture[imix]->ngroup;

  nvalue = narg - 4;
  which = new int[nvalue];

  nvalue = 0;
  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"n") == 0) which[nvalue++] = NUM;
    else if (strcmp(arg[iarg],"nwt") == 0) which[nvalue++] = NUMWT;
    else if (strcmp(arg[iarg],"mflux") == 0) which[nvalue++] = MFLUX;
    else if (strcmp(arg[iarg],"fx") == 0) which[nvalue++] = FX;
    else if (strcmp(arg[iarg],"fy") == 0) which[nvalue++] = FY;
    else if (strcmp(arg[iarg],"fz") == 0) which[nvalue++] = FZ;
    else if (strcmp(arg[iarg],"press") == 0) which[nvalue++] = PRESS;
    else if (strcmp(arg[iarg],"px") == 0) which[nvalue++] = XPRESS;
    else if (strcmp(arg[iarg],"py") == 0) which[nvalue++] = YPRESS;
    else if (strcmp(arg[iarg],"pz") == 0) which[nvalue++] = ZPRESS;
    else if (strcmp(arg[iarg],"shx") == 0) which[nvalue++] = XSHEAR;
    else if (strcmp(arg[iarg],"shy") == 0) which[nvalue++] = YSHEAR;
    else if (strcmp(arg[iarg],"shz") == 0) which[nvalue++] = ZSHEAR;
    else if (strcmp(arg[iarg],"ke") == 0) which[nvalue++] = KE;
    else if (strcmp(arg[iarg],"erot") == 0) which[nvalue++] = EROT;
    else if (strcmp(arg[iarg],"evib") == 0) which[nvalue++] = EVIB;
    else if (strcmp(arg[iarg],"etot") == 0) which[nvalue++] = ETOT;
    else error->all(FLERR,"Illegal compute surf command");
    iarg++;
  }

  // NOTE: maybe per_surf_flag should be 0, since this compute
  //       does not create array_surf, but only array_surf_tally

  surf_tally_flag = 1;
  timeflag = 1;
  per_surf_flag = 1;
  ntotal = ngroup*nvalue;
  size_per_surf_cols = ntotal;

  nlocal = maxlocal = 0;
  glob2loc = loc2glob = NULL;
  array_surf_tally = NULL;
  normflux = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeSurf::~ComputeSurf()
{
  delete [] which;
  memory->destroy(glob2loc);
  memory->destroy(loc2glob);
  memory->destroy(array_surf_tally);
  memory->destroy(normflux);
}

/* ---------------------------------------------------------------------- */

void ComputeSurf::init()
{
  if (ngroup != particle->mixture[imix]->ngroup)
    error->all(FLERR,"Number of groups in compute surf mixture has changed");

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

  // normalization nfactor = dt/fnum

  nfactor = update->dt/update->fnum;
  nfactor_inverse = 1.0/nfactor;

  // normflux for all surface elements, based on area and timestep size
  // store inverse, so can multipy by scale factor when tally
  // store for all surf elements, b/c don't know which ones I need to normalize
  // one-time only initialization

  if (normflux == NULL)  {
    memory->create(normflux,nsurf,"surf:normflux");

    int dimension = domain->dimension;
    int axisymmetric = domain->axisymmetric;
    double tmp;

    for (int i = 0; i < nsurf; i++) {
      if (dimension == 3) normflux[i] = surf->tri_size(i,tmp);
      else if (axisymmetric) normflux[i] = surf->axi_line_size(i);
      else normflux[i] = surf->line_size(i);
      normflux[i] *= nfactor;
      normflux[i] = 1.0/normflux[i];
    }
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

void ComputeSurf::compute_per_surf()
{
  invoked_per_surf = update->ntimestep;

  // no operation to perform, local tallies are already normalized
  // NOTE: this does not produce vector_surf or array_surf
  //       which might be what some callers expect
}

/* ---------------------------------------------------------------------- */

void ComputeSurf::clear()
{
  // reset all set glob2loc values to -1
  // called by Update at beginning of timesteps surf tallying is done

  for (int i = 0; i < nlocal; i++)
    glob2loc[loc2glob[i]] = -1;
  nlocal = 0;
}

/* ----------------------------------------------------------------------
   tally values for a single particle colliding with surface element isurf
   iorig = particle ip before collision
   ip,jp = particles after collision
   ip = NULL means no particles after collision
   jp = NULL means one particle after collision
   jp != NULL means two particles after collision
------------------------------------------------------------------------- */

void ComputeSurf::surf_tally(int isurf, Particle::OnePart *iorig, 
                             Particle::OnePart *ip, Particle::OnePart *jp)
{
  // skip if isurf not in surface group

  if (dimension == 2) {
    if (!(lines[isurf].mask & groupbit)) return;
  } else {
    if (!(tris[isurf].mask & groupbit)) return;
  }

  // skip if species not in mixture group

  int origspecies = iorig->ispecies;
  int igroup = particle->mixture[imix]->species2group[origspecies];
  if (igroup < 0) return;

  // ilocal = local index of global isurf
  // if 1st particle hitting isurf, add isurf to local list
  // grow local list if needed

  double *vec;

  int ilocal = glob2loc[isurf];
  if (ilocal < 0) {
    if (nlocal == maxlocal) grow();
    ilocal = nlocal++;
    loc2glob[ilocal] = isurf;
    glob2loc[isurf] = ilocal;
    vec = array_surf_tally[ilocal];
    for (int i = 0; i < ntotal; i++) vec[i] = 0.0;
  }

  double fluxscale = normflux[isurf];

  // tally all values associated with group into array
  // set nflag and tflag after normal and tangent computation is done once
  // particle weight used for all keywords except NUM
  // forcescale factor applied for keywords FX,FY,FZ
  // fluxscale factor applied for all keywords except NUM,FX,FY,FZ

  double vsqpre,ivsqpost,jvsqpost;
  double ierot,jerot,ievib,jevib,iother,jother,otherpre,etot;
  double pdelta[3],pnorm[3],ptang[3],pdelta_force[3];

  double *norm;
  if (dimension == 2) norm = lines[isurf].norm;
  else norm = tris[isurf].norm;

  double origmass,imass,jmass;
  if (weightflag) weight = iorig->weight;
  origmass = particle->species[origspecies].mass * weight;
  if (ip) imass = particle->species[ip->ispecies].mass * weight;
  if (jp) jmass = particle->species[jp->ispecies].mass * weight;

  double *vorig = iorig->v;
  double mvv2e = update->mvv2e;

  vec = array_surf_tally[ilocal];
  int k = igroup*nvalue;
  int fflag = 0;
  int nflag = 0;
  int tflag = 0;

  for (int m = 0; m < nvalue; m++) {
    switch (which[m]) {
    case NUM:
      vec[k++] += 1.0;
      break;
    case NUMWT:
      vec[k++] += weight;
      break;
    case MFLUX:
      vec[k++] += origmass;
      if (ip) vec[k++] -= imass;
      if (jp) vec[k++] -= jmass;
      break;
    case FX:
      if (!fflag) {
        fflag = 1;
        MathExtra::scale3(-origmass,vorig,pdelta_force);
        if (ip) MathExtra::axpy3(imass,ip->v,pdelta_force);
        if (jp) MathExtra::axpy3(jmass,jp->v,pdelta_force);
      }
      vec[k++] -= pdelta_force[0] * nfactor_inverse;
      break;
    case FY:
      if (!fflag) {
        fflag = 1;
        MathExtra::scale3(-origmass,vorig,pdelta_force);
        if (ip) MathExtra::axpy3(imass,ip->v,pdelta_force);
        if (jp) MathExtra::axpy3(jmass,jp->v,pdelta_force);
      }
      vec[k++] -= pdelta_force[1] * nfactor_inverse;
      break;
    case FZ:
      if (!fflag) {
        fflag = 1;
        MathExtra::scale3(-origmass,vorig,pdelta_force);
        if (ip) MathExtra::axpy3(imass,ip->v,pdelta_force);
        if (jp) MathExtra::axpy3(jmass,jp->v,pdelta_force);
      }
      vec[k++] -= pdelta_force[2] * nfactor_inverse;
      break;
    case PRESS:
      MathExtra::scale3(-origmass,vorig,pdelta);
      if (ip) MathExtra::axpy3(imass,ip->v,pdelta);
      if (jp) MathExtra::axpy3(jmass,jp->v,pdelta);
      vec[k++] += MathExtra::dot3(pdelta,norm) * fluxscale;
      break;
    case XPRESS:
      if (!nflag) {
        nflag = 1;
        MathExtra::scale3(-origmass,vorig,pdelta);
        if (ip) MathExtra::axpy3(imass,ip->v,pdelta);
        if (jp) MathExtra::axpy3(jmass,jp->v,pdelta);
        MathExtra::scale3(MathExtra::dot3(pdelta,norm),norm,pnorm);
      }
      vec[k++] -= pnorm[0] * fluxscale;
      break;
    case YPRESS:
      if (!nflag) {
        nflag = 1;
        MathExtra::scale3(-origmass,vorig,pdelta);
        if (ip) MathExtra::axpy3(imass,ip->v,pdelta);
        if (jp) MathExtra::axpy3(jmass,jp->v,pdelta);
        MathExtra::scale3(MathExtra::dot3(pdelta,norm),norm,pnorm);
      }
      vec[k++] -= pnorm[1] * fluxscale;
      break;
    case ZPRESS:
      if (!nflag) {
        nflag = 1;
        MathExtra::scale3(-origmass,vorig,pdelta);
        if (ip) MathExtra::axpy3(imass,ip->v,pdelta);
        if (jp) MathExtra::axpy3(jmass,jp->v,pdelta);
        MathExtra::scale3(MathExtra::dot3(pdelta,norm),norm,pnorm);
      }
      vec[k++] -= pnorm[2] * fluxscale;
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
      vec[k++] -= ptang[0] * fluxscale;
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
      vec[k++] -= ptang[1] * fluxscale;
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
      vec[k++] -= ptang[2] * fluxscale;
      break;
    case KE:
      vsqpre = origmass * MathExtra::lensq3(vorig);
      if (ip) ivsqpost = imass * MathExtra::lensq3(ip->v);
      else ivsqpost = 0.0;
      if (jp) jvsqpost = jmass * MathExtra::lensq3(jp->v);
      else jvsqpost = 0.0;
      vec[k++] -= 0.5*mvv2e * (ivsqpost + jvsqpost - vsqpre) * fluxscale;
      break;
    case EROT:
      if (ip) ierot = ip->erot;
      else ierot = 0.0;
      if (jp) jerot = jp->erot;
      else jerot = 0.0;
      vec[k++] -= weight * (ierot + jerot - iorig->erot) * fluxscale;
      break;
    case EVIB:
      if (ip) ievib = ip->evib;
      else ievib = 0.0;
      if (jp) jevib = jp->evib;
      else jevib = 0.0;
      vec[k++] -= weight * (ievib + jevib - iorig->evib) * fluxscale;
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
      etot = 0.5*mvv2e*(ivsqpost + jvsqpost - vsqpre) + 
        weight * (iother + jother - otherpre);
      vec[k++] -= etot * fluxscale;
      break;
    }
  }
}

/* ----------------------------------------------------------------------
   return ptr to norm vector used by column N
   input N is value from 1 to Ncols
------------------------------------------------------------------------- */

double *ComputeSurf::normptr(int n)
{
  return NULL;
}

/* ----------------------------------------------------------------------
   return ptr to norm vector used by column N
------------------------------------------------------------------------- */

int ComputeSurf::surfinfo(int *&locptr)
{
  locptr = loc2glob;
  return nlocal;
}

/* ---------------------------------------------------------------------- */

void ComputeSurf::grow()
{
  maxlocal += DELTA;
  memory->grow(loc2glob,maxlocal,"surf:loc2glob");
  memory->grow(array_surf_tally,maxlocal,ntotal,"surf:array_surf_tally");
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
  return bytes;
}

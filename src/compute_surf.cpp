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

enum{NUM,NUMWT,NFLUX,NFLUXIN,MFLUX,MFLUXIN,FX,FY,FZ,PRESS,XPRESS,YPRESS,ZPRESS,
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

  // process input values

  nvalue = 0;
  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"n") == 0) which[nvalue++] = NUM;
    else if (strcmp(arg[iarg],"nwt") == 0) which[nvalue++] = NUMWT;
    else if (strcmp(arg[iarg],"nflux") == 0) which[nvalue++] = NFLUX;
    else if (strcmp(arg[iarg],"nflux_incident") == 0) which[nvalue++] = NFLUXIN;
    else if (strcmp(arg[iarg],"mflux") == 0) which[nvalue++] = MFLUX;
    else if (strcmp(arg[iarg],"mflux_incident") == 0) which[nvalue++] = MFLUXIN;
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
    else break;
    iarg++;
  }

  // process optional keywords

  normarea = 1;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"norm") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Invalid compute surf optional keyword");
      if (strcmp(arg[iarg+1],"flow") == 0) normarea = 0;
      else if (strcmp(arg[iarg+1],"flux") == 0) normarea = 1;
      else error->all(FLERR,"Invalid compute surf optional keyword");
      iarg += 2;
    } else error->all(FLERR,"Invalid compute surf value or optional keyword");
  }

  // setup

  ntotal = ngroup*nvalue;

  per_surf_flag = 1;
  size_per_surf_cols = ntotal;

  surf_tally_flag = 1;
  timeflag = 1;

  ntally = maxtally = 0;
  array_surf_tally = NULL;
  tally2surf = NULL;

  maxsurf = 0;
  array_surf = NULL;
  normflux = NULL;
  combined = 0;

  hash = new MyHash;

  dim = domain->dimension;
}

/* ---------------------------------------------------------------------- */

ComputeSurf::~ComputeSurf()
{
  if (copy || copymode) return;

  delete [] which;
  memory->destroy(array_surf_tally);
  memory->destroy(tally2surf);
  memory->destroy(array_surf);
  memory->destroy(normflux);
  delete hash;
}

/* ---------------------------------------------------------------------- */

void ComputeSurf::init()
{
  if (!surf->exist)
    error->all(FLERR,"Cannot use compute surf when surfs do not exist");
  if (surf->implicit)
    error->all(FLERR,"Cannot use compute surf with implicit surfs");

  if (ngroup != particle->mixture[imix]->ngroup)
    error->all(FLERR,"Number of groups in compute surf mixture has changed");

  // set normflux for all owned + ghost surfs

  init_normflux();

  // set weightflag if cell weighting is enabled
  // else weight = 1.0 for all particles

  weight = 1.0;
  if (grid->cellweightflag) weightflag = 1;
  else weightflag = 0;

  // initialize tally array in case accessed before a tally timestep

  clear();

  combined = 0;
}

/* ----------------------------------------------------------------------
   set normflux for all surfs I store
   all: just nlocal
   distributed: nlocal + nghost
   called by init before each run (in case dt or fnum has changed)
   called whenever grid changes
------------------------------------------------------------------------- */

void ComputeSurf::init_normflux()
{
  // normalization nfactor = dt/fnum

  double nfactor = update->dt/update->fnum;
  nfactor_inverse = 1.0/nfactor;

  // normflux for all surface elements, based on area and timestep size
  // if normarea = 0, area is not included in flux
  //   mass/eng fluxes (mass/area/time) becomes mass/eng flows (mass/time)
  // nsurf = all explicit surfs in this procs grid cells
  // store inverse, so can multipy by scale factor when tally
  // store for all surf elements, b/c don't know which ones I need to normalize

  int nsurf = surf->nlocal + surf->nghost;
  memory->destroy(normflux);
  memory->create(normflux,nsurf,"surf:normflux");

  int axisymmetric = domain->axisymmetric;
  double tmp;

  for (int i = 0; i < nsurf; i++) {
    if (!normarea) normflux[i] = 1.0;
    else if (dim == 3) normflux[i] = surf->tri_size(i,tmp);
    else if (axisymmetric) normflux[i] = surf->axi_line_size(i);
    else normflux[i] = surf->line_size(i);
    normflux[i] *= nfactor;
    normflux[i] = 1.0/normflux[i];
  }
}

/* ----------------------------------------------------------------------
   no operations here, since compute results are stored in tally array
   just used by callers to indicate compute was used
   enables prediction of next step when update needs to tally
------------------------------------------------------------------------- */

void ComputeSurf::compute_per_surf()
{
  invoked_per_surf = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void ComputeSurf::clear()
{
  lines = surf->lines;
  tris = surf->tris;

  // clear hash of tallied surf IDs
  // called by Update at beginning of timesteps surf tallying is done

  hash->clear();
  ntally = 0;
  combined = 0;
}

/* ----------------------------------------------------------------------
   tally values for a single particle in icell
     colliding with surface element isurf, performing reaction (1 to N)
   iorig = particle ip before collision
   ip,jp = particles after collision
   ip = NULL means no particles after collision
   jp = NULL means one particle after collision
   jp != NULL means two particles after collision
------------------------------------------------------------------------- */

void ComputeSurf::surf_tally(int isurf, int icell, int reaction,
                             Particle::OnePart *iorig,
                             Particle::OnePart *ip, Particle::OnePart *jp)
{
  // skip if isurf not in surface group

  if (dim == 2) {
    if (!(lines[isurf].mask & groupbit)) return;
  } else {
    if (!(tris[isurf].mask & groupbit)) return;
  }

  // skip if species not in mixture group

  int origspecies = iorig->ispecies;
  int igroup = particle->mixture[imix]->species2group[origspecies];
  if (igroup < 0) return;

  // itally = tally index of isurf
  // if 1st particle hitting isurf, add surf ID to hash
  // grow tally list if needed

  int itally,transparent;
  double *vec;

  surfint surfID;
  if (dim == 2) {
    surfID = lines[isurf].id;
    transparent = lines[isurf].transparent;
  } else {
    surfID = tris[isurf].id;
    transparent = tris[isurf].transparent;
  }

  if (hash->find(surfID) != hash->end()) itally = (*hash)[surfID];
  else {
    if (ntally == maxtally) grow_tally();
    itally = ntally;
    (*hash)[surfID] = itally;
    tally2surf[itally] = surfID;
    vec = array_surf_tally[itally];
    for (int i = 0; i < ntotal; i++) vec[i] = 0.0;
    ntally++;
  }

  double fluxscale = normflux[isurf];

  // tally all values associated with group into array
  // set nflag and tflag after normal and tangent computation is done once
  // particle weight used for all keywords except NUM
  // forcescale factor applied for keywords FX,FY,FZ
  // fluxscale factor applied for all keywords except NUM,FX,FY,FZ
  // if surf is transparent, all flux tallying is for incident particle only

  double vsqpre,ivsqpost,jvsqpost;
  double ierot,jerot,ievib,jevib,iother,jother,otherpre,etot;
  double pdelta[3],pnorm[3],ptang[3],pdelta_force[3];

  double *norm;
  if (dim == 2) norm = lines[isurf].norm;
  else norm = tris[isurf].norm;

  double origmass,imass,jmass;
  if (weightflag) weight = iorig->weight;
  origmass = particle->species[origspecies].mass * weight;
  if (ip) imass = particle->species[ip->ispecies].mass * weight;
  if (jp) jmass = particle->species[jp->ispecies].mass * weight;

  double *vorig = iorig->v;
  double mvv2e = update->mvv2e;

  vec = array_surf_tally[itally];
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
    case NFLUX:
      vec[k] += weight * fluxscale;
      if (!transparent) {
        if (ip) vec[k] -= weight * fluxscale;
        if (jp) vec[k] -= weight * fluxscale;
      }
      k++;
      break;
    case NFLUXIN:
      vec[k] += weight * fluxscale;
      k++;
      break;
    case MFLUX:
      vec[k] += origmass * fluxscale;
      if (!transparent) {
        if (ip) vec[k] -= imass * fluxscale;
        if (jp) vec[k] -= jmass * fluxscale;
      }
      k++;
      break;
    case MFLUXIN:
      vec[k] += origmass * fluxscale;
      k++;
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
      if (transparent)
        vec[k++] += 0.5*mvv2e * vsqpre * fluxscale;
      else
        vec[k++] -= 0.5*mvv2e * (ivsqpost + jvsqpost - vsqpre) * fluxscale;
      break;
    case EROT:
      if (ip) ierot = ip->erot;
      else ierot = 0.0;
      if (jp) jerot = jp->erot;
      else jerot = 0.0;
      if (transparent)
        vec[k++] += weight * iorig->erot * fluxscale;
      else
        vec[k++] -= weight * (ierot + jerot - iorig->erot) * fluxscale;
      break;
    case EVIB:
      if (ip) ievib = ip->evib;
      else ievib = 0.0;
      if (jp) jevib = jp->evib;
      else jevib = 0.0;
      if (transparent)
        vec[k++] += weight * iorig->evib * fluxscale;
      else
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
      if (transparent)
        etot = -0.5*mvv2e*vsqpre - weight*otherpre;
      else
        etot = 0.5*mvv2e*(ivsqpost + jvsqpost - vsqpre) +
          weight * (iother + jother - otherpre);
      vec[k++] -= etot * fluxscale;
      break;
    }
  }
}

/* ----------------------------------------------------------------------
   return # of tallies and their indices into my local surf list
------------------------------------------------------------------------- */

int ComputeSurf::tallyinfo(surfint *&ptr)
{
  ptr = tally2surf;
  return ntally;
}

/* ----------------------------------------------------------------------
   sum tally values to owning surfs via surf->collate()
------------------------------------------------------------------------- */

void ComputeSurf::post_process_surf()
{
  if (combined) return;
  combined = 1;

  // reallocate array_surf if necessary

  int nown = surf->nown;

  if (nown > maxsurf) {
    memory->destroy(array_surf);
    maxsurf = nown;
    memory->create(array_surf,maxsurf,ntotal,"surf:array_surf");
  }

  // zero array_surf
  // NOTE: is this needed if collate zeroes ?

  int i,j;
  for (i = 0; i < nown; i++)
    for (j = 0; j < ntotal; j++)
      array_surf[i][j] = 0.0;

  // collate entire array of results

  surf->collate_array(ntally,ntotal,tally2surf,array_surf_tally,array_surf);

  /*
  if (array_surf_tally)
    surf->collate_vector(ntally,tally2surf,
                         &array_surf_tally[0][index-1],ntotal,vector_surf);
  else
    surf->collate_vector(ntally,tally2surf,NULL,ntotal,vector_surf);
  */
}


/* ---------------------------------------------------------------------- */

void ComputeSurf::grow_tally()
{
  maxtally += DELTA;
  memory->grow(tally2surf,maxtally,"surf:tally2surf");
  memory->grow(array_surf_tally,maxtally,ntotal,"surf:array_surf_tally");
}

/* ----------------------------------------------------------------------
   reset normflux for my surfs
   called whenever grid changes
------------------------------------------------------------------------- */

void ComputeSurf::reallocate()
{
  init_normflux();
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

bigint ComputeSurf::memory_usage()
{
  bigint bytes = 0;
  bytes += ntotal*maxtally * sizeof(double);    // array_surf_tally
  bytes += maxtally * sizeof(surfint);          // tally2surf
  return bytes;
}

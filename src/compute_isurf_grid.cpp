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
#include "compute_isurf_grid.h"
#include "particle.h"
#include "mixture.h"
#include "surf.h"
#include "grid.h"
#include "update.h"
#include "domain.h"
#include "comm.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{NUM,NUMWT,MFLUX,FX,FY,FZ,PRESS,XPRESS,YPRESS,ZPRESS,
     XSHEAR,YSHEAR,ZSHEAR,KE,EROT,EVIB,ETOT};

#define DELTA 4096

/* ---------------------------------------------------------------------- */

ComputeISurfGrid::ComputeISurfGrid(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal compute isurf/grid command");

  int igroup = grid->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Compute isurf/grid group ID does not exist");
  groupbit = grid->bitmask[igroup];

  imix = particle->find_mixture(arg[3]);
  if (imix < 0) error->all(FLERR,"Compute isurf/grid mixture ID does not exist");
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
    else error->all(FLERR,"Illegal compute isurf/grid command");
    iarg++;
  }

  ntotal = ngroup*nvalue;

  per_grid_flag = 1;
  size_per_grid_cols = ntotal;
  post_process_isurf_grid_flag = 1;

  surf_tally_flag = 1;
  timeflag = 1;

  ntally = maxtally = 0;
  array_surf_tally = NULL;
  tally2surf = NULL;

  maxgrid = 0;
  array_grid = NULL;
  normflux = NULL;
  combined = 0;

  hash = new MyHash;

  dim = domain->dimension;
}

/* ---------------------------------------------------------------------- */

ComputeISurfGrid::~ComputeISurfGrid()
{
  delete [] which;
  memory->destroy(array_surf_tally);
  memory->destroy(tally2surf);
  memory->destroy(array_grid);
  memory->destroy(normflux);
  delete hash;
}

/* ---------------------------------------------------------------------- */

void ComputeISurfGrid::init()
{
  if (!surf->exist)
    error->all(FLERR,"Cannot use compute isurf/grid when surfs do not exist");
  if (!surf->implicit)
    error->all(FLERR,"Cannot use compute isurf/grid with explicit surfs");

  if (ngroup != particle->mixture[imix]->ngroup)
    error->all(FLERR,
               "Number of groups in compute isurf/grid mixture has changed");

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
  reallocate();
}

/* ----------------------------------------------------------------------
   set normflux for all surfs I store
   distributed: nlocal + nghost
   called by init before each run (in case dt or fnum has changed)
   called whenever grid changes
------------------------------------------------------------------------- */

void ComputeISurfGrid::init_normflux()
{
  // normalization nfactor = dt/fnum

  double nfactor = update->dt/update->fnum;
  nfactor_inverse = 1.0/nfactor;

  // normflux for all surface elements, based on area and timestep size
  // nsurf = all explicit surfs in this procs grid cells
  // store inverse, so can multipy by scale factor when tally
  // store for all surf elements, b/c don't know which ones I need to normalize

  int nsurf = surf->nlocal + surf->nghost;
  memory->destroy(normflux);
  memory->create(normflux,nsurf,"isurf/grid:normflux");

  int axisymmetric = domain->axisymmetric;
  double tmp;

  for (int i = 0; i < nsurf; i++) {
    if (dim == 3) normflux[i] = surf->tri_size(i,tmp);
    else if (axisymmetric) normflux[i] = surf->axi_line_size(i);
    else normflux[i] = surf->line_size(i);
    normflux[i] *= nfactor;
    normflux[i] = 1.0/normflux[i];
  }
}

/* ----------------------------------------------------------------------
   no operations here, since compute results are stored in array_grid
   just used by callers to indicate compute was used
   enables prediction of next step when update needs to tally
   NOTE: also need to mark invoked_per_surf?
------------------------------------------------------------------------- */

void ComputeISurfGrid::compute_per_grid()
{
  invoked_per_grid = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void ComputeISurfGrid::clear()
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
   this method is exactly like ComputeSurf::surf_tally()
     except sum tally to to per-grid-cell array_grid
------------------------------------------------------------------------- */

void ComputeISurfGrid::surf_tally(int isurf, int icell, int reaction,
                                   Particle::OnePart *iorig,
                                   Particle::OnePart *ip, Particle::OnePart *jp)
{
  // skip if species not in mixture group

  int origspecies = iorig->ispecies;
  int igroup = particle->mixture[imix]->species2group[origspecies];
  if (igroup < 0) return;

  // itally = tally index of isurf
  // if 1st particle hitting isurf, add surf ID to hash
  // grow tally list if needed
  // for implicit surfs, surfID is really a cellID

  int itally;
  double *vec;

  surfint surfID;
  if (dim == 2) surfID = lines[isurf].id;
  else surfID = tris[isurf].id;

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
    case MFLUX:
      vec[k] += origmass * fluxscale;
      if (ip) vec[k] -= imass * fluxscale;
      if (jp) vec[k] -= jmass * fluxscale;
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
   return # of tallies and their indices into my local surf list
------------------------------------------------------------------------- */

int ComputeISurfGrid::tallyinfo(surfint *&ptr)
{
  ptr = tally2surf;
  return ntally;
}

/* ----------------------------------------------------------------------
   sum surf tallies to owning cells via surf->collate()
   also copy split cell values to sub-cells for use by dump grid
------------------------------------------------------------------------- */

void ComputeISurfGrid::post_process_isurf_grid()
{
  if (combined) return;
  combined = 1;

  // reallocate array_grid if necessary

  int nglocal = grid->nlocal;

  if (nglocal > maxgrid) {
    memory->destroy(array_grid);
    maxgrid = nglocal;
    memory->create(array_grid,maxgrid,ntotal,"isurf/grid:array_grid");
  }

  // zero array_grid

  int i,j;
  for (i = 0; i < nglocal; i++)
    for (j = 0; j < ntotal; j++)
      array_grid[i][j] = 0.0;

  // perform rendezvous comm on tallies to sum them to my grid cells
  // array_surf_tally can be NULL if this proc has performed no tallies

  surf->collate_array_implicit(ntally,ntotal,tally2surf,
                               array_surf_tally,array_grid);

  // zero out result if icell not in grid group
  // can't apply until now, b/c tally included surfs in ghost cells and
  // cinfo does not have mask values for ghost cells

  Grid::ChildInfo *cinfo = grid->cinfo;

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) {
      for (j = 0; j < ntotal; j++)
        array_grid[icell][j] = 0.0;
    }
  }

  // copy split cell values to their sub cells, used by dump grid

  Grid::ChildCell *cells = grid->cells;
  Grid::SplitInfo *sinfo = grid->sinfo;

  int jcell,nsplit;
  int *csubs;

  for (int icell = 0; icell < nglocal; icell++) {
    if (cells[icell].nsplit <= 1) continue;
    nsplit = cells[icell].nsplit;
    csubs = sinfo[cells[icell].isplit].csubs;
    for (int j = 0; j < nsplit; j++) {
      jcell = csubs[j];
      memcpy(array_grid[jcell],array_grid[icell],ntotal*sizeof(double));
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeISurfGrid::grow_tally()
{
  maxtally += DELTA;
  memory->grow(tally2surf,maxtally,"isurf/grid:tally2surf");
  memory->grow(array_surf_tally,maxtally,ntotal,"isurf/grid:array_surf_tally");
}

/* ----------------------------------------------------------------------
   reset normflux for my surfs
   called whenever grid changes
------------------------------------------------------------------------- */

void ComputeISurfGrid::reallocate()
{
  init_normflux();
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

bigint ComputeISurfGrid::memory_usage()
{
  bigint bytes = 0;
  bytes += ntotal*maxgrid * sizeof(double);     // array_grid
  bytes += ntotal*maxtally * sizeof(double);    // array_surf_tally
  bytes += maxtally * sizeof(surfint);          // tally2surf
  return bytes;
}

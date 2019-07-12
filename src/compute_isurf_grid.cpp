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

/* ---------------------------------------------------------------------- */

ComputeISurfGrid::ComputeISurfGrid(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal compute isurf/grid command");

  if (!surf->implicit) 
    error->all(FLERR,"Cannot use compute isurf/grid with explicit surfs");

  int igroup = grid->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Compute isurf/grid group ID does not exist");
  groupbit = surf->bitmask[igroup];

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
  post_process_grid_flag = 1;

  surf_tally_flag = 1;
  timeflag = 1;

  ngtotal = 0;
  vector_grid = NULL;
  array_grid = NULL;
  normflux = NULL;
  nfactor_previous = 0.0;

  maxsend = 0;
  proclist = NULL;
  sbuf = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeISurfGrid::~ComputeISurfGrid()
{
  delete [] which;
  memory->destroy(vector_grid);
  memory->destroy(array_grid);
  memory->destroy(normflux);
  memory->destroy(proclist);
  memory->destroy(sbuf);
}

/* ---------------------------------------------------------------------- */

void ComputeISurfGrid::init()
{
  if (ngroup != particle->mixture[imix]->ngroup)
    error->all(FLERR,
               "Number of groups in compute isurf/grid mixture has changed");

  // local copies

  dimension = domain->dimension;
  lines = surf->lines;
  tris = surf->tris;

  nsurf = surf->nlocal + surf->nghost;

  // normalization nfactor = dt/fnum

  nfactor = update->dt/update->fnum;
  nfactor_inverse = 1.0/nfactor;

  // normflux for all surface elements, based on area and timestep size
  // store inverse, so can multipy by scale factor when tally
  // store for all surf elements, b/c don't know which ones I need to normalize
  // one-time only initialization unless timestep or fnum changes between runs
  // NOTE: needs to be done for owned and ghost surfs
  // NOTE: when grid cells migrate or adapt, may have new surfs

  if (nfactor != nfactor_previous)  {
    memory->destroy(normflux);
    memory->create(normflux,nsurf,"surf:normflux");

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

  nfactor_previous = nfactor;

  // set weightflag if cell weighting is enabled
  // else weight = 1.0 for all particles

  weight = 1.0;
  if (grid->cellweightflag) weightflag = 1;
  else weightflag = 0;

  // allocate per grid cell array
  // zero array in case dump invokes post_process_grid() before run

  reallocate();

  collapsed = 0;

  int i,j;

  int n = grid->nlocal + grid->nghost;
  for (i = 0; i < n; i++)
    for (j = 0; j < ntotal; j++)
      array_grid[i][j] = 0.0;
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
  // reset array_grid - could do this with memset()
  // called by Update at beginning of timesteps surf tallying is done
  // NOTE: how to decide when to do this
  //       should not do it until caller sums tallies with ghost cells
  // NOTE: how to detect when local grid has changed (balance, adapt)

  int i,j;

  int n = grid->nlocal + grid->nghost;
  printf("CLEAR me %d: %d %d\n",comm->me,grid->nlocal,grid->nghost);
  for (i = 0; i < n; i++)
    for (j = 0; j < ntotal; j++)
      array_grid[i][j] = 0.0;

  collapsed = 0;

  cells = grid->cells;
  cinfo = grid->cinfo;
  sinfo = grid->sinfo;
  lines = surf->lines;
  tris = surf->tris;
}

/* ----------------------------------------------------------------------
   tally values for a single particle in icell 
     colliding with surface element isurf
   iorig = particle ip before collision
   ip,jp = particles after collision
   ip = NULL means no particles after collision
   jp = NULL means one particle after collision
   jp != NULL means two particles after collision
   this method is almost exactly like ComputeSurf::surf_tally()
     except sum tally to to per-grid-cell array_grid
------------------------------------------------------------------------- */

void ComputeISurfGrid::surf_tally(int isurf, int icell, Particle::OnePart *iorig, 
                                  Particle::OnePart *ip, Particle::OnePart *jp)
{
  // skip if icell not in grid group

  if (!(cinfo[icell].mask & groupbit)) return;

  // skip if species not in mixture group

  int origspecies = iorig->ispecies;
  int igroup = particle->mixture[imix]->species2group[origspecies];
  if (igroup < 0) return;

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

  // if icell is a sub-cell, reset icell to its parent split cell

  if (cells[icell].nsplit <= 0)
    icell = sinfo[cells[icell].isplit].icell;

  // sum into grid cell array
  // icell may be an owned or ghost cell

  double *vec = array_grid[icell];      
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
   sum ghost cell tallies to owning cells via irregular communication
   also copy split cell values to sub-cells for use by dump grid
   NOTE: function args are currently ignored
------------------------------------------------------------------------- */

double ComputeISurfGrid::post_process_grid(int index, int onecell, int nsample,
                                           double **etally, int *emap,
                                           double *vec, int nstride)
{
  cells = grid->cells;
  sinfo = grid->sinfo;
  int nglocal = grid->nlocal;
  int ngghost = grid->nghost;
  int n = nglocal + ngghost;

  // one-time collapse of ghost cell values into owned cells 

  if (!collapsed) {
    collapsed = 1;

    if (ngghost > maxsend) {
      maxsend = ngghost;
      memory->destroy(proclist);
      memory->create(proclist,maxsend,"isurf/grid:proclist");
      memory->destroy(sbuf);
      memory->create(sbuf,maxsend*(ntotal+1),"isurf/grid:sbuf");
    }

    // pack ghost cell data into sbuf
    // skip cells with no surfs and sub-cells
    // data = ilocal (of icell on other proc) + Ntotal values
    
    int m = 0;
    int nsend = 0;
    for (int icell = nglocal; icell < n; icell++) {
      if (cells[icell].nsurf <= 0) continue;
      if (cells[icell].nsplit <= 0) continue;
      proclist[nsend] = cells[icell].proc;
      sbuf[m++] = cells[icell].ilocal;    // NOTE: need ubuf
      memcpy(&sbuf[m],array_grid[icell],ntotal*sizeof(double));
      m += ntotal;
      nsend++;
    }

    printf("NSEND %d %ld: %d %d\n",comm->me,update->ntimestep,nsend,n);

    // perform irregular comm

    double *rbuf;
    int nrecv = comm->irregular_uniform_neighs(nsend,proclist,(char *) sbuf,
                                               (ntotal+1)*sizeof(double),
                                               (char **) &rbuf);
    
    printf("NRECV %d %d\n",comm->me,nrecv);

    // unpack received data and sum Ntotal values into array_grid

    int j,ilocal;

    m = 0;
    for (int i = 0; i < nrecv; i++) {
      ilocal = static_cast<int> (rbuf[m++]);   // NOTE: need ubuf
      for (j = 0; j < ntotal; j++)
        array_grid[ilocal][j] += rbuf[m++];
    }

    // copy split cell values to each of its sub cells
    // for use by dump grid

    int jcell,nsplit;
    int *csubs;

    for (int icell = nglocal; icell < n; icell++) {
      if (cells[icell].nsplit <= 1) continue;
      nsplit = cells[icell].nsplit;
      csubs = sinfo[cells[icell].isplit].csubs;
      for (int j = 0; j < nsplit; j++) {
        jcell = csubs[j];
        memcpy(array_grid[jcell],array_grid[icell],ntotal*sizeof(double));
      }
    }
  }

  // copy values for array_grid column index into vector_grid

  index--;
  for (int icell = 0; icell < nglocal; icell++)
    vector_grid[icell] = array_grid[icell][index];

  return 0.0;
}

/* ----------------------------------------------------------------------
   reallocate data storage if nglocal has changed
   called by init() and whenever grid changes
------------------------------------------------------------------------- */

void ComputeISurfGrid::reallocate()
{
  if (grid->nlocal + grid->nghost == ngtotal) return;

  memory->destroy(vector_grid);
  memory->destroy(array_grid);
  ngtotal = grid->nlocal + grid->nghost;
  memory->create(vector_grid,ngtotal,"isurf/grid:vector_grid");
  memory->create(array_grid,ngtotal,ntotal,"isurf/grid:array_grid");
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

bigint ComputeISurfGrid::memory_usage()
{
  bigint bytes = 0;
  bytes += ntotal*ngtotal * sizeof(double);     // array_grid
  return bytes;
}

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
#include "compute_grid.h"
#include "particle.h"
#include "mixture.h"
#include "grid.h"
#include "update.h"
#include "modify.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{NUM,NRHO,MASS,U,V,W,USQ,VSQ,WSQ,KE,TEMPERATURE,EROT,TROT,EVIB,TVIB};
enum{NONE,COUNT,MASSWT,RDOF,VDOF};

/* ---------------------------------------------------------------------- */

ComputeGrid::ComputeGrid(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute grid command");

  imix = particle->find_mixture(arg[2]);
  if (imix < 0) error->all(FLERR,"Compute grid mixture ID does not exist");

  nvalue = narg - 3;
  which = new int[nvalue];
  norm_style = new int[nvalue];

  nvalue = 0;
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"n") == 0) which[nvalue++] = NUM;
    else if (strcmp(arg[iarg],"nrho") == 0) which[nvalue++] = NRHO;
    else if (strcmp(arg[iarg],"mass") == 0) which[nvalue++] = MASS;
    else if (strcmp(arg[iarg],"u") == 0) which[nvalue++] = U;
    else if (strcmp(arg[iarg],"v") == 0) which[nvalue++] = V;
    else if (strcmp(arg[iarg],"w") == 0) which[nvalue++] = W;
    else if (strcmp(arg[iarg],"usq") == 0) which[nvalue++] = USQ;
    else if (strcmp(arg[iarg],"vsq") == 0) which[nvalue++] = VSQ;
    else if (strcmp(arg[iarg],"wsq") == 0) which[nvalue++] = WSQ;
    else if (strcmp(arg[iarg],"ke") == 0) which[nvalue++] = KE;
    else if (strcmp(arg[iarg],"temp") == 0) which[nvalue++] = TEMPERATURE;
    else if (strcmp(arg[iarg],"erot") == 0) which[nvalue++] = EROT;
    else if (strcmp(arg[iarg],"trot") == 0) which[nvalue++] = TROT;
    else if (strcmp(arg[iarg],"evib") == 0) which[nvalue++] = EVIB;
    else if (strcmp(arg[iarg],"tvib") == 0) which[nvalue++] = TVIB;
    else error->all(FLERR,"Illegal compute grid command");
    iarg++;
  }

  ngroup = particle->mixture[imix]->ngroup;
  ntotal = ngroup*nvalue;

  per_grid_flag = 1;
  size_per_grid_cols = ntotal;
  post_process_grid_flag = 1;

  nglocal = 0;
  array_grid = NULL;

  // norm vectors

  for (int i = 0; i < nvalue; i++) {
    if (which[i] == NUM || which[i] == NRHO) 
      norm_style[i] = NONE;
    else if (which[i] == MASS) 
      norm_style[i] = COUNT;
    else if (which[i] == U || which[i] == V || which[i] == W) 
      norm_style[i] = MASSWT;
    else if (which[i] == USQ || which[i] == VSQ || which[i] == WSQ) 
      norm_style[i] = MASSWT;
    else if (which[i] == KE)
      norm_style[i] = COUNT;
    else if (which[i] == TEMPERATURE) 
      norm_style[i] = COUNT;
    else if (which[i] == EROT)
      norm_style[i] = COUNT;
    else if (which[i] == TROT) 
      norm_style[i] = RDOF;
    else if (which[i] == EVIB)
      norm_style[i] = COUNT;
    else if (which[i] == TVIB) 
      norm_style[i] = VDOF;
  }

  norm_count = new double*[ngroup];
  norm_mass = new double*[ngroup];
  norm_rdof = new double*[ngroup];
  norm_vdof = new double*[ngroup];
  for (int i = 0; i < ngroup; i++)
    norm_count[i] = norm_mass[i] = norm_rdof[i] = norm_vdof[i] = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeGrid::~ComputeGrid()
{
  delete [] which;
  delete [] norm_style;

  memory->destroy(array_grid);

  for (int i = 0; i < ngroup; i++) {
    memory->destroy(norm_count[i]);
    memory->destroy(norm_mass[i]);
    memory->destroy(norm_rdof[i]);
    memory->destroy(norm_vdof[i]);
  }
  delete [] norm_count;
  delete [] norm_mass;
  delete [] norm_rdof;
  delete [] norm_vdof;
}

/* ---------------------------------------------------------------------- */

void ComputeGrid::init()
{
  if (ngroup != particle->mixture[imix]->ngroup)
    error->all(FLERR,"Number of groups in compute grid mixture has changed");

  reallocate();
}

/* ---------------------------------------------------------------------- */

void ComputeGrid::compute_per_grid()
{
  invoked_per_grid = update->ntimestep;

  Grid::ChildInfo *cinfo = grid->cinfo;
  Particle::Species *species = particle->species;
  Particle::OnePart *particles = particle->particles;
  int *s2g = particle->mixture[imix]->species2group;
  int nlocal = particle->nlocal;

  double fnum = update->fnum;
  double mvv2e = update->mvv2e;
  double tprefactor = mvv2e / (3.0*update->boltz);
  double trotprefactor = mvv2e * 2.0 / update->boltz;
  double tvibprefactor = mvv2e * 2.0 / update->boltz;

  int i,j,k,m,n,ispecies,igroup,icell;
  double mass;
  double *norm,*v,*vec;
  double sum,sumv,sumu,sumw,sume;

  // zero accumulator array and norm vectors

  for (i = 0; i < nglocal; i++)
    for (j = 0; j < ntotal; j++)
      array_grid[i][j] = 0.0;

  for (j = 0; j < ngroup; j++) {
    if (norm = norm_count[j])
      for (i = 0; i < nglocal; i++) norm[i] = 0.0;
    if (norm = norm_mass[j])
      for (i = 0; i < nglocal; i++) norm[i] = 0.0;
    if (norm = norm_rdof[j])
      for (i = 0; i < nglocal; i++) norm[i] = 0.0;
    if (norm = norm_vdof[j])
      for (i = 0; i < nglocal; i++) norm[i] = 0.0;
  }

  // loop over all particles, skip species not in mixture group
  // tally any norm associated with group into norms
  // tally all values associated with group into array_grid

  for (i = 0; i < nlocal; i++) {
    ispecies = particles[i].ispecies;
    igroup = s2g[ispecies];
    if (igroup < 0) continue;
    icell = particles[i].icell;

    mass = species[ispecies].mass;
    if (norm_mass[igroup]) norm_mass[igroup][icell] += mass;
    if (norm_count[igroup]) norm_count[igroup][icell] += 1.0;

    v = particles[i].v;

    vec = array_grid[icell];
    k = igroup*nvalue;

    for (m = 0; m < nvalue; m++) {
      switch (which[m]) {
      case NUM:
        vec[k++] += 1.0;
        break;
      case NRHO:
        vec[k++] += fnum * cinfo[icell].weight / cinfo[icell].volume;
        break;
      case MASS:
        vec[k++] += mass;
        break;
      case U:
        vec[k++] += mass*v[0];
        break;
      case V:
        vec[k++] += mass*v[1];
        break;
      case W:
        vec[k++] += mass*v[2];
        break;
      case USQ:
        vec[k++] += mass*v[0]*v[0];
        break;
      case VSQ:
        vec[k++] += mass*v[1]*v[1];
        break;
      case WSQ:
        vec[k++] += mass*v[2]*v[2];
        break;
      case KE:
        vec[k++] += 0.5*mvv2e*mass * (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
        break;
      case TEMPERATURE:
        vec[k++] += tprefactor*mass * (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
        break;
      case EROT:
        vec[k++] += particles[i].erot;
        break;
      case TROT:
        vec[k++] += trotprefactor*particles[i].erot;
        norm_rdof[igroup][icell] += species[ispecies].rotdof;
        break;
      case EVIB:
        vec[k++] += particles[i].evib;
        break;
      case TVIB:
        vec[k++] += tvibprefactor*particles[i].evib;
        norm_vdof[igroup][icell] += species[ispecies].vibdof;
        break;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   use tallied info to compute normalized values
   icell = -1, return values for entire group = index
     store them in out vector with nstride
   icell >= 0, return value for single icell in group = index
     store it in out[0]
   index = which column of this compute's output is requested
     0 = vector, 1-N = columns of array
   called by dumps with NULL input arrays, so use internal array/norm as input
   called by fix ave/grid with arrays it accumulated over many timesteps
------------------------------------------------------------------------- */

void ComputeGrid::post_process_grid(void *innumer, void *indenom,
                                    int icell, int index,
                                    double *out, int nstride)
{
  // use internal storage for input
  // should never be called with non-NULL inputs
  // fix ave/grid does its own normalization for this compute

  double **array,*norm;

  if (innumer == NULL) {
    array = array_grid;
    norm = normptr(index);
  } else error->all(FLERR,"Invalid call to ComputeGrid::post_process_grid()");

  // request for either a single value or entire column of values
  // for single value, iterate thru outer loop just once

  index--;
  int istart = icell;
  if (icell < 0) istart = 0;

  int m = 0;

  for (int i = istart; i < nglocal; i++) {
    if (!norm) out[m] = array[i][index];
    else if (norm[i] == 0.0) out[m] = 0.0;
    else out[m] = array[i][index]/norm[i];
    m += nstride;
    if (icell >= 0) return;
  }
}

/* ----------------------------------------------------------------------
   reallocate arrays if nglocal has changed
   called by init() and load balancer
------------------------------------------------------------------------- */

void ComputeGrid::reallocate()
{
  if (grid->nlocal == nglocal) return;

  nglocal = grid->nlocal;
  memory->destroy(array_grid);
  memory->create(array_grid,nglocal,ntotal,"grid:array_grid");

  for (int i = 0; i < ngroup; i++) {
    memory->destroy(norm_count[i]);
    memory->destroy(norm_mass[i]);
    memory->destroy(norm_rdof[i]);
    memory->destroy(norm_vdof[i]);
    norm_count[i] = norm_mass[i] = norm_rdof[i] = norm_vdof[i] = NULL;
    for (int j = 0; j < nvalue; j++) {
      if (norm_style[j] == COUNT && norm_count[i] == NULL)
        memory->create(norm_count[i],nglocal,"grid:norm_count");
      else if (norm_style[j] == MASSWT && norm_mass[i] == NULL)
        memory->create(norm_mass[i],nglocal,"grid:norm_mass");
      else if (norm_style[j] == RDOF && norm_rdof[i] == NULL)
        memory->create(norm_rdof[i],nglocal,"grid:norm_dof");
      else if (norm_style[j] == VDOF && norm_vdof[i] == NULL)
        memory->create(norm_vdof[i],nglocal,"grid:norm_dof");
    }
  }
}

/* ----------------------------------------------------------------------
   return info for norm vector used by column N
   input N is value from 1 to Ncols
   output: istyle = NONE, COUNT, etc
   output: igroup = which group within style
------------------------------------------------------------------------- */

void ComputeGrid::normwhich(int n, int &istyle, int &igroup)
{
  igroup = (n-1) / nvalue;
  int ivalue = (n-1) % nvalue;
  if (norm_style[ivalue] == COUNT) istyle = COUNT;
  else if (norm_style[ivalue] == MASSWT) istyle = MASSWT;
  else if (norm_style[ivalue] == RDOF) istyle = RDOF;
  else if (norm_style[ivalue] == VDOF) istyle = VDOF;
  else istyle = NONE;
}

/* ----------------------------------------------------------------------
   return ptr to norm vector used by column N
   input N is value from 1 to Ncols
------------------------------------------------------------------------- */

double *ComputeGrid::normptr(int n)
{
  int igroup = (n-1) / nvalue;
  int ivalue = (n-1) % nvalue;
  if (norm_style[ivalue] == COUNT) return norm_count[igroup];
  if (norm_style[ivalue] == MASSWT) return norm_mass[igroup];
  if (norm_style[ivalue] == RDOF) return norm_rdof[igroup];
  if (norm_style[ivalue] == VDOF) return norm_vdof[igroup];
  return NULL;
}

/* ----------------------------------------------------------------------
   memory usage of local grid-based array
------------------------------------------------------------------------- */

bigint ComputeGrid::memory_usage()
{
  bigint bytes;
  bytes = ntotal*nglocal * sizeof(double);
  for (int i = 0; i < ngroup; i++) {
    if (norm_count[i]) bytes += nglocal * sizeof(double);
    if (norm_mass[i]) bytes += nglocal * sizeof(double);
    if (norm_vdof[i]) bytes += nglocal * sizeof(double);
    if (norm_rdof[i]) bytes += nglocal * sizeof(double);
  }
  return bytes;
}

/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate alyzer
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
#include "compute_grid.h"
#include "particle.h"
#include "mixture.h"
#include "grid.h"
#include "update.h"
#include "modify.h"
#include "memory.h"
#include "error.h"

// DEBUG
#include "comm.h"

using namespace SPARTA_NS;

enum{NUM,U,V,W,USQ,VSQ,WSQ,KE,MASS,TEMPERATURE};
enum{NONE,COUNT,MASSWT};

/* ---------------------------------------------------------------------- */

ComputeGrid::ComputeGrid(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute grid command");

  imix = particle->find_mixture(arg[2]);
  if (imix < 0) error->all(FLERR,"Compute grid mixture ID does not exist");

  nvalue = narg - 3;
  which = new int[nvalue];

  nvalue = 0;
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"n") == 0) which[nvalue++] = NUM;
    else if (strcmp(arg[iarg],"u") == 0) which[nvalue++] = U;
    else if (strcmp(arg[iarg],"v") == 0) which[nvalue++] = V;
    else if (strcmp(arg[iarg],"w") == 0) which[nvalue++] = W;
    else if (strcmp(arg[iarg],"usq") == 0) which[nvalue++] = USQ;
    else if (strcmp(arg[iarg],"vsq") == 0) which[nvalue++] = VSQ;
    else if (strcmp(arg[iarg],"wsq") == 0) which[nvalue++] = WSQ;
    else if (strcmp(arg[iarg],"ke") == 0) which[nvalue++] = KE;
    else if (strcmp(arg[iarg],"mass") == 0) which[nvalue++] = MASS;
    else if (strcmp(arg[iarg],"temp") == 0) which[nvalue++] = TEMPERATURE;
    else error->all(FLERR,"Illegal compute grid command");
    iarg++;
  }

  per_grid_flag = 1;
  ngroup = particle->mixture[imix]->ngroup;
  ntotal = ngroup*nvalue;
  size_per_grid_cols = ntotal;

  nglocal = 0;
  array_grid = NULL;

  // norm vectors

  memory->create(value_norm_style,ngroup,nvalue,"grid:value_norm_style");
  for (int i = 0; i < ngroup; i++)
    for (int j = 0; j < nvalue; j++) {
      if (which[j] == NUM) value_norm_style[i][j] = NONE;
      else if (which[j] == KE || which[j] == MASS || which[j] == TEMPERATURE)
        value_norm_style[i][j] = COUNT; 
      else value_norm_style[i][j] = MASSWT;
    }

  norm_count = new double*[ngroup];
  norm_mass = new double*[ngroup];
  for (int i = 0; i < ngroup; i++)
    norm_count[i] = norm_mass[i] = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeGrid::~ComputeGrid()
{
  delete [] which;

  memory->destroy(array_grid);

  memory->destroy(value_norm_style);

  for (int i = 0; i < ngroup; i++) {
    memory->destroy(norm_count[i]);
    memory->destroy(norm_mass[i]);
  }
  delete [] norm_count;
  delete [] norm_mass;
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

  // compute kinetic energies for each group in each grid cell

  Grid::ChildInfo *cinfo = grid->cinfo;
  Particle::Species *species = particle->species;
  Particle::OnePart *particles = particle->particles;
  int *s2g = particle->mixture[imix]->species2group;
  int nlocal = particle->nlocal;
  double mvv2e = update->mvv2e;
  double tprefactor = mvv2e / (3.0*update->boltz);

  int i,j,k,m,n,ispecies,igroup,icell;
  double mass;
  double *norm,*v,*vec;
  double sum,sumv,sumu,sumw,sume;

  // zero accumulator array and norm vectors

  for (i = 0; i < nglocal; i++)
    for (j = 0; j < ntotal; j++) array_grid[i][j] = 0.0;

  for (j = 0; j < ngroup; j++) {
    if (norm = norm_count[j])
      for (i = 0; i < nglocal; i++) norm[i] = 0.0;
    if (norm = norm_mass[j])
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
      case MASS:
        vec[k++] += mass;
        break;
      case TEMPERATURE:
        vec[k++] += tprefactor*mass * (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
        break;
      }
    }
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
    norm_count[i] = norm_mass[i] = NULL;
    for (int j = 0; j < nvalue; j++) {
      if (which[j] == KE || which[j] == MASS || which[j] == TEMPERATURE) {
        if (norm_count[i] == NULL)
          memory->create(norm_count[i],nglocal,"grid:norm_count");
      } else {
        if (norm_mass[i] == NULL)
          memory->create(norm_mass[i],nglocal,"grid:norm_mass");
      }
    }
  }
}

/* ----------------------------------------------------------------------
   return flag for style of norm vector used by column N
   input N is value from 1 to Ncols
------------------------------------------------------------------------- */

int ComputeGrid::normflag(int n)
{
  int igroup = (n-1) / nvalue;
  int ivalue = (n-1) % nvalue;
  if (value_norm_style[igroup][ivalue] == COUNT) return COUNT;
  if (value_norm_style[igroup][ivalue] == MASSWT) return MASSWT;
  return NONE;
}

/* ----------------------------------------------------------------------
   return ptr to norm vector used by column N
   input N is value from 1 to Ncols
------------------------------------------------------------------------- */

double *ComputeGrid::normptr(int n)
{
  int igroup = (n-1) / nvalue;
  int ivalue = (n-1) % nvalue;
  if (value_norm_style[igroup][ivalue] == COUNT) return norm_count[igroup];
  if (value_norm_style[igroup][ivalue] == MASSWT) return norm_mass[igroup];
  return NULL;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

bigint ComputeGrid::memory_usage()
{
  bigint bytes;
  bytes = ntotal*nglocal * sizeof(double);
  for (int i = 0; i < ngroup; i++) {
    if (norm_count[i]) bytes += nglocal * sizeof(double);
    if (norm_mass[i]) bytes += nglocal * sizeof(double);
  }
  return bytes;
}

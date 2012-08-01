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
#include "compute_grid.h"
#include "particle.h"
#include "mixture.h"
#include "grid.h"
#include "update.h"
#include "modify.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{NUM,U,V,W,USQ,VSQ,WSQ,KE,TEMP};

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
    else if (strcmp(arg[iarg],"temp") == 0) which[nvalue++] = TEMP;
    else error->all(FLERR,"Illegal compute grid command");
  }

  per_grid_flag = 1;
  ngroup = particle->mixture[imix]->ngroup;
  ntotal = ngroup*nvalue;
  if (ntotal == 0) size_per_grid_cols = 0;
  else size_per_grid_cols = ntotal;

  vector_grid = NULL;
  array_grid = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeGrid::~ComputeGrid()
{
  delete [] which;
  memory->destroy(vector_grid);
  memory->destroy(array_grid);
}

/* ---------------------------------------------------------------------- */

void ComputeGrid::init()
{
  if (ngroup != particle->mixture[imix]->ngroup)
    error->all(FLERR,"Number of groups in compute grid mixture has changed");

  // one-time allocation

  if (ntotal == 1) {
    if (vector_grid == NULL)
      memory->create(vector_grid,grid->nlocal,"grid:vector_grid");
  } else {
    if (array_grid == NULL)
      memory->create(array_grid,grid->nlocal,ntotal,"grid:array_grid");
  }
}

/* ---------------------------------------------------------------------- */

void ComputeGrid::compute_per_grid()
{
  invoked_per_grid = update->ntimestep;

  // compute kinetic energies for each group in each grid cell

  Grid::OneCell *cells = grid->cells;
  Particle::Species *species = particle->species;
  Particle::OnePart *particles = particle->particles;
  int *s2g = particle->mixture[imix]->species2group;
  int nlocal = particle->nlocal;
  int nglocal = grid->nlocal;
  double mvv2e = update->mvv2e;

  int i,j,k,m,n,ispecies,igroup,ilocal;
  double *v;

  if (ntotal == 1) {
    for (i = 0; i < nglocal; i++) vector_grid[i] = 0.0;
  
    for (i = 0; i < nlocal; i++) {
      ispecies = particles[i].ispecies;
      if (s2g[ispecies] < 0) continue;

      ilocal = cells[particles[i].icell].local;
      v = particles[i].v;

      switch (which[0]) {
      case NUM:
        vector_grid[ilocal] += 1.0;
        break;
      case U:
        vector_grid[ilocal] += v[0];
        break;
      case V:
        vector_grid[ilocal] += v[1];
        break;
      case W:
        vector_grid[ilocal] += v[2];
        break;
      case USQ:
        vector_grid[ilocal] += v[0]*v[0];
        break;
      case VSQ:
        vector_grid[ilocal] += v[1]*v[1];
        break;
      case WSQ:
        vector_grid[ilocal] += v[2]*v[2];
        break;
      case KE:
        vector_grid[ilocal] += 
          0.5 * mvv2e * species[ispecies].mass * 
          (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        break;
      case TEMP:
        vector_grid[ilocal] += 
          0.5 * mvv2e * species[ispecies].mass * 
          (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        break;
      }
    }

  } else {
    for (i = 0; i < nglocal; i++)
      for (j = 0; j < ntotal; j++) array_grid[i][j] = 0.0;
  
    for (i = 0; i < nlocal; i++) {
      ispecies = particles[i].ispecies;
      igroup = s2g[ispecies];
      if (igroup < 0) continue;

      ilocal = cells[particles[i].icell].local;
      v = particles[i].v;
      k = igroup*nvalue;

      for (m = 0; m < nvalue; m++) {
        switch (which[m]) {
        case NUM:
          array_grid[ilocal][k++] += 1.0;
          break;
        case U:
          array_grid[ilocal][k++] += v[0];
          break;
        case V:
          array_grid[ilocal][k++] += v[1];
          break;
        case W:
          array_grid[ilocal][k++] += v[2];
          break;
        case USQ:
          array_grid[ilocal][k++] += v[0]*v[0];
          break;
        case VSQ:
          array_grid[ilocal][k++] += v[1]*v[1];
          break;
        case WSQ:
          array_grid[ilocal][k++] += v[2]*v[2];
          break;
        case KE:
          array_grid[ilocal][k++] += 
            0.5 * mvv2e * species[ispecies].mass * 
            (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
          break;
        case TEMP:
          array_grid[ilocal][k++] += 
            0.5 * mvv2e * species[ispecies].mass * 
            (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
          break;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

bigint ComputeGrid::memory_usage()
{
  bigint bytes;
  if (ntotal == 1) bytes = grid->nlocal * sizeof(double);
  else bytes = ntotal*grid->nlocal * sizeof(double);
  // extra mem for norm factors
  return bytes;
}

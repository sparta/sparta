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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "compute_sonine_grid.h"
#include "particle.h"
#include "mixture.h"
#include "grid.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{THERMAL,AMOM,BMOM};
enum{X,Y,Z};

/* ---------------------------------------------------------------------- */

ComputeSonineGrid::ComputeSonineGrid(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute sonine/grid command");

  imix = particle->find_mixture(arg[2]);
  if (imix < 0) 
    error->all(FLERR,"Compute sonine/grid mixture ID does not exist");

  int nmax = narg-3;
  which = new int[nmax];
  moment = new int[nmax];
  order = new int[nmax];
  nvalues = 0;
  npergroup = 0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"thermal") == 0) {
      which[nvalues] = THERMAL;
      nvalues++;
      npergroup++;
      iarg++;
    } else if (strcmp(arg[iarg],"a") == 0) {
      if (iarg+3 > narg) 
        error->all(FLERR,"Illegal compute sonine/grid command");
      which[nvalues] = AMOM;
      if (strcmp(arg[iarg+1],"x") == 0) moment[nvalues] = X;
      else if (strcmp(arg[iarg+1],"y") == 0) moment[nvalues] = Y;
      else if (strcmp(arg[iarg+1],"z") == 0) moment[nvalues] = Z;
      else error->all(FLERR,"Illegal compute sonine/grid command");
      order[nvalues] = atoi(arg[iarg+2]);
      if (order[nvalues] < 1 || order[nvalues] > 5)
        error->all(FLERR,"Illegal compute sonine/grid command");
      nvalues++;
      npergroup += order[nvalues];
      iarg += 3;
    } else if (strcmp(arg[iarg],"b") == 0) {
      if (iarg+3 > narg) 
        error->all(FLERR,"Illegal compute sonine/grid command");
      which[nvalues] = BMOM;
      if (strcmp(arg[iarg+1],"xx") == 0) moment[nvalues] = 3*X + X;
      else if (strcmp(arg[iarg+1],"yy") == 0) moment[nvalues] = 3*Y + Y;
      else if (strcmp(arg[iarg+1],"zz") == 0) moment[nvalues] = 3*Z + Z;
      else if (strcmp(arg[iarg+1],"xy") == 0) moment[nvalues] = 3*X + Y;
      else if (strcmp(arg[iarg+1],"yz") == 0) moment[nvalues] = 3*Y + Z;
      else if (strcmp(arg[iarg+1],"xz") == 0) moment[nvalues] = 3*X + Z;
      else error->all(FLERR,"Illegal compute sonine/grid command");
      order[nvalues] = atoi(arg[iarg+2]);
      if (order[nvalues] < 1 || order[nvalues] > 5)
        error->all(FLERR,"Illegal compute sonine/grid command");
      nvalues++;
      npergroup += order[nvalues];
      iarg += 3;
    } else error->all(FLERR,"Illegal compute sonine/grid command");
  }

  per_grid_flag = 1;
  ngroups = particle->mixture[imix]->ngroups;
  ntotal = ngroups*npergroup;
  if (ntotal == 0) size_per_grid_cols = 0;
  else size_per_grid_cols = ntotal;

  vave = NULL;
  avave = NULL;
  count = NULL;
  acount = NULL;
  sonine_vector = NULL;
  sonine_array = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeSonineGrid::~ComputeSonineGrid()
{
  memory->destroy(vave);
  memory->destroy(avave);
  memory->destroy(count);
  memory->destroy(acount);
  memory->destroy(sonine_vector);
  memory->destroy(sonine_array);
}

/* ---------------------------------------------------------------------- */

void ComputeSonineGrid::init()
{
  if (ngroups != particle->mixture[imix]->ngroups)
    error->all(FLERR,"Number of groups in compute ke/grid mixture has changed");

  // one-time allocation

  if (ntotal == 1) {
    if (sonine_vector == NULL) {
      memory->create(vave,grid->nlocal,3,"sonine/grid:vave");
      memory->create(count,grid->nlocal,"sonine/grid:count");
      memory->create(sonine_vector,grid->nlocal,"sonine/grid:sonine_vector");
      vector_grid = sonine_vector;
    }
  } else {
    if (sonine_array == NULL) {
      memory->create(avave,grid->nlocal,ngroups,3,"sonine/grid:vave");
      memory->create(acount,grid->nlocal,ngroups,"sonine/grid:count");
      memory->create(sonine_array,grid->nlocal,ntotal,
                     "sonine/grid:sonine_array");
      array_grid = sonine_array;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSonineGrid::compute_per_grid()
{
  invoked_per_grid = update->ntimestep;

  // compute average velocity for each group in each grid cell
  // then thermal temperature and/or sonine moments

  Grid::OneCell *cells = grid->cells;
  Particle::Species *species = particle->species;
  Particle::OnePart *particles = particle->particles;
  int *s2g = particle->mixture[imix]->species2group;
  int nlocal = particle->nlocal;
  int nglocal = grid->nlocal;
  double mvv2e = update->mvv2e;

  int i,j,k,m,n,ispecies,igroup,ilocal;
  double prefactor,csq;
  double *v;
  double vthermal[3];

  if (ntotal == 1) {
    for (i = 0; i < nglocal; i++) {
      sonine_vector[i] = 0.0;
      vave[i][0] = 0.0;
      vave[i][1] = 0.0;
      vave[i][2] = 0.0;
      count[i] = 0;
    }

    for (i = 0; i < nlocal; i++) {
      ispecies = particles[i].ispecies;
      igroup = s2g[ispecies];
      if (igroup < 0) continue;

      ilocal = cells[particles[i].icell].local;
      v = particles[i].v;
      vave[ilocal][0] += v[0];
      vave[ilocal][1] += v[1];
      vave[ilocal][2] += v[2];
      count[ilocal]++;
    }

    for (i = 0; i < nlocal; i++) {
      vave[ilocal][0] /= count[ilocal];
      vave[ilocal][1] /= count[ilocal];
      vave[ilocal][2] /= count[ilocal];
    }

    for (i = 0; i < nlocal; i++) {
      ispecies = particles[i].ispecies;
      igroup = s2g[ispecies];
      if (igroup < 0) continue;

      ilocal = cells[particles[i].icell].local;
      v = particles[i].v;

      vthermal[0] = v[0] - vave[ilocal][0];
      vthermal[1] = v[1] - vave[ilocal][1];
      vthermal[2] = v[2] - vave[ilocal][2];
      csq = vthermal[0]*vthermal[0] + vthermal[1]*vthermal[1] + 
        vthermal[2]*vthermal[2];

      if (which[0] == THERMAL)
        sonine_vector[ilocal] += 0.5*species[ispecies].mass*csq;
      else if (which[0] == AMOM)
        sonine_vector[ilocal] += vthermal[moment[0]]*csq;
      else if (which[0] == BMOM)
        sonine_vector[ilocal] += vthermal[moment[0] / 3] * 
          vthermal[moment[0] % 3] * csq;  
    }

  } else {
    for (i = 0; i < nglocal; i++) {
      for (j = 0; j < ntotal; j++) sonine_array[i][j] = 0.0;
      for (j = 0; j < ngroups; j++) {
        avave[i][j][0] = 0.0;
        avave[i][j][1] = 0.0;
        avave[i][j][2] = 0.0;
        acount[i][j] = 0;
      }
    }

    for (i = 0; i < nlocal; i++) {
      ispecies = particles[i].ispecies;
      igroup = s2g[ispecies];
      if (igroup < 0) continue;

      ilocal = cells[particles[i].icell].local;
      v = particles[i].v;
      avave[ilocal][igroup][0] += v[0];
      avave[ilocal][igroup][1] += v[1];
      avave[ilocal][igroup][2] += v[2];
      acount[ilocal][igroup]++;
    }

    for (i = 0; i < nlocal; i++)
      for (j = 0; j < ngroups; j++) {
        avave[ilocal][igroup][0] /= acount[ilocal][igroup];
        avave[ilocal][igroup][1] /= acount[ilocal][igroup];
        avave[ilocal][igroup][2] /= acount[ilocal][igroup];
      }

    for (i = 0; i < nlocal; i++) {
      ispecies = particles[i].ispecies;
      igroup = s2g[ispecies];
      if (igroup < 0) continue;

      ilocal = cells[particles[i].icell].local;
      v = particles[i].v;

      vthermal[0] = v[0] - avave[ilocal][igroup][0];
      vthermal[1] = v[1] - avave[ilocal][igroup][1];
      vthermal[2] = v[2] - avave[ilocal][igroup][2];
      csq = vthermal[0]*vthermal[0] + vthermal[1]*vthermal[1] + 
        vthermal[2]*vthermal[2];

      k = igroup*npergroup;
      for (m = 0; m < nvalues; m++) {
        if (which[m] == THERMAL) {
          sonine_array[ilocal][k++] += 0.5*species[ispecies].mass*csq;
        } else if (which[m] == AMOM) {
          prefactor = vthermal[moment[m]];
          for (n = 0; n < order[m]; n++)
            sonine_array[ilocal][k++] += prefactor*pow(csq,n);
        } else if (which[m] == BMOM) {
          prefactor = vthermal[moment[m] / 3] * vthermal[moment[m] % 3];
          for (n = 0; n < order[m]; n++)
            sonine_array[ilocal][k++] += prefactor*pow(csq,n);
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

bigint ComputeSonineGrid::memory_usage()
{
  bigint bytes = 0;
  if (ntotal == 1) {
    bytes += grid->nlocal*3 * sizeof(double);
    bytes += grid->nlocal * sizeof(int);
    bytes += grid->nlocal * sizeof(double);
  } else {
    bytes += grid->nlocal*ngroups*3 * sizeof(double);
    bytes += grid->nlocal*ngroups * sizeof(int);
    bytes += grid->nlocal*ntotal * sizeof(double);
  }
  return bytes;
}

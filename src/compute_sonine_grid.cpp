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
enum{NONE,COUNT,MASSWT,TEMPWT};

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
  nvalue = 0;
  npergroup = 0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"thermal") == 0) {
      which[nvalue] = THERMAL;
      nvalue++;
      npergroup++;
      iarg++;
    } else if (strcmp(arg[iarg],"a") == 0) {
      if (iarg+3 > narg) 
        error->all(FLERR,"Illegal compute sonine/grid command");
      which[nvalue] = AMOM;
      if (strcmp(arg[iarg+1],"x") == 0) moment[nvalue] = X;
      else if (strcmp(arg[iarg+1],"y") == 0) moment[nvalue] = Y;
      else if (strcmp(arg[iarg+1],"z") == 0) moment[nvalue] = Z;
      else error->all(FLERR,"Illegal compute sonine/grid command");
      order[nvalue] = atoi(arg[iarg+2]);
      if (order[nvalue] < 1 || order[nvalue] > 5)
        error->all(FLERR,"Illegal compute sonine/grid command");
      npergroup += order[nvalue];
      nvalue++;
      iarg += 3;
    } else if (strcmp(arg[iarg],"b") == 0) {
      if (iarg+3 > narg) 
        error->all(FLERR,"Illegal compute sonine/grid command");
      which[nvalue] = BMOM;
      if (strcmp(arg[iarg+1],"xx") == 0) moment[nvalue] = 3*X + X;
      else if (strcmp(arg[iarg+1],"yy") == 0) moment[nvalue] = 3*Y + Y;
      else if (strcmp(arg[iarg+1],"zz") == 0) moment[nvalue] = 3*Z + Z;
      else if (strcmp(arg[iarg+1],"xy") == 0) moment[nvalue] = 3*X + Y;
      else if (strcmp(arg[iarg+1],"yz") == 0) moment[nvalue] = 3*Y + Z;
      else if (strcmp(arg[iarg+1],"xz") == 0) moment[nvalue] = 3*X + Z;
      else error->all(FLERR,"Illegal compute sonine/grid command");
      order[nvalue] = atoi(arg[iarg+2]);
      if (order[nvalue] < 1 || order[nvalue] > 5)
        error->all(FLERR,"Illegal compute sonine/grid command");
      npergroup += order[nvalue];
      nvalue++;
      iarg += 3;
    } else error->all(FLERR,"Illegal compute sonine/grid command");
  }

  per_grid_flag = 1;
  ngroup = particle->mixture[imix]->ngroup;
  ntotal = ngroup*npergroup;
  size_per_grid_cols = ntotal;

  vcom = NULL;
  masstot = NULL;
  sonine = NULL;

  memory->create(value_norm_style,ngroup,npergroup,
                 "sonine/grid:value_norm_style");
  norm_count = new double*[ngroup];
  norm_mass = new double*[ngroup];
  norm_temp = new double*[ngroup];
}

/* ---------------------------------------------------------------------- */

ComputeSonineGrid::~ComputeSonineGrid()
{
  delete [] which;
  delete [] moment;
  delete [] order;

  memory->destroy(vcom);
  memory->destroy(masstot);
  memory->destroy(sonine);

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

void ComputeSonineGrid::init()
{
  if (ngroup != particle->mixture[imix]->ngroup)
    error->all(FLERR,"Number of groups in compute ke/grid mixture has changed");

  // one-time allocation of accumulators and norms
  // cannot allocate norms until now since depends on group sizes

  if (sonine == NULL) {
    memory->create(vcom,grid->nlocal,ngroup,3,"sonine/grid:vcom");
    memory->create(masstot,grid->nlocal,ngroup,"sonine/grid:masstot");
    memory->create(sonine,grid->nlocal,ntotal,"sonine/grid:sonine");
    array_grid = sonine;

    for (int i = 0; i < ngroup; i++) {
      int m = 0;
      for (int j = 0; j < nvalue; j++) {
        if (which[j] == THERMAL) value_norm_style[i][m++] = TEMPWT; 
        else {
          for (int k = 0; k < order[j]; k++)
            if (particle->mixture[imix]->groupsize[i] == 1)
              value_norm_style[i][m++] = COUNT; 
            else value_norm_style[i][m++] = MASSWT;
        }
      }
    }

    for (int i = 0; i < ngroup; i++) {
      norm_count[i] = norm_mass[i] = norm_temp[i] = NULL;
      for (int j = 0; j < npergroup; j++) {
        if (value_norm_style[i][j] == NONE) continue;
        if (value_norm_style[i][j] == COUNT && norm_count[i] == NULL)
          memory->create(norm_count[i],grid->nlocal,"sonine/grid:norm_count");
        if (value_norm_style[i][j] == MASSWT && norm_mass[i] == NULL)
          memory->create(norm_mass[i],grid->nlocal,"sonine/grid:norm_mass");
        if (value_norm_style[i][j] == TEMPWT && norm_temp[i] == NULL)
          memory->create(norm_temp[i],grid->nlocal,"sonine/grid:norm_temp");
      }
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
  int *groupsize = particle->mixture[imix]->groupsize;
  int nlocal = particle->nlocal;
  int nglocal = grid->nlocal;
  double mvv2e = update->mvv2e;
  double kbwt = 3.0*update->boltz;

  int i,j,k,m,n,ispecies,igroup,ilocal;
  double prefactor,csq,wt;
  double *norm,*v;
  double vthermal[3];

  // compute COM velocity for each cell and group

  for (i = 0; i < nglocal; i++) {
    for (j = 0; j < ngroup; j++) {
      vcom[i][j][0] = 0.0;
      vcom[i][j][1] = 0.0;
      vcom[i][j][2] = 0.0;
      masstot[i][j] = 0.0;
    }
  }

  for (i = 0; i < nlocal; i++) {
    ispecies = particles[i].ispecies;
    igroup = s2g[ispecies];
    if (igroup < 0) continue;
    ilocal = cells[particles[i].icell].local;

    if (groupsize[igroup] == 1) wt = 1.0;
    else wt = species[ispecies].mass;
    masstot[ilocal][igroup] += wt;

    v = particles[i].v;
    vcom[ilocal][igroup][0] += wt * v[0];
    vcom[ilocal][igroup][1] += wt * v[1];
    vcom[ilocal][igroup][2] += wt * v[2];
  }

  for (i = 0; i < nglocal; i++)
    for (j = 0; j < ngroup; j++) {
      if (masstot[i][j] == 0.0) continue;
      vcom[i][j][0] /= masstot[i][j];
      vcom[i][j][1] /= masstot[i][j];
      vcom[i][j][2] /= masstot[i][j];
    }
  
  // compute thermal temperature and sonine moments

  for (i = 0; i < nglocal; i++)
    for (j = 0; j < ntotal; j++) sonine[i][j] = 0.0;

  for (j = 0; j < ngroup; j++) {
    if (norm = norm_count[j])
      for (i = 0; i < nglocal; i++) norm[i] = 0.0;
    if (norm = norm_mass[j])
      for (i = 0; i < nglocal; i++) norm[i] = 0.0;
    if (norm = norm_temp[j])
      for (i = 0; i < nglocal; i++) norm[i] = 0.0;
  }

  for (i = 0; i < nlocal; i++) {
    ispecies = particles[i].ispecies;
    igroup = s2g[ispecies];
    if (igroup < 0) continue;
    ilocal = cells[particles[i].icell].local;

    if (norm_mass[igroup]) wt = species[ispecies].mass;
    else wt = 1.0;
    if (norm_count[igroup]) norm_count[igroup][ilocal] += 1.0;
    if (norm_mass[igroup]) norm_mass[igroup][ilocal] += wt;

    v = particles[i].v;
    vthermal[0] = v[0] - vcom[ilocal][igroup][0];
    vthermal[1] = v[1] - vcom[ilocal][igroup][1];
    vthermal[2] = v[2] - vcom[ilocal][igroup][2];
    csq = vthermal[0]*vthermal[0] + vthermal[1]*vthermal[1] + 
      vthermal[2]*vthermal[2];

    k = igroup*npergroup;
    
    for (m = 0; m < nvalue; m++) {
      switch (which[m]) {
      case THERMAL:
        sonine[ilocal][k++] += 0.5*species[ispecies].mass*csq;
        norm_temp[igroup][ilocal] += kbwt;
        break;
      case AMOM:
        sonine[ilocal][k++] = wt * vthermal[moment[m]];
        for (n = 1; n < order[m]; n++) {
          sonine[ilocal][k] = csq*sonine[ilocal][k-1];
          k++;
        }
        break;
      case BMOM:
        sonine[ilocal][k++] = wt * vthermal[moment[m] / 3] * 
          vthermal[moment[m] % 3];
        for (n = 1; n < order[m]; n++) {
          sonine[ilocal][k] = csq*sonine[ilocal][k-1];
          k++;
        }
        break;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   return ptr to norm vector used by column N
------------------------------------------------------------------------- */

double *ComputeSonineGrid::normptr(int n)
{
  int igroup = n / npergroup;
  int ivalue = n % npergroup;
  if (value_norm_style[igroup][ivalue] == COUNT) return norm_count[igroup];
  if (value_norm_style[igroup][ivalue] == MASSWT) return norm_mass[igroup];
  if (value_norm_style[igroup][ivalue] == TEMPWT) return norm_temp[igroup];
  return NULL;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

bigint ComputeSonineGrid::memory_usage()
{
  bigint bytes = 0;
  bytes += grid->nlocal*ngroup*3 * sizeof(double);
  bytes += grid->nlocal*ngroup * sizeof(int);
  bytes += grid->nlocal*ntotal * sizeof(double);
  for (int i = 0; i < ngroup; i++) {
    if (norm_count[i]) bytes += grid->nlocal * sizeof(double);
    if (norm_mass[i]) bytes += grid->nlocal * sizeof(double);
    if (norm_temp[i]) bytes += grid->nlocal * sizeof(double);
  }
  return bytes;
}

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
enum{NONE,COUNT,MASSWT};

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
      nvalue++;
      npergroup += order[nvalue];
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
      nvalue++;
      npergroup += order[nvalue];
      iarg += 3;
    } else error->all(FLERR,"Illegal compute sonine/grid command");
  }

  per_grid_flag = 1;
  ngroup = particle->mixture[imix]->ngroup;
  ntotal = ngroup*npergroup;
  size_per_grid_cols = ntotal;

  vave = NULL;
  sonine = NULL;

  memory->create(group_norm_style,ngroup,"sonine/grid:group_norm_style");
  norms = new double*[ngroup];
}

/* ---------------------------------------------------------------------- */

ComputeSonineGrid::~ComputeSonineGrid()
{
  memory->destroy(vave);
  memory->destroy(sonine);

  memory->destroy(group_norm_style);
  for (int i = 0; i < ngroup; i++)
    memory->destroy(norms[i]);
  delete [] norms;
}

/* ---------------------------------------------------------------------- */

void ComputeSonineGrid::init()
{
  if (ngroup != particle->mixture[imix]->ngroup)
    error->all(FLERR,"Number of groups in compute ke/grid mixture has changed");

  // one-time allocation of accumulators and norms
  // cannot allocate norms until now since depends on group sizes

  if (sonine == NULL) {
    memory->create(vave,grid->nlocal,ngroup,3,"sonine/grid:vave");
    memory->create(sonine,grid->nlocal,ntotal,"sonine/grid:sonine");
    array_grid = sonine;

    for (int i = 0; i < ngroup; i++) {
      if (particle->mixture[imix]->groupsize[i] == 1)
        group_norm_style[i] = COUNT; 
      else group_norm_style[i] = MASSWT;
      memory->create(norms[i],grid->nlocal,"sonine/grid:norms");
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
  double *norm,*v;
  double vthermal[3];

  for (i = 0; i < nglocal; i++) {
    for (j = 0; j < ntotal; j++) sonine[i][j] = 0.0;
    for (j = 0; j < ngroup; j++) {
      vave[i][j][0] = 0.0;
      vave[i][j][1] = 0.0;
      vave[i][j][2] = 0.0;
      count[i][j] = 0.0;
    }
  }

  for (j = 0; j < ngroup; j++) {
    norm = norms[j];
    for (i = 0; i < nglocal; i++) norm[i] = 0.0;
  }

  for (i = 0; i < nlocal; i++) {
    ispecies = particles[i].ispecies;
    igroup = s2g[ispecies];
    if (igroup < 0) continue;
    norm = norms[igroup];

    ilocal = cells[particles[i].icell].local;
    v = particles[i].v;
    vave[ilocal][igroup][0] += v[0];
    vave[ilocal][igroup][1] += v[1];
    vave[ilocal][igroup][2] += v[2];
    if (group_norm_style[igroup] == COUNT) norm[ilocal] += 1.0;
    else norm[ilocal] += species[particles[i].ispecies].mass;
  }

  for (i = 0; i < nlocal; i++)
    for (j = 0; j < ngroup; j++) {
      // check for div by 0.0
      vave[i][j][0] /= norms[i][j];
      vave[i][j][1] /= norms[i][j];
      vave[i][j][2] /= norms[i][j];
    }
  
  for (i = 0; i < nlocal; i++) {
    ispecies = particles[i].ispecies;
    igroup = s2g[ispecies];
    if (igroup < 0) continue;
    
    ilocal = cells[particles[i].icell].local;
    v = particles[i].v;
    
    // apply mass weighting

    vthermal[0] = v[0] - vave[ilocal][igroup][0];
    vthermal[1] = v[1] - vave[ilocal][igroup][1];
    vthermal[2] = v[2] - vave[ilocal][igroup][2];
    csq = vthermal[0]*vthermal[0] + vthermal[1]*vthermal[1] + 
      vthermal[2]*vthermal[2];
    
    k = igroup*npergroup;
    for (m = 0; m < nvalue; m++) {
      if (which[m] == THERMAL) {
        sonine[ilocal][k++] += 0.5*species[ispecies].mass*csq;
      } else if (which[m] == AMOM) {
        prefactor = vthermal[moment[m]];
        for (n = 0; n < order[m]; n++)
          sonine[ilocal][k++] += prefactor*pow(csq,n);
      } else if (which[m] == BMOM) {
        prefactor = vthermal[moment[m] / 3] * vthermal[moment[m] % 3];
        for (n = 0; n < order[m]; n++)
          sonine[ilocal][k++] += prefactor*pow(csq,n);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   return ptr to norm vector used by column N
------------------------------------------------------------------------- */

double *ComputeSonineGrid::normptr(int n)
{
  int igroup = n / nvalue;
  return norms[igroup];
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
  return bytes;
}

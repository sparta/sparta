/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate alyzer
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
#include "compute_tvib_grid.h"
#include "particle.h"
#include "mixture.h"
#include "grid.h"
#include "update.h"
#include "modify.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{NONE,COUNT,MASSWT,DOF};

/* ---------------------------------------------------------------------- */

ComputeTvibGrid::ComputeTvibGrid(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute grid command");

  imix = particle->find_mixture(arg[2]);
  if (imix < 0) error->all(FLERR,"Compute grid mixture ID does not exist");

  ngroup = particle->mixture[imix]->ngroup;
  nspecies = particle->mixture[imix]->nspecies;

  per_grid_flag = 1;
  size_per_grid_cols = ngroup;
  post_process_grid_flag = 1;
  size_per_grid_extra_cols = nspecies;

  nglocal = 0;
  array_grid = NULL;

  array_grid_extra = NULL;
  norm_grid_extra = NULL;

  double *tspecies = new double[nspecies];
  double *g2s = new double[nspecies];
}

/* ---------------------------------------------------------------------- */

ComputeTvibGrid::~ComputeTvibGrid()
{
  memory->destroy(array_grid);

  memory->destroy(array_grid_extra);
  memory->destroy(norm_grid_extra);

  delete [] tspecies;
  delete [] g2s;
}

/* ---------------------------------------------------------------------- */

void ComputeTvibGrid::init()
{
  if (ngroup != particle->mixture[imix]->ngroup)
    error->all(FLERR,"Number of groups in compute tvib/grid "
               "mixture has changed");
  if (nspecies != particle->mixture[imix]->nspecies)
    error->all(FLERR,"Number of species in compute tvib/grid "
               "mixture has changed");

  reallocate();
}

/* ---------------------------------------------------------------------- */

void ComputeTvibGrid::compute_per_grid()
{
  invoked_per_grid = update->ntimestep;

  Particle::Species *species = particle->species;
  Particle::OnePart *particles = particle->particles;
  int *s2g = particle->mixture[imix]->species2group;
  int *s2s = particle->mixture[imix]->species2species;
  int nlocal = particle->nlocal;

  int i,j,ispecies,mixspecies,igroup,icell;
  double *norm;

  // zero extra accumulator array and norm vectors

  for (i = 0; i < nglocal; i++)
    for (j = 0; j < nspecies; j++) {
      array_grid_extra[i][j] = 0.0;
      norm_grid_extra[i][j] = 0.0;
    }

  // loop over all particles, skip species not in mixture group
  // tally vibrational energy into array
  // tally particle count into norm

  for (i = 0; i < nlocal; i++) {
    ispecies = particles[i].ispecies;
    igroup = s2g[ispecies];
    if (igroup < 0) continue;
    mixspecies = s2s[ispecies];
    icell = particles[i].icell;

    array_grid_extra[icell][mixspecies] += particles[i].evib;
    norm_grid_extra[icell][mixspecies] += 1.0;
  }
}

/* ----------------------------------------------------------------------
   use tallied per-species info to compute a group vibrational temperature
   one temperature for each grid cell, stored in strided out vector
   index = which column of this compute's output is requested
     0 = vector, 1-N = columns of array
   called by dump with NULL input arrays, so use internal array/norm as input
   called by fix ave/grid with arrays it accumulated over many timesteps
------------------------------------------------------------------------- */

void ComputeTvibGrid::post_process_grid(double **inarray, double **innorm,
                                        int index, double *out, int nstride)
{
  // use internal storage for input

  if (inarray == NULL) {
    inarray = array_grid_extra;
    innorm = norm_grid_extra;
  }

  Particle::Species *species = particle->species;
  int *mixspecies = particle->mixture[imix]->species;
  int *m2g = particle->mixture[imix]->mix2group;

  // zero output vector

  int m = 0;
  for (int i = 0; i < nglocal; i++) {
    out[m] = 0.0;
    m += nstride;
  }

  // for each grid cell, compute per-species temperatures for species in group
  // then combine them into requested group temp as weighted average

  int i,j,k;
  double theta,ibar,numer,denom;

  index--;
  int gspecies = 0;
  for (i = 0; i < nspecies; i++)
    if (m2g[i] == index) g2s[gspecies++] = i;

  m = 0;
  for (i = 0; i < nglocal; i++) {
    for (j = 0; j < gspecies; j++) {
      k = g2s[j];
      theta = species[mixspecies[k]].vibtemp;
      if (theta == 0.0 || innorm[i][k] == 0.0) {
        tspecies[j] = 0.0;
        continue;
      }
      ibar = inarray[i][k] / innorm[i][k];
      ibar /= update->boltz * theta;
      if (ibar == 0.0) {
        tspecies[j] = 0.0;
        continue;
      }
      denom = update->boltz * innorm[i][k] * ibar * log(1.0 + 1.0/ibar);
      tspecies[j] = inarray[i][k] / denom;
    }

    numer = denom = 0.0;
    for (j = 0; j < gspecies; j++) {
      k = g2s[j];
      numer += tspecies[j]*innorm[i][k];
      denom += innorm[i][k];
    }

    if (denom == 0.0) out[m] = 0.0;
    else out[m] = numer/denom;
    m += nstride;
  }
}

/* ----------------------------------------------------------------------
   reallocate arrays if nglocal has changed
   called by init() and load balancer
------------------------------------------------------------------------- */

void ComputeTvibGrid::reallocate()
{
  if (grid->nlocal == nglocal) return;

  nglocal = grid->nlocal;
  memory->destroy(array_grid);
  memory->create(array_grid,nglocal,ngroup,"grid:array_grid");

  memory->destroy(array_grid_extra);
  memory->create(array_grid_extra,nglocal,nspecies,"grid:array_grid_extra");
  memory->destroy(norm_grid_extra);
  memory->create(norm_grid_extra,nglocal,nspecies,"grid:norm_grid_extra");
}

/* ----------------------------------------------------------------------
   return info for norm vector used by column N
   input N is value from 1 to Ncols
   output: istyle = NONE, COUNT, etc
   output: igroup = which group within style
------------------------------------------------------------------------- */

void ComputeTvibGrid::normwhich(int n, int &istyle, int &igroup)
{
  igroup = n-1;
  istyle = NONE;
}

/* ----------------------------------------------------------------------
   return ptr to norm vector used by column N
   input N is value from 1 to Ncols
------------------------------------------------------------------------- */

double *ComputeTvibGrid::normptr(int n)
{
  return NULL;
}

/* ----------------------------------------------------------------------
   memory usage of local grid-based arrays
------------------------------------------------------------------------- */

bigint ComputeTvibGrid::memory_usage()
{
  bigint bytes;
  bytes = ngroup*nglocal * sizeof(double);
  for (int i = 0; i < ngroup; i++) bytes += nglocal * sizeof(double);
  bytes += 2*nglocal*nspecies * sizeof(double);
  return bytes;
}

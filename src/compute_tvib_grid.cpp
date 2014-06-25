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

  per_grid_flag = 1;
  ngroup = particle->mixture[imix]->ngroup;
  nspecies = particle->mixture[imix]->nspecies;
  size_per_grid_cols = ngroup;

  nglocal = 0;
  array_grid = NULL;
  array_grid_extra = NULL;

  norm_count = new double*[ngroup];
  for (int i = 0; i < ngroup; i++) norm_count[i] = NULL;

  norm_count_extra = new double*[nspecies];
  for (int i = 0; i < nspecies; i++) norm_count_extra[i] = NULL;

  double *tspecies = new double[nspecies];
}

/* ---------------------------------------------------------------------- */

ComputeTvibGrid::~ComputeTvibGrid()
{
  memory->destroy(array_grid);
  memory->destroy(array_grid_extra);

  for (int i = 0; i < ngroup; i++) memory->destroy(norm_count[i]);
  delete [] norm_count;
  for (int i = 0; i < nspecies; i++) memory->destroy(norm_count_extra[i]);
  delete [] norm_count_extra;

  delete [] tspecies;
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
    for (j = 0; j < nspecies; j++) array_grid_extra[i][j] = 0.0;

  for (j = 0; j < nspecies; j++) {
    norm = norm_count_extra[j];
    for (i = 0; i < nglocal; i++) norm[i] = 0.0;
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

    array_grid_extra[mixspecies][icell] += particles[i].evib;
    norm_count_extra[mixspecies][icell] += 1.0;
  }
}

/* ----------------------------------------------------------------------
   user tallied per-species info to compute per-group vibrational temps
   called by dump with NULL arrays, so use internal arrays
   called by fix ave/grid with arrays it accumulated over many timesteps
------------------------------------------------------------------------- */

void ComputeTvibGrid::post_process_grid(double **one, double **onenorm)
{
  // setup ptrs for operation

  double **innumer,**indenom,**outnumer,**outdenom;

  if (one == NULL) {
    innumer = array_grid_extra;
    indenom = norm_count_extra;
    outnumer = array_grid;
    outdenom = norm_count;
  } else {
    innumer = one;
    indenom = onenorm;
  }

  Particle::Species *species = particle->species;
  int *mixspecies = particle->mixture[imix]->species;
  int *m2g = particle->mixture[imix]->mix2group;

  int i,j,igroup;
  double theta,ibar,denom;
  double *norm;

  // zero output array and norm vectors

  for (i = 0; i < nglocal; i++)
    for (j = 0; j < ngroup; j++) outnumer[i][j] = 0.0;

  for (j = 0; j < ngroup; j++) {
    norm = outdenom[j];
    for (i = 0; i < nglocal; i++) norm[i] = 0.0;
  }

  // for each grid cell, compute per-species temperatures
  // then combine per-species temperatures into groups as weighted averages

  for (i = 0; i < nglocal; i++) {
    for (j = 0; j < nspecies; j++) {
      theta = species[mixspecies[j]].vibtemp;
      if (theta == 0.0 || indenom[i][j] == 0.0) {
        tspecies[j] = 0.0;
        continue;
      }
      ibar = innumer[i][j] / indenom[i][j];
      ibar /= update->boltz * theta;
      if (ibar == 0.0) {
        tspecies[j] = 0.0;
        continue;
      }
      denom = update->boltz * indenom[i][j] * ibar * log(1.0 + 1.0/ibar);
      tspecies[j] = innumer[i][j] / denom;
    }

    for (j = 0; j < nspecies; j++) {
      igroup = m2g[j];
      outnumer[i][igroup] += tspecies[j]*norm_count_extra[i][j];
      outdenom[i][igroup] += norm_count_extra[i][j];
    }
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

  for (int i = 0; i < ngroup; i++) {
    memory->destroy(norm_count[i]);
    memory->create(norm_count[i],nglocal,"grid:norm_count");
  }
  for (int i = 0; i < nspecies; i++) {
    memory->destroy(norm_count_extra[i]);
    memory->create(norm_count_extra[i],nglocal,"grid:norm_count_extra");
  }
}

/* ----------------------------------------------------------------------
   return info for norm vector used by column N
   input N is value from 1 to Ncols
   output: istyle = NONE, COUNT, etc
   output: igroup = which group within style
------------------------------------------------------------------------- */

void ComputeTvibGrid::normwhich(int n, int &istyle, int &igroup)
{
  istyle = COUNT;
  igroup = n-1;
}

/* ----------------------------------------------------------------------
   return ptr to norm vector used by column N
   input N is value from 1 to Ncols
------------------------------------------------------------------------- */

double *ComputeTvibGrid::normptr(int n)
{
  return norm_count[n-1];
}

/* ----------------------------------------------------------------------
   memory usage of local grid-based arrays
------------------------------------------------------------------------- */

bigint ComputeTvibGrid::memory_usage()
{
  bigint bytes;
  bytes = ngroup*nglocal * sizeof(double);
  bytes = nspecies*nglocal * sizeof(double);
  for (int i = 0; i < ngroup; i++) bytes += nglocal * sizeof(double);
  for (int i = 0; i < nspecies; i++) bytes += nglocal * sizeof(double);
  return bytes;
}

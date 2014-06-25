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

  Grid::ChildInfo *cinfo = grid->cinfo;
  Particle::Species *species = particle->species;
  Particle::OnePart *particles = particle->particles;
  int *s2g = particle->mixture[imix]->species2group;
  int *s2s = particle->mixture[imix]->species2species;
  int nlocal = particle->nlocal;

  double tvibprefactor = update->mvv2e * 2.0 / update->boltz;

  int i,j,ispecies,mixspecies,igroup,icell;
  double *norm;

  // zero accumulator array and norm vectors

  for (i = 0; i < nglocal; i++)
    for (j = 0; j < nspecies; j++) array_grid_extra[i][j] = 0.0;

  for (j = 0; j < ngroup; j++) {
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

    array_grid_extra[mixspecies][icell] += tvibprefactor*particles[i].evib;
    norm_count_extra[mixspecies][icell] += 1.0;
  }
}

/* ----------------------------------------------------------------------
   user tallied per-species info to prouduce per-group temperatures
   called by dump with NULL arrays, so use own
   called by fix ave/grid with arrays it accumulated
------------------------------------------------------------------------- */

void ComputeTvibGrid::post_process_grid(double **one, double **onenorm)
{
  double **numerator,**denom;

  if (one == NULL) {
    numerator = array_grid_extra;
    denom = norm_count_extra;
  } else {
    numerator = one;
    denom = onenorm;
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

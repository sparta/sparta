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
  if (narg != 4) error->all(FLERR,"Illegal compute grid command");

  int igroup = grid->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Compute grid group ID does not exist");
  groupbit = grid->bitmask[igroup];

  imix = particle->find_mixture(arg[3]);
  if (imix < 0) error->all(FLERR,"Compute tvib/grid mixture ID does not exist");

  ngroup = particle->mixture[imix]->ngroup;
  mixspecies = particle->mixture[imix]->nspecies;
  nspecies = particle->nspecies;

  per_grid_flag = 1;
  size_per_grid_cols = ngroup;
  post_process_grid_flag = 1;

  // allocate and initialize nmap,map,s2t,t2s
  // must first set groupsize, groupspecies by mixture->init()
  // 2 tally quantities per species across all groups

  particle->mixture[imix]->init();
  int *groupsize = particle->mixture[imix]->groupsize;
  int **groupspecies = particle->mixture[imix]->groupspecies;

  int nmax = 0;
  for (int i = 0; i < ngroup; i++)
    nmax = MAX(nmax,groupsize[i]);

  nmap = new int[ngroup];
  memory->create(map,ngroup,2*nmax,"tvib/grid:map");

  tspecies = new double[nspecies];
  ntotal = 2*mixspecies;

  s2t = new int[nspecies];
  t2s = new int[ntotal];

  for (int isp = 0; isp < nspecies; isp++) s2t[isp] = -1;

  int itally = 0;
  for (int igroup = 0; igroup < ngroup; igroup++) {
    nmap[igroup] = 2*groupsize[igroup];
    for (int n = 0; n < groupsize[igroup]; n++) {
      s2t[groupspecies[igroup][n]] = itally;
      t2s[itally] = groupspecies[igroup][n];
      t2s[itally+1] = groupspecies[igroup][n];
      map[igroup][2*n] = itally;
      map[igroup][2*n+1] = itally+1;
      itally += 2;
    }
  }

  nglocal = 0;
  vector_grid = NULL;
  tally = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeTvibGrid::~ComputeTvibGrid()
{
  delete [] nmap;
  memory->destroy(map);

  delete [] tspecies;
  delete [] s2t;
  delete [] t2s;

  memory->destroy(vector_grid);
  memory->destroy(tally);
}

/* ---------------------------------------------------------------------- */

void ComputeTvibGrid::init()
{
  if (ngroup != particle->mixture[imix]->ngroup)
    error->all(FLERR,"Number of groups in compute tvib/grid "
               "mixture has changed");
  if (mixspecies != particle->mixture[imix]->nspecies)
    error->all(FLERR,"Number of species in compute tvib/grid "
               "mixture has changed");
  if (nspecies != particle->nspecies)
    error->all(FLERR,"Number of total species in compute tvib/grid "
               "has changed");

  int *groupsize = particle->mixture[imix]->groupsize;
  for (int i = 0; i < ngroup; i++)
    if (2*groupsize[i] != nmap[i])
      error->all(FLERR,"Number of species in compute tvib/grid "
                 "group has changed");

  reallocate();
}

/* ---------------------------------------------------------------------- */

void ComputeTvibGrid::compute_per_grid()
{
  invoked_per_grid = update->ntimestep;

  Grid::ChildInfo *cinfo = grid->cinfo;
  Particle::OnePart *particles = particle->particles;
  int *s2g = particle->mixture[imix]->species2group;
  int *s2s = particle->mixture[imix]->species2species;
  int nlocal = particle->nlocal;

  int i,j,ispecies,igroup,icell;

  // zero all accumulators - could do this with memset()

  for (i = 0; i < nglocal; i++)
    for (j = 0; j < ntotal; j++)
      tally[i][j] = 0.0;

  // loop over all particles, skip species not in mixture group
  // tally vibrational energy and particle count for each species

  for (i = 0; i < nlocal; i++) {
    ispecies = particles[i].ispecies;
    igroup = s2g[ispecies];
    if (igroup < 0) continue;
    icell = particles[i].icell;
    if (!(cinfo[icell].mask & groupbit)) continue;

    j = s2t[ispecies];
    tally[icell][j] += particles[i].evib;
    tally[icell][j+1] += 1.0;
  }
}

/* ----------------------------------------------------------------------
   query info about internal tally array for this compute
   index = which column of output (0 for vec, 1 to N for array)
   return # of tally quantities for this index
   also return array = ptr to tally array
   also return cols = ptr to list of columns in tally for this index
------------------------------------------------------------------------- */

int ComputeTvibGrid::query_tally_grid(int index, double **&array, int *&cols)
{
  index--;
  array = tally;
  cols = map[index];
  return nmap[index];
}

/* ----------------------------------------------------------------------
   tally accumulated info to compute final normalized values
   index = which column of output (0 for vec, 1 to N for array)
   for etally = NULL:
     use internal tallied info for single timestep, set nsample = 1
     if onecell = -1, compute values for all grid cells
       store results in vector_grid with nstride = 1 (single col of array_grid)
     if onecell >= 0, compute single value for onecell and return it
   for etaylly = ptr to caller array:
     use external tallied info for many timesteps
     nsample = additional normalization factor used by some values
     emap = list of etally columns to use, # of columns determined by index
     store results in caller's vec, spaced by nstride
   if norm = 0.0, set result to 0.0 directly so do not divide by 0.0
------------------------------------------------------------------------- */

double ComputeTvibGrid::post_process_grid(int index, int onecell, int nsample,
                                          double **etally, int *emap,
                                          double *vec, int nstride)
{
  index--;
  
  int lo = 0;
  int hi = nglocal;
  int k = 0;

  if (!etally) {
    nsample = 1;
    etally = tally;
    emap = map[index];
    vec = vector_grid;
    nstride = 1;
    if (onecell >= 0) {
      lo = onecell;
      hi = lo + 1;
      k = lo;
    }
  }

  // compute normalized final value for each grid cell
  // compute Ibar and Tspecies for each species in the requested group
  // loop over species in group to compute normalized Tgroup

  Particle::Species *species = particle->species;

  int isp,evib,count,ispecies;
  double theta,ibar,numer,denom;

  double boltz = update->boltz;
  int nsp = nmap[index] / 2;

  for (int icell = lo; icell < hi; icell++) {
    evib = emap[0];
    count = evib+1;
    for (isp = 0; isp < nsp; isp++) {
      ispecies = t2s[evib];
      theta = species[ispecies].vibtemp;
      if (theta == 0.0 || etally[icell][count] == 0.0) {
        tspecies[isp] = 0.0;
        continue;
      }
      ibar = etally[icell][evib] / (etally[icell][count] * boltz * theta);
      if (ibar == 0.0) {
        tspecies[isp] = 0.0;
        continue;
      }
      denom = boltz * etally[icell][count] * ibar * log(1.0 + 1.0/ibar);
      tspecies[isp] = etally[icell][evib] / denom;

      evib += 2;
      count = evib+1;
    }

    numer = denom = 0.0;
    count = emap[1];
    for (isp = 0; isp < nsp; isp++) {
      numer += tspecies[isp]*etally[icell][count];
      denom += etally[icell][count];
      count += 2;
    }

    if (denom == 0.0) vec[k] = 0.0;
    else vec[k] = numer/denom;
    k += nstride;
  }

  if (onecell < 0) return 0.0;
  return vec[onecell];
}

/* ----------------------------------------------------------------------
   reallocate arrays if nglocal has changed
   called by init() and whenever grid changes
------------------------------------------------------------------------- */

void ComputeTvibGrid::reallocate()
{
  if (grid->nlocal == nglocal) return;

  memory->destroy(vector_grid);
  memory->destroy(tally);
  nglocal = grid->nlocal;
  memory->create(vector_grid,nglocal,"tvib/grid:vector_grid");
  memory->create(tally,nglocal,ntotal,"tvib/grid:tally");
}

/* ----------------------------------------------------------------------
   memory usage of local grid-based data
------------------------------------------------------------------------- */

bigint ComputeTvibGrid::memory_usage()
{
  bigint bytes;
  bytes = nglocal * sizeof(double);
  bytes = ntotal*nglocal * sizeof(double);
  return bytes;
}

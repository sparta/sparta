/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "compute_thermal_grid.h"
#include "particle.h"
#include "mixture.h"
#include "grid.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

// user keywords

enum{TEMP,PRESS};

/* ---------------------------------------------------------------------- */

ComputeThermalGrid::ComputeThermalGrid(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal compute thermal/grid command");

  int igroup = grid->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Compute grid group ID does not exist");
  groupbit = grid->bitmask[igroup];

  imix = particle->find_mixture(arg[3]);
  if (imix < 0)
    error->all(FLERR,"Compute thermal/grid mixture ID does not exist");
  ngroup = particle->mixture[imix]->ngroup;

  nvalue = narg - 4;
  value = new int[nvalue];

  int ivalue = 0;
  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"temp") == 0) value[ivalue] = TEMP;
    else if (strcmp(arg[iarg],"press") == 0) value[ivalue] = PRESS;
    else error->all(FLERR,"Illegal compute thermal/grid command");
    ivalue++;
    iarg++;
  }

  per_grid_flag = 1;
  size_per_grid_cols = ngroup*nvalue;
  post_process_grid_flag = 1;

  // allocate and initialize nmap and map
  // npergroup = 6 tally quantities per group
  // same tally quantities for all user values

  npergroup = 6;
  ntotal = ngroup*npergroup;

  nmap = new int[nvalue];
  for (int i = 0; i < nvalue; i++) nmap[i] = npergroup;

  memory->create(map,ngroup*nvalue,npergroup,"thermal/grid:map");
  for (int i = 0; i < ngroup*nvalue; i++)
    for (int j = 0; j < npergroup; j++)
      map[i][j] = (i/nvalue)*npergroup + j;

  nglocal = 0;
  vector_grid = NULL;
  tally = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeThermalGrid::~ComputeThermalGrid()
{
  if (copymode) return;

  delete [] value;

  delete [] nmap;
  memory->destroy(map);

  memory->destroy(vector_grid);
  memory->destroy(tally);
}

/* ---------------------------------------------------------------------- */

void ComputeThermalGrid::init()
{
  if (ngroup != particle->mixture[imix]->ngroup)
    error->all(FLERR,"Number of groups in compute thermal/grid "
               "mixture has changed");

  tprefactor = update->mvv2e / (3.0*update->boltz);
  pprefactor = update->fnum * update->mvv2e / 3.0;

  reallocate();
}

/* ---------------------------------------------------------------------- */

void ComputeThermalGrid::compute_per_grid()
{
  invoked_per_grid = update->ntimestep;

  Grid::ChildInfo *cinfo = grid->cinfo;
  Particle::Species *species = particle->species;
  Particle::OnePart *particles = particle->particles;
  int *s2g = particle->mixture[imix]->species2group;
  int nlocal = particle->nlocal;

  int i,j,k,ispecies,igroup,icell;
  double mass;
  double *v,*vec;

  // zero all accumulators - could do this with memset()

  for (i = 0; i < nglocal; i++)
    for (j = 0; j < ntotal; j++)
      tally[i][j] = 0.0;

  // loop over all particles, skip species not in mixture group

  for (i = 0; i < nlocal; i++) {
    ispecies = particles[i].ispecies;
    igroup = s2g[ispecies];
    if (igroup < 0) continue;
    icell = particles[i].icell;
    if (!(cinfo[icell].mask & groupbit)) continue;

    mass = species[ispecies].mass;
    v = particles[i].v;

    // 6 tallies per particle: N, Mass, mVx, mVy, mVz, mV^2

    vec = tally[icell];
    k = igroup*npergroup;

    vec[k++] += 1.0;
    vec[k++] += mass;
    vec[k++] += mass*v[0];
    vec[k++] += mass*v[1];
    vec[k++] += mass*v[2];
    vec[k++] += mass * (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  }
}

/* ----------------------------------------------------------------------
   query info about internal tally array for this compute
   index = which column of output (0 for vec, 1 to N for array)
   return # of tally quantities for this index
   also return array = ptr to tally array
   also return cols = ptr to list of columns in tally for this index
------------------------------------------------------------------------- */

int ComputeThermalGrid::query_tally_grid(int index, double **&array, int *&cols)
{
  index--;
  int ivalue = index % nvalue;
  array = tally;
  cols = map[index];
  return nmap[ivalue];
}

/* ----------------------------------------------------------------------
   tally accumulated info to compute final normalized values
   index = which column of output (0 for vec, 1 to N for array)
   for etally = NULL:
     use internal tallied info for single timestep, set nsample = 1
     compute values for all grid cells
       store results in vector_grid with nstride = 1 (single col of array_grid)
   for etally = ptr to caller array:
     use external tallied info for many timesteps
     nsample = additional normalization factor used by some values
     emap = list of etally columns to use, # of columns determined by index
     store results in caller's vec, spaced by nstride
   if norm = 0.0, set result to 0.0 directly so do not divide by 0.0
------------------------------------------------------------------------- */

void ComputeThermalGrid::
post_process_grid(int index, int nsample,
                  double **etally, int *emap, double *vec, int nstride)
{
  index--;
  int ivalue = index % nvalue;

  int lo = 0;
  int hi = nglocal;
  int k = 0;

  if (!etally) {
    nsample = 1;
    etally = tally;
    emap = map[index];
    vec = vector_grid;
    nstride = 1;
  }

  // compute normalized final value for each grid cell
  // Vcm = Sum mv / Sum m
  // total KE = 0.5 * Sum m(v - Vcm)^2
  // KE = 0.5 * Sum(i=xyz) [(Sum mVi^2) - M Vcmx^2]
  // KE = 0.5 * Sum(i=xyz) [(Sum mVi^2) - (Sum mVx)^2 / M]
  // KE = 0.5 * Sum (mv^2) - [(Sum mVx)^2 + (Sum mVy)^2 + (Sum mVz)^2] / M
  // thermal temp = (2/(3NkB)) * KE
  // press = (2/(3V)) * KE

  int tflag = 1;
  if (value[ivalue] == PRESS) tflag = 0;

  double prefactor;
  if (value[ivalue] == TEMP) prefactor = tprefactor;
  else if (value[ivalue] == PRESS) prefactor = pprefactor;

  Grid::ChildInfo *cinfo = grid->cinfo;
  double ncount,mass,mvx,mvy,mvz,mvsq;
  double *values;

  int n = emap[0];

  for (int icell = lo; icell < hi; icell++) {
    values = etally[icell];
    ncount = values[n];
    if (ncount <= 1.0) vec[k] = 0.0;
    else {
      mass = values[n+1];
      mvx = values[n+2];
      mvy = values[n+3];
      mvz = values[n+4];
      mvsq = values[n+5];
      vec[k] = mvsq - (mvx*mvx + mvy*mvy + mvz*mvz)/mass;
      vec[k] *= prefactor;
      if (tflag) vec[k] /= ncount;
      else vec[k] *= cinfo[icell].weight / cinfo[icell].volume / nsample;
    }
    k += nstride;
  }
}

/* ----------------------------------------------------------------------
   reallocate arrays if nglocal has changed
   called by init() and whenever grid changes
------------------------------------------------------------------------- */

void ComputeThermalGrid::reallocate()
{
  if (grid->nlocal == nglocal) return;

  memory->destroy(vector_grid);
  memory->destroy(tally);
  nglocal = grid->nlocal;
  memory->create(vector_grid,nglocal,"thermal/grid:vector_grid");
  memory->create(tally,nglocal,ntotal,"thermal/grid:tally");
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based data
------------------------------------------------------------------------- */

bigint ComputeThermalGrid::memory_usage()
{
  bigint bytes = 0;
  bytes = nglocal * sizeof(double);
  bytes = ntotal*nglocal * sizeof(double);
  return bytes;
}

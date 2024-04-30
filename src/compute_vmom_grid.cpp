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

#include "string.h"
#include "compute_vmom_grid.h"
#include "particle.h"
#include "mixture.h"
#include "grid.h"
#include "update.h"
#include "modify.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

// user keywords

enum{V2,VX3,VY3,VZ3,V4};

// internal accumulators

enum{MASSSUM,mV2,mV4,mV3x,mV3y,mV3z,LASTSIZE};

// max # of quantities to accumulate for any user value

#define MAXACCUMULATE 12

/* ---------------------------------------------------------------------- */

ComputeVmomGrid::ComputeVmomGrid(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal compute vmom/grid command");

  int igroup = grid->find_group(arg[2]);
  if (igroup < 0)
    error->all(FLERR,"Compute vmom/grid group ID does not exist");
  groupbit = grid->bitmask[igroup];

  imix = particle->find_mixture(arg[3]);
  if (imix < 0) error->all(FLERR,"Compute vmom/grid mixture ID "
                           "does not exist");
  ngroup = particle->mixture[imix]->ngroup;

  nvalue = narg - 4;
  value = new int[nvalue];

  npergroup = 0;
  unique = new int[LASTSIZE];
  nmap = new int[nvalue];
  memory->create(map,ngroup*nvalue,MAXACCUMULATE,"vmom/grid:map");
  for (int i = 0; i < nvalue; i++) nmap[i] = 0;

  int ivalue = 0;
  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"v2") == 0) {
      value[ivalue] = V2;
      set_map(ivalue,MASSSUM);
      set_map(ivalue,mV2);
    } else if (strcmp(arg[iarg],"v4") == 0) {
      value[ivalue] = V4;
      set_map(ivalue,MASSSUM);
      set_map(ivalue,mV4);
    } else if (strcmp(arg[iarg],"vx3") == 0) {
      value[ivalue] = VX3;
      set_map(ivalue,MASSSUM);
      set_map(ivalue,mV3x);
    } else if (strcmp(arg[iarg],"vy3") == 0) {
      value[ivalue] = VY3;
      set_map(ivalue,MASSSUM);
      set_map(ivalue,mV3y);
    } else if (strcmp(arg[iarg],"vz3") == 0) {
      value[ivalue] = VZ3;
      set_map(ivalue,MASSSUM);
      set_map(ivalue,mV3z);
    } else error->all(FLERR,"Illegal compute vmom/grid command");

    ivalue++;
    iarg++;
  }

  // ntotal = total # of columns in tally array
  // reset_map() adjusts indices in initial map() using final npergroup
  // also adds columns to tally array for CELLCOUNT/CELLMASS

  ntotal = ngroup*npergroup;
  reset_map();

  per_grid_flag = 1;
  size_per_grid_cols = ngroup*nvalue;
  post_process_grid_flag = 1;

  nglocal = 0;
  vector_grid = NULL;
  tally = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeVmomGrid::~ComputeVmomGrid()
{
  if (copymode) return;

  delete [] value;
  delete [] unique;

  delete [] nmap;
  memory->destroy(map);

  memory->destroy(vector_grid);
  memory->destroy(tally);
}

/* ---------------------------------------------------------------------- */

void ComputeVmomGrid::init()
{
  if (ngroup != particle->mixture[imix]->ngroup)
    error->all(FLERR,"Number of groups in compute vmom/grid mixture "
               "has changed");

  reallocate();
}

/* ---------------------------------------------------------------------- */

void ComputeVmomGrid::compute_per_grid()
{
  invoked_per_grid = update->ntimestep;

  Grid::ChildInfo *cinfo = grid->cinfo;
  Particle::Species *species = particle->species;
  Particle::OnePart *particles = particle->particles;
  int *s2g = particle->mixture[imix]->species2group;
  int nlocal = particle->nlocal;

  int i,j,k,m,ispecies,igroup,icell;
  double mass, vsq;
  double *v,*vec;

  double *sweights;
  int index_sweight = particle->find_custom((char *) "sweight");
  if(index_sweight >= 0)
    sweights = particle->edvec[particle->ewhich[index_sweight]];

  // zero all accumulators - could do this with memset()

  for (i = 0; i < nglocal; i++)
    for (j = 0; j < ntotal; j++)
      tally[i][j] = 0.0;

  // loop over all particles, skip species not in mixture group
  // perform all tallies needed for each particle
  // depends on its species group and the user-requested values

  for (i = 0; i < nlocal; i++) {
    ispecies = particles[i].ispecies;
    igroup = s2g[ispecies];
    if (igroup < 0) continue;
    icell = particles[i].icell;
    if (!(cinfo[icell].mask & groupbit)) continue;

    mass = species[ispecies].mass;
    v = particles[i].v;
    if(index_sweight >= 0) mass *= sweights[i]/update->fnum;

    vec = tally[icell];

    // loop has all possible values particle needs to accumulate
    // subset defined by user values are indexed by accumulate vector
    // NOTE: at some point may need prefactors v,v^2,v^3 converted to p,eng

    k = igroup*npergroup;
    vsq = (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

    for (m = 0; m < npergroup; m++) {
      switch (unique[m]) {
      case MASSSUM:
        vec[k++] += mass;
        break;
      case mV2:
        vec[k++] += mass*vsq;
        break;
      case mV3x:
        vec[k++] += mass*vsq*v[0];
        break;
      case mV3y:
        vec[k++] += mass*vsq*v[1];
        break;
      case mV3z:
        vec[k++] += mass*vsq*v[2];
        break;
      case mV4:
        vec[k++] += mass*vsq*vsq/3.0;
        break;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   query info about internal tally array for this compute
   index = which column of output (0 for vec, 1 to N for array)
   return # of tally quantities for this index
   also return array = ptr to tally array
   also return cols = ptr to list of columns in tally for this index
------------------------------------------------------------------------- */

int ComputeVmomGrid::query_tally_grid(int index, double **&array, int *&cols)
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

void ComputeVmomGrid::post_process_grid(int index, int nsample,
                                         double **etally, int *emap,
                                         double *vec, int nstride)
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


  switch (value[ivalue]) {

  case V2:
  case VX3:
  case VY3:
  case VZ3:
  case V4:
    {
      int summass = emap[0];
      int mvn = emap[1];
      for (int icell = lo; icell < hi; icell++) {
        vec[k] = etally[icell][mvn] / etally[icell][summass];
        k += nstride;
      }
      break;
    }

  }

}

/* ----------------------------------------------------------------------
   add a tally quantity to all groups for ivalue
   also add it to unique list if first time this name is used
   name = name of tally quantity from enum{} at top of file
   nmap[i] = # of tally quantities for user value I
   map[i][k] = index of Kth tally quantity for output value I
   npergroup = length of unique list
------------------------------------------------------------------------- */

void ComputeVmomGrid::set_map(int ivalue, int name)
{
  // index = loc of name in current unique list if there, else npergroup

  int index = 0;
  for (index = 0; index < npergroup; index++)
    if (unique[index] == name) break;

  // if name is not already in unique, add it and increment npergroup

  if (index == npergroup) {
    index = npergroup;
    unique[npergroup++] = name;
  }

  // add index to map and nmap for all groups
  // will add group offset in reset_map()

  for (int igroup = 0; igroup < ngroup; igroup++)
    map[igroup*nvalue+ivalue][nmap[ivalue]] = index;
  nmap[ivalue]++;
}

/* ----------------------------------------------------------------------
   reset map indices to reflect final npergroup = unique quantities/group
------------------------------------------------------------------------- */

void ComputeVmomGrid::reset_map()
{
  for (int i = 0; i < ngroup*nvalue; i++) {
    int igroup = i / nvalue;
    int ivalue = i % nvalue;
    for (int k = 0; k < nmap[ivalue]; k++)
      map[i][k] += igroup*npergroup;
  }
}

/* ----------------------------------------------------------------------
   reallocate data storage if nglocal has changed
   called by init() and whenever grid changes
------------------------------------------------------------------------- */

void ComputeVmomGrid::reallocate()
{
  if (grid->nlocal == nglocal) return;

  memory->destroy(vector_grid);
  memory->destroy(tally);
  nglocal = grid->nlocal;
  memory->create(vector_grid,nglocal,"grid:vector_grid");
  memory->create(tally,nglocal,ntotal,"grid:tally");
}

/* ----------------------------------------------------------------------
   memory usage of local grid-based data
------------------------------------------------------------------------- */

bigint ComputeVmomGrid::memory_usage()
{
  bigint bytes;
  bytes = nglocal * sizeof(double);
  bytes = ntotal*nglocal * sizeof(double);
  return bytes;
}

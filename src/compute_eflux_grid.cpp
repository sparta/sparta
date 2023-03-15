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
#include "compute_eflux_grid.h"
#include "particle.h"
#include "mixture.h"
#include "grid.h"
#include "update.h"
#include "modify.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

// user keywords

enum{HEATX,HEATY,HEATZ};

// internal accumulators

enum{MASSSUM,mVx,mVy,mVz,mVxVx,mVyVy,mVzVz,mVxVy,mVyVz,mVxVz,
     mVxVxVx,mVyVyVy,mVzVzVz,
     mVxVyVy,mVxVzVz,mVyVxVx,mVyVzVz,mVzVxVx,mVzVyVy,LASTSIZE};

// max # of quantities to accumulate for any user value

#define MAXACCUMULATE 12

/* ---------------------------------------------------------------------- */

ComputeEFluxGrid::ComputeEFluxGrid(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal compute eflux/grid command");

  int igroup = grid->find_group(arg[2]);
  if (igroup < 0)
    error->all(FLERR,"Compute eflux/grid group ID does not exist");
  groupbit = grid->bitmask[igroup];

  imix = particle->find_mixture(arg[3]);
  if (imix < 0) error->all(FLERR,"Compute eflux/grid mixture ID "
                           "does not exist");
  ngroup = particle->mixture[imix]->ngroup;

  nvalue = narg - 4;
  value = new int[nvalue];

  npergroup = 0;
  unique = new int[LASTSIZE];
  nmap = new int[nvalue];
  memory->create(map,ngroup*nvalue,MAXACCUMULATE,"eflux/grid:map");
  for (int i = 0; i < nvalue; i++) nmap[i] = 0;

  int ivalue = 0;
  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"heatx") == 0) {
      value[ivalue] = HEATX;
      set_map(ivalue,MASSSUM);
      set_map(ivalue,mVx);
      set_map(ivalue,mVy);
      set_map(ivalue,mVz);
      set_map(ivalue,mVxVx);
      set_map(ivalue,mVyVy);
      set_map(ivalue,mVzVz);
      set_map(ivalue,mVxVy);
      set_map(ivalue,mVxVz);
      set_map(ivalue,mVxVxVx);
      set_map(ivalue,mVxVyVy);
      set_map(ivalue,mVxVzVz);
    } else if (strcmp(arg[iarg],"heaty") == 0) {
      value[ivalue] = HEATX;
      set_map(ivalue,MASSSUM);
      set_map(ivalue,mVy);
      set_map(ivalue,mVx);
      set_map(ivalue,mVz);
      set_map(ivalue,mVyVy);
      set_map(ivalue,mVxVx);
      set_map(ivalue,mVzVz);
      set_map(ivalue,mVxVy);
      set_map(ivalue,mVyVz);
      set_map(ivalue,mVyVyVy);
      set_map(ivalue,mVyVxVx);
      set_map(ivalue,mVyVzVz);
    } else if (strcmp(arg[iarg],"heatz") == 0) {
      value[ivalue] = HEATX;
      set_map(ivalue,MASSSUM);
      set_map(ivalue,mVz);
      set_map(ivalue,mVx);
      set_map(ivalue,mVy);
      set_map(ivalue,mVzVz);
      set_map(ivalue,mVxVx);
      set_map(ivalue,mVyVy);
      set_map(ivalue,mVxVz);
      set_map(ivalue,mVyVz);
      set_map(ivalue,mVzVzVz);
      set_map(ivalue,mVzVxVx);
      set_map(ivalue,mVzVyVy);
    } else error->all(FLERR,"Illegal compute eflux/grid command");

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

ComputeEFluxGrid::~ComputeEFluxGrid()
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

void ComputeEFluxGrid::init()
{
  if (ngroup != particle->mixture[imix]->ngroup)
    error->all(FLERR,"Number of groups in compute eflux/grid mixture "
               "has changed");

  reallocate();
}

/* ---------------------------------------------------------------------- */

void ComputeEFluxGrid::compute_per_grid()
{
  invoked_per_grid = update->ntimestep;

  Grid::ChildInfo *cinfo = grid->cinfo;
  Particle::Species *species = particle->species;
  Particle::OnePart *particles = particle->particles;
  int *s2g = particle->mixture[imix]->species2group;
  int nlocal = particle->nlocal;

  int i,j,k,m,ispecies,igroup,icell;
  double mass;
  double *v,*vec;

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

    vec = tally[icell];

    // loop has all possible values particle needs to accumulate
    // subset defined by user values are indexed by accumulate vector
    // NOTE: at some point may need prefactors v,v^2,v^3 converted to p,eng

    k = igroup*npergroup;

    for (m = 0; m < npergroup; m++) {
      switch (unique[m]) {
      case MASSSUM:
        vec[k++] += mass;
        break;
      case mVx:
        vec[k++] += mass*v[0];
        break;
      case mVy:
        vec[k++] += mass*v[1];
        break;
      case mVz:
        vec[k++] += mass*v[2];
        break;
      case mVxVx:
        vec[k++] += mass*v[0]*v[0];
        break;
      case mVyVy:
        vec[k++] += mass*v[1]*v[1];
        break;
      case mVzVz:
        vec[k++] += mass*v[2]*v[2];
        break;
      case mVxVy:
        vec[k++] += mass*v[0]*v[1];
        break;
      case mVyVz:
        vec[k++] += mass*v[1]*v[2];
        break;
      case mVxVz:
        vec[k++] += mass*v[0]*v[2];
        break;
      case mVxVxVx:
        vec[k++] += mass*v[0]*v[0]*v[0];
        break;
      case mVyVyVy:
        vec[k++] += mass*v[1]*v[1]*v[1];
        break;
      case mVzVzVz:
        vec[k++] += mass*v[2]*v[2]*v[2];
        break;
      case mVxVyVy:
        vec[k++] += mass*v[0]*v[1]*v[1];
        break;
      case mVxVzVz:
        vec[k++] += mass*v[0]*v[2]*v[2];
        break;
      case mVyVxVx:
        vec[k++] += mass*v[1]*v[0]*v[0];
        break;
      case mVyVzVz:
        vec[k++] += mass*v[1]*v[2]*v[2];
        break;
      case mVzVxVx:
        vec[k++] += mass*v[2]*v[0]*v[0];
        break;
      case mVzVyVy:
        vec[k++] += mass*v[2]*v[1]*v[1];
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

int ComputeEFluxGrid::query_tally_grid(int index, double **&array, int *&cols)
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

void ComputeEFluxGrid::post_process_grid(int index, int nsample,
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

  // compute normalized final value for each grid cell
  // Vcm = Sum mv / Sum m = (Wx,Wy,Wz)
  // Wi = Sum mVi / M
  // heati = 0.5 * F/V Sum m (Vi - Wi) (V - W)^2
  // (Vi - Wi) (V - W)^2 = (Vi - Wi)
  //                       [ (Vi - Wi)^2 + (V1 - W1)^2 + (V2 - W2)^2 ]
  // i = xyz and 1,2 = xyx indices different than i
  // heati = 3 terms = h+h1+h2 where h2 is same as h1 with V1 replaced by V2
  // h = F/V Sum m (Vi-Wi)^3
  //   (Vi-Wi)^3 = Vi^3 - 3WiVi^2 + 3Wi^2Vi - Wi^3
  //   Sum m (Vi-Wi)^3 = Sum(mVi^3) - 3 Sum(mVi^2) Sum(mVi) / M +
  //                     2 Sum(mVi)^3 / M^2
  // h1 = F/V Sum m (Vi-Wi) (V1-W1)^2
  //   (Vi-Wi) (V1-W1)^2 = ViV1^2 - 2ViV1W1 + ViW1^2 - WiV1^2 + 2V1WiW1 - WiW1^2
  //   3 terms in previous equation combine to 1 term in next equation
  //   Sum m (Vi-Wi) (V1-W1)^2 = Sum(mViV1^2) - 2 Sum(mViV1) Sum(mV1) / M -
  //                             Sum(mVi) Sum(mV1^2) / M +
  //                             2 Sum(mVi) Sum(mV1)^2 / M^2

  double summass,h,h1,h2,wt;
  double *t;

  double fnum = update->fnum;
  Grid::ChildInfo *cinfo = grid->cinfo;

  int mass = emap[0];
  int mv = emap[1];
  int mv1 = emap[2];
  int mv2 = emap[3];
  int mvv = emap[4];
  int mv1v1 = emap[5];
  int mv2v2 = emap[6];
  int mvv1 = emap[7];
  int mvv2 = emap[8];
  int mvvv = emap[9];
  int mvv1v1 = emap[10];
  int mvv2v2 = emap[11];

  for (int icell = lo; icell < hi; icell++) {
    t = etally[icell];
    summass = t[mass];
    if (summass == 0.0) vec[k] = 0.0;
    else {
      h = t[mvvv] - 3.0*t[mv]*t[mvv]/summass +
        2.0*t[mv]*t[mv]*t[mv]/summass/summass;
      h1 = t[mvv1v1] - 2.0*t[mvv1]*t[mv1]/summass - t[mv]*t[mv1v1]/summass +
        2.0*t[mv]*t[mv1]*t[mv1]/summass/summass;
      h2 = t[mvv2v2] - 2.0*t[mvv2]*t[mv2]/summass - t[mv]*t[mv2v2]/summass +
        2.0*t[mv]*t[mv2]*t[mv2]/summass/summass;
      wt = 0.5 * fnum * cinfo[icell].weight / cinfo[icell].volume;
      vec[k] = wt/nsample * (h + h1 + h2);
    }
    k += nstride;
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

void ComputeEFluxGrid::set_map(int ivalue, int name)
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

void ComputeEFluxGrid::reset_map()
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

void ComputeEFluxGrid::reallocate()
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

bigint ComputeEFluxGrid::memory_usage()
{
  bigint bytes;
  bytes = nglocal * sizeof(double);
  bytes = ntotal*nglocal * sizeof(double);
  return bytes;
}

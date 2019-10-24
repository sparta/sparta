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
#include "compute_react_surf.h"
#include "update.h"
#include "domain.h"
#include "surf_react.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

#define DELTA 4096

/* ---------------------------------------------------------------------- */

ComputeReactSurf::ComputeReactSurf(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg != 4) error->all(FLERR,"Illegal compute react/surf command");

  int igroup = surf->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Compute react/surf group ID does not exist");
  groupbit = surf->bitmask[igroup];

  isr = surf->find_react(arg[3]);
  if (isr < 0) error->all(FLERR,"Compute react/surf reaction ID does not exist");

  ntotal = surf->sr[isr]->nlist;

  per_surf_flag = 1;
  size_per_surf_cols = ntotal;

  surf_tally_flag = 1;
  timeflag = 1;

  ntally = maxtally = 0;
  array_surf_tally = NULL;
  tally2surf = NULL;

  maxsurf = 0;
  array_surf = NULL;
  combined = 0;

  hash = new MyHash;

  dim = domain->dimension;
}

/* ---------------------------------------------------------------------- */

ComputeReactSurf::~ComputeReactSurf()
{
  memory->destroy(array_surf_tally);
  memory->destroy(tally2surf);
  delete hash;
}

/* ---------------------------------------------------------------------- */

void ComputeReactSurf::init()
{
  if (!surf->exist)
    error->all(FLERR,"Cannot use compute react/surf when surfs do not exist");
  if (surf->implicit) 
    error->all(FLERR,"Cannot use compute react/surf with implicit surfs");

  // NOTE: warn if some surfs are assigned to different surf react model





  // initialize tally array in case accessed before a tally timestep

  clear();

  combined = 0;
}

/* ----------------------------------------------------------------------
   no operations here, since compute results are stored in tally array
   just used by callers to indicate compute was used
   enables prediction of next step when update needs to tally
------------------------------------------------------------------------- */

void ComputeReactSurf::compute_per_surf()
{
  invoked_per_surf = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void ComputeReactSurf::clear()
{
  lines = surf->lines;
  tris = surf->tris;

  // clear hash of tallied surf IDs
  // called by Update at beginning of timesteps surf tallying is done

  hash->clear();
  ntally = 0;
  combined = 0;
}

/* ----------------------------------------------------------------------
   tally values for a single particle in icell
     colliding with surface element isurf, performing reaction (1 to N)
   iorig = particle ip before collision
   ip,jp = particles after collision
   ip = NULL means no particles after collision
   jp = NULL means one particle after collision
   jp != NULL means two particles after collision
------------------------------------------------------------------------- */

void ComputeReactSurf::surf_tally(int isurf, int icell, int reaction,
                                  Particle::OnePart *iorig, 
                                  Particle::OnePart *ip, Particle::OnePart *jp)
{
  // skip if no reaction

  if (reaction == 0) return;

  // skip if isurf not in surface group
  // or if this surf's reaction model is not a match

  if (dim == 2) {
    if (!(lines[isurf].mask & groupbit)) return;
    if (lines[isurf].isr != isr) return;
  } else {
    if (!(tris[isurf].mask & groupbit)) return;
    if (tris[isurf].isr != isr) return;
  }

  // itally = tally index of isurf
  // if 1st reaction on this isurf, add surf ID to hash
  // grow tally list if needed

  int itally;
  double *vec;

  surfint surfID;
  if (dim == 2) surfID = lines[isurf].id;
  else surfID = tris[isurf].id;

  if (hash->find(surfID) != hash->end()) itally = (*hash)[surfID];
  else {
    if (ntally == maxtally) grow_tally();
    itally = ntally;
    (*hash)[surfID] = itally;
    tally2surf[itally] = surfID;
    vec = array_surf_tally[itally];
    for (int i = 0; i < ntotal; i++) vec[i] = 0.0;
    ntally++;
  }

  // tally the reaction

  vec = array_surf_tally[itally];
  vec[reaction-1] += 1.0;
}

/* ----------------------------------------------------------------------
   return # of tallies and their indices into my local surf list
------------------------------------------------------------------------- */

int ComputeReactSurf::tallyinfo(surfint *&ptr)
{
  ptr = tally2surf;
  return ntally;
}

/* ----------------------------------------------------------------------
   sum tally values to owning surfs via surf->collate()
------------------------------------------------------------------------- */

void ComputeReactSurf::post_process_surf()
{
  if (combined) return;
  combined = 1;

  // reallocate array_surf if necessary

  int nown = surf->nown;

  if (nown > maxsurf) {
    memory->destroy(array_surf);
    maxsurf = nown;
    memory->create(array_surf,maxsurf,ntotal,"surf:array_surf");
  }

  // zero array_surf

  int i,j;
  for (i = 0; i < nown; i++)
    for (j = 0; j < ntotal; j++)
      array_surf[i][j] = 0.0;

  // collate entire array of results

  surf->collate_array(ntally,ntotal,tally2surf,array_surf_tally,array_surf);
}


/* ---------------------------------------------------------------------- */

void ComputeReactSurf::grow_tally()
{
  maxtally += DELTA;
  memory->grow(tally2surf,maxtally,"surf:tally2surf");
  memory->grow(array_surf_tally,maxtally,ntotal,"surf:array_surf_tally");
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

bigint ComputeReactSurf::memory_usage()
{
  bigint bytes = 0;
  bytes += ntotal*maxtally * sizeof(double);    // array_surf_tally
  bytes += maxtally * sizeof(surfint);          // tally2surf
  return bytes;
}

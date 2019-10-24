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
#include "compute_react_isurf_grid.h"
#include "particle.h"
#include "mixture.h"
#include "surf.h"
#include "surf_react.h"
#include "grid.h"
#include "update.h"
#include "domain.h"
#include "comm.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

#define DELTA 4096

/* ---------------------------------------------------------------------- */

ComputeReactISurfGrid::
ComputeReactISurfGrid(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg != 4) error->all(FLERR,"Illegal compute react/isurf/grid command");

  int igroup = grid->find_group(arg[2]);
  if (igroup < 0) 
    error->all(FLERR,"Compute react/isurf/grid group ID does not exist");
  groupbit = grid->bitmask[igroup];

  isr = surf->find_react(arg[3]);
  if (isr < 0) 
    error->all(FLERR,"Compute react/isurf/grid reaction ID does not exist");

  ntotal = surf->sr[isr]->nlist;

  per_grid_flag = 1;
  size_per_grid_cols = ntotal;
  post_process_isurf_grid_flag = 1;

  surf_tally_flag = 1;
  timeflag = 1;

  ntally = maxtally = 0;
  array_surf_tally = NULL;
  tally2surf = NULL;

  maxgrid = 0;
  array_grid = NULL;
  combined = 0;

  hash = new MyHash;

  dim = domain->dimension;
}

/* ---------------------------------------------------------------------- */

ComputeReactISurfGrid::~ComputeReactISurfGrid()
{
  memory->destroy(array_surf_tally);
  memory->destroy(tally2surf);
  memory->destroy(array_grid);
  delete hash;
}

/* ---------------------------------------------------------------------- */

void ComputeReactISurfGrid::init()
{
  if (!surf->exist)
    error->all(FLERR,"Cannot use compute react/isurf/grid when surfs do not exist");
  if (!surf->implicit) 
    error->all(FLERR,"Cannot use compute react/isurf/grid with explicit surfs");

  // NOTE: warn if some surfs are assigned to different surf react model




  // initialize tally array in case accessed before a tally timestep

  clear();

  combined = 0;
}

/* ----------------------------------------------------------------------
   no operations here, since compute results are stored in array_grid
   just used by callers to indicate compute was used
   enables prediction of next step when update needs to tally
   NOTE: also need to mark invoked_per_surf?
------------------------------------------------------------------------- */

void ComputeReactISurfGrid::compute_per_grid()
{
  invoked_per_grid = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void ComputeReactISurfGrid::clear()
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
   this method is exactly like ComputeSurf::surf_tally()
     except sum tally to to per-grid-cell array_grid
------------------------------------------------------------------------- */

void ComputeReactISurfGrid::surf_tally(int isurf, int icell, int reaction, 
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
  // if 1st particle hitting isurf, add surf ID to hash
  // grow tally list if needed
  // for implicit surfs, surfID is really a cellID

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

int ComputeReactISurfGrid::tallyinfo(surfint *&ptr)
{
  ptr = tally2surf;
  return ntally;
}

/* ----------------------------------------------------------------------
   sum surf tallies to owning cells via surf->collate()
   also copy split cell values to sub-cells for use by dump grid
------------------------------------------------------------------------- */

void ComputeReactISurfGrid::post_process_isurf_grid()
{
  if (combined) return;
  combined = 1;

  // reallocate array_grid if necessary

  int nglocal = grid->nlocal;

  if (nglocal > maxgrid) {
    memory->destroy(array_grid);
    maxgrid = nglocal;
    memory->create(array_grid,maxgrid,ntotal,"isurf/grid:array_grid");
  }

  // zero array_grid

  int i,j;
  for (i = 0; i < nglocal; i++)
    for (j = 0; j < ntotal; j++)
      array_grid[i][j] = 0.0;

  // perform rendezvous comm on tallies to sum them to my grid cells
  // array_surf_tally can be NULL if this proc has performed no tallies

  surf->collate_array_implicit(ntally,ntotal,tally2surf,
                               array_surf_tally,array_grid);
  
  // zero out result if icell not in grid group
  // can't apply until now, b/c tally included surfs in ghost cells and
  // cinfo does not have mask values for ghost cells

  Grid::ChildInfo *cinfo = grid->cinfo;

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) {
      for (j = 0; j < ntotal; j++)
        array_grid[icell][j] = 0.0;
    }
  }

  // copy split cell values to their sub cells, used by dump grid

  Grid::ChildCell *cells = grid->cells;
  Grid::SplitInfo *sinfo = grid->sinfo;

  int jcell,nsplit;
  int *csubs;

  for (int icell = 0; icell < nglocal; icell++) {
    if (cells[icell].nsplit <= 1) continue;
    nsplit = cells[icell].nsplit;
    csubs = sinfo[cells[icell].isplit].csubs;
    for (int j = 0; j < nsplit; j++) {
      jcell = csubs[j];
      memcpy(array_grid[jcell],array_grid[icell],ntotal*sizeof(double));
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeReactISurfGrid::grow_tally()
{
  maxtally += DELTA;
  memory->grow(tally2surf,maxtally,"isurf/grid:tally2surf");
  memory->grow(array_surf_tally,maxtally,ntotal,"isurf/grid:array_surf_tally");
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

bigint ComputeReactISurfGrid::memory_usage()
{
  bigint bytes = 0;
  bytes += ntotal*maxgrid * sizeof(double);     // array_grid
  bytes += ntotal*maxtally * sizeof(double);    // array_surf_tally
  bytes += maxtally * sizeof(surfint);          // tally2surf
  return bytes;
}

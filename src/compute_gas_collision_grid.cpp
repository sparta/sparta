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
#include "compute_gas_collision_grid.h"
#include "particle.h"
#include "mixture.h"
#include "grid.h"
#include "update.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

ComputeGasCollisionGrid::ComputeGasCollisionGrid(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg != 4) error->all(FLERR,"Illegal compute gas/collision/grid command");

  int igroup = grid->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Compute gas/collision/grid group ID does not exist");
  groupbit = grid->bitmask[igroup];

  imix = particle->find_mixture(arg[3]);
  if (imix < 0) error->all(FLERR,"Compute gas/collision/grid mixture ID does not exist");

  // setup

  per_grid_flag = 1;
  size_per_grid_cols = 0;

  gas_tally_flag = 1;         // triggers Collide to invoke gas_tally() for each collision
  timeflag = 1;               // tells Collide which timesteps to invoke gas_tally()

  nglocal = 0;
  vector_grid = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeGasCollisionGrid::~ComputeGasCollisionGrid()
{
  if (copy || copymode) return;

  memory->destroy(vector_grid);
}

/* ---------------------------------------------------------------------- */

void ComputeGasCollisionGrid::init()
{
  reallocate();
}

/* ----------------------------------------------------------------------
   no operations here, since compute results are stored in tally array
   just used by callers to indicate compute was used
   enables prediction of next step when update needs to tally
------------------------------------------------------------------------- */

void ComputeGasCollisionGrid::compute_per_grid()
{
  invoked_per_grid = update->ntimestep;
}

/* ----------------------------------------------------------------------
   called by Update before timesteps which invoke gas_tally()
---------------------------------------------------------------------- */

void ComputeGasCollisionGrid::clear()
{
  cinfo = grid->cinfo;
  memset(vector_grid,0,nglocal*sizeof(double));
}

/* ----------------------------------------------------------------------
   tally values for a single gas collision in icell
   reaction = 0 for collision only
   reaction = 1 to N for which reaction
   iorig,jorig = particles before collision
   ip,jp = particles after collision
   jp = NULL means one particle after collision
   kp != NULL means three particles after collision
------------------------------------------------------------------------- */

void ComputeGasCollisionGrid::gas_tally(int icell, int reaction,
                                        Particle::OnePart *iorig,
                                        Particle::OnePart *jorig,
                                        Particle::OnePart *ip,
                                        Particle::OnePart *jp,
                                        Particle::OnePart *kp)
{
  // skip if a reaction
  // this compute only tallies collisions that do not induce a reaction
  // reactions can be tallied by compute gas/reaction/grid command

  if (reaction) return;

  // skip if icell not in grid group

  if (!(cinfo[icell].mask & groupbit)) return;

  // skip if either particle species not in mixture group

  int igroup = particle->mixture[imix]->species2group[iorig->ispecies];
  int jgroup = particle->mixture[imix]->species2group[jorig->ispecies];
  if (igroup < 0 || jgroup < 0) return;

  // simply tally the collision to its grid cell

  vector_grid[icell] += 1.0;
}

/* ----------------------------------------------------------------------
   reallocate data storage if nglocal has changed
   called by init() and whenever grid changes
------------------------------------------------------------------------- */

void ComputeGasCollisionGrid::reallocate()
{
  if (grid->nlocal == nglocal) return;

  memory->destroy(vector_grid);
  nglocal = grid->nlocal;
  memory->create(vector_grid,nglocal,"gas/collision/grid:vector_grid");

  // clear counts b/c may be accessed before tallying is done
  //   e.g. on initial timestep of a new run, e.g. by dump grid
  //   this is different than compute_grid.cpp b/c compute_per_grid() is a no-op
  // also note if load-balancing is done, tallies will be lost
  //   would need to implement (un)pack_grid_one() to avoid this

  memset(vector_grid,0,nglocal*sizeof(double));
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

bigint ComputeGasCollisionGrid::memory_usage()
{
  bigint bytes = 0;
  bytes += nglocal * sizeof(double);    // vector_grid
  return bytes;
}

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
#include "compute_gas_reaction_grid.h"
#include "particle.h"
#include "mixture.h"
#include "grid.h"
#include "react.h"
#include "update.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{ALL,EVERY,SELECT};

/* ---------------------------------------------------------------------- */

ComputeGasReactionGrid::ComputeGasReactionGrid(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (!react || react->nlist == 0)
    error->all(FLERR,"Compute gas/reaction/grid cannot be used if gas-phase reactions not defined");

  if (narg < 5) error->all(FLERR,"Illegal compute gas/reaction/grid command");

  int igroup = grid->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Compute gas/reaction/grid group ID does not exist");
  groupbit = grid->bitmask[igroup];

  imix = particle->find_mixture(arg[3]);
  if (imix < 0) error->all(FLERR,"Compute gas/reaction/grid mixture ID does not exist");

  selectlist = NULL;
  reaction2col = NULL;
  
  if (strcmp(arg[4],"all") == 0) {
    if (narg != 5) error->all(FLERR,"Illegal compute gas/reaction/grid command");
    mode = ALL;
    ncol = 0;
    
  } else if (strcmp(arg[4],"every") == 0) {
    if (narg != 5) error->all(FLERR,"Illegal compute gas/reaction/grid command");
    mode = EVERY;
    ncol = react->nlist;
    
  } else {
    if (narg < 6) error->all(FLERR,"Illegal compute gas/reaction/grid command");
    mode = SELECT;
    ncol = narg-5;
    selectlist = new int[ncol];
    
    int index;
    int iarg = 5;
    for (int icol = 0; icol < ncol; icol++) {
      index = atoi(arg[iarg]);
      if (index <= 0 || index > react->nlist)
        error->all(FLERR,"Compute gas/reaction/grid reaction index is invalid");
      selectlist[icol] = index;
      iarg++;
    }
  }

  // convert selectlist to reaction2col
  // reaction2col[I] = column index (0 to Ncol-1) in array_grid for reaction I (1 to M)
  
  if (mode == SELECT) {
    reaction2col = new int[react->nlist + 1];
    for (int i = 0; i <= react->nlist; i++) reaction2col[i] = -1;
    for (int icol = 0; icol < ncol; icol++)
      reaction2col[selectlist[icol]] = icol;
  }
  
  // setup

  per_grid_flag = 1;
  size_per_grid_cols = 0;

  gas_tally_flag = 1;         // triggers Collide to invoke gas_tally() for each collision
  timeflag = 1;               // tells Collide which timesteps to invoke gas_tally()

  nglocal = 0;
  if (ncol == 0) vector_grid = NULL;
  else array_grid = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeGasReactionGrid::~ComputeGasReactionGrid()
{
  if (copy || copymode) return;

  delete [] selectlist;
  delete [] reaction2col;
  
  if (ncol == 0) memory->destroy(vector_grid);
  else memory->destroy(array_grid);
}

/* ---------------------------------------------------------------------- */

void ComputeGasReactionGrid::init()
{
  reallocate();
}

/* ----------------------------------------------------------------------
   no operations here, since compute results are stored in tally array
   just used by callers to indicate compute was used
   enables prediction of next step when update needs to tally
------------------------------------------------------------------------- */

void ComputeGasReactionGrid::compute_per_grid()
{
  invoked_per_grid = update->ntimestep;
}

/* ----------------------------------------------------------------------
   called by Update before timesteps which invoke gas_tally()
---------------------------------------------------------------------- */

void ComputeGasReactionGrid::clear()
{
  cinfo = grid->cinfo;
  if (ncol == 0) memset(vector_grid,0,nglocal*sizeof(double));
  else if (nglocal) memset(&array_grid[0][0],0,nglocal*ncol*sizeof(double));
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

void ComputeGasReactionGrid::gas_tally(int icell, int reaction,
                                       Particle::OnePart *iorig,
                                       Particle::OnePart *jorig,
                                       Particle::OnePart *ip,
                                       Particle::OnePart *jp,
                                       Particle::OnePart *kp)
{
  // skip if not a reaction
  // this compute only tallies collisions that induce a reaction
  // collisions can be tallied by compute gas/collision/grid command

  if (!reaction) return;
  
  // skip if icell not in grid group

  if (!(cinfo[icell].mask & groupbit)) return;

  // skip if either particle species not in mixture group

  int igroup = particle->mixture[imix]->species2group[iorig->ispecies];
  int jgroup = particle->mixture[imix]->species2group[jorig->ispecies];
  if (igroup < 0 || jgroup < 0) return;

  // simply tally the reaction to its grid cell
  // for EVERY and SELECT mode, reaction index determines column of array_grid
  
  if (mode == ALL) vector_grid[icell] += 1.0;
  else if (mode == EVERY) {
    int icol = reaction - 1;
    array_grid[icell][icol] += 1.0;
  } else if (mode == SELECT) {
    int icol = reaction2col[reaction];
    if (icol >= 0) array_grid[icell][icol] += 1.0;
  }
}

/* ----------------------------------------------------------------------
   reallocate data storage if nglocal has changed
   called by init() and whenever grid changes
------------------------------------------------------------------------- */

void ComputeGasReactionGrid::reallocate()
{
  if (grid->nlocal == nglocal) return;

  if (ncol == 0) memory->destroy(vector_grid);
  else memory->destroy(array_grid);
  nglocal = grid->nlocal;
  if (ncol == 0)
    memory->create(vector_grid,nglocal,"gas/collision/grid:vector_grid");
  else
    memory->create(array_grid,nglocal,ncol,"gas/collision/grid:vector_grid");

  // clear counts b/c may be accessed before tallying is done
  //   e.g. on initial timestep of a new run, e.g. by dump grid
  //   this is different than compute_grid.cpp b/c compute_per_grid() is a no-op
  // also note if load-balancing is done, tallies will be lost
  //   would need to implement (un)pack_grid_one() to avoid this

  if (ncol == 0) memset(vector_grid,0,nglocal*sizeof(double));
  else if (nglocal) memset(&array_grid[0][0],0,nglocal*ncol*sizeof(double));
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

bigint ComputeGasReactionGrid::memory_usage()
{
  bigint bytes = 0;
  if (ncol == 0) bytes += nglocal * sizeof(double);    // vector_grid
  else bytes += nglocal * ncol * sizeof(double);       // array_grid
  return bytes;
}

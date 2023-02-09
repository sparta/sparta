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

#include "stdlib.h"
#include "string.h"
#include "ctype.h"
#include "compute.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

#define DELTA 4

/* ---------------------------------------------------------------------- */

Compute::Compute(SPARTA *sparta, int narg, char **arg) : Pointers(sparta)
{
  if (narg < 2) error->all(FLERR,"Illegal compute command");

  // compute ID and style
  // ID must be all alphanumeric chars or underscores

  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  for (int i = 0; i < n-1; i++)
    if (!isalnum(id[i]) && id[i] != '_')
      error->all(FLERR,
                 "Compute ID must be alphanumeric or underscore characters");

  n = strlen(arg[1]) + 1;
  style = new char[n];
  strcpy(style,arg[1]);

  // set child class defaults

  scalar_flag = vector_flag = array_flag = 0;
  per_particle_flag = per_grid_flag = per_surf_flag = 0;
  post_process_grid_flag = post_process_isurf_grid_flag = 0;
  surf_tally_flag = boundary_tally_flag = 0;

  timeflag = 0;
  ntime = maxtime = 0;
  tlist = NULL;

  invoked_scalar = invoked_vector = invoked_array = -1;
  invoked_per_particle = invoked_per_grid = invoked_per_surf = -1;

  kokkos_flag = 0;
  copy = copymode = 0;
}

/* ---------------------------------------------------------------------- */

Compute::~Compute()
{
  if (copy || copymode) return;

  delete [] id;
  delete [] style;
  memory->destroy(tlist);
}

/* ----------------------------------------------------------------------
   add ntimestep to list of timesteps the compute will be called on
   do not add if already in list
   search from top downward, since list of times is in decreasing order
------------------------------------------------------------------------- */

void Compute::addstep(bigint ntimestep)
{
  // i = location in list to insert ntimestep

  int i;
  for (i = ntime-1; i >= 0; i--) {
    if (ntimestep == tlist[i]) return;
    if (ntimestep < tlist[i]) break;
  }
  i++;

  // extend list as needed

  if (ntime == maxtime) {
    maxtime += DELTA;
    memory->grow(tlist,maxtime,"compute:tlist");
  }

  // move remainder of list upward and insert ntimestep

  for (int j = ntime-1; j >= i; j--) tlist[j+1] = tlist[j];
  tlist[i] = ntimestep;
  ntime++;
}

/* ----------------------------------------------------------------------
   return 1/0 if ntimestep is or is not in list of calling timesteps
   if value(s) on top of list are less than ntimestep, delete them
   search from top downward, since list of times is in decreasing order
------------------------------------------------------------------------- */

int Compute::matchstep(bigint ntimestep)
{
  for (int i = ntime-1; i >= 0; i--) {
    if (ntimestep < tlist[i]) return 0;
    if (ntimestep == tlist[i]) return 1;
    if (ntimestep > tlist[i]) ntime--;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   clean out list of timesteps to call the compute on
------------------------------------------------------------------------- */

void Compute::clearstep()
{
  ntime = 0;
}

/* ---------------------------------------------------------------------- */

bigint Compute::memory_usage()
{
  bigint bytes = 0;
  return bytes;
}

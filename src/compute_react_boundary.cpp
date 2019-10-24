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
#include "compute_react_boundary.h"
#include "update.h"
#include "domain.h"
#include "comm.h"
#include "surf.h"
#include "surf_react.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

ComputeReactBoundary::
ComputeReactBoundary(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute react/boundary command");

  isr = surf->find_react(arg[2]);
  if (isr < 0) error->all(FLERR,"Compute react/boundary reaction ID "
                          "does not exist");

  ntotal = surf->sr[isr]->nlist;

  boundary_tally_flag = 1;
  timeflag = 1;
  array_flag = 1;
  nrow = 2 * domain->dimension;
  size_array_rows = nrow;
  size_array_cols = ntotal;

  memory->create(array,size_array_rows,size_array_cols,"react/boundary:array");
  memory->create(myarray,size_array_rows,size_array_cols,"react/boundary:array");
}

/* ---------------------------------------------------------------------- */

ComputeReactBoundary::~ComputeReactBoundary()
{
  memory->destroy(array);
  memory->destroy(myarray);
}

/* ---------------------------------------------------------------------- */

void ComputeReactBoundary::init()
{
  if (!domain->surfreactany && comm->me == 0)
    error->warning(FLERR,"Using compute react/boundary "
                   "when no box faces are assigned a reaction model");

  // initialize tally array in case accessed before a tally timestep

  clear();
}

/* ---------------------------------------------------------------------- */

void ComputeReactBoundary::compute_array()
{
  invoked_array = update->ntimestep;

  // sum tallies across processors

  MPI_Allreduce(&myarray[0][0],&array[0][0],nrow*ntotal,
                MPI_DOUBLE,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void ComputeReactBoundary::clear()
{
  surf_react = domain->surf_react;

  // reset tally values to zero
  // called by Update at beginning of timesteps boundary tallying is done

  for (int i = 0; i < size_array_rows; i++)
    for (int j = 0; j < ntotal; j++)
      myarray[i][j] = 0.0;
}

/* ----------------------------------------------------------------------
   tally values for a single particle colliding with boundary iface/istyle,
     performing reaction (1 to N)
   iorig = particle ip before collision
   ip,jp = particles after collision
   ip = NULL means no particles after collision
   jp = NULL means one particle after collision
   jp != NULL means two particles after collision
------------------------------------------------------------------------- */

void ComputeReactBoundary::boundary_tally(int iface, int istyle, int reaction,
                                     Particle::OnePart *iorig, 
                                     Particle::OnePart *ip, 
                                     Particle::OnePart *jp)
{
  // skip if no reaction

  if (reaction == 0) return;

  // skip if this face's reaction model is not a match

  if (surf_react[iface] != isr) return;

  // tally the reaction

  double *vec = myarray[iface];
  vec[reaction-1] += 1.0;
}

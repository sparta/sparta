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
#include "compute_react_boundary.h"
#include "update.h"
#include "domain.h"
#include "comm.h"
#include "surf.h"
#include "surf_react.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{REACTANT,PRODUCT};

/* ---------------------------------------------------------------------- */

ComputeReactBoundary::
ComputeReactBoundary(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal compute react/boundary command");

  isr = surf->find_react(arg[2]);
  if (isr < 0) error->all(FLERR,"Compute react/boundary reaction ID "
                          "does not exist");

  ntotal = surf->sr[isr]->nlist;

  rpflag = 0;
  reaction2col = NULL;

  // parse per-column reactant/product args
  // reset rpflag = 1 and ntotal = # of args

  if (narg > 3) {
    rpflag = 1;
    int ncol = narg - 3;
    memory->create(reaction2col,ntotal,ncol,"react/surf:reaction2col");
    for (int i = 0; i < ntotal; i++)
      for (int j = 0; j < ncol; j++)
        reaction2col[i][j] = 0;
    int which;
    int icol = 0;
    int iarg = 3;
    while (iarg < narg) {
      if (strncmp(arg[iarg],"r:",2) == 0) which = REACTANT;
      else if (strncmp(arg[iarg],"p:",2) == 0) which = PRODUCT;
      else error->all(FLERR,"Illegal compute react/surf command");
      int n = strlen(&arg[iarg][2]) + 1;
      char *copy = new char[n];
      strcpy(copy,&arg[iarg][2]);
      char *ptr = copy;
      while ((ptr = strtok(ptr,"/")) != (char *) NULL) {
        for (int ireaction = 0; ireaction < ntotal; ireaction++) {
          reaction2col[ireaction][icol] = 0;
          if (which == REACTANT) {
            if (surf->sr[isr]->match_reactant(ptr,ireaction))
              reaction2col[ireaction][icol] = 1;
          } else if (which == PRODUCT) {
            if (surf->sr[isr]->match_product(ptr,ireaction))
              reaction2col[ireaction][icol] = 1;
          }
        }
        ptr = NULL;
      }
      delete [] copy;
      icol++;
      iarg++;
    }
    ntotal = narg - 3;
  }

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
  memory->destroy(reaction2col);
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
  reaction--;

  // skip if this face's reaction model is not a match

  if (surf_react[iface] != isr) return;

  // tally the reaction
  // for rpflag, tally each column if r2c is 1 for this reaction
  // for rpflag = 0, tally the reaction directly

  double *vec = myarray[iface];

  if (rpflag) {
    int *r2c = reaction2col[reaction];
    for (int i = 0; i < ntotal; i++)
      if (r2c[i]) vec[i] += 1.0;
  } else vec[reaction] += 1.0;
}

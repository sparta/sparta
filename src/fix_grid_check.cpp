/* ----------------------------------------------------------------------
   DSMC - Sandia parallel DSMC code
   www.sandia.gov/~sjplimp/dsmc.html
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2011) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level DSMC directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "fix_grid_check.h"
#include "update.h"
#include "particle.h"
#include "grid.h"
#include "error.h"

using namespace DSMC_NS;

/* ---------------------------------------------------------------------- */

FixGridCheck::FixGridCheck(DSMC *dsmc, int narg, char **arg) : 
  Fix(dsmc, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal fix grid/check command");

  nevery = atoi(arg[2]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix grid/check command");

  scalar_flag = 1;
  global_freq = 1;
}

/* ---------------------------------------------------------------------- */

int FixGridCheck::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGridCheck::init()
{
  nflag = 0;
}

/* ---------------------------------------------------------------------- */

void FixGridCheck::end_of_step()
{
  if (update->ntimestep % nevery) return;

  Particle::OnePart *particles = particle->particles;
  Grid::OneCell *cells = grid->cells;
  int *mycells = grid->mycells;

  int nlocal = particle->nlocal;

  int icell;
  double *x,*lo,*hi;

  for (int i = 0; i < nlocal; i++) {
    x = particles[i].x;
    //icell = mycells[particles[i].icell];
    icell = particles[i].icell;
    lo = cells[icell].lo;
    hi = cells[icell].hi;
    if (x[0] < lo[0] || x[0] > hi[0] ||
	x[1] < lo[1] || x[1] > hi[1] ||
	x[2] < lo[2] || x[2] > hi[2]) nflag++;
  }
}

/* ----------------------------------------------------------------------
   return total count of out-of-cell particles across all procs
------------------------------------------------------------------------- */

double FixGridCheck::compute_scalar()
{
  double one = nflag;
  double all;
  MPI_Allreduce(&one,&all,1,MPI_DOUBLE,MPI_SUM,world);
  return all;
}

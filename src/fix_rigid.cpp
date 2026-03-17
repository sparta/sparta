/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "fix_rigid.h"
#include "update.h"
#include "domain.h"
#include "particle.h"
#include "grid.h"
#include "comm.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

FixRigid::FixRigid(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix rigid command");

  // optional args

  //outside_check = 0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"outside") == 0) {
    }
  }

  // setup

  //dim = domain->dimension;

  vector_flag = 1;
  global_freq = 1;
}

/* ---------------------------------------------------------------------- */

int FixRigid::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRigid::init()
{
  //ntotal = 0;
}

/* ---------------------------------------------------------------------- */

void FixRigid::setup()
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixRigid::start_of_step()
{
  // update RB props
  // change COM within compute surf
  // remap RB surfs to grid cells
}

/* ---------------------------------------------------------------------- */

void FixRigid::end_of_step()
{
  // query compute surf
  // sum forces/torques on RB
}

/* ----------------------------------------------------------------------
   return properties of the single rigid body
------------------------------------------------------------------------- */

double FixRigid::compute_vector(int index)
{
  return 0.0;
}

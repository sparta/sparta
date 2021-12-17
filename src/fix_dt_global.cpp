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

#include "stdlib.h"
#include "string.h"
#include "fix_dt_global.h"
#include "update.h"
#include "input.h"
#include "modify.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include <iostream>

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

FixDtGlobal::FixDtGlobal(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal fix dt global command");
  nevery = atoi(arg[2]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix dt global command");

  MPI_Comm_rank(world,&me);
}

/* ---------------------------------------------------------------------- */

FixDtGlobal::~FixDtGlobal()
{
  ;
}

/* ---------------------------------------------------------------------- */

void FixDtGlobal::init()
{
  std::cout << "greetings from FixDtGlobal::init()\n";
}

/* ---------------------------------------------------------------------- */

int FixDtGlobal::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDtGlobal::end_of_step()
{
  std::cout << "greetings from FixDtGlobal::end_of_step()\n";
}

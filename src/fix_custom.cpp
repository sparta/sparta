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

#include "fix_custom.h"
#include "custom.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

FixCustom::FixCustom(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal fix custom command");

  nevery = atoi(arg[2]);

  // instantiate Custom class for use by this fix
  // use Custom class to parse mode and list of actions
  // final arg = 1 for calling from FixCustom
  
  custom = new Custom(sparta);
  custom->parse_actions(narg-3,&arg[3],1);
}

/* ---------------------------------------------------------------------- */

FixCustom::~FixCustom()
{
  delete custom;
}

/* ---------------------------------------------------------------------- */

int FixCustom::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ----------------------------------------------------------------------
   invoked once every Nevery steps
------------------------------------------------------------------------- */

void FixCustom::end_of_step()
{
  // use Custom class to parse mode and list of actions
  // final arg = 1 for calling from FixCustom
  // NOTE: could make count a scalar output of fix
  
  bigint count = custom->process_actions(1);
}

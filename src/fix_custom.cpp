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
#include "modify.h"
#include "error.h"

using namespace SPARTA_NS;

enum{EQUAL,PARTICLE,GRID,SURF};

/* ---------------------------------------------------------------------- */

FixCustom::FixCustom(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal fix custom command");

  nevery = atoi(arg[2]);

  // instantiate Custom class for use by this fix
  // use Custom class to parse mode and list of SET and FILESTYLE actions
  // final arg = 1 for calling from FixCustom
  
  custom = new Custom(sparta);
  custom->process_actions(narg-3,&arg[3],1);
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
  // use Custom class to invoke list of SET and FILESTYLE actions
  
  bigint count = custom->process_actions();

  // if per-surf custom attributes were changed,
  //   need to invoke fixes which use them
  // example: fix emit/surf
  // fix surf/temp is NOT an example
  //   since it sets per-surf custom values, not uses them
  
  int mode = custom->mode;
  if (mode == SURF) modify->custom_surf_changed();

  // sum counts across all procs
  // NOTE: could make count a scalar output of fix
  
  bigint countall;
  MPI_Allreduce(&count,&countall,1,MPI_SPARTA_BIGINT,MPI_SUM,world);

  // DEBUG
  // printf("Fix custom attributes set = " BIGINT_FORMAT "\n",countall);
}

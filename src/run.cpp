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

#include "dsmctype.h"
#include "stdlib.h"
#include "string.h"
#include "run.h"
#include "domain.h"
#include "update.h"
#include "finish.h"
#include "timer.h"
#include "error.h"

using namespace DSMC_NS;

/* ---------------------------------------------------------------------- */

Run::Run(DSMC *dsmc) : Pointers(dsmc) {}

/* ---------------------------------------------------------------------- */

void Run::command(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal run command");

  if (domain->box_exist == 0)
    error->all(FLERR,"Run command before simulation box is defined");

  bigint nsteps_input = ATOBIGINT(arg[0]);
  if (nsteps_input < 0 || nsteps_input > MAXSMALLINT)
    error->all(FLERR,"Invalid run command N value");
  int nsteps = static_cast<int> (nsteps_input);

  // perform a single run

  update->nsteps = nsteps;
  update->firststep = update->ntimestep;
  update->laststep = update->ntimestep + nsteps;
  if (update->laststep < 0 || update->laststep > MAXBIGINT)
    error->all(FLERR,"Too many timesteps");
  
  dsmc->init();
  
  timer->barrier_start(TIME_LOOP);
  update->run(nsteps);
  timer->barrier_stop(TIME_LOOP);
  
  Finish finish(dsmc);
  finish.end();
  
  update->firststep = update->laststep = 0;
}

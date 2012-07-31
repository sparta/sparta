/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   www.sandia.gov/sparta.html
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2012) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "sptype.h"
#include "stdlib.h"
#include "string.h"
#include "run.h"
#include "domain.h"
#include "update.h"
#include "finish.h"
#include "timer.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

Run::Run(SPARTA *sparta) : Pointers(sparta) {}

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

  update->runflag = 1;
  update->nsteps = nsteps;
  update->firststep = update->ntimestep;
  update->laststep = update->ntimestep + nsteps;
  if (update->laststep < 0 || update->laststep > MAXBIGINT)
    error->all(FLERR,"Too many timesteps");
  
  sparta->init();
  update->setup();

  timer->barrier_start(TIME_LOOP);
  update->run(nsteps);
  timer->barrier_stop(TIME_LOOP);
  
  Finish finish(sparta);
  finish.end();
  
  update->runflag = 0;
  update->firststep = update->laststep = 0;
}

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
#include "fix_adapt.h"
#include "adapt_grid.h"
#include "grid.h"
#include "surf.h"
#include "comm.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "input.h"
#include "output.h"
#include "dump.h"
#include "error.h"

using namespace SPARTA_NS;

enum{NONE,REFINE,COARSEN};              // also in AdaptGrid

/* ---------------------------------------------------------------------- */

FixAdapt::FixAdapt(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal fix adapt command");

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;

  // create instance of AdaptGrid class

  me = comm->me;
  nprocs = comm->nprocs;

  adapt = new AdaptGrid(sparta);
  adapt->mode = 1;

  // parse and check arguments using AdaptGrid class

  nevery = input->inumeric(FLERR,arg[2]);
  if (nevery < 0) error->all(FLERR,"Illegal fix adapt command");

  adapt->process_args(narg-3,&arg[3]);
  adapt->check_args(nevery);

  action1 = adapt->action1;
  action2 = adapt->action2;
  file = adapt->file;

  coarsen_flag = 0;
  if (action1 == COARSEN || action2 == COARSEN)
    coarsen_flag = 1;

  if (file && strchr(file,'*') == NULL)
    error->all(FLERR,"Fix adapt filename must contain '*' character");

  // compute initial outputs

  last_adapt = 0;
  nrefine = ncoarsen = 0;
}

/* ---------------------------------------------------------------------- */

FixAdapt::~FixAdapt()
{
  delete adapt;
}

/* ---------------------------------------------------------------------- */

int FixAdapt::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAdapt::init()
{
  // re-check args in case computes or fixes changed

  adapt->check_args(nevery);

  // if any fix ave/grid exists, insure it comes before this fix
  // so that its output values are up-to-date on timesteps adaptation occurs

  int fixme = modify->find_fix(id);
  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style,"ave/grid") == 0) {
      if (i > fixme)
        error->all(FLERR,"Fix adapt must come after fix ave/grid");
    }
  }
}

/* ----------------------------------------------------------------------
   perform grid adaptation via AdaptGrid class
------------------------------------------------------------------------- */

void FixAdapt::end_of_step()
{
  // wrap adaptivity with clearstep/addstep since it may invoke computes

  modify->clearstep_compute();

  // same operations as in AdaptGrid single invocation

  grid->remove_ghosts();

  // memory allocation in AdaptGrid class

  adapt->setup(0);

  // perform adaptation

  if (action1 == REFINE || action2 == REFINE) grid->maxlevel++;

  if (action1 == REFINE) nrefine = adapt->refine();
  else if (action1 == COARSEN) ncoarsen = adapt->coarsen();

  if (action2 == REFINE) nrefine = adapt->refine();
  else if (action2 == COARSEN) ncoarsen = adapt->coarsen();

  grid->set_maxlevel();
  grid->rehash();

  // if no refinement or coarsening, just reghost/reneighbor and return

  if (nrefine == 0 && ncoarsen == 0) {
    last_adapt = 0;
    adapt->cleanup();
    grid->acquire_ghosts();
    grid->find_neighbors();
    return;
  }

  // memory deallocation in AdaptGrid class

  adapt->cleanup();

  // reset all attributes of adapted grid
  // same steps as in adapt_grid

  grid->setup_owned();
  grid->acquire_ghosts();
  grid->find_neighbors();
  grid->check_uniform();
  comm->reset_neighbors();

  if (surf->exist) {
    grid->set_inout();
    grid->type_check(0);
  }

  // notify all classes that store per-grid data that grid may have changed

  grid->notify_changed();

  // write out new grid file

  if (file) adapt->write_file();

  // wrap adaptivity with clearstep/addstep since it may invoke computes

  modify->addstep_compute(update->ntimestep + nevery);

  // outputs

  last_adapt = 1;
}

/* ----------------------------------------------------------------------
   return 0/1 for whether last adaptation changed grid
------------------------------------------------------------------------- */

double FixAdapt::compute_scalar()
{
  return (double) last_adapt;
}

/* ----------------------------------------------------------------------
   return stats for last adaptation
------------------------------------------------------------------------- */

double FixAdapt::compute_vector(int i)
{
  if (i == 0) return (double) nrefine;
  return (double) ncoarsen;
}

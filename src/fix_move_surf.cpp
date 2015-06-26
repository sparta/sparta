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
#include "fix_move_surf.h"
#include "move_surf.h"
#include "comm.h"
#include "input.h"
#include "grid.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

FixMoveSurf::FixMoveSurf(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal fix move/surf command");

  // create instance of MoveSurf class

  me = comm->me;
  nprocs = comm->nprocs;

  movesurf = new MoveSurf(sparta);
  //adapt->mode = 1;

  // parse and check arguments using AdaptGrid class

  nevery = input->inumeric(FLERR,arg[2]);
  if (nevery < 0) error->all(FLERR,"Illegal fix adapt command");

  //adapt->process_args(narg-3,&arg[3]);
  //adapt->check_args(nevery);

  // NOTE: check that file has * char in it
}

/* ---------------------------------------------------------------------- */

FixMoveSurf::~FixMoveSurf()
{
  delete movesurf;
}

/* ---------------------------------------------------------------------- */

int FixMoveSurf::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMoveSurf::init()
{
  // re-check args in case computes or fixes changed

  //adapt->check_args(nevery);
}

/* ----------------------------------------------------------------------
   perform grid adaptation via AdaptGrid class
------------------------------------------------------------------------- */

void FixMoveSurf::end_of_step()
{
  grid->remove_ghosts();

  // memory allocation in AdaptGrid class

  //adapt->setup();

  // perform surface move

  // memory deallocation in AdaptGrid class

  //adapt->cleanup();

  // reset all attributes of adapted grid
  // same steps as in adapt_grid

  /*
  grid->setup_owned();
  grid->acquire_ghosts();
  grid->find_neighbors();
  grid->check_uniform();
  comm->reset_neighbors();

  if (surf->exist) {
    grid->set_inout();
    grid->type_check(0);
  }
  */
}

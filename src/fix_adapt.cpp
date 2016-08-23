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
  // DEBUG
  //if (update->ntimestep > 140) return;

  // wrap adaptivity with clearstep/addstep since it may invoke computes

  modify->clearstep_compute();

  // same operations as in AdaptGrid single invocation

  grid->remove_ghosts();

  // memory allocation in AdaptGrid class

  adapt->setup(0);

  // perform adaptation

  int pstop = grid->nparent;

  if (action1 == REFINE) nrefine = adapt->refine();
  else if (action1 == COARSEN) ncoarsen = adapt->coarsen(pstop);

  if (action2 == REFINE) nrefine = adapt->refine();
  else if (action2 == COARSEN) ncoarsen = adapt->coarsen(pstop);

  // if no refine or coarsen, just reghost/reneighbor and return

  //if (comm->me == 0) printf("NREF COARSE %d %d nlocal %d\n",nrefine,ncoarsen,
  //                          grid->nlocal);

  if (nrefine == 0 && ncoarsen == 0) {
    last_adapt = 0;
    grid->acquire_ghosts();
    grid->find_neighbors();
    adapt->cleanup();
    return;
  }

  // memory deallocation in AdaptGrid class

  adapt->cleanup();

  // reset all attributes of adapted grid
  // same steps as in adapt_grid

  grid->setup_owned();
  grid->acquire_ghosts();
  grid->find_neighbors();

  /*
  int flag = 0;
  if (update->ntimestep == 50) flag = 1;

  if (flag) printf("AG NCELLS %d: %d %d\n",comm->me,grid->nlocal,grid->nghost);
  MPI_Barrier(world);
  if (flag) {
    //if (comm->me == 5) {
  for (int i = 0; i < grid->nlocal+grid->nghost; i++) {
    printf("  neighs %d %d %d %d: %d %d: %d %d: %d %d: %d %d\n",
           comm->me,i,grid->cells[i].id,grid->pcells[grid->cells[i].iparent].id,
           grid->neigh_decode(grid->cells[i].nmask,0),grid->cells[i].neigh[0],
           grid->neigh_decode(grid->cells[i].nmask,1),grid->cells[i].neigh[1],
           grid->neigh_decode(grid->cells[i].nmask,2),grid->cells[i].neigh[2],
           grid->neigh_decode(grid->cells[i].nmask,3),grid->cells[i].neigh[3]);
  }
  //}
  }
  */

  grid->check_uniform();
  comm->reset_neighbors();

  if (surf->exist) {
    grid->set_inout();
    grid->type_check(0);
  }

  // final update of any per grid fixes for all new child cells
  
  if (modify->n_pergrid) adapt->add_grid_fixes();

  // reallocate per grid cell arrays in per grid computes

  Compute **compute = modify->compute;
  for (int i = 0; i < modify->ncompute; i++)
    if (compute[i]->per_grid_flag) compute[i]->reallocate();

  // reallocate per grid arrays in per grid dumps

  for (int i = 0; i < output->ndump; i++)
    output->dump[i]->reset_grid();

  // write out new parent grid file

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

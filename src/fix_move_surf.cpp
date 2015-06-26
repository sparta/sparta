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
#include "domain.h"
#include "error.h"

using namespace SPARTA_NS;

enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};           // several files

/* ---------------------------------------------------------------------- */

FixMoveSurf::FixMoveSurf(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg < 6) error->all(FLERR,"Illegal fix move/surf command");

  if (!surf->exist)
    error->all(FLERR,"Cannot fix move/surf with no surf elements are defined");

  // create instance of MoveSurf class

  me = comm->me;
  nprocs = comm->nprocs;

  movesurf = new MoveSurf(sparta);
  movesurf->mode = 1;

  // parse and check arguments using MoveSurf class

  int igroup = surf->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Compute surf group ID does not exist");
  movesurf->groupbit = surf->bitmask[igroup];

  nevery = input->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix move/surf command");

  nlarge = input->inumeric(FLERR,arg[4]);
  if (nlarge < 0) error->all(FLERR,"Illegal fix move/surf command");
  if (nlarge % nevery) 
    error->all(FLERR,"Fix move/surf nlarge must be multiple of nevery");

  movesurf->process_args(narg-5,&arg[5]);
  action = movesurf->action;

  fraction = 1.0*nevery / nlarge;
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
}

/* ----------------------------------------------------------------------
   perform grid adaptation via AdaptGrid class
------------------------------------------------------------------------- */

void FixMoveSurf::end_of_step()
{
  int dim = domain->dimension;

  // sort particles

  if (particle->exist) particle->sort();

  // move points via chosen action by fractional amount

  movesurf->move_points(fraction);
  int *pflags = movesurf->pflags;

  // remake list of surf elements I own
  // assign split cell particles to parent split cell
  // assign surfs to grid cells

  surf->setup_surf();

  grid->unset_neighbors();
  grid->remove_ghosts();

  if (grid->nsplitlocal) {
    Grid::ChildCell *cells = grid->cells;
    int nglocal = grid->nlocal;
    for (int icell = 0; icell < nglocal; icell++)
      if (cells[icell].nsplit > 1)
	grid->combine_split_cell_particles(icell,1);
  }

  grid->clear_surf();
  grid->surf2grid(1,0);

  // re-setup owned and ghost cell info

  grid->setup_owned();
  grid->acquire_ghosts();
  grid->reset_neighbors();
  comm->reset_neighbors();

  // flag cells and corners as OUTSIDE or INSIDE

  grid->set_inout();
  grid->type_check(0);

  // remove particles in any cell that is now INSIDE or contains moved surfs
  // reassign particles in split cells to sub cell owner
  // compress particles if any flagged for deletion
  // NOTE: doc this logic better, here and in ReadSurf

  bigint ndeleted;
  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  int nglocal = grid->nlocal;
  int delflag = 0;

  for (int icell = 0; icell < nglocal; icell++) {
    if (cinfo[icell].type == INSIDE) {
      if (cinfo[icell].count) delflag = 1;
      particle->remove_all_from_cell(cinfo[icell].first);
      cinfo[icell].count = 0;
      cinfo[icell].first = -1;
      continue;
    }
    if (cells[icell].nsurf && cells[icell].nsplit >= 1) {
      int nsurf = cells[icell].nsurf;
      int *csurfs = cells[icell].csurfs;
      int m;
      if (dim == 2) {
        for (m = 0; m < nsurf; m++) {
          if (pflags[lines[csurfs[m]].p1]) break;
          if (pflags[lines[csurfs[m]].p2]) break;
        }
      } else {
        for (m = 0; m < nsurf; m++) {
          if (pflags[tris[csurfs[m]].p1]) break;
          if (pflags[tris[csurfs[m]].p2]) break;
          if (pflags[tris[csurfs[m]].p3]) break;
        }
      }
      if (m < nsurf) {
        if (cinfo[icell].count) delflag = 1;
        particle->remove_all_from_cell(cinfo[icell].first);
        cinfo[icell].count = 0;
        cinfo[icell].first = -1;
      }
    }
    if (cells[icell].nsplit > 1)
      grid->assign_split_cell_particles(icell);
  }
  int nlocal_old = particle->nlocal;
  if (delflag) particle->compress_rebalance();
  bigint delta = nlocal_old - particle->nlocal;
  MPI_Allreduce(&delta,&ndeleted,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
}

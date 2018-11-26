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
#include "remove_surf.h"
#include "surf.h"
#include "grid.h"
#include "comm.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

RemoveSurf::RemoveSurf(SPARTA *sparta) : Pointers(sparta) {}

/* ---------------------------------------------------------------------- */

RemoveSurf::~RemoveSurf() {}

/* ---------------------------------------------------------------------- */

void RemoveSurf::command(int narg, char **arg)
{
  if (!surf->exist) 
    error->all(FLERR,"Cannot remove_surf before surf elements are defined");

  if (narg < 1) error->all(FLERR,"Illegal remove_surf command");

  int igroup = surf->find_group(arg[0]);
  if (igroup < 0) error->all(FLERR,"Remove surf group ID does not exist");
  int groupbit = surf->bitmask[igroup];

  if (comm->me == 0)
    if (screen) fprintf(screen,"Removing surfs ...\n");

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  // sort particles

  if (particle->exist) particle->sort();

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  // remove all surfs in group
  // check that remaining surfs are still watertight
  // remake list of surf elements I own
  // assign split cell particles to parent split cell
  // assign surfs to grid cells

  if (domain->dimension == 2) {
    remove_2d(groupbit);
    surf->check_watertight_2d(0);
    if (surf->nline == 0) surf->exist = 0 ;
  } else {
    remove_3d(groupbit);
    surf->check_watertight_3d(0);
    if (surf->ntri == 0) surf->exist = 0 ;
  }

  surf->setup_surf();

  grid->unset_neighbors();
  grid->remove_ghosts();

  if (particle->exist && grid->nsplitlocal) {
    Grid::ChildCell *cells = grid->cells;
    int nglocal = grid->nlocal;
    for (int icell = 0; icell < nglocal; icell++)
      if (cells[icell].nsplit > 1)
	grid->combine_split_cell_particles(icell,1);
  }

  grid->clear_surf();
  if (surf->exist) grid->surf2grid(1);

  MPI_Barrier(world);
  double time3 = MPI_Wtime();

  // reassign particles in split cells to sub cell owner

  if (particle->exist && grid->nsplitlocal) {
    Grid::ChildCell *cells = grid->cells;
    int nglocal = grid->nlocal;
    for (int icell = 0; icell < nglocal; icell++)
      if (cells[icell].nsplit > 1)
	grid->assign_split_cell_particles(icell);
  }

  MPI_Barrier(world);
  double time4 = MPI_Wtime();

  // re-setup owned and ghost cell info

  grid->setup_owned();
  grid->acquire_ghosts();
  grid->reset_neighbors();
  comm->reset_neighbors();

  grid->set_inout();
  grid->type_check();

  // DEBUG
  //grid->debug();

  MPI_Barrier(world);
  double time5 = MPI_Wtime();

  double time_total = time3-time1;

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"  CPU time = %g secs\n",time_total);
      fprintf(screen,"  sort/surf2grid/particle/ghost percent = %g %g %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total,
              100.0*(time4-time3)/time_total,100.0*(time5-time4)/time_total);
    }
    if (logfile) {
      fprintf(logfile,"  CPU time = %g secs\n",time_total);
      fprintf(logfile,"  sort/surf2grid/particle/ghost percent = %g %g %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total,
              100.0*(time4-time3)/time_total,100.0*(time5-time4)/time_total);
    }
  }
}

/* ----------------------------------------------------------------------
   remove all lines in surf group
   condense data structures by removing deleted points & lines
------------------------------------------------------------------------- */

void RemoveSurf::remove_2d(int groupbit)
{
  int i;

  Surf::Line *lines = surf->lines;
  int nline_old = surf->nline;

  // remove lines not in group

  int nline = surf->nline;
  int nbytes = sizeof(Surf::Line);

  int n = 0;
  for (i = 0; i < nline; i++) {
    if (lines[i].mask & groupbit) continue;
    if (i != n) memcpy(&lines[n],&lines[i],nbytes);
    n++;
  }

  surf->nline = nline = n;

  // print stats after removal

  int nline_remove = nline_old - surf->nline;

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"  removed %d lines\n",nline_remove);
      fprintf(screen,"  %d lines remain\n",surf->nline);
    }
    if (logfile) {
      fprintf(logfile,"  removed %d lines\n",nline_remove);
      fprintf(logfile,"  %d lines remain\n",surf->nline);
    }
  }
}

/* ----------------------------------------------------------------------
   remove all triangels in surf group
   condense data structures by removing deleted points & triangles
------------------------------------------------------------------------- */

void RemoveSurf::remove_3d(int groupbit)
{
  int i;

  Surf::Tri *tris = surf->tris;
  int ntri_old = surf->ntri;

  // remove triangles not in group

  int ntri = surf->ntri;
  int nbytes = sizeof(Surf::Tri);

  int n = 0;
  for (i = 0; i < ntri; i++) {
    if (tris[i].mask & groupbit) continue;
    if (i != n) memcpy(&tris[n],&tris[i],nbytes);
    n++;
  }

  surf->ntri = ntri = n;

  // print stats after removal

  int ntri_remove = ntri_old - surf->ntri;

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"  removed %d tris\n",ntri_remove);
      fprintf(screen,"  %d tris remain\n",surf->ntri);
    }
    if (logfile) {
      fprintf(logfile,"  removed %d tris\n",ntri_remove);
      fprintf(logfile,"  %d tris remain\n",surf->ntri);
    }
  }
}

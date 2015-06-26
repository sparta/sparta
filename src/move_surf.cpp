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

#include "move_surf.h"
#include "surf.h"
#include "grid.h"
#include "comm.h"
#include "update.h"
#include "domain.h"
#include "random_mars.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};           // several files

/* ---------------------------------------------------------------------- */

int *pflag;

MoveSurf::MoveSurf(SPARTA *sparta) : Pointers(sparta)
{
  if (!surf->exist) 
    error->all(FLERR,"Cannot remove_surf before surf elements are defined");

  // pflag = flags for which points have moded

  memory->create(pflags,surf->npoint,"move_surf:pflags");

  // create RNG for style = RANDOM

  random = new RanPark(update->ranmaster->uniform());
}

/* ---------------------------------------------------------------------- */

MoveSurf::~MoveSurf()
{
  memory->destroy(pflags);
  delete random;
}

/* ---------------------------------------------------------------------- */

void MoveSurf::command(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal move_surf command");

  int dim = domain->dimension;

  if (comm->me == 0)
    if (screen) fprintf(screen,"Removing surfs ...\n");

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  // sort particles

  if (particle->exist) particle->sort();

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  // move surfs to new positions via surf points
  // remake list of surf elements I own
  // assign split cell particles to parent split cell
  // assign surfs to grid cells

  if (dim == 2) move_2d();
  else move_3d();

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
  grid->surf2grid(1);

  // NOTE: is this needed - move method to grid to avoid code duplication
  //if (dim == 2) check_point_near_surf_2d();
  //else check_point_near_surf_3d();

  MPI_Barrier(world);
  double time3 = MPI_Wtime();

  // re-setup owned and ghost cell info

  grid->setup_owned();
  grid->acquire_ghosts();
  grid->reset_neighbors();
  comm->reset_neighbors();

  MPI_Barrier(world);
  double time4 = MPI_Wtime();

  // flag cells and corners as OUTSIDE or INSIDE

  grid->set_inout();
  grid->type_check();

  // DEBUG
  //grid->debug();

  MPI_Barrier(world);
  double time5 = MPI_Wtime();

  // remove particles in any cell that is now INSIDE or contains moved surfs
  // reassign particles in split cells to sub cell owner
  // compress particles if any flagged for deletion
  // NOTE: doc this logic better, here and in ReadSurf

  bigint ndeleted;
  if (particle->exist) {
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

  MPI_Barrier(world);
  double time6 = MPI_Wtime();

  double time_total = time6-time1;

  if (comm->me == 0) {
    if (screen) {
      if (particle->exist)
	fprintf(screen,"  " BIGINT_FORMAT " deleted particles\n",ndeleted);
      fprintf(screen,"  CPU time = %g secs\n",time_total);
      fprintf(screen,"  sort/surf2grid/ghost/inout/particle percent = "
	      "%g %g %g %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total,
              100.0*(time4-time3)/time_total,100.0*(time5-time4)/time_total,
	      100.0*(time6-time5)/time_total);
    }
    if (logfile) {
      if (particle->exist)
	fprintf(logfile,"  " BIGINT_FORMAT " deleted particles\n",ndeleted);
      fprintf(logfile,"  CPU time = %g secs\n",time_total);
      fprintf(logfile,"  sort/surf2grid/ghost/inout/particle percent = "
	      "%g %g %g %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total,
              100.0*(time4-time3)/time_total,100.0*(time5-time4)/time_total,
	      100.0*(time6-time5)/time_total);
    }
  }
}

/* ----------------------------------------------------------------------
   move 2d surf points
------------------------------------------------------------------------- */

void MoveSurf::move_2d()
{
  int nline = surf->nline;

  Surf::Point *pts = surf->pts;
  int npoint = surf->npoint;

  double lo = 0.50;
  double hi = 1.5;
   
  for (int i = 0; i < npoint; i++) {
    double rad = lo + (hi-lo)*random->uniform();
    double delx = pts[i].x[0] - 5.0;
    double dely = pts[i].x[1] - 7.0;
    delx *= rad;
    dely *= rad;
    //delx *= 1.5;
    //dely *= 1.5;
    pts[i].x[0] = 5.0 + delx;
    pts[i].x[1] = 7.0 + dely;
    pflags[i] = 1;
  }

  surf->compute_line_normal(0,nline);
}

/* ----------------------------------------------------------------------
   move 3d surf points
------------------------------------------------------------------------- */

void MoveSurf::move_3d()
{
  //surf->compute_tri_normal(0,ntri);
}

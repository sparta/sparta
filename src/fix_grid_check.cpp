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

#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "fix_grid_check.h"
#include "update.h"
#include "domain.h"
#include "particle.h"
#include "grid.h"
#include "comm.h"
#include "cut2d.h"
#include "cut3d.h"
#include "error.h"

using namespace SPARTA_NS;

enum{ERROR,WARNING,SILENT};
enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};   // same as Grid

/* ---------------------------------------------------------------------- */

FixGridCheck::FixGridCheck(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix grid/check command");

  nevery = atoi(arg[2]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix grid/check command");

  if (strcmp(arg[3],"error") == 0) outflag = ERROR;
  else if (strcmp(arg[3],"warn") == 0) outflag = WARNING;
  else if (strcmp(arg[3],"silent") == 0) outflag = SILENT;
  else error->all(FLERR,"Illegal fix grid/check command");

  // optional args

  outside_check = 0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"outside") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix grid/check command");
      if (strcmp(arg[iarg+1],"no") == 0) outside_check = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) outside_check = 1;
      else error->all(FLERR,"Illegal fix grid/check command");
      iarg += 2;
    }
  }

  if (outside_check && !surf->implicit)
    error->all(FLERR,"Fix grid/check outside yes requires implicit surfs");

  // setup

  dim = domain->dimension;
  cut2d = NULL;
  cut3d = NULL;

  if (outside_check) {
    if (dim == 3) cut3d = new Cut3d(sparta);
    else cut2d = new Cut2d(sparta,domain->axisymmetric);
  }

  scalar_flag = 1;
  global_freq = 1;
}

/* ---------------------------------------------------------------------- */

FixGridCheck::~FixGridCheck()
{
  delete cut2d;
  delete cut3d;
}

/* ---------------------------------------------------------------------- */

int FixGridCheck::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGridCheck::init()
{
  ntotal = 0;
}

/* ---------------------------------------------------------------------- */

void FixGridCheck::setup()
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixGridCheck::end_of_step()
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  Grid::SplitInfo *sinfo = grid->sinfo;
  Particle::OnePart *particles = particle->particles;
  int nglocal = grid->nlocal;
  int nlocal = particle->nlocal;

  int icell;
  double *x,*lo,*hi;

  // check if icell is a valid cell for owning particles
  // check for split cell is whether particle is inside parent cell

  int nflag = 0;

  for (int i = 0; i < nlocal; i++) {
    icell = particles[i].icell;

    // is icell a valid index

    if (icell < 0 || icell >= nglocal) {
      if (outflag == ERROR) {
        char str[128];
        sprintf(str,
                "Particle %d,%d on proc %d is in invalid cell " CELLINT_FORMAT
                " on timestep " BIGINT_FORMAT,
                i,particles[i].id,comm->me,cells[icell].id,update->ntimestep);
        error->one(FLERR,str);
      }
      nflag++;
    }

    // does particle coord match icell bounds

    lo = cells[icell].lo;
    hi = cells[icell].hi;
    x = particles[i].x;
    if (x[0] < lo[0] || x[0] > hi[0] ||
        x[1] < lo[1] || x[1] > hi[1] ||
        x[2] < lo[2] || x[2] > hi[2]) {
      if (outflag == ERROR) {
        //printf("BAD %d %d " CELLINT_FORMAT ": "
        //    "x %g %g %g: lo %g %g %g hi: %g %g %g\n",
        //    i,icell,cells[icell].id,x[0],x[1],x[2],
        //    lo[0],lo[1],lo[2],hi[0],hi[1],hi[2]);
        char str[128];
        sprintf(str,
                "Particle %d,%d on proc %d is outside cell " CELLINT_FORMAT
                " on timestep " BIGINT_FORMAT,
                i,particles[i].id,comm->me,cells[icell].id,
                update->ntimestep);
        error->one(FLERR,str);
      }
      nflag++;
    }

    // error if icell is a split cell, since should be a sub cell

    if (cells[icell].nsplit > 1) {
      if (outflag == ERROR) {
        char str[128];
        sprintf(str,
                "Particle %d,%d on proc %d is in split cell " CELLINT_FORMAT
                " on timestep " BIGINT_FORMAT,
                i,particles[i].id,comm->me,cells[icell].id,
                update->ntimestep);
        error->one(FLERR,str);
      }
      nflag++;
    }

    // error if icell is an interior cell, since particle is inside surfs

    if (cinfo[icell].type == INSIDE) {
      if (outflag == ERROR) {
        char str[128];
        sprintf(str,
                "Particle %d,%d on proc %d is in interior cell " CELLINT_FORMAT
                " on timestep " BIGINT_FORMAT,
                i,particles[i].id,comm->me,cells[icell].id,update->ntimestep);
        error->one(FLERR,str);
      }
      nflag++;
    }

    // error if icell has zero volume, since collision attempt freq will blow up

    if (cinfo[icell].volume == 0.0) {
      if (outflag == ERROR) {
        char str[128];
        sprintf(str,
                "Particle %d,%d on proc %d is in volume=0 cell " CELLINT_FORMAT
                " on timestep " BIGINT_FORMAT,
                i,particles[i].id,comm->me,cells[icell].id,update->ntimestep);
        error->one(FLERR,str);
      }
      nflag++;
    }

    // check if particle in a cell with surfs is outside the surfs
    // for split cell, also verify particle is in correct sub cell
    // expensive, so only do this check if requested

    if (!outside_check) continue;
    if (cells[icell].nsurf == 0) continue;

    int splitcell,subcell,flag;

    if (cells[icell].nsplit <= 0) {
      splitcell = sinfo[cells[icell].isplit].icell;
      flag = grid->outside_surfs(splitcell,x,cut3d,cut2d);
    } else flag = grid->outside_surfs(icell,x,cut3d,cut2d);

    if (!flag) {
      if (outflag == ERROR) {
        char str[128];
        sprintf(str,
                "Particle %d,%d on proc %d is inside surfs in cell "
                CELLINT_FORMAT " on timestep " BIGINT_FORMAT,
                i,particles[i].id,comm->me,cells[icell].id,
                update->ntimestep);
        error->one(FLERR,str);
      }
      nflag++;
    }

    if (cells[icell].nsplit <= 0) {
      int subcell;
      if (dim == 2) subcell = update->split2d(splitcell,x);
      else subcell = update->split3d(splitcell,x);

      if (subcell != icell) {
        if (outflag == ERROR) {
          char str[128];
          sprintf(str,
                  "Particle %d,%d on proc %d is in wrong sub cell %d not %d"
                  " on timestep " BIGINT_FORMAT,
                  i,particles[i].id,comm->me,icell,subcell,
                  update->ntimestep);
          error->one(FLERR,str);
        }
        nflag++;
      }
    }
  }

  // -------------------------------------
  // done with all tests
  // warning message instead of error

  if (outflag == WARNING) {
    int all;
    MPI_Allreduce(&nflag,&all,1,MPI_INT,MPI_SUM,world);
    if (all && comm->me == 0) {
      char str[128];
      sprintf(str,"%d particles were in wrong cells on timestep "
              BIGINT_FORMAT,all,update->ntimestep);
      error->warning(FLERR,str);
    }
  }

  ntotal += nflag;
}

/* ----------------------------------------------------------------------
   return cummulative total count of out-of-cell particles across all procs
------------------------------------------------------------------------- */

double FixGridCheck::compute_scalar()
{
  double one = ntotal;
  double all;
  MPI_Allreduce(&one,&all,1,MPI_DOUBLE,MPI_SUM,world);
  return all;
}

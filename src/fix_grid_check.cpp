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
#include "particle.h"
#include "grid.h"
#include "comm.h"
#include "error.h"

using namespace SPARTA_NS;

enum{ERROR,WARNING,SILENT};

/* ---------------------------------------------------------------------- */

FixGridCheck::FixGridCheck(SPARTA *sparta, int narg, char **arg) : 
  Fix(sparta, narg, arg)
{
  if (narg != 4) error->all(FLERR,"Illegal fix grid/check command");

  nevery = atoi(arg[2]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix grid/check command");

  if (strcmp(arg[3],"error") == 0) outflag = ERROR;
  else if (strcmp(arg[3],"warn") == 0) outflag = WARNING;
  else if (strcmp(arg[3],"silent") == 0) outflag = SILENT;
  else error->all(FLERR,"Illegal fix grid/check command");

  scalar_flag = 1;
  global_freq = 1;
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
  nflag = 0;
}

/* ---------------------------------------------------------------------- */

void FixGridCheck::end_of_step()
{
  if (update->ntimestep % nevery) return;

  Grid::ChildCell *cells = grid->cells;
  Particle::OnePart *particles = particle->particles;
  int nglocal = grid->nlocal;
  int nlocal = particle->nlocal;

  int icell;
  double *x,*lo,*hi;

  // check if icell is a valid cell for owning particles
  // check for split cell is whether particle is inside parent cell

  for (int i = 0; i < nlocal; i++) {
    icell = particles[i].icell;

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

    lo = cells[icell].lo;
    hi = cells[icell].hi;
    x = particles[i].x;
    if (x[0] < lo[0] || x[0] > hi[0] ||
	x[1] < lo[1] || x[1] > hi[1] ||
	x[2] < lo[2] || x[2] > hi[2]) {
      if (outflag == ERROR) {
        //printf("BAD %d %d " CELLINT_FORMAT ": "
        //      "x %g %g %g: lo %g %g %g hi: %g %g %g\n",
        //      i,icell,cells[icell].id,x[0],x[1],x[2],
        //     lo[0],lo[1],lo[2],hi[0],hi[1],hi[2]);
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
  }

  // warning message

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
}

/* ----------------------------------------------------------------------
   return total count of out-of-cell particles across all procs
------------------------------------------------------------------------- */

double FixGridCheck::compute_scalar()
{
  double one = nflag;
  double all;
  MPI_Allreduce(&one,&all,1,MPI_DOUBLE,MPI_SUM,world);
  return all;
}

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

#include "stdlib.h"
#include "string.h"
#include "create_particles.h"
#include "particle.h"
#include "grid.h"
#include "comm.h"
#include "domain.h"
#include "random_park.h"
#include "error.h"

using namespace DSMC_NS;

enum{ALL,LOCAL};

/* ---------------------------------------------------------------------- */

CreateParticles::CreateParticles(DSMC *dsmc) : Pointers(dsmc) {}

/* ---------------------------------------------------------------------- */

void CreateParticles::command(int narg, char **arg)
{
  if (!domain->box_exist) 
    error->all(FLERR,
	       "Cannot create_particles before simulation box is defined");

  if (narg < 2) error->all(FLERR,"Illegal create_particles command");

  bigint n = ATOBIGINT(arg[0]);
  if (n < 0) error->all(FLERR,"Illegal create_particles command");
  seed = atoi(arg[1]);
  if (seed <= 0) error->all(FLERR,"Illegal create_particles command");

  // optional args

  int loop = ALL;

  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"loop") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_particles command");
      if (strcmp(arg[iarg+1],"all") == 0) loop = ALL;
      else if (strcmp(arg[iarg+1],"local") == 0) loop = LOCAL;
      else error->all(FLERR,"Illegal create_particles command");
      iarg += 2;
    } else error->all(FLERR,"Illegal create_particles command");
  }

  // generate particles

  bigint nprevious = particle->nglobal;

  if (loop == ALL) create_all(n);
  else if (loop == LOCAL) create_local(n);

  // error check

  bigint nglobal;
  bigint nme = particle->nlocal;
  MPI_Allreduce(&nme,&nglobal,1,MPI_DSMC_BIGINT,MPI_SUM,world);
  if (nglobal - nprevious != n) {
    char str[128];
    sprintf(str,"Created incorrect # of particles = " BIGINT_FORMAT,
	    nglobal-nprevious);
    error->all(FLERR,str);
  }
  particle->nglobal = nglobal;

  // print stats

  if (comm->me == 0) {
    if (screen) fprintf(screen,"Created " BIGINT_FORMAT " particles\n",n);
    if (logfile) fprintf(logfile,"Created " BIGINT_FORMAT " particles\n",n);
  }
}

/* ----------------------------------------------------------------------
   create N particles in serial
   every proc generates all N coords, only keeps those in cells it owns
------------------------------------------------------------------------- */

void CreateParticles::create_all(bigint n)
{
  int dimension = domain->dimension;
  double xlo = domain->boxlo[0];
  double ylo = domain->boxlo[1];
  double zlo = domain->boxlo[2];
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  int me = comm->me;
  RanPark *random = new RanPark(dsmc,seed);

  int icell;
  double x,y,z;

  // loop over all N particles

  for (bigint m = 0; m < n; m++) {
    x = xlo + random->uniform()*xprd;
    y = ylo + random->uniform()*yprd;
    z = zlo + random->uniform()*zprd;
    if (dimension == 2) z = 0.0;

    // which_cell() returns global grid cell index the particle is in
    // if I own that grid cell, add particle

    icell = grid->which_cell(x,y,z);
    if (grid->cells[icell].proc == me)
      particle->add_particle(0,1,icell,x,y,z);
  }

  delete random;
}

/* ----------------------------------------------------------------------
   create N particles in serial
   every proc generates all N coords, only keeps those in cells it owns
------------------------------------------------------------------------- */

void CreateParticles::create_local(bigint n)
{
  int dimension = domain->dimension;

  int me = comm->me;
  RanPark *random = new RanPark(dsmc,seed+me);

  Grid::OneCell *cells = grid->cells;
  int *mycells = grid->mycells;
  int nglocal = grid->nlocal;

  bigint nme = n/comm->nprocs;
  if (me < n % comm->nprocs) nme++;

  // loop only over Nme particles I own

  int ilocal,icell;
  double x,y,z;
  double *lo,*hi;

  for (bigint m = 0; m < nme; m++) {
    ilocal = static_cast<int> (random->uniform()*nglocal);
    icell = mycells[ilocal];
    lo = cells[icell].lo;
    hi = cells[icell].hi;

    x = lo[0] + random->uniform() * (hi[0]-lo[0]);
    y = lo[1] + random->uniform() * (hi[1]-lo[1]);
    z = lo[2] + random->uniform() * (hi[2]-lo[2]);
    if (dimension == 2) z = 0.0;

    particle->add_particle(0,1,icell,x,y,z);
  }

  delete random;
}


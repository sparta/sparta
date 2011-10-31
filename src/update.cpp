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

#include "update.h"
#include "particle.h"
#include "grid.h"
#include "comm.h"
#include "timer.h"
#include "error.h"

using namespace DSMC_NS;

/* ---------------------------------------------------------------------- */

Update::Update(DSMC *dsmc) : Pointers(dsmc)
{
  dt = 1.0;
}

/* ---------------------------------------------------------------------- */

Update::~Update() {}

/* ---------------------------------------------------------------------- */

void Update::setup()
{
  bigint pbytes,gbytes,bytes;
  pbytes = particle->memory_usage();
  gbytes = grid->memory_usage();
  bytes = pbytes + gbytes;

  double scale = 1.0/1024.0/1024.0;

  bigint ave,min,max;

  MPI_Allreduce(&pbytes,&ave,1,MPI_DSMC_BIGINT,MPI_SUM,world);
  double pave = scale * ave/comm->nprocs;
  MPI_Allreduce(&pbytes,&min,1,MPI_DSMC_BIGINT,MPI_MIN,world);
  double pmin = scale * min;
  MPI_Allreduce(&pbytes,&max,1,MPI_DSMC_BIGINT,MPI_MAX,world);
  double pmax = scale * max;

  MPI_Allreduce(&gbytes,&ave,1,MPI_DSMC_BIGINT,MPI_SUM,world);
  double gave = scale * ave/comm->nprocs;
  MPI_Allreduce(&gbytes,&min,1,MPI_DSMC_BIGINT,MPI_MIN,world);
  double gmin = scale * min;
  MPI_Allreduce(&gbytes,&max,1,MPI_DSMC_BIGINT,MPI_MAX,world);
  double gmax = scale * max;

  MPI_Allreduce(&bytes,&ave,1,MPI_DSMC_BIGINT,MPI_SUM,world);
  double tave = scale * ave/comm->nprocs;
  MPI_Allreduce(&bytes,&min,1,MPI_DSMC_BIGINT,MPI_MIN,world);
  double tmin = scale * min;
  MPI_Allreduce(&bytes,&max,1,MPI_DSMC_BIGINT,MPI_MAX,world);
  double tmax = scale * max;

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"Memory usage per proc in Mbytes:\n");
      fprintf(screen,"  particles (ave,min,max) = %g %g %g\n",
	      pave,pmin,pmax);
      fprintf(screen,"  grid      (ave,min,max) = %g %g %g\n",
	      gave,gmin,gmax);
      fprintf(screen,"  total     (ave,min,max) = %g %g %g\n",
	      tave,tmin,tmax);
    }
    if (logfile) {
      fprintf(logfile,"Memory usage per proc in Mbytes:\n");
      fprintf(logfile,"  particles (ave,min,max) = %g %g %g\n",
	      pave,pmin,pmax);
      fprintf(logfile,"  grid      (ave,min,max) = %g %g %g\n",
	      gave,gmin,gmax);
      fprintf(logfile,"  total     (ave,min,max) = %g %g %g\n",
	      tave,tmin,tmax);
    }
  }
}

/* ---------------------------------------------------------------------- */

void Update::run(int nsteps)
{
  if (comm->me == 0) {
    if (screen)
      fprintf(screen,"Performing run ...\n");
    if (logfile)
      fprintf(logfile,"Performing run ...\n");
  }

  // loop over timesteps

  for (int i = 0; i < nsteps; i++) {

    ntimestep++;

    // move particles

    timer->stamp();
    particle->move();
    timer->stamp(TIME_MOVE);

    // communicate particles

    timer->stamp();
    comm->migrate();
    timer->stamp(TIME_COMM);

    // sanity check on particles in correct cells

    //check();
  }
}

/* ---------------------------------------------------------------------- */

void Update::check()
{
  int icell;
  double *x,*lo,*hi;

  Particle::OnePart *particles = particle->particles;
  Grid::OneCell *cells = grid->cells;

  int nlocal = particle->nlocal;

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    x = particles[i].x;
    icell = particles[i].icell;
    lo = cells[icell].lo;
    hi = cells[icell].hi;
    if (x[0] < lo[0] || x[0] > hi[0] ||
	x[1] < lo[1] || x[1] > hi[1] ||
	x[2] < lo[2] || x[2] > hi[2]) flag++;
  }

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) {
    char str[128];
    sprintf(str,"%d particles are not in correct cell",flagall);
    error->all(FLERR,str);
  }
}


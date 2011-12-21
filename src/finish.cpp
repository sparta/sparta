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

#include "dsmctype.h"
#include "mpi.h"
#include "math.h"
#include "string.h"
#include "stdio.h"
#include "finish.h"
#include "update.h"
#include "particle.h"
#include "collide.h"
#include "comm.h"
#include "timer.h"
#include "memory.h"

using namespace DSMC_NS;

/* ---------------------------------------------------------------------- */

Finish::Finish(DSMC *dsmc) : Pointers(dsmc) {}

/* ---------------------------------------------------------------------- */

void Finish::end()
{
  int i,m;
  int histo[10];
  int loopflag,timeflag,histoflag;
  double time,tmp,ave,max,min;
  double time_loop,time_other;

  int me,nprocs;
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  loopflag = 1;
  timeflag = histoflag = 1;

  // loop stats

  if (loopflag) {
    time_other = timer->array[TIME_LOOP] -
      (timer->array[TIME_MOVE] + timer->array[TIME_COLLIDE] + 
       timer->array[TIME_SORT] + timer->array[TIME_COMM] +
       timer->array[TIME_OUTPUT]);
    
    time_loop = timer->array[TIME_LOOP];
    MPI_Allreduce(&time_loop,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time_loop = tmp/nprocs;
  }

  // recalculate nglobal

  bigint n = particle->nlocal;
  MPI_Allreduce(&n,&particle->nglobal,1,MPI_DSMC_BIGINT,MPI_SUM,world);

  // overall loop time

  if (me == 0) {
    if (screen) fprintf(screen,
			"Loop time of %g on %d procs for %d steps with " 
			BIGINT_FORMAT " atoms\n",
			time_loop,nprocs,update->nsteps,particle->nglobal);
    if (logfile) fprintf(logfile,
			 "Loop time of %g on %d procs for %d steps with " 
			 BIGINT_FORMAT " atoms\n",
			 time_loop,nprocs,update->nsteps,particle->nglobal);
  }

  // dummy stats for now

  bigint nlocal = particle->nlocal;
  bigint nptotal,nmtotal,ncctotal,ncmtotal;
  bigint nclatotal = 0;
  bigint ncltotal = 0;

  MPI_Allreduce(&nlocal,&nptotal,1,MPI_DSMC_BIGINT,MPI_SUM,world);
  MPI_Allreduce(&update->nmove,&nmtotal,1,MPI_DSMC_BIGINT,MPI_SUM,world);
  MPI_Allreduce(&update->ncellcross,&ncctotal,1,MPI_DSMC_BIGINT,
		MPI_SUM,world);
  MPI_Allreduce(&comm->ncomm,&ncmtotal,1,MPI_DSMC_BIGINT,MPI_SUM,world);
  if (collide) {
    MPI_Allreduce(&collide->ncollattempt,&nclatotal,1,MPI_DSMC_BIGINT,
		  MPI_SUM,world);
    MPI_Allreduce(&collide->ncollision,&ncltotal,1,MPI_DSMC_BIGINT,
		  MPI_SUM,world);
  }

  if (me == 0) {
    if (screen) {
      fprintf(screen,"\n");
      fprintf(screen,"Particle count = " BIGINT_FORMAT "\n",nptotal);
      fprintf(screen,"Particle moves = " BIGINT_FORMAT "\n",nmtotal);
      fprintf(screen,"Cells touched  = " BIGINT_FORMAT "\n",ncctotal);
      fprintf(screen,"Particle comms = " BIGINT_FORMAT "\n",ncmtotal);
      fprintf(screen,"Coll attempt   = " BIGINT_FORMAT "\n",nclatotal);
      fprintf(screen,"Coll performed = " BIGINT_FORMAT "\n",ncltotal);
      fprintf(screen,"Cell-touches/particle/step: %g\n",
	      1.0*ncctotal/particle->nglobal/update->nsteps);
      fprintf(screen,"Particle fraction migrating: %g\n",
	      1.0*ncmtotal/particle->nglobal/update->nsteps);
      fprintf(screen,"Collisions/particle/step: %g\n",
	      1.0*ncltotal/particle->nglobal/update->nsteps);
      fprintf(screen,"CPU/particle/step per proc: %g\n",
	      time_loop/particle->nglobal/update->nsteps * comm->nprocs);
      fprintf(screen,"CPU/particle/step in aggregate: %g\n",
	      time_loop/particle->nglobal/update->nsteps);
    }
    if (logfile) {
      fprintf(logfile,"\n");
      fprintf(logfile,"Particle count = " BIGINT_FORMAT "\n",nptotal);
      fprintf(logfile,"Particle moves = " BIGINT_FORMAT "\n",nmtotal);
      fprintf(logfile,"Cells touched  = " BIGINT_FORMAT "\n",ncctotal);
      fprintf(logfile,"Particle comms = " BIGINT_FORMAT "\n",ncmtotal);
      fprintf(logfile,"Coll attempt   = " BIGINT_FORMAT "\n",nclatotal);
      fprintf(logfile,"Coll performed = " BIGINT_FORMAT "\n",ncltotal);
      fprintf(logfile,"Cell-touches/particle/step: %g\n",
	      1.0*ncctotal/particle->nglobal/update->nsteps);
      fprintf(logfile,"Particle fraction migrating: %g\n",
	      1.0*ncmtotal/particle->nglobal/update->nsteps);
      fprintf(logfile,"Collisions/particle/step: %g\n",
	      1.0*ncltotal/particle->nglobal/update->nsteps);
      fprintf(logfile,"CPU/particle/step per proc: %g\n",
	      time_loop/particle->nglobal/update->nsteps * comm->nprocs);
      fprintf(logfile,"CPU/particle/step in aggregate: %g\n",
	      time_loop/particle->nglobal/update->nsteps);
    }
  }
  
  // timing breakdowns

  if (timeflag) {
    if (me == 0) {
      if (screen) fprintf(screen,"\n");
      if (logfile) fprintf(logfile,"\n");
    }

    time = timer->array[TIME_MOVE];
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) {
      if (screen) 
	fprintf(screen,"Move  time (%%) = %g (%g)\n",
		time,time/time_loop*100.0);
      if (logfile) 
	fprintf(logfile,"Move  time (%%) = %g (%g)\n",
		time,time/time_loop*100.0);
    }

    time = timer->array[TIME_COLLIDE];
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) {
      if (screen) 
	fprintf(screen,"Coll  time (%%) = %g (%g)\n",
		time,time/time_loop*100.0);
      if (logfile)
	fprintf(logfile,"Coll  time (%%) = %g (%g)\n",
		time,time/time_loop*100.0);
    }
    
    time = timer->array[TIME_SORT];
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) {
      if (screen) 
	fprintf(screen,"Sort  time (%%) = %g (%g)\n",
		time,time/time_loop*100.0);
      if (logfile)
	fprintf(logfile,"Sort  time (%%) = %g (%g)\n",
		time,time/time_loop*100.0);
    }
    
    time = timer->array[TIME_COMM];
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) {
      if (screen) 
	fprintf(screen,"Comm  time (%%) = %g (%g)\n",
		time,time/time_loop*100.0);
      if (logfile) 
	fprintf(logfile,"Comm  time (%%) = %g (%g)\n",
		time,time/time_loop*100.0);
    }
    
    time = timer->array[TIME_OUTPUT];
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) {
      if (screen) 
	fprintf(screen,"Outpt time (%%) = %g (%g)\n",
		time,time/time_loop*100.0);
      if (logfile) 
	fprintf(logfile,"Outpt time (%%) = %g (%g)\n",
		time,time/time_loop*100.0);
    }
    
    time = time_other;
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;
    if (me == 0) {
      if (screen) 
	fprintf(screen,"Other time (%%) = %g (%g)\n",
		time,time/time_loop*100.0);
      if (logfile) 
	fprintf(logfile,"Other time (%%) = %g (%g)\n",
		time,time/time_loop*100.0);
    }
  }
       
  // histograms

  if (histoflag) {
    if (me == 0) {
      if (screen) fprintf(screen,"\n");
      if (logfile) fprintf(logfile,"\n");
    }
    
    tmp = particle->nlocal;
    stats(1,&tmp,&ave,&max,&min,10,histo);
    if (me == 0) {
      if (screen) {
	fprintf(screen,"Nlocal:    %g ave %g max %g min\n",ave,max,min);
	fprintf(screen,"Histogram:");
	for (i = 0; i < 10; i++) fprintf(screen," %d",histo[i]);
	fprintf(screen,"\n");
      }
      if (logfile) {
	fprintf(logfile,"Nlocal:    %g ave %g max %g min\n",ave,max,min);
	fprintf(logfile,"Histogram:");
	for (i = 0; i < 10; i++) fprintf(logfile," %d",histo[i]);
	fprintf(logfile,"\n");
      }
    }
  }
    
  if (logfile) fflush(logfile);
}

/* ---------------------------------------------------------------------- */

void Finish::stats(int n, double *data, 
		   double *pave, double *pmax, double *pmin,
		   int nhisto, int *histo)
{
  int i,m;
  int *histotmp;

  double min = 1.0e20;
  double max = -1.0e20;
  double ave = 0.0;
  for (i = 0; i < n; i++) {
    ave += data[i];
    if (data[i] < min) min = data[i];
    if (data[i] > max) max = data[i];
  }

  int ntotal;
  MPI_Allreduce(&n,&ntotal,1,MPI_INT,MPI_SUM,world);
  double tmp;
  MPI_Allreduce(&ave,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  ave = tmp/ntotal;
  MPI_Allreduce(&min,&tmp,1,MPI_DOUBLE,MPI_MIN,world);
  min = tmp;
  MPI_Allreduce(&max,&tmp,1,MPI_DOUBLE,MPI_MAX,world);
  max = tmp;

  for (i = 0; i < nhisto; i++) histo[i] = 0;

  double del = max - min;
  for (i = 0; i < n; i++) {
    if (del == 0.0) m = 0;
    else m = static_cast<int> ((data[i]-min)/del * nhisto);
    if (m > nhisto-1) m = nhisto-1;
    histo[m]++;
  }

  memory->create(histotmp,nhisto,"finish:histotmp");
  MPI_Allreduce(histo,histotmp,nhisto,MPI_INT,MPI_SUM,world);
  for (i = 0; i < nhisto; i++) histo[i] = histotmp[i];
  memory->destroy(histotmp);

  *pave = ave;
  *pmax = max;
  *pmin = min;
}

/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2012) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "spatype.h"
#include "mpi.h"
#include "math.h"
#include "string.h"
#include "stdio.h"
#include "finish.h"
#include "update.h"
#include "particle.h"
#include "collide.h"
#include "comm.h"
#include "math_extra.h"
#include "timer.h"
#include "memory.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

Finish::Finish(SPARTA *sparta) : Pointers(sparta) {}

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
  MPI_Allreduce(&n,&particle->nglobal,1,MPI_SPARTA_BIGINT,MPI_SUM,world);

  // overall loop time

  if (me == 0) {
    if (screen) fprintf(screen,
			"Loop time of %g on %d procs for %d steps with " 
			BIGINT_FORMAT " molecules\n",
			time_loop,nprocs,update->nsteps,particle->nglobal);
    if (logfile) fprintf(logfile,
			 "Loop time of %g on %d procs for %d steps with " 
			 BIGINT_FORMAT " molecules\n",
			 time_loop,nprocs,update->nsteps,particle->nglobal);
  }

  // cummulative stats over entire run

  bigint nmove_total,ntouch_total,ncomm_total;
  bigint nboundary_total,nexit_total;
  bigint nscheck_total,nscollide_total;
  bigint nattempt_total = 0;
  bigint ncollide_total = 0;

  MPI_Allreduce(&update->nmove_running,&nmove_total,1,
		MPI_SPARTA_BIGINT,MPI_SUM,world);
  MPI_Allreduce(&update->ntouch_running,&ntouch_total,1,
		MPI_SPARTA_BIGINT,MPI_SUM,world);
  MPI_Allreduce(&update->ncomm_running,&ncomm_total,1,
		MPI_SPARTA_BIGINT,MPI_SUM,world);
  MPI_Allreduce(&update->nboundary_running,&nboundary_total,1,
		MPI_SPARTA_BIGINT,MPI_SUM,world);
  MPI_Allreduce(&update->nexit_running,&nexit_total,1,
		MPI_SPARTA_BIGINT,MPI_SUM,world);
  MPI_Allreduce(&update->nscheck_running,&nscheck_total,1,
		MPI_SPARTA_BIGINT,MPI_SUM,world);
  MPI_Allreduce(&update->nscollide_running,&nscollide_total,1,
		MPI_SPARTA_BIGINT,MPI_SUM,world);
  if (collide) {
    MPI_Allreduce(&collide->nattempt_running,&nattempt_total,1,
		  MPI_SPARTA_BIGINT,MPI_SUM,world);
    MPI_Allreduce(&collide->ncollide_running,&ncollide_total,1,
		  MPI_SPARTA_BIGINT,MPI_SUM,world);
  }

  double pms,ctps,pfc,pfcwb,pfeb,schps,sclps,caps,cps;
  pms = ctps = pfc = pfcwb = pfeb = schps = sclps = caps = cps = 0.0;
  if (update->nsteps) pms = 1.0*nmove_total/update->nsteps;
  if (nmove_total) {
    ctps = 1.0*ntouch_total/nmove_total;
    pfc = 1.0*ncomm_total/nmove_total;
    pfcwb = 1.0*nboundary_total/nmove_total;
    pfeb = 1.0*nexit_total/nmove_total;
    schps = 1.0*nscheck_total/nmove_total;
    sclps = 1.0*nscollide_total/nmove_total;
    caps = 1.0*nattempt_total/nmove_total;
    cps = 1.0*ncollide_total/nmove_total;
  }

  char str[32];

  if (me == 0) {
    if (screen) {
      fprintf(screen,"\n");
      fprintf(screen,"Particle moves = " BIGINT_FORMAT " %s\n",
	      nmove_total,MathExtra::num2str(nmove_total,str));
      fprintf(screen,"Cells touched  = " BIGINT_FORMAT " %s\n",
	      ntouch_total,MathExtra::num2str(ntouch_total,str));
      fprintf(screen,"Particle comms = " BIGINT_FORMAT " %s\n",
	      ncomm_total,MathExtra::num2str(ncomm_total,str));
      fprintf(screen,"Bound collides = " BIGINT_FORMAT " %s\n",
	      nboundary_total,MathExtra::num2str(nboundary_total,str));
      fprintf(screen,"Bound exits    = " BIGINT_FORMAT " %s\n",
	      nexit_total,MathExtra::num2str(nexit_total,str));
      fprintf(screen,"SurfColl check = " BIGINT_FORMAT " %s\n",
	      nscheck_total,MathExtra::num2str(nscheck_total,str));
      fprintf(screen,"SurfColl occur = " BIGINT_FORMAT " %s\n",
	      nscollide_total,MathExtra::num2str(nscollide_total,str));
      fprintf(screen,"Collide attmpt = " BIGINT_FORMAT " %s\n",
	      nattempt_total,MathExtra::num2str(nattempt_total,str));
      fprintf(screen,"Collide occurs = " BIGINT_FORMAT " %s\n",
	      ncollide_total,MathExtra::num2str(ncollide_total,str));
      fprintf(screen,"\n");
      fprintf(screen,"Particle-moves/step: %g\n",pms);
      fprintf(screen,"Cell-touches/particle/step: %g\n",ctps);
      fprintf(screen,"Particle fraction communicated: %g\n",pfc);
      fprintf(screen,"Particle fraction colliding with boundary: %g\n",pfcwb);
      fprintf(screen,"Particle fraction exiting boundary: %g\n",pfeb);
      fprintf(screen,"Surface-checks/particle/step: %g\n",schps);
      fprintf(screen,"Surface-collisions/particle/step: %g\n",sclps);
      fprintf(screen,"Collision-attempts/particle/step: %g\n",caps);
      fprintf(screen,"Collisions/particle/step: %g\n",cps);
    }
    if (logfile) {
      fprintf(logfile,"\n");
      fprintf(logfile,"Particle moves = " BIGINT_FORMAT " %s\n",
	      nmove_total,MathExtra::num2str(nmove_total,str));
      fprintf(logfile,"Cells touched  = " BIGINT_FORMAT " %s\n",
	      ntouch_total,MathExtra::num2str(ntouch_total,str));
      fprintf(logfile,"Particle comms = " BIGINT_FORMAT " %s\n",
	      ncomm_total,MathExtra::num2str(ncomm_total,str));
      fprintf(logfile,"Bound collides = " BIGINT_FORMAT " %s\n",
	      nboundary_total,MathExtra::num2str(nboundary_total,str));
      fprintf(logfile,"Bound exits    = " BIGINT_FORMAT " %s\n",
	      nexit_total,MathExtra::num2str(nexit_total,str));
      fprintf(logfile,"SurfColl check = " BIGINT_FORMAT " %s\n",
	      nscheck_total,MathExtra::num2str(nscheck_total,str));
      fprintf(logfile,"SurfColl occur = " BIGINT_FORMAT " %s\n",
	      nscollide_total,MathExtra::num2str(nscollide_total,str));
      fprintf(logfile,"Collide attmpt = " BIGINT_FORMAT " %s\n",
	      nattempt_total,MathExtra::num2str(nattempt_total,str));
      fprintf(logfile,"Collide occurs = " BIGINT_FORMAT " %s\n",
	      ncollide_total,MathExtra::num2str(ncollide_total,str));
      fprintf(logfile,"\n");
      fprintf(logfile,"Particle-moves/step: %g\n",pms);
      fprintf(logfile,"Cell-touches/particle/step: %g\n",ctps);
      fprintf(logfile,"Particle fraction communicated: %g\n",pfc);
      fprintf(logfile,"Particle fraction colliding with boundary: %g\n",pfcwb);
      fprintf(logfile,"Particle fraction exiting boundary: %g\n",pfeb);
      fprintf(logfile,"Surface-checks/particle/step: %g\n",schps);
      fprintf(logfile,"Surface-collisions/particle/step: %g\n",sclps);
      fprintf(logfile,"Collision-attempts/particle/step: %g\n",caps);
      fprintf(logfile,"Collisions/particle/step: %g\n",cps);
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

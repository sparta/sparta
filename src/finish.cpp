/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
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
#include "grid.h"
#include "collide.h"
#include "react.h"
#include "surf.h"
#include "surf_react.h"
#include "comm.h"
#include "math_extra.h"
#include "timer.h"
#include "memory.h"

using namespace SPARTA_NS;

// local function prototypes, code at end of file

static void mpi_timings(const char *label, Timer *t, int tt,
                        MPI_Comm world, const int nprocs,
                        const int me, double time_loop, FILE *scr, FILE *log);

/* ---------------------------------------------------------------------- */

Finish::Finish(SPARTA *sparta) : Pointers(sparta) {}

/* ---------------------------------------------------------------------- */

void Finish::end(int flag, double time_multiple_runs)
{
  int i;
  int histo[10];
  int loopflag,statsflag,timeflag,histoflag;
  double time,tmp,ave,max,min;
  double time_loop,time_other;

  int me,nprocs;
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  // choose flavors of statistical output
  // flag = 0 = just loop summary
  // flag = 1 = dynamics or minimization

  loopflag = 1;
  statsflag = timeflag = histoflag = 0;
  if (flag == 1) statsflag = timeflag = histoflag = 1;

  // loop stats
  // time_multiple_runs used for moves/CPU/proc statistic below

  if (loopflag) {
    time_other = timer->array[TIME_LOOP] -
      (timer->array[TIME_MOVE] + timer->array[TIME_COLLIDE] +
       timer->array[TIME_SORT] + timer->array[TIME_COMM] +
       timer->array[TIME_MODIFY] + timer->array[TIME_OUTPUT]);

    time_loop = timer->array[TIME_LOOP];
    MPI_Allreduce(&time_loop,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time_loop = tmp/nprocs;

    if (time_multiple_runs == 0.0) time_multiple_runs = time_loop;
    else {
      tmp = time_multiple_runs;
      MPI_Allreduce(&tmp,&time_multiple_runs,1,MPI_DOUBLE,MPI_SUM,world);
      time_multiple_runs /= nprocs;
    }
  }

  // recalculate nglobal

  bigint n = particle->nlocal;
  MPI_Allreduce(&n,&particle->nglobal,1,MPI_SPARTA_BIGINT,MPI_SUM,world);

  // overall loop time

  if (me == 0) {
    if (screen) fprintf(screen,
                        "Loop time of %g on %d procs for %d steps with "
                        BIGINT_FORMAT " particles\n",
                        time_loop,nprocs,update->nsteps,particle->nglobal);
    if (logfile) fprintf(logfile,
                         "Loop time of %g on %d procs for %d steps with "
                         BIGINT_FORMAT " particles\n",
                         time_loop,nprocs,update->nsteps,particle->nglobal);
  }

  // timing breakdowns

  if (timeflag) {
    const char hdr[] = "\nMPI task timing breakdown:\n"
      "Section |  min time  |  avg time  |  max time  |%varavg| %total\n"
      "---------------------------------------------------------------\n";
    if (me == 0) {
      if (screen)  fputs(hdr,screen);
      if (logfile) fputs(hdr,logfile);
    }

    mpi_timings("Move",timer,TIME_MOVE,world,nprocs,
                me,time_loop,screen,logfile);
    mpi_timings("Coll",timer,TIME_COLLIDE,world,nprocs,
                me,time_loop,screen,logfile);
    mpi_timings("Sort",timer,TIME_SORT,world,nprocs,
                me,time_loop,screen,logfile);
    mpi_timings("Comm",timer,TIME_COMM,world,nprocs,
                me,time_loop,screen,logfile);
    mpi_timings("Modify",timer,TIME_MODIFY,world,nprocs,
                me,time_loop,screen,logfile);
    mpi_timings("Output",timer,TIME_OUTPUT,world,nprocs,
                me,time_loop,screen,logfile);

    time = time_other;
    MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    time = tmp/nprocs;

    const char *fmt;
    fmt = "Other   |            |%- 12.4g|            |       |%6.2f\n";

    if (me == 0) {
      if (screen) fprintf(screen,fmt,time,time/time_loop*100.0);
      if (logfile) fprintf(logfile,fmt,time,time/time_loop*100.0);
    }
  }

  // cummulative stats over entire run

  if (statsflag) {
    bigint nmove_total,ntouch_total,ncomm_total;
    bigint nboundary_total,nexit_total;
    bigint nscheck_total,nscollide_total,nsreact_total;
    bigint nattempt_total = 0;
    bigint ncollide_total = 0;
    bigint nreact_total = 0;
    int stuck_total,axibad_total;

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
    MPI_Allreduce(&surf->nreact_running,&nsreact_total,1,
                  MPI_SPARTA_BIGINT,MPI_SUM,world);
    if (collide) {
      MPI_Allreduce(&collide->nattempt_running,&nattempt_total,1,
                    MPI_SPARTA_BIGINT,MPI_SUM,world);
      MPI_Allreduce(&collide->ncollide_running,&ncollide_total,1,
                    MPI_SPARTA_BIGINT,MPI_SUM,world);
      MPI_Allreduce(&collide->nreact_running,&nreact_total,1,
                    MPI_SPARTA_BIGINT,MPI_SUM,world);
    }
    MPI_Allreduce(&update->nstuck,&stuck_total,1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&update->naxibad,&axibad_total,1,MPI_INT,MPI_SUM,world);

    double pms,pmsp,ctps,cis,pfc,pfcwb,pfeb,schps,sclps,srps,caps,cps,rps;
    pms = pmsp = ctps = cis = pfc = pfcwb = pfeb =
      schps = sclps = srps = caps = cps = rps = 0.0;

    bigint elapsed = update->ntimestep - update->first_running_step;
    if (elapsed) pms = 1.0*nmove_total/elapsed;

    if (nmove_total) {
      pmsp = 1.0*nmove_total/time_multiple_runs/nprocs;
      ctps = 1.0*ntouch_total/nmove_total;
      cis = 1.0*update->niterate_running/elapsed;
      pfc = 1.0*ncomm_total/nmove_total;
      pfcwb = 1.0*nboundary_total/nmove_total;
      pfeb = 1.0*nexit_total/nmove_total;
      schps = 1.0*nscheck_total/nmove_total;
      sclps = 1.0*nscollide_total/nmove_total;
      srps = 1.0*nsreact_total/nmove_total;
      caps = 1.0*nattempt_total/nmove_total;
      cps = 1.0*ncollide_total/nmove_total;
      rps = 1.0*nreact_total/nmove_total;
    }

    char str[32];

    if (me == 0) {
      if (screen) {
        fprintf(screen,"\n");
        fprintf(screen,"Particle moves    = " BIGINT_FORMAT " %s\n",
                nmove_total,MathExtra::num2str(nmove_total,str));
        fprintf(screen,"Cells touched     = " BIGINT_FORMAT " %s\n",
                ntouch_total,MathExtra::num2str(ntouch_total,str));
        fprintf(screen,"Particle comms    = " BIGINT_FORMAT " %s\n",
                ncomm_total,MathExtra::num2str(ncomm_total,str));
        fprintf(screen,"Boundary collides = " BIGINT_FORMAT " %s\n",
                nboundary_total,MathExtra::num2str(nboundary_total,str));
        fprintf(screen,"Boundary exits    = " BIGINT_FORMAT " %s\n",
                nexit_total,MathExtra::num2str(nexit_total,str));
        fprintf(screen,"SurfColl checks   = " BIGINT_FORMAT " %s\n",
                nscheck_total,MathExtra::num2str(nscheck_total,str));
        fprintf(screen,"SurfColl occurs   = " BIGINT_FORMAT " %s\n",
                nscollide_total,MathExtra::num2str(nscollide_total,str));
        fprintf(screen,"Surf reactions    = " BIGINT_FORMAT " %s\n",
                nsreact_total,MathExtra::num2str(nsreact_total,str));
        fprintf(screen,"Collide attempts  = " BIGINT_FORMAT " %s\n",
                nattempt_total,MathExtra::num2str(nattempt_total,str));
        fprintf(screen,"Collide occurs    = " BIGINT_FORMAT " %s\n",
                ncollide_total,MathExtra::num2str(ncollide_total,str));
        fprintf(screen,"Gas reactions     = " BIGINT_FORMAT " %s\n",
                nreact_total,MathExtra::num2str(nreact_total,str));
        fprintf(screen,"Particles stuck   = %d\n",stuck_total);
        fprintf(screen,"Axisymm bad moves = %d\n",axibad_total);

        fprintf(screen,"\n");
        fprintf(screen,"Particle-moves/CPUsec/proc: %g\n",pmsp);
        fprintf(screen,"Particle-moves/step: %g\n",pms);
        fprintf(screen,"Cell-touches/particle/step: %g\n",ctps);
        fprintf(screen,"Particle comm iterations/step: %g\n",cis);
        fprintf(screen,"Particle fraction communicated: %g\n",pfc);
        fprintf(screen,"Particle fraction colliding with boundary: %g\n",pfcwb);
        fprintf(screen,"Particle fraction exiting boundary: %g\n",pfeb);
        fprintf(screen,"Surface-checks/particle/step: %g\n",schps);
        fprintf(screen,"Surface-collisions/particle/step: %g\n",sclps);
        fprintf(screen,"Surface-reactions/particle/step: %g\n",srps);
        fprintf(screen,"Collision-attempts/particle/step: %g\n",caps);
        fprintf(screen,"Collisions/particle/step: %g\n",cps);
        fprintf(screen,"Gas-reactions/particle/step: %g\n",rps);
      }
      if (logfile) {
        fprintf(logfile,"\n");
        fprintf(logfile,"Particle moves    = " BIGINT_FORMAT " %s\n",
                nmove_total,MathExtra::num2str(nmove_total,str));
        fprintf(logfile,"Cells touched     = " BIGINT_FORMAT " %s\n",
                ntouch_total,MathExtra::num2str(ntouch_total,str));
        fprintf(logfile,"Particle comms    = " BIGINT_FORMAT " %s\n",
                ncomm_total,MathExtra::num2str(ncomm_total,str));
        fprintf(logfile,"Boundary collides = " BIGINT_FORMAT " %s\n",
                nboundary_total,MathExtra::num2str(nboundary_total,str));
        fprintf(logfile,"Boundary exits    = " BIGINT_FORMAT " %s\n",
                nexit_total,MathExtra::num2str(nexit_total,str));
        fprintf(logfile,"SurfColl checks   = " BIGINT_FORMAT " %s\n",
                nscheck_total,MathExtra::num2str(nscheck_total,str));
        fprintf(logfile,"SurfColl occurs   = " BIGINT_FORMAT " %s\n",
                nscollide_total,MathExtra::num2str(nscollide_total,str));
        fprintf(logfile,"Surf reactions    = " BIGINT_FORMAT " %s\n",
                nsreact_total,MathExtra::num2str(nsreact_total,str));
        fprintf(logfile,"Collide attempts  = " BIGINT_FORMAT " %s\n",
                nattempt_total,MathExtra::num2str(nattempt_total,str));
        fprintf(logfile,"Collide occurs    = " BIGINT_FORMAT " %s\n",
                ncollide_total,MathExtra::num2str(ncollide_total,str));
        fprintf(logfile,"Reactions         = " BIGINT_FORMAT " %s\n",
                nreact_total,MathExtra::num2str(nreact_total,str));
        fprintf(logfile,"Particles stuck   = %d\n",stuck_total);
        fprintf(logfile,"Axisymm bad moves = %d\n",axibad_total);

        fprintf(logfile,"\n");
        fprintf(logfile,"Particle-moves/CPUsec/proc: %g\n",pmsp);
        fprintf(logfile,"Particle-moves/step: %g\n",pms);
        fprintf(logfile,"Cell-touches/particle/step: %g\n",ctps);
        fprintf(logfile,"Particle comm iterations/step: %g\n",cis);
        fprintf(logfile,"Particle fraction communicated: %g\n",pfc);
        fprintf(logfile,"Particle fraction colliding with boundary: %g\n",
                pfcwb);
        fprintf(logfile,"Particle fraction exiting boundary: %g\n",pfeb);
        fprintf(logfile,"Surface-checks/particle/step: %g\n",schps);
        fprintf(logfile,"Surface-collisions/particle/step: %g\n",sclps);
        fprintf(logfile,"Surf-reactions/particle/step: %g\n",srps);
        fprintf(logfile,"Collision-attempts/particle/step: %g\n",caps);
        fprintf(logfile,"Collisions/particle/step: %g\n",cps);
        fprintf(logfile,"Reactions/particle/step: %g\n",rps);
      }
    }
  }

  // gas per-reaction stats

  if (react) {
    if (me == 0) {
      if (screen) fprintf(screen,"\nGas reaction tallies:\n");
      if (logfile) fprintf(logfile,"\nGas reaction tallies:\n");
    }

    double tally;
    char *rID;
    int nlist = react->nlist;
    if (me == 0) {
      if (screen) fprintf(screen,"  style %s #-of-reactions %d\n",
                          react->style,nlist);
      if (logfile) fprintf(logfile,"  style %s #-of-reactions %d\n",
                           react->style,nlist);
    }
    for (int m = 0; m < nlist; m++) {
      tally = react->extract_tally(m);
      if (tally == 0.0) continue;
      rID = react->reactionID(m);
      if (me == 0) {
        if (screen) fprintf(screen,"  reaction %s: %g\n",rID,tally);
        if (logfile) fprintf(logfile,"  reaction %s: %g\n",rID,tally);
      }
    }
  }

  // surface per-reaction stats

  if (surf->nsr) {
    if (me == 0) {
      if (screen) fprintf(screen,"\nSurface reaction tallies:\n");
      if (logfile) fprintf(logfile,"\nSurface reaction tallies:\n");
    }

    double tally;
    char *rID;
    for (int i = 0; i < surf->nsr; i++) {
      SurfReact *sr = surf->sr[i];
      int nlist = sr->nlist;
      if (me == 0) {
        if (screen) fprintf(screen,"  id %s style %s #-of-reactions %d\n",
                            sr->id,sr->style,nlist);
        if (logfile) fprintf(logfile,"  id %s style %s #-of-reactions %d\n",
                             sr->id,sr->style,nlist);
      }
      tally = sr->compute_vector(1);
      if (me == 0) {
        if (screen) fprintf(screen,"    reaction all: %g\n",tally);
        if (logfile) fprintf(logfile,"    reaction all: %g\n",tally);
      }
      for (int m = 0; m < nlist; m++) {
        tally = sr->compute_vector(2+nlist+m);
        if (tally == 0.0) continue;
        rID = sr->reactionID(m);
        if (me == 0) {
          if (screen) fprintf(screen,"    reaction %s: %g\n",rID,tally);
          if (logfile) fprintf(logfile,"    reaction %s: %g\n",rID,tally);
        }
      }
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
        fprintf(screen,"Particles: %g ave %g max %g min\n",ave,max,min);
        fprintf(screen,"Histogram:");
        for (i = 0; i < 10; i++) fprintf(screen," %d",histo[i]);
        fprintf(screen,"\n");
      }
      if (logfile) {
        fprintf(logfile,"Particles: %g ave %g max %g min\n",ave,max,min);
        fprintf(logfile,"Histogram:");
        for (i = 0; i < 10; i++) fprintf(logfile," %d",histo[i]);
        fprintf(logfile,"\n");
      }
    }

    tmp = grid->nlocal;
    stats(1,&tmp,&ave,&max,&min,10,histo);
    if (me == 0) {
      if (screen) {
        fprintf(screen,"Cells:     %g ave %g max %g min\n",ave,max,min);
        fprintf(screen,"Histogram:");
        for (i = 0; i < 10; i++) fprintf(screen," %d",histo[i]);
        fprintf(screen,"\n");
      }
      if (logfile) {
        fprintf(logfile,"Cells:      %g ave %g max %g min\n",ave,max,min);
        fprintf(logfile,"Histogram:");
        for (i = 0; i < 10; i++) fprintf(logfile," %d",histo[i]);
        fprintf(logfile,"\n");
      }
    }

    tmp = grid->nghost;
    stats(1,&tmp,&ave,&max,&min,10,histo);
    if (me == 0) {
      if (screen) {
        fprintf(screen,"GhostCell: %g ave %g max %g min\n",ave,max,min);
        fprintf(screen,"Histogram:");
        for (i = 0; i < 10; i++) fprintf(screen," %d",histo[i]);
        fprintf(screen,"\n");
      }
      if (logfile) {
        fprintf(logfile,"GhostCell: %g ave %g max %g min\n",ave,max,min);
        fprintf(logfile,"Histogram:");
        for (i = 0; i < 10; i++) fprintf(logfile," %d",histo[i]);
        fprintf(logfile,"\n");
      }
    }

    tmp = grid->nempty;
    stats(1,&tmp,&ave,&max,&min,10,histo);
    if (me == 0) {
      if (screen) {
        fprintf(screen,"EmptyCell: %g ave %g max %g min\n",ave,max,min);
        fprintf(screen,"Histogram:");
        for (i = 0; i < 10; i++) fprintf(screen," %d",histo[i]);
        fprintf(screen,"\n");
      }
      if (logfile) {
        fprintf(logfile,"EmptyCell: %g ave %g max %g min\n",ave,max,min);
        fprintf(logfile,"Histogram:");
        for (i = 0; i < 10; i++) fprintf(logfile," %d",histo[i]);
        fprintf(logfile,"\n");
      }
    }

    if (surf->exist) {
      tmp = surf->nlocal;
      stats(1,&tmp,&ave,&max,&min,10,histo);
      if (me == 0) {
        if (screen) {
          fprintf(screen,"Surfs:     %g ave %g max %g min\n",ave,max,min);
          fprintf(screen,"Histogram:");
          for (i = 0; i < 10; i++) fprintf(screen," %d",histo[i]);
          fprintf(screen,"\n");
        }
        if (logfile) {
          fprintf(logfile,"Surfs:    %g ave %g max %g min\n",ave,max,min);
          fprintf(logfile,"Histogram:");
          for (i = 0; i < 10; i++) fprintf(logfile," %d",histo[i]);
          fprintf(logfile,"\n");
        }
      }

      tmp = surf->nghost;
      stats(1,&tmp,&ave,&max,&min,10,histo);
      if (me == 0) {
        if (screen) {
          fprintf(screen,"GhostSurf: %g ave %g max %g min\n",ave,max,min);
          fprintf(screen,"Histogram:");
          for (i = 0; i < 10; i++) fprintf(screen," %d",histo[i]);
          fprintf(screen,"\n");
        }
        if (logfile) {
          fprintf(logfile,"GhostSurf: %g ave %g max %g min\n",ave,max,min);
          fprintf(logfile,"Histogram:");
          for (i = 0; i < 10; i++) fprintf(logfile," %d",histo[i]);
          fprintf(logfile,"\n");
        }
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

/* ---------------------------------------------------------------------- */

void mpi_timings(const char *label, Timer *t, int tt,
                        MPI_Comm world, const int nprocs,
                        const int me, double time_loop, FILE *scr, FILE *log)
{
  double tmp, time_max, time_min, time_sq;
  double time = t->array[tt];

  MPI_Allreduce(&time,&time_min,1,MPI_DOUBLE,MPI_MIN,world);
  MPI_Allreduce(&time,&time_max,1,MPI_DOUBLE,MPI_MAX,world);
  time_sq = time*time;
  MPI_Allreduce(&time,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  time = tmp/nprocs;
  MPI_Allreduce(&time_sq,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  time_sq = tmp/nprocs;

  // % variance from the average as measure of load imbalance
  if ((time > 0.001) && ((time_sq/time - time) > 1.0e-10))
    time_sq = sqrt(time_sq/time - time)*100.0;
  else
    time_sq = 0.0;


  if (me == 0) {
    tmp = time/time_loop*100.0;
    const char fmt[] = "%-8s|%- 12.5g|%- 12.5g|%- 12.5g|%6.1f |%6.2f\n";
    if (scr)
      fprintf(scr,fmt,label,time_min,time,time_max,time_sq,tmp);
    if (log)
      fprintf(log,fmt,label,time_min,time,time_max,time_sq,tmp);
  }
}

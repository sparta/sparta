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
#include "update.h"
#include "particle.h"
#include "grid.h"
#include "domain.h"
#include "collide.h"
#include "comm.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace DSMC_NS;

enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};

/* ---------------------------------------------------------------------- */

Update::Update(DSMC *dsmc) : Pointers(dsmc)
{
  dt = 1.0;

  fnum = 1.0;
  nrho = 1.0;
  vstream[0] = vstream[1] = vstream[2] = 0.0;
  temp_thermal = 300.0;

  maxmigrate = 0;
  mlist = NULL;
}

/* ---------------------------------------------------------------------- */

Update::~Update()
{
  memory->destroy(mlist);
}

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

  nmove = 0;
  ncellcross = 0;
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

    // create particles

    // move particles

    timer->stamp();
    move();
    timer->stamp(TIME_MOVE);

    // communicate particles

    timer->stamp();
    comm->migrate(nmigrate,mlist);
    timer->stamp(TIME_COMM);

    if (collide) {
      timer->stamp();
      particle->sort();
      timer->stamp(TIME_SORT);

      timer->stamp();
      collide->collisions();
      timer->stamp(TIME_COLLIDE);
    }

    // output

    // sanity check on particles in correct cells

    //check();
  }
}

/* ----------------------------------------------------------------------
   advect particles thru grid
------------------------------------------------------------------------- */

void Update::move()
{
  int icell,inface,outface;
  double xnew[3];
  double *x,*v,*lo,*hi;
  int *neigh;
  double frac,newfrac;

  // extend migration list if necessary

  int nlocal = particle->nlocal;
  int maxlocal = particle->maxlocal;

  if (nlocal > maxmigrate) {
    maxmigrate = maxlocal;
    memory->destroy(mlist);
    memory->create(mlist,maxmigrate,"particle:mlist");
  }

  int dimension = domain->dimension;
  Particle::OnePart *particles = particle->particles;
  Grid::OneCell *cells = grid->cells;
  double dt = update->dt;
  int me = comm->me;

  int count = 0;
  nmigrate = 0;

  for (int i = 0; i < nlocal; i++) {

    x = particles[i].x;
    v = particles[i].v;

    xnew[0] = x[0] + dt*v[0];
    xnew[1] = x[1] + dt*v[1];
    xnew[2] = x[2] + dt*v[2];

    icell = particles[i].icell;
    lo = cells[icell].lo;
    hi = cells[icell].hi;
    neigh = cells[icell].neigh;
    inface = INTERIOR;
    count++;

    // advect particle from cell to cell until the move is done

    while (1) {

      // check if particle crosses any cell face
      // frac = fraction of move completed before hitting cell face
      // this section should be as efficient as possible,
      // since most particles won't do anything else

      outface = INTERIOR;
      frac = 1.0;
      
      if (xnew[0] < lo[0] && inface != XLO) {
	frac = (lo[0]-x[0]) / (xnew[0]-x[0]);
	outface = XLO;
      } else if (xnew[0] >= hi[0] && inface != XHI) {
	frac = (hi[0]-x[0]) / (xnew[0]-x[0]);
	outface = XHI;
      }

      if (xnew[1] < lo[1] && inface != YLO) {
	newfrac = (lo[1]-x[1]) / (xnew[1]-x[1]);
	if (newfrac < frac) {
	  frac = newfrac;
	  outface = YLO;
	}
      } else if (xnew[1] >= hi[1] && inface != YHI) {
	newfrac = (hi[1]-x[1]) / (xnew[1]-x[1]);
	if (newfrac < frac) {
	  frac = newfrac;
	  outface = YHI;
	}
      }
      
      if (xnew[2] < lo[2] && inface != ZLO) {
	newfrac = (lo[2]-x[2]) / (xnew[2]-x[2]);
	if (newfrac < frac) {
	  frac = newfrac;
	  outface = ZLO;
	}
      } else if (xnew[2] >= hi[2] && inface != ZHI) {
	newfrac = (hi[2]-x[2]) / (xnew[2]-x[2]);
	if (newfrac < frac) {
	  frac = newfrac;
	  outface = ZHI;
	}
      }

      // particle is interior to cell

      if (outface == INTERIOR) break;

      // set particle position exactly on face of cell
      
      x[0] += frac * (xnew[0]-x[0]);
      x[1] += frac * (xnew[1]-x[1]);
      x[2] += frac * (xnew[2]-x[2]);
      
      if (outface == XLO) x[0] = lo[0];
      else if (outface == XHI) x[0] = hi[0]; 
      else if (outface == YLO) x[1] = lo[1];
      else if (outface == YHI) x[1] = hi[1]; 
      else if (outface == ZLO) x[2] = lo[2];
      else if (outface == ZHI) x[2] = hi[2]; 
      
      // assign new grid cell and new inface
      // if particle crosses global boundary:
      // reflect velocity and final position, remain in same cell
      
      if (neigh[outface] < 0) {
	inface = outface;
	if (outface == XLO) {
	  xnew[0] = lo[0] + (lo[0]-xnew[0]);
	  v[0] = -v[0];
	} else if (outface == XHI) {
	  xnew[0] = hi[0] - (xnew[0]-hi[0]);
	  v[0] = -v[0];
	} else if (outface == YLO) {
	  xnew[1] = lo[1] + (lo[1]-xnew[1]);
	  v[1] = -v[1];
	} else if (outface == YHI) {
	  xnew[1] = hi[1] - (xnew[1]-hi[1]);
	  v[1] = -v[1];
	} else if (outface == ZLO) {
	  xnew[2] = lo[2] + (lo[2]-xnew[2]);
	  v[2] = -v[2];
	} else if (outface == ZHI) {
	  xnew[2] = hi[2] - (xnew[2]-hi[2]);
	  v[2] = -v[2];
	}
      } else {
	icell = neigh[outface];
	lo = cells[icell].lo;
	hi = cells[icell].hi;
	neigh = cells[icell].neigh;
	
	if (outface == XLO) inface = XHI;
	else if (outface == XHI) inface = XLO;
	else if (outface == YLO) inface = YHI;
	else if (outface == YHI) inface = YLO;
	else if (outface == ZLO) inface = ZHI;
	else if (outface == ZHI) inface = ZLO;
	count++;
      }
    }

    // update final particle position and assign to new grid cell
    // add to migrate list if I don't own new cell

    x[0] = xnew[0];
    x[1] = xnew[1];
    x[2] = xnew[2];
    particles[i].icell = icell;
    if (cells[icell].proc != me) mlist[nmigrate++] = i;
  }

  nmove += nlocal;
  ncellcross += count;
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

/* ----------------------------------------------------------------------
   set global properites via global command in input script
------------------------------------------------------------------------- */

void Update::global(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal global command");

  if (strcmp(arg[0],"fnum") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal global command");
    fnum = atof(arg[1]);
    if (fnum <= 0.0) error->all(FLERR,"Illegal global command");
  } else if (strcmp(arg[0],"nrho") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal global command");
    nrho = atof(arg[1]);
    if (nrho <= 0.0) error->all(FLERR,"Illegal global command");
  } else if (strcmp(arg[0],"vstream") == 0) {
    if (narg != 4) error->all(FLERR,"Illegal global command");
    vstream[0] = atof(arg[1]);
    vstream[1] = atof(arg[2]);
    vstream[2] = atof(arg[3]);
  } else if (strcmp(arg[0],"temp") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal global command");
    temp_thermal = atof(arg[1]);
    if (temp_thermal <= 0.0) error->all(FLERR,"Illegal global command");
  } else error->all(FLERR,"Illegal global command");
}

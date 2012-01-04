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
#include "modify.h"
#include "domain.h"
#include "comm.h"
#include "collide.h"
#include "grid.h"
#include "output.h"
#include "random_mars.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace DSMC_NS;

enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};       // same as in Domain
enum{PERIODIC,OUTFLOW,SPECULAR};              // same as in Domain

/* ---------------------------------------------------------------------- */

Update::Update(DSMC *dsmc) : Pointers(dsmc)
{
  dt = 1.0;

  kboltz = 1.38066e-23;
  mvv2e = 1.0;

  fnum = 1.0;
  nrho = 1.0;
  vstream[0] = vstream[1] = vstream[2] = 0.0;
  temp_thermal = 300.0;

  maxmigrate = 0;
  mlist = NULL;

  ranmaster = new RanMars(dsmc);

  faceflip[XLO] = XHI;
  faceflip[XHI] = XLO;
  faceflip[YLO] = YHI;
  faceflip[YHI] = YLO;
  faceflip[ZLO] = ZHI;
  faceflip[ZHI] = ZLO;
}

/* ---------------------------------------------------------------------- */

Update::~Update()
{
  memory->destroy(mlist);
  delete ranmaster;
}

/* ---------------------------------------------------------------------- */

void Update::init()
{
  if (domain->dimension == 3) move = &Update::move3d;
  else if (domain->dimension == 2) move = &Update::move2d;
}

/* ---------------------------------------------------------------------- */

void Update::setup()
{
  nmove = 0;
  ncellcross = 0;

  output->setup(1);
}

/* ---------------------------------------------------------------------- */

void Update::run(int nsteps)
{
  int n_start_of_step = modify->n_start_of_step;
  int n_end_of_step = modify->n_end_of_step;

  // loop over timesteps

  for (int i = 0; i < nsteps; i++) {

    ntimestep++;

    // start of step fixes

    if (n_start_of_step) modify->start_of_step();

    // move particles

    timer->stamp();
    (this->*move)();
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

    // diagnostic fixes

    if (n_end_of_step) modify->end_of_step();

    // all output

    if (ntimestep == output->next) {
      timer->stamp();
      output->write(ntimestep);
      timer->stamp(TIME_OUTPUT);
    }
  }
}

/* ----------------------------------------------------------------------
   advect particles thru grid in 3d manner
------------------------------------------------------------------------- */

void Update::move3d()
{
  int icell,inface,outface,outflag;
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
    outflag = 0;
    count++;

    // advect particle from cell to cell until single-step move is done

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

      // particle stays interior to cell

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
      
      // if cell has neighbor cell, move into that cell
      // else enforce global boundary conditions

      if (neigh[outface] >= 0) {
	icell = neigh[outface];
	lo = cells[icell].lo;
	hi = cells[icell].hi;
	neigh = cells[icell].neigh;
	inface = faceflip[outface];
	count++;

      } else {
        outflag = domain->boundary(outface,icell,x,xnew,v);
	if (outflag == OUTFLOW) break;
	else if (outflag == SPECULAR) inface = outface;
	else if (outflag == PERIODIC) {
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
    }

    // update final particle position
    // if OUTFLOW particle, set icell to -1 and add to migrate list,
    //   which will delete it
    // else reset icell and add to migrate list if I don't own new cell

    x[0] = xnew[0];
    x[1] = xnew[1];
    x[2] = xnew[2];

    if (outflag == OUTFLOW) {
      particles[i].icell = -1;
      mlist[nmigrate++] = i;
    } else {
      particles[i].icell = icell;
      if (cells[icell].proc != me) mlist[nmigrate++] = i;
    }
  }

  nmove += nlocal;
  ncellcross += count;
}

/* ----------------------------------------------------------------------
   advect particles thru grid in 2d manner
------------------------------------------------------------------------- */

void Update::move2d()
{
  int icell,inface,outface,outflag;
  double xnew[2];
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

    icell = particles[i].icell;
    lo = cells[icell].lo;
    hi = cells[icell].hi;
    neigh = cells[icell].neigh;
    inface = INTERIOR;
    outflag = 0;
    count++;

    // advect particle from cell to cell until single-step move is done

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
      
      // particle stays interior to cell

      if (outface == INTERIOR) break;

      // set particle position exactly on edge of cell
      
      x[0] += frac * (xnew[0]-x[0]);
      x[1] += frac * (xnew[1]-x[1]);
      
      if (outface == XLO) x[0] = lo[0];
      else if (outface == XHI) x[0] = hi[0]; 
      else if (outface == YLO) x[1] = lo[1];
      else if (outface == YHI) x[1] = hi[1]; 
      
      // if cell has neighbor cell, move into that cell
      // else enforce global boundary conditions

      if (neigh[outface] >= 0) {
	icell = neigh[outface];
	lo = cells[icell].lo;
	hi = cells[icell].hi;
	neigh = cells[icell].neigh;
	inface = faceflip[outface];
	count++;

      } else {
        outflag = domain->boundary(outface,icell,x,xnew,v);
	if (outflag == OUTFLOW) break;
	else if (outflag == SPECULAR) inface = outface;
	else if (outflag == PERIODIC) {
	  lo = cells[icell].lo;
	  hi = cells[icell].hi;
	  neigh = cells[icell].neigh;

	  if (outface == XLO) inface = XHI;
	  else if (outface == XHI) inface = XLO;
	  else if (outface == YLO) inface = YHI;
	  else if (outface == YHI) inface = YLO;
	  count++;
	}
      }
    }

    // update final particle position
    // if OUTFLOW particle, set icell to -1 and add to migrate list,
    //   which will delete it
    // else reset icell and add to migrate list if I don't own new cell

    x[0] = xnew[0];
    x[1] = xnew[1];

    if (outflag == OUTFLOW) {
      particles[i].icell = -1;
      mlist[nmigrate++] = i;
    } else {
      particles[i].icell = icell;
      if (cells[icell].proc != me) mlist[nmigrate++] = i;
    }
  }

  nmove += nlocal;
  ncellcross += count;
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

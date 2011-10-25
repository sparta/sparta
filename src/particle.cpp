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
#include "stdlib.h"
#include "particle.h"
#include "domain.h"
#include "grid.h"
#include "update.h"
#include "comm.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace DSMC_NS;

#define DELTA 10
enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};

/* ---------------------------------------------------------------------- */

Particle::Particle(DSMC *dsmc) : Pointers(dsmc)
{
  nglobal = 0;
  nlocal = maxlocal = 0;
  particles = NULL;
}

/* ---------------------------------------------------------------------- */

Particle::~Particle()
{
  memory->sfree(particles);
}

/* ----------------------------------------------------------------------
   create particles
   called from input script
------------------------------------------------------------------------- */

void Particle::create(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal create_particles command");

  int n = atoi(arg[0]);
  if (n < 0) error->all(FLERR,"Illegal create_particles command");
  int seed = atoi(arg[1]);
  RanPark *random = new RanPark(dsmc,seed);

  // add N particles within simulation box
  // grid->which_cell() returns global grid cell index the particle is in
  // if I own that grid cell, store particle

  int dimension = domain->dimension;
  double xlo = domain->boxlo[0];
  double ylo = domain->boxlo[1];
  double zlo = domain->boxlo[2];
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  int icell;
  double x,y,z;
  int me = comm->me;

  int id = nglobal;
  nglobal += n;

  while (id < nglobal) {
    id++;
    x = xlo + random->uniform()*xprd;
    y = ylo + random->uniform()*yprd;
    z = zlo + random->uniform()*zprd;
    if (dimension == 2) z = 0.0;

    icell = grid->which_cell(x,y,z);

    if (grid->cells[icell].proc == me) {
      if (nlocal == maxlocal) {
	maxlocal += DELTA;
	particles = (OnePart *)
	  memory->srealloc(particles,maxlocal*sizeof(OnePart),
			   "particle:particles");
      }

      particles[nlocal].id = id;
      particles[nlocal].type = 1;
      particles[nlocal].icell = icell;
      particles[nlocal].x[0] = x;
      particles[nlocal].x[1] = y;
      particles[nlocal].x[2] = z;
      particles[nlocal].v[0] = 0.0;
      particles[nlocal].v[1] = 0.0;
      particles[nlocal].v[2] = 0.0;
      nlocal++;
    }
  }

  delete random;

  // print stats

  if (me == 0) {
    if (screen) fprintf(screen,"Created %d particles\n",n);
    if (logfile) fprintf(logfile,"Created %d particles\n",n);
  }
}

/* ---------------------------------------------------------------------- */

void Particle::move()
{
  int icell,inface,outface;
  double xnew[3];
  double *x,*v,*lo,*hi;
  int *neigh;
  double frac,newfrac;

  int dimension = domain->dimension;
  Grid::OneCell *cells = grid->cells;
  double dt = update->dt;

  int count = 0;

  for (int i = 0; i < nlocal; i++) {

    x = particles[i].x;
    v = particles[i].v;

    xnew[0] = x[0] + dt*v[0];
    xnew[1] = x[1] + dt*v[1];
    if (dimension == 3) xnew[2] = x[2] + dt*v[2];

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
    
    x[0] = xnew[0];
    x[1] = xnew[1];
    x[2] = xnew[2];
    particles[i].icell = icell;
  }

  cellcount += count;
}

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

#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "update.h"
#include "math_const.h"
#include "particle.h"
#include "modify.h"
#include "domain.h"
#include "comm.h"
#include "collide.h"
#include "grid.h"
#include "surf.h"
#include "surf_collide.h"
#include "output.h"
#include "geometry.h"
#include "random_mars.h"
#include "random_park.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace DSMC_NS;

enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};       // same as Domain
enum{PERIODIC,OUTFLOW,REFLECT,SURFACE};       // same as Domain
enum{OUTSIDE,INSIDE,ONSURF2OUT,ONSURF2IN};    // same as Geometry

//#define MOVE_DEBUG 1            // un-comment to debug one particle's motion
#define MOVE_DEBUG_PROC 0       // proc that owns particle
#define MOVE_DEBUG_PARTICLE 0   // index of particle on owning proc

/* ---------------------------------------------------------------------- */

Update::Update(DSMC *dsmc) : Pointers(dsmc)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  ntimestep = 0;
  runflag = 0;

  unit_style = NULL;
  set_units("si");

  fnum = 1.0;
  nrho = 1.0;
  vstream[0] = vstream[1] = vstream[2] = 0.0;
  temp_thermal = 300.0;

  maxmigrate = 0;
  mlist = NULL;

  ranmaster = new RanMars(dsmc);
  random = NULL;

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
  delete [] unit_style;
  memory->destroy(mlist);
  delete ranmaster;
  delete random;
}

/* ---------------------------------------------------------------------- */

void Update::init()
{
  if (domain->dimension == 3) {
    if (surf->surf_exist) move = &Update::move3d_surface;
    else move = &Update::move3d;
  } else if (domain->dimension == 2) {
    if (surf->surf_exist) move = &Update::move2d_surface;
    else move = &Update::move2d;
  }

  if (random == NULL) {
    random = new RanPark(ranmaster->uniform());
    double seed = ranmaster->uniform();
    random->reset(seed,comm->me,100);
  }
}

/* ---------------------------------------------------------------------- */

void Update::set_units(const char *style)
{
  // physical constants from:
  // http://physics.nist.gov/cuu/Constants/Table/allascii.txt
  
  if (strcmp(style,"cgs") == 0) {
    boltz = 1.3806504e-16;
    mvv2e = 1.0;

    dt = 1.0;

  } else if (strcmp(style,"si") == 0) {
    boltz = 1.380658e-23;
    mvv2e = 1.0;

    dt = 1.0;
    
  } else error->all(FLERR,"Illegal units command");

  delete [] unit_style;
  int n = strlen(style) + 1;
  unit_style = new char[n];
  strcpy(unit_style,style);
}

/* ---------------------------------------------------------------------- */

void Update::setup()
{
  // initialize counters in case stats outputs them
  // initialize running stats

  ntouch_one = ncomm_one = 0;
  nboundary_one = nexit_one = 0;
  nscheck_one = nscollide_one = 0;

  nmove_running = ntouch_running = ncomm_running = 0;
  nboundary_running = nexit_running = 0;
  nscheck_running = nscollide_running = 0;

  // initial output

  output->setup(1);
}

/* ---------------------------------------------------------------------- */

void Update::run(int nsteps)
{
  int n_start_of_step = modify->n_start_of_step;
  int n_end_of_step = modify->n_end_of_step;
  int dynamic = 0;
  
  // loop over timesteps

  for (int i = 0; i < nsteps; i++) {

    ntimestep++;

    // start of step fixes

    ncurrent = particle->nlocal;
    if (n_start_of_step) modify->start_of_step();

    // move particles

    timer->stamp();
    //if (dynamic) domain->dynamic();
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
   check for surface collisions
------------------------------------------------------------------------- */

void Update::move3d_surface()
{
}

/* ----------------------------------------------------------------------
   advect particles thru grid in 3d manner
   no check for surface collisions
------------------------------------------------------------------------- */

void Update::move3d()
{
  int icell,inface,outface,outflag;
  double xnew[3];
  double *x,*v,*lo,*hi;
  int *neigh;
  double dtfrac,frac,newfrac;

  Particle::Species *species = particle->species;

  // extend migration list if necessary

  int nlocal = particle->nlocal;
  int maxlocal = particle->maxlocal;

  if (nlocal > maxmigrate) {
    maxmigrate = maxlocal;
    memory->destroy(mlist);
    memory->create(mlist,maxmigrate,"particle:mlist");
  }

  // counters

  ntouch_one = ncomm_one = 0;
  nboundary_one = nexit_one = 0;
  nscheck_one = nscollide_one = 0;
  nmigrate = 0;

  // loop over all my particles

  Particle::OnePart *particles = particle->particles;
  Grid::OneCell *cells = grid->cells;
  double dt = update->dt;

  for (int i = 0; i < nlocal; i++) {

    x = particles[i].x;
    v = particles[i].v;

    if (i < ncurrent) {
      xnew[0] = x[0] + dt*v[0];
      xnew[1] = x[1] + dt*v[1];
      xnew[2] = x[2] + dt*v[2];
    } else {
      dtfrac = dt*random->uniform();
      xnew[0] = x[0] + dtfrac*v[0];
      xnew[1] = x[1] + dtfrac*v[1];
      xnew[2] = x[2] + dtfrac*v[2];
    }

    icell = particles[i].icell;
    lo = cells[icell].lo;
    hi = cells[icell].hi;
    neigh = cells[icell].neigh;
    inface = INTERIOR;
    outflag = -1;
    ntouch_one++;

    // advect particle from cell to cell until single-step move is done

    while (1) {

      // check if particle crosses any cell face
      // frac = fraction of move completed before hitting cell face
      // this section should be as efficient as possible,
      // since most particles won't do anything else
      // used to do this, but don't now think is necessary:
      // if (xnew[0] < lo[0] && inface != XLO) {

      outface = INTERIOR;
      frac = 1.0;
      
      if (xnew[0] < lo[0]) {
	frac = (lo[0]-x[0]) / (xnew[0]-x[0]);
	outface = XLO;
      } else if (xnew[0] >= hi[0]) {
	frac = (hi[0]-x[0]) / (xnew[0]-x[0]);
	outface = XHI;
      }

      if (xnew[1] < lo[1]) {
	newfrac = (lo[1]-x[1]) / (xnew[1]-x[1]);
	if (newfrac < frac) {
	  frac = newfrac;
	  outface = YLO;
	}
      } else if (xnew[1] >= hi[1]) {
	newfrac = (hi[1]-x[1]) / (xnew[1]-x[1]);
	if (newfrac < frac) {
	  frac = newfrac;
	  outface = YHI;
	}
      }
      
      if (xnew[2] < lo[2]) {
	newfrac = (lo[2]-x[2]) / (xnew[2]-x[2]);
	if (newfrac < frac) {
	  frac = newfrac;
	  outface = ZLO;
	}
      } else if (xnew[2] >= hi[2]) {
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
	ntouch_one++;

      } else {
        outflag = domain->collide(&particles[i],outface,icell,xnew);
	if (outflag == OUTFLOW) {
	  nexit_one++;
	  break;
	} else if (outflag == PERIODIC) {
	  lo = cells[icell].lo;
	  hi = cells[icell].hi;
	  neigh = cells[icell].neigh;
	  inface = faceflip[outface];
	  ntouch_one++;
	} else {
	  inface = outface;
	  nboundary_one++;
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
      if (cells[icell].proc != me) {
	mlist[nmigrate++] = i;
	ncomm_one++;
      }
    }
  }

  // accumulate running totals

  nmove_running += nlocal;
  ntouch_running += ntouch_one;
  ncomm_running += ncomm_one;
  nboundary_running += nboundary_one;
  nexit_running += nexit_one;
  nscheck_running += nscheck_one;
  nscollide_running += nscollide_one;
}

/* ----------------------------------------------------------------------
   advect particles thru grid in 2d manner
   check for surface collisions
------------------------------------------------------------------------- */

void Update::move2d_surface()
{
  bool hitflag;
  int i,m,icell,inface,outface,outflag,isurf,exclude;
  int side,minside,minsurf,nsurf,cflag;
  int *neigh;
  double dtremain,dtfrac,frac,newfrac,param,minparam;
  double *x,*v,*lo,*hi;
  Surf::Line *line;

  // xnew,xc passed to geometry routines which use or set 3rd component
  // xhold,minxc used locally w/out 3rd component

  double xnew[3],xc[3];      
  double xhold[2],minxc[2];
  xnew[2] = 0.0;

  // extend migration list if necessary

  int nlocal = particle->nlocal;
  int maxlocal = particle->maxlocal;

  if (nlocal > maxmigrate) {
    maxmigrate = maxlocal;
    memory->destroy(mlist);
    memory->create(mlist,maxmigrate,"particle:mlist");
  }

  // counters

  ntouch_one = ncomm_one = 0;
  nboundary_one = nexit_one = 0;
  nscheck_one = nscollide_one = 0;
  nmigrate = 0;

  // loop over all my particles

  Particle::OnePart *particles = particle->particles;
  Grid::OneCell *cells = grid->cells;
  int **csurfs = grid->csurfs;
  Surf::Line *lines = surf->lines;
  Surf::Point *pts = surf->pts;

  double dt = update->dt;

  for (i = 0; i < nlocal; i++) {

    x = particles[i].x;
    v = particles[i].v;

    if (i < ncurrent) {
      xnew[0] = x[0] + dt*v[0];
      xnew[1] = x[1] + dt*v[1];
    } else {
      dtfrac = dt*random->uniform();
      xnew[0] = x[0] + dtfrac*v[0];
      xnew[1] = x[1] + dtfrac*v[1];
    }

    icell = particles[i].icell;
    nsurf = cells[icell].nsurf;
    lo = cells[icell].lo;
    hi = cells[icell].hi;
    neigh = cells[icell].neigh;
    dtremain = dt;
    inface = INTERIOR;
    outflag = -1;
    exclude = -1;
    ntouch_one++;

    // advect particle from cell to cell until move is done

    while (1) {

#ifdef MOVE_DEBUG
      if (i == MOVE_DEBUG_PARTICLE && me == MOVE_DEBUG_PROC) 
	printf("PARTICLE %d: %d %d: %g %g: %g %g\n",MOVE_DEBUG_PARTICLE,
	       update->ntimestep,nsurf,x[0],x[1],xnew[0],xnew[1]);
#endif

      // check if particle crosses any cell face
      // frac = fraction of move completed before hitting cell face
      // this section should be as efficient as possible,
      // since most particles won't do anything else

      outface = INTERIOR;
      frac = 1.0;
      
      if (xnew[0] < lo[0]) {
	frac = (lo[0]-x[0]) / (xnew[0]-x[0]);
	outface = XLO;
      } else if (xnew[0] >= hi[0]) {
	frac = (hi[0]-x[0]) / (xnew[0]-x[0]);
	outface = XHI;
      }

      if (xnew[1] < lo[1]) {
	newfrac = (lo[1]-x[1]) / (xnew[1]-x[1]);
	if (newfrac < frac) {
	  frac = newfrac;
	  outface = YLO;
	}
      } else if (xnew[1] >= hi[1]) {
	newfrac = (hi[1]-x[1]) / (xnew[1]-x[1]);
	if (newfrac < frac) {
	  frac = newfrac;
	  outface = YHI;
	}
      }
      
      // particle crosses cell face

      if (outface != INTERIOR) {
	xhold[0] = xnew[0];
	xhold[1] = xnew[1];

	xnew[0] = x[0] + frac*(xnew[0]-x[0]);
	xnew[1] = x[1] + frac*(xnew[1]-x[1]);
      
	if (outface == XLO) xnew[0] = lo[0];
	else if (outface == XHI) xnew[0] = hi[0]; 
	else if (outface == YLO) xnew[1] = lo[1];
	else if (outface == YHI) xnew[1] = hi[1]; 
      }

      // check for collisions with lines in cell
      // find 1st surface hit via minparam
      // if collision occurs, perform collision with surface model
      // reset x,v,xnew,dtremain and continue particle trajectory

      if (nsurf) {
	cflag = 0;
	minparam = 2.0;
	for (m = 0; m < nsurf; m++) {
	  isurf = csurfs[icell][m];
	  if (isurf == exclude) continue;
	  line = &lines[isurf];
	  hitflag = Geometry::
	    line_line_intersect(x,xnew,
				pts[line->p1].x,pts[line->p2].x,line->norm,
				xc,param,side);

#ifdef MOVE_DEBUG
	  if (hitflag && i == MOVE_DEBUG_PARTICLE && me == MOVE_DEBUG_PROC)
	    printf("SURF COLLIDE: %d %d %d %d: P1 %g %g: P2 %g %g: L1 %g %g: "
		   "L2 %g %g: LN %g %g: XC %g %g: Param %g: Side %d: DTR %g\n",
		   MOVE_DEBUG_PARTICLE,icell,nsurf,isurf,
		   x[0],x[1],xnew[0],xnew[1],
		   pts[line->p1].x[0],pts[line->p1].x[1],
		   pts[line->p2].x[0],pts[line->p2].x[1],
		   line->norm[0],line->norm[1],xc[0],xc[1],param,side,
		   dtremain);
#endif

	  if (hitflag && side != ONSURF2OUT && param <= minparam) {
	    if (param == minparam && side == INSIDE) continue;
	    cflag = 1;
	    minparam = param;
	    minside = side;
	    minsurf = isurf;
	    minxc[0] = xc[0];
	    minxc[1] = xc[1];
	  }
	}
	nscheck_one += nsurf;

	if (cflag) {
	  if (minside == INSIDE) 
	    error->one(FLERR,"Particle hit inside of surface element");
	  x[0] = minxc[0];
	  x[1] = minxc[1];
	  line = &lines[minsurf];
	  surf->sc[line->isc]->collide(&particles[i],line->norm);
	  dtremain *= 1.0 - minparam*frac;
	  xnew[0] = x[0] + dtremain*v[0];
	  xnew[1] = x[1] + dtremain*v[1];
	  exclude = minsurf;
	  nscollide_one++;

#ifdef MOVE_DEBUG
	  if (i == MOVE_DEBUG_PARTICLE && me == MOVE_DEBUG_PROC)
	    printf("POST COLLISION %d: %g %g: %g %g: %g %g %g\n",
		   MOVE_DEBUG_PARTICLE,
		   x[0],x[1],xnew[0],xnew[1],minparam,frac,dtremain);
#endif
	  continue;
	}
      }

      // no surface collision & no cell crossing, so done

      if (outface == INTERIOR) break;

      // cell crossing: reset dtremain, restore xnew

      exclude = -1;
      dtremain *= 1.0-frac;
      xnew[0] = xhold[0];
      xnew[1] = xhold[1];

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
	nsurf = cells[icell].nsurf;
	lo = cells[icell].lo;
	hi = cells[icell].hi;
	neigh = cells[icell].neigh;
	inface = faceflip[outface];
	ntouch_one++;
 
      } else {
        outflag = domain->collide(&particles[i],outface,icell,xnew);
	if (outflag == OUTFLOW) {
	  nexit_one++;
	  break;
	} else if (outflag == PERIODIC) {
	  nsurf = cells[icell].nsurf;
	  lo = cells[icell].lo;
	  hi = cells[icell].hi;
	  neigh = cells[icell].neigh;
	  inface = faceflip[outface];
	  ntouch_one++;
	} else {
	  inface = outface;
	  nboundary_one++;
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
      if (cells[icell].proc != me) {
	mlist[nmigrate++] = i;
	ncomm_one++;
      }
    }
  }

  // accumulate running totals

  nmove_running += nlocal;
  ntouch_running += ntouch_one;
  ncomm_running += ncomm_one;
  nboundary_running += nboundary_one;
  nexit_running += nexit_one;
  nscheck_running += nscheck_one;
  nscollide_running += nscollide_one;
}

/* ----------------------------------------------------------------------
   advect particles thru grid in 2d manner
   no check for surface collisions
------------------------------------------------------------------------- */

void Update::move2d()
{
  int icell,inface,outface,outflag;
  double xnew[3];
  double *x,*v,*lo,*hi;
  int *neigh;
  double dtfrac,frac,newfrac;

  // needed for MathExtra calls

  xnew[2] = 0.0;

  // extend migration list if necessary

  int nlocal = particle->nlocal;
  int maxlocal = particle->maxlocal;

  if (nlocal > maxmigrate) {
    maxmigrate = maxlocal;
    memory->destroy(mlist);
    memory->create(mlist,maxmigrate,"particle:mlist");
  }

  // counters

  ntouch_one = ncomm_one = 0;
  nboundary_one = nexit_one = 0;
  nscheck_one = nscollide_one = 0;
  nmigrate = 0;

  // loop over all my particles

  Particle::OnePart *particles = particle->particles;
  Grid::OneCell *cells = grid->cells;
  double dt = update->dt;

  for (int i = 0; i < nlocal; i++) {

    x = particles[i].x;
    v = particles[i].v;

    if (i < ncurrent) {
      xnew[0] = x[0] + dt*v[0];
      xnew[1] = x[1] + dt*v[1];
    } else {
      dtfrac = dt*random->uniform();
      xnew[0] = x[0] + dtfrac*v[0];
      xnew[1] = x[1] + dtfrac*v[1];
    }

    icell = particles[i].icell;
    lo = cells[icell].lo;
    hi = cells[icell].hi;
    neigh = cells[icell].neigh;
    inface = INTERIOR;
    outflag = -1;
    ntouch_one++;

    // advect particle from cell to cell until single-step move is done

    while (1) {

      // check if particle crosses any cell face
      // frac = fraction of move completed before hitting cell face
      // this section should be as efficient as possible,
      // since most particles won't do anything else

      outface = INTERIOR;
      frac = 1.0;
      
      if (xnew[0] < lo[0]) {
	frac = (lo[0]-x[0]) / (xnew[0]-x[0]);
	outface = XLO;
      } else if (xnew[0] >= hi[0]) {
	frac = (hi[0]-x[0]) / (xnew[0]-x[0]);
	outface = XHI;
      }

      if (xnew[1] < lo[1]) {
	newfrac = (lo[1]-x[1]) / (xnew[1]-x[1]);
	if (newfrac < frac) {
	  frac = newfrac;
	  outface = YLO;
	}
      } else if (xnew[1] >= hi[1]) {
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
	ntouch_one++;
 
      } else {
        outflag = domain->collide(&particles[i],outface,icell,xnew);
	if (outflag == OUTFLOW) {
	  nexit_one++;
	  break;
	} else if (outflag == PERIODIC) {
	  lo = cells[icell].lo;
	  hi = cells[icell].hi;
	  neigh = cells[icell].neigh;
	  inface = faceflip[outface];
	  ntouch_one++;
	} else {
	  inface = outface;
	  nboundary_one++;
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
      if (cells[icell].proc != me) {
	mlist[nmigrate++] = i;
	ncomm_one++;
      }
    }
  }

  // accumulate running totals

  nmove_running += nlocal;
  ntouch_running += ntouch_one;
  ncomm_running += ncomm_one;
  nboundary_running += nboundary_one;
  nexit_running += nexit_one;
  nscheck_running += nscheck_one;
  nscollide_running += nscollide_one;
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

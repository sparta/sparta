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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "domain.h"
#include "update.h"
#include "grid.h"
#include "comm.h"
#include "particle.h"
#include "surf.h"
#include "surf_collide.h"
#include "comm.h"
#include "error.h"

using namespace DSMC_NS;

enum{PERIODIC,OUTFLOW,REFLECT,SURFACE};     // same as Update,
                                            // DumpMolecule, FixInflow
enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};     // same as Update, FixInflow

/* ---------------------------------------------------------------------- */

Domain::Domain(DSMC *dsmc) : Pointers(dsmc)
{
  box_exist = 0;
  dimension = 3;

  // surface normals of 6 box faces pointed inward towards particles

  norm[XLO][0] =  1.0; norm[XLO][1] =  0.0; norm[XLO][2] =  0.0;
  norm[XHI][0] = -1.0; norm[XHI][1] =  0.0; norm[XHI][2] =  0.0;
  norm[YLO][0] =  0.0; norm[YLO][1] =  1.0; norm[YLO][2] =  0.0;
  norm[YHI][0] =  0.0; norm[YHI][1] = -1.0; norm[YHI][2] =  0.0;
  norm[ZLO][0] =  0.0; norm[ZLO][1] =  0.0; norm[ZLO][2] =  1.0;
  norm[ZHI][0] =  0.0; norm[ZHI][1] =  0.0; norm[ZHI][2] = -1.0;
}

/* ----------------------------------------------------------------------
   set initial global box
   assumes boxlo/hi already set
------------------------------------------------------------------------- */

void Domain::set_initial_box()
{
  if (boxlo[0] >= boxhi[0] || boxlo[1] >= boxhi[1] || boxlo[2] >= boxhi[2])
    error->one(FLERR,"Box bounds are invalid");
}

/* ----------------------------------------------------------------------
   set global box params
   assumes boxlo/hi already set
------------------------------------------------------------------------- */

void Domain::set_global_box()
{
  prd[0] = xprd = boxhi[0] - boxlo[0];
  prd[1] = yprd = boxhi[1] - boxlo[1];
  prd[2] = zprd = boxhi[2] - boxlo[2];
}

/* ----------------------------------------------------------------------
   boundary settings from input script 
------------------------------------------------------------------------- */

void Domain::set_boundary(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Illegal boundary command");

  char c;
  int m = 0;
  for (int idim = 0; idim < 3; idim++)
    for (int iside = 0; iside < 2; iside++) {
      if (iside == 0) c = arg[idim][0];
      else if (iside == 1 && strlen(arg[idim]) == 1) c = arg[idim][0];
      else c = arg[idim][1];

      if (c == 'o') bflag[m] = OUTFLOW;
      else if (c == 'p') bflag[m] = PERIODIC;
      else if (c == 'r') bflag[m] = REFLECT;
      else if (c == 's') bflag[m] = SURFACE;
      else error->all(FLERR,"Illegal boundary command");

      m++;
    }

  if (dimension == 2 && (bflag[ZLO] != PERIODIC || bflag[ZHI] != PERIODIC))
    error->all(FLERR,"Z dimension must be periodic for 2d simulation");
      
  for (m = 0; m < 6; m += 2)
    if (bflag[m] == PERIODIC || bflag[m+1] == PERIODIC) {
      if (bflag[m] != PERIODIC || bflag[m+1] != PERIODIC)
	error->all(FLERR,"Both sides of boundary must be periodic");
    }
}

/* ----------------------------------------------------------------------
   boundary modifications from input script 
------------------------------------------------------------------------- */

void Domain::boundary_modify(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal bound_modify command");

  int face;
  if (strcmp(arg[0],"xlo") == 0) face = XLO;
  else if (strcmp(arg[0],"xhi") == 0) face = XHI;
  else if (strcmp(arg[0],"ylo") == 0) face = YLO;
  else if (strcmp(arg[0],"yhi") == 0) face = YHI;
  else if (strcmp(arg[0],"zlo") == 0) face = ZLO;
  else if (strcmp(arg[0],"zhi") == 0) face = ZHI;
  else error->all(FLERR,"Illegal bound_modify command");

  if (dimension == 2 && (face == ZLO || face == ZHI))
    error->all(FLERR,
	       "Cannot use bound_modify command "
	       "on Z dimension for 2d simulation");

  // setting surf_collide[] index here vs init()
  // assumes SurfCollide styles are static (never deleted)

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"surf") == 0) {
      if (iarg + 2 > narg) error->all(FLERR,"Illegal bound_modify command");
      surf_collide[face] = surf->find_collide(arg[iarg+1]);
      if (surf_collide[face] < 0) 
	error->all(FLERR,"Bound_modify surf_collide ID is unknown");
      iarg += 2;
    } else error->all(FLERR,"Illegal bound_modify command");
  }
}

/* ----------------------------------------------------------------------
   particle P hits global boundary on face
   called by Update::move()
   icell = current global grid cell particle is in
   xnew = final position of particle at end of move
   return boundary type of global boundary
   also update icell, xnew, particle velocity due to collision
   NOTE: tally stats?
   NOTE: code below depends on grid structure = Nx by Ny by Nz cells
------------------------------------------------------------------------- */

int Domain::collide(Particle::OnePart *p, int face, int &icell, double *xnew)
{
  switch (bflag[face]) {

  // outflow boundary, particle deleted by caller

  case OUTFLOW:
    return OUTFLOW;

  // periodic boundary
  // adjust x and xnew by periodic box length
  // assign new icell on other side of box

  case PERIODIC:
    {
      double *x = p->x;
      switch (face) {
      case XLO:
	x[0] += xprd;
	xnew[0] += xprd;
	icell += grid->nx - 1;
	break;
      case XHI:
	x[0] -= xprd;
	xnew[0] -= xprd;
	icell -= grid->nx - 1;
	break;
      case YLO:
	x[1] += yprd;
	xnew[1] += yprd;
	icell += (grid->ny-1) * grid->nx;
	break;
      case YHI:
	x[1] -= yprd;
	xnew[1] -= yprd;
	icell -= (grid->ny-1) * grid->nx;
	break;
      case ZLO:
	x[2] += zprd;
	xnew[2] += zprd;
	icell += (grid->nz-1) * grid->ny * grid->nx;
	break;
      case ZHI:
	x[2] -= zprd;
	xnew[2] -= zprd;
	icell -= (grid->nz-1) * grid->ny * grid->nx;
	break;
      }
    }

    return PERIODIC;

  // specular reflection boundary
  // adjust xnew and velocity

  case REFLECT:
    {
      double *v = p->v;
      double *lo = grid->cells[icell].lo;
      double *hi = grid->cells[icell].hi;
      int dim = face / 2;
      
      if (face % 2 == 0) {
	xnew[dim] = lo[dim] + (lo[dim]-xnew[dim]);
	v[dim] = -v[dim];
      } else {
	xnew[dim] = hi[dim] - (xnew[dim]-hi[dim]);
	v[dim] = -v[dim];
      }
    }
    
    return REFLECT;

  // treat global boundary as a surface
  // dtr = time remaining after collision
  // particle velocity is changed by surface collision model
  // reset one component of xnew due to new velocity

  case SURFACE: 
    {
      double *v = p->v;
      double *lo = grid->cells[icell].lo;
      double *hi = grid->cells[icell].hi;
      int dim = face / 2;

      if (face % 2 == 0) {
	double dtr = fabs((lo[dim]-xnew[dim])/v[dim]);
	surf->sc[surf_collide[face]]->collide(p,norm[face]);
	xnew[dim] = lo[dim] + v[dim]*dtr;
      } else {
	double dtr = fabs((xnew[dim]-hi[dim])/v[dim]);
	surf->sc[surf_collide[face]]->collide(p,norm[face]);
	xnew[dim] = hi[dim] - v[dim]*dtr;
      }
      
      return SURFACE;
    }

  }

  return 0;
}

/* ----------------------------------------------------------------------
   print box info
------------------------------------------------------------------------- */

void Domain::print_box(const char *str)
{
  if (comm->me == 0) {
    if (screen)
      fprintf(screen,"%sorthogonal box = (%g %g %g) to (%g %g %g)\n",
	      str,boxlo[0],boxlo[1],boxlo[2],boxhi[0],boxhi[1],boxhi[2]);
    if (logfile)
      fprintf(logfile,"%sorthogonal box = (%g %g %g) to (%g %g %g)\n",
	      str,boxlo[0],boxlo[1],boxlo[2],boxhi[0],boxhi[1],boxhi[2]);
  }
}

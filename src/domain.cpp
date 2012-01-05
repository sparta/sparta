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

#include "string.h"
#include "domain.h"
#include "grid.h"
#include "comm.h"
#include "error.h"

using namespace DSMC_NS;

enum{PERIODIC,OUTFLOW,SPECULAR};            // same as CreateGrid, Update,
                                            // DumpParticle, FixInflow
enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};     // same as Update, FixInflow

/* ---------------------------------------------------------------------- */

Domain::Domain(DSMC *dsmc) : Pointers(dsmc)
{
  box_exist = 0;
  dimension = 3;

  bflag[0] = bflag[1] = bflag[2] = bflag[3] = bflag[4] = bflag[5] = PERIODIC;
}

/* ---------------------------------------------------------------------- */

Domain::~Domain() {}

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
   boundary settings from the input script 
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
      else if (c == 's') bflag[m] = SPECULAR;
      else error->all(FLERR,"Illegal boundary command");

      m++;
    }

  for (m = 0; m < 6; m += 2)
    if (bflag[m] == PERIODIC || bflag[m+1] == PERIODIC) {
      if (bflag[m] != PERIODIC || bflag[m+1] != PERIODIC)
	error->all(FLERR,"Both sides of boundary must be periodic");
    }
}

/* ----------------------------------------------------------------------
   particle hits global boundary
   called by Update::move()
   periodic case is handled in Update::move() as moving into neighbor cell
   return 1 if OUTFLOW and particle is lost
------------------------------------------------------------------------- */

int Domain::boundary(int face, int &icell, double *x, double *xnew, double *v) 
{
  // tally stats?

  // outflow boundary, particle will be deleted if return 1

  if (bflag[face] == OUTFLOW) return OUTFLOW;

  // periodic boundary
  // adjust x and xnew by periodic box length
  // assign new icell on other side of box

  if (bflag[face] == PERIODIC) {
    if (face == XLO) {
      x[0] += xprd;
      xnew[0] += xprd;
      icell += grid->nx - 1;
    } else if (face == XHI) {
      x[0] -= xprd;
      xnew[0] -= xprd;
      icell -= grid->nx - 1;
    } else if (face == YLO) {
      x[1] += yprd;
      xnew[1] += yprd;
      icell += (grid->ny-1) * grid->nx;
    } else if (face == YHI) {
      x[1] -= yprd;
      xnew[1] -= yprd;
      icell -= (grid->ny-1) * grid->nx;
    } else if (face == ZLO) {
      x[2] += zprd;
      xnew[2] += zprd;
      icell += (grid->nz-1) * grid->ny * grid->nx;
    } else if (face == ZHI) {
      x[2] -= zprd;
      xnew[2] -= zprd;
      icell -= (grid->nz-1) * grid->ny * grid->nx;
    }

    return PERIODIC;
  }

  // specular reflection boundary
  // adjust xnew and velocity

  if (bflag[face] == SPECULAR) {
    double *lo = grid->cells[icell].lo;
    double *hi = grid->cells[icell].hi;

    if (face == XLO) {
      xnew[0] = lo[0] + (lo[0]-xnew[0]);
      v[0] = -v[0];
    } else if (face == XHI) {
      xnew[0] = hi[0] - (xnew[0]-hi[0]);
      v[0] = -v[0];
    } else if (face == YLO) {
      xnew[1] = lo[1] + (lo[1]-xnew[1]);
      v[1] = -v[1];
    } else if (face == YHI) {
      xnew[1] = hi[1] - (xnew[1]-hi[1]);
      v[1] = -v[1];
    } else if (face == ZLO) {
      xnew[2] = lo[2] + (lo[2]-xnew[2]);
      v[2] = -v[2];
    } else if (face == ZHI) {
      xnew[2] = hi[2] - (xnew[2]-hi[2]);
      v[2] = -v[2];
    }

    return SPECULAR;
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

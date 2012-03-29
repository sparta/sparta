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
#include "stdlib.h"
#include "domain.h"
#include "grid.h"
#include "comm.h"
#include "error.h"
#include "math.h"
#include "particle.h"
#include "update.h"
#include "comm.h"
#include "random_mars.h"
#include "random_park.h"
#include "math_const.h"

using namespace DSMC_NS;
using namespace MathConst;

enum{PERIODIC,OUTFLOW,SPECULAR,DIFFUSE};    // same as Update,
                                            // DumpMolecule, FixInflow
enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};     // same as Update, FixInflow

/* ---------------------------------------------------------------------- */

Domain::Domain(DSMC *dsmc) : Pointers(dsmc)
{

  box_exist = 0;
  dimension = 3;

  bflag[0] = bflag[1] = bflag[2] = bflag[3] = bflag[4] = bflag[5] = PERIODIC;

  // default values

  for (int i = 0; i < 6; i++) {
    acccoeff[i] = 1.0;
    twall[i] = 273.15;
  }

  // RNG for particle reflection off global boundaries

  random = NULL;
}

/* ---------------------------------------------------------------------- */

Domain::~Domain()
{
  delete random;
}

/* ---------------------------------------------------------------------- */

void Domain::init()
{
  if (random == NULL) {
    random = new RanPark(update->ranmaster->uniform());
    double seed = update->ranmaster->uniform();
    random->reset(seed,comm->me,100);
  }
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
      else if (c == 's') bflag[m] = SPECULAR;
      else if (c == 'd') bflag[m] = DIFFUSE;
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
  if (narg < 2) error->all(FLERR,"Illegal boundary_modify command");

  int face;
  if (strcmp(arg[0],"xlo") == 0) face = XLO;
  else if (strcmp(arg[0],"xhi") == 0) face = XHI;
  else if (strcmp(arg[0],"ylo") == 0) face = YLO;
  else if (strcmp(arg[0],"yhi") == 0) face = YHI;
  else if (strcmp(arg[0],"zlo") == 0) face = ZLO;
  else if (strcmp(arg[0],"zhi") == 0) face = ZHI;
  else error->all(FLERR,"Illegal boundary_modify command");

  if (dimension == 2 && (face == ZLO || face == ZHI))
    error->all(FLERR,"Cannot use boundary_modify command "
	       "on Z dimension for 2d simulation");

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"acc") == 0) {
      if (iarg + 2 > narg) error->all(FLERR,"Illegal boundary_modify command");
      acccoeff[face] = atof(arg[iarg+1]);
      if (acccoeff[face] < 0.0 || acccoeff[face] > 1.0) 
	error->all(FLERR,"Illegal boundary_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"temp") == 0) {
      if (iarg + 2 > narg) error->all(FLERR,"Illegal boundary_modify command");
      twall[face] = atof(arg[iarg+1]);
      if (twall[face] < 0.0)
	error->all(FLERR,"Illegal boundary_modify command");
      iarg += 2;
    } else error->all(FLERR,"Illegal boundary_modify command");
  }
}

/* ----------------------------------------------------------------------
   particle hits global boundary
   called by Update::move()
   periodic case is handled in Update::move() as moving into neighbor cell
   return boundary type
------------------------------------------------------------------------- */

int Domain::boundary(int face, int &icell, double *x, double *xnew, 
		     double *v, int isp) 
{
  // tally stats?

  // outflow boundary, particle deleted by caller

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
    int dim = face / 2;

    if (face % 2 == 0) {
      xnew[dim] = lo[dim] + (lo[dim]-xnew[dim]);
      v[dim] = -v[dim];
    } else {
      xnew[dim] = hi[dim] - (xnew[dim]-hi[dim]);
      v[dim] = -v[dim];
    }

    return SPECULAR;
  }

  // dtr = time remaining after collision
  // wall temprature Needs to come in through input
  // accommodation coefficient needs to come in through input

  if (bflag[face] == DIFFUSE) {
    double *lo = grid->cells[icell].lo;
    double *hi = grid->cells[icell].hi;
    int dim = face / 2;

    if (face % 2 == 0) {
      double dtr = fabs((lo[dim]-xnew[dim])/v[dim]);
      reflect(face,isp,v);
      xnew[dim] = lo[dim] + v[dim]*dtr;
    } else {
      double dtr = fabs((xnew[dim]-hi[dim])/v[dim]);
      reflect(face,isp,v);
      xnew[dim] = hi[dim] - v[dim]*dtr;
    }

    return DIFFUSE;
  }

  return 0;
}

/* ----------------------------------------------------------------------
   particle reflection off simulation box wall for DIFFUSE model
------------------------------------------------------------------------- */

void Domain::reflect(int face, int isp, double *v)
{
  // specular reflection

  if (random->uniform() > acccoeff[face]) {
    int dim = face / 2;
    v[dim] = -v[dim];

  // diffuse reflection
  // vrm = most probable speed of species isp, eqns (4.1) and (4.7)
  // generate normal velocity component, eqn (12.3)

  } else {
    Particle::Species *species = particle->species;
    double vrm = sqrt(2.0*update->boltz*twall[face]/species[isp].mass);
    v[0] = vrm * random->uniform();
    double theta1 = MY_2PI * random->uniform();
    double theta2 = MY_2PI * random->uniform();
    v[1] = vrm*sin(theta2);
    v[2] = vrm*cos(theta2);
    /*
      erot(isp);
      evib(isp);
    */
  }
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

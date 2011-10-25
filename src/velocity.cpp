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
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "velocity.h"
#include "domain.h"
#include "particle.h"
#include "comm.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace DSMC_NS;

enum{CREATE,SET};
enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

Velocity::Velocity(DSMC *dsmc) : Pointers(dsmc) {}

/* ---------------------------------------------------------------------- */

void Velocity::command(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal velocity command");

  if (domain->box_exist == 0) 
    error->all(FLERR,"Velocity command before simulation box is defined");
  if (particle->nglobal == 0)
    error->all(FLERR,"Velocity command with no particles existing");

  // parse args

  int style;
  sum_flag = 0;

  if (strcmp(arg[0],"create") == 0) {
    if (narg < 3) error->all(FLERR,"Illegal velocity command");
    style = CREATE;
    t_desired = atof(arg[1]);
    seed = atoi(arg[2]);
    options(narg-3,&arg[3]);

  } else if (strcmp(arg[0],"set") == 0) {
    if (narg < 4) error->all(FLERR,"Illegal velocity command");
    style = SET;
    vx = atof(arg[1]);
    vy = atof(arg[2]);
    vz = atof(arg[3]);
    if (domain->dimension == 2 && vz != 0.0)
      error->all(FLERR,
		 "Cannot set velocity to non-zero z value for 2d simulation");
    options(narg-4,&arg[4]);

  } else error->all(FLERR,"Illegal velocity command");

  // initialize velocities based on style

  if (style == CREATE) create();
  else if (style == SET) set();
}

/* ----------------------------------------------------------------------
   create velocities
   for now, use t_desired as a uniform delta spread
------------------------------------------------------------------------- */

void Velocity::create()
{
  if (seed <= 0) error->all(FLERR,"Illegal velocity create command");

  // different RNG for each proc

  RanPark *random = new RanPark(dsmc,seed+comm->me);

  // create randomized delta spread in velocity of each particle

  Particle::OnePart *particles = particle->particles;
  int nlocal = particle->nlocal;
  int dimension = domain->dimension;

  double *v;
  double xdelta,ydelta,zdelta;

  if (sum_flag) {
    for (int i = 0; i < nlocal; i++) {
      xdelta = 2.0*t_desired * (random->uniform() - 0.5);
      ydelta = 2.0*t_desired * (random->uniform() - 0.5);
      zdelta = 2.0*t_desired * (random->uniform() - 0.5);
      if (dimension == 2) zdelta = 0.0;
      v = particles[i].v;
      v[0] += xdelta;
      v[1] += ydelta;
      v[2] += zdelta;
    }

  } else {
    for (int i = 0; i < nlocal; i++) {
      xdelta = 2.0*t_desired * (random->uniform() - 0.5);
      ydelta = 2.0*t_desired * (random->uniform() - 0.5);
      zdelta = 2.0*t_desired * (random->uniform() - 0.5);
      if (dimension == 2) zdelta = 0.0;
      v = particles[i].v;
      v[0] = xdelta;
      v[1] = ydelta;
      v[2] = zdelta;
    }
  }

  delete random;
}

/* ----------------------------------------------------------------------
   set velocities to specified vx,vy,vz
------------------------------------------------------------------------- */

void Velocity::set()
{
  Particle::OnePart *particles = particle->particles;
  int nlocal = particle->nlocal;

  double *v;

  if (sum_flag) {
    for (int i = 0; i < nlocal; i++) {
      v = particles[i].v;
      v[0] += vx;
      v[1] += vy;
      v[2] += vz;
    }

  } else {
    for (int i = 0; i < nlocal; i++) {
      v = particles[i].v;
      v[0] = vx;
      v[1] = vy;
      v[2] = vz;
    }
  }
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of velocity input line 
------------------------------------------------------------------------- */

void Velocity::options(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"sum") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal velocity command");
      if (strcmp(arg[iarg+1],"no") == 0) sum_flag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) sum_flag = 1;
      else error->all(FLERR,"Illegal velocity command");
      iarg += 2;
    } else error->all(FLERR,"Illegal velocity command");
  }
}

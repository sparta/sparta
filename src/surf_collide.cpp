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
#include "surf_collide.h"

using namespace DSMC_NS;

/* ---------------------------------------------------------------------- */

SurfCollide::SurfCollide(DSMC *dsmc, int narg, char **arg) : Pointers(dsmc)
{
  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  n = strlen(arg[0]) + 1;
  style = new char[n];
  strcpy(style,arg[0]);
}

/* ---------------------------------------------------------------------- */

SurfCollide::~SurfCollide()
{
  delete [] id;
  delete [] style;
}

/* ---------------------------------------------------------------------- */

void SurfCollide::init()
{
}

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
#include "ctype.h"
#include "fix.h"
#include "memory.h"
#include "error.h"

using namespace DSMC_NS;

/* ---------------------------------------------------------------------- */

Fix::Fix(DSMC *dsmc, int narg, char **arg) : Pointers(dsmc)
{
  // fix ID and style
  // ID must be all alphanumeric chars or underscores

  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  for (int i = 0; i < n-1; i++)
    if (!isalnum(id[i]) && id[i] != '_')
      error->all(FLERR,"Fix ID must be alphanumeric or underscore characters");

  n = strlen(arg[1]) + 1;
  style = new char[n];
  strcpy(style,arg[1]);

  scalar_flag = vector_flag = array_flag = 0;
  per_molecule_flag = per_grid_flag = 0;

  // mask settings - same as in modify.cpp

  START_OF_STEP = 1;
  END_OF_STEP = 2;
}

/* ---------------------------------------------------------------------- */

Fix::~Fix()
{
  delete [] id;
  delete [] style;
}

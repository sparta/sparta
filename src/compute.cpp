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
#include "ctype.h"
#include "compute.h"
#include "error.h"

using namespace DSMC_NS;

/* ---------------------------------------------------------------------- */

Compute::Compute(DSMC *dsmc, int narg, char **arg) : Pointers(dsmc)
{
  if (narg < 2) error->all(FLERR,"Illegal compute command");

  // compute ID and style
  // ID must be all alphanumeric chars or underscores

  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  for (int i = 0; i < n-1; i++)
    if (!isalnum(id[i]) && id[i] != '_')
      error->all(FLERR,
		 "Compute ID must be alphanumeric or underscore characters");

  n = strlen(arg[1]) + 1;
  style = new char[n];
  strcpy(style,arg[1]);

  // set child class defaults

  scalar_flag = vector_flag = array_flag = 0;
  per_particle_flag = per_cell_flag = 0;

  invoked_scalar = invoked_vector = invoked_array = -1;
  invoked_per_particle = invoked_per_cell = -1;
}

/* ---------------------------------------------------------------------- */

Compute::~Compute()
{
  delete [] id;
  delete [] style;
}

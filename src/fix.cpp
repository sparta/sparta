/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "string.h"
#include "ctype.h"
#include "fix.h"
#include "memory.h"
#include "error.h"
#include "sparta_masks.h"

using namespace SPARTA_NS;

#define DELTA 4

/* ---------------------------------------------------------------------- */

Fix::Fix(SPARTA *sparta, int, char **arg) : Pointers(sparta)
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

  // set child class defaults

  time_depend = 0;
  gridmigrate = 0;
  flag_update_custom = flag_gas_react = flag_surf_react = 0;

  scalar_flag = vector_flag = array_flag = 0;
  per_particle_flag = per_grid_flag = per_surf_flag = 0;
  per_particle_field = per_grid_field = 0;

  // mask settings - same as in modify.cpp

  START_OF_STEP = 1;
  END_OF_STEP = 2;

  kokkos_flag = 0;
  copy = copymode = 0;

  execution_space = Host;
  datamask_read = ALL_MASK;
  datamask_modify = ALL_MASK;
}

/* ---------------------------------------------------------------------- */

Fix::~Fix()
{
  if (copy || copymode) return;

  delete [] id;
  delete [] style;
}

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

#include "mpi.h"
#include "ctype.h"
#include "string.h"
#include "surf_collide.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

SurfCollide::SurfCollide(SPARTA *sparta, int, char **arg) :
  Pointers(sparta)
{
  // ID and style
  // ID must be all alphanumeric chars or underscores

  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  for (int i = 0; i < n-1; i++)
    if (!isalnum(id[i]) && id[i] != '_')
      error->all(FLERR,"Surf_collide ID must be alphanumeric or "
                 "underscore characters");

  n = strlen(arg[1]) + 1;
  style = new char[n];
  strcpy(style,arg[1]);

  dynamicflag = 0;
  allowreact = 1;
  transparent = 0;
  vector_flag = 1;
  size_vector = 2;

  nsingle = ntotal = 0;

  kokkosable = copy = copymode = 0;
}

/* ---------------------------------------------------------------------- */

SurfCollide::~SurfCollide()
{
  if (copy) return;

  delete [] id;
  delete [] style;
}

/* ---------------------------------------------------------------------- */

void SurfCollide::init()
{
  nsingle = ntotal = 0;
}

/* ---------------------------------------------------------------------- */

void SurfCollide::tally_reset()
{
  nsingle = 0;
}

/* ---------------------------------------------------------------------- */

void SurfCollide::tally_update()
{
  ntotal += nsingle;
}

/* ---------------------------------------------------------------------- */

double SurfCollide::compute_vector(int i)
{
  one[0] = nsingle;
  one[1] = ntotal;
  MPI_Allreduce(one,all,2,MPI_DOUBLE,MPI_SUM,world);

  return all[i];
}

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

#include "math.h"
#include "ctype.h"
#include "stdlib.h"
#include "string.h"
#include "surf_react.h"
#include "input.h"
#include "update.h"
#include "comm.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

SurfReact::SurfReact(SPARTA *sparta, int, char **arg) :
  Pointers(sparta)
{
  // ID and style
  // ID must be all alphanumeric chars or underscores

  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  for (int i = 0; i < n-1; i++)
    if (!isalnum(id[i]) && id[i] != '_')
      error->all(FLERR,"Surf_react ID must be alphanumeric or "
                 "underscore characters");

  n = strlen(arg[1]) + 1;
  style = new char[n];
  strcpy(style,arg[1]);

  vector_flag = 1;
  size_vector = 2;

  // tallies

  nsingle = ntotal = 0;
  tally_two_flag = tally_single_flag = tally_total_flag = 0;

  kokkosable = copy = copymode = 0;
}

/* ---------------------------------------------------------------------- */

SurfReact::~SurfReact()
{
  if (copy) return;

  delete [] id;
  delete [] style;

  delete [] tally_single;
  delete [] tally_total;
  delete [] tally_single_all;
  delete [] tally_total_all;
}

/* ---------------------------------------------------------------------- */

void SurfReact::init()
{
  nsingle = ntotal = 0;
  for (int i = 0; i < nlist; i++)
    tally_single[i] = tally_total[i] = 0;
}

/* ---------------------------------------------------------------------- */

void SurfReact::tally_reset()
{
  nsingle = 0;
  for (int i = 0; i < nlist; i++) tally_single[i] = 0;
  tally_two_flag = tally_single_flag = tally_total_flag = 0;
}

/* ---------------------------------------------------------------------- */

void SurfReact::tally_update()
{
  ntotal += nsingle;
  for (int i = 0; i < nlist; i++) tally_total[i] += tally_single[i];
}

/* ---------------------------------------------------------------------- */

double SurfReact::compute_vector(int i)
{
  if (i < 2) {
    if (!tally_two_flag) {
      tally_two_flag = 1;
      one[0] = nsingle;
      one[1] = ntotal;
      MPI_Allreduce(one,all,2,MPI_DOUBLE,MPI_SUM,world);
    }
    return all[i];
  }

  if (i < 2+nlist) {
    if (!tally_single_flag) {
      tally_single_flag = 1;
      MPI_Allreduce(tally_single,tally_single_all,nlist,MPI_INT,MPI_SUM,world);
    }
    return 1.0*tally_single_all[i-2];
  }

  if (!tally_total_flag) {
    tally_total_flag = 1;
    MPI_Allreduce(tally_total,tally_total_all,nlist,MPI_INT,MPI_SUM,world);
  }
  return 1.0*tally_total_all[i-nlist-2];
}

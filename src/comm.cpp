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

#include "comm.h"
#include "error.h"

using namespace DSMC_NS;

/* ---------------------------------------------------------------------- */

Comm::Comm(DSMC *dsmc) : Pointers(dsmc)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  if (nprocs > 1) error->all(FLERR,"Cannot yet run in parallel");
}

/* ---------------------------------------------------------------------- */

Comm::~Comm()
{
}

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

#include "update.h"
#include "particle.h"
#include "comm.h"

using namespace DSMC_NS;

/* ---------------------------------------------------------------------- */

Update::Update(DSMC *dsmc) : Pointers(dsmc)
{
  dt = 0.1;
}

/* ---------------------------------------------------------------------- */

Update::~Update() {}

/* ---------------------------------------------------------------------- */

void Update::run(int nsteps)
{
  if (comm->me == 0) {
    if (screen)
      fprintf(screen,"Performing run ...\n");
    if (logfile)
      fprintf(logfile,"Performing run ...\n");
  }

  // loop over timesteps

  for (int i = 0; i < nsteps; i++) {

    ntimestep++;

    // move and communicate particles

    particle->move();
  }
}

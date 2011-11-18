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
#include "create_particles.h"
#include "domain.h"
#include "error.h"

using namespace DSMC_NS;

/* ---------------------------------------------------------------------- */

CreateParticles::CreateParticles(DSMC *dsmc) : Pointers(dsmc) {}

/* ---------------------------------------------------------------------- */

void CreateParticles::command(int narg, char **arg)
{
  if (narg != 6) error->all(FLERR,"Illegal create_particles command");

  if (!domain->box_exist) 
    error->all(FLERR,
	       "Cannot create_particles before simulation box is defined");;

}

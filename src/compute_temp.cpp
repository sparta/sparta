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
#include "string.h"
#include "compute_temp.h"
#include "update.h"
#include "particle.h"
#include "domain.h"
#include "error.h"

using namespace DSMC_NS;

/* ---------------------------------------------------------------------- */

ComputeTemp::ComputeTemp(DSMC *dsmc, int narg, char **arg) : 
  Compute(dsmc, narg, arg)
{
  if (narg != 2) error->all(FLERR,"Illegal compute temp command");

  scalar_flag = 1;
}

/* ---------------------------------------------------------------------- */

double ComputeTemp::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  Particle::Species *species = particle->species;
  Particle::OnePart *particles = particle->particles;
  int nlocal = particle->nlocal;

  double *v;
  double t = 0.0;

  for (int i = 0; i < nlocal; i++) {
    v = particles[i].v;
    t += (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) * 
      species[particles[i].ispecies].mass;
  }

  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);

  bigint n = particle->nlocal;
  MPI_Allreduce(&n,&particle->nglobal,1,MPI_DSMC_BIGINT,MPI_SUM,world);
  if (particle->nglobal == 0) return 0.0;

  double factor = update->mvv2e / 
    (domain->dimension * particle->nglobal * update->kboltz);
  scalar *= factor;
  return scalar;
}

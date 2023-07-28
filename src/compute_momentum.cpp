/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "string.h"
#include "compute_momentum.h"
#include "update.h"
#include "particle.h"
#include "domain.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

ComputeMomentum::ComputeMomentum(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg != 2) error->all(FLERR,"Illegal compute temp command");

  vector_flag = 1;
  size_vector = 3;
  vector = new double[size_vector];
}

/* ---------------------------------------------------------------------- */

void ComputeMomentum::compute_vector()
{
  invoked_scalar = update->ntimestep;

  Particle::Species *species = particle->species;
  Particle::OnePart *particles = particle->particles;
  int nlocal = particle->nlocal;

  double *v;
  double mom[3] = {0.0, 0.0, 0.0};

  for (int i = 0; i < nlocal; i++) {
    v = particles[i].v;
      for(int j = 0; j < 3; ++j){
          mom[j] += species[particles[i].ispecies].mass* v[j];
      }

  }
    MPI_Allreduce(&mom, vector, 3, MPI_DOUBLE, MPI_SUM, world);
}

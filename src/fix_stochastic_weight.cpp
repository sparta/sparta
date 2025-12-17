/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_stochastic_weight.h"
#include "update.h"
#include "particle.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{INT,DOUBLE};                      // several files

/* ---------------------------------------------------------------------- */

FixStochasticWeight::FixStochasticWeight(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg != 2) error->all(FLERR,"Illegal fix stochastic_weight command");

  flag_update_custom = 1;

  // check if custom attribute already exists, due to restart file
  // else create per-particle vector

  stochastic_wt_index = particle->find_custom((char *) "stochastic_wt");

  if (stochastic_wt_index < 0)
    stochastic_wt_index = particle->add_custom((char *) "stochastic_wt",DOUBLE,0);
}

/* ---------------------------------------------------------------------- */

FixStochasticWeight::~FixStochasticWeight()
{
  if (copy || copymode) return;

  particle->remove_custom(stochastic_wt_index);
}

/* ---------------------------------------------------------------------- */

int FixStochasticWeight::setmask()
{
  int mask = 0;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixStochasticWeight::init()
{
  // initialize custom weight attribute to 1
  // subsequent updates will store weight/fnum

  if (stochastic_wt_index >= 0) {
    double *stochastic_weights = particle->edvec[particle->ewhich[stochastic_wt_index]];
    for (int i = 0; i < particle->nlocal; i++)
      stochastic_weights[i] = 1.0;
  }
}

/* ----------------------------------------------------------------------
   called when a particle with index is created
   or when temperature dependent properties need to be updated
   initialize custom weight attribute with particle weight
------------------------------------------------------------------------- */

void FixStochasticWeight::update_custom(int index, double temp_thermal,
                            double, double,
                            double *vstream)
{
  double *stochastic_weights = particle->edvec[particle->ewhich[stochastic_wt_index]];
  stochastic_weights[index] = particle->particles[index].weight / update->fnum;
}



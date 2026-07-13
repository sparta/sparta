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

  // reuse the custom attribute if it already exists (e.g. restored from a
  // restart file, where the stored weights must be preserved); otherwise
  // create it and initialize any pre-existing particles to a weight of 1.0
  // (relative to fnum).  newly created particles are set by update_custom().

  stochastic_wt_index = particle->find_custom((char *) "stochastic_wt");

  if (stochastic_wt_index < 0) {
    stochastic_wt_index = particle->add_custom((char *) "stochastic_wt",DOUBLE,0);
    double *stochastic_weights = particle->edvec[particle->ewhich[stochastic_wt_index]];
    for (int i = 0; i < particle->nlocal; i++)
      stochastic_weights[i] = 1.0;
  }
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
  // intentionally empty: stochastic weights persist across runs and restarts.
  // they are initialized once at fix creation and for each new particle in
  // update_custom(); resetting them here would discard accumulated weights.
}

/* ----------------------------------------------------------------------
   called when a particle with index is created (emit, create_particles)
   set its stochastic weight to 1.0 (relative to fnum)
------------------------------------------------------------------------- */

void FixStochasticWeight::update_custom(int index, double temp_thermal,
                            double, double,
                            double *vstream)
{
  // a standard particle represents fnum real molecules, so its weight
  // relative to fnum is 1.0 (grid-based particle->weight is mutually
  // exclusive with stochastic weighting and is not used here)

  double *stochastic_weights = particle->edvec[particle->ewhich[stochastic_wt_index]];
  stochastic_weights[index] = 1.0;
}



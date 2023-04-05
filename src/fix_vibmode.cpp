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
#include "fix_vibmode.h"
#include "update.h"
#include "particle.h"
#include "collide.h"
#include "comm.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "math_const.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{INT,DOUBLE};                      // several files
enum{NONE,DISCRETE,SMOOTH};            // several files

/* ---------------------------------------------------------------------- */

FixVibmode::FixVibmode(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg != 2) error->all(FLERR,"Illegal fix vibmode command");

  flag_update_custom = 1;

  // random = RNG for vibrational mode initialization

  random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,comm->me,100);

  // create per-particle array

  if (!collide)
    error->all(FLERR,"Cannot use fix vibmode without "
               "collide style defined");

  if (collide->vibstyle != DISCRETE)
    error->all(FLERR,"Cannot use fix vibmode without "
               "collide_modify vibrate discrete");

  maxmode = particle->maxvibmode;
  if (maxmode <= 1)
    error->all(FLERR,"No multiple vibrational modes in fix vibmode "
               "for any species");

  vibmodeindex = particle->add_custom((char *) "vibmode",INT,maxmode);
}

/* ---------------------------------------------------------------------- */

FixVibmode::~FixVibmode()
{
  if (copy) return;

  delete random;
  particle->remove_custom(vibmodeindex);
}

/* ---------------------------------------------------------------------- */

int FixVibmode::setmask()
{
  int mask = 0;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixVibmode::init()
{
  if (collide->vibstyle != DISCRETE)
    error->all(FLERR,"Cannot use fix vibmode without "
               "collide_modify vibrate discrete");
  if (maxmode != particle->maxvibmode)
    error->all(FLERR,"Number of maximum vibrational modes has changed since "
               "fix vibmode was specified");
}

/* ----------------------------------------------------------------------
   called when a particle with index is created
    or when temperature dependent properties need to be updated
   populate all vibrational modes and set evib = sum of mode energies
------------------------------------------------------------------------- */

void FixVibmode::update_custom(int index, double temp_thermal,
                               double temp_rot, double temp_vib,
                               double *vstream)
{
  int **vibmode = particle->eiarray[particle->ewhich[vibmodeindex]];

  int isp = particle->particles[index].ispecies;
  int nmode = particle->species[isp].nvibmode;

  // no modes, just return

  if (nmode == 0) return;

  // single mode, evib already set by Particle::evib()
  // just convert evib back to mode level

  if (nmode == 1) {
    vibmode[index][0] = static_cast<int>
      (particle->particles[index].evib / update->boltz /
       particle->species[isp].vibtemp[0]);
    return;
  }

  // loop over modes and populate each
  // accumlate new total evib

  int ivib;
  double evib = 0.0;

  for (int imode = 0; imode < nmode; imode++) {
    ivib = static_cast<int> (-log(random->uniform()) * temp_vib /
                             particle->species[isp].vibtemp[imode]);
    vibmode[index][imode] = ivib;
    evib += ivib * update->boltz * particle->species[isp].vibtemp[imode];
  }

  particle->particles[index].evib = evib;
}

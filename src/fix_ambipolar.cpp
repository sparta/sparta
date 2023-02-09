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
#include "fix_ambipolar.h"
#include "update.h"
#include "particle.h"
#include "comm.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{INT,DOUBLE};                      // several files

/* ---------------------------------------------------------------------- */

FixAmbipolar::FixAmbipolar(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix ambipolar command");

  flag_update_custom = 1;
  flag_surf_react = 1;

  especies = particle->find_species(arg[2]);
  if (especies < 0) error->all(FLERR,"Fix ambipolar species does not exist");
  if (particle->species[especies].charge >= 0.0)
    error->all(FLERR,"Fix ambipolar electron species has charge >= 0.0");

  // ions[i] = 0 if species I is not a listed ion, 1 if it is

  maxion = particle->nspecies;
  memory->create(ions,maxion,"ambipolar:ions");
  for (int i = 0; i < maxion; i++) ions[i] = 0;

  int ispecies;
  for (int iarg = 3; iarg < narg; iarg++) {
    ispecies = particle->find_species(arg[iarg]);
    if (ispecies < 0) error->all(FLERR,"Fix ambipolar species does not exist");
    if (particle->species[ispecies].charge <= 0.0)
    error->all(FLERR,"Fix ambipolar ion species has charge <= 0.0");
    ions[ispecies] = 1;
  }

  // random = RNG for electron velocity creation

  random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,comm->me,100);

  // create per-particle vector and array

  ionindex = particle->add_custom((char *) "ionambi",INT,0);
  velindex = particle->add_custom((char *) "velambi",DOUBLE,3);
}

/* ---------------------------------------------------------------------- */

FixAmbipolar::~FixAmbipolar()
{
  if (copy || copymode) return;

  memory->destroy(ions);
  delete random;
  particle->remove_custom(ionindex);
  particle->remove_custom(velindex);
}

/* ---------------------------------------------------------------------- */

int FixAmbipolar::setmask()
{
  int mask = 0;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAmbipolar::init()
{
  if (maxion != particle->nspecies)
    error->all(FLERR,"Number of particle species has changed since "
               "fix ambipolar was specified");
}

/* ----------------------------------------------------------------------
   called when a particle with index is created
    or when temperature dependent properties need to be updated
   creation used temp_thermal and vstream to set particle velocity
   if an ion, set ionambi and velambi for particle
------------------------------------------------------------------------- */

void FixAmbipolar::update_custom(int index, double temp_thermal,
                                double, double,
                                double *vstream)
{
  int *ionambi = particle->eivec[particle->ewhich[ionindex]];
  double **velambi = particle->edarray[particle->ewhich[velindex]];

  // if species is not ambipolar ion, set ionambi off and return

  int ispecies = particle->particles[index].ispecies;

  if (ions[ispecies] == 0) {
    ionambi[index] = 0;
    return;
  }

  // set velocity of electron
  // based on electron mass, thermal temperature, and streaming velocity

  ionambi[index] = 1;

  double vscale = sqrt(2.0 * update->boltz * temp_thermal /
                       particle->species[especies].mass);

  double vn = vscale * sqrt(-log(random->uniform()));
  double vr = vscale * sqrt(-log(random->uniform()));
  double theta1 = MY_2PI * random->uniform();
  double theta2 = MY_2PI * random->uniform();

  velambi[index][0] = vstream[0] + vn*cos(theta1);
  velambi[index][1] = vstream[1] + vr*cos(theta2);
  velambi[index][2] = vstream[2] + vr*sin(theta2);
}

/* ----------------------------------------------------------------------
   called when a surface reaction occurs
   iorig = particle I before reaction
   I,J = indices of two particles after reaction
         either can be -1, meaning particle does not exist
------------------------------------------------------------------------- */

void FixAmbipolar::surf_react(Particle::OnePart *iorig, int &i, int &j)
{
  int ispecies = iorig->ispecies;

  // recombination reaction, just return

  if (i < 0) return;

  // exchange reaction
  // if ion -> non-ion, unset ionambi flag

  if (j < 0) {
    if (ions[ispecies] == 0) return;
    Particle::OnePart *particles = particle->particles;
    if (ions[particles[i].ispecies] == 1) return;
    int *ionambi = particle->eivec[particle->ewhich[ionindex]];
    ionambi[i] = 0;
  }

  // dissociation reaction
  // if non-ion -> ion + electron, create an ambipolar ion
  // use global temp_thermal and vstream for electron creation
  // set j = -1 to delete electron that was just created by caller

  else {
    if (ions[ispecies] == 1) return;
    Particle::OnePart *particles = particle->particles;
    if (ions[particles[i].ispecies] == 0) return;
    if (particles[j].ispecies != especies) return;
    update_custom(i,update->temp_thermal,update->temp_thermal,
                 update->temp_thermal,update->vstream);
    j = -1;
  }
}

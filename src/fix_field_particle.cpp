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

#include "stdlib.h"
#include "string.h"
#include "fix_field_particle.h"
#include "particle.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

FixFieldParticle::FixFieldParticle(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Illegal fix field/particle command");

  int ncols = 0;

  if (strcmp(arg[2],"NULL") == 0) axstr = NULL;
  else {
    int n = strlen(arg[2]) + 1;
    axstr = new char[n];
    strcpy(axstr,arg[2]);
    ncols++;
  }
  if (strcmp(arg[3],"NULL") == 0) aystr = NULL;
  else {
    int n = strlen(arg[3]) + 1;
    aystr = new char[n];
    strcpy(aystr,arg[3]);
    ncols++;
  }
  if (strcmp(arg[4],"NULL") == 0) azstr = NULL;
  else {
    int n = strlen(arg[4]) + 1;
    azstr = new char[n];
    strcpy(azstr,arg[4]);
    ncols++;
  }

  // fix settings

  per_particle_flag = 1;
  size_per_particle_cols = ncols;
  per_particle_freq = 1;
  per_particle_field = 1;

  field_active[0] = field_active[1] = field_active[2] = 0;
  if (axstr) field_active[0] = 1;
  if (aystr) field_active[1] = 1;
  if (azstr) field_active[2] = 1;

  // per-particle memory initialization

  maxparticle = 0;
  array_particle = NULL;
}

/* ---------------------------------------------------------------------- */

FixFieldParticle::~FixFieldParticle()
{
  delete [] axstr;
  delete [] aystr;
  delete [] azstr;

  memory->destroy(array_particle);
}

/* ---------------------------------------------------------------------- */

int FixFieldParticle::setmask()
{
  int mask = 0;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixFieldParticle::init()
{
  // check if all variables exist and are particle-style vars

  if (axstr) {
    axvar = input->variable->find(axstr);
    if (axvar < 0)
      error->all(FLERR,"Variable name for fix field/particle does not exist");
    if (!input->variable->particle_style(axvar))
      error->all(FLERR,"Variable for fix field/particle is invalid style");
  }
  if (aystr) {
    ayvar = input->variable->find(aystr);
    if (ayvar < 0)
      error->all(FLERR,"Variable name for fix field/particle does not exist");
    if (!input->variable->particle_style(ayvar))
      error->all(FLERR,"Variable for fix field/particle is invalid style");
  }
  if (azstr) {
    azvar = input->variable->find(azstr);
    if (azvar < 0)
      error->all(FLERR,"Variable name for fix field/particle does not exist");
    if (!input->variable->particle_style(azvar))
      error->all(FLERR,"Variable for fix field/particle is invalid style");
  }

  // set initial particle values to zero in case dump is performed at step 0

  if (particle->nlocal > maxparticle) {
    maxparticle = particle->maxlocal;
    memory->destroy(array_particle);
    memory->create(array_particle,maxparticle,size_per_particle_cols,
                   "array_particle");
  }

  bigint nbytes = (bigint) particle->nlocal * size_per_particle_cols;
  memset(&array_particle[0][0],0,nbytes*sizeof(double));
}

/* ---------------------------------------------------------------------- */

void FixFieldParticle::compute_field()
{
  // reallocate array_particle if necessary

  if (particle->nlocal > maxparticle) {
    maxparticle = particle->maxlocal;
    memory->destroy(array_particle);
    memory->create(array_particle,maxparticle,size_per_particle_cols,
                   "array_particle");
  }

  // evaluate each particle-style variable
  // results are put into strided array_particle

  int stride = size_per_particle_cols;
  int icol = 0;

  if (axstr) {
    input->variable->compute_particle(axvar,&array_particle[0][icol],stride,0);
    icol++;
  }

  if (aystr) {
    input->variable->compute_particle(ayvar,&array_particle[0][icol],stride,0);
    icol++;
  }

  if (azstr) {
    input->variable->compute_particle(azvar,&array_particle[0][icol],stride,0);
    icol++;
  }
}

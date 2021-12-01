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
  if (narg != 8) error->all(FLERR,"Illegal fix field/particle command");

  int ncols = 0;

  if (strcmp(arg[2],"NULL") == 0) xstr = NULL;
  else {
    int n = strlen(arg[2]) + 1;
    xstr = new char[n];
    strcpy(xstr,arg[2]);
    ncols++;
  }
  if (strcmp(arg[3],"NULL") == 0) ystr = NULL;
  else {
    int n = strlen(arg[3]) + 1;
    ystr = new char[n];
    strcpy(ystr,arg[3]);
    ncols++;
  }
  if (strcmp(arg[4],"NULL") == 0) zstr = NULL;
  else {
    int n = strlen(arg[4]) + 1;
    zstr = new char[n];
    strcpy(zstr,arg[4]);
    ncols++;
  }

  if (strcmp(arg[5],"NULL") == 0) vxstr = NULL;
  else {
    int n = strlen(arg[5]) + 1;
    vxstr = new char[n];
    strcpy(vxstr,arg[5]);
    ncols++;
  }
  if (strcmp(arg[6],"NULL") == 0) vystr = NULL;
  else {
    int n = strlen(arg[6]) + 1;
    vystr = new char[n];
    strcpy(vystr,arg[6]);
    ncols++;
  }
  if (strcmp(arg[7],"NULL") == 0) vzstr = NULL;
  else {
    int n = strlen(arg[7]) + 1;
    vzstr = new char[n];
    strcpy(vzstr,arg[7]);
    ncols++;
  }

  // fix settings

  per_particle_flag = 1;
  size_per_particle_cols = ncols;
  per_particle_freq = 1;
  per_particle_field = 1;

  field_active[0] = field_active[1] = field_active[2] = 0;
  field_active[3] = field_active[4] = field_active[5] = 0;
  if (xstr) field_active[0] = 1;
  if (ystr) field_active[1] = 1;
  if (zstr) field_active[2] = 1;
  if (vxstr) field_active[3] = 1;
  if (vystr) field_active[4] = 1;
  if (vzstr) field_active[5] = 1;

  // per-particle memory initialization

  maxparticle = 0;
  array_particle = NULL;
}

/* ---------------------------------------------------------------------- */

FixFieldParticle::~FixFieldParticle()
{  
  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  delete [] vxstr;
  delete [] vystr;
  delete [] vzstr;
  
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

  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0)
      error->all(FLERR,"Variable name for fix field/particles does not exist");
    if (!input->variable->particle_style(xvar))
      error->all(FLERR,"Variable for fix field/particles is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (xvar < 0)
      error->all(FLERR,"Variable name for fix field/particles does not exist");
    if (!input->variable->particle_style(yvar))
      error->all(FLERR,"Variable for fix field/particles is invalid style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0)
      error->all(FLERR,"Variable name for fix field/particles does not exist");
    if (!input->variable->particle_style(zvar))
      error->all(FLERR,"Variable for fix field/particles is invalid style");
  }

  if (vxstr) {
    vxvar = input->variable->find(vxstr);
    if (vxvar < 0)
      error->all(FLERR,"Variable name for fix field/particles does not exist");
    if (!input->variable->particle_style(vxvar))
      error->all(FLERR,"Variable for fix field/particles is invalid style");
  }
  if (vystr) {
    vyvar = input->variable->find(vystr);
    if (vxvar < 0)
      error->all(FLERR,"Variable name for fix field/particles does not exist");
    if (!input->variable->particle_style(vyvar))
      error->all(FLERR,"Variable for fix field/particles is invalid style");
  }
  if (vzstr) {
    vzvar = input->variable->find(vzstr);
    if (vzvar < 0)
      error->all(FLERR,"Variable name for fix field/particles does not exist");
    if (!input->variable->particle_style(vzvar))
      error->all(FLERR,"Variable for fix field/particles is invalid style");
  }
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

  if (xstr) {
    input->variable->compute_particle(xvar,&array_particle[0][icol],stride,0);
    icol++;
  }

  if (ystr) {
    input->variable->compute_particle(yvar,&array_particle[0][icol],stride,0);
    icol++;
  }

  if (zstr) {
    input->variable->compute_particle(zvar,&array_particle[0][icol],stride,0);
    icol++;
  }

  if (vxstr) {
    input->variable->compute_particle(vxvar,&array_particle[0][icol],stride,0);
    icol++;
  }

  if (vystr) {
    input->variable->compute_particle(vyvar,&array_particle[0][icol],stride,0);
    icol++;
  }

  if (vzstr) {
    input->variable->compute_particle(vzvar,&array_particle[0][icol],stride,0);
    icol++;
  }
}

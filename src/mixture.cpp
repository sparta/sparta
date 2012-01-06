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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "ctype.h"
#include "mixture.h"
#include "update.h"
#include "particle.h"
#include "memory.h"
#include "error.h"

using namespace DSMC_NS;

#define DELTA 8

/* ---------------------------------------------------------------------- */

Mixture::Mixture(DSMC *dsmc, char *userid) : Pointers(dsmc)
{
  // mixture ID must be all alphanumeric chars or underscores

  int n = strlen(userid) + 1;
  id = new char[n];
  strcpy(id,userid);

  for (int i = 0; i < n-1; i++)
    if (!isalnum(id[i]) && id[i] != '_')
      error->all(FLERR,
		 "Mixture ID must be alphanumeric or underscore characters");

  // initialize mixture values

  species = NULL;

  nspecies = maxspecies = 0;

  fraction = NULL;
  fraction_user = NULL;
  fraction_flag = NULL;

  nrho_flag = 0;
  vstream_flag = 0;
  temp_thermal_flag = 0;

  cummulative = NULL;
  vscale = NULL;
  active = NULL;

  allocate(DELTA);
}

/* ---------------------------------------------------------------------- */

Mixture::~Mixture()
{
  delete [] id;

  memory->destroy(species);
  memory->destroy(fraction);
  memory->destroy(fraction_user);
  memory->destroy(fraction_flag);

  memory->destroy(cummulative);
  memory->destroy(vscale);
  memory->destroy(active);
}

/* ----------------------------------------------------------------------
   set attributes of each mixture species to either user settings or defaults
   defaults = background gas properties stored in Update class
------------------------------------------------------------------------- */

void Mixture::init()
{
  if (nrho_flag) nrho = nrho_user;
  else nrho = update->nrho;
  if (vstream_flag) {
    vstream[0] = vstream_user[0];
    vstream[1] = vstream_user[1];
    vstream[2] = vstream_user[2];
  } else {
    vstream[0] = update->vstream[0];
    vstream[1] = update->vstream[1];
    vstream[2] = update->vstream[2];
  }
  if (temp_thermal_flag) temp_thermal = temp_thermal_user;
  else temp_thermal = update->temp_thermal;

  // frac_explicit = sum for species with explicity set fractions
  // frac_implicit = number of unset species

  double fraction_explicit = 0.0;
  int fraction_implicit = 0;
  for (int i = 0; i < nspecies; i++) {
    if (fraction_flag[i]) fraction_explicit += fraction_user[i];
    else fraction_implicit++;
  }

  if (fraction_explicit > 1.0) {
    char str[128];
    sprintf(str,"Mixture %s fractions exceed 1.0",id);
    error->all(FLERR,str);
  }

  // fraction for each unset species = equal portion of unset remainder
  // cummulative = cummulative fraction across species

  for (int i = 0; i < nspecies; i++) {
    if (fraction_flag[i]) fraction[i] = fraction_user[i];
    else fraction[i] = (1.0-fraction_explicit) / fraction_implicit;
    if (i) cummulative[i] = cummulative[i-1] + fraction[i];
    else cummulative[i] = fraction[i];
  }  
  cummulative[nspecies-1] = 1.0;

  // vscale = factor to scale Gaussian unit variance by
  //          to get thermal distribution of velocities
  // per-species value since includes species mass

  for (int i = 0; i < nspecies; i++) {
    int index = species[i];
    vscale[i] = sqrt(update->boltz * temp_thermal /
		     particle->species[index].mass);
  }
}

/* ---------------------------------------------------------------------- */

void Mixture::add_species(int narg, char **arg)
{
  int i,j,index,exist;

  // nnew = # of species not in mixture list

  int nnew = 0;
  for (i = 0; i < narg; i++) {
    index = particle->find_species(arg[i]);
    if (index < 0) error->all(FLERR,"Mixture species is not defined");
    exist = 0;
    for (j = 0; j < nspecies; j++)
      if (species[j] == index) exist = 1;
    if (!exist) nnew++;
  }

  // grow per-species arrays if necessary

  if (nspecies + nnew > maxspecies) allocate(nspecies+nnew);

  // add new species to mixture list
  // set their user flags to 0

  for (i = 0; i < narg; i++) {
    index = particle->find_species(arg[i]);
    exist = 0;
    for (j = 0; j < nspecies; j++)
      if (species[j] == index) exist = 1;
    if (!exist) species[nspecies++] = index;
  }

  // active[i] = 1 if mixture species is active in this input script command

  for (i = 0; i < nspecies; i++) active[i] = 0;
  for (i = 0; i < narg; i++) {
    index = particle->find_species(arg[i]);
    for (j = 0; j < nspecies; j++)
      if (species[j] == index) active[j] = 1;
  }
}

/* ---------------------------------------------------------------------- */

void Mixture::params(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"frac") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal mixture command");
      double value = atof(arg[iarg+1]);
      if (value < 0.0 || value > 1.0) 
	error->all(FLERR,"Illegal mixture command");
      for (int i = 0; i < nspecies; i++)
	if (active[i]) {
	  fraction_flag[i] = 1;
	  fraction_user[i] = value;
	}
      iarg += 2;
    } else if (strcmp(arg[iarg],"nrho") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal mixture command");
      nrho_user = atof(arg[iarg+1]);
      if (nrho_user <= 0.0) error->all(FLERR,"Illegal mixture command");
      nrho_flag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"vstream") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal mixture command");
      vstream_user[0] = atof(arg[iarg+1]);
      vstream_user[1] = atof(arg[iarg+2]);
      vstream_user[2] = atof(arg[iarg+3]);
      vstream_flag = 1;
      iarg += 4;
    } else if (strcmp(arg[iarg],"temp") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal mixture command");
      temp_thermal_user = atof(arg[iarg+1]);
      if (temp_thermal_user <= 0.0) error->all(FLERR,"Illegal mixture command");
      temp_thermal_flag = 1;
      iarg += 2;
    } else error->all(FLERR,"Illegal mixture command");
  }
}

/* ----------------------------------------------------------------------
   (re)allocate all mixture arrays to hold N species
   set new user flags to 0
------------------------------------------------------------------------- */

void Mixture::allocate(int n)
{
  int old = maxspecies;
  while (n > maxspecies) maxspecies += DELTA;

  memory->grow(species,maxspecies,"mixture:species");
  memory->grow(fraction,maxspecies,"mixture:species");
  memory->grow(fraction_user,maxspecies,"mixture:fraction_user");
  memory->grow(fraction_flag,maxspecies,"mixture:fraction_flag");

  memory->grow(cummulative,maxspecies,"mixture:cummulative");
  memory->grow(vscale,maxspecies,"mixture:vscale");
  memory->grow(active,maxspecies,"mixture:active");

  for (int i = old; i < maxspecies; i++) fraction_flag[i] = 0;
}

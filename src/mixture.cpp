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
#include "comm.h"
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

  nrho_flag = 0;
  vstream_flag = 0;
  temp_thermal_flag = 0;

  fraction = NULL;
  fraction_user = NULL;
  fraction_flag = NULL;
  cummulative = NULL;

  ngroups = maxgroup = 0;
  groups = NULL;
  mix2group = NULL;
  species2group = NULL;

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

  delete_groups();
  memory->sfree(groups);
  memory->destroy(mix2group);
  memory->destroy(species2group);

  memory->destroy(vscale);
  memory->destroy(active);
}

/* ----------------------------------------------------------------------
   set attributes of each mixture species to either user settings or defaults
   defaults = background gas properties stored in Update class
------------------------------------------------------------------------- */

void Mixture::init()
{
  // global attributes

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

  // sum = frac sum for species with explicity set fractions
  // nimplicit = number of unset species

  double sum = 0.0;
  int nimplicit = 0;
  for (int i = 0; i < nspecies; i++) {
    if (fraction_flag[i]) sum += fraction_user[i];
    else nimplicit++;
  }

  if (sum > 1.0) {
    char str[128];
    sprintf(str,"Mixture %s fractions exceed 1.0",id);
    error->all(FLERR,str);
  }

  // fraction for each unset species = equal portion of unset remainder
  // cummulative = cummulative fraction across species

  for (int i = 0; i < nspecies; i++) {
    if (fraction_flag[i]) fraction[i] = fraction_user[i];
    else fraction[i] = (1.0-sum) / nimplicit;
    if (i) cummulative[i] = cummulative[i-1] + fraction[i];
    else cummulative[i] = fraction[i];
  }  
  cummulative[nspecies-1] = 1.0;

  // warn if any group has no species assigned to it

  int count = 0;
  for (int i = 0; i < ngroups; i++) {
    int flag = 0;
    for (int j = 0; j < nspecies; j++)
      if (mix2group[j] == i) flag = 1;
    if (!flag) count++;
  }

  if (count && comm->me == 0) {
    char str[128];
    sprintf(str,"Mixture %s has %d groups with no species",id,count);
    error->warning(FLERR,str);
  }

  // setup species2group

  memory->destroy(species2group);
  memory->create(species2group,particle->nspecies,"mixture:species2group");
  for (int i = 0; i < particle->nspecies; i++) species2group[i] = -1;
  for (int i = 0; i < nspecies; i++) species2group[species[i]] = mix2group[i];

  // vscale = factor to scale Gaussian unit variance by
  //          to get thermal distribution of velocities
  // per-species value since includes species mass

  for (int i = 0; i < nspecies; i++) {
    int index = species[i];
    vscale[i] = sqrt(2.*update->boltz * temp_thermal /
		     particle->species[index].mass);
  }
}

/* ----------------------------------------------------------------------
   process list of species appearing in a mixture command
------------------------------------------------------------------------- */

void Mixture::add_species(int narg, char **arg)
{
  int i,j,index,exist;

  // nnew = # of species not already in mixture list

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

  for (i = 0; i < narg; i++) {
    index = particle->find_species(arg[i]);
    exist = 0;
    for (j = 0; j < nspecies; j++)
      if (species[j] == index) exist = 1;
    if (!exist) species[nspecies++] = index;
  }

  // active[i] = 1 if mixture species appears in this mixture command
  // else 0

  nactive = 0;
  for (i = 0; i < nspecies; i++) active[i] = 0;
  for (i = 0; i < narg; i++) {
    index = particle->find_species(arg[i]);
    for (j = 0; j < nspecies; j++)
      if (species[j] == index) {
	active[j] = 1;
	nactive++;
      }
  }
}

/* ----------------------------------------------------------------------
   process list of global or per-species params appearing in a mixture command
------------------------------------------------------------------------- */

void Mixture::params(int narg, char **arg)
{
  // for global attributes, set immediately
  // for per-species attributes, store flags

  int fracflag = 0;
  int groupflag = 0;

  double fracvalue;
  int grouparg;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"nrho") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal mixture command");
      nrho_flag = 1;
      nrho_user = atof(arg[iarg+1]);
      if (nrho_user <= 0.0) error->all(FLERR,"Illegal mixture command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"vstream") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal mixture command");
      vstream_flag = 1;
      vstream_user[0] = atof(arg[iarg+1]);
      vstream_user[1] = atof(arg[iarg+2]);
      vstream_user[2] = atof(arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"temp") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal mixture command");
      temp_thermal_flag = 1;
      temp_thermal_user = atof(arg[iarg+1]);
      if (temp_thermal_user <= 0.0) error->all(FLERR,"Illegal mixture command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"frac") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal mixture command");
      fracflag = 1;
      fracvalue = atof(arg[iarg+1]);
      if (fracvalue < 0.0 || fracvalue > 1.0) 
	error->all(FLERR,"Illegal mixture command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"group") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal mixture command");
      groupflag = 1;
      grouparg = iarg+1;
      iarg += 2;

    } else error->all(FLERR,"Illegal mixture command");
  }

  // apply per-species attributes

  if (fracflag) {
    for (int i = 0; i < nspecies; i++)
      if (active[i]) {
	fraction_flag[i] = 1;
	fraction_user[i] = fracvalue;
      }
  }

  // assign species to groups
  // if group-ID = all:
  ///  no listed species: delete groups, assign all species to group "all"
  //   yes listed species: assign listed species to group "all"
  // if group-ID = species 
  //   no listed species: delete groups, assign each species to own group
  //   yes listed species: assign each listed species to own group
  // if group-ID is user specified: assign listed species to that group
  // if no groupflag and species listed: assign listed species to group "all"

  if (groupflag) {
    if (strcmp(arg[grouparg],"all") == 0) {
      if (nactive == 0) {
	delete_groups();
	add_group("all");
	for (int i = 0; i < nspecies; i++)
	  mix2group[i] = 0;
      } else {
	int igroup = find_group("all");
	if (igroup < 0) {
	  add_group("all");
	  igroup = 0;
	}
	for (int i = 0; i < nspecies; i++)
	  if (active[i]) mix2group[i] = 0;
      }

    } else if (strcmp(arg[grouparg],"species") == 0) {
      if (nactive == 0) {
	delete_groups();
	for (int i = 0; i < nspecies; i++) {
	  add_group(particle->species[species[i]].id);
	  mix2group[i] = ngroups-1;
	}
      } else {
	for (int i = 0; i < nspecies; i++) {
	  if (!active[i]) continue;
	  int igroup = find_group(particle->species[species[i]].id);
	  if (igroup < 0) {
	    add_group(particle->species[species[i]].id);
	    igroup = ngroups-1;
	  }
	  mix2group[i] = igroup;
	}
      }

    } else {
      for (int i = 0; i < nspecies; i++) {
	if (!active[i]) continue;
	int igroup = find_group(arg[grouparg]);
	if (igroup < 0) {
	  add_group(arg[grouparg]);
	  igroup = ngroups-1;
	}
	mix2group[i] = igroup;
      }
    }

  } else if (nactive) {
    int igroup = find_group("all");
    if (igroup < 0) {
      add_group("all");
      igroup = 0;
    }
    for (int i = 0; i < nspecies; i++)
      if (active[i]) mix2group[i] = 0;
  }
}

/* ----------------------------------------------------------------------
   (re)allocate all mixture arrays to hold N species
   set new per-species flags to 0
------------------------------------------------------------------------- */

void Mixture::allocate(int n)
{
  int old = maxspecies;
  while (n > maxspecies) maxspecies += DELTA;

  memory->grow(species,maxspecies,"mixture:species");
  memory->grow(fraction,maxspecies,"mixture:species");
  memory->grow(fraction_flag,maxspecies,"mixture:fraction_flag");
  memory->grow(fraction_user,maxspecies,"mixture:fraction_user");
  memory->grow(cummulative,maxspecies,"mixture:cummulative");
  memory->grow(mix2group,maxspecies,"mixture:cummulative");
  memory->grow(vscale,maxspecies,"mixture:vscale");
  memory->grow(active,maxspecies,"mixture:active");

  for (int i = old; i < maxspecies; i++) fraction_flag[i] = 0;
}

/* ----------------------------------------------------------------------
   delete all groups
------------------------------------------------------------------------- */

void Mixture::delete_groups()
{
  for (int i = 0; i < ngroups; i++) delete [] groups[i];
  ngroups = 0;
}

/* ----------------------------------------------------------------------
   add a group with ID, assumed to not exist
------------------------------------------------------------------------- */

void Mixture::add_group(const char *idgroup)
{
  if (ngroups == maxgroup) {
    maxgroup += DELTA;
    groups = (char **) memory->srealloc(groups,maxgroup*sizeof(char *),
					"mixture:groups");
  }

  int n = strlen(idgroup) + 1;
  groups[ngroups] = new char[n];
  strcpy(groups[ngroups],idgroup);
  ngroups++;
}

/* ----------------------------------------------------------------------
   find group with ID and return its index
   return -1 if group does not exist
------------------------------------------------------------------------- */

int Mixture::find_group(const char *idgroup)
{
  for (int i = 0; i < ngroups; i++)
    if (strcmp(groups[i],idgroup) == 0) return i;
  return -1;
}

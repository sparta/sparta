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

using namespace SPARTA_NS;

#define DELTA 8

/* ---------------------------------------------------------------------- */

Mixture::Mixture(SPARTA *sparta, char *userid) : Pointers(sparta)
{
  // mixture ID must be all alphanumeric chars or underscores

  int n = strlen(userid) + 1;
  id = new char[n];
  strcpy(id,userid);

  for (int i = 0; i < n-1; i++)
    if (!isalnum(id[i]) && id[i] != '_')
      error->all(FLERR,
		 "Mixture ID must be alphanumeric or underscore characters");

  // special default mixtures

  all_default = species_default = 0;
  if (strcmp(id,"all") == 0) all_default = 1;
  if (strcmp(id,"species") == 0) species_default = 1;

  // initialize mixture values

  nspecies = maxspecies = 0;
  species = NULL;

  nrho_flag = 0;
  vstream_flag = 0;
  temp_thermal_flag = 0;

  fraction = NULL;
  fraction_user = NULL;
  fraction_flag = NULL;
  cummulative = NULL;

  ngroup = maxgroup = 0;
  groups = NULL;
  groupsize = NULL;
  mix2group = NULL;
  species2group = NULL;

  vscale = NULL;
  active = NULL;

  allocate();
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
  delete [] groupsize;
  memory->destroy(mix2group);
  memory->destroy(species2group);

  memory->destroy(vscale);
  memory->destroy(active);
}

/* ----------------------------------------------------------------------
   process args of a mixture command
   zero of more species IDs followed by optional args
------------------------------------------------------------------------- */

void Mixture::command(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal mixture command");

  // nsp = # of listed species before optional keywords
  // iarg = start of optional keywords

  int iarg;
  for (iarg = 1; iarg < narg; iarg++) {
    if (strcmp(arg[iarg],"nrho") == 0) break;
    if (strcmp(arg[iarg],"vstream") == 0) break;
    if (strcmp(arg[iarg],"temp") == 0) break;
    if (strcmp(arg[iarg],"frac") == 0) break;
    if (strcmp(arg[iarg],"group") == 0) break;
  }
  int nsp = iarg - 1;

  // add_species() for list of species
  // params() for reamining optional keywords

  add_species(nsp,&arg[1]);
  params(narg-iarg,&arg[iarg]);
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
  if (nspecies) cummulative[nspecies-1] = 1.0;

  // vscale = factor to scale Gaussian unit variance by
  //          to get thermal distribution of velocities
  // per-species value since includes species mass

  for (int i = 0; i < nspecies; i++) {
    int index = species[i];
    vscale[i] = sqrt(2.0 * update->boltz * temp_thermal /
		     particle->species[index].mass);
  }

  // setup species2group

  memory->destroy(species2group);
  memory->create(species2group,particle->nspecies,"mixture:species2group");
  for (int i = 0; i < particle->nspecies; i++) species2group[i] = -1;
  for (int i = 0; i < nspecies; i++) species2group[species[i]] = mix2group[i];

  // setup groupsize

  delete [] groupsize;
  groupsize = new int[ngroup];
  for (int i = 0; i < ngroup; i++) groupsize[i] = 0;
  for (int i = 0; i < nspecies; i++) groupsize[mix2group[i]]++;
}

/* ----------------------------------------------------------------------
   process list of species appearing in a mixture command
------------------------------------------------------------------------- */

void Mixture::add_species(int narg, char **arg)
{
  int i,j,index;

  // active[i] = 0 if species I not in list
  // active[i] = 1 if species I is in list and already exists in mixture
  // active[i] = 2 if species I is in list and added to mixture
  // nactive = # of species in list

  nactive = narg;
  for (i = 0; i < nspecies; i++) active[i] = 0;

  for (i = 0; i < narg; i++) {
    index = particle->find_species(arg[i]);
    if (index < 0) error->all(FLERR,"Mixture species is not defined");
    for (j = 0; j < nspecies; j++)
      if (species[j] == index) break;
    if (j < nspecies) active[j] = 1;
    else {
      if (all_default || species_default)
        error->all(FLERR,"Cannot add new species to mixture all or species");
      if (nspecies == maxspecies) allocate();
      active[nspecies] = 2;
      species[nspecies++] = index;
    }
  }
}

/* ----------------------------------------------------------------------
   add a species to the default mixture "all" or "species"
   set group assignment accordingly
------------------------------------------------------------------------- */

void Mixture::add_species_default(char *name)
{
  int index = particle->find_species(name);
  if (nspecies == maxspecies) allocate();
  species[nspecies] = index;

  if (all_default && ngroup == 0) add_group("all");
  if (species_default) add_group(name);
  mix2group[nspecies] = ngroup-1;

  nspecies++;
}

/* ----------------------------------------------------------------------
   process optional keywords appearing in a mixture command
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
      if (temp_thermal_user <= 0.0) 
        error->all(FLERR,"Illegal mixture command");
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
      int n = strlen(arg[grouparg]);
      for (int i = 0; i < n; i++)
        if (!isalnum(arg[grouparg][i]) && arg[grouparg][i] != '_')
          error->all(FLERR,"Mixture group ID must be "
                     "alphanumeric or underscore characters");
      if (all_default || species_default)
        error->all(FLERR,
                   "Cannot use group keyword with mixture all or species");
      iarg += 2;

    } else error->all(FLERR,"Illegal mixture command");
  }

  // assign per-species attributes

  if (fracflag) {
    for (int i = 0; i < nspecies; i++)
      if (active[i]) {
	fraction_flag[i] = 1;
	fraction_user[i] = fracvalue;
      }
  }

  // assign species to groups
  // end up with:
  //   every species assigned to exactly one group via mix2group
  //   no empty groups via shrink_groups()
  // if group-ID = SELF:
  //   no listed species: delete groups, assign each species to own group
  //   listed species: assign each listed species to own group
  // else if group-ID = user name:
  //   no listed species: delete groups, assign all species to group-ID
  //   listed species: assign listed species to group-ID
  // else if group keyword not specified:
  //   assign listed species that are new each to own group

  if (groupflag) {
    if (strcmp(arg[grouparg],"SELF") == 0) {
      if (nactive == 0) {
	delete_groups();
	for (int i = 0; i < nspecies; i++) {
	  add_group(particle->species[species[i]].id);
	  mix2group[i] = ngroup-1;
	}
      } else {
	for (int i = 0; i < nspecies; i++) {
	  if (!active[i]) continue;
	  int igroup = find_group(particle->species[species[i]].id);
	  if (igroup < 0) {
	    add_group(particle->species[species[i]].id);
	    igroup = ngroup-1;
	  }
	  mix2group[i] = igroup;
	}
      }

    } else {
      if (nactive == 0) {
	delete_groups();
        add_group(arg[grouparg]);
	for (int i = 0; i < nspecies; i++) mix2group[i] = ngroup-1;
      } else {
        int igroup = find_group(arg[grouparg]);
        if (igroup < 0) {
          add_group(arg[grouparg]);
          igroup = ngroup-1;
        }
        for (int i = 0; i < nspecies; i++)
          if (active[i]) mix2group[i] = igroup;
      }
    }

  } else if (nactive) {
    for (int i = 0; i < nspecies; i++) {
      if (active[i] != 2) continue;
      int igroup = find_group(particle->species[species[i]].id);
      if (igroup < 0) {
        add_group(particle->species[species[i]].id);
        igroup = ngroup-1;
      }
      mix2group[i] = igroup;
    }
  }

  shrink_groups();
}

/* ----------------------------------------------------------------------
   grow all mixture arrays
   set new per-species flags to 0
------------------------------------------------------------------------- */

void Mixture::allocate()
{
  int old = maxspecies;
  maxspecies += DELTA;

  memory->grow(species,maxspecies,"mixture:species");
  memory->grow(fraction,maxspecies,"mixture:species");
  memory->grow(fraction_flag,maxspecies,"mixture:fraction_flag");
  memory->grow(fraction_user,maxspecies,"mixture:fraction_user");
  memory->grow(cummulative,maxspecies,"mixture:cummulative");
  memory->grow(mix2group,maxspecies,"mixture:cummulative");
  memory->grow(vscale,maxspecies,"mixture:vscale");
  memory->grow(active,maxspecies,"mixture:active");

  for (int i = old; i < maxspecies; i++) {
    fraction_flag[i] = 0;
    fraction_user[i] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   delete all groups, leave maxgroup unchanged
------------------------------------------------------------------------- */

void Mixture::delete_groups()
{
  for (int i = 0; i < ngroup; i++) delete [] groups[i];
  ngroup = 0;
}

/* ----------------------------------------------------------------------
   remove empty groups, leave maxgroup unchanged
   compress group IDs, reset mix2group indices
------------------------------------------------------------------------- */

void Mixture::shrink_groups()
{
  int i,nsp;

  int igroup = 0;
  while (igroup < ngroup) {
    nsp = 0;
    for (i = 0; i < nspecies; i++)
      if (mix2group[i] == igroup) nsp++;
    if (nsp == 0) {
      delete [] groups[igroup];
      for (i = igroup; i < ngroup-1; i++)
        groups[i] = groups[i+1];
      for (i = 0; i < nspecies; i++)
        if (mix2group[i] > igroup) mix2group[i]--;
      ngroup--;
    } else igroup++;
  }
}

/* ----------------------------------------------------------------------
   add a group with ID, assumed to not exist
------------------------------------------------------------------------- */

void Mixture::add_group(const char *idgroup)
{
  if (ngroup == maxgroup) {
    maxgroup += DELTA;
    groups = (char **) memory->srealloc(groups,maxgroup*sizeof(char *),
					"mixture:groups");
  }

  int n = strlen(idgroup) + 1;
  groups[ngroup] = new char[n];
  strcpy(groups[ngroup],idgroup);
  ngroup++;
}

/* ----------------------------------------------------------------------
   find group with ID and return its index
   return -1 if group does not exist
------------------------------------------------------------------------- */

int Mixture::find_group(const char *idgroup)
{
  for (int i = 0; i < ngroup; i++)
    if (strcmp(groups[i],idgroup) == 0) return i;
  return -1;
}

/* ----------------------------------------------------------------------
   proc 0 writes mixture info to restart file
------------------------------------------------------------------------- */

void Mixture::write_restart(FILE *fp)
{
  fwrite(&nspecies,sizeof(int),1,fp);
  fwrite(species,sizeof(int),nspecies,fp);

  fwrite(&nrho_flag,sizeof(int),1,fp);
  if (nrho_flag) fwrite(&nrho_user,sizeof(double),1,fp);
  fwrite(&vstream_flag,sizeof(int),1,fp);
  if (vstream_flag) fwrite(vstream_user,sizeof(double),3,fp);
  fwrite(&temp_thermal_flag,sizeof(int),1,fp);
  if (temp_thermal_flag) fwrite(&temp_thermal_user,sizeof(double),1,fp);

  fwrite(fraction_flag,sizeof(int),nspecies,fp);
  fwrite(fraction_user,sizeof(double),nspecies,fp);

  fwrite(&ngroup,sizeof(int),1,fp);
  for (int i = 0; i < ngroup; i++) {
    int n = strlen(groups[i]) + 1;
    fwrite(&n,sizeof(int),1,fp);
    fwrite(groups[i],sizeof(char),n,fp);
  }
  fwrite(mix2group,sizeof(int),nspecies,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads mixture info from restart file
   bcast to other procs and all procs populate mixture with info
------------------------------------------------------------------------- */

void Mixture::read_restart(FILE *fp)
{
}

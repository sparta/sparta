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
  temp_rot_flag = 0;
  temp_vib_flag = 0;

  fraction = NULL;
  fraction_user = NULL;
  fraction_flag = NULL;
  cummulative = NULL;

  ngroup = maxgroup = 0;
  groups = NULL;
  groupsize = NULL;
  groupspecies = NULL;

  mix2group = NULL;
  species2group = NULL;
  species2species = NULL;

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
  memory->destroy(groupspecies);

  memory->destroy(mix2group);
  memory->destroy(species2group);
  memory->destroy(species2species);

  memory->destroy(vscale);
  memory->destroy(active);
}

/* ----------------------------------------------------------------------
   this empty mixture becomes clone of mixture old
------------------------------------------------------------------------- */

void Mixture::copy(Mixture *old)
{
  nrho = old->nrho;
  nrho_flag = old->nrho_flag;
  nrho_user = old->nrho_user;
  vstream[0] = old->vstream[0];
  vstream[1] = old->vstream[1];
  vstream[2] = old->vstream[2];
  vstream_flag = old->vstream_flag;
  vstream_user[0] = old->vstream_user[0];
  vstream_user[1] = old->vstream_user[1];
  vstream_user[2] = old->vstream_user[2];
  temp_thermal = old->temp_thermal;
  temp_thermal_flag = old->temp_thermal_flag;
  temp_thermal_user = old->temp_thermal_user;

  nspecies = maxspecies = old->nspecies;
  allocate();

  for (int i = 0; i < nspecies; i++) {
    species[i] = old->species[i];
    fraction[i] = old->fraction[i];
    fraction_flag[i] = old->fraction_flag[i];
    fraction_user[i] = old->fraction_user[i];
    mix2group[i] = old->mix2group[i];
  }

  for (int i = 0; i < old->ngroup; i++)
    add_group(old->groups[i]);
}

/* ----------------------------------------------------------------------
   process args of a mixture command
   zero of more species IDs followed by optional args
------------------------------------------------------------------------- */

void Mixture::command(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal mixture command");

  // nsp = # of species listed before optional keywords
  // iarg = start of optional keywords

  int iarg;
  for (iarg = 1; iarg < narg; iarg++) {
    if (strcmp(arg[iarg],"nrho") == 0) break;
    if (strcmp(arg[iarg],"vstream") == 0) break;
    if (strcmp(arg[iarg],"temp") == 0) break;
    if (strcmp(arg[iarg],"trot") == 0) break;
    if (strcmp(arg[iarg],"tvib") == 0) break;
    if (strcmp(arg[iarg],"frac") == 0) break;
    if (strcmp(arg[iarg],"group") == 0) break;
    if (strcmp(arg[iarg],"copy") == 0) break;
    if (strcmp(arg[iarg],"delete") == 0) break;
  }
  int nsp = iarg - 1;

  // add_species() processes list of species
  // params() processes remaining optional keywords

  add_species(nsp,&arg[1]);
  params(narg-iarg,&arg[iarg]);

  // if copy keyword was used, create a new mixture via add_mixture()
  // then invoke its copy() method, passing it this mixture

  if (copyflag) {
    particle->add_mixture(1,&arg[iarg+copyarg]);
    particle->mixture[particle->nmixture-1]->copy(this);
  }
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
  if (temp_rot_flag) temp_rot = temp_rot_user;
  else temp_rot = temp_thermal;
  if (temp_vib_flag) temp_vib = temp_vib_user;
  else temp_vib = temp_thermal;

  // mixture temperarate cannot be 0.0 if streaming velocity is 0.0

  if (temp_thermal == 0.0 &&
      vstream[0] == 0.0 && vstream[1] == 0.0 && vstream[2] == 0.0)
    error->all(FLERR,"Mixture streaming velocity and "
               "temperature cannot both be zero");

  // initialize all per-species fraction and cummulative values
  // account for both explicitly and implicitly set fractions

  int err = init_fraction(fraction_flag,fraction_user,fraction,cummulative);

  if (err) {
    char str[128];
    sprintf(str,"Mixture %s fractions exceed 1.0",id);
    error->all(FLERR,str);
  }

  // vscale = factor to scale Gaussian unit variance by
  //          to get thermal distribution of velocities
  // per-species value since includes species mass

  for (int i = 0; i < nspecies; i++) {
    int index = species[i];
    vscale[i] = sqrt(2.0 * update->boltz * temp_thermal /
                     particle->species[index].mass);
  }

  // setup species2group and species2species

  memory->destroy(species2group);
  memory->create(species2group,particle->nspecies,"mixture:species2group");
  for (int i = 0; i < particle->nspecies; i++) species2group[i] = -1;
  for (int i = 0; i < nspecies; i++) species2group[species[i]] = mix2group[i];

  memory->destroy(species2species);
  memory->create(species2species,particle->nspecies,"mixture:species2group");
  for (int i = 0; i < particle->nspecies; i++) species2species[i] = -1;
  for (int i = 0; i < nspecies; i++) species2species[species[i]] = i;

  // setup groupsize

  delete [] groupsize;
  groupsize = new int[ngroup];
  for (int i = 0; i < ngroup; i++) groupsize[i] = 0;
  for (int i = 0; i < nspecies; i++) groupsize[mix2group[i]]++;

  // setup groupspecies

  memory->destroy(groupspecies);
  memory->create_ragged(groupspecies,ngroup,groupsize,"mixture:groupspecies");

  for (int i = 0; i < ngroup; i++) groupsize[i] = 0;
  for (int i = 0; i < nspecies; i++)
    groupspecies[mix2group[i]][groupsize[mix2group[i]]++] = species[i];
}

/* ----------------------------------------------------------------------
   set f = fraction and c = cummulative fraction for each species in mixture
   fflag[I] = 1 if species I fraction is set by fuser[I]
   otherwise species I fraction is an implicit value
   implicit value = 1/nimplicit of unset remainder,
     where nimplicit = # of unset species
   return 0 for success
   return 1 for error if sum of specified fractions > 1.0
   called by init() and also by FixInflowFile::interpolate()
------------------------------------------------------------------------- */

int Mixture::init_fraction(int *fflag, double *fuser, double *f, double *c)
{
  // sum = total frac for species with explicity set fractions
  // nimplicit = number of unset species

  double sum = 0.0;
  int nimplicit = 0;
  for (int i = 0; i < nspecies; i++) {
    if (fflag[i]) sum += fuser[i];
    else nimplicit++;
  }

  if (sum > 1.0) return 1;

  // fraction for each unset species = equal portion of unset remainder
  // cummulative = cummulative fraction across species

  for (int i = 0; i < nspecies; i++) {
    if (fflag[i]) f[i] = fuser[i];
    else f[i] = (1.0-sum) / nimplicit;
    if (i) c[i] = c[i-1] + f[i];
    else c[i] = f[i];
  }
  if (nspecies) c[nspecies-1] = 1.0;
  return 0;
}

/* ----------------------------------------------------------------------
   process list of species appearing in a mixture command
------------------------------------------------------------------------- */

void Mixture::add_species(int narg, char **arg)
{
  int i,j,index;

  // activeflag = 1 if species are listed
  // active[i] = 0 if current mixture species I is not in the list
  // active[i] = 1 if species I is in the list but already existed in mixture
  // active[i] = 2 if species I was just added b/c it was in the list
  // active[i] = 3 if species is being removed from mixture via delete keyword
  //             this flag is set in params()

  if (narg) activeflag = 1;
  else activeflag = 0;
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

  copyflag = 0;
  int deleteflag = 0;
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
      if (temp_thermal_user < 0.0)
        error->all(FLERR,"Illegal mixture command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"trot") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal mixture command");
      temp_rot_flag = 1;
      temp_rot_user = atof(arg[iarg+1]);
      if (temp_rot_user < 0.0)
        error->all(FLERR,"Illegal mixture command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"tvib") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal mixture command");
      temp_vib_flag = 1;
      temp_vib_user = atof(arg[iarg+1]);
      if (temp_vib_user < 0.0)
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

    } else if (strcmp(arg[iarg],"copy") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal mixture command");
      if (particle->find_mixture(arg[iarg+1]) >= 0)
        error->all(FLERR,"New mixture copy mixture already exists");
      copyflag = 1;
      copyarg = iarg+1;
      iarg += 2;

    } else if (strcmp(arg[iarg],"delete") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal mixture command");
      deleteflag = 1;
      if (activeflag)
        error->all(FLERR,"Mixture delete cannot list species before keyword");
      if (iarg != 0)
        error->all(FLERR,"Mixture delete must be only keyword");
      iarg = narg;

    } else error->all(FLERR,"Illegal mixture command");
  }

  // process delete keyword

  if (deleteflag) {
    int m,ispecies;
    for (iarg = 1; iarg < narg; iarg++) {
      ispecies = particle->find_species(arg[iarg]);
      if (ispecies < 0)
        error->all(FLERR,"Mixture delete species is not recognized");
      for (m = 0; m < nspecies; m++)
        if (species[m] == ispecies) break;
      if (m == nspecies)
        error->all(FLERR,"Mixture delete species is not in mixture");
      active[m] = 1;
    }

    int nspecies_original = nspecies;
    m = 0;
    for (int i = 0; i < nspecies_original; i++) {
      species[m] = species[i];
      mix2group[m] = mix2group[i];
      if (active[i]) nspecies--;
      else m++;
    }
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
  //   assign any listed species that are new each to group "default"

  if (groupflag) {
    if (strcmp(arg[grouparg],"SELF") == 0) {
      if (!activeflag) {
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
      if (!activeflag) {
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

  } else if (activeflag) {
    for (int i = 0; i < nspecies; i++) {
      if (active[i] != 2) continue;
      int igroup = find_group("default");
      if (igroup < 0) {
        add_group("default");
        igroup = ngroup-1;
      }
      mix2group[i] = igroup;
    }
  }

  // remove empty groups due to deleteflag or groupflag operations

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
  fwrite(&temp_rot_flag,sizeof(int),1,fp);
  if (temp_rot_flag) fwrite(&temp_rot_user,sizeof(double),1,fp);
  fwrite(&temp_vib_flag,sizeof(int),1,fp);
  if (temp_vib_flag) fwrite(&temp_vib_user,sizeof(double),1,fp);

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
  int tmp;

  int me = comm->me;

  if (me == 0) tmp = fread(&nspecies,sizeof(int),1,fp);
  MPI_Bcast(&nspecies,1,MPI_INT,0,world);

  while (nspecies > maxspecies) allocate();

  if (me == 0) tmp = fread(species,sizeof(int),nspecies,fp);
  MPI_Bcast(species,nspecies,MPI_INT,0,world);

  if (me == 0) tmp = fread(&nrho_flag,sizeof(int),1,fp);
  MPI_Bcast(&nrho_flag,1,MPI_INT,0,world);
  if (nrho_flag) {
    if (me == 0) tmp = fread(&nrho_user,sizeof(double),1,fp);
    MPI_Bcast(&nrho_user,1,MPI_DOUBLE,0,world);
  }
  if (me == 0) tmp = fread(&vstream_flag,sizeof(int),1,fp);
  MPI_Bcast(&vstream_flag,1,MPI_INT,0,world);
  if (vstream_flag) {
    if (me == 0) tmp = fread(vstream_user,sizeof(double),3,fp);
    MPI_Bcast(vstream_user,3,MPI_DOUBLE,0,world);
  }
  if (me == 0) tmp = fread(&temp_thermal_flag,sizeof(int),1,fp);
  MPI_Bcast(&temp_thermal_flag,1,MPI_INT,0,world);
  if (temp_thermal_flag) {
    if (me == 0) tmp = fread(&temp_thermal_user,sizeof(double),1,fp);
    MPI_Bcast(&temp_thermal_user,1,MPI_DOUBLE,0,world);
  }
  if (me == 0) tmp = fread(&temp_rot_flag,sizeof(int),1,fp);
  MPI_Bcast(&temp_rot_flag,1,MPI_INT,0,world);
  if (temp_rot_flag) {
    if (me == 0) tmp = fread(&temp_rot_user,sizeof(double),1,fp);
    MPI_Bcast(&temp_rot_user,1,MPI_DOUBLE,0,world);
  }
  if (me == 0) tmp = fread(&temp_vib_flag,sizeof(int),1,fp);
  MPI_Bcast(&temp_vib_flag,1,MPI_INT,0,world);
  if (temp_vib_flag) {
    if (me == 0) tmp = fread(&temp_vib_user,sizeof(double),1,fp);
    MPI_Bcast(&temp_vib_user,1,MPI_DOUBLE,0,world);
  }

  if (me == 0) tmp = fread(fraction_flag,sizeof(int),nspecies,fp);
  MPI_Bcast(fraction_flag,nspecies,MPI_INT,0,world);
  if (me == 0) tmp = fread(fraction_user,sizeof(double),nspecies,fp);
  MPI_Bcast(fraction_user,nspecies,MPI_DOUBLE,0,world);

  int ngroup_file;
  if (me == 0) tmp = fread(&ngroup_file,sizeof(int),1,fp);
  MPI_Bcast(&ngroup_file,1,MPI_INT,0,world);

  int n;
  char *id;

  for (int i = 0; i < ngroup_file; i++) {
    if (me == 0) tmp = fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    id = new char[n];
    if (me == 0) tmp = fread(id,sizeof(char),n,fp);
    MPI_Bcast(id,n,MPI_CHAR,0,world);
    add_group(id);
    delete [] id;
  }

  if (me == 0) tmp = fread(mix2group,sizeof(int),nspecies,fp);
  MPI_Bcast(mix2group,nspecies,MPI_INT,0,world);
}

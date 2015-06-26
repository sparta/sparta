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

#include "mpi.h"
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "ctype.h"
#include "particle.h"
#include "grid.h"
#include "update.h"
#include "comm.h"
#include "mixture.h"
#include "collide.h"
#include "random_mars.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{PKEEP,PINSERT,PDONE,PDISCARD,PENTRY,PEXIT,PSURF};  // several files
enum{NONE,DISCRETE,SMOOTH};            // several files
enum{INT,DOUBLE};                      // several files

#define DELTA 16384
#define DELTASPECIES 16
#define DELTAMIXTURE 8
#define MAXLINE 1024

// customize by adding an abbreviation string
// also add a check for the keyword in 2 places in add_species()

#define AIR "N O NO"

/* ---------------------------------------------------------------------- */

Particle::Particle(SPARTA *sparta) : Pointers(sparta)
{
  MPI_Comm_rank(world,&me);

  exist = sorted = 0;
  nglobal = 0;
  nlocal = maxlocal = 0;
  particles = NULL;

  nspecies = maxspecies  = 0;
  species = NULL;

  //maxgrid = 0;
  //cellcount = NULL;
  //first = NULL;
  maxsort = 0;
  next = NULL;

  // create two default mixtures

  nmixture = maxmixture = 0;
  mixture = NULL;

  char **newarg = new char*[1];
  newarg[0] = (char *) "all";
  add_mixture(1,newarg);
  newarg[0] = (char *) "species";
  add_mixture(1,newarg);
  delete [] newarg;

  // custom per-particle vectors/arrays

  ncustom = 0;
  ename = NULL;
  etype = esize = ewhich = NULL;

  ncustom_ivec = ncustom_iarray = 0;
  icustom_ivec = icustom_iarray = NULL;
  eivec = NULL;
  eiarray = NULL;
  eicol = NULL;

  ncustom_dvec = ncustom_darray = 0;
  icustom_dvec = icustom_darray = NULL;
  edvec = NULL;
  edarray = NULL;
  edcol = NULL;

  custom_restart_flag = NULL;

  // RNG for particle weighting

  wrandom = NULL;
}

/* ---------------------------------------------------------------------- */

Particle::~Particle()
{
  memory->sfree(species);
  for (int i = 0; i < nmixture; i++) delete mixture[i];
  memory->sfree(mixture);

  memory->sfree(particles);
  //memory->destroy(cellcount);
  //memory->destroy(first);
  memory->destroy(next);

  for (int i = 0; i < ncustom; i++) delete [] ename[i];
  memory->sfree(ename);
  memory->destroy(etype);
  memory->destroy(esize);
  memory->destroy(ewhich);

  for (int i = 0; i < ncustom_ivec; i++) memory->destroy(eivec[i]);
  for (int i = 0; i < ncustom_iarray; i++) memory->destroy(eiarray[i]);
  for (int i = 0; i < ncustom_dvec; i++) memory->destroy(edvec[i]);
  for (int i = 0; i < ncustom_darray; i++) memory->destroy(edarray[i]);

  memory->destroy(icustom_ivec);
  memory->destroy(icustom_iarray);
  memory->sfree(eivec);
  memory->sfree(eiarray);
  memory->destroy(eicol);
  memory->destroy(icustom_dvec);
  memory->destroy(icustom_darray);
  memory->sfree(edvec);
  memory->sfree(edarray);
  memory->destroy(edcol);

  delete wrandom;
}

/* ---------------------------------------------------------------------- */

void Particle::init()
{
  for (int i = 0; i < nmixture; i++) mixture[i]->init();

  // RNG for particle weighting

  if (!wrandom) {
    wrandom = new RanPark(update->ranmaster->uniform());
    double seed = update->ranmaster->uniform();
    wrandom->reset(seed,me,100);
  }

  // if first run after reading a restart file,
  // delete any custom particle attributes that have not been re-defined
  // use nactive since remove_custom() may alter ncustom

  if (custom_restart_flag) {
    int nactive = ncustom;
    for (int i = 0; i < nactive; i++)
      if (custom_restart_flag[i] == 0) remove_custom(i);
    delete [] custom_restart_flag;
    custom_restart_flag = NULL;
  }

  // reallocate cellcount and first lists as needed
  // NOTE: when grid becomes dynamic, will need to do this in sort()

  //if (maxgrid < grid->nlocal) {
  //  maxgrid = grid->nlocal;
    //    memory->destroy(cellcount);
    //memory->destroy(first);
    //memory->create(first,maxgrid,"particle:first");
    //memory->create(cellcount,maxgrid,"particle:cellcount");
  // }
}

/* ----------------------------------------------------------------------
   compress particle list to remove particles with indices in mlist
   mlist indices MUST be in ascending order
   overwrite deleted particle with particle from end of nlocal list
   inner while loop avoids overwrite with deleted particle at end of mlist
   called every step from Comm::migrate_particles() when particles migrate
   this is similar to compress_reactions(), but does not need
     an auxiliary vector b/c indices are in ascending order
------------------------------------------------------------------------- */

void Particle::compress_migrate(int nmigrate, int *mlist)
{
  int i,j,k;
  int nbytes = sizeof(OnePart);

  if (!ncustom) {
    for (i = 0; i < nmigrate; i++) {
      j = mlist[i];
      k = nlocal - 1;
      while (k == mlist[nmigrate-1] && k > j) {
        nmigrate--;
        nlocal--;
        k--;
      }
      nlocal--;
      if (j == k) continue;
      memcpy(&particles[j],&particles[k],nbytes);
    }
    
  } else {
    for (i = 0; i < nmigrate; i++) {
      j = mlist[i];
      k = nlocal - 1;
      while (k == mlist[nmigrate-1] && k > j) {
        nmigrate--;
        nlocal--;
        k--;
      }
      nlocal--;
      if (j == k) continue;
      memcpy(&particles[j],&particles[k],nbytes);
      copy_custom(j,k);
    }
  }

  sorted = 0;
}

/* ----------------------------------------------------------------------
   compress particle list to remove particles with icell < 0
   all particles MUST be in owned cells
   overwrite deleted particle with particle from end of nlocal list
   called from Comm::migrate_cells() when cells+particles migrate on rebalance
   called from AdaptGrid when particles are sent to other procs
   called from ReadSurf to remove particles from cells with surfs
   this does NOT preserve particle sorting
------------------------------------------------------------------------- */

void Particle::compress_rebalance()
{
  int nbytes = sizeof(OnePart);

  if (!ncustom) {
    int i = 0;
    while (i < nlocal) {
      if (particles[i].icell < 0) {
	memcpy(&particles[i],&particles[nlocal-1],nbytes);
	nlocal--;
      } else i++;
    }

  } else {
    int i = 0;
    while (i < nlocal) {
      if (particles[i].icell < 0) {
	memcpy(&particles[i],&particles[nlocal-1],nbytes);
	copy_custom(i,nlocal-1);
	nlocal--;
      } else i++;
    }
  }

  sorted = 0;
}

/* ----------------------------------------------------------------------
   compress particle list to remove particles with indices in dellist
   dellist indices can be in ANY order
   overwrite deleted particle with particle from end of nlocal list
   use of next vector does bookkeeping for particles
     that are moved from their original location before they are deleted
   called from Collide::migrate_particles() each timestep
     if any particles were deleted by gas-phase collision reactions
   this is similar to compress_migrate(), but needs to use
     an auxiliary vector b/c indices are in random order
------------------------------------------------------------------------- */

void Particle::compress_reactions(int ndelete, int *dellist)
{
  int i;

  int nbytes = sizeof(OnePart);

  // use next as a scratch vector
  // is always defined when performing collisions and thus gas reactions
  // next is only used for upper locs from nlocal-ndelete to nlocal
  // next[i] = current index of atom originally at index i, when i >= nlocal
  // next[i] = original index of atom currently at index i, when i <= nlocal

  int upper = nlocal-ndelete;
  for (i = upper; i < nlocal; i++) next[i] = i;

  // i = current index of atom to remove, even if it previously moved

  if (!ncustom) {
    for (int m = 0; m < ndelete; m++) {
      i = dellist[m];
      if (i >= nlocal) i = next[i];
      nlocal--;
      if (i == nlocal) continue;
      memcpy(&particles[i],&particles[nlocal],nbytes);
      if (i >= upper) next[i] = next[nlocal];
      next[next[nlocal]] = i;
    }

  } else {
    for (int m = 0; m < ndelete; m++) {
      i = dellist[m];
      if (i >= nlocal) i = next[i];
      nlocal--;
      if (i == nlocal) continue;
      memcpy(&particles[i],&particles[nlocal],nbytes);
      copy_custom(i,nlocal);
      if (i >= upper) next[i] = next[nlocal];
      next[next[nlocal]] = i;
    }
  }
}

/* ----------------------------------------------------------------------
   sort particles into grid cells
   set cinfo.first = index of first particle in cell
   set cinfo.count = # of particles in cell
   next[] = index of next particle in same cell, -1 for no more
------------------------------------------------------------------------- */

void Particle::sort()
{
  sorted = 1;

  // reallocate next list as needed

  if (maxsort < maxlocal) {
    maxsort = maxlocal;
    memory->destroy(next);
    memory->create(next,maxsort,"particle:next");
  }

  // initialize linked list of particles in cells I own

  Grid::ChildInfo *cinfo = grid->cinfo;
  int nglocal = grid->nlocal;

  for (int icell = 0; icell < nglocal; icell++) {
    cinfo[icell].first = -1;
    cinfo[icell].count = 0;
    //cellcount[i] = 0;
    //first[i] = -1;
  }

  // reverse loop stores linked list in forward order
  // icell = global cell the particle is in

  int icell;
  for (int i = nlocal-1; i >= 0; i--) {
    icell = particles[i].icell;
    next[i] = cinfo[icell].first;
    cinfo[icell].first = i;
    cinfo[icell].count++;

    // NOTE: this method seems much slower for some reason
    // uses separate, smaller vectors for first & cellcount
    //icell = cells[particles[i].icell].local;
    //next[i] = first[icell];
    //first[icell] = i;
    //cellcount[icell]++;
  }
}

/* ----------------------------------------------------------------------
   mark all particles in a grid cell, stating with index ip, with icell = -1
   so can be deleted by compress_rebalance()
   assume particles are already sorted
   called from ReadSurf to remove particles from a grid cell with surfaces
------------------------------------------------------------------------- */

void Particle::remove_all_from_cell(int ip)
{
  while (ip >= 0) {
    particles[ip].icell = -1;
    ip = next[ip];
  }
}

/* ----------------------------------------------------------------------
   set the initial weight of each particle
   called by Update before particle move
   only called if particle weighting is enabled
   only grid-based weighting is currently implemented
------------------------------------------------------------------------- */

void Particle::pre_weight()
{
  int icell;
  Grid::ChildInfo *cinfo = grid->cinfo;

  for (int i = 0; i < nlocal; i++) {
    icell = particles[i].icell;
    particles[i].weight = cinfo[icell].weight;
  }
}

/* ----------------------------------------------------------------------
   clone/delete each particle based on ratio of its initial/final weights
   called by Update after particle move and migration
   only called if particle weighting is enabled
   only grid-based weighting is currently implemented
------------------------------------------------------------------------- */

void Particle::post_weight()
{
  int m,icell,nclone;
  double ratio,fraction;

  int nbytes = sizeof(OnePart);
  Grid::ChildInfo *cinfo = grid->cinfo;

  // nlocal_original-1 = index of last original particle

  int nlocal_original = nlocal;
  int i = 0;

  while (i < nlocal_original) {
    icell = particles[i].icell;

    // next particle will be an original particle
    // skip it if no weight change

    if (particles[i].weight == cinfo[icell].weight) {
      i++;
      continue;
    }

    // ratio < 1.0 is candidate for deletion
    // if deleted and particle that takes its place is cloned (Nloc > Norig)
    //   then skip it via i++, else will examine it on next iteration

    ratio = particles[i].weight / cinfo[icell].weight;

    if (ratio < 1.0) {
      if (wrandom->uniform() > ratio) {
        memcpy(&particles[i],&particles[nlocal-1],nbytes);
	if (ncustom) copy_custom(i,nlocal-1);
        if (nlocal > nlocal_original) i++;
        else nlocal_original--;
        nlocal--;
      } else i++;
      continue;
    }

    // ratio > 1.0 is candidate for cloning
    // create Nclone new particles each with unique ID 

    nclone = static_cast<int> (ratio);
    fraction = ratio - nclone;
    nclone--;
    if (wrandom->uniform() < fraction) nclone++;

    for (m = 0; m < nclone; m++) {
      clone_particle(i);
      particles[nlocal-1].id = MAXSMALLINT*wrandom->uniform();
    }
    i++;
  }
}

/* ----------------------------------------------------------------------
   insure particle list can hold nextra new particles
   if defined, also grow custom particle arrays and initialize with zeroes
------------------------------------------------------------------------- */

void Particle::grow(int nextra)
{
  bigint target = (bigint) nlocal + nextra;
  if (target <= maxlocal) return;
  
  int oldmax = maxlocal;
  bigint newmax = maxlocal;
  while (newmax < target) newmax += DELTA;
  
  if (newmax > MAXSMALLINT) 
    error->one(FLERR,"Per-processor particle count is too big");

  maxlocal = newmax;
  particles = (OnePart *)
    memory->srealloc(particles,maxlocal*sizeof(OnePart),
		     "particle:particles");
  memset(&particles[oldmax],0,(maxlocal-oldmax)*sizeof(OnePart));

  if (ncustom == 0) return;

  for (int i = 0; i < ncustom; i++) {
    if (ename[i] == NULL) continue;
    grow_custom(i,oldmax,maxlocal);
  }
}

/* ----------------------------------------------------------------------
   add a particle to particle list
   return 1 if particle array was reallocated, else 0
------------------------------------------------------------------------- */

int Particle::add_particle(int id, int ispecies, int icell,
                           double *x, double *v, double erot, double evib)
{
  int reallocflag = 0;
  if (nlocal == maxlocal) {
    grow(1);
    reallocflag = 1;
  }

  OnePart *p = &particles[nlocal];
  
  p->id = id;
  p->ispecies = ispecies;
  p->icell = icell;
  p->x[0] = x[0];
  p->x[1] = x[1];
  p->x[2] = x[2];
  p->v[0] = v[0];
  p->v[1] = v[1];
  p->v[2] = v[2];
  p->erot = erot;
  p->evib = evib;
  p->flag = PKEEP;

  //p->dtremain = 0.0;    not needed due to memset in grow() ??
  //p->weight = 1.0;      not needed due to memset in grow() ??

  nlocal++;
  return reallocflag;
}

/* ----------------------------------------------------------------------
   clone particle Index and add to end of particle list
   return 1 if particle array was reallocated, else 0
------------------------------------------------------------------------- */

int Particle::clone_particle(int index)
{
  int reallocflag = 0;
  if (nlocal == maxlocal) {
    grow(1);
    reallocflag = 1;
  }

  memcpy(&particles[nlocal],&particles[index],sizeof(OnePart));
  if (ncustom) copy_custom(nlocal,index);

  nlocal++;
  return reallocflag;
}

/* ----------------------------------------------------------------------
   add one or more species to species list
------------------------------------------------------------------------- */

void Particle::add_species(int narg, char **arg)
{
  int i,j,n;;

  if (narg < 2) error->all(FLERR,"Illegal species command");

  if (me == 0) {
    fp = fopen(arg[0],"r");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open species file %s",arg[0]);
      error->one(FLERR,str);
    }
  }

  // filespecies = list of species defined in file

  nfilespecies = maxfilespecies = 0;
  filespecies = NULL;

  if (me == 0) read_species_file();
  MPI_Bcast(&nfilespecies,1,MPI_INT,0,world);
  if (nfilespecies >= maxfilespecies) {
    memory->destroy(filespecies);
    maxfilespecies = nfilespecies;
    filespecies = (Species *) 
      memory->smalloc(maxfilespecies*sizeof(Species),
		      "particle:filespecies");
  }
  MPI_Bcast(filespecies,nfilespecies*sizeof(Species),MPI_BYTE,0,world);

  // newspecies = # of new user-requested species
  // names = list of new species IDs
  // customize abbreviations by adding new keyword in 2 places

  char line[MAXLINE];

  int newspecies = 0;
  for (int iarg = 1; iarg < narg; iarg++) {
    if (strcmp(arg[iarg],"air") == 0) {
      strcpy(line,AIR);
      newspecies += wordcount(line,NULL);
    } else newspecies++;
  }

  char **names = new char*[newspecies];
  newspecies = 0;

  for (int iarg = 1; iarg < narg; iarg++) {
    if (strcmp(arg[iarg],"air") == 0) {
      strcpy(line,AIR);
      newspecies += wordcount(line,&names[newspecies]);
    } else names[newspecies++] = arg[iarg];
  }

  // species ID must be all alphanumeric chars, underscore, plus/minus

  for (i = 0; i < newspecies; i++) {
    n = strlen(names[i]);
    for (j = 0; j < n-1; j++)
      if (!isalnum(names[i][j]) && names[i][j] != '_' && 
          names[i][j] != '+' && names[i][j] != '-')
        error->all(FLERR,"Invalid character in species ID");
  }

  // extend species list if necessary

  if (nspecies + newspecies > maxspecies) {
    while (nspecies+newspecies > maxspecies) maxspecies += DELTASPECIES;
    species = (Species *) 
      memory->srealloc(species,maxspecies*sizeof(Species),"particle:species");
  }

  // extract info on user-requested species from file species list
  // add new species to default mixtures "all" and "species"

  int imix_all = find_mixture((char *) "all");
  int imix_species = find_mixture((char *) "species");

  for (i = 0; i < newspecies; i++) {
    for (j = 0; j < nspecies; j++)
      if (strcmp(names[i],species[j].id) == 0) break;
    if (j < nspecies) error->all(FLERR,"Species ID is already defined");
    for (j = 0; j < nfilespecies; j++)
      if (strcmp(names[i],filespecies[j].id) == 0) break;
    if (j == nfilespecies)
      error->all(FLERR,"Species ID does not appear in species file");
    memcpy(&species[nspecies],&filespecies[j],sizeof(Species));
    nspecies++;

    mixture[imix_all]->add_species_default(species[nspecies-1].id);
    mixture[imix_species]->add_species_default(species[nspecies-1].id);
  }

  memory->sfree(filespecies);
  delete [] names;
}

/* ----------------------------------------------------------------------
   create or augment a mixture of species
------------------------------------------------------------------------- */

void Particle::add_mixture(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal mixture command");

  // imix = index if mixture ID already exists
  // else instantiate a new mixture

  int imix = find_mixture(arg[0]);

  if (imix < 0) {
    if (nmixture == maxmixture) {
      maxmixture += DELTAMIXTURE;
      mixture = (Mixture **) memory->srealloc(mixture,
                                              maxmixture*sizeof(Mixture *),
                                              "particle:mixture");
    }
    imix = nmixture;
    mixture[nmixture++] = new Mixture(sparta,arg[0]);
  }

  mixture[imix]->command(narg,arg);
}

/* ----------------------------------------------------------------------
   return index of ID in list of species IDs
   return -1 if not found
------------------------------------------------------------------------- */

int Particle::find_species(char *id)
{
  for (int i = 0; i < nspecies; i++)
    if (strcmp(id,species[i].id) == 0) return i;
  return -1;
}

/* ----------------------------------------------------------------------
   return index of ID in list of mixture IDs
   return -1 if not found
------------------------------------------------------------------------- */

int Particle::find_mixture(char *id)
{
  for (int i = 0; i < nmixture; i++)
    if (strcmp(id,mixture[i]->id) == 0) return i;
  return -1;
}

/* ----------------------------------------------------------------------
   generate random rotational energy for a particle
   only a function of species index and species properties
------------------------------------------------------------------------- */

double Particle::erot(int isp, double temp_thermal, RanPark *erandom)
{
 double eng,a,erm,b;

 if (!collide || collide->rotstyle == NONE) return 0.0;
 if (species[isp].rotdof < 2) return 0.0;

 if (species[isp].rotdof == 2)
   eng = -log(erandom->uniform()) * update->boltz * temp_thermal;
 else {
   a = 0.5*particle->species[isp].rotdof-1.0;
   while (1) {
     // energy cut-off at 10 kT
     erm = 10.0*erandom->uniform();
     b = pow(erm/a,a) * exp(a-erm);
     if (b > erandom->uniform()) break;
   }
   eng = erm * update->boltz * temp_thermal;
 }

 return eng;
}

/* ----------------------------------------------------------------------
   generate random vibrational energy for a particle
   only a function of species index and species properties
------------------------------------------------------------------------- */

double Particle::evib(int isp, double temp_thermal, RanPark *erandom)
{
  double eng,a,erm,b;

  int vibstyle = NONE;
  if (collide) vibstyle = collide->vibstyle;
  if (vibstyle == NONE || species[isp].vibdof < 2) return 0.0;

  eng = 0.0;
  if (vibstyle == DISCRETE && species[isp].vibdof == 2) {
    int ivib = -log(erandom->uniform()) * temp_thermal / 
      particle->species[isp].vibtemp;
    eng = ivib * update->boltz * particle->species[isp].vibtemp;
  } else if (vibstyle == SMOOTH || species[isp].vibdof >= 2) {
    if (species[isp].vibdof == 2)
      eng = -log(erandom->uniform()) * update->boltz * temp_thermal;
    else if (species[isp].vibdof > 2) {
      a = 0.5*particle->species[isp].vibdof-1.;
      while (1) {
        // energy cut-off at 10 kT
        erm = 10.0*erandom->uniform();
        b = pow(erm/a,a) * exp(a-erm);
        if (b > erandom->uniform()) break;
      }
      eng = erm * update->boltz * temp_thermal;
    }
  }

  return eng;
}

/* ----------------------------------------------------------------------
   read list of species defined in species file
   store info in filespecies and nfilespecies
   only invoked by proc 0
------------------------------------------------------------------------- */

void Particle::read_species_file()
{
  nfilespecies = maxfilespecies = 0;
  filespecies = NULL;

  // read file line by line
  // skip blank lines or comment lines starting with '#'
  // all other lines must have NWORDS 

  int NWORDS = 10;
  char **words = new char*[NWORDS];
  char line[MAXLINE],copy[MAXLINE];

  while (fgets(line,MAXLINE,fp)) {
    int pre = strspn(line," \t\n");
    if (pre == strlen(line) || line[pre] == '#') continue;

    strcpy(copy,line);
    int nwords = wordcount(copy,NULL);
    if (nwords != NWORDS)
      error->one(FLERR,"Incorrect line format in species file");

    if (nfilespecies == maxfilespecies) {
      maxfilespecies += DELTASPECIES;
      filespecies = (Species *) 
	memory->srealloc(filespecies,maxfilespecies*sizeof(Species),
			 "particle:filespecies");
      memset(&filespecies[nfilespecies],0,
             (maxfilespecies-nfilespecies)*sizeof(Species));
    }

    nwords = wordcount(line,words);
    Species *fsp = &filespecies[nfilespecies];

    if (strlen(words[0]) + 1 > 16)
      error->one(FLERR,"Invalid species ID in species file");
    strcpy(fsp->id,words[0]);

    fsp->molwt = atof(words[1]);
    fsp->mass = atof(words[2]);
    fsp->rotdof = atoi(words[3]);
    fsp->rotrel = atof(words[4]);
    fsp->vibdof = atoi(words[5]);
    fsp->vibrel = atof(words[6]);
    fsp->vibtemp = atof(words[7]);
    fsp->specwt = atof(words[8]);
    fsp->charge = atof(words[9]);

    if (fsp->rotdof > 0 || fsp->vibdof > 0) fsp->internaldof = 1;
    else fsp->internaldof = 0;

    nfilespecies++;
  }

  delete [] words;

  fclose(fp);
}

/* ----------------------------------------------------------------------
   count whitespace-delimited words in line
   line will be modified, since strtok() inserts NULLs
   if words is non-NULL, store ptr to each word
------------------------------------------------------------------------- */

int Particle::wordcount(char *line, char **words)
{
  int nwords = 0;
  char *word = strtok(line," \t");

  while (word) {
    if (words) words[nwords] = word;
    nwords++;
    word = strtok(NULL," \t");
  }

  return nwords;
}

/* ----------------------------------------------------------------------
   proc 0 writes species info to restart file
------------------------------------------------------------------------- */

void Particle::write_restart_species(FILE *fp)
{
  fwrite(&nspecies,sizeof(int),1,fp);
  fwrite(species,sizeof(Species),nspecies,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads species info from restart file
   bcast to other procs
------------------------------------------------------------------------- */

void Particle::read_restart_species(FILE *fp)
{
  if (me == 0) fread(&nspecies,sizeof(int),1,fp);
  MPI_Bcast(&nspecies,1,MPI_INT,0,world);

  if (nspecies > maxspecies) {
    while (nspecies > maxspecies) maxspecies += DELTASPECIES;
    species = (Species *) 
      memory->srealloc(species,maxspecies*sizeof(Species),"particle:species");
  }

  if (me == 0) fread(species,sizeof(Species),nspecies,fp);
  MPI_Bcast(species,nspecies*sizeof(Species),MPI_CHAR,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes mixture info to restart file
------------------------------------------------------------------------- */

void Particle::write_restart_mixture(FILE *fp)
{
  fwrite(&nmixture,sizeof(int),1,fp);
  for (int i = 0; i < nmixture; i++) {
    int n = strlen(mixture[i]->id) + 1;
    fwrite(&n,sizeof(int),1,fp);
    fwrite(mixture[i]->id,sizeof(char),n,fp);
    mixture[i]->write_restart(fp);
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads mixture info from restart file
   bcast to other procs and all procs instantiate series of Mixtures
------------------------------------------------------------------------- */

void Particle::read_restart_mixture(FILE *fp)
{
  // must first clear existing default mixtures

  for (int i = 0; i < nmixture; i++) delete mixture[i];
  nmixture = 0;

  // now process restart file data

  if (me == 0) fread(&nmixture,sizeof(int),1,fp);
  MPI_Bcast(&nmixture,1,MPI_INT,0,world);

  if (nmixture > maxmixture) {
    while (nmixture > maxmixture) maxmixture += DELTAMIXTURE;
    mixture = (Mixture **) 
      memory->srealloc(mixture,maxmixture*sizeof(Mixture *),"particle:mixture");
  }

  int n;
  char *id;

  for (int i = 0; i < nmixture; i++) {
    if (me == 0) fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    id = new char[n];
    if (me == 0) fread(id,sizeof(char),n,fp);
    MPI_Bcast(id,n,MPI_CHAR,0,world);
    mixture[i] = new Mixture(sparta,id);
    mixture[i]->read_restart(fp);
    delete [] id;
  }
}

/* ----------------------------------------------------------------------
   return size of particle restart info for this proc
   NOTE: worry about N overflowing int and IROUNDUP ???
------------------------------------------------------------------------- */

int Particle::size_restart()
{
  int n = sizeof(int);
  n = IROUNDUP(n);
  n += nlocal * sizeof(OnePartRestart);
  n += nlocal * sizeof_custom();
  n = IROUNDUP(n);
  return n;
}

/* ----------------------------------------------------------------------
   pack my particle info into buf
   use OnePartRestart data struct for permanent info and to encode cell ID
   include per-particle custom attributes if defined
------------------------------------------------------------------------- */

int Particle::pack_restart(char *buf)
{
  Grid::ChildCell *cells = grid->cells;
  OnePart *p;
  OnePartRestart *pr;
  int nbytes_custom = sizeof_custom();

  char *ptr = buf;
  int *ibuf = (int *) ptr;
  *ibuf = nlocal;
  ptr += sizeof(int);
  ptr = ROUNDUP(ptr);

  for (int i = 0; i < nlocal; i++) {
    p = &particles[i];
    pr = (OnePartRestart *) ptr;
    pr->id = p->id;
    pr->ispecies = p->ispecies;
    pr->icell = cells[p->icell].id;
    pr->nsplit = cells[p->icell].nsplit;
    pr->x[0] = p->x[0];
    pr->x[1] = p->x[1];
    pr->x[2] = p->x[2];
    pr->v[0] = p->v[0];
    pr->v[1] = p->v[1];
    pr->v[2] = p->v[2];
    pr->erot = p->erot;
    pr->evib = p->evib;

    ptr += sizeof(OnePartRestart);
    if (!ncustom) continue;

    pack_custom(i,ptr);
    ptr += nbytes_custom;
  }

  ptr = ROUNDUP(ptr);
  return ptr - buf;
}

/* ----------------------------------------------------------------------
   unpack particle info into restart storage
   allocate data structure here, will be deallocated by ReadRestart
   include per-particle custom attributes if defined
------------------------------------------------------------------------- */

int Particle::unpack_restart(char *buf)
{
  int nbytes_particle = sizeof(OnePartRestart);
  int nbytes_custom = sizeof_custom();
  int nbytes = nbytes_particle + nbytes_custom;

  char *ptr = buf;
  int *ibuf = (int *) buf;
  nlocal_restart = *ibuf;
  ptr += sizeof(int);
  ptr = ROUNDUP(ptr);

  particle_restart = (char *) 
    memory->smalloc(nlocal_restart*nbytes,"grid:particle_restart");

  memcpy(particle_restart,ptr,nlocal_restart*nbytes);
  ptr += nlocal_restart * sizeof(OnePartRestart);
  ptr = ROUNDUP(ptr);

  return ptr - buf;
}

// ----------------------------------------------------------------------
// methods for per-particle custom attributes
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   find custom per-atom vector/array with name
   return index if found
   return -1 if not found
------------------------------------------------------------------------- */

int Particle::find_custom(char *name)
{
  for (int i = 0; i < ncustom; i++)
    if (ename[i] && strcmp(ename[i],name) == 0) return i;
  return -1;
}

/* ----------------------------------------------------------------------
   add a custom attribute with name
   assumes name does not already exist, except in case of restart
   type = 0/1 for int/double
   size = 0 for vector, size > 0 for array with size columns
   allocate the vector or array to current maxlocal via grow_custom()
   return index of its location;
------------------------------------------------------------------------- */

int Particle::add_custom(char *name, int type, int size)
{
  int index;

  // if name already exists
  // just return index if a restart script and re-defining the name
  // else error

  index = find_custom(name);
  if (index >= 0) {
    if (custom_restart_flag == NULL || custom_restart_flag[index] == 1)
      error->all(FLERR,"Custom particle attribute name already exists");
    custom_restart_flag[index] = 1;
    return index;
  }

  // use first available NULL entry or allocate a new one

  for (index = 0; index < ncustom; index++)
    if (ename[index] == NULL) break;

  if (index == ncustom) {
    ncustom++;
    ename = (char **) memory->srealloc(ename,ncustom*sizeof(char *),
                                       "particle:ename");
    memory->grow(etype,ncustom,"particle:etype");
    memory->grow(esize,ncustom,"particle:etype");
    memory->grow(ewhich,ncustom,"particle:etype");
  }

  int n = strlen(name) + 1;
  ename[index] = new char[n];
  strcpy(ename[index],name);
  etype[index] = type;
  esize[index] = size;

  if (type == INT) {
    if (size == 0) {
      ewhich[index] = ncustom_ivec++;
      eivec = (int **) 
        memory->srealloc(eivec,ncustom_ivec*sizeof(int *),"particle:eivec");
      eivec[ncustom_ivec-1] = NULL;
      memory->grow(icustom_ivec,ncustom_ivec,"particle:icustom_ivec");
      icustom_ivec[ncustom_ivec-1] = index;
    } else {
      ewhich[index] = ncustom_iarray++;
      eiarray = (int ***) 
        memory->srealloc(eiarray,ncustom_iarray*sizeof(int **),
                         "particle:eivec");
      eiarray[ncustom_iarray-1] = NULL;
      memory->grow(icustom_iarray,ncustom_iarray,"particle:icustom_iarray");
      icustom_iarray[ncustom_iarray-1] = index;
      memory->grow(eicol,ncustom_iarray,"particle:eicol");
      eicol[ncustom_iarray-1] = size;
    }
  } else if (type == DOUBLE) {
    if (size == 0) {
      ewhich[index] = ncustom_dvec++;
      edvec = (double **) 
        memory->srealloc(edvec,ncustom_dvec*sizeof(double *),"particle:edvec");
      edvec[ncustom_dvec-1] = NULL;
      memory->grow(icustom_dvec,ncustom_dvec,"particle:icustom_dvec");
      icustom_dvec[ncustom_dvec-1] = index;
    } else {
      ewhich[index] = ncustom_darray++;
      edarray = (double ***) 
        memory->srealloc(edarray,ncustom_darray*sizeof(double **),
                         "particle:edvec");
      edarray[ncustom_darray-1] = NULL;
      memory->grow(icustom_darray,ncustom_darray,"particle:icustom_darray");
      icustom_darray[ncustom_darray-1] = index;
      memory->grow(edcol,ncustom_darray,"particle:edcol");
      edcol[ncustom_darray-1] = size;
    }
  }

  grow_custom(index,0,maxlocal);

  return index;
}

/* ----------------------------------------------------------------------
   grow the vector/array associated with custom attribute with index
   nold = old length, nnew = new length (typically maxlocal)
   set new values to 0 via memset()
------------------------------------------------------------------------- */

void Particle::grow_custom(int index, int nold, int nnew)
{
  if (etype[index] == INT) {
    if (esize[index] == 0) {
      int *ivector = eivec[ewhich[index]];
      memory->grow(ivector,nnew,"particle:eivec");
      if (ivector) memset(&ivector[nold],0,(nnew-nold)*sizeof(int));
      eivec[ewhich[index]] = ivector;
    } else {
      int **iarray = eiarray[ewhich[index]];
      memory->grow(iarray,nnew,esize[index],"particle:eiarray");
      if (iarray) 
        memset(&iarray[nold][0],0,(nnew-nold)*esize[index]*sizeof(int));
      eiarray[ewhich[index]] = iarray;
    }

  } else {
    if (esize[index] == 0) {
      double *dvector = edvec[ewhich[index]];
      memory->grow(dvector,nnew,"particle:edvec");
      if (dvector) memset(&dvector[nold],0,(nnew-nold)*sizeof(double));
      edvec[ewhich[index]] = dvector;
    } else {
      double **darray = edarray[ewhich[index]];
      memory->grow(darray,nnew,esize[index],"particle:edarray");
      if (darray)
        memset(&darray[nold][0],0,(nnew-nold)*esize[index]*sizeof(double));
      edarray[ewhich[index]] = darray;
    }
  }
}

/* ----------------------------------------------------------------------
   remove a custom attribute at location index
   free memory for name and vector/array and set ptrs to NULL
   ncustom lists never shrink, but indices stored between
     the ncustom list and the dense vector/array lists must be reset
------------------------------------------------------------------------- */

void Particle::remove_custom(int index)
{
  delete [] ename[index];
  ename[index] = NULL;

  if (etype[index] == INT) {
    if (esize[index] == 0) {
      memory->destroy(eivec[ewhich[index]]);
      ncustom_ivec--;
      for (int i = ewhich[index]; i < ncustom_ivec; i++) {
        icustom_ivec[i] = icustom_ivec[i+1];
        ewhich[icustom_ivec[i]] = i;
        eivec[i] = eivec[i+1];
      }
    } else{
      memory->destroy(eiarray[ewhich[index]]);
      ncustom_iarray--;
      for (int i = ewhich[index]; i < ncustom_iarray; i++) {
        icustom_iarray[i] = icustom_iarray[i+1];
        ewhich[icustom_iarray[i]] = i;
        eiarray[i] = eiarray[i+1];
        eicol[i] = eicol[i+1];
      }
    }
  } else if (etype[index] == DOUBLE) {
    if (esize[index] == 0) {
      memory->destroy(edvec[ewhich[index]]);
      ncustom_dvec--;
      for (int i = ewhich[index]; i < ncustom_dvec; i++) {
        icustom_dvec[i] = icustom_dvec[i+1];
        ewhich[icustom_dvec[i]] = i;
        edvec[i] = edvec[i+1];
      }
    } else{
      memory->destroy(edarray[ewhich[index]]);
      ncustom_darray--;
      for (int i = ewhich[index]; i < ncustom_darray; i++) {
        icustom_darray[i] = icustom_darray[i+1];
        ewhich[icustom_darray[i]] = i;
        edarray[i] = edarray[i+1];
        edcol[i] = edcol[i+1];
      }
    }
  }

  // set ncustom = 0 if custom list is now entirely empty

  int empty = 1;
  for (int i = 0; i < ncustom; i++) 
    if (ename[i]) empty = 0;
  if (empty) ncustom = 0;
}

/* ----------------------------------------------------------------------
   copy info for one particle in custom attribute vectors/arrays
   into location I from location J
------------------------------------------------------------------------- */

void Particle::copy_custom(int i, int j)
{
  int m;

  // caller does not always check this
  // shouldn't be a problem, but valgrind can complain if memcpy to self
  // oddly memcpy(&particles[i],&particles[j],sizeof(OnePart)) seems OK

  if (i == j) return;

  // 4 flavors of vectors/arrays

  if (ncustom_ivec) {
    for (m = 0; m < ncustom_ivec; m++) eivec[m][i] = eivec[m][j];
  }
  if (ncustom_iarray) {
    for (m = 0; m < ncustom_iarray; m++)
      memcpy(eiarray[m][i],eiarray[m][j],eicol[m]*sizeof(int));
  }
  if (ncustom_dvec) {
    for (m = 0; m < ncustom_dvec; m++) edvec[m][i] = edvec[m][j];
  }
  if (ncustom_darray) {
    for (m = 0; m < ncustom_darray; m++)
      memcpy(edarray[m][i],edarray[m][j],edcol[m]*sizeof(double));
  }
}

/* ----------------------------------------------------------------------
   return size of all custom attributes in bytes for one particle
   used by callers to allocate buffer memory for particles
   assume integer attributes can be put at start of buffer
   only alignment needed is between integers and doubles
------------------------------------------------------------------------- */

int Particle::sizeof_custom()
{
  int n = 0;

  n += ncustom_ivec*sizeof(int);
  if (ncustom_iarray)
    for (int i = 0; i < ncustom_iarray; i++)
      n += eicol[i]*sizeof(int);

  n = IROUNDUP(n);

  n += ncustom_dvec*sizeof(double);
  if (ncustom_darray)
    for (int i = 0; i < ncustom_darray; i++)
      n += edcol[i]*sizeof(double);

  return n;
}

/* ----------------------------------------------------------------------
   proc 0 writes custom attribute definition info to restart file
------------------------------------------------------------------------- */

void Particle::write_restart_custom(FILE *fp)
{
  int m,index;

  // nactive = # of ncustom that have active vectors/arrays

  int nactive = 0;
  for (int i = 0; i < ncustom; i++)
    if (ename[i]) nactive++;

  fwrite(&nactive,sizeof(int),1,fp);

  // must write custom info in same order 
  //   the per-particle custom values will be written into file
  // not necessarily the same as ncustom list, due to deletions & additions

  for (m = 0; m < ncustom_ivec; m++) {
    index = icustom_ivec[m];
    int n = strlen(ename[index]) + 1;
    fwrite(&n,sizeof(int),1,fp);
    fwrite(ename[index],sizeof(char),n,fp);
    fwrite(&etype[index],sizeof(int),1,fp);
    fwrite(&esize[index],sizeof(int),1,fp);
  }
  for (m = 0; m < ncustom_iarray; m++) {
    index = icustom_iarray[m];
    int n = strlen(ename[index]) + 1;
    fwrite(&n,sizeof(int),1,fp);
    fwrite(ename[index],sizeof(char),n,fp);
    fwrite(&etype[index],sizeof(int),1,fp);
    fwrite(&esize[index],sizeof(int),1,fp);
  }
  for (m = 0; m < ncustom_dvec; m++) {
    index = icustom_dvec[m];
    int n = strlen(ename[index]) + 1;
    fwrite(&n,sizeof(int),1,fp);
    fwrite(ename[index],sizeof(char),n,fp);
    fwrite(&etype[index],sizeof(int),1,fp);
    fwrite(&esize[index],sizeof(int),1,fp);
  }
  for (m = 0; m < ncustom_darray; m++) {
    index = icustom_darray[m];
    int n = strlen(ename[index]) + 1;
    fwrite(&n,sizeof(int),1,fp);
    fwrite(ename[index],sizeof(char),n,fp);
    fwrite(&etype[index],sizeof(int),1,fp);
    fwrite(&esize[index],sizeof(int),1,fp);
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads custom attribute definition info from restart file
   bcast to other procs and all procs instantiate series of Mixtures
------------------------------------------------------------------------- */

void Particle::read_restart_custom(FILE *fp)
{
  // ncustom is 0 at time restart file is read
  // will be incremented as add_custom() for each nactive

  int nactive;
  if (me == 0) fread(&nactive,sizeof(int),1,fp);
  MPI_Bcast(&nactive,1,MPI_INT,0,world);
  if (nactive == 0) return;

  // order that custom vectors/arrays are in restart file
  //   matches order the per-particle custom values will be read from file

  int n,type,size;
  char *name;

  for (int i = 0; i < nactive; i++) {
    if (me == 0) fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    name = new char[n];
    if (me == 0) fread(name,sizeof(char),n,fp);
    MPI_Bcast(name,n,MPI_CHAR,0,world);
    if (me == 0) fread(&type,sizeof(int),1,fp);
    MPI_Bcast(&type,n,MPI_CHAR,0,world);
    if (me == 0) fread(&size,sizeof(int),1,fp);
    MPI_Bcast(&size,n,MPI_CHAR,0,world);

    // create the custom attribute

    add_custom(name,type,size);
    delete [] name;
  }

  // set flag for each newly created custom attribute to 0
  // will be reset to 1 if restart script redefines attribute with same name

  custom_restart_flag = new int[ncustom];
  for (int i = 0; i < ncustom; i++) custom_restart_flag[i] = 0;
}

/* ----------------------------------------------------------------------
   pack a custom attributes for a single particle N into buf
   this is done in order of 4 styles of vectors/arrays, not in ncustom order
------------------------------------------------------------------------- */

void Particle::pack_custom(int n, char *buf)
{
  int i;
  char *ptr = buf;

  if (ncustom_ivec) {
    for (i = 0; i < ncustom_ivec; i++) {
      memcpy(ptr,&eivec[i][n],sizeof(int));
      ptr += sizeof(int);
    }
  }
  if (ncustom_iarray) {
    for (i = 0; i < ncustom_iarray; i++) {
      memcpy(ptr,eiarray[i][n],eicol[i]*sizeof(int));
      ptr += eicol[i]*sizeof(int);
    }
  }

  ptr = ROUNDUP(ptr);

  if (ncustom_dvec) {
    for (i = 0; i < ncustom_dvec; i++) {
      memcpy(ptr,&edvec[i][n],sizeof(double));
      ptr += sizeof(double);
    }
  }
  if (ncustom_darray) {
    for (i = 0; i < ncustom_darray; i++) {
      memcpy(ptr,edarray[i][n],edcol[i]*sizeof(double));
      ptr += edcol[i]*sizeof(double);
    }
  }
}

/* ----------------------------------------------------------------------
   unpack custom attributes for a single particle N from buf
   this is done in order of 4 styles of vectors/arrays, not in ncustom order
------------------------------------------------------------------------- */

void Particle::unpack_custom(char *buf, int n)
{
  int i;
  char *ptr = buf;

  if (ncustom_ivec) {
    for (i = 0; i < ncustom_ivec; i++) {
      memcpy(&eivec[i][n],ptr,sizeof(int));
      ptr += sizeof(int);
    }
  }
  if (ncustom_iarray) {
    for (i = 0; i < ncustom_iarray; i++) {
      memcpy(eiarray[i][n],ptr,eicol[i]*sizeof(int));
      ptr += eicol[i]*sizeof(int);
    }
  }

  ptr = ROUNDUP(ptr);

  if (ncustom_dvec) {
    for (i = 0; i < ncustom_dvec; i++) {
      memcpy(&edvec[i][n],ptr,sizeof(double));
      ptr += sizeof(double);
    }
  }
  if (ncustom_darray) {
    for (i = 0; i < ncustom_darray; i++) {
      memcpy(edarray[i][n],ptr,edcol[i]*sizeof(double));
      ptr += edcol[i]*sizeof(double);
    }
  }
}

/* ---------------------------------------------------------------------- */

bigint Particle::memory_usage()
{
  bigint bytes = (bigint) maxlocal * sizeof(OnePart);
  bytes += (bigint) maxlocal * sizeof(int);
  for (int i = 0; i < ncustom_ivec; i++)
    bytes += (bigint) maxlocal * sizeof(int);
  for (int i = 0; i < ncustom_iarray; i++)
    bytes += (bigint) maxlocal*eicol[i] * sizeof(int);
  for (int i = 0; i < ncustom_dvec; i++)
    bytes += (bigint) maxlocal * sizeof(double);
  for (int i = 0; i < ncustom_darray; i++)
    bytes += (bigint) maxlocal*edcol[i] * sizeof(double);
  return bytes;
}

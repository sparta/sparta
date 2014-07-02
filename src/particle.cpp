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

enum{PKEEP,PINSERT,PDONE,PDISCARD,PENTRY,PEXIT};  // several files
enum{NONE,DISCRETE,SMOOTH};            // several files

#define DELTA 10000
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

  exist = 0;
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
   called from Comm::migrate_particles() when particles migrate every step
------------------------------------------------------------------------- */

void Particle::compress(int nmigrate, int *mlist)
{
  int j,k;
  int nbytes = sizeof(OnePart);

  for (int i = 0; i < nmigrate; i++) {
    j = mlist[i];
    k = nlocal - 1;
    while (k == mlist[nmigrate-1] && k > j) {
      nmigrate--;
      nlocal--;
      k--;
    }
    memcpy(&particles[j],&particles[k],nbytes);
    nlocal--;
  }
}

/* ----------------------------------------------------------------------
   compress particle list to remove particles with icell < 0
   all particles MUST be in owned cells
   overwrite deleted particle with particle from end of nlocal list
   called from Comm::migrate_cells() when cells+particles migrate on rebalance
------------------------------------------------------------------------- */

void Particle::compress()
{
  int nbytes = sizeof(OnePart);

  int i = 0;
  while (i < nlocal) {
    if (particles[i].icell < 0) {
      memcpy(&particles[i],&particles[nlocal-1],nbytes);
      nlocal--;
    } else i++;
  }
}

/* ----------------------------------------------------------------------
   sort particles into grid cells
------------------------------------------------------------------------- */

void Particle::sort()
{
  // reallocate next list as needed

  if (maxsort < maxlocal) {
    maxsort = maxlocal;
    memory->destroy(next);
    memory->create(next,maxsort,"particle:next");
  }

  // initialize linked list of particles in cells I own

  Grid::ChildCell *cells = grid->cells;
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

  nlocal++;
  return reallocflag;
}

/* ----------------------------------------------------------------------
   add one or more species to species list
------------------------------------------------------------------------- */

void Particle::add_species(int narg, char **arg)
{
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

  int j;

  for (int i = 0; i < newspecies; i++) {
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

double Particle::erot(int isp, RanPark *erandom)
{
 double eng,a,erm,b;

 if (!collide || collide->rotstyle == NONE) return 0.0;
 if (species[isp].rotdof < 2) return 0.0;

 // NOTE: is temp_thermal always set?

 if (species[isp].rotdof == 2)
   eng = -log(erandom->uniform()) * update->boltz * update->temp_thermal;
 else {
   a = 0.5*particle->species[isp].rotdof-1.;
   while (1) {
     // energy cut-off at 10 kT
     erm = 10.0*erandom->uniform();
     b = pow(erm/a,a) * exp(a-erm);
     if (b > erandom->uniform()) break;
   }
   eng = erm * update->boltz * update->temp_thermal;
 }

 return eng;
}

/* ----------------------------------------------------------------------
   generate random vibrational energy for a particle
   only a function of species index and species properties
------------------------------------------------------------------------- */

double Particle::evib(int isp, RanPark *erandom)
{
  double eng,a,erm,b;

  int vibstyle = NONE;
  if (collide) vibstyle = collide->vibstyle;

  // NOTE: is temp_thermal always set?
  
  eng = 0.0;
  if (vibstyle == DISCRETE && species[isp].vibdof == 2) {
    int ivib = -log(erandom->uniform()) * 
      update->temp_thermal / particle->species[isp].vibtemp;
    eng = ivib * update->boltz * particle->species[isp].vibtemp;
  } else if (vibstyle == SMOOTH || species[isp].vibdof >= 2) {
    if (species[isp].vibdof == 2)
      eng = -log(erandom->uniform()) * update->boltz * update->temp_thermal;
    else {
      a = 0.5*particle->species[isp].vibdof-1.;
      while (1) {
        // energy cut-off at 10 kT
        erm = 10.0*erandom->uniform();
        b = pow(erm/a,a) * exp(a-erm);
        if (b > erandom->uniform()) break;
      }
      eng = erm * update->boltz * update->temp_thermal;
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
    fsp->rotrel = atoi(words[4]);
    fsp->vibdof = atoi(words[5]);
    fsp->vibrel = atoi(words[6]);
    fsp->vibtemp = atof(words[7]);
    fsp->specwt = atof(words[8]);
    fsp->charge = atof(words[9]);

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
   // NOTE: worry about N overflowing int, and in IROUNDUP ???
------------------------------------------------------------------------- */

int Particle::size_restart()
{
  int n = sizeof(int);
  n = IROUNDUP(n);
  n += nlocal * sizeof(OnePartRestart);
  n = IROUNDUP(n);
  return n;
}

/* ----------------------------------------------------------------------
   pack my particle info into buf
   use OnePartRestart data struct for permanent info and to encode cell ID
   // NOTE: worry about N overflowing int, and in IROUNDUP ???
------------------------------------------------------------------------- */

int Particle::pack_restart(char *buf)
{
  int n,nsplit;

  Grid::ChildCell *cells = grid->cells;
  OnePart *p;
  OnePartRestart *pr;

  int *ibuf = (int *) buf;
  *ibuf = nlocal;
  n = sizeof(int);
  n = IROUNDUP(n);

  OnePartRestart *pbuf = (OnePartRestart *) &buf[n];

  for (int i = 0; i < nlocal; i++) {
    p = &particles[i];
    pr = &pbuf[i];
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
  }

  n += nlocal * sizeof(OnePartRestart);
  n = IROUNDUP(n);

  return n;
}

/* ----------------------------------------------------------------------
   unpack particle info into restart storage
   allocate data structure here, will be deallocated by ReadRestart
------------------------------------------------------------------------- */

int Particle::unpack_restart(char *buf)
{
  int n;
  Grid::ChildCell *cells = grid->cells;
  OnePartRestart *p;

  int *ibuf = (int *) buf;
  nlocal_restart = *ibuf;
  n = sizeof(int);
  n = IROUNDUP(n);

  particle_restart = (OnePartRestart *) 
    memory->smalloc(nlocal_restart*sizeof(OnePartRestart),
                    "grid:particle_restart");

  memcpy(particle_restart,&buf[n],nlocal_restart*sizeof(OnePartRestart));
  n += nlocal_restart * sizeof(OnePartRestart);
  n = IROUNDUP(n);

  return n;
}

/* ---------------------------------------------------------------------- */

bigint Particle::memory_usage()
{
  bigint bytes = (bigint) maxlocal * sizeof(OnePart);
  bytes += (bigint) maxlocal * sizeof(int);
  return bytes;
}

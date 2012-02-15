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

#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "particle.h"
#include "grid.h"
#include "update.h"
#include "comm.h"
#include "mixture.h"
#include "memory.h"
#include "error.h"

using namespace DSMC_NS;

#define DELTA 10000
#define DELTASMALL 16
#define MAXLINE 1024

// customize by adding an abbreviation string
// also add a check for the keyword in 2 places in add_species()

#define AIR "O N NO"

/* ---------------------------------------------------------------------- */

Particle::Particle(DSMC *dsmc) : Pointers(dsmc)
{
  nglobal = 0;
  nlocal = maxlocal = 0;
  particles = NULL;

  nspecies = maxspecies  = 0;
  species = NULL;

  nmixture = 0;
  mixture = NULL;

  maxgrid = 0;
  //cellcount = NULL;
  //first = NULL;
  maxsort = 0;
  next = NULL;
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
}

/* ---------------------------------------------------------------------- */

void Particle::init()
{
  for (int i = 0; i < nmixture; i++) mixture[i]->init();

  // reallocate cellcount and first lists as needed
  // NOTE: when grid becomes dynamic, will need to do this in sort()

  if (maxgrid < grid->nlocal) {
    maxgrid = grid->nlocal;
    //    memory->destroy(cellcount);
    //memory->destroy(first);
    //memory->create(first,maxgrid,"particle:first");
    //memory->create(cellcount,maxgrid,"particle:cellcount");
  }
}

/* ----------------------------------------------------------------------
   compress particle list to remove mlist of migrating particles
   overwrite deleted particle with particle from end of nlocal list
   j = mlist loop avoids overwrite with deleted particle at end of mlist
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
   sort particles into grid cells
------------------------------------------------------------------------- */

void Particle::sort()
{
  int i,icell;

  // reallocate next list as needed

  if (maxsort < maxlocal) {
    maxsort = maxlocal;
    memory->destroy(next);
    memory->create(next,maxsort,"particle:next");
  }

  // initialize linked list of particles in cells I own

  Grid::OneCell *cells = grid->cells;
  int *mycells = grid->mycells;
  int nglocal = grid->nlocal;

  for (i = 0; i < nglocal; i++) {
    icell = mycells[i];
    cells[icell].first = -1;
    cells[icell].count = 0;

    //cellcount[i] = 0;
    //first[i] = -1;
  }

  // reverse loop stores linked list in forward order
  // icell = global cell the particle is in

  for (i = nlocal-1; i >= 0; i--) {
    icell = particles[i].icell;
    next[i] = cells[icell].first;
    cells[icell].first = i;
    cells[icell].count++;

    // NOTE: this method seems much slower for some reason
    // uses separate, smaller vectors for first & cellcount
    //icell = cells[particles[i].icell].local;
    //next[i] = first[icell];
    //first[icell] = i;
    //cellcount[icell]++;
  }
}

/* ----------------------------------------------------------------------
   insure particle list can hold nextra new particles
------------------------------------------------------------------------- */

void Particle::grow(int nextra)
{
  bigint target = (bigint) nlocal + nextra;
  if (target <= maxlocal) return;
  
  bigint newmax = maxlocal;
  while (newmax < target) newmax += DELTA;
  
  if (newmax > MAXSMALLINT) 
    error->one(FLERR,"Per-processor grid count is too big");

  maxlocal = newmax;
  particles = (OnePart *)
    memory->srealloc(particles,maxlocal*sizeof(OnePart),
		     "particle:particles");
}

/* ----------------------------------------------------------------------
   add a particle to particle list
------------------------------------------------------------------------- */

void Particle::add_particle(int id, int ispecies, int icell,
			    double *x, double *v)
{
  if (nlocal == maxlocal) grow(1);

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

  nlocal++;
}

/* ----------------------------------------------------------------------
   add one of more species to species list
------------------------------------------------------------------------- */

void Particle::add_species(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal species command");

  if (comm->me == 0) {
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

  if (comm->me == 0) read_species_file();
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
    while (nspecies+newspecies > maxspecies) maxspecies += DELTASMALL;
    species = (Species *) 
      memory->srealloc(species,maxspecies*sizeof(Species),"particle:species");
  }

  // extract info on user-requested species from file species list

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
  }

  memory->sfree(filespecies);
  delete [] names;
}

/* ----------------------------------------------------------------------
   add or augment a mixture of species
------------------------------------------------------------------------- */

void Particle::add_mixture(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal mixture command");

  // imix = index if mixture ID already exists
  // else instantiate a new mixture

  int imix = find_mixture(arg[0]);

  if (imix < 0) {
    mixture = (Mixture **) memory->srealloc(mixture,
					    (nmixture+1)*sizeof(Mixture *),
					    "particle:mixture");
    imix = nmixture;
    mixture[nmixture++] = new Mixture(dsmc,arg[0]);
  }

  // nsp = # of species args before optional keywords
  // iarg = start of optional keywords

  int iarg = 1;
  for (iarg = 2; iarg < narg; iarg++) {
    if (strcmp(arg[iarg],"nrho") == 0) break;
    if (strcmp(arg[iarg],"vstream") == 0) break;
    if (strcmp(arg[iarg],"temp") == 0) break;
    if (strcmp(arg[iarg],"frac") == 0) break;
    if (strcmp(arg[iarg],"group") == 0) break;
  }
  int nsp = iarg - 1;

  // pass list of species and list of params to Mixture

  mixture[imix]->add_species(nsp,&arg[1]);
  mixture[imix]->params(narg-iarg,&arg[iarg]);
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

  int NWORDS = 14;
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
      maxfilespecies += DELTASMALL;
      filespecies = (Species *) 
	memory->srealloc(filespecies,maxfilespecies*sizeof(Species),
			 "particle:filespecies");
    }

    nwords = wordcount(line,words);
    Species *fsp = &filespecies[nfilespecies];

    if (strlen(words[0]) + 1 > 16) error->one(FLERR,"");
    strcpy(fsp->id,words[0]);

    fsp->molwt = atof(words[1]);
    fsp->mass = atof(words[2]);
    fsp->diam = atof(words[3]);
    fsp->rotdof = atoi(words[4]);
    fsp->rotrel = atoi(words[5]);
    fsp->vibdof = atoi(words[6]);
    fsp->vibrel = atoi(words[7]);
    fsp->vibtemp = atof(words[8]);
    fsp->specwt = atof(words[9]);
    fsp->charge = atof(words[10]);
    fsp->omega = atof(words[11]);
    fsp->tref = atof(words[12]);
    fsp->alpha = atof(words[13]);

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

/* ---------------------------------------------------------------------- */

bigint Particle::memory_usage()
{
  bigint bytes = (bigint) maxlocal * sizeof(OnePart);
  bytes += (bigint) maxlocal * sizeof(int);
  return bytes;
}

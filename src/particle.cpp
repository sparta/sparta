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

#include "string.h"
#include "stdlib.h"
#include "particle.h"
#include "domain.h"
#include "grid.h"
#include "update.h"
#include "comm.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace DSMC_NS;

#define DELTA 10000
#define DELTASPECIES 16
#define MAXLINE 1024

// customize by adding an abbreviation string
// also add a check for the keyword in 2 places in add_species()

#define AIR "O N NO"

enum{ALL,LOCAL};

/* ---------------------------------------------------------------------- */

Particle::Particle(DSMC *dsmc) : Pointers(dsmc)
{
  nglobal = 0;
  nlocal = maxlocal = 0;
  particles = NULL;

  nspecies = maxspecies = 0;
  species = NULL;
}

/* ---------------------------------------------------------------------- */

Particle::~Particle()
{
  memory->sfree(particles);
  memory->sfree(species);
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
      maxfilespecies += DELTASPECIES;
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
   create particles
   called from input script
------------------------------------------------------------------------- */

void Particle::create(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal create_particles command");

  bigint n = ATOBIGINT(arg[0]);
  if (n < 0) error->all(FLERR,"Illegal create_particles command");
  seed = atoi(arg[1]);
  if (seed <= 0) error->all(FLERR,"Illegal create_particles command");

  // optional args

  int loop = ALL;

  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"loop") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_particles command");
      if (strcmp(arg[iarg+1],"all") == 0) loop = ALL;
      else if (strcmp(arg[iarg+1],"local") == 0) loop = LOCAL;
      else error->all(FLERR,"Illegal create_particles command");
      iarg += 2;
    } else error->all(FLERR,"Illegal create_particles command");
  }

  // generate particles

  bigint nprevious = nglobal;

  if (loop == ALL) create_all(n);
  else if (loop == LOCAL) create_local(n);

  // error check

  bigint nme = nlocal;
  MPI_Allreduce(&nme,&nglobal,1,MPI_DSMC_BIGINT,MPI_SUM,world);
  if (nglobal - nprevious != n) {
    char str[128];
    sprintf(str,"Created incorrect # of particles = " BIGINT_FORMAT,
	    nglobal-nprevious);
    error->all(FLERR,str);
  }

  // print stats

  if (comm->me == 0) {
    if (screen) fprintf(screen,"Created " BIGINT_FORMAT " particles\n",n);
    if (logfile) fprintf(logfile,"Created " BIGINT_FORMAT " particles\n",n);
  }
}

/* ----------------------------------------------------------------------
   create N particles in serial
   every proc generates all N coords, only keeps those in cells it owns
------------------------------------------------------------------------- */

void Particle::create_all(bigint n)
{
  int dimension = domain->dimension;
  double xlo = domain->boxlo[0];
  double ylo = domain->boxlo[1];
  double zlo = domain->boxlo[2];
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  int me = comm->me;
  RanPark *random = new RanPark(dsmc,seed);

  int icell;
  double x,y,z;

  // loop over all N particles

  for (bigint m = 0; m < n; m++) {
    x = xlo + random->uniform()*xprd;
    y = ylo + random->uniform()*yprd;
    z = zlo + random->uniform()*zprd;
    if (dimension == 2) z = 0.0;

    // which_cell() returns global grid cell index the particle is in
    // if I own that grid cell, store particle

    icell = grid->which_cell(x,y,z);

    if (grid->cells[icell].proc == me) {
      if (nlocal == maxlocal) grow(1);
      particles[nlocal].id = 0;
      particles[nlocal].type = 1;
      particles[nlocal].icell = icell;
      particles[nlocal].x[0] = x;
      particles[nlocal].x[1] = y;
      particles[nlocal].x[2] = z;
      particles[nlocal].v[0] = 0.0;
      particles[nlocal].v[1] = 0.0;
      particles[nlocal].v[2] = 0.0;
      nlocal++;
    }
  }

  delete random;
}

/* ----------------------------------------------------------------------
   create N particles in serial
   every proc generates all N coords, only keeps those in cells it owns
------------------------------------------------------------------------- */

void Particle::create_local(bigint n)
{
  int dimension = domain->dimension;

  int me = comm->me;
  RanPark *random = new RanPark(dsmc,seed+me);

  Grid::OneCell *cells = grid->cells;
  int nglocal = grid->nlocal;

  bigint nme = n/comm->nprocs;
  if (me < n % comm->nprocs) nme++;

  // list = list of indices of grid cells I own

  int *list = new int[nglocal];
  int ncell = grid->ncell;

  int j = 0;
  for (int i = 0; i < ncell; i++)
    if (cells[i].proc == me) list[j++] = i;

  // loop only over Nme particles I own

  int icell;
  double x,y,z;
  double *lo,*hi;

  for (bigint m = 0; m < nme; m++) {
    j = static_cast<int> (random->uniform()*nglocal);
    icell = list[j];
    lo = cells[icell].lo;
    hi = cells[icell].hi;

    x = lo[0] + random->uniform() * (hi[0]-lo[0]);
    y = lo[1] + random->uniform() * (hi[1]-lo[1]);
    z = lo[2] + random->uniform() * (hi[2]-lo[2]);
    if (dimension == 2) z = 0.0;

    if (nlocal == maxlocal) grow(1);
    particles[nlocal].id = 0;
    particles[nlocal].type = 1;
    particles[nlocal].icell = icell;
    particles[nlocal].x[0] = x;
    particles[nlocal].x[1] = y;
    particles[nlocal].x[2] = z;
    particles[nlocal].v[0] = 0.0;
    particles[nlocal].v[1] = 0.0;
    particles[nlocal].v[2] = 0.0;
    nlocal++;
  }

  delete random;
  delete [] list;
}

/* ----------------------------------------------------------------------
   compress particle list to remove mlist of migrating particles
   overwrite deleted particle with particle from end of nlocal list
   j = mlist loop avoids overwrite with deleted particle at end of mlist
------------------------------------------------------------------------- */

void Particle::compress()
{
  int j,k;
  int nbytes = sizeof(OnePart);
  int *mlist = update->mlist;
  int nmigrate = update->nmigrate;

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

  // reallocate sort list as needed

  if (maxsortparticle < maxlocal) {
    maxsortparticle = maxlocal;
    memory->destroy(next);
    memory->create(next,maxsortparticle,"sort:next");
  }

  // build linked list for each cell

  Grid::OneCell *cells = grid->cells;
  int ncelllocal = grid->nlocal;

  for (icell = 0; icell < ncelllocal; icell++) {
    cells[icell].nparticles = 0;
    cells[icell].first = -1;
  }

  for (i = 0; i < nlocal; i++) {
    icell = particles[i].icell;
    if (cells[icell].first < 0) {
      cells[icell].first = i;
      next[i] = -1;
    } else {


      next[i] = -1;
    }
    grid->cells[icell].nparticles++;
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

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
#include "fix_vibmode.h"
#include "random_mars.h"
#include "random_knuth.h"
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

  nspecies = maxspecies = 0;
  species = NULL;
  maxvibmode = 0;
  maxelecstate = 0;
  cumulative_probabilities = NULL;

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

  // RNG for particle weighting

  wrandom = NULL;

  copy = uncopy = copymode = 0;
}

/* ---------------------------------------------------------------------- */

Particle::~Particle()
{
  if (!uncopy && (copy || copymode)) return;

  for (int i = 0; i < nspecies; i++) {
    if (species && species[i].elecdat != NULL) {
      memory->sfree(species[i].elecdat->states);
      memory->destroy(species[i].elecdat->default_rel);
      for (int j = 0; j < nspecies; j++) {
        if (species[i].elecdat->species_rel[j] != NULL) {
          memory->destroy(species[i].elecdat->species_rel[j]);
        }
      }
      memory->sfree(species[i].elecdat->species_rel);
      memory->destroy(species[i].elecdat->enforce_spin_conservation);
      delete species[i].elecdat;
    }
  }
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
  memory->destroy(cumulative_probabilities);

  delete wrandom;
}

/* ---------------------------------------------------------------------- */

void Particle::init()
{
  // check for errors in custom particle vectors/arrays

  error_custom();

  // initialize mixtures

  for (int i = 0; i < nmixture; i++) mixture[i]->init();

  // RNG for particle weighting

  if (!wrandom) {
    wrandom = new RanKnuth(update->ranmaster->uniform());
    double seed = update->ranmaster->uniform();
    wrandom->reset(seed,me,100);
  }

  // if vibstyle = DISCRETE,
  // all species with vibdof > 2 must have info read from a species.vib file

  if (collide && collide->vibstyle == DISCRETE) {
    for (int isp = 0; isp < nspecies; isp++) {
      if (species[isp].vibdof <= 2) continue;
      if (species[isp].vibdiscrete_read == 0) {
        char str[128];
        sprintf(str,"Discrete vibrational info for species %s not read in",
                species[isp].id);
        error->all(FLERR,str);
      }
    }
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
   this does NOT preserve particle sorting
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
   called from AdaptGrid when coarsening occurs
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
   compress particle list to remove particles with icell < 0
   same as compress_rebalance() except this DOES preserve particle sorting
   invoked by balance migrate_cells_less_memory()
------------------------------------------------------------------------- */

void Particle::compress_rebalance_sorted()
{
  int nbytes = sizeof(OnePart);

  Grid::ChildInfo *cinfo = grid->cinfo;

  if (!ncustom) {
    int i = 0;
    while (i < nlocal) {
      if (particles[i].icell < 0) {
        int icell = particles[nlocal-1].icell;
        if (icell >= 0) {
          if (cinfo[icell].first == nlocal-1) cinfo[icell].first = i;
          else {
            int ip = cinfo[icell].first;
            while (ip >= 0) {
              if (next[ip] == nlocal-1) {
                next[ip] = i;
                break;
              }
              ip = next[ip];
            }
          }
          next[i] = next[nlocal-1];
        }
        memcpy(&particles[i],&particles[nlocal-1],nbytes);
        nlocal--;
      } else i++;
    }

  } else {
    int i = 0;
    while (i < nlocal) {
      if (particles[i].icell < 0) {
        int icell = particles[nlocal-1].icell;
        if (icell >= 0) {
          if (cinfo[icell].first == nlocal-1) cinfo[icell].first = i;
          else {
            int ip = cinfo[icell].first;
            while (ip >= 0) {
              if (next[ip] == nlocal-1) {
                next[ip] = i;
                break;
              }
              ip = next[ip];
            }
          }
          next[i] = next[nlocal-1];
        }
        memcpy(&particles[i],&particles[nlocal-1],nbytes);
        copy_custom(i,nlocal-1);
        nlocal--;
      } else i++;
    }
  }
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

  // reallocate next list as needed

  if (maxsort < maxlocal) {
    maxsort = maxlocal;
    memory->destroy(next);
    memory->create(next,maxsort,"particle:next");
  }

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
  // NOTE: why not just compare maxsort to nlocal?
  //       then could realloc less often?
  //       ditto for all reallocs of next in related methods

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
  }

  // reverse loop over partlcles to store linked lists in forward order
  // icell = global cell the particle is in

  int icell;
  for (int i = nlocal-1; i >= 0; i--) {
    icell = particles[i].icell;
    next[i] = cinfo[icell].first;
    cinfo[icell].first = i;
    cinfo[icell].count++;
  }
}

/* ----------------------------------------------------------------------
   reallocate next list if necessary
   called before partial sort by FixEmit classes in subsonic case
------------------------------------------------------------------------- */

void Particle::sort_allocate()
{
  if (maxsort < maxlocal) {
    maxsort = maxlocal;
    memory->destroy(next);
    memory->create(next,maxsort,"particle:next");
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
                     "particle:particles",SPARTA_GET_ALIGN(OnePart));
  memset(&particles[oldmax],0,(maxlocal-oldmax)*sizeof(OnePart));

  if (ncustom == 0) return;

  for (int i = 0; i < ncustom; i++) {
    if (ename[i] == NULL) continue;
    grow_custom(i,oldmax,maxlocal);
  }
}

/* ----------------------------------------------------------------------
   insure species list can hold maxspecies species
   assumes that maxspecies has already been increased
------------------------------------------------------------------------- */

void Particle::grow_species()
{
  species = (Species *)
    memory->srealloc(species,maxspecies*sizeof(Species),"particle:species");

  for (int isp = 0; isp < nspecies; ++isp) {
    if (species[isp].elecdat != NULL) {
      memory->srealloc(species[isp].elecdat->species_rel, maxspecies*sizeof(double*),"elecdat:species_rel");
      memory->grow(species[isp].elecdat->enforce_spin_conservation, maxspecies, "elecdat:enforce_spin_conservation");
    }
  }
}

/* ----------------------------------------------------------------------
   grow next list if more particles now exist than there is room for
   called from Grid::unpack_particles_adapt() when grid adaptation
     takes place and acquire particles from other procs due to coarsening
   unlike sort(), this requires next list be grown, not destroy/create
     b/c sorted particle list is maintained during adaptation
------------------------------------------------------------------------- */

void Particle::grow_next()
{
  // compare maxsort (length of next) to new particle count (nlocal)
  // grow to maxlocal (max length of particles) to avoid frequent re-grow

  if (maxsort < nlocal) {
    maxsort = maxlocal;
    memory->grow(next,maxsort,"particle:next");
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

  if (ncustom) zero_custom(nlocal);

  nlocal++;
  return reallocflag;
}

/* ----------------------------------------------------------------------
   add an empty particle to particle list, caller will fill it
   return 1 if particle array was reallocated, else 0
------------------------------------------------------------------------- */

int Particle::add_particle()
{
  int reallocflag = 0;
  if (nlocal == maxlocal) {
    grow(1);
    reallocflag = 1;
  }

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
  int i,j,k,n;;

  if (narg < 2) error->all(FLERR,"Illegal species command");

  if (me == 0) {
    fp = fopen(arg[0],"r");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open species file %s",arg[0]);
      error->one(FLERR,str);
    }
  }

  // nfile = # of species defined in file
  // filespecies = list of species defined in file

  nfile = maxfile = 0;
  filespecies = NULL;

  if (me == 0) read_species_file();
  MPI_Bcast(&nfile,1,MPI_INT,0,world);
  if (comm->me) {
    filespecies = (Species *)
      memory->smalloc(nfile*sizeof(Species),"particle:filespecies");
  }
  MPI_Bcast(filespecies,nfile*sizeof(Species),MPI_BYTE,0,world);

  // newspecies = # of new user-requested species,
  //   not including trailing optional keywords
  // names = list of new species IDs
  // customize abbreviations by adding new keyword in 2 places

  char line[MAXLINE];

  int newspecies = 0;
  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"air") == 0) {
      strcpy(line,AIR);
      newspecies += wordcount(line,NULL);
    } else if (strcmp(arg[iarg],"rotfile") == 0) {
      break;
    } else if (strcmp(arg[iarg],"vibfile") == 0) {
      break;
    } else {
      newspecies++;
    }
    iarg++;
  }

  char **names = new char*[newspecies];
  newspecies = 0;

  iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"air") == 0) {
      strcpy(line,AIR);
      newspecies += wordcount(line,&names[newspecies]);
    } else if (strcmp(arg[iarg],"rotfile") == 0) {
      break;
    } else if (strcmp(arg[iarg],"vibfile") == 0) {
      break;
    } else if (strcmp(arg[iarg],"elecfile") == 0) {
      break;
    } else {
      names[newspecies++] = arg[iarg];
    }
    iarg++;
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
    grow_species();
  }

  // extract info on user-requested species from file species list
  // add new species to default mixtures "all" and "species"

  int nspecies_original = nspecies;

  int imix_all = find_mixture((char *) "all");
  int imix_species = find_mixture((char *) "species");

  for (i = 0; i < newspecies; i++) {
    for (j = 0; j < nspecies; j++)
      if (strcmp(names[i],species[j].id) == 0) break;
    if (j < nspecies) error->all(FLERR,"Species ID is already defined");
    for (j = 0; j < nfile; j++)
      if (strcmp(names[i],filespecies[j].id) == 0) break;
    if (j == nfile)
      error->all(FLERR,"Species ID does not appear in species file");
    memcpy(&species[nspecies],&filespecies[j],sizeof(Species));
    nspecies++;

    mixture[imix_all]->add_species_default(species[nspecies-1].id);
    mixture[imix_species]->add_species_default(species[nspecies-1].id);
  }

  memory->sfree(filespecies);

  // process any optional keywords
  // NOTE: rotfile is not yet supported

  int rotindex = 0;
  int vibindex = 0;
  int elecindex = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"rotfile") == 0) {
      // not yet supported
      error->all(FLERR,"Illegal species command");
      if (iarg+2 > narg) error->all(FLERR,"Illegal species command");
      if (rotindex)
        error->all(FLERR,"Species command can only use a single rotfile");
      rotindex = iarg+1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"vibfile") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal species command");
      if (vibindex)
        error->all(FLERR,"Species command can only use a single vibfile");
      vibindex = iarg+1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"elecfile") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal species command");
      if (elecindex)
        error->all(FLERR,"Species command can only use a single elecfile");
      elecindex = iarg+1;
      iarg += 2;
    } else error->all(FLERR,"Illegal species command");
  }

  // read rotational species file and setup per-species params

  if (rotindex) {
    if (me == 0) {
      fp = fopen(arg[rotindex],"r");
      if (fp == NULL) {
        char str[128];
        sprintf(str,"Cannot open rotation file %s",arg[rotindex]);
        error->one(FLERR,str);
      }
    }

    nfile = maxfile = 0;
    filerot = NULL;

    if (me == 0) read_rotation_file();
    MPI_Bcast(&nfile,1,MPI_INT,0,world);
    if (comm->me) {
      filerot = (RotFile *)
        memory->smalloc(nfile*sizeof(RotFile),"particle:filerot");
    }
    MPI_Bcast(filerot,nfile*sizeof(RotFile),MPI_BYTE,0,world);

    for (i = 0; i < newspecies; i++) {
      int ii = nspecies_original + i;

      for (j = 0; j < nfile; j++)
        if (strcmp(names[i],filerot[j].id) == 0) break;
      if (j == nfile) {
        if (species[ii].rotdof == 0) continue;
        error->all(FLERR,"Species ID does not appear in rotation file");
      }

      int ntemp = filerot[j].ntemp;
      if ((species[ii].rotdof == 0) ||
          (species[ii].rotdof == 2 && ntemp != 1) ||
          (species[ii].rotdof == 3 && ntemp != 3))
          error->all(FLERR,"Mismatch between species rotdof "
                     "and rotation file entry");

      species[ii].nrottemp = ntemp;
      for (k = 0; k < ntemp; k++) {
        species[ii].rottemp[k] = filerot[j].rottemp[k];
      }
    }

    memory->sfree(filerot);
  }

  // read vibrational species file and setup per-species params

  if (vibindex) {
    if (me == 0) {
      fp = fopen(arg[vibindex],"r");
      if (fp == NULL) {
        char str[128];
        sprintf(str,"Cannot open vibration file %s",arg[vibindex]);
        error->one(FLERR,str);
      }
    }

    nfile = maxfile = 0;
    filevib = NULL;

    if (me == 0) read_vibration_file();
    MPI_Bcast(&nfile,1,MPI_INT,0,world);
    if (comm->me) {
      filevib = (VibFile *)
        memory->smalloc(nfile*sizeof(VibFile),"particle:filevib");
    }
    MPI_Bcast(filevib,nfile*sizeof(VibFile),MPI_BYTE,0,world);

    for (i = 0; i < newspecies; i++) {
      int ii = nspecies_original + i;

      for (j = 0; j < nfile; j++)
        if (strcmp(names[i],filevib[j].id) == 0) break;
      if (j == nfile) {
        if (species[ii].vibdof <= 2) continue;
        error->all(FLERR,"Species ID does not appear in vibration file");
      }

      int nmode = filevib[j].nmode;
      if (species[ii].nvibmode != nmode)
        error->all(FLERR,"Mismatch between species vibdof "
                   "and vibration file entry");

      species[ii].nvibmode = nmode;
      for (k = 0; k < nmode; k++) {
        species[ii].vibtemp[k] = filevib[j].vibtemp[k];
        species[ii].vibrel[k] = filevib[j].vibrel[k];
        species[ii].vibdegen[k] = filevib[j].vibdegen[k];
      }

      maxvibmode = MAX(maxvibmode,species[ii].nvibmode);
      species[ii].vibdiscrete_read = 1;
    }

    memory->sfree(filevib);
  }

  if (elecindex) {
    if (me == 0) {
      fp = fopen(arg[elecindex],"r");
      if (fp == NULL) {
        char str[128];
        sprintf(str,"Cannot open electronic file %s",arg[elecindex]);
        error->one(FLERR,str);
      }
    }

    nfile = maxfile = 0;
    fileelec = NULL;

    if (me == 0) read_electronic_file();
    MPI_Bcast(&nfile,1,MPI_INT,0,world);

    if (comm->me) {
      fileelec = (ElecFile *)
        memory->smalloc(nfile*sizeof(ElecFile),"particle:fileelec");
    }
    MPI_Bcast(fileelec,nfile*sizeof(ElecFile),MPI_BYTE,0,world);

    if (comm->me) {
      for (i = 0; i < nfile; ++i) {
        fileelec[i].default_rel = (double*) memory->smalloc(
          fileelec[i].nmode*sizeof(double), "vsp:default_rel");
        fileelec[i].elecrel = (double**) memory->smalloc(
          nspecies*sizeof(double*), "vsp:elecrel");
        for (int j = 0; j < nspecies; ++j) {
          fileelec[i].elecrel[j] = (double*) memory->smalloc(
            fileelec[i].nmode*sizeof(double), "vsp:elecrel[]");
        }
        fileelec[i].enforce_spin_conservation = (bool*) memory->smalloc(
          nspecies*sizeof(bool), "vsp:enforce_spin_conservation");
        fileelec[i].electemp = (double*) memory->smalloc(
          fileelec[i].nmode*sizeof(double), "vsp:electemp");
        fileelec[i].elecdegen = (int*) memory->smalloc(
          fileelec[i].nmode*sizeof(int), "vsp:elecdegen");
        fileelec[i].elecspin = (int*) memory->smalloc(
          fileelec[i].nmode*sizeof(int), "vsp:elecspin");
      }
    }

    for (i = 0; i < nfile; ++i) {
      MPI_Bcast(fileelec[i].default_rel,fileelec[i].nmode*sizeof(double),MPI_BYTE,0,world);
      for (int j = 0; j < nspecies; ++j) {
        MPI_Bcast(fileelec[i].elecrel[j],fileelec[i].nmode*sizeof(double),MPI_BYTE,0,world);
      }
      MPI_Bcast(fileelec[i].enforce_spin_conservation,nspecies*sizeof(bool),MPI_BYTE,0,world);
      MPI_Bcast(fileelec[i].electemp,fileelec[i].nmode*sizeof(double),MPI_BYTE,0,world);
      MPI_Bcast(fileelec[i].elecdegen,fileelec[i].nmode*sizeof(int),MPI_BYTE,0,world);
      MPI_Bcast(fileelec[i].elecspin,fileelec[i].nmode*sizeof(int),MPI_BYTE,0,world);
    }

    for (i = 0; i < newspecies; i++) {
      int ii = nspecies_original + i;

      // Allocate memory for the electronic mode only for species with
      // such data defined in the elecfile.
      for (j = 0; j < nfile; j++)
        if (strcmp(names[i],fileelec[j].id) == 0) break;
      if (j == nfile) {
        species[ii].elecdat = NULL;
        continue;
      }

      int nmode = fileelec[j].nmode;
      maxelecstate = MAX(maxelecstate, nmode);
      species[ii].elecdat = new ElectronicData();
      species[ii].elecdat->nelecstate = nmode;
      species[ii].elecdat->states = (ElecState*) memory->smalloc(nmode*sizeof(ElecState),"elecdat:elecstate");
      memory->create(species[ii].elecdat->default_rel, nmode, "elecdat:default_rel");
      species[ii].elecdat->species_rel = (double**) memory->smalloc(maxspecies*sizeof(double*),"elecdat:species_rel");
      memory->create(species[ii].elecdat->enforce_spin_conservation, maxspecies, "elecdat:enforce_spin_conservation");
      for (int isp = 0; isp < particle->nspecies; ++isp) {
        species[ii].elecdat->species_rel[isp] = NULL;
      }
      for (int isp = 0; isp < particle->nspecies; ++isp) {
        species[ii].elecdat->enforce_spin_conservation[isp] = fileelec[j].enforce_spin_conservation[isp];
      }
      for (k = 0; k < nmode; k++) {
        species[ii].elecdat->states[k].temp = fileelec[j].electemp[k];
        species[ii].elecdat->states[k].degen = fileelec[j].elecdegen[k];
        species[ii].elecdat->states[k].spin = fileelec[j].elecspin[k];
        species[ii].elecdat->default_rel[k] = fileelec[j].default_rel[k];

      }
      for (int isp = 0; isp < particle->nspecies; ++isp) {
        if (isp == ii) continue;
        if (fileelec[j].elecrel[isp][0] >= 0) {
          memory->create(species[ii].elecdat->species_rel[isp], nmode, "elecdat:species_rel[]");
          for (k = 0; k < nmode; k++) {
            species[ii].elecdat->species_rel[isp][k] = fileelec[j].elecrel[isp][k];
          }
        }
      }
      species[ii].elecdiscrete_read = 1;
    }
    for (i = 0; i < nfile; ++i) {
      memory->sfree(fileelec[i].default_rel);
      for (int j = 0; j < nspecies; ++j) {
        memory->sfree(fileelec[i].elecrel[j]);
      }
      memory->sfree(fileelec[i].elecrel);
      memory->sfree(fileelec[i].enforce_spin_conservation);
      memory->sfree(fileelec[i].electemp);
      memory->sfree(fileelec[i].elecdegen);
      memory->sfree(fileelec[i].elecspin);
    }
    memory->sfree(fileelec);
  }
  // clean up
  if (cumulative_probabilities) memory->destroy(cumulative_probabilities);
  if (maxelecstate)
    memory->create(cumulative_probabilities, maxelecstate+1, "particle:cumulative_probabilities");

  delete [] names;
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
   return index of ID in list of species IDs
   return -1 if not found
------------------------------------------------------------------------- */

void Particle::species_modify(int narg, char **arg)
{
  if (narg < 3) error->all(FLERR,"Illegal species_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (iarg+3 > narg) error->all(FLERR,"Illegal species_modify command");

    int ispecies = find_species(arg[iarg]);
    if (ispecies < 0) error->all(FLERR,"Species_modify species does not exist");

    if (strcmp(arg[iarg+1],"mu") == 0)
      species[ispecies].magmoment = atof(arg[iarg+2]);
    else error->all(FLERR,"Unrecognized species_modify property");

    iarg += 3;
  }
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

double Particle::erot(int isp, double temp_thermal, RanKnuth *erandom)
{
  double eng,a,erm,b;
  int rotstyle = NONE;
  if (collide) rotstyle = collide->rotstyle;

  if (!collide || collide->rotstyle == NONE) return 0.0;
  if (species[isp].rotdof < 2) return 0.0;

  if (rotstyle == DISCRETE && species[isp].rotdof == 2) {
    int irot = -log(erandom->uniform()) * temp_thermal /
      particle->species[isp].rottemp[0];
    eng = irot * update->boltz * particle->species[isp].rottemp[0];
  } else if (rotstyle == SMOOTH && species[isp].rotdof == 2) {
    eng = -log(erandom->uniform()) * update->boltz * temp_thermal;
  } else {
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
   index_vibmode = index of extra per-particle vibrational mode storage
     -1 if not defined for this model
------------------------------------------------------------------------- */

double Particle::evib(int isp, double temp_thermal, RanKnuth *erandom)
{
  double eng,a,erm,b;

  int vibstyle = NONE;
  if (collide) vibstyle = collide->vibstyle;
  if (vibstyle == NONE || species[isp].vibdof < 2) return 0.0;

  // for DISCRETE, only need set evib for vibdof = 2
  // mode levels and evib will be set by FixVibmode::update_custom()

  eng = 0.0;

  if (vibstyle == DISCRETE && species[isp].vibdof == 2) {
    int ivib = -log(erandom->uniform()) * temp_thermal /
      particle->species[isp].vibtemp[0];
    eng = ivib * update->boltz * particle->species[isp].vibtemp[0];
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

/* ---------------------------------------------------------------------- */

int Particle::ielec(int isp, double temp_elec, RanKnuth *erandom)
{
  int ielec = 0;

  int elecstyle = NONE;
  if (collide) elecstyle = collide->elecstyle;
  if (elecstyle == DISCRETE) {
    Species species = particle->species[isp];
    if (species.elecdat == NULL) return 0.0;

    electronic_distribution_func(isp, temp_elec);

    for (int i = 1; i < species.elecdat->nelecstate; ++i)
      cumulative_probabilities[i] += cumulative_probabilities[i-1];

    double ran = erandom->uniform();
    ielec = 0;
    while (ran > cumulative_probabilities[ielec])
      ++ielec;
  }
  return ielec;
}

/* ---------------------------------------------------------------------- */

double* Particle::electronic_distribution_func(int isp, double temp_elec) {
  double* distribution = cumulative_probabilities;
  int elecstyle = NONE;
  if (collide) elecstyle = collide->elecstyle;
  if (elecstyle == DISCRETE) {
    Particle::Species species = particle->species[isp];
    double partition_function = 0.0;

    for (int i = 0; i < species.elecdat->nelecstate; ++i) {

      // Calculate boltzmann fractions

      distribution[i] = species.elecdat->states[i].degen*exp(-species.elecdat->states[i].temp/temp_elec);

      // Calculate partition function

      partition_function += distribution[i];
    }

    for (int i = 0; i < species.elecdat->nelecstate; ++i)
      distribution[i] /= partition_function;
  }

  return distribution;
}

/* ----------------------------------------------------------------------
   read list of species defined in species file
   store info in filespecies and nfile
   only invoked by proc 0
------------------------------------------------------------------------- */

void Particle::read_species_file()
{
  // read file line by line
  // skip blank lines or comment lines starting with '#'
  // all other lines must have NWORDS

  int NWORDS = 10;
  char **words = new char*[NWORDS];
  char line[MAXLINE],copy[MAXLINE];

  while (fgets(line,MAXLINE,fp)) {
    int pre = strspn(line," \t\n\r");
    if (pre == strlen(line) || line[pre] == '#') continue;

    strcpy(copy,line);
    int nwords = wordcount(copy,NULL);
    if (nwords != NWORDS)
      error->one(FLERR,"Incorrect line format in species file");

    if (nfile == maxfile) {
      maxfile += DELTASPECIES;
      filespecies = (Species *)
        memory->srealloc(filespecies,maxfile*sizeof(Species),
                         "particle:filespecies");
      memset(&filespecies[nfile],0,(maxfile-nfile)*sizeof(Species));
    }

    nwords = wordcount(line,words);
    Species *fsp = &filespecies[nfile];

    if (strlen(words[0]) + 1 > 16)
      error->one(FLERR,"Invalid species ID in species file");
    strcpy(fsp->id,words[0]);

    fsp->molwt = atof(words[1]);
    fsp->mass = atof(words[2]);
    fsp->rotdof = atoi(words[3]);
    fsp->rotrel = atof(words[4]);
    fsp->vibdof = atoi(words[5]);
    fsp->vibrel[0] = atof(words[6]);
    fsp->vibtemp[0] = atof(words[7]);
    fsp->specwt = atof(words[8]);
    fsp->charge = atof(words[9]);

    if (fsp->rotdof > 0 || fsp->vibdof > 0) fsp->internaldof = 1;
    else fsp->internaldof = 0;

    // error checks

    if (fsp->rotdof != 0 && fsp->rotdof != 2 && fsp->rotdof != 3)
      error->all(FLERR,"Invalid rotational DOF in species file");

    if (fsp->vibdof < 0 || fsp->vibdof > 2*MAXVIBMODE || fsp->vibdof % 2)
      error->all(FLERR,"Invalid vibrational DOF in species file");

    // initialize additional rotation/vibration fields
    // may be overwritten by rotfile or vibfile

    fsp->nrottemp = 0;
    fsp->nvibmode = fsp->vibdof / 2;

    fsp->rottemp[0] = fsp->rottemp[1] = fsp->rottemp[2] = 0.0;

    fsp->vibdegen[0] = 0;
    for (int m = 1; m < MAXVIBMODE; m++) {
      fsp->vibtemp[m] = 0.0;
      fsp->vibrel[m] = 0.0;
      fsp->vibdegen[m] = 0;
    }

    fsp->vibdiscrete_read = 0;

    nfile++;
  }

  delete [] words;

  fclose(fp);
}

/* ----------------------------------------------------------------------
   read list of extra rotation info in rotation file
   store info in filerot and nfile
   only invoked by proc 0
------------------------------------------------------------------------- */

void Particle::read_rotation_file()
{
  // read file line by line
  // skip blank lines or comment lines starting with '#'
  // all other lines can have up to NWORDS

  int NWORDS = 5;
  char **words = new char*[NWORDS];
  char line[MAXLINE],copy[MAXLINE];

  while (fgets(line,MAXLINE,fp)) {
    int pre = strspn(line," \t\n\r");
    if (pre == strlen(line) || line[pre] == '#') continue;

    strcpy(copy,line);
    int nwords = wordcount(copy,NULL);
    if (nwords > NWORDS)
      error->one(FLERR,"Incorrect line format in rotation file");

    if (nfile == maxfile) {
      maxfile += DELTASPECIES;
      filerot = (RotFile *)
        memory->srealloc(filerot,maxfile*sizeof(RotFile),
                         "particle:filerot");
      memset(&filerot[nfile],0,(maxfile-nfile)*sizeof(RotFile));
    }

    nwords = wordcount(line,words);
    RotFile *rsp = &filerot[nfile];

    if (strlen(words[0]) + 1 > 16)
      error->one(FLERR,"Invalid species ID in rotation file");
    strcpy(rsp->id,words[0]);

    rsp->ntemp = atoi(words[1]);
    if (rsp->ntemp != 1 && rsp->ntemp != 3)
      error->one(FLERR,"Invalid N count in rotation file");
    if (nwords != 2 + rsp->ntemp)
      error->one(FLERR,"Incorrect line format in rotation file");

    int j = 2;
    for (int i = 0; i < rsp->ntemp; i++)
      rsp->rottemp[i] = atof(words[j++]);

    nfile++;
  }

  delete [] words;

  fclose(fp);
}

/* ----------------------------------------------------------------------
   read list of extra rotation info in vibration file
   store info in filevib and nfile
   only invoked by proc 0
------------------------------------------------------------------------- */

void Particle::read_vibration_file()
{
  // read file line by line
  // skip blank lines or comment lines starting with '#'
  // all other lines can have up to NWORDS

  int NWORDS = 2 + 3*MAXVIBMODE;
  char **words = new char*[NWORDS];
  char line[MAXLINE],copy[MAXLINE];

  while (fgets(line,MAXLINE,fp)) {
    int pre = strspn(line," \t\n\r");
    if (pre == strlen(line) || line[pre] == '#') continue;

    strcpy(copy,line);
    int nwords = wordcount(copy,NULL);
    if (nwords > NWORDS)
      error->one(FLERR,"Incorrect line format in vibration file");

    if (nfile == maxfile) {
      maxfile += DELTASPECIES;
      filevib = (VibFile *)
        memory->srealloc(filevib,maxfile*sizeof(VibFile),
                         "particle:filevib");
      memset(&filevib[nfile],0,(maxfile-nfile)*sizeof(VibFile));
    }

    nwords = wordcount(line,words);
    VibFile *vsp = &filevib[nfile];

    if (strlen(words[0]) + 1 > 16)
      error->one(FLERR,"Invalid species ID in vibration file");
    strcpy(vsp->id,words[0]);

    vsp->nmode = atoi(words[1]);
    if (vsp->nmode < 2 || vsp->nmode > MAXVIBMODE)
      error->one(FLERR,"Invalid N count in vibration file");
    if (nwords != 2 + 3*vsp->nmode)
      error->one(FLERR,"Incorrect line format in vibration file");

    int j = 2;
    for (int i = 0; i < vsp->nmode; i++) {
      vsp->vibtemp[i] = atof(words[j++]);
      vsp->vibrel[i] = atof(words[j++]);
      vsp->vibdegen[i] = atoi(words[j++]);
    }
    nfile++;
  }

  delete [] words;

  fclose(fp);
}

/* ----------------------------------------------------------------------
   read list of extra electronic info in electronic file
   store info in fileelec and nfile
   only invoked by proc 0
------------------------------------------------------------------------- */

void Particle::read_electronic_file()
{
  // read file line by line
  // skip blank lines or comment lines starting with '#'

  char **words = new char*[128];
  char line[MAXLINE],copy[MAXLINE];

  while (fgets(line,MAXLINE,fp)) {
    int pre = strspn(line," \t\n\r");
    if (pre == strlen(line) || line[pre] == '#') continue;

    strcpy(copy,line);
    int nwords = wordcount(copy,NULL);
    if (nwords < 3)
      error->one(FLERR,"Incorrect line format in electronic file");

    if (nfile == maxfile) {
      maxfile += DELTASPECIES;
      fileelec = (ElecFile *)
      memory->srealloc(fileelec,maxfile*sizeof(ElecFile),
                       "particle:fileelec");
      memset(&fileelec[nfile],0,(maxfile-nfile)*sizeof(ElecFile));
    }

    nwords = wordcount(line,words);
    ElecFile *vsp = &fileelec[nfile];

    if (strlen(words[0]) + 1 > 16)
      error->one(FLERR,"Invalid species ID in electronic file");
    strcpy(vsp->id,words[0]);

    int isp = particle->find_species(words[0]);
    if (isp < 0) continue;

    int test_nmode = atoi(words[1]);
    if (test_nmode > 0) {
      vsp->nmode = test_nmode;

      // Line defining electronic states for one species
      if (vsp->nmode < 2)
        error->one(FLERR,"Invalid N count in electronic file");
      if (nwords != 2 + 4*vsp->nmode)
        error->one(FLERR,"Incorrect line format in electronic file");

      vsp->default_rel = (double*) memory->smalloc(vsp->nmode*sizeof(double), "vsp:default_rel");
      vsp->elecrel = (double**) memory->smalloc(nspecies*sizeof(double*), "vsp:elecrel");
      for (int i = 0; i < nspecies; ++i) {
        vsp->elecrel[i] = (double*) memory->smalloc(vsp->nmode*sizeof(double), "vsp:elecrel[]");
      }
      vsp->enforce_spin_conservation = (bool*) memory->smalloc(nspecies*sizeof(bool), "vsp:enforce_spin_conservation");
      vsp->electemp = (double*) memory->smalloc(vsp->nmode*sizeof(double), "vsp:electemp");
      vsp->elecdegen = (int*) memory->smalloc(vsp->nmode*sizeof(int), "vsp:elecdegen");
      vsp->elecspin = (int*) memory->smalloc(vsp->nmode*sizeof(int), "vsp:elecspin");

      for (int j = 0; j < nspecies; ++j) {
        for (int i = 0; i < vsp->nmode; ++i) {
          vsp->elecrel[j][i] = -1.0;
        }
        vsp->enforce_spin_conservation[j] = true;
      }

      int j = 2;
      for (int i = 0; i < vsp->nmode; ++i) {
        vsp->electemp[i] = atof(words[j++]);
        vsp->default_rel[i] = atof(words[j++]);
        vsp->elecdegen[i] = atoi(words[j++]);
        vsp->elecspin[i] = atoi(words[j++]);
      }
      nfile++;
    } else {
      // Cross-species line defining species-specific relaxation collision numbers
      for (int i = 0; i < nfile; ++i) {
        if (strcmp(words[0],fileelec[i].id) == 0) {
          vsp = &fileelec[i];
          break;
        }
      }
      if (nwords != 3 + vsp->nmode)
        error->one(FLERR,"Incorrect line format in electronic file");

      int jsp = particle->find_species(words[1]);
      if (jsp < 0) continue;

      if (strcmp(words[2],"T") == 0) {
        vsp->enforce_spin_conservation[jsp] = true;
      } else if (strcmp(words[2],"F") == 0) {
        vsp->enforce_spin_conservation[jsp] = false;
      } else {
        error->one(FLERR,"Must specify 'T' or 'F' for spin conservation flag");
      }
      int j = 3;
      for (int i = 0; i < vsp->nmode; ++i) {
        vsp->elecrel[jsp][i] = atof(words[j++]);
      }
    }
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
  char *word = strtok(line," \t\n");

  while (word) {
    if (words) words[nwords] = word;
    nwords++;
    word = strtok(NULL," \t\n");
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
  int tmp;

  if (me == 0) tmp = fread(&nspecies,sizeof(int),1,fp);
  MPI_Bcast(&nspecies,1,MPI_INT,0,world);

  if (nspecies > maxspecies) {
    while (nspecies > maxspecies) maxspecies += DELTASPECIES;
    grow_species();
  }

  if (me == 0) tmp = fread(species,sizeof(Species),nspecies,fp);
  MPI_Bcast(species,nspecies*sizeof(Species),MPI_CHAR,0,world);

  maxvibmode = 0;
  for (int isp = 0; isp < nspecies; isp++)
    maxvibmode = MAX(maxvibmode,species[isp].nvibmode);
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
  int tmp;

  // must first clear existing default mixtures

  for (int i = 0; i < nmixture; i++) delete mixture[i];
  nmixture = 0;

  // now process restart file data

  if (me == 0) tmp = fread(&nmixture,sizeof(int),1,fp);
  MPI_Bcast(&nmixture,1,MPI_INT,0,world);

  if (nmixture > maxmixture) {
    while (nmixture > maxmixture) maxmixture += DELTAMIXTURE;
    mixture = (Mixture **)
      memory->srealloc(mixture,maxmixture*sizeof(Mixture *),"particle:mixture");
  }

  int n;
  char *id;

  for (int i = 0; i < nmixture; i++) {
    if (me == 0) tmp = fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    id = new char[n];
    if (me == 0) tmp = fread(id,sizeof(char),n,fp);
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
   return size of particle restart info for this proc
------------------------------------------------------------------------- */

bigint Particle::size_restart_big()
{
  bigint n = sizeof(int);
  n = BIROUNDUP(n);
  n += nlocal * sizeof(OnePartRestart);
  n += nlocal * sizeof_custom();
  n = BIROUNDUP(n);
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
   pack my particle info into buf
   use multiple passes to reduce memory use
   use OnePartRestart data struct for permanent info and to encode cell ID
   include per-particle custom attributes if defined
   NOTE: does not ROUNDUP(ptr) at the end, this is done by caller
------------------------------------------------------------------------- */

void Particle::pack_restart(char *buf, int step, int pass)
{
  Grid::ChildCell *cells = grid->cells;
  OnePart *p;
  OnePartRestart *pr;
  int nbytes_custom = sizeof_custom();

  char *ptr = buf;
  if (pass == 0) {
    int *ibuf = (int *) ptr;
    *ibuf = nlocal;
    ptr += sizeof(int);
    ptr = ROUNDUP(ptr);
  }

  int start = step*pass;
  int end = start+step;
  end = MIN(nlocal,end);
  for (int i = start; i < end; i++) {
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
    memory->smalloc(nlocal_restart*nbytes,"particle:particle_restart");

  memcpy(particle_restart,ptr,nlocal_restart*nbytes);
  ptr += nlocal_restart * sizeof(OnePartRestart);
  ptr = ROUNDUP(ptr);

  return ptr - buf;
}

/* ----------------------------------------------------------------------
   unpack particle info into restart storage
   use multiple passes to reduce memory use
   allocate data structure here, will be deallocated by ReadRestart
   include per-particle custom attributes if defined
   NOTE: does not ROUNDUP(ptr) at the end, this is done by caller
------------------------------------------------------------------------- */

void Particle::unpack_restart(char *buf, int &nlocal_restart, int step, int pass)
{
  int nbytes_particle = sizeof(OnePartRestart);
  int nbytes_custom = sizeof_custom();
  int nbytes = nbytes_particle + nbytes_custom;

  char *ptr = buf;
  if (pass == 0) {
    int *ibuf = (int *) buf;
    nlocal_restart = *ibuf;
    ptr += sizeof(int);
    ptr = ROUNDUP(ptr);
  }
  int start = step*pass;
  int end = start+step;
  end = MIN(nlocal_restart,end);
  step = end - start;

  particle_restart = (char *)
    memory->smalloc(step*nbytes,"particle:particle_restart");

  memcpy(particle_restart,ptr,step*nbytes);

  this->nlocal_restart = step;
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

/* ---------------------------------------------------------------------- */

double Particle::elec_energy(int isp, double temp_elec) {
  Particle::Species species = particle->species[isp];

  double* state_probabilities = particle->electronic_distribution_func(isp, temp_elec);

  double total_energy = 0.0;
  for (int i = 0; i < species.elecdat->nelecstate; ++i)
    total_energy += state_probabilities[i]*species.elecdat->states[i].temp*update->boltz;

  return total_energy;
}

/* ---------------------------------------------------------------------- */

double Particle::bisectTelec(int isp, double eelec, int count)
{
  double first_elec_eng, t_elec, degen0, degen1, numer, denom;
  double boltz = update->boltz;

  // We calculate a first guess at the temp assuming
  //   all the electronic energy is stored in the first
  //   excited state

  double target_energy_per_part = eelec/count;
  first_elec_eng = species[isp].elecdat->states[1].temp*boltz;
  degen0 = species[isp].elecdat->states[0].degen;
  degen1 = species[isp].elecdat->states[1].degen;
  t_elec = first_elec_eng / (boltz*(
      - log( target_energy_per_part*degen0 /
          (first_elec_eng*degen1)
         ))
    );
  // If the electronic excitation is very high, the above calculation will
  //   give negative numbers, including -inf. Negative temperatures are physically
  //   meaningful if over 50% of particles are in excited states. However, in the
  //   case of a truly broken value (-inf) we initialize the guess at something
  //   more sane

  if (isinf(t_elec)) t_elec = -species[isp].elecdat->states[1].temp;

  // Bisection method to find T accurate to 1%

  // Find initial bounds based on our first guess

  double T_low = 0.9*t_elec;
  while (elec_energy(isp, T_low) > target_energy_per_part) {
    T_low /= 2.0;
  }

  double T_high = t_elec*1.1;
  while (elec_energy(isp, T_high) < target_energy_per_part) {
    T_high *= 2.0;
  }

  // Bisect

  double T_mid = t_elec;
  double e_mid = elec_energy(isp, T_mid);
  while ((T_high - T_low) > 0.01) {
    if (e_mid > target_energy_per_part) {
      T_high = T_mid;
    } else {
      T_low = T_mid;
    }
    T_mid = (T_high - T_low)/2.0 + T_low;
    e_mid = elec_energy(isp, T_mid);
  }

  return T_mid;
}

/* ---------------------------------------------------------------------- */

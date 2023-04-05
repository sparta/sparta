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
#include "string.h"
#include "compute_count.h"
#include "update.h"
#include "particle.h"
#include "mixture.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{SPECIES,MIXTURE};

#define DELTAVALUES 4

/* ---------------------------------------------------------------------- */

ComputeCount::ComputeCount(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal compute count command");

  nvalues = maxvalues =  0;
  spmix = index = indexgroup = mixgroups = NULL;
  allocate(narg-2);

  for (int iarg = 2; iarg < narg; iarg++) {
    if (strchr(arg[iarg],'/')) {
      char *ptr = strchr(arg[iarg],'/');
      *ptr = '\0';
      int imix = particle->find_mixture(arg[iarg]);
      if (imix < 0)
        error->all(FLERR,"Unknown species/mixture in compute count command");
      int igroup = particle->mixture[imix]->find_group(ptr+1);
      if (imix < 0)
        error->all(FLERR,"Unknown mixture group in compute count command");
      *ptr = '/';
      if (nvalues == maxvalues) allocate(maxvalues+DELTAVALUES);
      spmix[nvalues] = MIXTURE;
      index[nvalues] = imix;
      indexgroup[nvalues] = igroup;
      mixgroups[nvalues] = particle->mixture[imix]->ngroup;
      nvalues++;

    } else {
      int isp = particle->find_species(arg[iarg]);
      int imix = particle->find_mixture(arg[iarg]);
      if (isp < 0 && imix < 0)
        error->all(FLERR,"Unknown species/mixture in compute count command");
      if (isp >= 0 && imix >= 0)
        error->all(FLERR,"Ambiguous species/mixture in compute count command");
      if (isp >= 0) {
        if (nvalues == maxvalues) allocate(maxvalues+DELTAVALUES);
        spmix[nvalues] = SPECIES;
        index[nvalues] = isp;
        nvalues++;
      } else {
        int ngroup = particle->mixture[imix]->ngroup;
        if (nvalues+ngroup > maxvalues) allocate(nvalues+ngroup);
        for (int i = 0; i < ngroup; i++) {
          spmix[nvalues] = MIXTURE;
          index[nvalues] = imix;
          indexgroup[nvalues] = i;
          mixgroups[nvalues] = ngroup;
          nvalues++;
        }
      }
    }
  }

  // initialization

  maxspecies = 0;
  count = NULL;

  lasttally = -1;
  vector = NULL;
  onevec = sumvec = NULL;

  if (nvalues == 1) scalar_flag = 1;
  else {
    vector_flag = 1;
    size_vector = nvalues;
    vector = new double[nvalues];
    onevec = new bigint[nvalues];
    sumvec = new bigint[nvalues];
  }
}

/* ---------------------------------------------------------------------- */

ComputeCount::~ComputeCount()
{
  if (copymode) return;

  memory->destroy(spmix);
  memory->destroy(index);
  memory->destroy(indexgroup);
  memory->destroy(mixgroups);

  delete [] count;
  delete [] vector;
  delete [] onevec;
  delete [] sumvec;
}


/* ---------------------------------------------------------------------- */

void ComputeCount::init()
{
  // insure count vector is long enough for species count

  if (maxspecies < particle->nspecies) {
    maxspecies = particle->nspecies;
    delete [] count;
    count = new int[maxspecies];
  }

  // check if the group count in any accessed mixtures has changed

  int warn = 0;
  int err = 0;
  for (int i = 0; i < nvalues; i++) {
    if (spmix[i] == SPECIES) continue;
    int imix = index[i];
    int igroup = indexgroup[i];
    int ngroup = mixgroups[i];
    if (ngroup != particle->mixture[imix]->ngroup) warn = 1;
    if (igroup >= particle->mixture[imix]->ngroup) err = 1;
  }

  if (err)
    error->all(FLERR,
               "Group in mixture used by compute count no longer valid");
  if (warn && comm->me == 0)
    error->warning(FLERR,
                   "Group count in mixture used by compute count has changed");
}

/* ---------------------------------------------------------------------- */

double ComputeCount::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  per_species_tally();

  bigint one;
  if (spmix[0] == SPECIES) one = count[index[0]];
  else {
    int nspecies = particle->mixture[index[0]]->groupsize[indexgroup[0]];
    one = 0;
    for (int m = 0; m < nspecies; m++) {
      int isp = particle->mixture[index[0]]->groupspecies[indexgroup[0]][m];
      one += count[isp];
    }
  }

  bigint sum;
  MPI_Allreduce(&one,&sum,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  scalar = sum;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeCount::compute_vector()
{
  int i,m;

  invoked_scalar = update->ntimestep;

  per_species_tally();

  for (i = 0; i < nvalues; i++) {
    onevec[i] = 0;
    if (spmix[i] == SPECIES) onevec[i] = count[index[i]];
    else {
      int nspecies = particle->mixture[index[i]]->groupsize[indexgroup[i]];
      for (m = 0; m < nspecies; m++) {
        int isp = particle->mixture[index[i]]->groupspecies[indexgroup[i]][m];
        onevec[i] += count[isp];
      }
    }
  }

  MPI_Allreduce(onevec,sumvec,nvalues,MPI_SPARTA_BIGINT,MPI_SUM,world);
  for (i = 0; i < nvalues; i++) vector[i] = sumvec[i];
}

/* ---------------------------------------------------------------------- */

void ComputeCount::per_species_tally()
{
  if (lasttally == update->ntimestep) return;
  lasttally = update->ntimestep;

  int nspecies = particle->nspecies;
  for (int i = 0; i < nspecies; i++) count[i] = 0;

  Particle::OnePart *particles = particle->particles;
  int nlocal = particle->nlocal;

  for (int i = 0; i < nlocal; i++)
    count[particles[i].ispecies]++;
}

/* ---------------------------------------------------------------------- */

void ComputeCount::allocate(int n)
{
  maxvalues = n;
  memory->grow(spmix,maxvalues,"count:spmix");
  memory->grow(index,maxvalues,"count:index");
  memory->grow(indexgroup,maxvalues,"count:indexgroup");
  memory->grow(mixgroups,maxvalues,"count:mixgroups");
}

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
#include "collide.h"
#include "particle.h"
#include "mixture.h"
#include "update.h"
#include "grid.h"
#include "comm.h"
#include "random_mars.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

#define DELTAGRID 1000            // must be bigger than split cells per cell

/* ---------------------------------------------------------------------- */

Collide::Collide(SPARTA *sparta, int narg, char **arg) : Pointers(sparta)
{
  int n = strlen(arg[0]) + 1;
  style = new char[n];
  strcpy(style,arg[0]);

  n = strlen(arg[1]) + 1;
  mixID = new char[n];
  strcpy(mixID,arg[1]);

  random = new RanPark(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,comm->me,100);

  ngroups = 0;

  npmax = 0;
  plist = NULL;

  nglocal = nglocalmax = 0;

  ngroup = NULL;
  maxgroup = NULL;
  glist = NULL;
  gpair = NULL;

  vre_first = 1;
  vre_start = 1;
  vre_every = 0;
  remainflag = 1;
  vremax = NULL;
  vremax_initial = NULL;
  remain = NULL;

  // initialize counters in case stats outputs them

  ncollide_one = nattempt_one = nreact_one = 0;
  ncollide_running = nattempt_running = nreact_running = 0;
}

/* ---------------------------------------------------------------------- */

Collide::~Collide()
{
  delete [] style;
  delete [] mixID;
  delete random;

  if (ngroups == 1) memory->destroy(plist);
  if (ngroups > 1) {
    delete [] ngroup;
    delete [] maxgroup;
    for (int i = 0; i < ngroups; i++) memory->destroy(glist[i]);
    delete [] glist;
    memory->destroy(gpair);
  }

  memory->destroy(vremax);
  memory->destroy(vremax_initial);
  memory->destroy(remain);
}

/* ---------------------------------------------------------------------- */

void Collide::init()
{
  // require mixture to contain all species

  int imix = particle->find_mixture(mixID);
  if (imix < 0) error->all(FLERR,"Collision mixture does not exist");
  mixture = particle->mixture[imix];

  if (mixture->nspecies != particle->nspecies)
    error->all(FLERR,"Collision mixture does not contain all species");

  // reallocate one-cell data structs for one or many groups

  int oldgroups = ngroups;
  ngroups = mixture->ngroup;

  if (ngroups != oldgroups) {
    if (oldgroups == 1) {
      memory->destroy(plist);
      plist = NULL;
    }
    if (oldgroups > 1) {
      delete [] ngroup;
      delete [] maxgroup;
      for (int i = 0; i < oldgroups; i++) memory->destroy(glist[i]);
      delete [] glist;
      memory->destroy(gpair);
      ngroup = NULL;
      maxgroup = NULL;
      glist = NULL;
      gpair = NULL;
    }

    if (ngroups == 1) {
      npmax = DELTAPART;
      memory->create(plist,npmax,"collide:plist");
    }
    if (ngroups > 1) {
      ngroup = new int[ngroups];
      maxgroup = new int[ngroups];
      glist = new int*[ngroups];
      for (int i = 0; i < ngroups; i++) {
        maxgroup[i] = DELTAPART;
        memory->create(glist[i],DELTAPART,"collide:glist");
      }
      memory->create(gpair,ngroups*ngroups,3,"collide:gpair");
    }
  }

  // allocate vremax,remain if group count changed
  // will always be allocated on first run since oldgroups = 0
  // set vremax_intitial via values calculated by collide style

  if (ngroups != oldgroups) {
    memory->destroy(vremax);
    memory->destroy(vremax_initial);
    memory->destroy(remain);
    nglocal = grid->nlocal;
    nglocalmax = nglocal;
    memory->create(vremax,nglocalmax,ngroups,ngroups,"collide:vremax");
    memory->create(vremax_initial,ngroups,ngroups,"collide:vremax_initial");
    if (remainflag)
      memory->create(remain,nglocalmax,ngroups,ngroups,"collide:remain");

    for (int igroup = 0; igroup < ngroups; igroup++)
      for (int jgroup = 0; jgroup < ngroups; jgroup++)
        vremax_initial[igroup][jgroup] = vremax_init(igroup,jgroup);
  }

  // vre_next = next timestep to zero vremax & remain, based on vre_every

  if (vre_every) vre_next = (update->ntimestep/vre_every)*vre_every + vre_every;
  else vre_next = update->laststep + 1;

  // if requested reset vremax & remain
  // must be after per-species vremax_initial is setup

  if (vre_first || vre_start) {
    reset_vremax();
    vre_first = 0;
  }

  // initialize running stats before each run

  ncollide_running = nattempt_running = nreact_running = 0;
}

/* ---------------------------------------------------------------------- */

void Collide::modify_params(int narg, char **arg)
{
  if (narg == 0) error->all(FLERR,"Illegal collide_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"vremax") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal collide_modify command");
      vre_every = atoi(arg[iarg+1]);
      if (vre_every < 0) error->all(FLERR,"Illegal collide_modify command");
      if (strcmp(arg[iarg+2],"yes") == 0) vre_start = 1;
      else if (strcmp(arg[iarg+2],"no") == 0) vre_start = 0;
      else error->all(FLERR,"Illegal collide_modify command");
      iarg += 3;
    } else if (strcmp(arg[iarg],"remain") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal collide_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) remainflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) remainflag = 0;
      else error->all(FLERR,"Illegal collide_modify command");
      iarg += 2;
    } else error->all(FLERR,"Illegal collide_modify command");
  }
}

/* ----------------------------------------------------------------------
   reset vremax to initial species-based values
   reset remain to 0.0
------------------------------------------------------------------------- */

void Collide::reset_vremax()
{
  int ngroups = mixture->ngroup;

  for (int icell = 0; icell < nglocal; icell++)
    for (int igroup = 0; igroup < ngroups; igroup++)
      for (int jgroup = 0; jgroup < ngroups; jgroup++) {
        vremax[icell][igroup][jgroup] = vremax_initial[igroup][jgroup];
        if (remainflag) remain[icell][igroup][jgroup] = 0.0;
      }
}

/* ----------------------------------------------------------------------
  NTC algorithm
------------------------------------------------------------------------- */

void Collide::collisions()
{
  // if requested, reset vrwmax & remain

  if (update->ntimestep == vre_next) {
    reset_vremax();
    vre_next += vre_every;
  }

  // counters

  ncollide_one = nattempt_one = nreact_one = 0;

  // perform collisions
  // one variant is optimized for a single group

  if (ngroups == 1) collisions_one();
  else collisions_group();

  // accumulate running totals

  nattempt_running += nattempt_one;
  ncollide_running += ncollide_one;
  nreact_running += nreact_one;
}

/* ----------------------------------------------------------------------
  NTC algorithm on a single group
------------------------------------------------------------------------- */

void Collide::collisions_one()
{
  int i,j,k,n,ip,jp,np;
  int nattempt;
  double attempt,volume;
  Particle::OnePart *ipart,*jpart,*kpart;

  // loop over cells I own

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  Particle::OnePart *particles = particle->particles;
  int *next = particle->next;

  for (int icell = 0; icell < nglocal; icell++) {
    np = cinfo[icell].count;
    if (np <= 1) continue;
    ip = cinfo[icell].first;
    volume = cinfo[icell].volume;

    // setup particle list for this cell

    if (np > npmax) {
      npmax = np + DELTAPART;
      memory->destroy(plist);
      memory->create(plist,npmax,"collide:plist");
    }

    n = 0;
    while (ip >= 0) {
      plist[n++] = ip;
      ip = next[ip];
    }

    // attempt = exact collision attempt count for a pair of groups
    // nattempt = rounded attempt with RN

    attempt = attempt_collision(icell,np,volume);
    nattempt = static_cast<int> (attempt);

    if (!nattempt) continue;
    nattempt_one += nattempt;

    // perform collisions
    // select random particles, cannot be same
    // test if collision actually occurs
    // if chemistry occurs, move output I,J,K particles to new group lists
    // if chemistry occurs, exit attempt loop if group count goes to 0

    for (k = 0; k < nattempt; k++) {
      i = np * random->uniform();
      j = np * random->uniform();
      while (i == j) j = np * random->uniform();

      ipart = &particles[plist[i]];
      jpart = &particles[plist[j]];

      if (!test_collision(icell,0,0,ipart,jpart)) continue;
      setup_collision(ipart,jpart);
      kpart = perform_collision(ipart,jpart);
      ncollide_one++;

      // jpart destroyed, delete from plist
      
      if (!jpart) {
        plist[j] = plist[np-1];
        np--;
        if (np <= 1) break;
      }
      
      // if kpart created, add to plist
      // kpart was just added to particle list, so index = nlocal-1
      
      if (kpart) {
        if (np == npmax) {
          npmax = np + DELTAPART;
          memory->grow(plist,npmax,"collide:plist");
        }
        plist[np++] = particle->nlocal-1;
      }
    }
  }
}

/* ----------------------------------------------------------------------
  NTC algorithm on pairs of groups
  pre-compute # of attempts per group pair
------------------------------------------------------------------------- */

void Collide::collisions_group()
{
  int i,j,k,ip,jp,np,isp,ipair,igroup,jgroup,newgroup;
  int nattempt;
  int *ni,*nj,*ilist,*jlist;
  double attempt,volume;
  Particle::OnePart *ipart,*jpart,*kpart;

  // counters

  ncollide_one = nattempt_one = nreact_one = 0;

  // loop over cells I own

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  Particle::OnePart *particles = particle->particles;
  int *next = particle->next;

  int *species2group = mixture->species2group;

  for (int icell = 0; icell < nglocal; icell++) {
    np = cinfo[icell].count;
    if (np <= 1) continue;
    ip = cinfo[icell].first;
    volume = cinfo[icell].volume;

    // setup per-group particle lists for this cell

    for (i = 0; i < ngroups; i++) ngroup[i] = 0;

    while (ip >= 0) {
      isp = particles[ip].ispecies;
      igroup = species2group[isp];
      if (ngroup[igroup] == maxgroup[igroup]) {
	maxgroup[igroup] += DELTAPART;
	memory->grow(glist[igroup],maxgroup[igroup],"collide:grouplist");
      }
      glist[igroup][ngroup[igroup]++] = ip;
      ip = next[ip];
    }

    // attempt = exact collision attempt count for a pair of groups
    // nattempt = rounded attempt with RN
    // add pairing to gpair

    npair = 0;
    for (igroup = 0; igroup < ngroups; igroup++)
      for (jgroup = 0; jgroup < ngroups; jgroup++) {
	attempt = attempt_collision(icell,igroup,jgroup,volume);
	nattempt = static_cast<int> (attempt);

	if (nattempt) {
	  gpair[npair][0] = igroup;
	  gpair[npair][1] = jgroup;
	  gpair[npair][2] = nattempt;
	  nattempt_one += nattempt;
	  npair++;
	}
      }

    // perform collisions for each pair of groups in gpair list
    // select random particle in each group
    // if igroup = jgroup, cannot be same particle
    // test if collision actually occurs
    // if chemistry occurs, move output I,J,K particles to new group lists
    // if chemistry occurs, exit attempt loop if group count goes to 0
    // NOTE: need to reset vremax ?
    // NOTE: OK to use pre-computed nattempt when Ngroup may have changed?

    for (ipair = 0; ipair < npair; ipair++) {
      igroup = gpair[ipair][0];
      jgroup = gpair[ipair][1];
      nattempt = gpair[ipair][2];

      ni = &ngroup[igroup];
      nj = &ngroup[jgroup];
      ilist = glist[igroup];
      jlist = glist[jgroup];

      if (*ni == 0 || *nj == 0) continue;
      if (igroup == jgroup && *ni == 1) continue;

      for (k = 0; k < nattempt; k++) {
	i = *ni * random->uniform();
	j = *nj * random->uniform();
	if (igroup == jgroup)
	  while (i == j) j = *nj * random->uniform();

	ipart = &particles[ilist[i]];
	jpart = &particles[jlist[j]];

	if (!test_collision(icell,igroup,jgroup,ipart,jpart)) continue;
	setup_collision(ipart,jpart);
	kpart = perform_collision(ipart,jpart);
	ncollide_one++;

	// ipart may now be in different group

	newgroup = species2group[ipart->ispecies];
	if (newgroup != igroup) {
	  addgroup(newgroup,ilist[i]);
	  ilist[i] = ilist[*ni-1];
	  (*ni)--;
	  if (*ni <= 1) {
	    if (*ni == 0) break;
	    if (igroup == jgroup) break;
	  }
	}

	// jpart may be in different group or 

	if (jpart) {
	  newgroup = species2group[jpart->ispecies];
	  if (newgroup != jgroup) {
	    addgroup(newgroup,jlist[j]);
	    jlist[j] = jlist[*nj-1];
	    (*nj)--;
	    if (*nj <= 1) {
	      if (*nj == 0) break;
	      if (igroup == jgroup) break;
	    }
	  }
	} else {
	  jlist[j] = jlist[*nj-1];
	  (*nj)--;
	  if (*nj <= 1) {
	    if (*nj == 0) break;
	    if (igroup == jgroup) break;
	  }
	}

        // if kpart created, add to group list
	// kpart was just added to particle list, so index = nlocal-1

	if (kpart) {
	  newgroup = species2group[kpart->ispecies];
	  addgroup(newgroup,particle->nlocal-1);
	}
      }
    }
  }

  // accumulate running totals

  nattempt_running += nattempt_one;
  ncollide_running += ncollide_one;
  nreact_running += nreact_one;
}

/* ----------------------------------------------------------------------
   pack icell values for per-cell arrays into buf
   if icell is a split cell, also pack all sub cell values 
   return byte count of amount packed
   if memflag, only return count, do not fill buf
   NOTE: why packing/unpacking parent cell if a split cell?
------------------------------------------------------------------------- */

int Collide::pack_grid_one(int icell, char *buf, int memflag)
{
  if (!vremax) return 0;
  int nbytes = ngroups*ngroups*sizeof(double);

  Grid::ChildCell *cells = grid->cells;

  int n;
  if (remainflag) {
    if (memflag) {
      memcpy(buf,&vremax[icell][0][0],nbytes);
      memcpy(&buf[nbytes],&remain[icell][0][0],nbytes);
    }
    n = 2*nbytes;
  } else {
    if (memflag) memcpy(buf,&vremax[icell][0][0],nbytes);
    n = nbytes;
  }

  if (cells[icell].nsplit > 1) {
    int isplit = cells[icell].isplit;
    int nsplit = cells[icell].nsplit;
    for (int i = 0; i < nsplit; i++) {
      int m = grid->sinfo[isplit].csubs[i];
      if (remainflag) {
        if (memflag) {
          memcpy(&buf[n],&vremax[m][0][0],nbytes);
          n += nbytes;
          memcpy(&buf[n],&remain[m][0][0],nbytes);
          n += nbytes;
        } else n += 2*nbytes;
      } else {
        if (memflag) memcpy(&buf[n],&vremax[m][0][0],nbytes);
        n += nbytes;
      }
    }
  }
  
  return n;
}

/* ----------------------------------------------------------------------
   unpack icell values for per-cell arrays from buf
   if icell is a split cell, also unpack all sub cell values 
   return byte count of amount unpacked
------------------------------------------------------------------------- */

int Collide::unpack_grid_one(int icell, char *buf)
{
  if (!vremax) return 0;
  int nbytes = ngroups*ngroups*sizeof(double);

  Grid::ChildCell *cells = grid->cells;
  Grid::SplitInfo *sinfo = grid->sinfo;

  grow_percell(1);
  memcpy(&vremax[icell][0][0],buf,nbytes);
  int n = nbytes;
  if (remainflag) {
    memcpy(&remain[icell][0][0],&buf[n],nbytes);
    n += nbytes;
  }
  nglocal++;

  if (cells[icell].nsplit > 1) {
    int isplit = cells[icell].isplit;
    int nsplit = cells[icell].nsplit;
    grow_percell(nsplit);
    for (int i = 0; i < nsplit; i++) {
      int m = sinfo[isplit].csubs[i];
      memcpy(&vremax[m][0][0],&buf[n],nbytes);
      n += nbytes;
      if (remainflag) {
        memcpy(&remain[m][0][0],&buf[n],nbytes);
        n += nbytes;
      }
    }
    nglocal += nsplit;
  }
  
  return n;
}

/* ----------------------------------------------------------------------
   compress per-cell arrays due to cells migrating to new procs
   criteria for keeping/discarding a cell is same as in Grid::compress()
   this keeps final ordering of per-cell arrays consistent with Grid class
------------------------------------------------------------------------- */

void Collide::compress_grid()
{
  if (!vremax) return;
  int nbytes = ngroups*ngroups*sizeof(double);

  int me = comm->me;
  Grid::ChildCell *cells = grid->cells;

  // keep an unsplit or split cell if staying on this proc
  // keep a sub cell if its split cell is staying on this proc

  int ncurrent = nglocal;
  nglocal = 0;
  for (int icell = 0; icell < ncurrent; icell++) {
    if (cells[icell].nsplit >= 1) {
      if (cells[icell].proc != me) continue;
    } else {
      int isplit = cells[icell].isplit;
      if (cells[grid->sinfo[isplit].icell].proc != me) continue;
    }

    if (nglocal != icell) { 
      memcpy(&vremax[nglocal][0][0],&vremax[icell][0][0],nbytes);
      if (remainflag) 
        memcpy(&remain[nglocal][0][0],&remain[icell][0][0],nbytes);
    }
    nglocal++;
  }
}

/* ----------------------------------------------------------------------
   insure per-cell arrays are allocated long enough for N new cells
------------------------------------------------------------------------- */

void Collide::grow_percell(int n)
{
  if (nglocal+n < nglocalmax) return;
  nglocalmax += DELTAGRID;
  memory->grow(vremax,nglocalmax,ngroups,ngroups,"collide:vremax");
  if (remainflag) 
    memory->grow(remain,nglocalmax,ngroups,ngroups,"collide:remain");
}

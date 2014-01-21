/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2012) Sandia Corporation.  Under the terms of Contract
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

  ngroup = NULL;
  maxgroup = NULL;
  glist = NULL;
  gpair = NULL;

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

  // initialize running stats before each run

  ncollide_running = nattempt_running = nreact_running = 0;
}

/* ---------------------------------------------------------------------- */

void Collide::modify_params(int narg, char **arg)
{
}

/* ----------------------------------------------------------------------
  NTC algorithm
------------------------------------------------------------------------- */

void Collide::collisions()
{
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
  int nglocal = grid->nlocal;

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
  int nglocal = grid->nlocal;

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

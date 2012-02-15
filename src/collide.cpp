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

using namespace DSMC_NS;

#define DELTAPART 128
#define SEED 12345

/* ---------------------------------------------------------------------- */

Collide::Collide(DSMC *dsmc, int narg, char **arg) : Pointers(dsmc)
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
  ngroup = NULL;
  maxgroup = NULL;
  glist = NULL;
  gpair = NULL;
}

/* ---------------------------------------------------------------------- */

Collide::~Collide()
{
  delete [] style;
  delete [] mixID;
  delete random;

  delete [] ngroup;
  delete [] maxgroup;
  for (int i = 0; i < ngroups; i++) memory->destroy(glist[i]);
  delete [] glist;
  memory->destroy(gpair);
}

/* ---------------------------------------------------------------------- */

void Collide::init()
{
  // require mixture to contain all species

  int imix = particle->find_mixture(mixID);
  if (imix < 0) error->all(FLERR,"Collision mixture does not exist");
  mixture = particle->mixture[imix];

  if (mixture->allspecies == 0)
    error->all(FLERR,"Collision mixture does not exist");

  // reallocate one-cell data structs that depend on # of groups

  int oldgroups = ngroups;
  ngroups = mixture->ngroups;

  if (ngroups != oldgroups) {
    delete [] ngroup;
    delete [] maxgroup;
    for (int i = 0; i < oldgroups; i++) memory->destroy(glist[i]);
    delete [] glist;
    memory->destroy(gpair);

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

/* ----------------------------------------------------------------------
  NTC algorithm on pairs of groups
  pre-compute # of attempts per group pair
------------------------------------------------------------------------- */

void Collide::collisions()
{
  int i,j,k,ip,jp,np,icell,itype,igroup,jgroup,nattempt;
  double attempt,volume;
  Particle::OnePart *ipart,*jpart,*kpart;

  Grid::OneCell *cells = grid->cells;
  int *mycells = grid->mycells;
  int nglocal = grid->nlocal;

  Particle::OnePart *particles = particle->particles;
  int *cellcount = particle->cellcount;
  int *first = particle->first;
  int *next = particle->next;

  int *species2group = mixture->species2group;

  // loop over cells I own

  ncoll_attempt = 0;
  ncoll = 0;

  for (int m = 0; m < nglocal; m++) {
    icell = mycells[m];
    volume = cells[icell].volume;
    np = cellcount[m];
    ip = first[m];

    // setup per-group particle lists for this cell

    for (i = 0; i < ngroups; i++) ngroup[i] = 0;

    while (ip >= 0) {
      itype = particles[ip].ispecies;
      igroup = species2group[itype];
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
	if (attempt-nattempt > random->uniform()) nattempt++;

	if (nattempt) {
	  gpair[npair][0] = igroup;
	  gpair[npair][1] = jgroup;
	  gpair[npair][2] = nattempt;
	  ncoll_attempt += nattempt;
	  npair++;
	}
      }

    // perform collisions for each group pairing in gpair list
    // select random particle in each group
    // if igroup = jgroup, cannot be same particle
    // test if collision actually occurs
    // NOTE: will need to reset vremax
    // NOTE: will need to account for chemistry & update group lists

    for (i = 0; i < npair; i++) {
      igroup = gpair[i][0];
      jgroup = gpair[i][1];
      nattempt = gpair[i][2];

      for (k = 0; k < nattempt; k++) {
	i = ngroup[igroup] * random->uniform();
	j = ngroup[jgroup] * random->uniform();
	if (igroup == jgroup)
	  while (i == j) j = ngroup[jgroup] * random->uniform();

	ipart = &particles[glist[igroup][i]];
	jpart = &particles[glist[jgroup][j]];

	if (!test_collision(icell,igroup,jgroup,ipart,jpart)) continue;
	setup_collision(ipart,jpart);
	kpart = perform_collision(ipart,jpart);
	ncoll++;
      }
    }
  }
}

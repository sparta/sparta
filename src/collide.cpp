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
#include "grid.h"
#include "comm.h"
#include "random_park.h"
#include "memory.h"

using namespace DSMC_NS;

#define DELTAPART 128
#define SEED 12345

/* ---------------------------------------------------------------------- */

Collide::Collide(DSMC *dsmc, int narg, char **arg) : Pointers(dsmc)
{
  int n = strlen(arg[0]) + 1;
  style = new char[n];
  strcpy(style,arg[0]);

  nspecies = 0;
  nsp = NULL;
  maxsp = NULL;
  splist = NULL;

  sscoll = NULL;

  // NOTE: need to handle RNG and seeds better, like kMC
  random = new RanPark(dsmc,SEED+comm->me);
}

/* ---------------------------------------------------------------------- */

Collide::~Collide()
{
  delete [] style;
  delete random;

  delete [] nsp;
  delete [] maxsp;
  for (int i = 0; i < nspecies; i++) memory->destroy(splist[i]);
  delete [] splist;
  memory->destroy(sscoll);
}

/* ---------------------------------------------------------------------- */

void Collide::init()
{
  oldspecies = nspecies;
  nspecies = particle->nspecies;

  // reallocate one-cell data structures that depend on # of species
  // NOTE: this assumes nglocal is static

  if (nspecies != oldspecies) {
    delete [] nsp;
    delete [] maxsp;
    for (int i = 0; i < oldspecies; i++) memory->destroy(splist[i]);
    delete [] splist;
    memory->destroy(sscoll);

    nsp = new int[nspecies];
    maxsp = new int[nspecies];
    splist = new int*[nspecies];
    for (int i = 0; i < nspecies; i++) {
      maxsp[i] = DELTAPART;
      memory->create(splist[i],DELTAPART,"collide:splist");
    }
    memory->create(sscoll,nspecies*nspecies,3,"collide:sscoll");
  }

  ncollattempt = 0;
  ncollision = 0;
}

/* ---------------------------------------------------------------------- */

void Collide::collisions()
{
  int i,j,k,icell,itype,isp,jsp,ip,jp,np;
  int nattempt;
  double volume;
  Particle::OnePart *ipart,*jpart,*kpart;

  Grid::OneCell *cells = grid->cells;
  int *mycells = grid->mycells;
  Particle::OnePart *particles = particle->particles;
  int *cellcount = particle->cellcount;
  int *first = particle->first;
  int *next = particle->next;
  int nglocal = grid->nlocal;
  int nspecies = particle->nspecies;




  return;



  // loop over cells I own
  // NOTE: this is the NTC algorithm, could do simpler TC instead

  for (int m = 0; m < nglocal; m++) {

    icell = mycells[m];
    volume = cells[icell].volume;
    np = cellcount[m];;
    ip = first[m];

    // setup per-species particle lists for this cell

    for (i = 0; i < nspecies; i++) nsp[i] = 0;

    while (ip >= 0) {
      itype = particles[ip].ispecies;
      if (nsp[itype] == maxsp[itype]) {
	maxsp[itype] += DELTAPART;
	memory->grow(splist[itype],maxsp[itype],"collide:splist");
      }
      splist[itype][nsp[itype]++] = ip;
      ip = next[ip];
    }

    // compute collision attempt count for each species pair

    nsspair = 0;
    for (isp = 0; isp < nspecies; isp++)
      for (jsp = 0; jsp < nspecies; jsp++) {
	// NOTE: comment out for now to exercise list buidling
	//nattempt = attempt_collision(icell,isp,jsp,volume);






	// species is off by 1 to Ntypes
	//if (nsp[0]




	nattempt = 10;
	if (nattempt) {
	  sscoll[nsspair][0] = isp;
	  sscoll[nsspair][1] = jsp;
	  sscoll[nsspair][2] = nattempt;
	  nsspair++;
	}
      }

    // perform collisions for each species pair in sscoll list

    for (i = 0; i < nsspair; i++) {
      isp = sscoll[i][0];
      jsp = sscoll[i][1];
      nattempt = sscoll[i][2];
      ncollattempt += nattempt;

      for (k = 0; k < nattempt; k++) {

	// select random particle of each species
	// if isp = jsp, cannot be same particle

	i = nsp[isp] * random->uniform();
	j = nsp[jsp] * random->uniform();
	if (isp == jsp)
	  while (i == j) j = nsp[jsp] * random->uniform();

	ipart = &particles[splist[isp][i]];
	jpart = &particles[splist[jsp][j]];

	// test if collision actually occurs
	// if so, perform a collision
	// kpart = NULL means no new particle
        // NOTE: worry about created kpart when add chemistry

	if (!test_collision(icell,isp,jsp,ipart,jpart)) continue;
	setup_collision(ipart,jpart);
	kpart = perform_collision(ipart,jpart);
	ncollision++;
      }
    }
  }
}

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
#include "react.h"
#include "modify.h"
#include "fix.h"
#include "fix_ambipolar.h"
#include "random_mars.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{NONE,DISCRETE,SMOOTH};       // several files
enum{PKEEP,PINSERT,PDONE,PDISCARD,PENTRY,PEXIT,PSURF};   // several files

#define DELTAGRID 1000            // must be bigger than split cells per cell
#define DELTADELETE 1024
#define DELTAELECTRON 128

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

  maxdelete = 0;
  dellist = NULL;

  vre_first = 1;
  vre_start = 1;
  vre_every = 0;
  remainflag = 1;
  vremax = NULL;
  vremax_initial = NULL;
  remain = NULL;
  rotstyle = SMOOTH;
  vibstyle = NONE;

  ambiflag = 0;
  maxelectron = 0;
  elist = NULL;

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

  memory->destroy(dellist);
  memory->sfree(elist);
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

  // find ambipolar fix
  // set ambipolar vector/array indices
  // if reactions defined, check that they are valid ambipolar reactions

  if (ambiflag) {
    index_ionambi = particle->find_custom((char *) "ionambi");
    index_velambi = particle->find_custom((char *) "velambi");
    if (index_ionambi < 0 || index_velambi < 0) 
      error->all(FLERR,"Collision ambipolar without fix ambipolar");
    if (react) react->ambi_check();

    int ifix;
    for (ifix = 0; ifix < modify->nfix; ifix++)
      if (strcmp(modify->fix[ifix]->style,"ambipolar") == 0) break;
    FixAmbipolar *afix = (FixAmbipolar *) modify->fix[ifix];
    ambispecies = afix->especies;
    ions = afix->ions;
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
    } else if (strcmp(arg[iarg],"rotate") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal collide_modify command");
      if (strcmp(arg[iarg+1],"no") == 0) rotstyle = NONE;
      else if (strcmp(arg[iarg+1],"yes") == 0) rotstyle = SMOOTH;
      else error->all(FLERR,"Illegal collide_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"vibrate") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal collide_modify command");
      if (strcmp(arg[iarg+1],"no") == 0) vibstyle = NONE;
      else if (strcmp(arg[iarg+1],"discrete") == 0) vibstyle = DISCRETE;
      else if (strcmp(arg[iarg+1],"smooth") == 0) vibstyle = SMOOTH;
      else error->all(FLERR,"Illegal collide_modify command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"ambipolar") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal collide_modify command");
      if (strcmp(arg[iarg+1],"no") == 0) ambiflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) ambiflag = 1;
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
  ndelete = 0;

  // perform collisions without or with ambipolar approximation
  // one variant is optimized for a single group

  if (!ambiflag) {
    if (ngroups == 1) collisions_one();
    else collisions_group();
  } else {
    if (ngroups == 1) collisions_one_ambipolar();
    else collisions_group_ambipolar();
  }

  // remove any particles deleted in chemistry reactions
  // if particles deleted/created by chemistry, particles are no longer sorted
  // NOTE: not sure if this is the case, need to check
  //       important if grid adapt will happen at end of timsestepx

  if (ndelete) particle->compress_reactions(ndelete,dellist);
  if (react) particle->sorted = 0;

  // accumulate running totals

  nattempt_running += nattempt_one;
  ncollide_running += ncollide_one;
  nreact_running += nreact_one;
}

/* ----------------------------------------------------------------------
   NTC algorithm for a single group
------------------------------------------------------------------------- */

void Collide::collisions_one()
{
  int i,j,k,n,ip,np;
  int nattempt,reactflag;
  double attempt,volume;
  Particle::OnePart *ipart,*jpart,*kpart;

  // loop over cells I own

  Grid::ChildInfo *cinfo = grid->cinfo;

  Particle::OnePart *particles = particle->particles;
  int *next = particle->next;

  for (int icell = 0; icell < nglocal; icell++) {
    np = cinfo[icell].count;
    if (np <= 1) continue;
    ip = cinfo[icell].first;
    volume = cinfo[icell].volume / cinfo[icell].weight;

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
    // select random pair of particles, cannot be same
    // test if collision actually occurs

    for (k = 0; k < nattempt; k++) {
      i = np * random->uniform();
      j = np * random->uniform();
      while (i == j) j = np * random->uniform();

      ipart = &particles[plist[i]];
      jpart = &particles[plist[j]];

      // test if collision actually occurs, then perform it
      // ijspecies = species before collision chemistry
      // continue to next collision if no reaction

      if (!test_collision(icell,0,0,ipart,jpart)) continue;
      setup_collision(ipart,jpart);
      reactflag = perform_collision(ipart,jpart,kpart);
      ncollide_one++;
      if (reactflag) nreact_one++;
      else continue;

      // jpart destroyed, delete from plist
      // also add particle to deletion list
      // exit attempt loop if only single particle left

      if (!jpart) {
        if (ndelete == maxdelete) {
          maxdelete += DELTADELETE;
          memory->grow(dellist,maxdelete,"collide:dellist");
        }
        dellist[ndelete++] = plist[j];
        plist[j] = plist[np-1];
        np--;
        if (np < 2) break;
      }
      
      // if kpart created, add to plist
      // kpart was just added to particle list, so index = nlocal-1
      // particle data structs may have been realloced by kpart
      
      if (kpart) {
        if (np == npmax) {
          npmax = np + DELTAPART;
          memory->grow(plist,npmax,"collide:plist");
        }
        plist[np++] = particle->nlocal-1;
        particles = particle->particles;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   NTC algorithm for multiple groups, loop over pairs of groups
   pre-compute # of attempts per group pair
------------------------------------------------------------------------- */

void Collide::collisions_group()
{
  int i,j,k,ip,np,isp,ipair,igroup,jgroup,newgroup;
  int nattempt,reactflag;
  int *ni,*nj,*ilist,*jlist;
  double attempt,volume;
  Particle::OnePart *ipart,*jpart,*kpart;

  // loop over cells I own

  Grid::ChildInfo *cinfo = grid->cinfo;

  Particle::OnePart *particles = particle->particles;
  int *next = particle->next;
  int *species2group = mixture->species2group;

  for (int icell = 0; icell < nglocal; icell++) {
    np = cinfo[icell].count;
    if (np <= 1) continue;
    ip = cinfo[icell].first;
    volume = cinfo[icell].volume / cinfo[icell].weight;

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
    // double loop over N^2 / 2 pairs of groups
    // nattempt = rounded attempt with RN
    // add pair of groups to gpair

    npair = 0;
    for (igroup = 0; igroup < ngroups; igroup++)
      for (jgroup = igroup; jgroup < ngroups; jgroup++) {
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

        // test if collision actually occurs, then perform it
        // ijspecies = species before collision chemistry
        // continue to next collision if no reaction

	if (!test_collision(icell,igroup,jgroup,ipart,jpart)) continue;
	setup_collision(ipart,jpart);
	reactflag = perform_collision(ipart,jpart,kpart);
	ncollide_one++;
        if (reactflag) nreact_one++;
        else continue;

	// ipart may now be in different group
        // reset jlist after addgroup() b/c may have realloced if igroup=jgroup

	newgroup = species2group[ipart->ispecies];
	if (newgroup != igroup) {
	  addgroup(newgroup,ilist[i]);
          jlist = glist[jgroup];
	  ilist[i] = ilist[*ni-1];
	  (*ni)--;
          // this line needed if jgroup=igroup and just moved jlist[j]
          if (jlist == ilist && j == *ni) j = i;
	  if (*ni <= 1) {
	    if (*ni == 0) break;
	    if (igroup == jgroup) break;
	  }
	}

	// jpart may now be in different group or destroyed
        // reset ilist after addgroup() b/c may have realloced if igroup=jgroup

	if (jpart) {
	  newgroup = species2group[jpart->ispecies];
	  if (newgroup != jgroup) {
	    addgroup(newgroup,jlist[j]);
            ilist = glist[igroup];
	    jlist[j] = jlist[*nj-1];
	    (*nj)--;
	    if (*nj <= 1) {
	      if (*nj == 0) break;
	      if (igroup == jgroup) break;
	    }
	  }
	} else {
          if (ndelete == maxdelete) {
            maxdelete += DELTADELETE;
            memory->grow(dellist,maxdelete,"collide:dellist");
          }
          dellist[ndelete++] = jlist[j];
	  jlist[j] = jlist[*nj-1];
	  (*nj)--;
	  if (*nj <= 1) {
	    if (*nj == 0) break;
	    if (igroup == jgroup) break;
	  }
	}

        // if kpart created, add to group list
	// kpart was just added to particle list, so index = nlocal-1
        // reset ilist,jlist after addgroup() b/c may have been realloced
        // particles data struct may also have been realloced

	if (kpart) {
	  newgroup = species2group[kpart->ispecies];
	  addgroup(newgroup,particle->nlocal-1);
          ilist = glist[igroup];
          jlist = glist[jgroup];
          particles = particle->particles;
	}
      }
    }
  }
}

/* ----------------------------------------------------------------------
   NTC algorithm for a single group with ambipolar approximation
------------------------------------------------------------------------- */

void Collide::collisions_one_ambipolar()
{
  int i,j,k,n,ip,np,ispecies,jspecies,tmp;
  int nattempt,reactflag;
  double attempt,volume;
  Particle::OnePart *ipart,*jpart,*kpart,*p,*ep;

  // ambipolar vectors

  int *ionambi = particle->eivec[particle->ewhich[index_ionambi]];
  double **velambi = particle->edarray[particle->ewhich[index_velambi]];

  // loop over cells I own

  Grid::ChildInfo *cinfo = grid->cinfo;

  Particle::OnePart *particles = particle->particles;
  int *next = particle->next;
  int nbytes = sizeof(Particle::OnePart);

  for (int icell = 0; icell < nglocal; icell++) {
    np = cinfo[icell].count;
    if (np <= 1) continue;
    ip = cinfo[icell].first;
    volume = cinfo[icell].volume / cinfo[icell].weight;

    // DEBUG test that there are no electrons
    // can remove at some point

    /*
    while (ip >= 0) {
      if (particles[ip].ispecies == ambispecies)
        error->one(FLERR,"Pre-collision particle is ambipolar electron");
      ip = next[ip];
    }
    ip = cinfo[icell].first;
    */

    // setup particle list for this cell
    // allow for up to Np extra electrons

    if (2*np > npmax) {
      npmax = 2*np + DELTAPART;
      memory->destroy(plist);
      memory->create(plist,npmax,"collide:plist");
    }

    n = 0;
    while (ip >= 0) {
      plist[n++] = ip;
      ip = next[ip];
    }

    // grow electron array as needed

    if (np >= maxelectron) {
      while (maxelectron < np) maxelectron += DELTAELECTRON;
      memory->sfree(elist);
      elist = (Particle::OnePart *) 
        memory->smalloc(maxelectron*nbytes,"collide:elist");
    }
    
    // create electrons for ambipolar ions, then increment np
    // create them in separate array since will never become real particles
    // plist indexes them with negative indices (-1 to -Nelectron)
    // ion->flag stores same negative index of its matching electron

    nelectron = 0;
    for (i = 0; i < np; i++) {
      if (ionambi[plist[i]]) {
        p = &particles[plist[i]];
        ep = &elist[nelectron];
        memcpy(ep,p,nbytes);
        memcpy(ep->v,velambi[plist[i]],3*sizeof(double));
        ep->ispecies = ambispecies;
        nelectron++;
        plist[n++] = -nelectron;
        p->flag = -nelectron;
      }
    }

    np += nelectron;

    // attempt = exact collision attempt count for a pair of groups
    // nattempt = rounded attempt with RN
    // if no attempts, reset particle flags to PKEEP before continuing

    attempt = attempt_collision(icell,np,volume);
    nattempt = static_cast<int> (attempt);

    if (!nattempt) {
      np -= nelectron;
      for (i = 0; i < np; i++)
        if (ionambi[plist[i]]) particles[plist[i]].flag = PKEEP;
      continue;
    }
    nattempt_one += nattempt;

    // perform collisions
    // select random pair of particles, cannot be same
    // test if collision actually occurs
    // if chemistry occurs, exit attempt loop if group count goes to 0

    for (k = 0; k < nattempt; k++) {
      i = np * random->uniform();
      j = np * random->uniform();
      while (i == j) j = np * random->uniform();

      // plist index >= 0 for particles array
      // plist index < 0 for electron array

      if (plist[i] >= 0) ipart = &particles[plist[i]];
      else ipart = &elist[-plist[i]-1];
      if (plist[j] >= 0) jpart = &particles[plist[j]];
      else jpart = &elist[-plist[j]-1];

      // check for e/e pair
      // count as collision, but do not perform it

      if (ipart->ispecies == ambispecies && jpart->ispecies == ambispecies) {
        ncollide_one++;
        continue;
      }

      // if particle I is electron
      // swap with J, since electron must be 2nd in any ambipolar reaction
      // just need to swap i/j, ipart/jpart
      // don't have to worry if an ambipolar ion is I or J

      if (ipart->ispecies == ambispecies) {
        tmp = i;
        i = j;
        j = tmp;
        p = ipart;
        ipart = jpart;
        jpart = p;
      }

      // test if collision actually occurs, then perform it
      // ijspecies = species before collision chemistry
      // continue to next collision if no reaction

      if (!test_collision(icell,0,0,ipart,jpart)) continue;
      ispecies = ipart->ispecies;
      jspecies = jpart->ispecies;
      setup_collision(ipart,jpart);
      reactflag = perform_collision(ipart,jpart,kpart);
      ncollide_one++;
      if (reactflag) nreact_one++;
      else continue;


      // reset ambipolar ions and ion/electron pairings due to reaction
      // must do now before group reset below can break out of loop
      // first reset ionambi if added kpart since ambi_reset uses it

      if (kpart) ionambi = particle->eivec[particle->ewhich[index_ionambi]];
      ambi_reset(plist[i],plist[j],ispecies,jspecies,
                 ipart,jpart,kpart,ionambi);

      // jpart destroyed, delete from plist
      // also add particle to deletion list
      // exit attempt loop if only single particle left

      if (!jpart) {
        if (ndelete == maxdelete) {
          maxdelete += DELTADELETE;
          memory->grow(dellist,maxdelete,"collide:dellist");
        }
        dellist[ndelete++] = plist[j];
        plist[j] = plist[np-1];
        np--;
        if (np < 2) break;
      }
      
      // if kpart created, add to plist
      // kpart was just added to particle list, so index = nlocal-1
      // particles and custom data structs may have been realloced by kpart
      
      if (kpart) {
        if (np == npmax) {
          npmax = np + DELTAPART;
          memory->grow(plist,npmax,"collide:plist");
        }
        plist[np++] = particle->nlocal-1;
        particles = particle->particles;
        ionambi = particle->eivec[particle->ewhich[index_ionambi]];
        velambi = particle->edarray[particle->ewhich[index_velambi]];
      }
    }

    // DEBUG test that ions and electrons are still matched one-to-one
    // flag each electron with -1 as find ions
    // can remove at some point

    /*
    int nelec = 0;
    int nion = 0;

    for (n = 0; n < np; n++) {
      i = plist[n];
      if (i < 0 || particles[i].ispecies == ambispecies) {
        nelec++;
        continue;
      }
      p = &particles[i];
      if (ionambi[i]) {
        nion++;
        if (p->flag >= 0) ep = &particles[p->flag];
        else ep = &elist[-(p->flag)-1];
        if (ep->ispecies != ambispecies)
          error->one(FLERR,"Ambipolar ion is not coupled to electron");
        if (ep->flag < 0) 
          error->one(FLERR,"Ambipolar electron is coupled to multiple ions");
        ep->flag = -1;
      }
    }

    if (nion != nelec)
      error->one(FLERR,"Ambipolar ion and electron counts do not match");
    */

    // done with collisions/chemistry for one grid cell
    // recombine ambipolar ions with their matching electrons
    //   by copying electron velocity into velambi
    // reset all flags to PKEEP
    // add any newly created electrons (not in elist) to delete list

    for (n = 0; n < np; n++) {
      i = plist[n];
      if (i < 0) continue;
      p = &particles[i];
      p->flag = PKEEP;
      if (ionambi[i]) {
        if (p->flag >= 0) ep = &particles[p->flag];
        else ep = &elist[-(p->flag)-1];
        memcpy(velambi[i],ep->v,3*sizeof(double));
      } else if (p->ispecies == ambispecies) {
        if (ndelete == maxdelete) {
          maxdelete += DELTADELETE;
          memory->grow(dellist,maxdelete,"collide:dellist");
        }
        dellist[ndelete++] = i;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   NTC algorithm for multiple groups with ambipolar approximation
   loop over pairs of groups, pre-compute # of attempts per group pair
------------------------------------------------------------------------- */

void Collide::collisions_group_ambipolar()
{
  int i,j,k,n,ip,np,isp,ipair,igroup,jgroup,newgroup,ispecies,jspecies,tmp;
  int nattempt,reactflag;
  int *ni,*nj,*ilist,*jlist,*tmpvec;
  double attempt,volume;
  Particle::OnePart *ipart,*jpart,*kpart,*p,*ep;

  // ambipolar vectors

  int *ionambi = particle->eivec[particle->ewhich[index_ionambi]];
  double **velambi = particle->edarray[particle->ewhich[index_velambi]];

  // loop over cells I own

  Grid::ChildInfo *cinfo = grid->cinfo;

  Particle::OnePart *particles = particle->particles;
  int *next = particle->next;
  int nbytes = sizeof(Particle::OnePart);
  int *species2group = mixture->species2group;

  int egroup = species2group[ambispecies];

  for (int icell = 0; icell < nglocal; icell++) {
    np = cinfo[icell].count;
    if (np <= 1) continue;
    ip = cinfo[icell].first;
    volume = cinfo[icell].volume / cinfo[icell].weight;

    // DEBUG test that there are no electrons
    // can remove at some point

    /*
    while (ip >= 0) {
      if (particles[ip].ispecies == ambispecies)
        error->one(FLERR,"Pre-collision particle is ambipolar electron");
      ip = next[ip];
    }
    ip = cinfo[icell].first;
    */

    // grow electron array as needed

    if (np >= maxelectron) {
      while (maxelectron < np) maxelectron += DELTAELECTRON;
      memory->sfree(elist);
      elist = (Particle::OnePart *) 
        memory->smalloc(maxelectron*nbytes,"collide:elist");
    }

    // setup per-group particle lists for this cell
    // create electrons for ambipolar ions, and put in egroup
    // create them in separate array since will never become real particles
    // group lists index them with negative indices (-1 to -Nelectron)
    // ion->flag stores same negative index of its matching electron

    for (i = 0; i < ngroups; i++) ngroup[i] = 0;
    nelectron = 0;

    while (ip >= 0) {
      isp = particles[ip].ispecies;
      igroup = species2group[isp];
      if (ngroup[igroup] == maxgroup[igroup]) {
	maxgroup[igroup] += DELTAPART;
	memory->grow(glist[igroup],maxgroup[igroup],"collide:grouplist");
      }
      glist[igroup][ngroup[igroup]++] = ip;

      if (ionambi[ip]) {
        p = &particles[ip];
        ep = &elist[nelectron];
        memcpy(ep,p,nbytes);
        memcpy(ep->v,velambi[ip],3*sizeof(double));
        ep->ispecies = ambispecies;
        nelectron++;
        if (ngroup[egroup] == maxgroup[egroup]) {
          maxgroup[egroup] += DELTAPART;
          memory->grow(glist[egroup],maxgroup[egroup],"collide:grouplist");
        }
        glist[egroup][ngroup[egroup]++] = -nelectron;
        p->flag = -nelectron;
      }

      ip = next[ip];
    }

    // attempt = exact collision attempt count for a pair of groups
    // double loop over N^2 / 2 pairs of groups
    // nattempt = rounded attempt with RN
    // add pair of groups to gpair

    npair = 0;
    for (igroup = 0; igroup < ngroups; igroup++)
      for (jgroup = igroup; jgroup < ngroups; jgroup++) {
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

        // ilist/jlist indices >= 0 for particles array
        // ilist/jlist indices < 0 for electron array

        if (ilist[i] >= 0) ipart = &particles[ilist[i]];
        else ipart = &elist[-ilist[i]-1];
        if (jlist[j] >= 0) jpart = &particles[jlist[j]];
        else jpart = &elist[-jlist[j]-1];

        // check for e/e pair
        // count as collision, but do not perform it
        
        if (ipart->ispecies == ambispecies && jpart->ispecies == ambispecies) {
          ncollide_one++;
          continue;
        }

        // if particle I is electron
        // swap with J, since electron must be 2nd in any ambipolar reaction
        // need to swap i/j, igroup/jgroup, ni/nj, ilist/jlist, ipart/jpart
        // don't have to worry if an ambipolar ion is I or J

        if (ipart->ispecies == ambispecies) {
          tmp = i;
          i = j;
          j = tmp;
          tmp = igroup;
          igroup = jgroup;
          jgroup = tmp;
          tmpvec = ni;
          ni = nj;
          nj = tmpvec;
          tmpvec = ilist;
          ilist = jlist;
          jlist = tmpvec;
          p = ipart;
          ipart = jpart;
          jpart = p;
        }

        // test if collision actually occurs, then perform it
        // ijspecies = species before collision chemistry
        // continue to next collision if no reaction

	if (!test_collision(icell,igroup,jgroup,ipart,jpart)) continue;
        ispecies = ipart->ispecies;
        jspecies = jpart->ispecies;
	setup_collision(ipart,jpart);
	reactflag = perform_collision(ipart,jpart,kpart);
	ncollide_one++;
        if (reactflag) nreact_one++;
        else continue;

        // reset ambipolar ions and ion/electron pairings due to reaction
        // must do now before group reset below can break out of loop
        // first reset ionambi if added kpart since ambi_reset uses it

        if (kpart) ionambi = particle->eivec[particle->ewhich[index_ionambi]];
        ambi_reset(ilist[i],jlist[j],ispecies,jspecies,
                   ipart,jpart,kpart,ionambi);

	// ipart may now be in different group
        // reset jlist after addgroup() b/c may have realloced if igroup=jgroup

	newgroup = species2group[ipart->ispecies];
	if (newgroup != igroup) {
	  addgroup(newgroup,ilist[i]);
          jlist = glist[jgroup];
	  ilist[i] = ilist[*ni-1];
	  (*ni)--;
          // this line needed if jgroup=igroup and just moved jlist[j]
          if (jlist == ilist && j == *ni) j = i;
	  if (*ni <= 1) {
	    if (*ni == 0) break;
	    if (igroup == jgroup) break;
	  }
	}

	// jpart may now be in different group or destroyed
        // reset ilist after addgroup() b/c may have realloced if igroup=jgroup

	if (jpart) {
	  newgroup = species2group[jpart->ispecies];
	  if (newgroup != jgroup) {
	    addgroup(newgroup,jlist[j]);
            ilist = glist[igroup];
	    jlist[j] = jlist[*nj-1];
	    (*nj)--;
            if (*nj <= 1) {
	      if (*nj == 0) break;
	      if (igroup == jgroup) break;
	    }
	  }

	} else {
          if (ndelete == maxdelete) {
            maxdelete += DELTADELETE;
            memory->grow(dellist,maxdelete,"collide:dellist");
          }
          dellist[ndelete++] = jlist[j];
	  jlist[j] = jlist[*nj-1];
	  (*nj)--;
	  if (*nj <= 1) {
	    if (*nj == 0) break;
	    if (igroup == jgroup) break;
	  }
	}

        // if kpart created, add to group list
	// kpart was just added to particle list, so index = nlocal-1
        // reset ilist,jlist after addgroup() b/c may have been realloced
        // particle and custom data structs may have been realloced by kpart

	if (kpart) {
	  newgroup = species2group[kpart->ispecies];
	  addgroup(newgroup,particle->nlocal-1);
          ilist = glist[igroup];
          jlist = glist[jgroup];
          particles = particle->particles;
          ionambi = particle->eivec[particle->ewhich[index_ionambi]];
          velambi = particle->edarray[particle->ewhich[index_velambi]];
	}
      }
    }

    // DEBUG test that ions and electrons are still matched one-to-one
    // flag each electron with -1 as find ions
    // can remove at some point

    /*
    int nelec = 0;
    int nion = 0;

    for (k = 0; k < ngroups; k++) {
      n = ngroup[k];
      for (j = 0; j < n; j++) {
        i = glist[k][j];
        if (i < 0 || particles[i].ispecies == ambispecies) {
          nelec++;
          continue;
        }
        p = &particles[i];
        if (ionambi[i]) {
          nion++;
          if (p->flag >= 0) ep = &particles[p->flag];
          else ep = &elist[-(p->flag)-1];
          if (ep->ispecies != ambispecies)
            error->one(FLERR,"Ambipolar ion is not coupled to electron");
          if (ep->flag < 0) 
            error->one(FLERR,"Ambipolar electron is coupled to multiple ions");
          ep->flag = -1;
        }
      }
    }

    if (nion != nelec)
      error->one(FLERR,"Ambipolar ion and electron counts do not match");
    */

    // done with collisions/chemistry for one grid cell
    // recombine ambipolar ions with their matching electrons
    //   by copying electron velocity into velambi
    // reset all flags to PKEEP
    // add any newly created electrons (not in elist) to delete list

    for (k = 0; k < ngroups; k++) {
      n = ngroup[k];
      for (j = 0; j < n; j++) {
        i = glist[k][j];
        if (i < 0) continue;
        p = &particles[i];
        p->flag = PKEEP;
        if (ionambi[i]) {
          if (p->flag >= 0) ep = &particles[p->flag];
          else ep = &elist[-(p->flag)-1];
          memcpy(velambi[i],ep->v,3*sizeof(double));
        } else if (p->ispecies == ambispecies) {
          if (ndelete == maxdelete) {
            maxdelete += DELTADELETE;
            memory->grow(dellist,maxdelete,"collide:dellist");
          }
          dellist[ndelete++] = i;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   reset ionambi and ion/electron coupling if ambipolar reaction occurred
   i/j = indices of I,J reactants
   isp/jsp = pre-reaction species of I,J
     both will not be electrons, if one is electron it will be jsp
   reactants i,j and isp/jsp will always be in order listed below
   products ip,jp,kp will always be in order listed below
   logic must be valid for all ambipolar AND non-ambipolar reactions
   check for 3 versions of 2 -> 3: dissociation or ionization
     all have J product = electron
     D: AB + e -> A + e + B
        if I reactant = neutral and K product not electron:
        set K product = neutral
     D: AB+ + e -> A+ + e + B
        if I reactant = ion:
        set K product = neutral
     I: A + e -> A+ + e + e
        if I reactant = neutral and K product = electron:
        set I product = ion, couple I ion to K electron
     all other 2 -> 3 cases, set K product = neutral
   check for 3 versions of 2 -> 2: ionization or exchange
     I: A + B -> AB+ + e
        if J product = electron:
        set I product to ion, couple I ion to J electron
     E: AB+ + C -> A + BC+
        if I reactant = ion:
        set I/J products to neutral/ion, couple J ion to I's original electron
     E: C + AB+ -> A + BC+
        if J reactant = ion:
        nothing to change for products
     all other 2 -> 2 cases, no changes
   check for one version of 2 -> 1: recombination
     R: A+ + e -> A
        if ej = elec, set I product to not ion
     all other 2 -> 1 cases, no changes
   WARNING:
     do not index by I,J if could be e, since may be negative I,J index
     do not access ionambi if could be e, since e may be in elist
------------------------------------------------------------------------- */

void Collide::ambi_reset(int i, int j, int isp, int jsp, 
                         Particle::OnePart *ip, Particle::OnePart *jp, 
                         Particle::OnePart *kp, int *ionambi)
{
  int e = ambispecies;

  // 2 reactants become 3 products
  // in all ambi reactions, J reactant is electron

  if (kp) {
    int k = particle->nlocal-1;
    ionambi[k] = 0;
    if (jsp != e) return;

    if (ionambi[i]) {                // nothing more to change
    } else if (kp->ispecies == e) {
      ionambi[i] = 1;                // 1st reactant is now ion
      ip->flag = k;                  // couple I ion to K electron
    } else {                         // nothing more to change
    }

  // 2 reactants become 2 products
  // ambi reaction if J product is electron or either reactant is ion

  } else if (jp) {
    if (jp->ispecies == e) {
      ionambi[i] = 1;         // 1st reactant is now ion
      ip->flag = j;           // couple I ion to J electron
    } else if (ionambi[i]) {
      ionambi[i] = 0;         // 1st reactant is now neutral
      ionambi[j] = 1;         // 2nd reactant is now ion
      jp->flag = ip->flag;    // couple J ion to I reactant's electron
    } else if (ionambi[j]) {  // nothing to change
    }

  // 2 reactants become 1 product
  // ambi reaction if J reactant is electron

  } else if (!jp) {
    if (jsp == e) ionambi[i] = 0;
  }
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
   reinitialize per-cell arrays due to grid cell adaptation
   count of owned grid cells has changed
   called from adapt_grid
------------------------------------------------------------------------- */

void Collide::adapt_grid()
{
  int nglocal_old = nglocal;
  nglocal = grid->nlocal;

  // reallocate vremax and remain
  // initialize only new added locations
  // this leaves vremax/remain for non-adapted cells the same

  nglocalmax = nglocal;
  memory->grow(vremax,nglocalmax,ngroups,ngroups,"collide:vremax");
  if (remainflag)
    memory->grow(remain,nglocalmax,ngroups,ngroups,"collide:remain");

  for (int icell = nglocal_old; icell < nglocal; icell++)
    for (int igroup = 0; igroup < ngroups; igroup++)
      for (int jgroup = 0; jgroup < ngroups; jgroup++) {
        vremax[icell][igroup][jgroup] = vremax_initial[igroup][jgroup];
        if (remainflag) remain[icell][igroup][jgroup] = 0.0;
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

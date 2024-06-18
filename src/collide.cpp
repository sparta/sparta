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

#include "math.h"
#include "math_extra.h"
#include "math_eigen.h"
#include "math_eigen_impl.h"
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
#include "fix_swpm.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{NONE,DISCRETE,SMOOTH};       // several files  (NOTE: change order)
enum{PKEEP,PINSERT,PDONE,PDISCARD,PENTRY,PEXIT,PSURF};   // several files
enum{ENERGY,HEAT,STRESS};   // particle reduction choices
enum{BINARY,WEIGHT,OCTREE,OPTIMIZE}; // grouping choices

#define DELTAGRID 1000            // must be bigger than split cells per cell
#define DELTADELETE 1024
#define DELTAELECTRON 128

#define BIG 1.0e20
#define SMALLISH 1.0e-12
#define SMALL 1.0e-16

/* ---------------------------------------------------------------------- */

Collide::Collide(SPARTA *sparta, int, char **arg) : Pointers(sparta)
{
  int n = strlen(arg[0]) + 1;
  style = new char[n];
  strcpy(style,arg[0]);

  n = strlen(arg[1]) + 1;
  mixID = new char[n];
  strcpy(mixID,arg[1]);

  random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,comm->me,100);

  ngroups = 0;

  npmax = 0;
  plist = NULL;
  p2g = NULL;

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
  nearcp = 0;
  nearlimit = 10;

  recomb_ijflag = NULL;

  ambiflag = 0;
  maxelectron = 0;
  elist = NULL;

  // stochastic weighted particle method
  sweight_max = update->fnum;
  reduceflag = 0;
  Ncmin = Ncmax = Ngmin = Ngmax = 0;
  pL = NULL;
  pLU = NULL;

  // used if near-neighbor model is invoked

  max_nn = 1;
  memory->create(nn_last_partner,max_nn,"collide:nn_last_partner");
  memory->create(nn_last_partner_igroup,max_nn,"collide:nn_last_partner");
  memory->create(nn_last_partner_jgroup,max_nn,"collide:nn_last_partner");

  // initialize counters in case stats outputs them

  ncollide_one = nattempt_one = nreact_one = 0;
  ncollide_running = nattempt_running = nreact_running = 0;

  copymode = kokkos_flag = 0;
}

/* ---------------------------------------------------------------------- */

Collide::~Collide()
{
  if (copymode) return;

  delete [] style;
  delete [] mixID;
  delete random;

  memory->destroy(plist);
  memory->destroy(p2g);

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
  memory->destroy(nn_last_partner);
  memory->destroy(nn_last_partner_igroup);
  memory->destroy(nn_last_partner_jgroup);

  memory->destroy(recomb_ijflag);

  memory->destroy(pL);
  memory->destroy(pLU);
}

/* ---------------------------------------------------------------------- */

void Collide::init()
{
  // error check

  if (ambiflag && nearcp)
    error->all(FLERR,"Ambipolar collision model does not yet support "
               "near-neighbor collisions");

  // require mixture to contain all species

  int imix = particle->find_mixture(mixID);
  if (imix < 0) error->all(FLERR,"Collision mixture does not exist");
  mixture = particle->mixture[imix];

  if (mixture->nspecies != particle->nspecies)
    error->all(FLERR,"Collision mixture does not contain all species");

  if (sparta->kokkos && !kokkos_flag)
    error->all(FLERR,"Must use Kokkos-supported collision style if "
               "Kokkos is enabled");

  // if rotstyle or vibstyle = DISCRETE,
  // check that extra rotation/vibration info is defined
  // for species that require it

  if (rotstyle == DISCRETE) {
    Particle::Species *species = particle->species;
    int nspecies = particle->nspecies;

    int flag = 0;
    for (int isp = 0; isp < nspecies; isp++) {
      if (species[isp].rotdof == 0) continue;
      if (species[isp].rotdof == 2 && species[isp].nrottemp != 1) flag++;
      if (species[isp].rotdof == 3 && species[isp].nrottemp != 3) flag++;
    }
    if (flag) {
      char str[128];
      sprintf(str,"%d species do not define correct rotational "
              "temps for discrete model",flag);
      error->all(FLERR,str);
    }
  }

  if (vibstyle == DISCRETE) {
    index_vibmode = particle->find_custom((char *) "vibmode");

    Particle::Species *species = particle->species;
    int nspecies = particle->nspecies;

    int flag = 0;
    for (int isp = 0; isp < nspecies; isp++) {
      if (species[isp].vibdof <= 2) continue;
      if (index_vibmode < 0)
        error->all(FLERR,
                   "Fix vibmode must be used with discrete vibrational modes");
      if (species[isp].nvibmode != species[isp].vibdof / 2) flag++;
    }
    if (flag) {
      char str[128];
      sprintf(str,"%d species do not define correct vibrational "
              "modes for discrete model",flag);
      error->all(FLERR,str);
    }
  }

  // reallocate one-cell data structs for one or many groups

  oldgroups = ngroups;
  ngroups = mixture->ngroup;

  if (ngroups != oldgroups) {
    if (oldgroups == 1) {
      memory->destroy(plist);
      memory->destroy(pL);
      memory->destroy(pLU);
      npmax = 0;
      plist = NULL;
      pL = NULL;
      pLU = NULL;
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
      if(swpmflag) {
        memory->create(pL,npmax,"collide:pL");
        memory->create(pLU,npmax,"collide:pLU");
      }
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

  // if recombination reactions exist, set flags per species pair

  recombflag = 0;
  if (react) {
    recombflag = react->recombflag;
    recomb_boost_inverse = react->recomb_boost_inverse;
  }

  if (recombflag) {
    int nspecies = particle->nspecies;
    memory->destroy(recomb_ijflag);
    memory->create(recomb_ijflag,nspecies,nspecies,"collide:recomb_ijflag");
    for (int i = 0; i < nspecies; i++)
      for (int j = 0; j < nspecies; j++)
        recomb_ijflag[i][j] = react->recomb_exist(i,j);
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
  }

  // if ambipolar and multiple groups in mixture, ambispecies must be its own group

  if (ambiflag && mixture->ngroup > 1) {
    int *species2group = mixture->species2group;
    int egroup = species2group[ambispecies];
    if (mixture->groupsize[egroup] != 1)
      error->all(FLERR,"Multigroup ambipolar collisions require "
                 "electrons be their own group");
  }

  // find swpm fix

  if (swpmflag) {
    index_sweight = particle->find_custom((char *) "sweight");
    if (index_sweight < 0)
      error->all(FLERR,"Collision swpm without fix swpm");
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
      // not yet supported
      //else if (strcmp(arg[iarg+1],"discrete") == 0) rotstyle = DISCRETE;
      else if (strcmp(arg[iarg+1],"smooth") == 0) rotstyle = SMOOTH;
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
    } else if (strcmp(arg[iarg],"nearcp") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal collide_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) nearcp = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) nearcp = 0;
      else error->all(FLERR,"Illegal collide_modify command");
      nearlimit = atoi(arg[iarg+2]);
      if (nearcp && nearlimit <= 0)
        error->all(FLERR,"Illegal collide_modify command");
      iarg += 3;
    } else if (strcmp(arg[iarg],"swpm") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal collide_modify command");
      if (strcmp(arg[iarg+1],"no") == 0) swpmflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) swpmflag = 1;
      else error->all(FLERR,"Illegal collide_modify command");
      Ncmin = atoi(arg[iarg+2]);
      wtf = atof(arg[iarg+3]);
      if (wtf < 0) error->all(FLERR,"Illegal collide_modify command");
      iarg += 4;
    } else if (strcmp(arg[iarg],"reduce") == 0) {
      if (!swpmflag) error->all(FLERR,"Must have swpm enabled first");
      if (iarg+5 > narg) error->all(FLERR,"Illegal collide_modify command");
      reduceflag = 1;
      if (strcmp(arg[iarg+1],"energy") == 0) reduction_type = ENERGY;
      else if (strcmp(arg[iarg+1],"heat") == 0) reduction_type = HEAT;
      else if (strcmp(arg[iarg+1],"stress") == 0) reduction_type = STRESS;
      else error->all(FLERR,"Requested reduction scheme not available");
      if (reduction_type == ENERGY || reduction_type == HEAT) Ngmin = 2;
      else Ngmin = 6;
      Ncmax = atoi(arg[iarg+2]);
      if (strcmp(arg[iarg+3],"binary") == 0) group_type = BINARY;
      else if (strcmp(arg[iarg+3],"weight") == 0) group_type = WEIGHT;
      else if (strcmp(arg[iarg+3],"octree") == 0) group_type = OCTREE;
      else error->all(FLERR,"Requested grouping scheme not available");
      Ngmax = atoi(arg[iarg+4]);
      if(Ngmax < Ngmin) error->all(FLERR,"Max group size too small");
      iarg += 5;
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

  // perform collisions:
  // variant for single group or multiple groups
  // variant for nearcp flag or not
  // variant for ambipolar approximation or not
  // variant for stochastic weighted collisions or not

  if (swpmflag) {
    collisions_one_sw();
    particle->sort();
    if (reduceflag) group_reduce();
  } else if (!ambiflag) {
    if (nearcp == 0) {
      if (ngroups == 1) collisions_one<0>();
      else collisions_group<0>();
    } else {
      if (ngroups == 1) collisions_one<1>();
      else collisions_group<1>();
    }
  } else {
    if (ngroups == 1) collisions_one_ambipolar();
    else collisions_group_ambipolar();
  }

  // remove any particles deleted in chemistry reactions
  // if reactions occurred, particles are no longer sorted
  // e.g. compress_reactions may have reallocated particle->next vector

  if (ndelete) particle->compress_reactions(ndelete,dellist);
  if (react) particle->sorted = 0;
  if (swpmflag) particle->sorted = 0;

  // accumulate running totals

  nattempt_running += nattempt_one;
  ncollide_running += ncollide_one;
  nreact_running += nreact_one;
}

/* ----------------------------------------------------------------------
   NTC algorithm for a single group
------------------------------------------------------------------------- */

template < int NEARCP > void Collide::collisions_one()
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

    if (NEARCP) {
      if (np > max_nn) realloc_nn(np,nn_last_partner);
      memset(nn_last_partner,0,np*sizeof(int));
    }

    ip = cinfo[icell].first;
    volume = cinfo[icell].volume / cinfo[icell].weight;
    if (volume == 0.0) error->one(FLERR,"Collision cell volume is zero");

    // setup particle list for this cell

    if (np > npmax) {
      while (np > npmax) npmax += DELTAPART;
      memory->destroy(plist);
      memory->create(plist,npmax,"collide:plist");
    }

    n = 0;
    while (ip >= 0) {
      plist[n++] = ip;
      ip = next[ip];
    }

    // attempt = exact collision attempt count for all particles in cell
    // nattempt = rounded attempt with RN
    // if no attempts, continue to next grid cell

    attempt = attempt_collision(icell,np,volume);
    nattempt = static_cast<int> (attempt);

    if (!nattempt) continue;
    nattempt_one += nattempt;

    // perform collisions
    // select random pair of particles, cannot be same
    // test if collision actually occurs

    for (int iattempt = 0; iattempt < nattempt; iattempt++) {
      i = np * random->uniform();
      if (NEARCP) j = find_nn(i,np);
      else {
        j = np * random->uniform();
        while (i == j) j = np * random->uniform();
      }

      ipart = &particles[plist[i]];
      jpart = &particles[plist[j]];

      // test if collision actually occurs
      // continue to next collision if no reaction

      if (!test_collision(icell,0,0,ipart,jpart)) continue;

      if (NEARCP) {
        nn_last_partner[i] = j+1;
        nn_last_partner[j] = i+1;
      }

      // if recombination reaction is possible for this IJ pair
      // pick a 3rd particle to participate and set cell number density
      // unless boost factor turns it off, or there is no 3rd particle

      if (recombflag && recomb_ijflag[ipart->ispecies][jpart->ispecies]) {
        if (random->uniform() > react->recomb_boost_inverse)
          react->recomb_species = -1;
        else if (np <= 2)
          react->recomb_species = -1;
        else {
          k = np * random->uniform();
          while (k == i || k == j) k = np * random->uniform();
          react->recomb_part3 = &particles[plist[k]];
          react->recomb_species = react->recomb_part3->ispecies;
          react->recomb_density = np * update->fnum / volume;
        }
      }

      // perform collision and possible reaction

      setup_collision(ipart,jpart);
      reactflag = perform_collision(ipart,jpart,kpart);
      ncollide_one++;
      if (reactflag) nreact_one++;
      else continue;

      // if jpart destroyed: delete from plist, add particle to deletion list
      // exit attempt loop if only single particle left

      if (!jpart) {
        if (ndelete == maxdelete) {
          maxdelete += DELTADELETE;
          memory->grow(dellist,maxdelete,"collide:dellist");
        }
        dellist[ndelete++] = plist[j];
        np--;
        plist[j] = plist[np];
        if (NEARCP) nn_last_partner[j] = nn_last_partner[np];
        if (np < 2) break;
      }

      // if kpart created, add to plist
      // kpart was just added to particle list, so index = nlocal-1
      // particle data structs may have been realloced by kpart

      if (kpart) {
        if (np == npmax) {
          npmax += DELTAPART;
          memory->grow(plist,npmax,"collide:plist");
        }
        if (NEARCP) set_nn(np);
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

template < int NEARCP > void Collide::collisions_group()
{
  int i,j,k,n,ii,jj,ip,np,isp,ng;
  int pindex,ipair,igroup,jgroup,newgroup,ngmax;
  int nattempt,reactflag;
  int *ni,*nj,*ilist,*jlist;
  int *nn_igroup,*nn_jgroup;
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
    if (volume == 0.0) error->one(FLERR,"Collision cell volume is zero");

    // reallocate plist and p2g if necessary

    if (np > npmax) {
      while (np > npmax) npmax += DELTAPART;
      memory->destroy(plist);
      memory->create(plist,npmax,"collide:plist");
      memory->destroy(p2g);
      memory->create(p2g,npmax,2,"collide:p2g");
    }

    // plist = particle list for entire cell
    // glist[igroup][i] = index in plist of Ith particle in Igroup
    // ngroup[igroup] = particle count in Igroup
    // p2g[i][0] = Igroup for Ith particle in plist
    // p2g[i][1] = index within glist[igroup] of Ith particle in plist

    for (i = 0; i < ngroups; i++) ngroup[i] = 0;
    n = 0;

    while (ip >= 0) {
      isp = particles[ip].ispecies;
      igroup = species2group[isp];
      if (ngroup[igroup] == maxgroup[igroup]) {
        maxgroup[igroup] += DELTAPART;
        memory->grow(glist[igroup],maxgroup[igroup],"collide:glist");
      }
      ng = ngroup[igroup];
      glist[igroup][ng] = n;
      p2g[n][0] = igroup;
      p2g[n][1] = ng;
      plist[n] = ip;
      ngroup[igroup]++;
      n++;
      ip = next[ip];
    }

    if (NEARCP) {
      ngmax = 0;
      for (i = 0; i < ngroups; i++) ngmax = MAX(ngmax,ngroup[i]);
      if (ngmax > max_nn) {
        realloc_nn(ngmax,nn_last_partner_igroup);
        realloc_nn(ngmax,nn_last_partner_jgroup);
      }
    }

    // attempt = exact collision attempt count for a pair of groups
    // double loop over N^2 / 2 pairs of groups
    // nattempt = rounded attempt with RN
    // NOTE: not using RN for rounding of nattempt
    // gpair = list of group pairs when nattempt > 0

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
    // if chemistry occurs, exit attempt loop if group counts become too small
    // Ni and Nj are pointers to value in ngroup vector
    //   b/c need to stay current as chemistry occurs
    // NOTE: OK to use pre-computed nattempt when Ngroup may have changed via react?

    for (ipair = 0; ipair < npair; ipair++) {
      igroup = gpair[ipair][0];
      jgroup = gpair[ipair][1];
      nattempt = gpair[ipair][2];

      ni = &ngroup[igroup];
      nj = &ngroup[jgroup];
      ilist = glist[igroup];
      jlist = glist[jgroup];

      // re-test for no possible attempts
      // could have changed due to reactions in previous group pairs

      if (*ni == 0 || *nj == 0) continue;
      if (igroup == jgroup && *ni == 1) continue;

      if (NEARCP) {
        nn_igroup = nn_last_partner_igroup;
        if (igroup == jgroup) nn_jgroup = nn_last_partner_igroup;
        else nn_jgroup = nn_last_partner_jgroup;
        memset(nn_igroup,0,(*ni)*sizeof(int));
        if (igroup != jgroup) memset(nn_jgroup,0,(*nj)*sizeof(int));
      }

      for (int iattempt = 0; iattempt < nattempt; iattempt++) {
        i = *ni * random->uniform();
        if (NEARCP) j = find_nn_group(i,ilist,*nj,jlist,plist,nn_igroup,nn_jgroup);
        else {
          j = *nj * random->uniform();
          if (igroup == jgroup)
            while (i == j) j = *nj * random->uniform();
        }

        ipart = &particles[plist[ilist[i]]];
        jpart = &particles[plist[jlist[j]]];

        // test if collision actually occurs
        // continue to next collision if no reaction

        if (!test_collision(icell,igroup,jgroup,ipart,jpart)) continue;

        if (NEARCP) {
          nn_igroup[i] = j+1;
          nn_jgroup[j] = i+1;
        }

        // if recombination reaction is possible for this IJ pair
        // pick a 3rd particle to participate and set cell number density
        // unless boost factor turns it off, or there is no 3rd particle

        if (recombflag && recomb_ijflag[ipart->ispecies][jpart->ispecies]) {
          if (random->uniform() > react->recomb_boost_inverse)
            react->recomb_species = -1;
          else if (np <= 2)
            react->recomb_species = -1;
          else {
            ii = ilist[i];
            jj = jlist[j];
            k = np * random->uniform();
            while (k == ii || k == jj) k = np * random->uniform();
            react->recomb_part3 = &particles[plist[k]];
            react->recomb_species = react->recomb_part3->ispecies;
            react->recomb_density = np * update->fnum / volume;
          }
        }

        // perform collision and possible reaction

        setup_collision(ipart,jpart);
        reactflag = perform_collision(ipart,jpart,kpart);
        ncollide_one++;
        if (reactflag) nreact_one++;
        else continue;

        // ipart may now be in different group
        // reset ilist,jlist after addgroup() in case it realloced glist

        newgroup = species2group[ipart->ispecies];
        if (newgroup != igroup) {
          addgroup(newgroup,ilist[i]);
          delgroup(igroup,i);
          ilist = glist[igroup];
          jlist = glist[jgroup];
          // this line needed if jgroup=igroup and delgroup() moved J particle
          if (jgroup == igroup && j == *ni) j = i;
        }

        // jpart may now be in different group or destroyed
        // if new group: reset ilist,jlist after addgroup() in case it realloced glist
        // if destroyed: delete from plist and group, add particle to deletion list

        if (jpart) {
          newgroup = species2group[jpart->ispecies];
          if (newgroup != jgroup) {
            addgroup(newgroup,jlist[j]);
            delgroup(jgroup,j);
            ilist = glist[igroup];
            jlist = glist[jgroup];
          }

        } else {
          if (ndelete == maxdelete) {
            maxdelete += DELTADELETE;
            memory->grow(dellist,maxdelete,"collide:dellist");
          }
          pindex = jlist[j];
          dellist[ndelete++] = plist[pindex];

          delgroup(jgroup,j);

          plist[pindex] = plist[np-1];
          p2g[pindex][0] = p2g[np-1][0];
          p2g[pindex][1] = p2g[np-1][1];
          if (pindex < np-1) glist[p2g[pindex][0]][p2g[pindex][1]] = pindex;
          np--;

          if (NEARCP) nn_jgroup[j] = nn_jgroup[*nj];
        }

        // if kpart created, add to plist and group list
        // kpart was just added to particle list, so index = nlocal-1
        // reset ilist,jlist after addgroup() in case it realloced
        // particles data struct may also have been realloced

        if (kpart) {
          newgroup = species2group[kpart->ispecies];

          if (NEARCP) {
            if (newgroup == igroup || newgroup == jgroup) {
              n = ngroup[newgroup];
              set_nn_group(n);
              nn_igroup = nn_last_partner_igroup;
              if (igroup == jgroup) nn_jgroup = nn_last_partner_igroup;
              else nn_jgroup = nn_last_partner_jgroup;
              nn_igroup[n] = 0;
              nn_jgroup[n] = 0;
            }
          }

          if (np == npmax) {
            npmax += DELTAPART;
            memory->grow(plist,npmax,"collide:plist");
            memory->grow(p2g,npmax,2,"collide:p2g");
          }
          plist[np++] = particle->nlocal-1;

          addgroup(newgroup,np-1);
          ilist = glist[igroup];
          jlist = glist[jgroup];
          particles = particle->particles;
        }

        // test to exit attempt loop due to groups becoming too small

        if (*ni <= 1) {
          if (*ni == 0) break;
          if (igroup == jgroup) break;
        }
        if (*nj <= 1) {
          if (*nj == 0) break;
          if (igroup == jgroup) break;
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
  int i,j,k,n,ip,np,nelectron,nptotal,jspecies,tmp;
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
    if (volume == 0.0) error->one(FLERR,"Collision cell volume is zero");

    // setup particle list for this cell

    if (np > npmax) {
      while (np > npmax) npmax += DELTAPART;
      memory->destroy(plist);
      memory->create(plist,npmax,"collide:plist");
    }

    n = 0;
    while (ip >= 0) {
      plist[n++] = ip;
      ip = next[ip];
    }

    // setup elist of ionized electrons for this cell
    // create them in separate array since will never become real particles

    if (np >= maxelectron) {
      while (maxelectron < np) maxelectron += DELTAELECTRON;
      memory->sfree(elist);
      elist = (Particle::OnePart *)
        memory->smalloc(maxelectron*nbytes,"collide:elist");
    }

    // create electrons for ambipolar ions

    nelectron = 0;
    for (i = 0; i < np; i++) {
      if (ionambi[plist[i]]) {
        p = &particles[plist[i]];
        ep = &elist[nelectron];
        memcpy(ep,p,nbytes);
        memcpy(ep->v,velambi[plist[i]],3*sizeof(double));
        ep->ispecies = ambispecies;
        nelectron++;
      }
    }

    // attempt = exact collision attempt count for all particles in cell
    // nptotal = includes neutrals, ions, electrons
    // nattempt = rounded attempt with RN

    nptotal = np + nelectron;
    attempt = attempt_collision(icell,nptotal,volume);
    nattempt = static_cast<int> (attempt);

    if (!nattempt) continue;
    nattempt_one += nattempt;

    // perform collisions
    // select random pair of particles, cannot be same
    // test if collision actually occurs
    // if chemistry occurs, exit attempt loop if group count goes to 0

    for (int iattempt = 0; iattempt < nattempt; iattempt++) {
      i = nptotal * random->uniform();
      j = nptotal * random->uniform();
      while (i == j) j = nptotal * random->uniform();

      // ipart,jpart = heavy particles or electrons

      if (i < np) ipart = &particles[plist[i]];
      else ipart = &elist[i-np];
      if (j < np) jpart = &particles[plist[j]];
      else jpart = &elist[j-np];

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

      // test if collision actually occurs

      if (!test_collision(icell,0,0,ipart,jpart)) continue;

      // if recombination reaction is possible for this IJ pair
      // pick a 3rd particle to participate and set cell number density
      // unless boost factor turns it off, or there is no 3rd particle
      // 3rd particle cannot be an electron, so select from Np

      if (recombflag && recomb_ijflag[ipart->ispecies][jpart->ispecies]) {
        if (random->uniform() > react->recomb_boost_inverse)
          react->recomb_species = -1;
        else if (np == 1)
          react->recomb_species = -1;
        else if (np == 2 && jpart->ispecies != ambispecies)
          react->recomb_species = -1;
        else {
          k = np * random->uniform();
          while (k == i || k == j) k = np * random->uniform();
          react->recomb_part3 = &particles[plist[k]];
          react->recomb_species = react->recomb_part3->ispecies;
          react->recomb_density = np * update->fnum / volume;
        }
      }

      // perform collision
      // ijspecies = species before collision chemistry
      // continue to next collision if no reaction

      jspecies = jpart->ispecies;
      setup_collision(ipart,jpart);
      reactflag = perform_collision(ipart,jpart,kpart);
      ncollide_one++;
      if (reactflag) nreact_one++;
      else continue;

      // reset ambipolar ion flags due to collision
      // must do now before particle count reset below can break out of loop
      // first reset ionambi if kpart was added since ambi_reset() uses it

      if (kpart) ionambi = particle->eivec[particle->ewhich[index_ionambi]];
      if (jspecies == ambispecies)
        ambi_reset(plist[i],-1,jspecies,ipart,jpart,kpart,ionambi);
      else
        ambi_reset(plist[i],plist[j],jspecies,ipart,jpart,kpart,ionambi);

      // if kpart created:
      // particles and custom data structs may have been realloced by kpart
      // add kpart to plist or elist
      // kpart was just added to particle list, so index = nlocal-1
      // must come before jpart code below since it modifies nlocal

      if (kpart) {
        particles = particle->particles;
        ionambi = particle->eivec[particle->ewhich[index_ionambi]];
        velambi = particle->edarray[particle->ewhich[index_velambi]];

        if (kpart->ispecies != ambispecies) {
          if (np == npmax) {
            npmax += DELTAPART;
            memory->grow(plist,npmax,"collide:plist");
          }
          plist[np++] = particle->nlocal-1;

        } else {
          if (nelectron == maxelectron) {
            maxelectron += DELTAELECTRON;
            elist = (Particle::OnePart *)
              memory->srealloc(elist,maxelectron*nbytes,"collide:elist");
          }
          ep = &elist[nelectron];
          memcpy(ep,kpart,nbytes);
          ep->ispecies = ambispecies;
          nelectron++;
          particle->nlocal--;
        }
      }

      // if jpart exists, was originally not an electron, now is an electron:
      //   ionization reaction converted 2 neutrals to one ion
      //   add to elist, remove from plist, flag J for deletion
      // if jpart exists, was originally an electron, now is not an electron:
      //   exchange reaction converted ion + electron to two neutrals
      //   add neutral J to master particle list, remove from elist, add to plist
      // if jpart destroyed, was an electron:
      //   recombination reaction converted ion + electron to one neutral
      //   remove electron from elist
      // else if jpart destroyed:
      //   non-ambipolar recombination reaction
      //   remove from plist, flag J for deletion

      if (jpart) {
        if (jspecies != ambispecies && jpart->ispecies == ambispecies) {
          if (nelectron == maxelectron) {
            maxelectron += DELTAELECTRON;
            elist = (Particle::OnePart *)
              memory->srealloc(elist,maxelectron*nbytes,"collide:elist");
          }
          ep = &elist[nelectron];
          memcpy(ep,jpart,nbytes);
          ep->ispecies = ambispecies;
          nelectron++;
          jpart = NULL;

        } else if (jspecies == ambispecies && jpart->ispecies != ambispecies) {
          int reallocflag = particle->add_particle();
          if (reallocflag) {
            particles = particle->particles;
            ionambi = particle->eivec[particle->ewhich[index_ionambi]];
            velambi = particle->edarray[particle->ewhich[index_velambi]];
          }

          int index = particle->nlocal-1;
          memcpy(&particles[index],jpart,nbytes);
          particles[index].id = MAXSMALLINT*random->uniform();
          ionambi[index] = 0;

          if (nelectron-1 != j-np) memcpy(&elist[j-np],&elist[nelectron-1],nbytes);
          nelectron--;

          if (np == npmax) {
            npmax += DELTAPART;
            memory->grow(plist,npmax,"collide:plist");
          }
          plist[np++] = index;
        }
      }

      if (!jpart && jspecies == ambispecies) {
        if (nelectron-1 != j-np) memcpy(&elist[j-np],&elist[nelectron-1],nbytes);
        nelectron--;

      } else if (!jpart) {
        if (ndelete == maxdelete) {
          maxdelete += DELTADELETE;
          memory->grow(dellist,maxdelete,"collide:dellist");
        }
        dellist[ndelete++] = plist[j];
        plist[j] = plist[np-1];
        np--;
      }

      // update particle counts
      // quit if no longer enough particles for another collision

      nptotal = np + nelectron;
      if (nptotal < 2) break;
    }

    // done with collisions/chemistry for one grid cell
    // recombine ambipolar ions with their matching electrons
    //   by copying electron velocity into velambi
    // which ion is combined with which electron does not matter
    // error if ion count does not match electron count

    int melectron = 0;
    for (n = 0; n < np; n++) {
      i = plist[n];
      if (ionambi[i]) {
        if (melectron < nelectron) {
          ep = &elist[melectron];
          memcpy(velambi[i],ep->v,3*sizeof(double));
        }
        melectron++;
      }
    }
    if (melectron != nelectron)
      error->one(FLERR,"Collisions in cell did not conserve electron count");
  }
}

/* ----------------------------------------------------------------------
   NTC algorithm for multiple groups with ambipolar approximation
   loop over pairs of groups, pre-compute # of attempts per group pair
------------------------------------------------------------------------- */

void Collide::collisions_group_ambipolar()
{
  int i,j,k,n,ii,jj,ip,np,isp,ng;
  int pindex,ipair,igroup,jgroup,newgroup,jspecies,tmp;
  int nattempt,reactflag,nelectron;
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
    if (volume == 0.0) error->one(FLERR,"Collision cell volume is zero");

    // reallocate plist and p2g if necessary

    if (np > npmax) {
      while (np > npmax) npmax += DELTAPART;
      memory->destroy(plist);
      memory->create(plist,npmax,"collide:plist");
      memory->destroy(p2g);
      memory->create(p2g,npmax,2,"collide:p2g");
    }

    // setup elist of ionized electrons for this cell
    // create them in separate array since will never become real particles

    if (np >= maxelectron) {
      while (maxelectron < np) maxelectron += DELTAELECTRON;
      memory->sfree(elist);
      elist = (Particle::OnePart *)
        memory->smalloc(maxelectron*nbytes,"collide:elist");
    }

    // plist = particle list for entire cell
    // glist[igroup][i] = index in plist of Ith particle in Igroup
    // ngroup[igroup] = particle count in Igroup
    // p2g[i][0] = Igroup for Ith particle in plist
    // p2g[i][1] = index within glist[igroup] of Ith particle in plist
    // also populate elist with ionized electrons, now separated from ions
    // ngroup[egroup] = nelectron

    for (i = 0; i < ngroups; i++) ngroup[i] = 0;
    n = 0;
    nelectron = 0;

    while (ip >= 0) {
      isp = particles[ip].ispecies;
      igroup = species2group[isp];
      if (ngroup[igroup] == maxgroup[igroup]) {
        maxgroup[igroup] += DELTAPART;
        memory->grow(glist[igroup],maxgroup[igroup],"collide:glist");
      }
      ng = ngroup[igroup];
      glist[igroup][ng] = n;
      p2g[n][0] = igroup;
      p2g[n][1] = ng;
      plist[n] = ip;
      ngroup[igroup]++;

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
        ng = ngroup[egroup];
        glist[egroup][ng] = nelectron-1;
        ngroup[egroup]++;
      }

      n++;
      ip = next[ip];
    }

    // attempt = exact collision attempt count for a pair of groups
    // double loop over N^2 / 2 pairs of groups
    // temporarily include nelectrons in count for egroup
    // nattempt = rounded attempt with RN
    // NOTE: not using RN for rounding of nattempt
    // gpair = list of group pairs when nattempt > 0
    //         flip igroup/jgroup if igroup = egroup
    // egroup/egroup collisions are not included in gpair

    npair = 0;
    for (igroup = 0; igroup < ngroups; igroup++)
      for (jgroup = igroup; jgroup < ngroups; jgroup++) {
        if (igroup == egroup && jgroup == egroup) continue;
        attempt = attempt_collision(icell,igroup,jgroup,volume);
        nattempt = static_cast<int> (attempt);

        if (nattempt) {
          if (igroup == egroup) {
              gpair[npair][0] = jgroup;
              gpair[npair][1] = igroup;
            } else {
              gpair[npair][0] = igroup;
              gpair[npair][1] = jgroup;
            }
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
    // if chemistry occurs, exit attempt loop if group counts become too small
    // Ni and Nj are pointers to value in ngroup vector
    //   b/c need to stay current as chemistry occurs
    // NOTE: OK to use pre-computed nattempt when Ngroup may have changed via react?

    for (ipair = 0; ipair < npair; ipair++) {
      igroup = gpair[ipair][0];
      jgroup = gpair[ipair][1];
      nattempt = gpair[ipair][2];

      ni = &ngroup[igroup];
      nj = &ngroup[jgroup];
      ilist = glist[igroup];
      jlist = glist[jgroup];

      // re-test for no possible attempts
      // could have changed due to reactions in previous group pairs

      if (*ni == 0 || *nj == 0) continue;
      if (igroup == jgroup && *ni == 1) continue;

      for (int iattempt = 0; iattempt < nattempt; iattempt++) {
        i = *ni * random->uniform();
        j = *nj * random->uniform();
        if (igroup == jgroup)
          while (i == j) j = *nj * random->uniform();

        // ipart/jpart can be from particles or elist

        if (igroup == egroup) ipart = &elist[i];
        else ipart = &particles[plist[ilist[i]]];
        if (jgroup == egroup) jpart = &elist[j];
        else jpart = &particles[plist[jlist[j]]];

        // NOTE: unlike single group, no possibility of e/e collision
        //       means collision stats may be different

        //if (ipart->ispecies == ambispecies && jpart->ispecies == ambispecies) {
        //  ncollide_one++;
        //  continue;
        //}

        // test if collision actually occurs

        if (!test_collision(icell,igroup,jgroup,ipart,jpart)) continue;

        // if recombination reaction is possible for this IJ pair
        // pick a 3rd particle to participate and set cell number density
        // unless boost factor turns it off, or there is no 3rd particle
        // 3rd particle will never be an electron since plist has no electrons
        // if jgroup == egroup, no need to check k for match to jj

        if (recombflag && recomb_ijflag[ipart->ispecies][jpart->ispecies]) {
          if (random->uniform() > react->recomb_boost_inverse)
            react->recomb_species = -1;
          else if (np <= 2)
            react->recomb_species = -1;
          else {
            ii = ilist[i];
            if (jgroup == egroup) jj = -1;
            else jj = jlist[j];
            k = np * random->uniform();
            while (k == ii || k == jj) k = np * random->uniform();
            react->recomb_part3 = &particles[plist[k]];
            react->recomb_species = react->recomb_part3->ispecies;
            react->recomb_density = np * update->fnum / volume;
          }
        }

        // perform collision
        // ijspecies = species before collision chemistry
        // continue to next collision if no reaction

        jspecies = jpart->ispecies;
        setup_collision(ipart,jpart);
        reactflag = perform_collision(ipart,jpart,kpart);
        ncollide_one++;
        if (reactflag) nreact_one++;
        else continue;

        // reset ambipolar ion flags due to reaction
        // must do now before group reset below can break out of loop
        // first reset ionambi if kpart was added since ambi_reset() uses it

        if (kpart) ionambi = particle->eivec[particle->ewhich[index_ionambi]];
        if (jgroup == egroup)
          ambi_reset(plist[ilist[i]],-1,jspecies,ipart,jpart,kpart,ionambi);
        else
          ambi_reset(plist[ilist[i]],plist[jlist[j]],jspecies,
                     ipart,jpart,kpart,ionambi);

        // ipart may now be in different group
        // reset ilist,jlist after addgroup() in case it realloced glist

        newgroup = species2group[ipart->ispecies];
        if (newgroup != igroup) {
          addgroup(newgroup,ilist[i]);
          delgroup(igroup,i);
          ilist = glist[igroup];
          jlist = glist[jgroup];
          // this line needed if jgroup=igroup and delgroup() moved J particle
          if (jlist == ilist && j == *ni) j = i;
        }

        // if kpart created:
        // particles and custom data structs may have been realloced by kpart
        // add kpart to plist or elist and to group
        // kpart was just added to particle list, so index = nlocal-1
        // must come before jpart code below since it modifies nlocal

        if (kpart) {
          particles = particle->particles;
          ionambi = particle->eivec[particle->ewhich[index_ionambi]];
          velambi = particle->edarray[particle->ewhich[index_velambi]];

          newgroup = species2group[kpart->ispecies];

          if (newgroup != egroup) {
            if (np == npmax) {
              npmax += DELTAPART;
              memory->grow(plist,npmax,"collide:plist");
              memory->grow(p2g,npmax,2,"collide:p2g");
            }
            plist[np++] = particle->nlocal-1;
            addgroup(newgroup,np-1);
            ilist = glist[igroup];
            jlist = glist[jgroup];

          } else {
            if (nelectron == maxelectron) {
              maxelectron += DELTAELECTRON;
              elist = (Particle::OnePart *)
                memory->srealloc(elist,maxelectron*nbytes,"collide:elist");
            }
            ep = &elist[nelectron];
            memcpy(ep,kpart,nbytes);
            ep->ispecies = ambispecies;
            nelectron++;
            particle->nlocal--;

            if (ngroup[egroup] == maxgroup[egroup]) {
              maxgroup[egroup] += DELTAPART;
              memory->grow(glist[egroup],maxgroup[egroup],"collide:grouplist");
            }
            ng = ngroup[egroup];
            glist[egroup][ng] = nelectron-1;
            ngroup[egroup]++;
          }
        }

        // jpart may now be in a different group or destroyed
        // if jpart exists, now in a different group, neither group is egroup:
        //   add/del group, reset ilist,jlist after addgroup() in case glist realloced
        // if jpart exists, was originally not an electron, now is an electron:
        //   ionization reaction converted 2 neutrals to one ion
        //   add to elist, remove from plist, flag J for deletion
        // if jpart exists, was originally an electron, now is not an electron:
        //   exchange reaction converted ion + electron to two neutrals
        //   add neutral J to master particle list, remove from elist, add to plist
        // if jpart destroyed, was an electron:
        //   recombination reaction converted ion + electron to one neutral
        //   remove electron from elist
        // else if jpart destroyed:
        //   non-ambipolar recombination reaction
        //   remove from plist and group, add particle to deletion list

        if (jpart) {
          newgroup = species2group[jpart->ispecies];

          if (newgroup == jgroup) {
            // nothing to do

          } else if (jgroup != egroup && newgroup != egroup) {
            addgroup(newgroup,jlist[j]);
            delgroup(jgroup,j);
            ilist = glist[igroup];
            jlist = glist[jgroup];

          } else if (jgroup != egroup && jpart->ispecies == ambispecies) {
            if (nelectron == maxelectron) {
              maxelectron += DELTAELECTRON;
              elist = (Particle::OnePart *)
                memory->srealloc(elist,maxelectron*nbytes,"collide:elist");
            }
            ep = &elist[nelectron];
            memcpy(ep,jpart,nbytes);
            ep->ispecies = ambispecies;
            nelectron++;

            if (ngroup[egroup] == maxgroup[egroup]) {
              maxgroup[egroup] += DELTAPART;
              memory->grow(glist[egroup],maxgroup[egroup],"collide:grouplist");
            }
            ng = ngroup[egroup];
            glist[egroup][ng] = nelectron-1;
            ngroup[egroup]++;

            jpart = NULL;

          } else if (jgroup == egroup && jpart->ispecies != ambispecies) {
            int reallocflag = particle->add_particle();
            if (reallocflag) {
              particles = particle->particles;
              ionambi = particle->eivec[particle->ewhich[index_ionambi]];
              velambi = particle->edarray[particle->ewhich[index_velambi]];
            }

            int index = particle->nlocal-1;
            memcpy(&particles[index],jpart,nbytes);
            particles[index].id = MAXSMALLINT*random->uniform();
            ionambi[index] = 0;

            if (nelectron-1 != j) memcpy(&elist[j],&elist[nelectron-1],nbytes);
            nelectron--;
            ngroup[egroup]--;

            if (np == npmax) {
              npmax += DELTAPART;
              memory->grow(plist,npmax,"collide:plist");
              memory->grow(p2g,npmax,2,"collide:p2g");
            }
            plist[np++] = index;
            addgroup(newgroup,np-1);
            ilist = glist[igroup];
            jlist = glist[jgroup];
          }
        }

        if (!jpart && jspecies == ambispecies) {
          if (nelectron-1 != j) memcpy(&elist[j],&elist[nelectron-1],nbytes);
          nelectron--;
          ngroup[egroup]--;

        } else if (!jpart) {
          if (ndelete == maxdelete) {
            maxdelete += DELTADELETE;
            memory->grow(dellist,maxdelete,"collide:dellist");
          }
          pindex = jlist[j];
          dellist[ndelete++] = plist[pindex];

          delgroup(jgroup,j);

          plist[pindex] = plist[np-1];
          p2g[pindex][0] = p2g[np-1][0];
          p2g[pindex][1] = p2g[np-1][1];
          if (pindex < np-1) glist[p2g[pindex][0]][p2g[pindex][1]] = pindex;
          np--;
        }

        // test to exit attempt loop due to groups becoming too small

        if (*ni <= 1) {
          if (*ni == 0) break;
          if (igroup == jgroup) break;
        }
        if (*nj <= 1) {
          if (*nj == 0) break;
          if (igroup == jgroup) break;
        }
      }
    }

    // done with collisions/chemistry for one grid cell
    // recombine ambipolar ions with their matching electrons
    //   by copying electron velocity into velambi
    // which ion is combined with which electron does not matter
    // error if do not use all nelectrons in cell

    int melectron = 0;
    for (n = 0; n < np; n++) {
      i = plist[n];
      if (ionambi[i]) {
        if (melectron < nelectron) {
          ep = &elist[melectron];
          memcpy(velambi[i],ep->v,3*sizeof(double));
        }
        melectron++;
      }
    }
    if (melectron != nelectron)
      error->one(FLERR,"Collisions in cell did not conserve electron count");
  }
}

/* ----------------------------------------------------------------------
   reset ionambi flags if ambipolar reaction occurred
   this operates independent of cell particle counts and plist/elist data structs
     caller will adjust those after this method returns
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
        set I product = ion
     all other 2 -> 3 cases, set K product = neutral
   check for 4 versions of 2 -> 2: ionization or exchange
     I: A + B -> AB+ + e
        if J product = electron:
        set I product to ion
     E: AB+ + e -> A + B
        if I reactant = ion and J reactant = elecrton
        set I/J products to neutral
     E: AB+ + C -> A + BC+
        if I reactant = ion:
        set I/J products to neutral/ion
     E: C + AB+ -> A + BC+
        if J reactant = ion:
        nothing to change for products
     all other 2 -> 2 cases, no changes
   check for one version of 2 -> 1: recombination
     R: A+ + e -> A
        if ej = elec, set I product to neutral
     all other 2 -> 1 cases, no changes
   WARNING:
     do not index by I,J if could be e, since may be negative I,J index
     do not access ionambi if could be e, since e may be in elist
------------------------------------------------------------------------- */

void Collide::ambi_reset(int i, int j, int jsp,
                         Particle::OnePart *ip, Particle::OnePart *jp,
                         Particle::OnePart *kp, int *ionambi)
{
  int e = ambispecies;

  // 2 reactants become 3 products
  // in all ambi reactions with an electron reactant, it is J

  if (kp) {
    int k = particle->nlocal-1;
    ionambi[k] = 0;
    if (jsp != e) return;

    if (ionambi[i]) {                // nothing to change
    } else if (kp->ispecies == e) {
      ionambi[i] = 1;                // 1st reactant is now 1st product ion
    }

  // 2 reactants become 2 products
  // ambi reaction if J product is electron or either reactant is ion

  } else if (jp) {
    if (jp->ispecies == e) {
      ionambi[i] = 1;         // 1st reactant is now 1st product ion
    } else if (ionambi[i] && jsp == e) {
      ionambi[i] = 0;         // 1st reactant is now 1st product neutral
    } else if (ionambi[i]) {
      ionambi[i] = 0;         // 1st reactant is now 1st product neutral
      ionambi[j] = 1;         // 2nd reactant is now 2nd product ion
    }

  // 2 reactants become 1 product
  // ambi reaction if J reactant is electron

  } else if (!jp) {
    if (jsp == e) ionambi[i] = 0;   // 1st reactant is now 1st product neutral
  }
}

/* ----------------------------------------------------------------------
   DEBUG : Checks particles
------------------------------------------------------------------------- */

/*template <int ORD> void Collide::check_particles()
{
  //if(ORD==0) printf("before collide\n");
  //else if(ORD==1) printf("after collide\n");
  //else if(ORD==2) printf("after sort\n");
  //else if(ORD==3) printf("after reduce\n");
  //else if(ORD==4) printf("after compress\n");

  int i,j,m,n,ip,np;
  int nattempt,reactflag;
  double attempt,volume;
  Particle::OnePart *ipart;
  Grid::ChildInfo *cinfo = grid->cinfo;
  Particle::OnePart *particles = particle->particles;
  int *next = particle->next;

  double rho;
  double E;
  double V[3], vsq;
  E = rho = 0.0;
  V[0] = V[1] = V[2] = 0.0;
  for (int icell = 0; icell < nglocal; icell++) {
    // create particle list
    ip = cinfo[icell].first;
    while (ip >= 0) {
      ipart = &particles[ip];
      if(ipart->sweight > 0) {
        rho += ipart->sweight;
        vsq = ipart->v[0]*ipart->v[0]+
              ipart->v[1]*ipart->v[1]+
              ipart->v[2]*ipart->v[2];
        for (int d = 0; d < 3; d++) V[d] += ipart->sweight*ipart->v[d];
        E += ipart->sweight*vsq;
      }
      ip = next[ip];
    }
  }

  if(update->ntimestep == 1) {
    rho0 = rho;
    for(int d = 0; d < 3; d++) V0[d] = V[d];
    E0 = E;
  }

  if(update->ntimestep % 1000 == 0)
    printf("drho: %2.8e; dE: %2.8e; dV: %2.8e,%2.8e,%2.8e\n",
          1-rho/rho0,
          1-E/E0,
          1-V[0]/V0[0],
          1-V[1]/V0[1],
          1-V[2]/V0[2]);
  return;
}*/

/* ----------------------------------------------------------------------
   Stochastic weighted algorithm
------------------------------------------------------------------------- */

void Collide::collisions_one_sw()
{
  int i,j,n,ip,np,newp;
  int nattempt;
  double attempt,volume;
  Particle::OnePart *ipart,*jpart,*kpart,*lpart,*mpart;

  int ilevel;
  double np_scale;

  // loop over cells I own

  Grid::ChildInfo *cinfo = grid->cinfo;
  Grid::ChildCell *cells = grid->cells;

  Particle::OnePart *particles = particle->particles;
  int *next = particle->next;

  double isw;
  double *sweights = particle->edvec[particle->ewhich[index_sweight]];

  for (int icell = 0; icell < nglocal; icell++) {
    np = cinfo[icell].count;
    if (np <= 1) continue;

    volume = cinfo[icell].volume / cinfo[icell].weight;
    if (volume == 0.0) error->one(FLERR,"Collision cell volume is zero");

    // setup particle list for this cell

    if (np > npmax) {
      while (np > npmax) npmax += DELTAPART;
      memory->destroy(plist);
      memory->create(plist,npmax,"collide:plist");
      memory->create(pL,npmax,"collide:pL");
      memory->create(pLU,npmax,"collide:pLU");
    }

    // get stochastic weights

    // build particle list and find maximum particle weight

    ip = cinfo[icell].first;
    n = 0;
    sweight_max = 0.0;
    while (ip >= 0) {
      plist[n++] = ip;
      ipart = &particles[ip];
      isw = sweights[ip];
      sweight_max = MAX(sweight_max,isw);

      if (isw != isw) error->all(FLERR,"Particle has NaN weight");
      if (isw <= 0.0) error->all(FLERR,"Particle has negative or zero weight");
      ip = next[ip];
    }

    ilevel = cells[icell].level;
    if(ilevel==1) np_scale = 1.0;
    else np_scale = pow(8,ilevel-1);

    // attempt = exact collision attempt count for all particles in cell
    // nattempt = rounded attempt with RN
    // if no attempts, continue to next grid cell

    attempt = attempt_collision(icell,np,volume);
    nattempt = static_cast<int> (attempt);

    if (!nattempt) continue;
    nattempt_one += nattempt;

    for (int iattempt = 0; iattempt < nattempt; iattempt++) {

      i = np * random->uniform();
      j = np * random->uniform();
      while (i == j) j = np * random->uniform();

      ipart = &particles[plist[i]];
      jpart = &particles[plist[j]];

      if (!test_collision(icell,0,0,ipart,jpart)) continue;

      // check if pair has zero or negative weight

      if(sweights[i] <= 0.0 || sweights[j] <= 0.0) {
        printf("i: %i -- %4.3e; j: %i -- %4.3e\n", i, sweights[i], j, sweights[j]);
        error->one(FLERR,"bad weight before split");
      }

      // split particles
      if (np >= Ncmin/np_scale && Ncmin > 0.0) pre_wtf = 0.0;
      else pre_wtf = 1.0;

      newp = split(ipart,jpart,kpart,lpart);

      // add new particles to particle list

      if (newp > 1) {
        particles = particle->particles;
        sweights = particle->edvec[particle->ewhich[index_sweight]];
        if (np+2 >= npmax) {
          npmax += DELTAPART;
          memory->grow(plist,npmax,"collide:plist");
          memory->grow(pL,npmax,"collide:pL");
          memory->grow(pLU,npmax,"collide:pLU");
        }
        plist[np++] = particle->nlocal-2;
        plist[np++] = particle->nlocal-1;
        particles = particle->particles;
      } else if (newp > 0) {
        particles = particle->particles;
        sweights = particle->edvec[particle->ewhich[index_sweight]];
        if (np+1 >= npmax) {
          npmax += DELTAPART;
          memory->grow(plist,npmax,"collide:plist");
          memory->grow(pL,npmax,"collide:pL");
          memory->grow(pLU,npmax,"collide:pLU");
        }
        plist[np++] = particle->nlocal-1;
        particles = particle->particles;
      }

      // since ipart and jpart have same weight, do not need
      // ... to account for weight during collision itself
      // also the splits are all handled beforehand

      mpart = NULL; // dummy particle
      setup_collision(ipart,jpart);
      perform_collision(ipart,jpart,mpart);
      ncollide_one++;

    } // end attempt loop

    // removes very small weighted particles
    //remove_small();

  } // loop for cells
}

/* ----------------------------------------------------------------------
   Splits particles and generates two new particles (for SWPM)
------------------------------------------------------------------------- */

int Collide::split(Particle::OnePart *&ip, Particle::OnePart *&jp,
                   Particle::OnePart *&kp, Particle::OnePart *&lp)
{
  double xk[3],vk[3];
  double xl[3],vl[3];
  double erotk, erotl;
  int ks, ls;
  int kcell, lcell;

  // checks if particles properly deleted

  int id;
  int reallocflag;
  Particle::OnePart *particles = particle->particles;

  kp = NULL;
  lp = NULL;

  // weight transfer function is assumed to be
  // ... MIN(ip->sweight,jp->sweight)/(1 + pre_wtf * wtf)

  double *sweights = particle->edvec[particle->ewhich[index_sweight]];
  double isw = sweights[ip-particle->particles];
  double jsw = sweights[jp-particle->particles];
  double Gwtf, ksw, lsw;

  if (isw <= 0.0 || jsw <= 0.0)
    error->one(FLERR,"Zero or negative weight before split");

  // particle ip has larger weight

  if(isw >= jsw) {
    Gwtf = jsw/(1.0+pre_wtf*wtf);
    ksw  = isw-Gwtf;
    lsw  = jsw-Gwtf;

    ks = ip->ispecies;
    ls = jp->ispecies;

    kcell = ip->icell;
    lcell = jp->icell;

    memcpy(xk,ip->x,3*sizeof(double));
    memcpy(vk,ip->v,3*sizeof(double));
    memcpy(xl,jp->x,3*sizeof(double));
    memcpy(vl,jp->v,3*sizeof(double));

    erotk = ip->erot;
    erotl = jp->erot;

  // particle jp has larger weight

  } else {
    Gwtf = isw/(1.0+pre_wtf*wtf);
    ksw  = jsw-Gwtf;
    lsw  = isw-Gwtf;

    ks = jp->ispecies;
    ls = ip->ispecies;

    kcell = jp->icell;
    lcell = ip->icell;

    memcpy(xk,jp->x,3*sizeof(double));
    memcpy(vk,jp->v,3*sizeof(double));
    memcpy(xl,ip->x,3*sizeof(double));
    memcpy(vl,ip->v,3*sizeof(double));

    erotk = jp->erot;
    erotl = ip->erot;
  }

  // Gwtf should never be negative or zero

  if (Gwtf <= 0.0)
    error->one(FLERR,"Negative weight assigned after split");

  if (Gwtf > 0.0 && pre_wtf > 0.0)
    if (ksw <= 0.0 || lsw <= 0.0)
      error->one(FLERR,"Zero or negative weight after split");

  // number of new particles

  int newp = 0;

  // gk is always the bigger of the two

  if(ksw > 0) {
    id = MAXSMALLINT*random->uniform();
    reallocflag = particle->add_particle(id,ks,kcell,xk,vk,erotk,0.0);
    if (reallocflag) {
      sweights = particle->edvec[particle->ewhich[index_sweight]];
      ip = particle->particles + (ip - particles);
      jp = particle->particles + (jp - particles);
    }
    kp = &particle->particles[particle->nlocal-1];
    sweights[particle->nlocal-1] = ksw;
    newp++;
  }

  // there should never be case where you add particle "l" if
  // ... you did not add particle "k"

  if(lsw > 0) {
    if(ksw <= 0) error->one(FLERR,"Bad addition to particle list");
    id = MAXSMALLINT*random->uniform();
    reallocflag = particle->add_particle(id,ls,lcell,xl,vl,erotl,0.0);
    if (reallocflag) {
      sweights = particle->edvec[particle->ewhich[index_sweight]];
      ip = particle->particles + (ip - particles);
      jp = particle->particles + (jp - particles);
      kp = particle->particles + (kp - particles);
    }
    sweights[particle->nlocal-1] = lsw;
    newp++;
  }

  // update weights

  sweights[ip - particle->particles] = Gwtf;
  sweights[jp - particle->particles] = Gwtf;

  return newp;
}

/* ----------------------------------------------------------------------
   Reorder plist depending on grouping strategy used and prepare for
   grouping and particle reduction
------------------------------------------------------------------------- */

void Collide::group_reduce()
{
  int n,nold,np,ip;
  double isw;

  Particle::OnePart *ipart;
  Grid::ChildInfo *cinfo = grid->cinfo;
  Grid::ChildCell *cells = grid->cells;
  Particle::OnePart *particles = particle->particles;
  int *next = particle->next;
  int ilevel, nthresh;
  double np_scale;

  double *sweights = particle->edvec[particle->ewhich[index_sweight]];
  double swmean, swvar, swstd;
  double d1, d2;
  double lLim, uLim;
  int npL, npLU;

  for (int icell = 0; icell < nglocal; icell++) {
    np = cinfo[icell].count;

    ilevel = cells[icell].level;
    if(ilevel == 1) nthresh = Ncmax;
    else {
      np_scale = pow(8,ilevel-1);
      // number of particles should at least exceed minimum group size
      nthresh = MAX(Ncmax/np_scale, Ngmax);
    }

    if (np <= nthresh) continue;

    // create particle list

    ip = cinfo[icell].first;
    n = 0;
    while (ip >= 0) {
      ipart = &particles[ip];
      isw = sweights[ip];
      if(isw > 0) plist[n++] = ip;
      ip = next[ip];
    }

    gbuf = 0;

    while (n > nthresh) {
      nold = n;

      // seems to be more stable than weighted

      if (group_type == BINARY) {

        // shuffle indices to choose random positions

        /*int j;
        for (int i = n-1; i > 0; --i) {
          j = random->uniform()*i;
          if (j < 0 || j >= n) error->one(FLERR,"bad index");
          std::swap(plist[i], plist[j]);
        }*/

        group_bt(plist,n);

      } else if (group_type == WEIGHT) {

        // find mean / standard deviation of weight

        ip = cinfo[icell].first;
        n = 0;
        swmean = swvar = 0.0;
        while (ip >= 0) {
          ipart = &particles[ip];
          isw = sweights[ip];

          // Incremental variance

          if(isw > 0) {
            n++;
            d1 = isw - swmean;
            swmean += (d1/n);
            swvar += (n-1.0)/n*d1*d1;
          }
          ip = next[ip];
        }
        swstd = sqrt(swvar/n);

        // weight limits to separate particles

        lLim = MAX(swmean-1.25*swstd,0);
        uLim = swmean+2.0*swstd;

        // recreate particle list and omit large weighted particles

        ip = cinfo[icell].first;
        npL = npLU = 0;
        while (ip >= 0) {
          ipart = &particles[ip];
          isw = sweights[ip];
          if(isw > 0 && isw < lLim) pL[npL++] = ip;
          else if(isw >= lLim && isw < uLim) pLU[npLU++] = ip;
          ip = next[ip];
        }

        // shuffle indices to choose random positions

        /*int j;
        for(int i = n-1; i > 0; --i) {
          j = random->uniform()*i;
          if(j < 0 || j >= n) error->one(FLERR,"bad index");
          std::swap(plist[i], plist[j]);
        }*/

        // rearrange so that small weighted particles in front

        /*int pmid = 0;
        for(int i = 0; i < n; i++) {
          ipart = &particles[plist[i]];
          isw = sweights[plist[i]];
          if(isw < lLim)
            std::swap(plist[pmid++],plist[i]);
        }*/

        // can reuse binary tree division here

        group_bt(pL,  npL);
        group_bt(pLU, npLU);

      } else if (group_type == OCTREE) {
        int j;
        for(int i = n-1; i > 0; --i) {
          j = random->uniform()*i;
          if(j < 0 || j >= n) error->one(FLERR,"bad index");
          std::swap(plist[i], plist[j]);
        }

        group_ot(0,n);

      }

      // recreate particle list after reduction

      ip = cinfo[icell].first;
      n = 0;
      while (ip >= 0) {
        ipart = &particles[ip];
        isw = sweights[ip];
        if(isw > 0) plist[n++] = ip;
        ip = next[ip];
      }

      // if no particles reduced, increase group size
      
      if (gbuf > n) error->one(FLERR,"too big");

      if (n == nold) gbuf += 2;

    }
  }// loop for cells
  return;
}

/* ----------------------------------------------------------------------
   Recursivley divides particles using the binary tree strategy
------------------------------------------------------------------------- */
void Collide::group_bt(int *plist_leaf, int np)
{
  Particle::OnePart *ipart;
  Particle::OnePart *particles = particle->particles;
  Particle::Species *species = particle->species;

  // ignore groups which have too few particles

  if (np <= Ngmin) return;

  // compute stress tensor since it's needed for
  // .. further dividing and reduction

  double *sweights = particle->edvec[particle->ewhich[index_sweight]];

  double gsum, msum, mV[3], mVV[3][3], mVVV[3][3];
  gsum = msum = 0.0;
  for (int i = 0; i < 3; i++) {
    mV[i] = 0.0;
    for (int j = 0; j < 3; j++) {
      mVV[i][j] = 0.0;
      mVVV[i][j] = 0.0;
    }
  }

  // find maximum particle weight

  int ispecies;
	double mass, psw, pmsw, vp[3];
  double Erot;
  for (int p = 0; p < np; p++) {
    ipart = &particles[plist_leaf[p]];
    ispecies = ipart->ispecies;
    mass = species[ispecies].mass;

    psw = sweights[plist_leaf[p]];
    pmsw = psw * mass;
    memcpy(vp, ipart->v, 3*sizeof(double));
   	gsum += psw;
    msum += pmsw;
    Erot += psw*ipart->erot;
    for (int i = 0; i < 3; i++) {
      mV[i] += (pmsw*vp[i]);
      for (int j = 0; j < 3; j++) {
        mVV[i][j] += (pmsw*vp[i]*vp[j]);
        mVVV[i][j] += (pmsw*vp[i]*vp[j]*vp[j]);
      }
    }
  }

  // mean velocity

	double V[3];
  for (int i = 0; i < 3; i++) V[i] = mV[i]/msum;

  // stress tensor

  double pij[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      pij[i][j] = mVV[i][j] - mV[i]*mV[j]/msum;

  // if group is small enough, merge the particles

  if (np <= Ngmax+gbuf) {

    // remove small stress tensor components

    //for(int i = 0; i < 3; i++)
    //  for(int j = 0; j < 3; j++)
    //    if(fabs(pij[i][j]/mass) < SMALLISH) pij[i][j] = 0.0;

    // temperature
    double T = (pij[0][0] + pij[1][1] + pij[2][2])/
      (3.0 * gsum * update->boltz);

    // heat flux
    double Vsq = V[0]*V[0] + V[1]*V[1] + V[2]*V[2];
    double h,h1,h2,q[3];
    int i1,i2;
    for (int i = 0; i < 3; i++) {
      if (i == 0) {
        i1 = 1;
        i2 = 2;
      } else if (i == 1) {
        i1 = 2;
        i2 = 0;
      } else {
        i1 = 0;
        i2 = 1;
      }

      h  = mVVV[i][i] - 3.0*mV[i]*mVV[i][i]/msum +
           2.0*mV[i]*mV[i]*mV[i]/msum/msum;
      h1 = mVVV[i][i1] - 2.0*mVV[i][i1]*mV[i1]/msum -
           mV[i]*mVV[i1][i1]/msum + 2.0*mV[i]*mV[i1]*mV[i1]/msum/msum;
      h2 = mVVV[i][i2] - 2.0*mVV[i][i2]*mV[i2]/msum -
           mV[i]*mVV[i2][i2]/msum + 2.0*mV[i]*mV[i2]*mV[i2]/msum/msum;
      q[i] = (h + h1 + h2) * 0.5;
    }

    // remove small heat fluxes
    //for(int i = 0; i < 3; i++)
    //  if(fabs(q[i]/mass) < SMALLISH) q[i] = 0.0;

    // scale values to be consistent with definitions in
    // .. stochastic numerics book

    T *= update->boltz/mass;
    for (int i = 0; i < 3; i++) {
      q[i] /= mass;
      for(int j = 0; j < 3; j++) pij[i][j] /= mass;
    }

    // reduce based on type
    if (reduction_type == ENERGY) {
      reduce(plist_leaf, np, gsum, V, T, Erot);
    } else if (reduction_type == HEAT) {
      reduce(plist_leaf, np, gsum, V, T, Erot, q);
    } else if (reduction_type == STRESS) {
      reduce(plist_leaf, np, gsum, V, T, Erot, q, pij);
    }

  // group still too large so divide further

  } else {

    // Compute covariance matrix

    double Rij[3][3];
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        Rij[i][j] = pij[i][j]/gsum;

    // Find eigenpairs

    double eval[3], evec[3][3];
    int ierror = MathEigen::jacobi3(Rij,eval,evec);

    // Find largest eigenpair

    double maxeval;
    double maxevec[3]; // normal of splitting plane

    maxeval = 0;
    for (int i = 0; i < 3; i++) {
      if (std::abs(eval[i]) > maxeval) {
        maxeval = std::abs(eval[i]);
        for (int j = 0; j < 3; j++) {
          maxevec[j] = evec[j][i];  
        }
      }
    }

    // Separate based on particle velocity

    double center = V[0]*maxevec[0] + V[1]*maxevec[1] + V[2]*maxevec[2];
    int pid, pidL[np], pidR[np];
    int npL, npR;
    npL = npR = 0;
    for (int i = 0; i < np; i++) {
      pid = plist_leaf[i];
      ipart = &particles[pid];
      if (MathExtra::dot3(ipart->v,maxevec) < center)
        pidL[npL++] = pid;
      else
        pidR[npR++] = pid;
    }

    if(npL > Ngmin) group_bt(pidL,npL);
    if(npR > Ngmin) group_bt(pidR,npR);
  }

  return;
}

/* ----------------------------------------------------------------------
   Recursivley divides particles using the octree strategy
------------------------------------------------------------------------- */
void Collide::group_ot(int pfirst, int plast)
{
  Particle::OnePart *ipart;
  Particle::OnePart *particles = particle->particles;
  Particle::Species *species = particle->species;
  int np = plast-pfirst;

  // ignore groups which have too few particles

  if (np <= Ngmin)
    return;

  // compute stress tensor since it's needed for
  // .. further dividing and reduction

  double *sweights = particle->edvec[particle->ewhich[index_sweight]];

  double gsum, msum, mV[3], mVV[3][3], mVVV[3][3];
  gsum = msum = 0.0;
  for (int i = 0; i < 3; i++) {
    mV[i] = 0.0;
    for (int j = 0; j < 3; j++) {
      mVV[i][j] = 0.0;
      mVVV[i][j] = 0.0;
    }
  }

  // find maximum particle weight

  int ispecies;
	double mass, psw, pmsw, vp[3];
  double Erot = 0.0;
  for (int p = pfirst; p < plast; p++) {
    ipart = &particles[plist[p]];
    ispecies = ipart->ispecies;
    mass = species[ispecies].mass;

    psw = sweights[plist[p]];
    pmsw = psw * mass;
    memcpy(vp, ipart->v, 3*sizeof(double));
   	gsum += psw;
    msum += pmsw;
    Erot += psw*ipart->erot;
    for (int i = 0; i < 3; i++) {
      mV[i] += (pmsw*vp[i]);
      for (int j = 0; j < 3; j++) {
        mVV[i][j] += (pmsw*vp[i]*vp[j]);
        mVVV[i][j] += (pmsw*vp[i]*vp[j]*vp[j]);
      }
    }
  }

  // mean velocity

	double V[3];
  for (int i = 0; i < 3; i++) V[i] = mV[i]/msum;

  // stress tensor

  double pij[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      pij[i][j] = mVV[i][j] - mV[i]*mV[j]/msum;

  // if group is small enough, merge the particles

  if (np <= Ngmax+gbuf) {

    // remove small stress tensor components

    //for(int i = 0; i < 3; i++)
    //  for(int j = 0; j < 3; j++)
    //    if(fabs(pij[i][j]/mass) < SMALLISH) pij[i][j] = 0.0;

    // temperature

    double T = (pij[0][0] + pij[1][1] + pij[2][2])/
      (3.0 * gsum * update->boltz);

    // heat flux

    double Vsq = V[0]*V[0] + V[1]*V[1] + V[2]*V[2];
    double h,h1,h2,q[3];
    int i1,i2;
    for (int i = 0; i < 3; i++) {
      if (i == 0) {
        i1 = 1;
        i2 = 2;
      } else if(i == 1) {
        i1 = 2;
        i2 = 0;
      } else {
        i1 = 0;
        i2 = 1;
      }

      h  = mVVV[i][i] - 3.0*mV[i]*mVV[i][i]/msum +
           2.0*mV[i]*mV[i]*mV[i]/msum/msum;
      h1 = mVVV[i][i1] - 2.0*mVV[i][i1]*mV[i1]/msum -
           mV[i]*mVV[i1][i1]/msum + 2.0*mV[i]*mV[i1]*mV[i1]/msum/msum;
      h2 = mVVV[i][i2] - 2.0*mVV[i][i2]*mV[i2]/msum -
           mV[i]*mVV[i2][i2]/msum + 2.0*mV[i]*mV[i2]*mV[i2]/msum/msum;
      q[i] = (h + h1 + h2) * 0.5;
    }

    // remove small heat fluxes

    //for(int i = 0; i < 3; i++)
    //  if(fabs(q[i]/mass) < SMALLISH) q[i] = 0.0;

    // scale values to be consistent with definitions in
    // .. stochastic numerics book

    T *= update->boltz/mass;
    for (int i = 0; i < 3; i++) {
      q[i] /= mass;
      for (int j = 0; j < 3; j++) pij[i][j] /= mass;
    }

    // reduce based on type

    if (reduction_type == ENERGY) {
      //reduce(plist_leaf, np, gsum, V, T, Erot);
    } else if (reduction_type == HEAT) {
      //reduce(plist_leaf, np, gsum, V, T, Erot, q);
    } else if (reduction_type == STRESS) {
      //reduce(plist_leaf, np, gsum, V, T, Erot, q, pij);
    }

  // group still too large so divide further

  } else {

    // sort particles into octants

    int temp[8][np];
    int ip, iquad, nquad, nq[8];
    for (int i = 0; i < 8; i++) nq[i] = 0;

    for (int i = pfirst; i < plast; i++) {
      ip = plist[i];
      ipart = &particles[ip];
      memcpy(vp, ipart->v, 3*sizeof(double));

      iquad = 0;
      if (vp[0] > V[0]) iquad += 1;
      if (vp[1] > V[1]) iquad += 2;
      if (vp[2] > V[2]) iquad += 4;

      nquad = nq[iquad];
      temp[iquad][nquad] = ip;
      nq[iquad] = nquad + 1;
    }

    // rebuild particle list

    ip = pfirst;
    for (int i = 0; i < 8; i++)
      for (int j = 0; j < nq[i]; j++)
        plist[ip++] = temp[i][j];

    // start next iteration

    int start, end;
    start = pfirst;
    end = pfirst + nq[0];
    group_ot(start,end);
    for (int i = 0; i < 7; i++) {
      start = end;
      end += nq[i+1];
      group_ot(start,end);
    }
  }

  return;
}

/* ----------------------------------------------------------------------
   Merge particles using energy scheme
------------------------------------------------------------------------- */
void Collide::reduce(int *pleaf, int np,
                     double rho, double *V, double T, double Erot)
{

  // reduced particles chosen as first two in group

  Particle::OnePart *particles = particle->particles;
  Particle::OnePart *ipart, *jpart;

  int ip, jp;
  ip = np * random->uniform();
  jp = np * random->uniform();
  while (ip == jp) jp = np * random->uniform();

  ipart = &particles[pleaf[ip]];
  jpart = &particles[pleaf[jp]];

  double *sweights = particle->edvec[particle->ewhich[index_sweight]];

  // find direction of velocity wrt CoM frame

  double theta = 2.0 * 3.14159 * random->uniform();
  double phi = acos(1.0 - 2.0 * random->uniform());
  double uvec[3];
  uvec[0] = sin(phi) * cos(theta);
  uvec[1] = sin(phi) * sin(theta);
  uvec[2] = cos(phi);

  // set reduced particle velocities

  double sqT = sqrt(3.0*T);
  for (int d = 0; d < 3; d++) {
    ipart->v[d] = V[d] + sqT*uvec[d];
    jpart->v[d] = V[d] - sqT*uvec[d];
  }

  // set reduced particle rotational energies

  ipart->erot = Erot/(rho*0.5)*0.5;
  jpart->erot = Erot/(rho*0.5)*0.5;

  // set reduced particle weights

  sweights[pleaf[ip]] = rho*0.5;
  sweights[pleaf[jp]] = rho*0.5;

  // delete other particles

  for (int i = 0; i < np; i++) {
    if (i == ip || i == jp) continue;
    if (ndelete == maxdelete) {
      maxdelete += DELTADELETE;
      memory->grow(dellist,maxdelete,"collide:dellist");
    }
    sweights[pleaf[i]] = -1.0;
    dellist[ndelete++] = pleaf[i];
  }

  return;
}

/* ----------------------------------------------------------------------
   Merge particles using heat flux scheme
------------------------------------------------------------------------- */
void Collide::reduce(int *pleaf, int np,
                     double rho, double *V, double T, double Erot, double *q)
{

  // reduced particles chosen as first two in group

  Particle::OnePart *particles = particle->particles;
  Particle::OnePart *ipart, *jpart;

  int ip, jp;
  ip = np * random->uniform();
  jp = np * random->uniform();
  while (ip == jp) jp = np * random->uniform();

  ipart = &particles[pleaf[ip]];
  jpart = &particles[pleaf[jp]];

  double *sweights = particle->edvec[particle->ewhich[index_sweight]];

  // precompute

  double sqT = sqrt(3.0*T);
  double qmag = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);
  double qge = qmag / (rho * pow(sqT,3.0));
  double itheta = qge + sqrt(1.0 + qge*qge);
  double alpha = sqT*itheta;
  double beta = sqT/itheta;

  // find direction of velocity wrt CoM frame

  double uvec[3];
  if (qmag < SMALL) {
    for (int d = 0; d < 3; d++) {
      double A = sqrt(-log(random->uniform()));
      double phi = 6.283185308 * random->uniform();
      if (random->uniform() < 0.5) uvec[d] = A * cos(phi);
      else uvec[d] = A * sin(phi);
    }
  } else 
    for (int d = 0; d < 3; d++) uvec[d] = q[d]/qmag;

  // set reduced particle velocities

  for (int d = 0; d < 3; d++) {
    ipart->v[d] = V[d] + alpha*uvec[d];
    jpart->v[d] = V[d] - beta*uvec[d];
  }

  // set reduced particle weights

  double isw = rho/(1.0+itheta*itheta);
  double jsw = rho - isw;

  // set reduced particle rotational energies

  ipart->erot = Erot/isw*0.5;
  jpart->erot = Erot/jsw*0.5;

  sweights[pleaf[ip]] = isw;
  sweights[pleaf[jp]] = jsw;

  // delete other particles
  for (int i = 0; i < np; i++) {
    if (i == ip || i == jp) continue;
    if (ndelete == maxdelete) {
      maxdelete += DELTADELETE;
      memory->grow(dellist,maxdelete,"collide:dellist");
    }
    sweights[pleaf[i]] = -1.0;
    dellist[ndelete++] = pleaf[i];
  }

  return;
}

/* ----------------------------------------------------------------------
   Merge particles using stress scheme
------------------------------------------------------------------------- */
void Collide::reduce(int *pleaf, int np,
                     double rho, double *V, double T, double Erot,
                     double *q, double pij[3][3])
{

  // find eigenpairs of stress tensor

  double eval[3], evec[3][3];
  int ierror = MathEigen::jacobi3(pij,eval,evec);

  // find number of non-zero eigenvalues

  int nK = 0;
  for (int i = 0; i < 3; i++) {
    if (fabs(eval[i]) >= SMALL && eval[i] > 0) {
      eval[nK] = eval[i];
      for (int d = 0; d < 3; d++) evec[nK][d] = evec[i][d];
      nK++;
    }
  }

  // iterate through nK pairs

  Particle::OnePart *particles = particle->particles;
  Particle::OnePart *ipart, *jpart;

  double *sweights = particle->edvec[particle->ewhich[index_sweight]];

  double qli, itheta;
  double isw, jsw;
  double uvec[3];

  for (int iK = 0; iK < nK; iK++) {

    // reduced particles chosen as first two

    ipart = &particles[pleaf[2*iK]];
    jpart = &particles[pleaf[2*iK+1]];

    qli = evec[0][iK]*q[0] + evec[1][iK]*q[1] + evec[2][iK]*q[2];
    if (qli < 0)
      for (int d = 0; d < 3; d++) evec[d][iK] *= -1.0;
    qli = fabs(qli);

    itheta = sqrt(rho) * qli / (sqrt(nK) * pow(eval[iK],1.5))
      + sqrt(1.0 + (rho*qli*qli)/(nK*pow(eval[iK],3.0)));

    // set reduced particle velocities

    for (int d = 0; d < 3; d++) {
      ipart->v[d] = V[d] + itheta*sqrt(nK*eval[iK]/rho)*evec[d][iK];
      jpart->v[d] = V[d] - 1.0/itheta*sqrt(nK*eval[iK]/rho)*evec[d][iK];
    }

    // set reduced particle weights

    isw = rho/(nK*(1.0+itheta*itheta));
    jsw = rho/nK - isw;

    // set reduced particle rotational energies

    ipart->erot = Erot/isw*0.5/nK;
    jpart->erot = Erot/jsw*0.5/nK;

    sweights[pleaf[2*iK]] = isw;
    sweights[pleaf[2*iK+1]] = jsw;
  } // end nK
  

  // delete other particles
  for (int i = 0; i < np; i++) {
    if (i < 2*nK) continue;
    if (ndelete == maxdelete) {
      maxdelete += DELTADELETE;
      memory->grow(dellist,maxdelete,"collide:dellist");
    }
    sweights[pleaf[i]] = -1.0;
    dellist[ndelete++] = pleaf[i];
  }

  return;
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
   copy per-cell collision info from Icell to Jcell
   called whenever a grid cell is removed from this processor's list
   caller checks that Icell != Jcell
------------------------------------------------------------------------- */

void Collide::copy_grid_one(int icell, int jcell)
{
  int nbytes = ngroups*ngroups*sizeof(double);

  memcpy(&vremax[jcell][0][0],&vremax[icell][0][0],nbytes);
  if (remainflag)
    memcpy(&remain[jcell][0][0],&remain[icell][0][0],nbytes);
}

/* ----------------------------------------------------------------------
   reset final grid cell count after grid cell removals
------------------------------------------------------------------------- */

void Collide::reset_grid_count(int nlocal)
{
  nglocal = nlocal;
}

/* ----------------------------------------------------------------------
   add a grid cell
   called when a grid cell is added to this processor's list
   initialize values to 0.0
------------------------------------------------------------------------- */

void Collide::add_grid_one()
{
  grow_percell(1);

  for (int igroup = 0; igroup < ngroups; igroup++)
    for (int jgroup = 0; jgroup < ngroups; jgroup++) {
      vremax[nglocal][igroup][jgroup] = vremax_initial[igroup][jgroup];
      if (remainflag) remain[nglocal][igroup][jgroup] = 0.0;
    }

  nglocal++;
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
  if (nglocal+n < nglocalmax || !ngroups) return;
  while (nglocal+n >= nglocalmax) nglocalmax += DELTAGRID;
  memory->grow(vremax,nglocalmax,ngroups,ngroups,"collide:vremax");
  if (remainflag)
    memory->grow(remain,nglocalmax,ngroups,ngroups,"collide:remain");
}

/* ----------------------------------------------------------------------
   for particle I, find collision partner J via near neighbor algorithm
   always returns a J neighbor, even if not that near
   near neighbor algorithm:
     check up to nearlimit particles, starting with random particle
     as soon as find one within distance moved by particle I, return it
     else return the closest one found
     also exclude an I,J pair if both most recently collided with each other
   this version is for single group collisions
------------------------------------------------------------------------- */

int Collide::find_nn(int i, int np)
{
  int jneigh;
  double dx,dy,dz,rsq;
  double *xj;

  // if np = 2, just return J = non-I particle
  // np is never < 2

  if (np == 2) return (i+1) % 2;

  Particle::OnePart *ipart,*jpart;
  Particle::OnePart *particles = particle->particles;
  double dt = update->dt;

  // thresh = distance particle I moves in this timestep

  ipart = &particles[plist[i]];
  double *vi = ipart->v;
  double *xi = ipart->x;
  double threshsq =  dt*dt * (vi[0]*vi[0]+vi[1]*vi[1]+vi[2]*vi[2]);
  double minrsq = BIG;

  // nlimit = max # of J candidates to consider

  int nlimit = MIN(nearlimit,np-1);
  int count = 0;

  // pick a random starting J
  // jneigh = collision partner when exit loop
  //   set to initial J as default in case no Nlimit J meets criteria

  int j = np * random->uniform();
  while (i == j) j = np * random->uniform();
  jneigh = j;

  while (count < nlimit) {
    count++;

    // skip this J if I,J last collided with each other

    if (nn_last_partner[i] == j+1 && nn_last_partner[j] == i+1) {
      j++;
      if (j == np) j = 0;
      continue;
    }

    // rsq = squared distance between particles I and J
    // if rsq = 0.0, skip this J
    //   could be I = J, or a cloned J at same position as I
    // if rsq <= threshsq, this J is collision partner
    // if rsq = smallest yet seen, this J is tentative collision partner

    jpart = &particles[plist[j]];
    xj = jpart->x;
    dx = xi[0] - xj[0];
    dy = xi[1] - xj[1];
    dz = xi[2] - xj[2];
    rsq = dx*dx + dy*dy + dz*dz;

    if (rsq > 0.0) {
      if (rsq <= threshsq) {
        jneigh = j;
        break;
      }
      if (rsq < minrsq) {
        minrsq = rsq;
        jneigh = j;
      }
    }

    j++;
    if (j == np) j = 0;
  }

  return jneigh;
}

/* ----------------------------------------------------------------------
   for particle I, find collision partner J via near neighbor algorithm
   always returns a J neighbor, even if not that near
   same near neighbor algorithm as in find_nn()
     looking for J particles in jlist of length Np = ngroup[jgroup]
     ilist = jlist when igroup = jgroup
   this version is for multi group collisions
------------------------------------------------------------------------- */

int Collide::find_nn_group(int i, int *ilist, int np, int *jlist, int *plist,
                           int *nn_igroup, int *nn_jgroup)
{
  int jneigh;
  double dx,dy,dz,rsq;
  double *xj;

  // if ilist = jlist and np = 2, just return J = non-I particle
  // np is never < 2 for ilist = jlist
  // np is never < 1 for ilist != jlist

  if (ilist == jlist && np == 2) return (i+1) % 2;

  Particle::OnePart *ipart,*jpart;
  Particle::OnePart *particles = particle->particles;
  double dt = update->dt;

  // thresh = distance particle I moves in this timestep

  ipart = &particles[plist[ilist[i]]];
  double *vi = ipart->v;
  double *xi = ipart->x;
  double threshsq =  dt*dt * (vi[0]*vi[0]+vi[1]*vi[1]+vi[2]*vi[2]);
  double minrsq = BIG;

  // nlimit = max # of J candidates to consider

  int nlimit = MIN(nearlimit,np-1);
  int count = 0;

  // pick a random starting J
  // jneigh = collision partner when exit loop
  //   set to initial J as default in case no Nlimit J meets criteria

  int j = np * random->uniform();
  if (ilist == jlist)
    while (i == j) j = np * random->uniform();
  jneigh = j;

  while (count < nlimit) {
    count++;

    // skip this J if I,J last collided with each other

    if (nn_igroup[i] == j+1 && nn_jgroup[j] == i+1) {
      j++;
      if (j == np) j = 0;
      continue;
    }

    // rsq = squared distance between particles I and J
    // if rsq = 0.0, skip this J
    //   could be I = J, or a cloned J at same position as I
    // if rsq <= threshsq, this J is collision partner
    // if rsq = smallest yet seen, this J is tentative collision partner

    jpart = &particles[plist[jlist[j]]];
    xj = jpart->x;
    dx = xi[0] - xj[0];
    dy = xi[1] - xj[1];
    dz = xi[2] - xj[2];
    rsq = dx*dx + dy*dy + dz*dz;

    if (rsq > 0.0) {
      if (rsq <= threshsq) {
        jneigh = j;
        break;
      }
      if (rsq < minrsq) {
        minrsq = rsq;
        jneigh = j;
      }
    }

    j++;
    if (j == np) j = 0;
  }

  return jneigh;
}

/* ----------------------------------------------------------------------
   reallocate a nn_last_partner vector to allow for N values
   increase size by multiples of 2x
------------------------------------------------------------------------- */

void Collide::realloc_nn(int n, int *&vec)
{
  while (n > max_nn) max_nn *= 2;
  memory->destroy(vec);
  memory->create(vec,max_nn,"collide:nn_last_partner");
}

/* ----------------------------------------------------------------------
   set nn_last_partner[N] = 0 for newly created particle
   grow the vector if necessary
------------------------------------------------------------------------- */

void Collide::set_nn(int n)
{
  if (n == max_nn) {
    max_nn *= 2;
    memory->grow(nn_last_partner,max_nn,"collide:nn_last_partner");
  }
  nn_last_partner[n] = 0;
}

/* ----------------------------------------------------------------------
   grow the group last partner vectors if necessary
------------------------------------------------------------------------- */

void Collide::set_nn_group(int n)
{
  if (n == max_nn) {
    max_nn *= 2;
    memory->grow(nn_last_partner_igroup,max_nn,"collide:nn_last_partner");
    memory->grow(nn_last_partner_jgroup,max_nn,"collide:nn_last_partner");
  }
}

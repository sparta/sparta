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

#include "string.h"
#include "compute_surf_collision_tally.h"
#include "particle.h"
#include "mixture.h"
#include "surf.h"
#include "update.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{IDSURF,ID,TYPE,TIME,XC,YC,ZC,VXPRE,VYPRE,VZPRE,VXPOST,VYPOST,VZPOST};
enum{DOUBLE,INT,BIGINT,UINT,BIGUINT,STRING};    // same as Dump

#define DELTA 4096

/* ---------------------------------------------------------------------- */

ComputeSurfCollisionTally::ComputeSurfCollisionTally(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal compute surf/collision/tally command");

  int igroup = surf->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Compute surf/collision/tally group ID does not exist");
  groupbit = surf->bitmask[igroup];

  imix = particle->find_mixture(arg[3]);
  if (imix < 0) error->all(FLERR,"Compute surf/collision/tally mixture ID does not exist");

  nvalue = narg - 4;
  which = new int[nvalue];

  // process input values

  nvalue = 0;
  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"id/surf") == 0) which[nvalue++] = IDSURF;
    else if (strcmp(arg[iarg],"id") == 0) which[nvalue++] = ID;
    else if (strcmp(arg[iarg],"type") == 0) which[nvalue++] = TYPE;
    else if (strcmp(arg[iarg],"time") == 0) which[nvalue++] = TIME;
    else if (strcmp(arg[iarg],"xc") == 0) which[nvalue++] = XC;
    else if (strcmp(arg[iarg],"yc") == 0) which[nvalue++] = YC;
    else if (strcmp(arg[iarg],"zc") == 0) which[nvalue++] = ZC;
    else if (strcmp(arg[iarg],"vx/pre") == 0) which[nvalue++] = VXPRE;
    else if (strcmp(arg[iarg],"vy/pre") == 0) which[nvalue++] = VYPRE;
    else if (strcmp(arg[iarg],"vz/pre") == 0) which[nvalue++] = VZPRE;
    else if (strcmp(arg[iarg],"vx/post") == 0) which[nvalue++] = VXPOST;
    else if (strcmp(arg[iarg],"vy/post") == 0) which[nvalue++] = VYPOST;
    else if (strcmp(arg[iarg],"vz/post") == 0) which[nvalue++] = VZPOST;
    else error->all(FLERR,"Invalid value for compute surf/collision/tally");
    iarg++;
  }

  // setup

  per_tally_flag = 1;
  size_per_tally_cols = nvalue;

  surf_tally_flag = 1;         // triggers Update to invoke surf_tally() for each collision
  timeflag = 1;                // tells Update which timesteps to invoke surf_tally()

  ntally = maxtally = 0;
  array_tally = NULL;
  
  dim = domain->dimension;
}

/* ---------------------------------------------------------------------- */

ComputeSurfCollisionTally::~ComputeSurfCollisionTally()
{
  if (copy || copymode) return;

  delete [] which;

  memory->destroy(array_tally);
}

/* ---------------------------------------------------------------------- */

void ComputeSurfCollisionTally::init()
{
  if (!surf->exist)
    error->all(FLERR,"Cannot use compute surf/collision/tally when surfs do not exist");
}

/* ----------------------------------------------------------------------
   no operations here, since compute results are stored in tally array
   just used by callers to indicate compute was used
   enables prediction of next step when update needs to tally
------------------------------------------------------------------------- */

void ComputeSurfCollisionTally::compute_per_tally()
{
  invoked_per_tally = update->ntimestep;
}

/* ----------------------------------------------------------------------
   called by Update before timesteps if will invoke surf_tally()
---------------------------------------------------------------------- */

void ComputeSurfCollisionTally::clear()
{
  lines = surf->lines;
  tris = surf->tris;

  ntally = 0;
}

/* ----------------------------------------------------------------------
   tally values for a single particle in icell
     colliding with surface element isurf, performing reaction (1 to N)
   iorig = particle before collision
   ip,jp = particles after collision
   ip = NULL means no particles after collision
   jp = NULL means one particle after collision
   jp != NULL means two particles after collision
------------------------------------------------------------------------- */

void ComputeSurfCollisionTally::surf_tally(double dtremain, int isurf,
                                           int icell, int reaction,
                                           Particle::OnePart *iorig,
                                           Particle::OnePart *ip, Particle::OnePart *jp)
{
  // skip if a reaction
  // this compute only tallies collisions that do not induce a reaction
  // reactions can be tallied by compute surf/reaction/tally command

  if (ip == NULL || jp) return;
  if (iorig->ispecies != ip->ispecies) return;
  
  // skip if isurf not in surface group

  if (dim == 2) {
    if (!(lines[isurf].mask & groupbit)) return;
  } else {
    if (!(tris[isurf].mask & groupbit)) return;
  }

  // skip if particle species not in mixture group

  int origspecies = iorig->ispecies;
  int igroup = particle->mixture[imix]->species2group[origspecies];
  if (igroup < 0) return;

  // grow tally array if necessary
  
  if (ntally == maxtally) grow_tally();

  // tally all values associated with group into array
  // particle iorig,ip have same collision point but before/after velocities
  
  double *vec = array_tally[ntally++];
  
  for (int m = 0; m < nvalue; m++) {
    switch (which[m]) {
    case IDSURF:
      if (dim == 2) vec[m] = ubuf(lines[isurf].id).d;
      else vec[m] = ubuf(tris[isurf].id).d;
      break;
    case ID:
      vec[m] = ubuf(ip->id).d;
      break;
    case TYPE:
      vec[m] = ubuf(ip->ispecies+1).d;
      break;
    case XC:
      vec[m] = iorig->x[0];
      break;
    case YC: 
      vec[m] = iorig->x[1];
      break;
    case ZC: 
      vec[m] = iorig->x[2];
      break;
    case TIME: 
      vec[m] = update->dt - dtremain;
      break;
    case VXPRE:
      vec[m] = iorig->v[0];
      break;
    case VYPRE: 
      vec[m] = iorig->v[1];
      break;
    case VZPRE: 
      vec[m] = iorig->v[2];
      break;
    case VXPOST:
      vec[m] = ip->v[0];
      break;
    case VYPOST: 
      vec[m] = ip->v[1];
      break;
    case VZPOST: 
      vec[m] = ip->v[2];
      break;
    }
  }
}

/* ----------------------------------------------------------------------
   return # of tallies
------------------------------------------------------------------------- */

int ComputeSurfCollisionTally::tallyinfo(surfint *&dummy)
{
  return ntally;
}

/* ----------------------------------------------------------------------
   return datatype of tally quantity
   icol = 0 for vector
   icol = 1 to N for array column
   datatype = INT,DOUBLE,BIGINT
------------------------------------------------------------------------- */

int ComputeSurfCollisionTally::datatype(int icol)
{
  if (which[icol-1] == IDSURF) {
    if (sizeof(surfint) == sizeof(smallint)) return INT;
    if (sizeof(surfint) == sizeof(bigint)) return BIGINT;
  }
  if (which[icol-1] == ID) return INT;
  if (which[icol-1] == TYPE) return INT;

  return DOUBLE;
}

/* ---------------------------------------------------------------------- */

void ComputeSurfCollisionTally::grow_tally()
{
  maxtally += DELTA;
  memory->grow(array_tally,maxtally,nvalue,"surf/collision/tally:array_tally");
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

bigint ComputeSurfCollisionTally::memory_usage()
{
  bigint bytes = 0;
  bytes += nvalue*maxtally * sizeof(double);    // array_tally
  return bytes;
}

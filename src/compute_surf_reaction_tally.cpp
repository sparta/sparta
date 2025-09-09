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
#include "compute_surf_reaction_tally.h"
#include "particle.h"
#include "mixture.h"
#include "surf.h"
#include "update.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{REACTION,IDSURF,IDPRE,ID1POST,ID2POST,TYPEPRE,TYPE1POST,TYPE2POST,TIME,XC,YC,ZC,
  VXPRE,VYPRE,VZPRE,VX1POST,VY1POST,VZ1POST,VX2POST,VY2POST,VZ2POST};
enum{INT,DOUBLE,BIGINT,UINT,BIGUINT,STRING};    // same as Dump

#define DELTA 4096

/* ---------------------------------------------------------------------- */

ComputeSurfReactionTally::ComputeSurfReactionTally(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal compute surf/reaction/tally command");

  int igroup = surf->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Compute surf/reaction/tally group ID does not exist");
  groupbit = surf->bitmask[igroup];

  imix = particle->find_mixture(arg[3]);
  if (imix < 0) error->all(FLERR,"Compute surf/reaction/tally mixture ID does not exist");

  nvalue = narg - 4;
  which = new int[nvalue];

  // process input values

  nvalue = 0;
  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"reaction") == 0) which[nvalue++] = REACTION;
    else if (strcmp(arg[iarg],"id/surf") == 0) which[nvalue++] = IDSURF;
    else if (strcmp(arg[iarg],"id/pre") == 0) which[nvalue++] = IDPRE;
    else if (strcmp(arg[iarg],"id1/post") == 0) which[nvalue++] = ID1POST;
    else if (strcmp(arg[iarg],"id2/post") == 0) which[nvalue++] = ID2POST;
    else if (strcmp(arg[iarg],"type/pre") == 0) which[nvalue++] = TYPEPRE;
    else if (strcmp(arg[iarg],"type1/post") == 0) which[nvalue++] = TYPE1POST;
    else if (strcmp(arg[iarg],"type2/post") == 0) which[nvalue++] = TYPE2POST;
    else if (strcmp(arg[iarg],"time") == 0) which[nvalue++] = TIME;
    else if (strcmp(arg[iarg],"xc") == 0) which[nvalue++] = XC;
    else if (strcmp(arg[iarg],"yc") == 0) which[nvalue++] = YC;
    else if (strcmp(arg[iarg],"zc") == 0) which[nvalue++] = ZC;
    else if (strcmp(arg[iarg],"vx/pre") == 0) which[nvalue++] = VXPRE;
    else if (strcmp(arg[iarg],"vy/pre") == 0) which[nvalue++] = VYPRE;
    else if (strcmp(arg[iarg],"vz/pre") == 0) which[nvalue++] = VZPRE;
    else if (strcmp(arg[iarg],"vx1/post") == 0) which[nvalue++] = VX1POST;
    else if (strcmp(arg[iarg],"vy1/post") == 0) which[nvalue++] = VY1POST;
    else if (strcmp(arg[iarg],"vz1/post") == 0) which[nvalue++] = VZ1POST;
    else if (strcmp(arg[iarg],"vx2/post") == 0) which[nvalue++] = VX2POST;
    else if (strcmp(arg[iarg],"vy2/post") == 0) which[nvalue++] = VY2POST;
    else if (strcmp(arg[iarg],"vz2/post") == 0) which[nvalue++] = VZ2POST;
    else error->all(FLERR,"Invalid value for compute surf/reaction/tally");
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

ComputeSurfReactionTally::~ComputeSurfReactionTally()
{
  if (copy || copymode) return;

  delete [] which;

  memory->destroy(array_tally);
}

/* ---------------------------------------------------------------------- */

void ComputeSurfReactionTally::init()
{
  if (!surf->exist)
    error->all(FLERR,"Cannot use compute surf/reaction/tally when surfs do not exist");
}

/* ----------------------------------------------------------------------
   no operations here, since compute results are stored in tally array
   just used by callers to indicate compute was used
   enables prediction of next step when update needs to tally
------------------------------------------------------------------------- */

void ComputeSurfReactionTally::compute_per_tally()
{
  invoked_per_tally = update->ntimestep;
}

/* ----------------------------------------------------------------------
   called by Update before timesteps if will invoke surf_tally()
---------------------------------------------------------------------- */

void ComputeSurfReactionTally::clear()
{
  lines = surf->lines;
  tris = surf->tris;

  ntally = 0;
}

/* ----------------------------------------------------------------------
   tally values for a single particle in icell colliding with
     surface element isurf
   reaction = 0 for collision only
   reaction = 1 to N for which reaction
   iorig = particle before collision
   ip,jp = particles after collision
   ip = NULL means no particles after collision
   jp = NULL means one particle after collision
   jp != NULL means two particles after collision
------------------------------------------------------------------------- */

void ComputeSurfReactionTally::surf_tally(double dtremain, int isurf,
                                          int icell, int reaction,
                                          Particle::OnePart *iorig,
                                          Particle::OnePart *ip,
                                          Particle::OnePart *jp)
{
  // skip if not a reaction
  // this compute only tallies collisions that induce a reaction
  // simple collisions can be tallied by compute surf/collision/tally command

  if (!reaction) return;

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
    case REACTION:
      vec[m] = ubuf(reaction).d;
      break;

    case IDSURF:
      if (dim == 2) vec[m] = ubuf(lines[isurf].id).d;
      else vec[m] = ubuf(tris[isurf].id).d;
      break;
    case IDPRE:
      vec[m] = ubuf(iorig->id).d;
      break;
    case ID1POST:
      if (ip == NULL) vec[m] = ubuf(0).d;
      else vec[m] = ubuf(ip->id).d;
      break;
    case ID2POST:
      if (jp == NULL) vec[m] = ubuf(0).d;
      else vec[m] = ubuf(jp->id).d;
      break;

    case TYPEPRE:
      vec[m] = ubuf(iorig->ispecies+1).d;
      break;
    case TYPE1POST:
      if (ip == NULL) vec[m] = ubuf(0).d;
      else vec[m] = ubuf(ip->ispecies+1).d;
      break;
    case TYPE2POST:
      if (jp == NULL) vec[m] = ubuf(0).d;
      else vec[m] = ubuf(jp->ispecies+1).d;
      break;

    case TIME:
      vec[m] = update->dt - dtremain;
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

    case VXPRE:
      vec[m] = iorig->v[0];
      break;
    case VYPRE:
      vec[m] = iorig->v[1];
      break;
    case VZPRE:
      vec[m] = iorig->v[2];
      break;
    case VX1POST:
      if (ip == NULL) vec[m] = 0.0;
      else vec[m] = ip->v[0];
      break;
    case VY1POST:
      if (ip == NULL) vec[m] = 0.0;
      else vec[m] = ip->v[1];
      break;
    case VZ1POST:
      if (ip == NULL) vec[m] = 0.0;
      else vec[m] = ip->v[2];
      break;
    case VX2POST:
      if (jp == NULL) vec[m] = 0.0;
      else vec[m] = jp->v[0];
      break;
    case VY2POST:
      if (jp == NULL) vec[m] = 0.0;
      else vec[m] = jp->v[1];
      break;
    case VZ2POST:
      if (jp == NULL) vec[m] = 0.0;
      else vec[m] = jp->v[2];
      break;
    }
  }
}

/* ----------------------------------------------------------------------
   return # of tallies
------------------------------------------------------------------------- */

int ComputeSurfReactionTally::tallyinfo(surfint *&dummy)
{
  return ntally;
}

/* ----------------------------------------------------------------------
   return datatype of tally quantity
   icol = 0 for vector
   icol = 1 to N for array column
   datatype = INT,DOUBLE,BIGINT
------------------------------------------------------------------------- */

int ComputeSurfReactionTally::datatype(int icol)
{
  if (which[icol-1] == REACTION) return INT;
  if (which[icol-1] == IDSURF) {
    if (sizeof(surfint) == sizeof(smallint)) return INT;
    if (sizeof(surfint) == sizeof(bigint)) return BIGINT;
  }
  if (which[icol-1] == IDPRE) return INT;
  if (which[icol-1] == ID2POST || which[icol-1] == ID2POST) return INT;
  if (which[icol-1] == TYPEPRE) return INT;
  if (which[icol-1] == TYPE1POST || which[icol-1] == TYPE2POST) return INT;

  return DOUBLE;
}

/* ---------------------------------------------------------------------- */

void ComputeSurfReactionTally::grow_tally()
{
  maxtally += DELTA;
  memory->grow(array_tally,maxtally,nvalue,"surf/reaction/tally:array_tally");
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

bigint ComputeSurfReactionTally::memory_usage()
{
  bigint bytes = 0;
  bytes += nvalue*maxtally * sizeof(double);    // array_tally
  return bytes;
}

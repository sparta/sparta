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
#include "compute_collide_tally.h"
#include "particle.h"
#include "mixture.h"
#include "surf.h"
#include "update.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{IDPART,IDSURF,XC,YC,ZC,VXOLD,VYOLD,VZOLD,VXNEW,VYNEW,VZNEW};

#define DELTA 4096

/* ---------------------------------------------------------------------- */

ComputeCollideTally::ComputeCollideTally(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal compute collide/tally command");

  int igroup = surf->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Compute collide/tally group ID does not exist");
  groupbit = surf->bitmask[igroup];

  imix = particle->find_mixture(arg[3]);
  if (imix < 0) error->all(FLERR,"Compute collide/tally mixture ID does not exist");

  nvalue = narg - 4;
  which = new int[nvalue];

  // process input values

  nvalue = 0;
  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"id/part") == 0) which[nvalue++] = IDPART;
    else if (strcmp(arg[iarg],"id/surf") == 0) which[nvalue++] = IDSURF;
    else if (strcmp(arg[iarg],"xc") == 0) which[nvalue++] = XC;
    else if (strcmp(arg[iarg],"yc") == 0) which[nvalue++] = YC;
    else if (strcmp(arg[iarg],"zc") == 0) which[nvalue++] = ZC;
    else if (strcmp(arg[iarg],"vxold") == 0) which[nvalue++] = VXOLD;
    else if (strcmp(arg[iarg],"vyold") == 0) which[nvalue++] = VYOLD;
    else if (strcmp(arg[iarg],"vzold") == 0) which[nvalue++] = VZOLD;
    else if (strcmp(arg[iarg],"vxnew") == 0) which[nvalue++] = VXNEW;
    else if (strcmp(arg[iarg],"vynew") == 0) which[nvalue++] = VYNEW;
    else if (strcmp(arg[iarg],"vznew") == 0) which[nvalue++] = VZNEW;
    else error->all(FLERR,"Invalid value for compute collide/tally");
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

ComputeCollideTally::~ComputeCollideTally()
{
  if (copy || copymode) return;

  delete [] which;

  memory->destroy(array_tally);
}

/* ---------------------------------------------------------------------- */

void ComputeCollideTally::init()
{
  if (!surf->exist)
    error->all(FLERR,"Cannot use compute collide/tally when surfs do not exist");
  if (surf->implicit)
    error->all(FLERR,"Cannot use compute collide/tally with implicit surfs");
}

/* ----------------------------------------------------------------------
   no operations here, since compute results are stored in tally array
   just used by callers to indicate compute was used
   enables prediction of next step when update needs to tally
------------------------------------------------------------------------- */

void ComputeCollideTally::compute_per_tally()
{
  invoked_per_tally = update->ntimestep;
}

/* ----------------------------------------------------------------------
   called by Update before timesteps if will invoke surf_tally()
---------------------------------------------------------------------- */

void ComputeCollideTally::clear()
{
  lines = surf->lines;
  tris = surf->tris;

  ntally = 0;
}

/* ----------------------------------------------------------------------
   tally values for a single particle in icell
     colliding with surface element isurf, performing reaction (1 to N)
   iorig = particle ip before collision
   ip,jp = particles after collision
   ip = NULL means no particles after collision
   jp = NULL means one particle after collision
   jp != NULL means two particles after collision
------------------------------------------------------------------------- */

void ComputeCollideTally::surf_tally(int isurf, int icell, int reaction,
                                     Particle::OnePart *iorig,
                                     Particle::OnePart *ip, Particle::OnePart *jp)
{
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
  // particle iorig and ip have same position = collision point
  
  double *vec = array_tally[ntally++];
  
  for (int m = 0; m < nvalue; m++) {
    switch (which[m]) {
    case IDPART:
      vec[m] = ip->id;
      break;
    case IDSURF:
      if (dim == 2) vec[m] = lines[isurf].id;
      else vec[m] = tris[isurf].id;
      break;
    case XC:
      vec[m] = ip->x[0];
      break;
    case YC: 
      vec[m] = ip->x[1];
      break;
    case ZC: 
      vec[m] = ip->x[2];
      break;
    case VXOLD:
      vec[m] = iorig->v[0];
      break;
    case VYOLD: 
      vec[m] = iorig->v[1];
      break;
    case VZOLD: 
      vec[m] = iorig->v[2];
      break;
    case VXNEW:
      vec[m] = ip->v[0];
      break;
    case VYNEW: 
      vec[m] = ip->v[1];
      break;
    case VZNEW: 
      vec[m] = ip->v[2];
      break;
    }
  }
}

/* ----------------------------------------------------------------------
   return # of tallies
------------------------------------------------------------------------- */

int ComputeCollideTally::tallyinfo(surfint *&dummy)
{
  return ntally;
}

/* ---------------------------------------------------------------------- */

void ComputeCollideTally::grow_tally()
{
  maxtally += DELTA;
  memory->grow(array_tally,maxtally,nvalue,"collide/tally:array_tally");
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

bigint ComputeCollideTally::memory_usage()
{
  bigint bytes = 0;
  bytes += nvalue*maxtally * sizeof(double);    // array_tally
  return bytes;
}

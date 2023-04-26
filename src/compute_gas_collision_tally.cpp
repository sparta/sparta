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
#include "compute_gas_collision_tally.h"
#include "particle.h"
#include "mixture.h"
#include "grid.h"
#include "update.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{ID1,ID2,IDCELL,TYPE1,TYPE2,VX1PRE,VY1PRE,VZ1PRE,VX2PRE,VY2PRE,VZ2PRE,
  VX1POST,VY1POST,VZ1POST,VX2POST,VY2POST,VZ2POST};
enum{DOUBLE,INT,BIGINT,UINT,BIGUINT,STRING};    // same as Dump

#define DELTA 4096

/* ---------------------------------------------------------------------- */

ComputeGasCollisionTally::ComputeGasCollisionTally(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal compute gas/collision/tally command");

  int igroup = grid->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Compute gas/collision/tally group ID does not exist");
  groupbit = grid->bitmask[igroup];

  imix = particle->find_mixture(arg[3]);
  if (imix < 0) error->all(FLERR,"Compute gas/collision/tally mixture ID does not exist");

  nvalue = narg - 4;
  which = new int[nvalue];

  // process input values

  nvalue = 0;
  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"id1") == 0) which[nvalue++] = ID1;
    else if (strcmp(arg[iarg],"id2") == 0) which[nvalue++] = ID2;
    else if (strcmp(arg[iarg],"id/cell") == 0) which[nvalue++] = IDCELL;
    else if (strcmp(arg[iarg],"type1") == 0) which[nvalue++] = TYPE1;
    else if (strcmp(arg[iarg],"type2") == 0) which[nvalue++] = TYPE1;
    else if (strcmp(arg[iarg],"vx1/pre") == 0) which[nvalue++] = VX1PRE;
    else if (strcmp(arg[iarg],"vy1/pre") == 0) which[nvalue++] = VY1PRE;
    else if (strcmp(arg[iarg],"vz1/pre") == 0) which[nvalue++] = VZ1PRE;
    else if (strcmp(arg[iarg],"vx2/pre") == 0) which[nvalue++] = VX2PRE;
    else if (strcmp(arg[iarg],"vy2/pre") == 0) which[nvalue++] = VY2PRE;
    else if (strcmp(arg[iarg],"vz2/pre") == 0) which[nvalue++] = VZ2PRE;
    else if (strcmp(arg[iarg],"vx1/post") == 0) which[nvalue++] = VX1POST;
    else if (strcmp(arg[iarg],"vy1/post") == 0) which[nvalue++] = VY1POST;
    else if (strcmp(arg[iarg],"vz1/post") == 0) which[nvalue++] = VZ1POST;
    else if (strcmp(arg[iarg],"vx2/post") == 0) which[nvalue++] = VX2POST;
    else if (strcmp(arg[iarg],"vy2/post") == 0) which[nvalue++] = VY2POST;
    else if (strcmp(arg[iarg],"vz2/post") == 0) which[nvalue++] = VZ2POST;
    else error->all(FLERR,"Invalid value for compute gas/collision/tally");
    iarg++;
  }

  // setup

  per_tally_flag = 1;
  size_per_tally_cols = nvalue;

  gas_tally_flag = 1;         // triggers Collide to invoke gas_tally() for each collision
  timeflag = 1;               // tells Collide which timesteps to invoke gas_tally()

  ntally = maxtally = 0;
  array_tally = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeGasCollisionTally::~ComputeGasCollisionTally()
{
  if (copy || copymode) return;

  delete [] which;

  memory->destroy(array_tally);
}

/* ----------------------------------------------------------------------
   no operations here, since compute results are stored in tally array
   just used by callers to indicate compute was used
   enables prediction of next step when update needs to tally
------------------------------------------------------------------------- */

void ComputeGasCollisionTally::compute_per_tally()
{
  invoked_per_tally = update->ntimestep;
}

/* ----------------------------------------------------------------------
   called by Collide before timesteps which invoke gas_tally()
---------------------------------------------------------------------- */

void ComputeGasCollisionTally::clear()
{
  cells = grid->cells;
  cinfo = grid->cinfo;
  
  ntally = 0;
}

/* ----------------------------------------------------------------------
   tally values for a single gas collision in icell,
     performing reaction (1 to N)
   iorig,jorig = particles before collision
   ip,jp = particles after collision
   jp = NULL means one particle after collision
   kp != NULL means three particles after collision
------------------------------------------------------------------------- */

void ComputeGasCollisionTally::gas_tally(int icell, int reaction,
                                         Particle::OnePart *iorig,
                                         Particle::OnePart *jorig,
                                         Particle::OnePart *ip,
                                         Particle::OnePart *jp,
                                         Particle::OnePart *kp)
{
  // skip if a reaction
  // this compute only tallies collisions that do not induce a reaction
  // reactions can be tallied by compute gas/reaction/tally command

  if (!reaction) return;
  
  // skip if icell not in grid group

  if (!(cinfo[icell].mask & groupbit)) return;

  // skip if either particle species not in mixture group

  int igroup = particle->mixture[imix]->species2group[iorig->ispecies];
  int jgroup = particle->mixture[imix]->species2group[jorig->ispecies];
  if (igroup < 0 || jgroup < 0) return;

  // grow tally array if necessary
  
  if (ntally == maxtally) grow_tally();

  // tally all values associated with group into array
  // particle iorig,ip have same collision point but before/after velocities
  
  double *vec = array_tally[ntally++];
  
  for (int m = 0; m < nvalue; m++) {
    switch (which[m]) {
    case ID1:
      vec[m] = ubuf(ip->id).d;
      break; 
    case ID2:
      vec[m] = ubuf(jp->id).d;
      break;
    case IDCELL:
      vec[m] = ubuf(cells[icell].id).d;
      break;
    case TYPE1:
      vec[m] = ubuf(ip->ispecies+1).d;
      break;
    case TYPE2:
      vec[m] = ubuf(jp->ispecies+1).d;
      break;
    case VX1PRE: 
      vec[m] = iorig->v[0];
      break;
    case VY1PRE: 
      vec[m] = iorig->v[1];
      break;
    case VZ1PRE: 
      vec[m] = iorig->v[2];
      break;
    case VX2PRE: 
      vec[m] = jorig->v[0];
      break;
    case VY2PRE: 
      vec[m] = jorig->v[1];
      break;
    case VZ2PRE: 
      vec[m] = jorig->v[2];
      break;
    case VX1POST:
      vec[m] = ip->v[0];
      break;
    case VY1POST: 
      vec[m] = ip->v[1];
      break;
    case VZ1POST: 
      vec[m] = ip->v[2];
      break;
    case VX2POST:
      vec[m] = jp->v[0];
      break;
    case VY2POST: 
      vec[m] = jp->v[1];
      break;
    case VZ2POST: 
      vec[m] = jp->v[2];
      break;
    }
  }
}

/* ----------------------------------------------------------------------
   return # of tallies
------------------------------------------------------------------------- */

int ComputeGasCollisionTally::tallyinfo(surfint *&dummy)
{
  return ntally;
}

/* ----------------------------------------------------------------------
   return datatype of tally quantity
   icol = 0 for vector
   icol = 1 to N for array column
   datatype = INT,DOUBLE,BIGINT
------------------------------------------------------------------------- */

int ComputeGasCollisionTally::datatype(int icol)
{
  if (which[icol-1] == ID1 || which[icol-1] == ID2) return INT;
  if (which[icol-1] == IDCELL) {
    if (sizeof(cellint) == sizeof(smallint)) return UINT;
    if (sizeof(cellint) == sizeof(bigint)) return BIGUINT;
  }
  if (which[icol-1] == TYPE1 || which[icol-1] == TYPE2) return INT;

  return DOUBLE;
}

/* ---------------------------------------------------------------------- */

void ComputeGasCollisionTally::grow_tally()
{
  maxtally += DELTA;
  memory->grow(array_tally,maxtally,nvalue,"gas/collision/tally:array_tally");
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

bigint ComputeGasCollisionTally::memory_usage()
{
  bigint bytes = 0;
  bytes += nvalue*maxtally * sizeof(double);    // array_tally
  return bytes;
}

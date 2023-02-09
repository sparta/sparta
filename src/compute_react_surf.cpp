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
#include "compute_react_surf.h"
#include "update.h"
#include "domain.h"
#include "comm.h"
#include "surf_react.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{REACTANT,PRODUCT};
#define DELTA 4096

/* ---------------------------------------------------------------------- */

ComputeReactSurf::ComputeReactSurf(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute react/surf command");

  int igroup = surf->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Compute react/surf group ID does not exist");
  groupbit = surf->bitmask[igroup];

  isr = surf->find_react(arg[3]);
  if (isr < 0) error->all(FLERR,"Compute react/surf reaction ID does not exist");

  ntotal = surf->sr[isr]->nlist;
  rpflag = 0;
  reaction2col = NULL;

  // parse per-column reactant/product args
  // reset rpflag = 1 and ntotal = # of args

  if (narg > 4) {
    rpflag = 1;
    int ncol = narg - 4;
    memory->create(reaction2col,ntotal,ncol,"react/surf:reaction2col");
    for (int i = 0; i < ntotal; i++)
      for (int j = 0; j < ncol; j++)
        reaction2col[i][j] = 0;
    int which;
    int icol = 0;
    int iarg = 4;
    while (iarg < narg) {
      if (strncmp(arg[iarg],"r:",2) == 0) which = REACTANT;
      else if (strncmp(arg[iarg],"p:",2) == 0) which = PRODUCT;
      else error->all(FLERR,"Illegal compute react/surf command");
      int n = strlen(&arg[iarg][2]) + 1;
      char *copy = new char[n];
      strcpy(copy,&arg[iarg][2]);
      char *ptr = copy;
      while ((ptr = strtok(ptr,"/")) != (char *) NULL) {
        for (int ireaction = 0; ireaction < ntotal; ireaction++) {
          reaction2col[ireaction][icol] = 0;
          if (which == REACTANT) {
            if (surf->sr[isr]->match_reactant(ptr,ireaction))
              reaction2col[ireaction][icol] = 1;
          } else if (which == PRODUCT) {
            if (surf->sr[isr]->match_product(ptr,ireaction))
              reaction2col[ireaction][icol] = 1;
          }
        }
        ptr = NULL;
      }
      delete [] copy;
      icol++;
      iarg++;
    }
    ntotal = narg - 4;
  }

  per_surf_flag = 1;
  size_per_surf_cols = ntotal;

  surf_tally_flag = 1;
  timeflag = 1;

  ntally = maxtally = 0;
  array_surf_tally = NULL;
  tally2surf = NULL;

  maxsurf = 0;
  array_surf = NULL;
  combined = 0;

  hash = new MyHash;

  dim = domain->dimension;
}

/* ---------------------------------------------------------------------- */

ComputeReactSurf::~ComputeReactSurf()
{
  memory->destroy(reaction2col);
  memory->destroy(array_surf_tally);
  memory->destroy(tally2surf);
  memory->destroy(array_surf);
  delete hash;
}

/* ---------------------------------------------------------------------- */

void ComputeReactSurf::init()
{
  if (!surf->exist)
    error->all(FLERR,"Cannot use compute react/surf when surfs do not exist");
  if (surf->implicit)
    error->all(FLERR,"Cannot use compute react/surf with implicit surfs");

  // warn if any surfs in group are assigned to different surf react model

  lines = surf->lines;
  tris = surf->tris;
  int nslocal = surf->nlocal;

  bigint flag = 0;
  if (dim == 2) {
    for (int i = 0; i < nslocal; i++) {
      if (!(lines[i].mask & groupbit)) return;
      if (lines[i].isr != isr) flag++;
    }
  } else {
    for (int i = 0; i < nslocal; i++) {
      if (!(tris[i].mask & groupbit)) return;
      if (tris[i].isr != isr) flag++;
    }
  }

  bigint flagall;
  if (surf->distributed) {
    MPI_Allreduce(&flag,&flagall,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  } else flagall = flag;

  if (flagall && comm->me == 0) {
    char str[128];
    sprintf(str,
            "Compute react/surf " BIGINT_FORMAT
            " surfs are not assigned to surf react model",flagall);
    error->warning(FLERR,str);
  }

  // initialize tally array in case accessed before a tally timestep

  clear();

  combined = 0;
}

/* ----------------------------------------------------------------------
   no operations here, since compute results are stored in tally array
   just used by callers to indicate compute was used
   enables prediction of next step when update needs to tally
------------------------------------------------------------------------- */

void ComputeReactSurf::compute_per_surf()
{
  invoked_per_surf = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void ComputeReactSurf::clear()
{
  lines = surf->lines;
  tris = surf->tris;

  // clear hash of tallied surf IDs
  // called by Update at beginning of timesteps surf tallying is done

  hash->clear();
  ntally = 0;
  combined = 0;
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

void ComputeReactSurf::surf_tally(int isurf, int icell, int reaction,
                                  Particle::OnePart *iorig,
                                  Particle::OnePart *ip, Particle::OnePart *jp)
{
  // skip if no reaction

  if (reaction == 0) return;
  reaction--;

  // skip if isurf not in surface group
  // or if this surf's reaction model is not a match

  if (dim == 2) {
    if (!(lines[isurf].mask & groupbit)) return;
    if (lines[isurf].isr != isr) return;
  } else {
    if (!(tris[isurf].mask & groupbit)) return;
    if (tris[isurf].isr != isr) return;
  }

  // itally = tally index of isurf
  // if 1st reaction on this isurf, add surf ID to hash
  // grow tally list if needed

  int itally;
  double *vec;

  surfint surfID;
  if (dim == 2) surfID = lines[isurf].id;
  else surfID = tris[isurf].id;

  if (hash->find(surfID) != hash->end()) itally = (*hash)[surfID];
  else {
    if (ntally == maxtally) grow_tally();
    itally = ntally;
    (*hash)[surfID] = itally;
    tally2surf[itally] = surfID;
    vec = array_surf_tally[itally];
    for (int i = 0; i < ntotal; i++) vec[i] = 0.0;
    ntally++;
  }

  // tally the reaction
  // for rpflag, tally each column if r2c is 1 for this reaction
  // for rpflag = 0, tally the reaction directly

  vec = array_surf_tally[itally];

  if (rpflag) {
    int *r2c = reaction2col[reaction];
    for (int i = 0; i < ntotal; i++)
      if (r2c[i]) vec[i] += 1.0;
  } else vec[reaction] += 1.0;
}

/* ----------------------------------------------------------------------
   return # of tallies and their indices into my local surf list
------------------------------------------------------------------------- */

int ComputeReactSurf::tallyinfo(surfint *&ptr)
{
  ptr = tally2surf;
  return ntally;
}

/* ----------------------------------------------------------------------
   sum tally values to owning surfs via surf->collate()
------------------------------------------------------------------------- */

void ComputeReactSurf::post_process_surf()
{
  if (combined) return;
  combined = 1;

  // reallocate array_surf if necessary

  int nown = surf->nown;

  if (nown > maxsurf) {
    memory->destroy(array_surf);
    maxsurf = nown;
    memory->create(array_surf,maxsurf,ntotal,"react/surf:array_surf");
  }

  // zero array_surf

  int i,j;
  for (i = 0; i < nown; i++)
    for (j = 0; j < ntotal; j++)
      array_surf[i][j] = 0.0;

  // collate entire array of results

  surf->collate_array(ntally,ntotal,tally2surf,array_surf_tally,array_surf);
}


/* ---------------------------------------------------------------------- */

void ComputeReactSurf::grow_tally()
{
  maxtally += DELTA;
  memory->grow(tally2surf,maxtally,"react/surf:tally2surf");
  memory->grow(array_surf_tally,maxtally,ntotal,"react/surf:array_surf_tally");
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

bigint ComputeReactSurf::memory_usage()
{
  bigint bytes = 0;
  bytes += ntotal*maxtally * sizeof(double);    // array_surf_tally
  bytes += maxtally * sizeof(surfint);          // tally2surf
  return bytes;
}

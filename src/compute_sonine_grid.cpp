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
#include "stdlib.h"
#include "string.h"
#include "compute_sonine_grid.h"
#include "particle.h"
#include "mixture.h"
#include "grid.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{AMOM,BMOM};
enum{X,Y,Z};

/* ---------------------------------------------------------------------- */

ComputeSonineGrid::ComputeSonineGrid(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal compute sonine/grid command");

  int igroup = grid->find_group(arg[2]);
  if (igroup < 0)
    error->all(FLERR,"Compute sonine/grid group ID does not exist");
  groupbit = grid->bitmask[igroup];

  imix = particle->find_mixture(arg[3]);
  if (imix < 0)
    error->all(FLERR,"Compute sonine/grid mixture ID does not exist");
  ngroup = particle->mixture[imix]->ngroup;

  // assume args are correct to infer nvalue, error check when process args

  nvalue = (narg-4) / 3;
  which = new int[nvalue];
  moment = new int[nvalue];
  order = new int[nvalue];

  noutpergroup = 0;

  int ivalue = 0;
  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"a") == 0) {
      if (iarg+3 > narg)
        error->all(FLERR,"Illegal compute sonine/grid command");
      which[ivalue] = AMOM;
      if (strcmp(arg[iarg+1],"x") == 0) moment[ivalue] = X;
      else if (strcmp(arg[iarg+1],"y") == 0) moment[ivalue] = Y;
      else if (strcmp(arg[iarg+1],"z") == 0) moment[ivalue] = Z;
      else error->all(FLERR,"Illegal compute sonine/grid command");
      order[ivalue] = atoi(arg[iarg+2]);
      if (order[ivalue] < 1 || order[ivalue] > 5)
        error->all(FLERR,"Illegal compute sonine/grid command");
      noutpergroup += order[ivalue];
      iarg += 3;
    } else if (strcmp(arg[iarg],"b") == 0) {
      if (iarg+3 > narg)
        error->all(FLERR,"Illegal compute sonine/grid command");
      which[ivalue] = BMOM;
      if (strcmp(arg[iarg+1],"xx") == 0) moment[ivalue] = 3*X + X;
      else if (strcmp(arg[iarg+1],"yy") == 0) moment[ivalue] = 3*Y + Y;
      else if (strcmp(arg[iarg+1],"zz") == 0) moment[ivalue] = 3*Z + Z;
      else if (strcmp(arg[iarg+1],"xy") == 0) moment[ivalue] = 3*X + Y;
      else if (strcmp(arg[iarg+1],"yz") == 0) moment[ivalue] = 3*Y + Z;
      else if (strcmp(arg[iarg+1],"xz") == 0) moment[ivalue] = 3*X + Z;
      else error->all(FLERR,"Illegal compute sonine/grid command");
      order[ivalue] = atoi(arg[iarg+2]);
      if (order[ivalue] < 1 || order[ivalue] > 5)
        error->all(FLERR,"Illegal compute sonine/grid command");
      noutpergroup += order[ivalue];
      iarg += 3;
    } else error->all(FLERR,"Illegal compute sonine/grid command");

    ivalue++;
  }

  // npergroup = # of outputs per group + 1 for mass tally
  // ntotal = total # of columns in tally array
  // reset_map() adjusts indices in initial map() using final npergroup

  npergroup = noutpergroup + 1;
  ntotal = ngroup*npergroup;

  per_grid_flag = 1;
  size_per_grid_cols = ngroup*noutpergroup;
  post_process_grid_flag = 1;

  // allocate and initialize nmap and map
  // map[1] = mass = first tally in each group
  // map[0] = numerator = tallies that follow mass, one per output

  nmap = 2;
  memory->create(map,ngroup*noutpergroup,nmap,"sonine/grid:map");
  for (int i = 0; i < ngroup; i++)
    for (int j = 0; j < noutpergroup; j++) {
      map[i*noutpergroup+j][0] = i*npergroup+j+1;
      map[i*noutpergroup+j][1] = i*npergroup;
    }

  nglocal = 0;
  vector_grid = NULL;
  tally = NULL;
  vcom = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeSonineGrid::~ComputeSonineGrid()
{
  if (copymode) return;

  delete [] which;
  delete [] moment;
  delete [] order;

  memory->destroy(map);

  memory->destroy(vector_grid);
  memory->destroy(tally);
  memory->destroy(vcom);
}

/* ---------------------------------------------------------------------- */

void ComputeSonineGrid::init()
{
  if (ngroup != particle->mixture[imix]->ngroup)
    error->all(FLERR,"Number of groups in compute sonine/grid "
               "mixture has changed");

  reallocate();
}

/* ---------------------------------------------------------------------- */

void ComputeSonineGrid::compute_per_grid()
{
  invoked_per_grid = update->ntimestep;

  Grid::ChildInfo *cinfo = grid->cinfo;
  Particle::Species *species = particle->species;
  Particle::OnePart *particles = particle->particles;
  int *s2g = particle->mixture[imix]->species2group;
  int nlocal = particle->nlocal;

  int i,j,k,m,n,ispecies,igroup,icell;
  double mass,norm,csq,value;
  double *v,*vec;
  double vthermal[3];

  // compute COM velocity on this timestep for each cell and group

  for (i = 0; i < nglocal; i++) {
    for (j = 0; j < ngroup; j++) {
      vcom[i][j][0] = 0.0;
      vcom[i][j][1] = 0.0;
      vcom[i][j][2] = 0.0;
      vcom[i][j][3] = 0.0;
    }
  }

  for (i = 0; i < nlocal; i++) {
    ispecies = particles[i].ispecies;
    igroup = s2g[ispecies];
    if (igroup < 0) continue;
    icell = particles[i].icell;
    if (!(cinfo[icell].mask & groupbit)) continue;

    mass = species[ispecies].mass;
    v = particles[i].v;

    vcom[icell][igroup][0] += mass * v[0];
    vcom[icell][igroup][1] += mass * v[1];
    vcom[icell][igroup][2] += mass * v[2];
    vcom[icell][igroup][3] += mass;
  }

  for (i = 0; i < nglocal; i++)
    for (j = 0; j < ngroup; j++) {
      norm = vcom[i][j][3];
      if (norm == 0.0) continue;
      vcom[i][j][0] /= norm;
      vcom[i][j][1] /= norm;
      vcom[i][j][2] /= norm;
    }

  // zero all accumulators - could do this with memset()

  for (i = 0; i < nglocal; i++)
    for (j = 0; j < ntotal; j++)
      tally[i][j] = 0.0;

  // loop over all particles, skip species not in mixture group
  // perform all tallies needed for each particle
  // mass is first tally of group
  // each output value adds one extra tally

  for (i = 0; i < nlocal; i++) {
    ispecies = particles[i].ispecies;
    igroup = s2g[ispecies];
    if (igroup < 0) continue;
    icell = particles[i].icell;
    if (!(cinfo[icell].mask & groupbit)) continue;

    vec = tally[icell];
    k = igroup*npergroup;

    mass = species[ispecies].mass;
    vec[k++] += mass;

    v = particles[i].v;
    vthermal[0] = v[0] - vcom[icell][igroup][0];
    vthermal[1] = v[1] - vcom[icell][igroup][1];
    vthermal[2] = v[2] - vcom[icell][igroup][2];
    csq = vthermal[0]*vthermal[0] + vthermal[1]*vthermal[1] +
      vthermal[2]*vthermal[2];

    for (m = 0; m < nvalue; m++) {
      if (which[m] == AMOM) {
        value = mass*vthermal[moment[m]] * csq;
        vec[k++] += value;
        for (n = 1; n < order[m]; n++) {
          value *= csq;
          vec[k++] += value;
        }
      } else {
        value = mass * vthermal[moment[m]/3] * vthermal[moment[m]%3] * csq;
        vec[k++] += value;
        for (n = 1; n < order[m]; n++) {
          value *= csq;
          vec[k++] += value;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   query info about internal tally array for this compute
   index = which column of output (0 for vec, 1 to N for array)
   return # of tally quantities for this index
   also return array = ptr to tally array
   also return cols = ptr to list of columns in tally for this index
------------------------------------------------------------------------- */

int ComputeSonineGrid::query_tally_grid(int index, double **&array, int *&cols)
{
  index--;
  array = tally;
  cols = map[index];
  return nmap;
}

/* ----------------------------------------------------------------------
   tally accumulated info to compute final normalized values
   index = which column of output (0 for vec, 1 to N for array)
   for etally = NULL:
     use internal tallied info for single timestep
     compute values for all grid cells
       store results in vector_grid with nstride = 1 (single col of array_grid)
   for etaylly = ptr to caller array:
     use external tallied info for many timesteps
     emap = list of etally columns to use, # of columns determined by index
     store results in caller's vec, spaced by nstride
   if norm = 0.0, set result to 0.0 directly so do not divide by 0.0
------------------------------------------------------------------------- */

void ComputeSonineGrid::
post_process_grid(int index, int nsample,
                  double **etally, int *emap, double *vec, int nstride)
{
  index--;

  int lo = 0;
  int hi = nglocal;
  int k = 0;

  if (!etally) {
    nsample = 1;
    etally = tally;
    emap = map[index];
    vec = vector_grid;
    nstride = 1;
  }

  // compute normalized final value for each grid cell

  double norm;
  int numerator = emap[0];
  int mass = emap[1];
  for (int icell = lo; icell < hi; icell++) {
    norm = etally[icell][mass];
    if (norm == 0.0) vec[k] = 0.0;
    else vec[k] = etally[icell][numerator] / norm;
    k += nstride;
  }
}

/* ----------------------------------------------------------------------
   reallocate arrays if nglocal has changed
   called by init() and load balancer
------------------------------------------------------------------------- */

void ComputeSonineGrid::reallocate()
{
  if (grid->nlocal == nglocal) return;

  memory->destroy(vector_grid);
  memory->destroy(tally);
  memory->destroy(vcom);
  nglocal = grid->nlocal;
  memory->create(vector_grid,nglocal,"sonine/grid:vector_grid");
  memory->create(tally,nglocal,ntotal,"sonine/grid:tally");
  memory->create(vcom,nglocal,ngroup,4,"sonine/grid:vcom");
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

bigint ComputeSonineGrid::memory_usage()
{
  bigint bytes;
  bytes = nglocal * sizeof(double);              // vector_grid
  bytes = ntotal*nglocal * sizeof(double);       // tally array
  bytes += nglocal*ngroup*4 * sizeof(double);    // vcom
  return bytes;
}

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

enum{THERMAL,AMOM,BMOM};
enum{X,Y,Z};
enum{NONE,COUNT,MASSWT,DOF};

/* ---------------------------------------------------------------------- */

ComputeSonineGrid::ComputeSonineGrid(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute sonine/grid command");

  imix = particle->find_mixture(arg[2]);
  if (imix < 0) 
    error->all(FLERR,"Compute sonine/grid mixture ID does not exist");

  int nmax = narg-3;
  which = new int[nmax];
  moment = new int[nmax];
  order = new int[nmax];
  nvalue = 0;
  npergroup = 0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"thermal") == 0) {
      which[nvalue] = THERMAL;
      nvalue++;
      npergroup++;
      iarg++;
    } else if (strcmp(arg[iarg],"a") == 0) {
      if (iarg+3 > narg) 
        error->all(FLERR,"Illegal compute sonine/grid command");
      which[nvalue] = AMOM;
      if (strcmp(arg[iarg+1],"x") == 0) moment[nvalue] = X;
      else if (strcmp(arg[iarg+1],"y") == 0) moment[nvalue] = Y;
      else if (strcmp(arg[iarg+1],"z") == 0) moment[nvalue] = Z;
      else error->all(FLERR,"Illegal compute sonine/grid command");
      order[nvalue] = atoi(arg[iarg+2]);
      if (order[nvalue] < 1 || order[nvalue] > 5)
        error->all(FLERR,"Illegal compute sonine/grid command");
      npergroup += order[nvalue];
      nvalue++;
      iarg += 3;
    } else if (strcmp(arg[iarg],"b") == 0) {
      if (iarg+3 > narg) 
        error->all(FLERR,"Illegal compute sonine/grid command");
      which[nvalue] = BMOM;
      if (strcmp(arg[iarg+1],"xx") == 0) moment[nvalue] = 3*X + X;
      else if (strcmp(arg[iarg+1],"yy") == 0) moment[nvalue] = 3*Y + Y;
      else if (strcmp(arg[iarg+1],"zz") == 0) moment[nvalue] = 3*Z + Z;
      else if (strcmp(arg[iarg+1],"xy") == 0) moment[nvalue] = 3*X + Y;
      else if (strcmp(arg[iarg+1],"yz") == 0) moment[nvalue] = 3*Y + Z;
      else if (strcmp(arg[iarg+1],"xz") == 0) moment[nvalue] = 3*X + Z;
      else error->all(FLERR,"Illegal compute sonine/grid command");
      order[nvalue] = atoi(arg[iarg+2]);
      if (order[nvalue] < 1 || order[nvalue] > 5)
        error->all(FLERR,"Illegal compute sonine/grid command");
      npergroup += order[nvalue];
      nvalue++;
      iarg += 3;
    } else error->all(FLERR,"Illegal compute sonine/grid command");
  }

  ngroup = particle->mixture[imix]->ngroup;
  ntotal = ngroup*npergroup;

  per_grid_flag = 1;
  size_per_grid_cols = ntotal;
  post_process_grid_flag = 1;

  nglocal = 0;
  vcom = NULL;
  masstot = NULL;
  sonine = array_grid = NULL;

  // norm vectors

  memory->create(value_norm_style,ngroup,npergroup,
                 "sonine/grid:value_norm_style");

  for (int i = 0; i < ngroup; i++) {
    int m = 0;
    for (int j = 0; j < nvalue; j++) {
      if (which[j] == THERMAL) value_norm_style[i][m++] = COUNT; 
      else
        for (int k = 0; k < order[j]; k++)
          value_norm_style[i][m++] = MASSWT;
    }
  }

  norm_count = new double*[ngroup];
  norm_mass = new double*[ngroup];
  for (int i = 0; i < ngroup; i++)
    norm_count[i] = norm_mass[i] = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeSonineGrid::~ComputeSonineGrid()
{
  delete [] which;
  delete [] moment;
  delete [] order;

  memory->destroy(vcom);
  memory->destroy(masstot);
  memory->destroy(sonine);

  memory->destroy(value_norm_style);

  for (int i = 0; i < ngroup; i++) {
    memory->destroy(norm_count[i]);
    memory->destroy(norm_mass[i]);
  }
  delete [] norm_count;
  delete [] norm_mass;
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
  double mvv2e = update->mvv2e;
  double tprefactor = mvv2e / (3.0*update->boltz);

  int i,j,k,m,n,ispecies,igroup,icell;
  double csq,mass,value;
  double *norm,*v,*vec;
  double vthermal[3];

  // compute COM velocity for each cell and group

  for (i = 0; i < nglocal; i++) {
    for (j = 0; j < ngroup; j++) {
      vcom[i][j][0] = 0.0;
      vcom[i][j][1] = 0.0;
      vcom[i][j][2] = 0.0;
      masstot[i][j] = 0.0;
    }
  }

  for (i = 0; i < nlocal; i++) {
    ispecies = particles[i].ispecies;
    igroup = s2g[ispecies];
    if (igroup < 0) continue;
    icell = particles[i].icell;

    mass = species[ispecies].mass;
    masstot[icell][igroup] += mass;

    v = particles[i].v;
    vcom[icell][igroup][0] += mass * v[0];
    vcom[icell][igroup][1] += mass * v[1];
    vcom[icell][igroup][2] += mass * v[2];
  }

  for (i = 0; i < nglocal; i++)
    for (j = 0; j < ngroup; j++) {
      if (masstot[i][j] == 0.0) continue;
      vcom[i][j][0] /= masstot[i][j];
      vcom[i][j][1] /= masstot[i][j];
      vcom[i][j][2] /= masstot[i][j];
    }
  
  // zero accumulator array and norm vectors

  for (i = 0; i < nglocal; i++)
    for (j = 0; j < ntotal; j++) sonine[i][j] = 0.0;

  for (j = 0; j < ngroup; j++) {
    if (norm = norm_count[j])
      for (i = 0; i < nglocal; i++) norm[i] = 0.0;
    if (norm = norm_mass[j])
      for (i = 0; i < nglocal; i++) norm[i] = 0.0;
  }

  // compute thermal temperature and sonine moments

  for (i = 0; i < nlocal; i++) {
    ispecies = particles[i].ispecies;
    igroup = s2g[ispecies];
    if (igroup < 0) continue;
    icell = particles[i].icell;

    mass = species[ispecies].mass;
    if (norm_mass[igroup]) norm_mass[igroup][icell] += mass;
    if (norm_count[igroup]) norm_count[igroup][icell] += 1.0;

    v = particles[i].v;
    vthermal[0] = v[0] - vcom[icell][igroup][0];
    vthermal[1] = v[1] - vcom[icell][igroup][1];
    vthermal[2] = v[2] - vcom[icell][igroup][2];
    csq = vthermal[0]*vthermal[0] + vthermal[1]*vthermal[1] + 
      vthermal[2]*vthermal[2];

    vec = sonine[icell];
    k = igroup*npergroup;
    
    for (m = 0; m < nvalue; m++) {
      switch (which[m]) {
      case THERMAL:
        vec[k++] += tprefactor*mass * csq;
        break;
      case AMOM:
        value = mass*vthermal[moment[m]] * csq;
        vec[k++] += value;
        for (n = 1; n < order[m]; n++) {
          value *= csq;
          vec[k++] += value;
        }
        break;
      case BMOM:
        value = mass * vthermal[moment[m]/3] * vthermal[moment[m]%3] * csq;
        vec[k++] += value;
        for (n = 1; n < order[m]; n++) {
          value *= csq;
          vec[k++] += value;
        }
        break;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   use tallied info to compute a normalized value
   icell = -1, return values for entire group = index
     store them in out vector with nstride
   icell >= 0, return value for single icell in group = index
     store it in out[0]
   index = which column of this compute's output is requested
     0 = vector, 1-N = columns of array
   called by dumps with NULL input arrays, so use internal array/norm as input
------------------------------------------------------------------------- */

void ComputeSonineGrid::post_process_grid(void *innumer, void *indenom,
                                          int icell, int index,
                                          double *out, int nstride)
{
  // use internal storage for input
  // should never be called with non-NULL inputs
  // fix ave/grid does its own normalization for this compute

  double **array,*norm;

  if (innumer == NULL) {
    array = array_grid;
    norm = normptr(index);
  } else error->all(FLERR,"Invalid call to "
                    "ComputeSonineGrid::post_process_grid()");

  // request for either a single value or entire column of values
  // for single value, iterate thru outer loop just once

  index--;
  int istart = icell;
  if (icell < 0) istart = 0;

  int m = 0;
  for (int i = istart; i < nglocal; i++) {
    if (norm[i] == 0.0) out[m] = 0.0;
    else out[m] = array[i][index]/norm[i];
    m += nstride;
    if (icell >= 0) return;
  }
}


/* ----------------------------------------------------------------------
   reallocate arrays if nglocal has changed
   called by init() and load balancer
------------------------------------------------------------------------- */

void ComputeSonineGrid::reallocate()
{
  if (grid->nlocal == nglocal) return;

  nglocal = grid->nlocal;
  memory->destroy(vcom);
  memory->destroy(masstot);
  memory->destroy(sonine);

  memory->create(vcom,nglocal,ngroup,3,"sonine/grid:vcom");
  memory->create(masstot,nglocal,ngroup,"sonine/grid:masstot");
  memory->create(sonine,nglocal,ntotal,"sonine/grid:sonine");
  array_grid = sonine;

  for (int i = 0; i < ngroup; i++) {
    memory->destroy(norm_count[i]);
    memory->destroy(norm_mass[i]);
    norm_count[i] = norm_mass[i] = NULL;
    for (int j = 0; j < nvalue; j++) {
      if (which[j] == THERMAL) {
        if (norm_count[i] == NULL)
          memory->create(norm_count[i],nglocal,"sonine/grid:norm_count");
      } else {
        if (norm_mass[i] == NULL)
          memory->create(norm_mass[i],nglocal,"sonine/grid:norm_mass");
      }
    }
  }
}

/* ----------------------------------------------------------------------
   return info for norm vector used by column N
   input N is value from 1 to Ncols
   output: istyle = NONE, COUNT, etc
   output: igroup = which group within style
------------------------------------------------------------------------- */

void ComputeSonineGrid::normwhich(int n, int &istyle, int &igroup)
{
  igroup = (n-1) / npergroup;
  int ivalue = (n-1) % npergroup;
  if (value_norm_style[igroup][ivalue] == COUNT) istyle = COUNT;
  else if (value_norm_style[igroup][ivalue] == MASSWT) istyle = MASSWT;
  else istyle = NONE;
}

/* ----------------------------------------------------------------------
   return ptr to norm vector used by column N
   input N is value from 1 to Ncols
------------------------------------------------------------------------- */

double *ComputeSonineGrid::normptr(int n)
{
  int igroup = (n-1) / npergroup;
  int ivalue = (n-1) % npergroup;
  if (value_norm_style[igroup][ivalue] == COUNT) return norm_count[igroup];
  if (value_norm_style[igroup][ivalue] == MASSWT) return norm_mass[igroup];
  return NULL;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

bigint ComputeSonineGrid::memory_usage()
{
  bigint bytes = 0;
  bytes += nglocal*ngroup*3 * sizeof(double);
  bytes += nglocal*ngroup * sizeof(int);
  bytes += nglocal*ntotal * sizeof(double);
  for (int i = 0; i < ngroup; i++) {
    if (norm_count[i]) bytes += nglocal * sizeof(double);
    if (norm_mass[i]) bytes += nglocal * sizeof(double);
  }
  return bytes;
}

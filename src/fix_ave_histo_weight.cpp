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

#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "fix_ave_histo_weight.h"
#include "particle.h"
#include "mixture.h"
#include "domain.h"
#include "region.h"
#include "grid.h"
#include "modify.h"
#include "compute.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{X,V,F,COMPUTE,FIX,VARIABLE};
enum{SCALAR,VECTOR,WINDOW};
enum{GLOBAL,PERPARTICLE,PERGRID};
enum{IGNORE,END,EXTRA};

#define INVOKED_SCALAR 1
#define INVOKED_VECTOR 2
#define INVOKED_ARRAY 4
#define INVOKED_PER_PARTICLE 8
#define INVOKED_PER_GRID 16

/* ---------------------------------------------------------------------- */

FixAveHistoWeight::FixAveHistoWeight(SPARTA *spa, int narg, char **arg) :
  FixAveHisto(spa, narg, arg)
{
  weightflag = 1;

  // nvalues = 2 required for histo/weight

  if (nvalues != 2) error->all(FLERR,"Illegal fix ave/histo/weight command");

  // check that length of 2 values is the same

  int size[2];

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == X || which[i] == V) {
      size[i] = particle->nlocal;
    } else if (which[i] == COMPUTE && kind == GLOBAL && mode == SCALAR) {
      int icompute = modify->find_compute(ids[i]);
      size[i] = modify->compute[icompute]->size_vector;
    } else if (which[i] == COMPUTE && kind == GLOBAL && mode == VECTOR) {
      int icompute = modify->find_compute(ids[i]);
      size[i] = modify->compute[icompute]->size_array_rows;
    } else if (which[i] == COMPUTE && kind == PERPARTICLE) {
      size[i] = particle->nlocal;
    } else if (which[i] == COMPUTE && kind == PERGRID) {
      size[i] = grid->nlocal;
    } else if (which[i] == FIX && kind == GLOBAL && mode == SCALAR) {
      int ifix = modify->find_fix(ids[i]);
      size[i] = modify->fix[ifix]->size_vector;
    } else if (which[i] == FIX && kind == GLOBAL && mode == VECTOR) {
      int ifix = modify->find_fix(ids[i]);
      size[i]= modify->fix[ifix]->size_array_rows;
    } else if (which[i] == FIX && kind == PERPARTICLE) {
      size[i] = particle->nlocal;
    } else if (which[i] == FIX && kind == PERGRID) {
      size[i] = grid->nlocal;
    } else if (which[i] == VARIABLE && kind == PERPARTICLE) {
      size[i] = particle->nlocal;
    } else if (which[i] == VARIABLE && kind == PERGRID) {
      size[i] = grid->nlocal;
    }
  }

  if (size[0] != size[1])
    error->all(FLERR,"Fix ave/histo/weight value and weight vector "
               "lengths do not match");

  vectorwt = NULL;
  maxvectorwt = 0;
}

/* ---------------------------------------------------------------------- */

FixAveHistoWeight::~FixAveHistoWeight()
{
  if (copymode) return;

  memory->destroy(vectorwt);
}

/* ---------------------------------------------------------------------- */

void FixAveHistoWeight::calculate_weights()
{
  // weight factors are 2nd value (i = 1)

  weight = 0.0;
  weights = NULL;
  stridewt = 0;

  int i = 1;
  int m = value2index[i];
  int j = argindex[i];

  // invoke compute if not previously invoked

  if (which[i] == COMPUTE) {

    Compute *compute = modify->compute[m];

    if (kind == GLOBAL && mode == SCALAR) {
      if (j == 0) {
        if (!(compute->invoked_flag & INVOKED_SCALAR)) {
          compute->compute_scalar();
          compute->invoked_flag |= INVOKED_SCALAR;
        }
        weight = compute->scalar;
      } else {
        if (!(compute->invoked_flag & INVOKED_VECTOR)) {
          compute->compute_vector();
          compute->invoked_flag |= INVOKED_VECTOR;
        }
        weight = compute->vector[j-1];
      }

    } else if (kind == GLOBAL && mode == VECTOR) {
      if (j == 0) {
        if (!(compute->invoked_flag & INVOKED_VECTOR)) {
          compute->compute_vector();
          compute->invoked_flag |= INVOKED_VECTOR;
        }
        weights = compute->vector;
        stridewt = 1;
      } else {
        if (!(compute->invoked_flag & INVOKED_ARRAY)) {
          compute->compute_array();
          compute->invoked_flag |= INVOKED_ARRAY;
        }
        if (compute->array) weights = &compute->array[0][j-1];
        stridewt = compute->size_array_cols;
      }

    } else if (kind == PERPARTICLE) {
      if (!(compute->invoked_flag & INVOKED_PER_PARTICLE)) {
        compute->compute_per_particle();
        compute->invoked_flag |= INVOKED_PER_PARTICLE;
      }
      if (j == 0) {
        weights = compute->vector_particle;
        stridewt = 1;
      } else if (compute->array_particle) {
        weights = &compute->array_particle[0][j-1];
        stridewt = compute->size_per_particle_cols;
      }

    } else if (kind == PERGRID) {
      if (!(compute->invoked_flag & INVOKED_PER_GRID)) {
        compute->compute_per_grid();
        compute->invoked_flag |= INVOKED_PER_GRID;
      }
      if (j == 0) {
        weights = compute->vector_grid;
        stridewt = 1;
      } else if (compute->array_grid) {
        weights = &compute->array_grid[0][j-1];
        stridewt = compute->size_per_grid_cols;
      }
    }

  // access fix fields, guaranteed to be ready

  } else if (which[i] == FIX) {

    Fix *fix = modify->fix[m];

    if (kind == GLOBAL && mode == SCALAR) {
      if (j == 0) weight = fix->compute_scalar();
      else weight = fix->compute_vector(j-1);

    } else if (kind == GLOBAL && mode == VECTOR) {
      error->all(FLERR,"Fix ave/histo/weight option not yet supported");
      // NOTE: need to allocate local storage
      if (j == 0) {
        int n = fix->size_vector;
        for (i = 0; i < n; i++) weights[n] = fix->compute_vector(i);
      } else {
        int n = fix->size_vector;
        for (i = 0; i < n; i++) weights[n] = fix->compute_array(i,j-1);
      }

    } else if (kind == PERPARTICLE) {
      if (j == 0) {
        weights = fix->vector_particle;
        stridewt = 1;
      } else if (fix->array_particle) {
        weights = fix->array_particle[j-1];
        stridewt = fix->size_per_particle_cols;
      }

    } else if (kind == PERGRID) {
      if (j == 0) {
        weights = fix->vector_grid;
        stridewt = 1;
      } else if (fix->array_grid) {
        weights = &fix->array_grid[0][j-1];
        stridewt = fix->size_per_grid_cols;
      }
    }

  // evaluate equal-style or particle-style or grid-style variable

  } else if (which[i] == VARIABLE) {
    if (kind == GLOBAL && mode == SCALAR) {
      weight = input->variable->compute_equal(m);

    } else if (which[i] == VARIABLE && kind == PERPARTICLE) {
      if (particle->maxlocal > maxvectorwt) {
        memory->destroy(vectorwt);
        maxvectorwt = particle->maxlocal;
        memory->create(vectorwt,maxvectorwt,"ave/histo/weight:vectorwt");
      }
      input->variable->compute_particle(m,vectorwt,1,0);
      weights = vectorwt;
      stridewt = 1;

    } else if (which[i] == VARIABLE && kind == PERGRID) {
      if (grid->maxlocal > maxvectorwt) {
        memory->destroy(vectorwt);
        maxvectorwt = grid->maxlocal;
        memory->create(vectorwt,maxvectorwt,"ave/histo/weight:vectorwt");
      }
      input->variable->compute_grid(m,vectorwt,1,0);
      weights = vectorwt;
      stridewt = 1;
    }

  // explicit per-particle attributes
  // NOTE: need to allocate local storage

  } else {
    error->all(FLERR,"Fix ave/histo/weight option not yet supported");
  }
}

/* ----------------------------------------------------------------------
   bin a single value with weight
------------------------------------------------------------------------- */

void FixAveHistoWeight::bin_one_weight(double value, double wt)
{
  stats[2] = MIN(stats[2],value);
  stats[3] = MAX(stats[3],value);

  if (value < lo) {
    if (beyond == IGNORE) {
      stats[1] += wt;
      return;
    } else bin[0] += wt;
  } else if (value > hi) {
    if (beyond == IGNORE) {
      stats[1] += wt;
      return;
    } else bin[nbins-1] += wt;
  } else {
    int ibin = static_cast<int> ((value-lo)*bininv);
    ibin = MIN(ibin,nbins-1);
    if (beyond == EXTRA) ibin++;
    bin[ibin] += wt;
  }

  stats[0] += wt;
}

/* ----------------------------------------------------------------------
   bin a single value with weight
------------------------------------------------------------------------- */

void FixAveHistoWeight::bin_one(double value)
{
  bin_one_weight(value,weight);
}

/* ----------------------------------------------------------------------
   bin a vector of values with weights
   values and weights each have a stride
------------------------------------------------------------------------- */

void FixAveHistoWeight::bin_vector(int n, double *values, int stride)
{
  int m = 0;
  int mwt = 0;
  for (int i = 0; i < n; i++) {
    bin_one_weight(values[m],weights[mwt]);
    m += stride;
    mwt += stridewt;
  }
}

/* ----------------------------------------------------------------------
   bin a per-particle attribute with weights
   index is 0,1,2 if attribute is X or V
   weights have a stridewt
------------------------------------------------------------------------- */

void FixAveHistoWeight::bin_particles(int attribute, int index)
{
  Particle::OnePart *particles = particle->particles;
  int *s2g = particle->mixture[imix]->species2group;
  int nlocal = particle->nlocal;

  Region *region;
  if (regionflag) region = domain->regions[iregion];

  int mwt = 0;

  if (attribute == X) {
    if (regionflag && mixflag) {
      for (int i = 0; i < nlocal; i++) {
        if (region->match(particles[i].x) &&
            s2g[particles[i].ispecies] < 0)
          bin_one_weight(particles[i].x[index],weights[mwt]);
        mwt += stridewt;
      }
    } else if (regionflag) {
      for (int i = 0; i < nlocal; i++) {
        if (region->match(particles[i].x))
          bin_one_weight(particles[i].x[index],weights[mwt]);
        mwt += stridewt;
      }
    } else if (mixflag) {
      for (int i = 0; i < nlocal; i++) {
        if (s2g[particles[i].ispecies] >= 0)
          bin_one_weight(particles[i].x[index],weights[mwt]);
        mwt += stridewt;
      }
    } else {
      for (int i = 0; i < nlocal; i++) {
        bin_one_weight(particles[i].x[index],weights[mwt]);
        mwt += stridewt;
      }
    }

  } else if (attribute == V) {
    if (regionflag && mixflag) {
      for (int i = 0; i < nlocal; i++) {
        if (region->match(particles[i].x) &&
            s2g[particles[i].ispecies] < 0)
          bin_one_weight(particles[i].v[index],weights[mwt]);
        mwt += stridewt;
      }
    } else if (regionflag) {
      for (int i = 0; i < nlocal; i++) {
        if (region->match(particles[i].x))
          bin_one_weight(particles[i].v[index],weights[mwt]);
        mwt += stridewt;
      }
    } else if (mixflag) {
      for (int i = 0; i < nlocal; i++) {
        if (s2g[particles[i].ispecies] >= 0)
          bin_one_weight(particles[i].v[index],weights[mwt]);
        mwt += stridewt;
      }
    } else {
      for (int i = 0; i < nlocal; i++) {
        bin_one_weight(particles[i].v[index],weights[mwt]);
        mwt += stridewt;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   bin a per-particle vector of values with weights
   values and weights each have a stride
------------------------------------------------------------------------- */

void FixAveHistoWeight::bin_particles(double *values, int stride)
{
  Particle::OnePart *particles = particle->particles;
  int *s2g = particle->mixture[imix]->species2group;
  int nlocal = particle->nlocal;

  Region *region;
  if (regionflag) region = domain->regions[iregion];

  int m = 0;
  int mwt = 0;

  if (regionflag && mixflag) {
    for (int i = 0; i < nlocal; i++) {
      if (region->match(particles[i].x) &&
          s2g[particles[i].ispecies] >= 0)
        bin_one_weight(values[m],weights[mwt]);
      m += stride;
      mwt += stridewt;
    }
  } else if (regionflag) {
    for (int i = 0; i < nlocal; i++) {
      if (region->match(particles[i].x))
        bin_one_weight(values[m],weights[mwt]);
      m += stride;
      mwt += stridewt;
    }
  } else if (mixflag) {
    for (int i = 0; i < nlocal; i++) {
      if (s2g[particles[i].ispecies] < 0)
        bin_one_weight(values[m],weights[mwt]);
      m += stride;
      mwt += stridewt;
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      bin_one_weight(values[m],weights[mwt]);
      m += stride;
      mwt += stridewt;
    }
  }
}

/* ----------------------------------------------------------------------
   bin a per-grid vector of values with weights
------------------------------------------------------------------------- */

void FixAveHistoWeight::bin_grid_cells(double *values, int stride)
{
  Grid::ChildInfo *cinfo = grid->cinfo;
  int nglocal = grid->nlocal;

  int m = 0;
  int mwt = 0;

  if (groupflag) {
    for (int i = 0; i < nglocal; i++) {
      if (cinfo[i].mask & groupbit)
        bin_one_weight(values[m],weights[mwt]);
      m += stride;
      mwt += stridewt;
    }
  } else {
    for (int i = 0; i < nglocal; i++) {
      bin_one_weight(values[m],weights[mwt]);
      m += stride;
      mwt += stridewt;
    }
  }
}

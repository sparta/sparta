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

/* ----------------------------------------------------------------------
   Contributing author: Tim Fuller (Sandia)
------------------------------------------------------------------------- */

#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "sparta_masks.h"
#include "kokkos_base.h"
#include "fix_ave_histo_weight_kokkos.h"
#include "particle_kokkos.h"
#include "update.h"
#include "particle.h"
#include "mixture.h"
#include "grid.h"
#include "grid_kokkos.h"
#include "domain.h"
#include "region.h"
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

FixAveHistoWeightKokkos::FixAveHistoWeightKokkos(SPARTA *spa, int narg, char **arg) :
  FixAveHistoKokkos(spa, narg, arg)
{

  weightflag = 1;

  // nvalues = 2 required for histo/weight

  if (nvalues != 2) error->all(FLERR,"Illegal fix ave/histo/weight/kokkos command");

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
    error->all(FLERR,"Fix ave/histo/weight/kokkos value and weight vector "
               "lengths do not match");

  vectorwt = NULL;
  maxvectorwt = 0;

  }

/* ---------------------------------------------------------------------- */

FixAveHistoWeightKokkos::~FixAveHistoWeightKokkos()
{
  if (copymode) return;
}

/* ---------------------------------------------------------------------- */
void FixAveHistoWeightKokkos::calculate_weights()
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
    if (!compute->kokkos_flag)
      error->all(FLERR,"Cannot (yet) use non-Kokkos computes with fix ave/histo/kk");
    KokkosBase* computeKKBase = dynamic_cast<KokkosBase*>(compute);

    if (kind == GLOBAL && mode == SCALAR) {
      if (j == 0) {
        if (!(compute->invoked_flag & INVOKED_SCALAR)) {
          compute->compute_scalar();
          compute->invoked_flag |= INVOKED_SCALAR;
        }
        weight = compute->scalar;
      } else {
        error->all(FLERR,"Compute not compatible with fix ave/histo/kk");
        if (!(compute->invoked_flag & INVOKED_VECTOR)) {
          compute->compute_vector();
          compute->invoked_flag |= INVOKED_VECTOR;
        }
        weight = compute->vector[j-1];
      }
    } else if (kind == GLOBAL && mode == VECTOR) {
      error->all(FLERR,"Compute not compatible with fix ave/histo/kk");
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
      error->all(FLERR,"Compute not compatible with fix ave/histo/kk");
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
        computeKKBase->compute_per_grid_kokkos();
        compute->invoked_flag |= INVOKED_PER_GRID;
      }
      if (j == 0) {
        d_weights = computeKKBase->d_vector;
      } else if (compute->array_grid) {
        d_weights = Kokkos::subview(computeKKBase->d_array_grid,Kokkos::ALL(),j-1);
      }
    }

  // access fix fields, guaranteed to be ready

  } else if (which[i] == FIX) {
    Fix *fix = modify->fix[m];
    if (!fix->kokkos_flag)
      error->all(FLERR,"Cannot (yet) use non-Kokkos fixes with fix ave/histo/weight/kk");
      KokkosBase* fixKKBase = dynamic_cast<KokkosBase*>(fix);

    if (kind == GLOBAL && mode == SCALAR) {
      if (j == 0) {
        weight = fix->compute_scalar();
      } else {
        error->all(FLERR,"Fix not compatible with fix ave/histo/weight/kk");
        weight = fix->compute_vector(j-1);
      }
    } else if (kind == GLOBAL && mode == VECTOR) {
      error->all(FLERR,"Fix not compatible with fix ave/histo/weight/kk");
      // NOTE: need to allocate local storage
      if (j == 0) {
        int n = fix->size_vector;
        for (i = 0; i < n; i++) weights[n] = fix->compute_vector(i);
      } else {
        int n = fix->size_vector;
        for (i = 0; i < n; i++) weights[n] = fix->compute_array(i,j-1);
      }

    } else if (kind == PERPARTICLE) {
      error->all(FLERR,"Fix not compatible with fix ave/histo/weight/kk");
      if (j == 0) {
        weights = fix->vector_particle;
        stridewt = 1;
      } else if (fix->array_particle) {
        weights = fix->array_particle[j-1];
        stridewt = fix->size_per_particle_cols;
      }

    } else if (kind == PERGRID) {
      if (j == 0) {
        d_weights = fixKKBase->d_vector;
      } else if (fixKKBase->d_array_grid.data()) {
        d_weights = Kokkos::subview(fixKKBase->d_array_grid,Kokkos::ALL(),j-1);
      }
    }

  // evaluate equal-style or particle-style or grid-style variable

  } else if (which[i] == VARIABLE) {
    error->all(FLERR,"Cannot (yet) use variables with fix ave/histo/weight/kk");
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
    printf("%d, %d\n", which[i] == VARIABLE, kind == PERGRID);
    error->all(FLERR,"Fix ave/histo/weight/kokkos option not yet supported");
  }
}

/* ----------------------------------------------------------------------
   bin a single value with weight
------------------------------------------------------------------------- */
void FixAveHistoWeightKokkos::bin_scalar(
    typename minmax_type::value_type& minmax,
    double value)
{
  bin_one(minmax, value, weight);
}

/* ----------------------------------------------------------------------
   bin a vector of values with stride
------------------------------------------------------------------------- */
void FixAveHistoWeightKokkos::bin_vector(
    minmax_type& reducer,
    int n, double *values, int stride)
{
  using FixKokkosDetails::mirror_view_from_raw_host_array;
  this->stride = stride;

  d_values = mirror_view_from_raw_host_array<double,DeviceType>(values, n, stride);

  auto policy = Kokkos::RangePolicy<TagFixAveHistoWeight_BinVector,DeviceType>(0, n);
  Kokkos::parallel_reduce(policy, *this, reducer);
}

/* ----------------------------------------------------------------------
   bin a per-particle attribute
   index is 0,1,2 if attribute is X or V
------------------------------------------------------------------------- */
void FixAveHistoWeightKokkos::bin_particles(
    minmax_type& reducer,
    int attribute, int index)
{
  using Kokkos::RangePolicy;
  using FixKokkosDetails::mirror_view_from_raw_host_array;

  this->index = index;
  int n = particle->nlocal;

  // FIXME: Kokkos version of region
  //Region *region;
  //if (regionflag) region = domain->regions[iregion];

  if (regionflag)
    error->all(FLERR,"Cannot (yet) use regionflag with fix ave/histo/kk");

  if (attribute == X) {

    if (regionflag && mixflag) {
      //auto policy = RangePolicy<TagFixAveHistoWeight_BinParticlesX1,DeviceType>(0, n);
      //Kokkos::parallel_reduce(policy, *this, reducer);
    } else if (regionflag) {
      //auto policy = RangePolicy<TagFixAveHistoWeight_BinParticlesX2,DeviceType>(0, n);
      //Kokkos::parallel_reduce(policy, *this, reducer);
    } else if (mixflag) {
      auto policy = RangePolicy<TagFixAveHistoWeight_BinParticlesX3,DeviceType>(0, n);
      Kokkos::parallel_reduce(policy, *this, reducer);
    } else {
      auto policy = RangePolicy<TagFixAveHistoWeight_BinParticlesX4,DeviceType>(0, n);
      Kokkos::parallel_reduce(policy, *this, reducer);
    }

  } else if (attribute == V) {

    if (regionflag && mixflag) {
      //auto policy = RangePolicy<TagFixAveHistoWeight_BinParticlesV1,DeviceType>(0, n);
      //Kokkos::parallel_reduce(policy, *this, reducer);
    } else if (regionflag) {
      //auto policy = RangePolicy<TagFixAveHistoWeight_BinParticlesV2,DeviceType>(0, n);
      //Kokkos::parallel_reduce(policy, *this, reducer);
    } else if (mixflag) {
      auto policy = RangePolicy<TagFixAveHistoWeight_BinParticlesV3,DeviceType>(0, n);
      Kokkos::parallel_reduce(policy, *this, reducer);
    } else {
      auto policy = RangePolicy<TagFixAveHistoWeight_BinParticlesV4,DeviceType>(0, n);
      Kokkos::parallel_reduce(policy, *this, reducer);
    }
  }
}

/* ----------------------------------------------------------------------
   bin a per-particle vector of values with stride
------------------------------------------------------------------------- */
void FixAveHistoWeightKokkos::bin_particles(
    minmax_type& reducer,
    double *values, int stride)
{
  using Kokkos::RangePolicy;
  using FixKokkosDetails::mirror_view_from_raw_host_array;

  this->stride = stride;
  int n = particle->nlocal;

  d_values = mirror_view_from_raw_host_array<double,DeviceType>(values, n, stride);

  // FIXME: Kokkos version of region
  // FIXME: Does values need to be made a view that lives on Device?
  //Region *region;
  //if (regionflag) region = domain->regions[iregion];

  if (regionflag)
    error->all(FLERR,"Cannot (yet) use regionflag with fix ave/histo/kk");

  if (regionflag && mixflag) {
    //auto policy = RangePolicy<TagFixAveHistoWeight_BinParticles1,DeviceType>(0, n);
    //Kokkos::parallel_reduce(policy, *this, reducer);
  } else if (regionflag) {
    //auto policy = RangePolicy<TagFixAveHistoWeight_BinParticles2,DeviceType>(0, n);
    //Kokkos::parallel_reduce(policy, *this, reducer);
  } else if (mixflag) {
    auto policy = RangePolicy<TagFixAveHistoWeight_BinParticles3,DeviceType>(0, n);
    Kokkos::parallel_reduce(policy, *this, reducer);
  } else {
    auto policy = RangePolicy<TagFixAveHistoWeight_BinParticles4,DeviceType>(0, n);
    Kokkos::parallel_reduce(policy, *this, reducer);
  }
}

/* ----------------------------------------------------------------------
   bin a per-grid vector of values with stride
------------------------------------------------------------------------- */
void FixAveHistoWeightKokkos::bin_grid_cells(
    minmax_type& reducer,
    DAT::t_float_1d_strided d_vec)
{
  using Kokkos::RangePolicy;
  using FixKokkosDetails::mirror_view_from_raw_host_array;

  int n = grid->nlocal;

  d_values = d_vec;

  if (groupflag) {
    GridKokkos* grid_kk = (GridKokkos*) grid;
    grid_kk->sync(Device, CINFO_MASK);
    auto policy = RangePolicy<TagFixAveHistoWeight_BinGridCells1,DeviceType>(0, n);
    Kokkos::parallel_reduce(policy, *this, reducer);
  } else {
    auto policy = RangePolicy<TagFixAveHistoWeight_BinGridCells2,DeviceType>(0, n);
    Kokkos::parallel_reduce(policy, *this, reducer);
  }
}

/* ------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void
FixAveHistoWeightKokkos::operator()(TagFixAveHistoWeight_BinVector, const int i,
                                    minmax_type::value_type& lminmax) const
{
  bin_one(lminmax, d_values(i), d_weights(i));
}

/* ------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void
FixAveHistoWeightKokkos::operator()(TagFixAveHistoWeight_BinParticles1, const int i,
                                    minmax_type::value_type& lminmax) const
{
  /*
   * region is not Kokkos compatible
   * If a Kokkos compatible region becomes available,
   * this code can be recommissioned.
   *
  const int ispecies = d_particles(i).ispecies;
  if (region_kk->match(d_particles(i).x) && d_s2g(imix, ispecies) >= 0)
  {
    bin_one(lminmax, d_values(i), d_weights(i));
  }
  */
}

/* ------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void
FixAveHistoWeightKokkos::operator()(TagFixAveHistoWeight_BinParticles2, const int i,
                                    minmax_type::value_type& lminmax) const
{
  /*
   * region is not Kokkos compatible.
   * If a Kokkos compatible region becomes available,
   * this code can be recommissioned.
   *
  if (region_kk->match(d_particles(i).x))
  {
    bin_one(lminmax, d_values(i), d_weights(i));
  }
  */
}

/* ------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void
FixAveHistoWeightKokkos::operator()(TagFixAveHistoWeight_BinParticles3, const int i,
                                    minmax_type::value_type& lminmax) const
{
  const int ispecies = d_particles(i).ispecies;
  if (d_s2g(imix, ispecies) < 0)
  {
    bin_one(lminmax, d_values(i), d_weights(i));
  }
}

/* ------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void
FixAveHistoWeightKokkos::operator()(TagFixAveHistoWeight_BinParticles4, const int i,
                                    minmax_type::value_type& lminmax) const
{
  bin_one(lminmax, d_values(i), d_weights(i));
}

/* ------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void
FixAveHistoWeightKokkos::operator()(TagFixAveHistoWeight_BinGridCells1, const int i,
                                    minmax_type::value_type& lminmax) const
{
  if (grid_kk->k_cinfo.d_view[i].mask & groupbit)
  {
    bin_one(lminmax, d_values(i), d_weights(i));
  }
}

/* ------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void
FixAveHistoWeightKokkos::operator()(TagFixAveHistoWeight_BinGridCells2, const int i,
                                    minmax_type::value_type& lminmax) const
{
  bin_one(lminmax, d_values(i), d_weights(i));
}

/* ------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void
FixAveHistoWeightKokkos::operator()(TagFixAveHistoWeight_BinParticlesX1, const int i,
                                    minmax_type::value_type& lminmax) const
{
  /*
   * region is not Kokkos compatible
   * If a Kokkos compatible region becomes available,
   * this code can be recommissioned.
   *
  const int ispecies = d_particles(i).ispecies;
  if (region_kk->match(d_particles(i).x) && d_s2g(imix, ispecies) < 0)
  {
    bin_one(lminmax, d_particles(i).x[index], d_weights(i));
  }
  */
}

/* ------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void
FixAveHistoWeightKokkos::operator()(TagFixAveHistoWeight_BinParticlesX2, const int i,
                                    minmax_type::value_type& lminmax) const
{
  /*
   * region is not Kokkos compatible
   * If a Kokkos compatible region becomes available,
   * this code can be recommissioned.
   *
  if (region_kk->match(d_particles(i).x))
  {
    bin_one(lminmax, d_particles(i).x[index], d_weights(i));
  }
  */
}

/* ------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void
FixAveHistoWeightKokkos::operator()(TagFixAveHistoWeight_BinParticlesX3, const int i,
                                    minmax_type::value_type& lminmax) const
{
  const int ispecies = d_particles(i).ispecies;
  if (d_s2g(imix, ispecies) >= 0)
  {
    bin_one(lminmax, d_particles(i).x[index], d_weights(i));
  }
}

/* ------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void
FixAveHistoWeightKokkos::operator()(TagFixAveHistoWeight_BinParticlesX4, const int i,
                                    minmax_type::value_type& lminmax) const
{
  bin_one(lminmax, d_particles(i).x[index], d_weights(i));
}

/* ------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void
FixAveHistoWeightKokkos::operator()(TagFixAveHistoWeight_BinParticlesV1, const int i,
                                    minmax_type::value_type& lminmax) const
{
  /*
   * region is not Kokkos compatible
   * If a Kokkos compatible region becomes available,
   * this code can be recommissioned.
   *
  const int ispecies = d_particles(i).ispecies;
  if (region_kk->match(d_particles(i).x) && d_s2g(imix, ispecies) < 0)
  {
    bin_one(lminmax, d_particles(i).v[index], d_weights(i));
  }
  */
}

/* ------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void
FixAveHistoWeightKokkos::operator()(TagFixAveHistoWeight_BinParticlesV2, const int i,
                                    minmax_type::value_type& lminmax) const
{
  /*
   * region is not Kokkos compatible
   * If a Kokkos compatible region becomes available,
   * this code can be recommissioned.
   *
  if (region_kk->match(d_particles(i).x))
  {
    bin_one(lminmax, d_particles(i).v[index], d_weights(i));
  }
  */
}

/* ------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void
FixAveHistoWeightKokkos::operator()(TagFixAveHistoWeight_BinParticlesV3, const int i,
                                    minmax_type::value_type& lminmax) const
{
  const int ispecies = d_particles(i).ispecies;
  if (d_s2g(imix, ispecies) >= 0)
  {
    bin_one(lminmax, d_particles(i).v[index], d_weights(i));
  }
}

/* ------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void
FixAveHistoWeightKokkos::operator()(TagFixAveHistoWeight_BinParticlesV4, const int i,
                                    minmax_type::value_type& lminmax) const
{
  bin_one(lminmax, d_particles(i).v[index], d_weights(i));
}

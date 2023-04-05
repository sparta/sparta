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
#include "fix_ave_histo_kokkos.h"
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
#include "memory_kokkos.h"
#include "error.h"

using namespace SPARTA_NS;

enum{X,V,F,COMPUTE,FIX,VARIABLE};
enum{ONE,RUNNING};
enum{SCALAR,VECTOR,WINDOW};
enum{GLOBAL,PERPARTICLE,PERGRID};
enum{IGNORE,END,EXTRA};

#define INVOKED_SCALAR 1
#define INVOKED_VECTOR 2
#define INVOKED_ARRAY 4
#define INVOKED_PER_PARTICLE 8
#define INVOKED_PER_GRID 16

/* ---------------------------------------------------------------------- */

FixAveHistoKokkos::FixAveHistoKokkos(SPARTA *spa, int narg, char **arg) :
  FixAveHisto(spa, narg, arg)
{
  kokkos_flag = 1;
  execution_space = Device;

  k_stats.resize(4);
  d_stats = k_stats.d_view;

  memory->destroy(bin);
  bin = NULL;
  memoryKK->grow_kokkos(k_bin, bin, nbins, "ave/histo:bin");
  d_bin = k_bin.d_view;

}

/* ---------------------------------------------------------------------- */

FixAveHistoKokkos::~FixAveHistoKokkos()
{
  if (copymode) return;

  k_stats = DAT::tdual_float_1d();
  memoryKK->destroy_kokkos(k_bin, bin);
}

/* ---------------------------------------------------------------------- */

void FixAveHistoKokkos::init()
{
  // set current indices for all computes,fixes,variables

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix ave/histo/kk does not exist");
      value2index[i] = icompute;

    } else if (which[i] == FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for fix ave/histo/kk does not exist");
      value2index[i] = ifix;

    } else if (which[i] == VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for fix ave/histo/kk does not exist");
      value2index[i] = ivariable;
    }
  }
}

/* ----------------------------------------------------------------------
   only does something if nvalid = current timestep
------------------------------------------------------------------------- */

void FixAveHistoKokkos::setup()
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveHistoKokkos::end_of_step()
{

  using FixKokkosDetails::mirror_view_from_raw_host_array;

  int j,m;

  // skip if not step which requires doing something

  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) {
    return;
  }

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device, PARTICLE_MASK|SPECIES_MASK);
  d_particles = particle_kk->k_particles.d_view;

  d_s2g = particle_kk->k_species2group.d_view;

  copymode = 1;

  // zero if first step
  if (irepeat == 0) {
    for (int i=0; i<4; i++) k_stats.h_view(i) = 0.0;
    k_stats.modify_host();
    k_stats.sync_device();

    for (int i=0; i<nbins; i++) k_bin.h_view(i) = 0.0;
    k_bin.modify_host();
    k_bin.sync_device();
  }

  minmax_type::value_type minmax;
  minmax_type reducer(minmax);

  // accumulate results of computes,fixes,variables to local copy
  // compute/fix/variable may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  // for fix ave/histo/weight, nvalues will be 2
  // first calculate weight factors, then histogram single value

  int ncount = nvalues;
  if (weightflag) {
    calculate_weights();
    ncount = 1;
  }

  for (int i = 0; i < ncount; i++) {
    m = value2index[i];
    j = argindex[i];

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
          bin_scalar(minmax, compute->scalar);
        }
        else {
          error->all(FLERR,"Compute kind not compatible with fix ave/histo/kk");
          if (!(compute->invoked_flag & INVOKED_VECTOR)) {
            compute->compute_vector();
            compute->invoked_flag |= INVOKED_VECTOR;
          }
          bin_scalar(minmax, compute->vector[j-1]);
        }
      } else if (kind == GLOBAL && mode == VECTOR) {
          error->all(FLERR,"Compute kind not compatible with fix ave/histo/kk");
        if (j == 0) {
          if (!(compute->invoked_flag & INVOKED_VECTOR)) {
            compute->compute_vector();
            compute->invoked_flag |= INVOKED_VECTOR;
          }
          bin_vector(reducer, compute->size_vector,compute->vector,1);
        } else {
          if (!(compute->invoked_flag & INVOKED_ARRAY)) {
            compute->compute_array();
            compute->invoked_flag |= INVOKED_ARRAY;
          }
          if (compute->array)
            bin_vector(reducer, compute->size_array_rows,&compute->array[0][j-1],
                       compute->size_array_cols);
        }
      } else if (kind == PERPARTICLE) {
          error->all(FLERR,"Compute kind not compatible with fix ave/histo/kk");
        if (!(compute->invoked_flag & INVOKED_PER_PARTICLE)) {
          compute->compute_per_particle();
          compute->invoked_flag |= INVOKED_PER_PARTICLE;
        }
        if (j == 0)
          bin_particles(reducer, compute->vector_particle,1);
        else if (compute->array_particle)
          bin_particles(reducer, &compute->array_particle[0][j-1],
                        compute->size_per_particle_cols);
      } else if (kind == PERGRID) {
        if (!(compute->invoked_flag & INVOKED_PER_GRID)) {
          computeKKBase->compute_per_grid_kokkos();
          compute->invoked_flag |= INVOKED_PER_GRID;
        }

        if (compute->post_process_grid_flag) {
          DAT::t_float_2d_lr d_etally;
          DAT::t_float_1d_strided d_vec;
          computeKKBase->post_process_grid_kokkos(j,1,d_etally,NULL,d_vec);
        }
        else if (compute->post_process_isurf_grid_flag)
          compute->post_process_isurf_grid();

        if (j == 0 || compute->post_process_grid_flag)
          bin_grid_cells(reducer, computeKKBase->d_vector);
        else if (computeKKBase->d_array_grid.data())
          // @stamoor: fix_ave_histo.cpp passes compute->array_grid[0][j-1],
          // @stamoor: so send subview of d_array_grid.
          bin_grid_cells(reducer,
                         Kokkos::subview(computeKKBase->d_array_grid,Kokkos::ALL(),j-1));
      }

    // access fix fields, guaranteed to be ready

    } else if (which[i] == FIX) {

      Fix *fix = modify->fix[m];
      if (!fix->kokkos_flag)
        error->all(FLERR,"Cannot (yet) use non-Kokkos fixes with fix ave/histo/kk");
      KokkosBase* fixKKBase = dynamic_cast<KokkosBase*>(fix);

      if (kind == GLOBAL && mode == SCALAR) {
        if (j == 0) {
          bin_scalar(minmax, fix->compute_scalar());
        }
        else {
          error->all(FLERR,"Fix not compatible with fix ave/histo/kk");
          bin_scalar(minmax, fix->compute_vector(j-1));
        }
      } else if (kind == GLOBAL && mode == VECTOR) {
        error->all(FLERR,"Fix not compatible with fix ave/histo/kk");
        if (j == 0) {
          int n = fix->size_vector;
          for (i = 0; i < n; i++) bin_scalar(minmax, fix->compute_vector(i));
        } else {
          int n = fix->size_vector;
          for (i = 0; i < n; i++) bin_scalar(minmax, fix->compute_array(i,j-1));
        }

      } else if (kind == PERPARTICLE) {
        error->all(FLERR,"Fix not compatible with fix ave/histo/kk");
        if (j == 0) bin_particles(reducer, fix->vector_particle,1);
        else if (fix->array_particle)
          bin_particles(reducer, fix->array_particle[j-1],fix->size_per_particle_cols);
      } else if (kind == PERGRID) {
        if (j == 0) {
          bin_grid_cells(reducer, fixKKBase->d_vector);
        } else if (fixKKBase->d_array_grid.data()) {
          // @stamoor: fix_ave_histo.cpp passes fix->array_grid[j-1], which is
          // not the same as what happens above with the compute object, it is
          // also inconsistent with fix_ave_histo_weight which uses
          // fix->array_grid[0][j-1] too.  Is this a type in fix_ave_histo?
          bin_grid_cells(reducer,
                         Kokkos::subview(fixKKBase->d_array_grid,Kokkos::ALL(),j-1));
        }
      }

    // evaluate equal-style or particle-style or grid-style variable
    } else if (which[i] == VARIABLE) {
      error->all(FLERR,"Cannot (yet) use variables with fix ave/histo/kk");
      if (kind == GLOBAL && mode == SCALAR) {
        bin_scalar(minmax, input->variable->compute_equal(m));

      } else if (which[i] == VARIABLE && kind == PERPARTICLE) {
        if (particle->maxlocal > maxvector) {
          memory->destroy(vector);
          maxvector = particle->maxlocal;
          memory->create(vector,maxvector,"ave/histo:vector");
        }
        input->variable->compute_particle(m,vector,1,0);
        bin_particles(reducer, vector,1);

      } else if (which[i] == VARIABLE && kind == PERGRID) {
        if (grid->maxlocal > maxvector) {
          memory->destroy(vector);
          maxvector = grid->maxlocal;
          memory->create(vector,maxvector,"ave/histo:vector");
        }
        input->variable->compute_grid(m,vector,1,0);
        //bin_grid_cells(reducer, vector);
      }
    } else {
      // explicit per-particle attributes
      bin_particles(reducer, which[i], j);
    }
  }

  k_stats.modify_device();
  k_stats.sync_host();

  k_bin.modify_device();
  k_bin.sync_host();

  // Copy data back
  stats[0] = k_stats.h_view(0);
  stats[1] = k_stats.h_view(1);
  stats[2] = minmax.min_val;
  stats[3] = minmax.max_val;

  // done if irepeat < nrepeat
  // else reset irepeat and nvalid

  irepeat++;
  if (irepeat < nrepeat) {
    nvalid += nevery;
    modify->addstep_compute(nvalid);
    copymode = 0;
    return;
  }

  irepeat = 0;
  nvalid = ntimestep + nfreq - (nrepeat-1)*nevery;
  modify->addstep_compute(nvalid);

  // merge histogram stats across procs if necessary
  if (kind == PERPARTICLE || kind == PERGRID) {
    MPI_Allreduce(stats,stats_all,2,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&stats[2],&stats_all[2],1,MPI_DOUBLE,MPI_MIN,world);
    MPI_Allreduce(&stats[3],&stats_all[3],1,MPI_DOUBLE,MPI_MAX,world);
    MPI_Allreduce(bin,bin_all,nbins,MPI_DOUBLE,MPI_SUM,world);

    stats[0] = stats_all[0];
    stats[1] = stats_all[1];
    stats[2] = stats_all[2];
    stats[3] = stats_all[3];
    for (int i = 0; i < nbins; i++) bin[i] = bin_all[i];
  }

  // if ave = ONE, only single Nfreq timestep value is needed
  // if ave = RUNNING, combine with all previous Nfreq timestep values
  // if ave = WINDOW, combine with nwindow most recent Nfreq timestep values

  if (ave == ONE) {
    stats_total[0] = stats[0];
    stats_total[1] = stats[1];
    stats_total[2] = stats[2];
    stats_total[3] = stats[3];
    for (int i = 0; i < nbins; i++) bin_total[i] = bin[i];

  } else if (ave == RUNNING) {
    stats_total[0] += stats[0];
    stats_total[1] += stats[1];
    stats_total[2] = MIN(stats_total[2],stats[2]);
    stats_total[3] = MAX(stats_total[3],stats[3]);
    for (int i = 0; i < nbins; i++) bin_total[i] += bin[i];

  } else if (ave == WINDOW) {
    stats_total[0] += stats[0];
    if (window_limit) stats_total[0] -= stats_list[iwindow][0];
    stats_list[iwindow][0] = stats[0];
    stats_total[1] += stats[1];
    if (window_limit) stats_total[1] -= stats_list[iwindow][1];
    stats_list[iwindow][1] = stats[1];

    if (window_limit) m = nwindow;
    else m = iwindow+1;

    stats_list[iwindow][2] = stats[2];
    stats_total[2] = stats_list[0][2];
    for (int i = 1; i < m; i++)
      stats_total[2] = MIN(stats_total[2],stats_list[i][2]);
    stats_list[iwindow][3] = stats[3];
    stats_total[3] = stats_list[0][3];
    for (int i = 1; i < m; i++)
      stats_total[3] = MAX(stats_total[3],stats_list[i][3]);

    for (int i = 0; i < nbins; i++) {
      bin_total[i] += bin[i];
      if (window_limit) bin_total[i] -= bin_list[iwindow][i];
      bin_list[iwindow][i] = bin[i];
    }

    iwindow++;
    if (iwindow == nwindow) {
      iwindow = 0;
      window_limit = 1;
    }
  }

  // output result to file

  if (fp && me == 0) {
    clearerr(fp);
    if (overwrite) fseek(fp,filepos,SEEK_SET);
    fprintf(fp,BIGINT_FORMAT " %d %g %g %g %g\n",ntimestep,nbins,
            stats_total[0],stats_total[1],stats_total[2],stats_total[3]);
    if (stats_total[0] != 0.0) {
      for (int i = 0; i < nbins; i++) {
        fprintf(fp,"%d %g %g %g\n",
                i+1,coord[i],bin_total[i],bin_total[i]/stats_total[0]);
      }
    } else {
      for (int i = 0; i < nbins; i++)
        fprintf(fp,"%d %g %g %g\n",i+1,coord[i],0.0,0.0);
    }

    if (ferror(fp)) {
      error->one(FLERR,"Error writing out histogram data");
    }

    fflush(fp);
    if (overwrite) {
      long fileend = ftell(fp);
      if (fileend > 0) ftruncate(fileno(fp),fileend);
    }
  }
  copymode = 0;
}

/* ----------------------------------------------------------------------
   return Ith vector value
------------------------------------------------------------------------- */

double FixAveHistoKokkos::compute_vector(int i)
{
  return stats_total[i];
}

/* ----------------------------------------------------------------------
   return I,J array value
------------------------------------------------------------------------- */

double FixAveHistoKokkos::compute_array(int i, int j)
{
  if (j == 0) return coord[i];
  else if (j == 1) return bin_total[i];
  else if (stats_total[0] != 0.0) return bin_total[i]/stats_total[0];
  return 0.0;
}

/* ----------------------------------------------------------------------
   bin a Scalar
------------------------------------------------------------------------- */
void FixAveHistoKokkos::bin_scalar(minmax_type::value_type& minmax, double val)
{
  bin_one(minmax, val);
}

/* ----------------------------------------------------------------------
   bin a vector of values with stride
------------------------------------------------------------------------- */
void FixAveHistoKokkos::bin_vector(
    minmax_type& reducer,
    int n, double *values, int stride)
{
  using FixKokkosDetails::mirror_view_from_raw_host_array;
  this->stride = stride;

  d_values = mirror_view_from_raw_host_array<double,DeviceType>(values, n, stride);

  auto policy = Kokkos::RangePolicy<TagFixAveHisto_BinVector,DeviceType>(0, n);
  Kokkos::parallel_reduce(policy, *this, reducer);
}

/* ----------------------------------------------------------------------
   bin a per-particle attribute
   index is 0,1,2 if attribute is X or V
------------------------------------------------------------------------- */
void FixAveHistoKokkos::bin_particles(
    minmax_type& reducer,
    int attribute, int index)
{
  using Kokkos::RangePolicy;

  this->index = index;
  int n = particle->nlocal;

  // FIXME: Kokkos version of region
  //Region *region;
  //if (regionflag) region = domain->regions[iregion];

  if (regionflag)
    error->all(FLERR,"Cannot (yet) use regionflag with fix ave/histo/kk");

  if (attribute == X) {

    if (regionflag && mixflag) {
      //auto policy = RangePolicy<TagFixAveHisto_BinParticlesX1,DeviceType>(0, n);
      //Kokkos::parallel_reduce(policy, *this, reducer);
    } else if (regionflag) {
      //auto policy = RangePolicy<TagFixAveHisto_BinParticlesX2,DeviceType>(0, n);
      //Kokkos::parallel_reduce(policy, *this, reducer);
    } else if (mixflag) {
      auto policy = RangePolicy<TagFixAveHisto_BinParticlesX3,DeviceType>(0, n);
      Kokkos::parallel_reduce(policy, *this, reducer);
    } else {
      auto policy = RangePolicy<TagFixAveHisto_BinParticlesX4,DeviceType>(0, n);
      Kokkos::parallel_reduce(policy, *this, reducer);
    }

  } else if (attribute == V) {

    if (regionflag && mixflag) {
      //auto policy = RangePolicy<TagFixAveHisto_BinParticlesV1,DeviceType>(0, n);
      //Kokkos::parallel_reduce(policy, *this, reducer);
    } else if (regionflag) {
      //auto policy = RangePolicy<TagFixAveHisto_BinParticlesV2,DeviceType>(0, n);
      //Kokkos::parallel_reduce(policy, *this, reducer);
    } else if (mixflag) {
      auto policy = RangePolicy<TagFixAveHisto_BinParticlesV3,DeviceType>(0, n);
      Kokkos::parallel_reduce(policy, *this, reducer);
    } else {
      auto policy = RangePolicy<TagFixAveHisto_BinParticlesV4,DeviceType>(0, n);
      Kokkos::parallel_reduce(policy, *this, reducer);
    }
  }
}

/* ----------------------------------------------------------------------
   bin a per-particle vector of values with stride
------------------------------------------------------------------------- */
void FixAveHistoKokkos::bin_particles(
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
    //auto policy = RangePolicy<TagFixAveHisto_BinParticles1,DeviceType>(0, n);
    //Kokkos::parallel_reduce(policy, *this, reducer);
  } else if (regionflag) {
    //auto policy = RangePolicy<TagFixAveHisto_BinParticles2,DeviceType>(0, n);
    //Kokkos::parallel_reduce(policy, *this, reducer);
  } else if (mixflag) {
    auto policy = RangePolicy<TagFixAveHisto_BinParticles3,DeviceType>(0, n);
    Kokkos::parallel_reduce(policy, *this, reducer);
  } else {
    auto policy = RangePolicy<TagFixAveHisto_BinParticles4,DeviceType>(0, n);
    Kokkos::parallel_reduce(policy, *this, reducer);
  }
}

/* ----------------------------------------------------------------------
   bin a per-grid vector of values with stride
------------------------------------------------------------------------- */
void FixAveHistoKokkos::bin_grid_cells(
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
    auto policy = RangePolicy<TagFixAveHisto_BinGridCells1,DeviceType>(0, n);
    Kokkos::parallel_reduce(policy, *this, reducer);
  } else {
    auto policy = RangePolicy<TagFixAveHisto_BinGridCells2,DeviceType>(0, n);
    Kokkos::parallel_reduce(policy, *this, reducer);
  }
}


/* ----------------------------------------------------------------------
   calculate nvalid = next step on which end_of_step does something
   can be this timestep if multiple of nfreq and nrepeat = 1
   else backup from next multiple of nfreq
   startstep is lower bound on nfreq multiple
------------------------------------------------------------------------- */

bigint FixAveHistoKokkos::nextvalid()
{
  bigint nvalid = (update->ntimestep/nfreq)*nfreq + nfreq;
  while (nvalid < startstep) nvalid += nfreq;
  if (nvalid-nfreq == update->ntimestep && nrepeat == 1)
    nvalid = update->ntimestep;
  else
    nvalid -= (nrepeat-1)*nevery;
  if (nvalid < update->ntimestep) nvalid += nfreq;
  return nvalid;
}

/* ------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void
FixAveHistoKokkos::operator()(TagFixAveHisto_BinVector, const int i,
                              minmax_type::value_type& lminmax) const
{
  bin_one(lminmax, d_values(i));
}

/* ------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void
FixAveHistoKokkos::operator()(TagFixAveHisto_BinParticles1, const int i,
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
    bin_one(lminmax, d_values(i));
  }
  */
}

/* ------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void
FixAveHistoKokkos::operator()(TagFixAveHisto_BinParticles2, const int i,
                              minmax_type::value_type& lminmax) const
{
  /*
   * region is not Kokkos compatible.
   * If a Kokkos compatible region becomes available,
   * this code can be recommissioned.
   *
  if (region_kk->match(d_particles(i).x))
  {
    bin_one(lminmax, d_values(i));
  }
  */
}

/* ------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void
FixAveHistoKokkos::operator()(TagFixAveHisto_BinParticles3, const int i,
                              minmax_type::value_type& lminmax) const
{
  const int ispecies = d_particles(i).ispecies;
  if (d_s2g(imix, ispecies) < 0)
  {
    bin_one(lminmax, d_values(i));
  }
}

/* ------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void
FixAveHistoKokkos::operator()(TagFixAveHisto_BinParticles4, const int i,
                              minmax_type::value_type& lminmax) const
{
  bin_one(lminmax, d_values(i));
}

/* ------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void
FixAveHistoKokkos::operator()(TagFixAveHisto_BinGridCells1, const int i,
                              minmax_type::value_type& lminmax) const
{
  if (grid_kk->k_cinfo.d_view[i].mask & groupbit)
  {
    bin_one(lminmax, d_values(i));
  }
}

/* ------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void
FixAveHistoKokkos::operator()(TagFixAveHisto_BinGridCells2, const int i,
                              minmax_type::value_type& lminmax) const
{
  bin_one(lminmax, d_values(i));
}

/* ------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void
FixAveHistoKokkos::operator()(TagFixAveHisto_BinParticlesX1, const int i,
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
    bin_one(lminmax, d_particles(i).x[index]);
  }
  */
}

/* ------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void
FixAveHistoKokkos::operator()(TagFixAveHisto_BinParticlesX2, const int i,
                              minmax_type::value_type& lminmax) const
{
  /*
   * region is not Kokkos compatible
   * If a Kokkos compatible region becomes available,
   * this code can be recommissioned.
   *
  if (region_kk->match(d_particles(i).x))
  {
    bin_one(lminmax, d_particles(i).x[index]);
  }
  */
}

/* ------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void
FixAveHistoKokkos::operator()(TagFixAveHisto_BinParticlesX3, const int i,
                              minmax_type::value_type& lminmax) const
{
  const int ispecies = d_particles(i).ispecies;
  if (d_s2g(imix, ispecies) >= 0)
  {
    bin_one(lminmax, d_particles(i).x[index]);
  }
}

/* ------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void
FixAveHistoKokkos::operator()(TagFixAveHisto_BinParticlesX4, const int i,
                              minmax_type::value_type& lminmax) const
{
  bin_one(lminmax, d_particles(i).x[index]);
}

/* ------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void
FixAveHistoKokkos::operator()(TagFixAveHisto_BinParticlesV1, const int i,
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
    bin_one(lminmax, d_particles(i).v[index]);
  }
  */
}

/* ------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void
FixAveHistoKokkos::operator()(TagFixAveHisto_BinParticlesV2, const int i,
                              minmax_type::value_type& lminmax) const
{
  /*
   * region is not Kokkos compatible
   * If a Kokkos compatible region becomes available,
   * this code can be recommissioned.
   *
  if (region_kk->match(d_particles(i).x))
  {
    bin_one(lminmax, d_particles(i).v[index]);
  }
  */
}

/* ------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void
FixAveHistoKokkos::operator()(TagFixAveHisto_BinParticlesV3, const int i,
                              minmax_type::value_type& lminmax) const
{
  const int ispecies = d_particles(i).ispecies;
  if (d_s2g(imix, ispecies) >= 0)
  {
    bin_one(lminmax, d_particles(i).v[index]);
  }
}

/* ------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void
FixAveHistoKokkos::operator()(TagFixAveHisto_BinParticlesV4, const int i,
                              minmax_type::value_type& lminmax) const
{
  bin_one(lminmax, d_particles(i).v[index]);
}

/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "compute_reduce_kokkos.h"
#include "update.h"
#include "domain.h"
#include "particle_kokkos.h"
#include "mixture.h"
#include "grid_kokkos.h"
#include "surf_kokkos.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "input.h"
#include "variable.h"
#include "memory_kokkos.h"
#include "error.h"
#include "sparta_masks.h"
#include "kokkos.h"

using namespace SPARTA_NS;

// these must stay in lock step with the enums in compute_reduce.cpp

enum{SUM,SUMSQ,MINN,MAXX,AVE,AVESQ,SUMAREA,AVEAREA};
enum{X,V,KE,EROT,EVIB,COMPUTE,FIX,VARIABLE,PCUSTOM,GCUSTOM,SCUSTOM};
enum{PARTICLE,GRID,SURF};

enum{INT,DOUBLE};                       // several files

#define INVOKED_PER_PARTICLE 8
#define INVOKED_PER_GRID 16
#define INVOKED_PER_SURF 32

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

ComputeReduceKokkos::ComputeReduceKokkos(SPARTA *sparta, int narg, char **arg) :
  ComputeReduce(sparta, narg, arg)
{
  kokkos_flag = 1;

  nelements = maxelements = 0;
}

/* ---------------------------------------------------------------------- */

ComputeReduceKokkos::~ComputeReduceKokkos()
{
  // no explicit deallocation: the scratch Kokkos views free themselves.
  // deliberately NOT guarded with "if (copymode) return;": nothing here is
  //   ever copied into a Kokkos functor (see the comment in the header), and
  //   a guard would only paper over the fact that the base class destructor
  //   would double-free anyway if a copy were ever made
}

/* ---------------------------------------------------------------------- */

void ComputeReduceKokkos::init()
{
  // everything this style needs is set up by the base class:
  //   s2g (mixture map), gridgroupbit, smasks/areasurf/area_total for SURF,
  //   and value2index for every input.  the device-side subset test reads a
  //   fresh copy of s2g in build_include(), so nothing extra is recorded here

  ComputeReduce::init();
}

/* ---------------------------------------------------------------------- */

double ComputeReduceKokkos::compute_scalar()
{
  if (sparta->kokkos->prewrap) return ComputeReduce::compute_scalar();

  invoked_scalar = update->ntimestep;

  double one = compute_one_kokkos(0,-1);

  // MPI reduction is identical to ComputeReduce::compute_scalar()

  if (mode == SUM || mode == SUMSQ || mode == SUMAREA) {
    MPI_Allreduce(&one,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  } else if (mode == MINN) {
    MPI_Allreduce(&one,&scalar,1,MPI_DOUBLE,MPI_MIN,world);
  } else if (mode == MAXX) {
    MPI_Allreduce(&one,&scalar,1,MPI_DOUBLE,MPI_MAX,world);
  } else if (mode == AVE || mode == AVESQ) {
    MPI_Allreduce(&one,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
    bigint n = count_included_kokkos();
    if (n) scalar /= n;
  } else if (mode == AVEAREA) {
    MPI_Allreduce(&one,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
    if (area_total > 0.0) scalar /= area_total;
  }

  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeReduceKokkos::compute_vector()
{
  if (sparta->kokkos->prewrap) {
    ComputeReduce::compute_vector();
    return;
  }

  invoked_vector = update->ntimestep;

  for (int m = 0; m < nvalues; m++)
    if (!replace || replace[m] < 0) {
      onevec[m] = compute_one_kokkos(m,-1);
      indices[m] = index;
    }

  // MPI reduction is identical to ComputeReduce::compute_vector()

  if (mode == SUM || mode == SUMSQ || mode == SUMAREA) {
    for (int m = 0; m < nvalues; m++)
      MPI_Allreduce(&onevec[m],&vector[m],1,MPI_DOUBLE,MPI_SUM,world);

  } else if (mode == MINN) {
    if (!replace) {
      for (int m = 0; m < nvalues; m++)
        MPI_Allreduce(&onevec[m],&vector[m],1,MPI_DOUBLE,MPI_MIN,world);

    } else {
      for (int m = 0; m < nvalues; m++)
        if (replace[m] < 0) {
          pairme.value = onevec[m];
          pairme.proc = me;
          MPI_Allreduce(&pairme,&pairall,1,MPI_DOUBLE_INT,MPI_MINLOC,world);
          vector[m] = pairall.value;
          owner[m] = pairall.proc;
        }
      for (int m = 0; m < nvalues; m++)
        if (replace[m] >= 0) {
          if (me == owner[replace[m]])
            vector[m] = compute_one_kokkos(m,indices[replace[m]]);
          MPI_Bcast(&vector[m],1,MPI_DOUBLE,owner[replace[m]],world);
        }
    }

  } else if (mode == MAXX) {
    if (!replace) {
      for (int m = 0; m < nvalues; m++)
        MPI_Allreduce(&onevec[m],&vector[m],1,MPI_DOUBLE,MPI_MAX,world);

    } else {
      for (int m = 0; m < nvalues; m++)
        if (replace[m] < 0) {
          pairme.value = onevec[m];
          pairme.proc = me;
          MPI_Allreduce(&pairme,&pairall,1,MPI_DOUBLE_INT,MPI_MAXLOC,world);
          vector[m] = pairall.value;
          owner[m] = pairall.proc;
        }
      for (int m = 0; m < nvalues; m++)
        if (replace[m] >= 0) {
          if (me == owner[replace[m]])
            vector[m] = compute_one_kokkos(m,indices[replace[m]]);
          MPI_Bcast(&vector[m],1,MPI_DOUBLE,owner[replace[m]],world);
        }
    }

  } else if (mode == AVE || mode == AVESQ) {
    bigint n = count_included_kokkos();
    for (int m = 0; m < nvalues; m++) {
      MPI_Allreduce(&onevec[m],&vector[m],1,MPI_DOUBLE,MPI_SUM,world);
      if (n) vector[m] /= n;
    }

  } else if (mode == AVEAREA) {
    for (int m = 0; m < nvalues; m++) {
      MPI_Allreduce(&onevec[m],&vector[m],1,MPI_DOUBLE,MPI_SUM,world);
      if (area_total > 0.0) vector[m] /= area_total;
    }
  }
}

/* ----------------------------------------------------------------------
   bring the host-side arrays that ComputeReduce::compute_one() reads back
     up to date before delegating one input to it.
   in a Kokkos run the authoritative copy of particles / cinfo / custom
     attributes is the device one, so a host fall-back that skipped this
     would silently reduce stale data (the subset test in particular reads
     particle->particles[i].ispecies and grid->cinfo[i].mask directly)
------------------------------------------------------------------------- */

void ComputeReduceKokkos::sync_host_for_fallback(int m)
{
  if (flavor[m] == PARTICLE) {
    ParticleKokkos *particle_kk = (ParticleKokkos *) particle;
    particle_kk->sync(Host,PARTICLE_MASK|SPECIES_MASK|CUSTOM_MASK);
  } else if (flavor[m] == GRID) {
    GridKokkos *grid_kk = (GridKokkos *) grid;
    grid_kk->sync(Host,CINFO_MASK|CUSTOM_MASK);
  } else if (surf->exist) {
    SurfKokkos *surf_kk = (SurfKokkos *) surf;
    surf_kk->sync(Host,LINE_MASK|TRI_MASK|MYLINE_MASK|MYTRI_MASK|CUSTOM_MASK);
  }
}

/* ----------------------------------------------------------------------
   calculate reduced value for one input M and return it
   same contract as ComputeReduce::compute_one():
     flag = -1: sum/min/max/ave all values, set index for MIN/MAX
     flag >= 0: simply return the value of element flag
   falls back to the host implementation whenever input M has no
     device-resident source (see setup_values)
------------------------------------------------------------------------- */

double ComputeReduceKokkos::compute_one_kokkos(int m, int flag)
{
  index = -1;

  // SURF inputs stay on the host.  KokkosBase exposes no d_vector_surf or
  //   d_array_surf, the surf->nown ownership decomposition and the smasks[]
  //   /areasurf[] weights that ComputeReduce::init() builds live only on the
  //   host, and SUMAREA/AVEAREA are per-surf by construction.  Reducing on
  //   the host is also trivially bit-for-bit identical

  if (flavor[m] == SURF) {
    sync_host_for_fallback(m);
    return ComputeReduce::compute_one(m,flag);
  }

  // gather element values for input m into d_values; 0 means no device source

  if (!setup_values(m)) {
    sync_host_for_fallback(m);
    return ComputeReduce::compute_one(m,flag);
  }

  // flag >= 0: return the single element, ignoring the subset mask,
  //   exactly as the host does

  if (flag >= 0) {
    double one = 0.0;
    if (flag < nelements) {

      // the scalar handed to deep_copy(value,View) is a non-deduced
      //   parameter, so it has to be spelled with the view's own value type,
      //   SPARTA_FLOAT, not double, or template deduction fails

      SPARTA_FLOAT tmp = 0.0;
      Kokkos::deep_copy(tmp,Kokkos::subview(d_values,flag));
      one = tmp;
    }
    return one;
  }

  return reduce_values();
}

/* ----------------------------------------------------------------------
   set up the device source for input m and gather it into d_values
   return 1 if the gather was done on the device, 0 to fall back to the host
   an early 0 return must not have invoked the source compute/fix, so that
     ComputeReduce::compute_one() can still invoke it itself
------------------------------------------------------------------------- */

int ComputeReduceKokkos::setup_values(int m)
{
  const int vidx = value2index[m];
  const int aidx = argindex[m];
  const int acol = aidx - 1;

  // particle-style / grid-style variables are evaluated by the host Variable
  //   class, there is no device evaluator.  Copying the host result to the
  //   device only to reduce it there would buy nothing and would give up the
  //   host's summation order, so reduce on the host instead.
  //   decided first, so that nothing at all has been touched on this path

  if (which[m] == VARIABLE) return 0;

  // explicit per-particle attributes: always resident on the device

  if (which[m] == X || which[m] == V || which[m] == KE ||
      which[m] == EROT || which[m] == EVIB) {

    ParticleKokkos *particle_kk = (ParticleKokkos *) particle;
    particle_kk->sync(Device,PARTICLE_MASK|SPECIES_MASK);
    auto l_particles = particle_kk->k_particles.view_device();
    auto l_species = particle_kk->k_species.view_device();

    build_include(m);
    if (nelements == 0) return 1;

    auto l_values = d_values;
    const int j = aidx;
    const double mvv2e = update->mvv2e;

    if (which[m] == X) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,nelements),
        KOKKOS_LAMBDA(const int i) {
          l_values(i) = l_particles(i).x[j];
        });

    } else if (which[m] == V) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,nelements),
        KOKKOS_LAMBDA(const int i) {
          l_values(i) = l_particles(i).v[j];
        });

    } else if (which[m] == KE) {

      // expression order is copied verbatim from ComputeReduce::compute_one()
      //   so the per-particle value is bit-identical to the host's

      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,nelements),
        KOKKOS_LAMBDA(const int i) {
          const double *v = l_particles(i).v;
          l_values(i) = mvv2e * 0.5 * l_species(l_particles(i).ispecies).mass *
            (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        });

    } else if (which[m] == EROT) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,nelements),
        KOKKOS_LAMBDA(const int i) {
          l_values(i) = l_particles(i).erot;
        });

    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,nelements),
        KOKKOS_LAMBDA(const int i) {
          l_values(i) = l_particles(i).evib;
        });
    }

    return 1;

  // per-particle or per-grid output of another compute

  } else if (which[m] == COMPUTE) {
    Compute *c = modify->compute[vidx];

    // a non-Kokkos compute, or one that does not derive from KokkosBase,
    //   has no device output at all.  neither test invokes anything

    if (!c->kokkos_flag) return 0;
    KokkosBase *ckk = dynamic_cast<KokkosBase*>(c);
    if (!ckk) return 0;

    DAT::t_float_1d_strided d_src;

    if (flavor[m] == PARTICLE) {

      // KokkosBase declares no compute_per_particle_kokkos() entry point, so
      //   the compute is invoked through the normal Compute API.  The /kk
      //   per-particle styles dispatch compute_per_particle() to their own
      //   device kernel whenever Kokkos is past its prewrap phase

      if (!(c->invoked_flag & INVOKED_PER_PARTICLE)) {
        c->compute_per_particle();
        c->invoked_flag |= INVOKED_PER_PARTICLE;
      }

      // CAVEAT: if the compute filled only its host vector_particle/
      //   array_particle we fall back, and ComputeReduce::compute_one() will
      //   see invoked_flag already set and read those host arrays without
      //   re-invoking.  That is correct only for a /kk style whose
      //   compute_per_particle() really did leave the HOST arrays current
      //   (e.g. one that just calls its non-Kokkos base).  A style that
      //   computes on the device and does not publish d_vector_particle /
      //   d_array_particle would be reduced from stale host data here.  No
      //   such style exists today (compute ke/particle/kk publishes
      //   d_vector_particle, fix field/particle/kk publishes
      //   d_array_particle) and KokkosBase gives us no way to force a
      //   sync_host() on someone else's DualView, so this is a documented
      //   limitation rather than a fixed problem

      if (aidx == 0) {
        if (!ckk->d_vector_particle.data()) return 0;
        if ((int) ckk->d_vector_particle.extent(0) < particle->nlocal) return 0;
        d_src = ckk->d_vector_particle;
      } else {
        if (!ckk->d_array_particle.data()) return 0;
        if ((int) ckk->d_array_particle.extent(0) < particle->nlocal ||
            (int) ckk->d_array_particle.extent(1) <= acol) return 0;
        d_src = Kokkos::subview(ckk->d_array_particle,Kokkos::ALL(),acol);
      }

    } else {                                            // GRID

      // post_process_isurf_grid() has no device counterpart, so this whole
      //   input is reduced on the host.  Checked before the invocation, so
      //   the host path is free to invoke the compute itself

      if (c->post_process_isurf_grid_flag) return 0;

      // canonical "consume another Kokkos compute's per-grid device output"
      //   sequence, same as FixAveGridKokkos::end_of_step() and
      //   FixAveHistoKokkos::end_of_step()

      if (!(c->invoked_flag & INVOKED_PER_GRID)) {
        ckk->compute_per_grid_kokkos();
        c->invoked_flag |= INVOKED_PER_GRID;
      }

      // must run for every input, not just the first: the column selects
      //   which post-processed quantity lands in d_vector_grid.  passing
      //   null views is the device spelling of the host's
      //   post_process_grid(aidx,1,NULL,NULL,NULL,1)

      if (c->post_process_grid_flag)
        ckk->post_process_grid_kokkos(aidx,1,DAT::t_float_2d_lr(),NULL,
                                      DAT::t_float_1d_strided());

      // a post-processing compute writes its answer to d_vector_grid
      //   regardless of aidx, matching the host's cvec/carray choice

      if (aidx == 0 || c->post_process_grid_flag) {
        if (!ckk->d_vector_grid.data()) return 0;
        if ((int) ckk->d_vector_grid.extent(0) < grid->nlocal) return 0;
        d_src = ckk->d_vector_grid;
      } else {
        if (!ckk->d_array_grid.data()) return 0;
        if ((int) ckk->d_array_grid.extent(0) < grid->nlocal ||
            (int) ckk->d_array_grid.extent(1) <= acol) return 0;
        d_src = Kokkos::subview(ckk->d_array_grid,Kokkos::ALL(),acol);
      }
    }

    build_include(m);
    gather_float(d_src);
    return 1;

  // per-particle or per-grid output of a fix
  // fixes are not invoked here, they are guaranteed to be up to date

  } else if (which[m] == FIX) {
    Fix *fix = modify->fix[vidx];

    // preserve the host's timestep compatibility check, before any fallback

    if (flavor[m] == PARTICLE) {
      if (update->ntimestep % fix->per_particle_freq)
        error->all(FLERR,"Fix used in compute reduce not "
                   "computed at compatible time");
    } else {
      if (update->ntimestep % fix->per_grid_freq)
        error->all(FLERR,"Fix used in compute reduce not "
                   "computed at compatible time");
    }

    if (!fix->kokkos_flag) return 0;
    KokkosBase *fkk = dynamic_cast<KokkosBase*>(fix);
    if (!fkk) return 0;

    DAT::t_float_1d_strided d_src;

    if (flavor[m] == PARTICLE) {
      if (aidx == 0) {
        if (!fkk->d_vector_particle.data()) return 0;
        if ((int) fkk->d_vector_particle.extent(0) < particle->nlocal) return 0;
        d_src = fkk->d_vector_particle;
      } else {
        if (!fkk->d_array_particle.data()) return 0;
        if ((int) fkk->d_array_particle.extent(0) < particle->nlocal ||
            (int) fkk->d_array_particle.extent(1) <= acol) return 0;
        d_src = Kokkos::subview(fkk->d_array_particle,Kokkos::ALL(),acol);
      }
    } else {
      if (aidx == 0) {
        if (!fkk->d_vector_grid.data()) return 0;
        if ((int) fkk->d_vector_grid.extent(0) < grid->nlocal) return 0;
        d_src = fkk->d_vector_grid;
      } else {
        if (!fkk->d_array_grid.data()) return 0;
        if ((int) fkk->d_array_grid.extent(0) < grid->nlocal ||
            (int) fkk->d_array_grid.extent(1) <= acol) return 0;
        d_src = Kokkos::subview(fkk->d_array_grid,Kokkos::ALL(),acol);
      }
    }

    build_include(m);
    gather_float(d_src);
    return 1;

  // per-particle custom attribute
  // index the plain host ewhich[]/etype[] arrays, exactly as the non-Kokkos
  //   path does, and take the device half of the matching DualView

  } else if (which[m] == PCUSTOM) {
    ParticleKokkos *particle_kk = (ParticleKokkos *) particle;
    particle_kk->sync(Device,CUSTOM_MASK);

    const int ew = particle->ewhich[vidx];
    if (ew < 0) return 0;

    if (particle->etype[vidx] == INT) {
      if (aidx == 0) {
        if (ew >= (int) particle_kk->k_eivec.view_host().extent(0)) return 0;
        auto d_src = particle_kk->k_eivec.view_host()[ew].k_view.view_device();
        if (!d_src.data()) return 0;
        if ((int) d_src.extent(0) < particle->nlocal) return 0;
        build_include(m);
        gather_int_vec(d_src);
      } else {
        if (ew >= (int) particle_kk->k_eiarray.view_host().extent(0)) return 0;
        auto d_src = particle_kk->k_eiarray.view_host()[ew].k_view.view_device();
        if (!d_src.data()) return 0;
        if ((int) d_src.extent(0) < particle->nlocal ||
            (int) d_src.extent(1) <= acol) return 0;
        build_include(m);
        gather_int_array(d_src,acol);
      }
    } else {
      if (aidx == 0) {
        if (ew >= (int) particle_kk->k_edvec.view_host().extent(0)) return 0;
        auto d_src = particle_kk->k_edvec.view_host()[ew].k_view.view_device();
        if (!d_src.data()) return 0;
        if ((int) d_src.extent(0) < particle->nlocal) return 0;
        build_include(m);
        gather_float(d_src);
      } else {
        if (ew >= (int) particle_kk->k_edarray.view_host().extent(0)) return 0;
        auto d_src = particle_kk->k_edarray.view_host()[ew].k_view.view_device();
        if (!d_src.data()) return 0;
        if ((int) d_src.extent(0) < particle->nlocal ||
            (int) d_src.extent(1) <= acol) return 0;
        build_include(m);
        gather_float(Kokkos::subview(d_src,Kokkos::ALL(),acol));
      }
    }

    return 1;

  // per-grid custom attribute

  } else if (which[m] == GCUSTOM) {
    GridKokkos *grid_kk = (GridKokkos *) grid;
    grid_kk->sync(Device,CUSTOM_MASK);

    const int ew = grid->ewhich[vidx];
    if (ew < 0) return 0;

    if (grid->etype[vidx] == INT) {
      if (aidx == 0) {
        if (ew >= (int) grid_kk->k_eivec.view_host().extent(0)) return 0;
        auto d_src = grid_kk->k_eivec.view_host()[ew].k_view.view_device();
        if (!d_src.data()) return 0;
        if ((int) d_src.extent(0) < grid->nlocal) return 0;
        build_include(m);
        gather_int_vec(d_src);
      } else {
        if (ew >= (int) grid_kk->k_eiarray.view_host().extent(0)) return 0;
        auto d_src = grid_kk->k_eiarray.view_host()[ew].k_view.view_device();
        if (!d_src.data()) return 0;
        if ((int) d_src.extent(0) < grid->nlocal ||
            (int) d_src.extent(1) <= acol) return 0;
        build_include(m);
        gather_int_array(d_src,acol);
      }
    } else {
      if (aidx == 0) {
        if (ew >= (int) grid_kk->k_edvec.view_host().extent(0)) return 0;
        auto d_src = grid_kk->k_edvec.view_host()[ew].k_view.view_device();
        if (!d_src.data()) return 0;
        if ((int) d_src.extent(0) < grid->nlocal) return 0;
        build_include(m);
        gather_float(d_src);
      } else {
        if (ew >= (int) grid_kk->k_edarray.view_host().extent(0)) return 0;
        auto d_src = grid_kk->k_edarray.view_host()[ew].k_view.view_device();
        if (!d_src.data()) return 0;
        if ((int) d_src.extent(0) < grid->nlocal ||
            (int) d_src.extent(1) <= acol) return 0;
        build_include(m);
        gather_float(Kokkos::subview(d_src,Kokkos::ALL(),acol));
      }
    }

    return 1;
  }

  // SCUSTOM and anything else: host.  (SCUSTOM is flavor SURF and has
  //   already been intercepted by compute_one_kokkos)

  return 0;
}

/* ----------------------------------------------------------------------
   copy the selected source into the contiguous d_values scratch view
   build_include() must have run first: it sets nelements and (re)sizes
     d_values, so taking the local copy afterwards is mandatory
------------------------------------------------------------------------- */

void ComputeReduceKokkos::gather_float(DAT::t_float_1d_strided d_src)
{
  if (nelements == 0) return;
  auto l_values = d_values;
  auto l_src = d_src;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,nelements),
    KOKKOS_LAMBDA(const int i) {
      l_values(i) = l_src(i);
    });
}

/* ---------------------------------------------------------------------- */

void ComputeReduceKokkos::gather_int_vec(DAT::t_int_1d d_src)
{
  if (nelements == 0) return;
  auto l_values = d_values;
  auto l_src = d_src;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,nelements),
    KOKKOS_LAMBDA(const int i) {
      l_values(i) = l_src(i);
    });
}

/* ---------------------------------------------------------------------- */

void ComputeReduceKokkos::gather_int_array(DAT::t_int_2d_lr d_src, int acol)
{
  if (nelements == 0) return;
  auto l_values = d_values;
  auto l_src = d_src;
  const int j = acol;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,nelements),
    KOKKOS_LAMBDA(const int i) {
      l_values(i) = l_src(i,j);
    });
}

/* ----------------------------------------------------------------------
   set nelements, size the scratch views, and fill d_include with the subset
     membership test.  the test only depends on flavor and subsetID, not on
     which input it is
------------------------------------------------------------------------- */

void ComputeReduceKokkos::build_include(int m)
{
  if (flavor[m] == PARTICLE) nelements = particle->nlocal;
  else nelements = grid->nlocal;

  grow_scratch(nelements);
  if (nelements == 0) return;

  // no subset: every element participates

  if (!subsetID) {
    Kokkos::deep_copy(Kokkos::subview(d_include,
                                      Kokkos::make_pair(0,nelements)),1);
    return;
  }

  auto l_include = d_include;

  if (flavor[m] == PARTICLE) {
    ParticleKokkos *particle_kk = (ParticleKokkos *) particle;
    particle_kk->sync(Device,PARTICLE_MASK);
    auto l_particles = particle_kk->k_particles.view_device();

    // refresh the device copy of the base class's s2g every time.
    //   ParticleKokkos::k_species2group is only built once, inside
    //   wrap_kokkos(), so it can be both stale and too short if species or
    //   mixtures change after the first setup -- indexing it on the device
    //   would then read out of bounds.  s2g always has particle->nspecies
    //   entries (Mixture::init) and is the exact array the non-Kokkos path
    //   reads, so copying it is both cheap and definitionally in agreement

    const int nsp = particle->nspecies;
    if ((int) k_s2g.view_host().extent(0) != nsp)
      MemKK::realloc_kokkos(k_s2g,"reduce/kk:s2g",nsp);
    auto hv_s2g = k_s2g.view_host();
    for (int i = 0; i < nsp; i++) hv_s2g(i) = s2g[i];
    k_s2g.modify_host();
    k_s2g.sync_device();

    auto l_s2g = k_s2g.view_device();
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,nelements),
      KOKKOS_LAMBDA(const int i) {
        l_include(i) = (l_s2g(l_particles(i).ispecies) >= 0) ? 1 : 0;
      });

  } else {
    GridKokkos *grid_kk = (GridKokkos *) grid;
    grid_kk->sync(Device,CINFO_MASK);
    auto l_cinfo = grid_kk->k_cinfo.view_device();
    const int l_groupbit = gridgroupbit;
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,nelements),
      KOKKOS_LAMBDA(const int i) {
        l_include(i) = (l_cinfo(i).mask & l_groupbit) ? 1 : 0;
      });
  }
}

/* ----------------------------------------------------------------------
   reduce d_values over the included elements according to mode
   also sets index for MIN/MAX when the replace option is in use
------------------------------------------------------------------------- */

double ComputeReduceKokkos::reduce_values()
{
  double one = 0.0;
  if (mode == MINN) one = BIG;
  else if (mode == MAXX) one = -BIG;

  if (nelements == 0) return one;

  const int n = nelements;

#ifdef SPARTA_KOKKOS_EXACT

  // SPARTA_KOKKOS_EXACT asks the Kokkos package to reproduce the non-Kokkos
  //   result bit for bit so the KOKKOS build can be regression tested
  //   against the existing gold-standard logs.  Floating point addition is
  //   not associative, so no tree reduction can be relied on to reproduce
  //   ComputeReduce::compute_one()'s left-to-right accumulation over
  //   ascending element index.  Bring the gathered values back to the host
  //   and run the base class's own combine() over them in exactly that
  //   order.  MIN/MAX would in principle survive a tree reduction, but they
  //   go through here as well: combine() then also sets "index" the way the
  //   host sets it, including its first-index-wins tie break and its
  //   treatment of +/-0.0, so there is nothing left to argue about

  // local mirrors on purpose: View::HostMirror is deprecation-gated in
  //   recent Kokkos, and this path only runs in a SPARTA_KOKKOS_EXACT
  //   regression build where the extra allocation does not matter

  auto h_values = Kokkos::create_mirror_view(d_values);
  auto h_include = Kokkos::create_mirror_view(d_include);
  Kokkos::deep_copy(h_values,d_values);
  Kokkos::deep_copy(h_include,d_include);

  for (int i = 0; i < n; i++) {
    if (!h_include(i)) continue;
    combine(one,h_values(i),i);
  }
  return one;

#else

  auto l_values = d_values;
  auto l_include = d_include;

  if (mode == SUM || mode == AVE || mode == SUMAREA || mode == AVEAREA) {

    // NOTE: summation order is the Kokkos backend's, not the host's.  For
    //   SUM/AVE that difference is not observable except in the last bits;
    //   see the SPARTA_KOKKOS_EXACT branch above for the bit-exact path

    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,n),
      KOKKOS_LAMBDA(const int i, double &lsum) {
        if (l_include(i)) lsum += l_values(i);
      },one);

  } else if (mode == SUMSQ || mode == AVESQ) {
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,n),
      KOKKOS_LAMBDA(const int i, double &lsum) {
        if (l_include(i)) lsum += l_values(i)*l_values(i);
      },one);

  } else if (mode == MINN) {

    // min/max are order independent, so the device result matches the host.
    //   Kokkos ignores the value a reducer is constructed with and starts
    //   from its own identity (+/- DBL_MAX), so the result has to be folded
    //   back into the host's BIG seed rather than used directly.  That also
    //   covers the empty / all-excluded case

    double vmin = BIG;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,n),
      KOKKOS_LAMBDA(const int i, double &lmin) {
        if (l_include(i) && l_values(i) < lmin) lmin = l_values(i);
      },Kokkos::Min<double>(vmin));
    if (vmin < one) one = vmin;

  } else if (mode == MAXX) {
    double vmax = -BIG;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,n),
      KOKKOS_LAMBDA(const int i, double &lmax) {
        if (l_include(i) && l_values(i) > lmax) lmax = l_values(i);
      },Kokkos::Max<double>(vmax));
    if (vmax > one) one = vmax;
  }

  // ComputeReduce::combine() also records the index of the winning element.
  //   It is only ever consumed by the "replace" option, which in turn is only
  //   legal for min/max, so the extra pass is skipped otherwise.
  //   The host's strict < / > leaves the FIRST winning index behind, so take
  //   the smallest index attaining the winning value.  Kokkos::MinLoc cannot
  //   be used for this: its join only breaks a tie against the identity, so
  //   the index it returns depends on the thread decomposition.
  //   "seeded" reproduces the host leaving index = -1 when nothing ever beat
  //   its BIG / -BIG seed

  if (replace && (mode == MINN || mode == MAXX)) {
    const double target = one;
    const int seeded = (mode == MINN) ? (one < BIG) : (one > -BIG);
    if (seeded) {
      int iwin = n;
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,n),
        KOKKOS_LAMBDA(const int i, int &lidx) {
          if (l_include(i) && l_values(i) == target && i < lidx) lidx = i;
        },Kokkos::Min<int>(iwin));
      if (iwin < n) index = iwin;
    }
  }

  return one;

#endif
}

/* ----------------------------------------------------------------------
   count elements included in the reduction, summed across procs
   device version of ComputeReduce::count_included()
------------------------------------------------------------------------- */

bigint ComputeReduceKokkos::count_included_kokkos()
{
  // the surf ownership decomposition and masks are host-only

  if (flavor[0] == SURF) {
    sync_host_for_fallback(0);
    return ComputeReduce::count_included();
  }

  bigint ncount = 0;
  bigint ncountall = 0;

  if (!subsetID) {
    if (flavor[0] == PARTICLE) ncount = particle->nlocal;
    else ncount = grid->nlocal;

  } else {
    build_include(0);

    // an integer count is order independent, so no SPARTA_KOKKOS_EXACT
    //   special case is needed here

    bigint nsum = 0;
    if (nelements) {
      auto l_include = d_include;
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,nelements),
        KOKKOS_LAMBDA(const int i, bigint &lsum) {
          lsum += l_include(i);
        },nsum);
    }
    ncount = nsum;
  }

  MPI_Allreduce(&ncount,&ncountall,1,MPI_SPARTA_BIGINT,MPI_SUM,world);

  return ncountall;
}

/* ---------------------------------------------------------------------- */

void ComputeReduceKokkos::grow_scratch(int n)
{
  if (n <= maxelements) return;
  maxelements = n;
  MemKK::realloc_kokkos(d_values,"reduce/kk:values",maxelements);
  MemKK::realloc_kokkos(d_include,"reduce/kk:include",maxelements);
}

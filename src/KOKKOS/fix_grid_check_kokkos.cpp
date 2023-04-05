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

#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "fix_grid_check_kokkos.h"
#include "update.h"
#include "particle_kokkos.h"
#include "grid_kokkos.h"
#include "comm.h"
#include "error.h"
#include "sparta_masks.h"

#include <cstdlib>

using namespace SPARTA_NS;

enum{ERROR,WARNING,SILENT};
enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};   // same as Grid

FixGridCheckKokkos::FixGridCheckKokkos(SPARTA *sparta, int narg, char **arg) :
  FixGridCheck(sparta, narg, arg)
{
  kokkos_flag = 1;
  execution_space = Device;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

void FixGridCheckKokkos::end_of_step()
{
  if (update->ntimestep % nevery) return;

  auto particleKK = dynamic_cast<ParticleKokkos*>(particle);
  particleKK->k_particles.sync_device();
  auto d_particles = particleKK->k_particles.d_view;
  auto gridKK = dynamic_cast<GridKokkos*>(grid);
  gridKK->k_cells.sync_device();
  auto d_cells = gridKK->k_cells.d_view;
  gridKK->k_cinfo.sync_device();
  auto d_cinfo = gridKK->k_cinfo.d_view;
  gridKK->k_sinfo.sync_device();
  auto d_sinfo = gridKK->k_sinfo.d_view;
  int nglocal = grid->nlocal;
  int nlocal = particle->nlocal;

  enum {
    NO_PROBLEM             = 0,
    IS_IN_INVALID_CELL     = 1,
    IS_OUTSIDE_CELL        = 2,
    IS_IN_SPLIT_CELL       = 4,
    IS_IN_INTERIOR_CELL    = 8,
    IS_IN_ZERO_VOLUME_CELL =16,
  };
  Kokkos::View<int*> d_particle_problems("d_particles_problems", nlocal);

  // check if icell is a valid cell for owning particles
  // check for split cell is whether particle is inside parent cell

  int nflag = 0;
  Kokkos::parallel_reduce(nlocal, KOKKOS_LAMBDA(int i, int& local_nflag) {
    auto icell = d_particles(i).icell;

    // is icell a valid index
    if (icell < 0 || icell >= nglocal) {
      d_particle_problems(i) |= IS_IN_INVALID_CELL;
      local_nflag++;
    }

    // does particle coord match icell bounds
    double* lo = d_cells[icell].lo;
    double* hi = d_cells[icell].hi;
    double* x = d_particles[i].x;
    if (x[0] < lo[0] || x[0] > hi[0] ||
        x[1] < lo[1] || x[1] > hi[1] ||
        x[2] < lo[2] || x[2] > hi[2]) {
      d_particle_problems(i) |= IS_OUTSIDE_CELL;
      local_nflag++;
    }

    // error if icell is a split cell, since should be a sub cell
    if (d_cells[icell].nsplit > 1) {
      d_particle_problems(i) |= IS_IN_SPLIT_CELL;
      local_nflag++;
    }

    // error if icell is an interior cell, since particle is inside surfs
    if (d_cinfo[icell].type == INSIDE) {
      d_particle_problems(i) |= IS_IN_INTERIOR_CELL;
      local_nflag++;
    }

    // error if icell has zero volume, since collision attempt freq will blow up
    if (d_cinfo[icell].volume == 0.0) {
      d_particle_problems(i) |= IS_IN_ZERO_VOLUME_CELL;
      local_nflag++;
    }

    // This check not yet supported

    //// check if particle in a cell with surfs is outside the surfs
    //// for split cell, also verify particle is in correct sub cell
    //// expensive, so only do this check if requested
    //
    //if (!outside_check) continue;
    //if (cells[icell].nsurf == 0) continue;
    //
    //int splitcell,subcell,flag;
    //
    //if (cells[icell].nsplit <= 0) {
    //  splitcell = sinfo[cells[icell].isplit].icell;
    //  flag = grid->outside_surfs(splitcell,x,cut3d,cut2d);
    //} else flag = grid->outside_surfs(icell,x,cut3d,cut2d);
    //
    //if (!flag) {
    //  if (outflag == ERROR) {
    //    char str[128];
    //    sprintf(str,
    //            "Particle %d,%d on proc %d is inside surfs in cell "
    //            CELLINT_FORMAT " on timestep " BIGINT_FORMAT,
    //            i,particles[i].id,comm->me,cells[icell].id,
    //            update->ntimestep);
    //    error->one(FLERR,str);
    //  }
    //  nflag++;
    //}
    //
    //if (cells[icell].nsplit <= 0) {
    //  int subcell;
    //  if (dim == 2) subcell = update->split2d(splitcell,x);
    //  else subcell = update->split3d(splitcell,x,particles[i].id);
    //
    //  if (subcell != icell) {
    //    if (outflag == ERROR) {
    //      char str[128];
    //      sprintf(str,
    //              "Particle %d,%d on proc %d is in wrong sub cell %d not %d"
    //              " on timestep " BIGINT_FORMAT,
    //              i,particles[i].id,comm->me,icell,subcell,
    //              update->ntimestep);
    //      error->one(FLERR,str);
    //    }
    //    nflag++;
    //  }
    //}

  }, nflag);

  // -------------------------------------
  // done with all tests
  // warning message instead of error

  if (outflag == WARNING) {
    int all;
    MPI_Allreduce(&nflag,&all,1,MPI_INT,MPI_SUM,world);
    if (all && comm->me == 0) {
      char str[128];
      sprintf(str,"%d particles were in wrong cells on timestep "
              BIGINT_FORMAT,all,update->ntimestep);
      error->warning(FLERR,str);
    }
  }
  if (outflag == ERROR && nflag) {
    auto h_particle_problems = Kokkos::create_mirror_view(d_particle_problems);
    Kokkos::deep_copy(h_particle_problems, d_particle_problems);
    auto cells = grid->cells;
    auto particles = particle->particles;
    char str[128];
    for (int i = 0; i < nlocal; ++i) {
      auto icell = particles[i].icell;
      if (h_particle_problems(i) & IS_IN_INVALID_CELL) {
        sprintf(str,
                "Particle %d,%d on proc %d is in invalid cell " CELLINT_FORMAT
                " on timestep " BIGINT_FORMAT,
                i,particles[i].id,comm->me,cells[icell].id,update->ntimestep);
        error->one(FLERR,str);
      }
      if (h_particle_problems(i) & IS_OUTSIDE_CELL) {
        sprintf(str,
                "Particle %d,%d on proc %d is outside cell " CELLINT_FORMAT
                " on timestep " BIGINT_FORMAT,
                i,particles[i].id,comm->me,cells[icell].id,
                update->ntimestep);
        error->one(FLERR,str);
      }
      if (h_particle_problems(i) & IS_IN_SPLIT_CELL) {
        sprintf(str,
                "Particle %d,%d on proc %d is in split cell " CELLINT_FORMAT
                " on timestep " BIGINT_FORMAT,
                i,particles[i].id,comm->me,cells[icell].id,
                update->ntimestep);
        error->one(FLERR,str);
      }
      if (h_particle_problems(i) & IS_IN_INTERIOR_CELL) {
        sprintf(str,
                "Particle %d,%d on proc %d is in interior cell " CELLINT_FORMAT
                " on timestep " BIGINT_FORMAT,
                i,particles[i].id,comm->me,cells[icell].id,update->ntimestep);
        error->one(FLERR,str);
      }
      if (h_particle_problems(i) & IS_IN_ZERO_VOLUME_CELL) {
        sprintf(str,
                "Particle %d,%d on proc %d is in volume=0 cell " CELLINT_FORMAT
                " on timestep " BIGINT_FORMAT,
                i,particles[i].id,comm->me,cells[icell].id,update->ntimestep);
        error->one(FLERR,str);
      }
    }
  }

  ntotal += nflag;
}

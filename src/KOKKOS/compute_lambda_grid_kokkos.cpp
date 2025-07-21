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

/* ----------------------------------------------------------------------
   Contributing author: Alan Stagg (Sandia)
------------------------------------------------------------------------- */

#include "string.h"
#include "compute_lambda_grid_kokkos.h"
#include "update.h"
#include "grid_kokkos.h"
#include "domain.h"
#include "collide.h"
#include "modify.h"
#include "particle_kokkos.h"
#include "fix.h"
#include "compute.h"
#include "math_const.h"
#include "memory.h"
#include "memory_kokkos.h"
#include "error.h"
#include "kokkos.h"
#include "sparta_masks.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{NONE,COMPUTE,FIX};
enum{LAMBDA,TAU,KNALL,KNX,KNY,KNZ};

#define MAXOUTPUT 6
#define INVOKED_PER_GRID 16
#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

ComputeLambdaGridKokkos::ComputeLambdaGridKokkos(SPARTA *sparta, int narg, char **arg) :
  ComputeLambdaGrid(sparta, narg, arg)
{
  kokkos_flag = 1;

  auto k_numap = DAT::tdual_float_1d("lambda/grid:numap",nrho_values);
  auto k_umap = DAT::tdual_float_2d("lambda/grid:umap",nrho_values,tmax);
  auto k_uomap = DAT::tdual_float_2d("lambda/grid:uomap",nrho_values,tmax);

  for (int i = 0; i < nrho_values; i++) {
    k_numap.h_view(i) = numap[i];
    for (int j = 0; j < tmax; j++) {
      k_umap.h_view(i,j) = umap[i][j];
      k_uomap.h_view(i,j) = uomap[i][j];
    }
  }
  k_numap.modify_host();
  k_numap.sync_device();
  d_numap = k_numap.d_view;

  k_umap.modify_host();
  k_umap.sync_device();
  d_umap = k_umap.d_view;

  k_uomap.modify_host();
  k_uomap.sync_device();
  d_uomap = k_uomap.d_view;

  auto k_output_order = DAT::tdual_int_1d("lambda/grid:output_order",MAXOUTPUT);

  for (int i = 0; i < MAXOUTPUT; i++)
    k_output_order.h_view(i) = output_order[i];

  k_output_order.modify_host();
  k_output_order.sync_device();
  d_output_order = k_output_order.d_view;
}

/* ---------------------------------------------------------------------- */

ComputeLambdaGridKokkos::~ComputeLambdaGridKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_vector_grid,vector_grid);
  memoryKK->destroy_kokkos(k_array_grid,array_grid);
  vector_grid = NULL;
  array_grid = NULL;
}

/* ---------------------------------------------------------------------- */

void ComputeLambdaGridKokkos::compute_per_grid()
{
  if (sparta->kokkos->prewrap)
    ComputeLambdaGrid::compute_per_grid();
  else
    compute_per_grid_kokkos();
}

/* ---------------------------------------------------------------------- */

void ComputeLambdaGridKokkos::compute_per_grid_kokkos()
{
  nspecies = ntotal;

  Kokkos::deep_copy(d_nrho,0.0);
  Kokkos::deep_copy(d_lambdainv,0.0);
  Kokkos::deep_copy(d_tauinv,0.0);

  invoked_per_grid = update->ntimestep;

  if (tempwhich == FIX && update->ntimestep % ftemp->per_grid_freq)
    error->all(FLERR,"Compute lambda/grid fix not computed at compatible time");

  // grab nrho and temp values from compute or fix
  // invoke nrho and temp computes as needed

  auto l_nrho = d_nrho;

  for (int m = 0; m < nrho_values; m++) {
    const int n = value2index[m];
    const int j = nrhoindex[m];

    if (nrhowhich[m] == FIX && update->ntimestep % modify->fix[n]->per_grid_freq)
      error->all(FLERR,"Compute lambda/grid fix not computed at compatible time");

    if (nrhowhich[m] == COMPUTE) {
      Compute *compute = modify->compute[n];
      if (!compute->kokkos_flag)
        error->all(FLERR,"Cannot (yet) use non-Kokkos computes with compute lambda/grid/kk");
      KokkosBase* cKKBase = dynamic_cast<KokkosBase*>(compute);

      if (!(compute->invoked_flag & INVOKED_PER_GRID)) {
        cKKBase->compute_per_grid_kokkos();
        compute->invoked_flag |= INVOKED_PER_GRID;
      }

      // accumulate one or more compute values to umap columns of tally array
      // if compute does not post-process, access its vec/array grid directly
      // else access uomap columns in its ctally array

      if (post_process[m]) {
        const int ntally_col = numap[m];
        DAT::t_float_2d_lr l_ctally;
        auto l_umap = d_umap;
        auto l_uomap = d_uomap;
        cKKBase->query_tally_grid_kokkos(l_ctally);
        Kokkos::parallel_for(nglocal, SPARTA_LAMBDA(int i) {
          for (int itally = 0; itally < ntally_col; itally++) {
            const int k = l_umap(m,itally);
            const int kk = l_uomap(m,itally);
            l_nrho(i,k) = l_ctally(i,kk);
          }
        });

        const int k = umap[m][0];
        const int jm1 = j - 1;
        if (nrho_values == 1) {
            cKKBase->post_process_grid_kokkos(j,1,d_nrho,map[0],d_vector_grid);
            auto l_vector_grid = d_vector_grid;
            Kokkos::parallel_for(nglocal, SPARTA_LAMBDA(int i) {
              l_nrho(i,k) = l_vector_grid[i];
            });
        } else {
            cKKBase->post_process_grid_kokkos(j,1,d_nrho,map[m],Kokkos::subview(d_array_grid1,Kokkos::ALL(),m));
            auto l_array_grid1 = d_array_grid1;
            Kokkos::parallel_for(nglocal, SPARTA_LAMBDA(int i) {
              l_nrho(i,k) = l_array_grid1(i,jm1);;
            });
        }
      } else {
        const int k = umap[m][0];
        if (j == 0) {
          auto l_compute_vector = cKKBase->d_vector_grid;
          Kokkos::parallel_for(nglocal, SPARTA_LAMBDA(int i) {
            l_nrho(i,k) = l_compute_vector[i];
          });
        } else {
          const int jm1 = j - 1;
          auto l_compute_array = cKKBase->d_array_grid;
          Kokkos::parallel_for(nglocal, SPARTA_LAMBDA(int i) {
            l_nrho(i,k) = l_compute_array(i,jm1);
          });
        }
      }

    // access fix fields, guaranteed to be ready

    } else if (nrhowhich[m] == FIX) {
      Fix *fix = modify->fix[n];
      if (!fix->kokkos_flag)
        error->all(FLERR,"Cannot (yet) use non-Kokkos fixes with compute lambda/grid/kk");
      KokkosBase* fKKBase = dynamic_cast<KokkosBase*>(fix);

      const int k = umap[m][0];
      if (j == 0) {
        auto l_fix_vector = fKKBase->d_vector_grid;
        Kokkos::parallel_for(nglocal, SPARTA_LAMBDA(int i) {
          l_nrho(i,k) = l_fix_vector[i];
        });
      } else {
        int jm1 = j - 1;
        auto l_fix_array = fKKBase->d_array_grid;
        Kokkos::parallel_for(nglocal, SPARTA_LAMBDA(int i) {
          l_nrho(i,k) = l_fix_array(i,jm1);
        });
      }
    }
  }

  auto l_temp = d_temp;

  if (tempwhich == COMPUTE) {
    if (!ctemp->kokkos_flag)
      error->all(FLERR,"Cannot (yet) use non-Kokkos computes with compute lambda/grid/kk");
    KokkosBase* ctempKKBase = dynamic_cast<KokkosBase*>(ctemp);

    if (!(ctemp->invoked_flag & INVOKED_PER_GRID)) {
      ctempKKBase->compute_per_grid_kokkos();
      ctemp->invoked_flag |= INVOKED_PER_GRID;
    }

    if (ctemp->post_process_grid_flag)
      ctempKKBase->post_process_grid_kokkos(tempindex,1,DAT::t_float_2d_lr(),NULL,DAT::t_float_1d_strided());

    if (tempindex == 0 || ctemp->post_process_grid_flag)
      Kokkos::deep_copy(d_temp,ctempKKBase->d_vector_grid);
    else {
      const int index = tempindex-1;
      auto l_array = ctempKKBase->d_array_grid;
      Kokkos::parallel_for(nglocal, SPARTA_LAMBDA(int i) {
        l_temp[i] = l_array(i,index);
      });
    }

  } else if (tempwhich == FIX) {
    if (!ftemp->kokkos_flag)
      error->all(FLERR,"Cannot (yet) use non-Kokkos fixes with compute lambda/grid/kk");
    KokkosBase* ftempKKBase = dynamic_cast<KokkosBase*>(ftemp);

    if (tempindex == 0)
      Kokkos::deep_copy(d_temp,ftempKKBase->d_vector_grid);
    else {
      auto l_array = ftempKKBase->d_array_grid;
      const int index = tempindex-1;
      Kokkos::parallel_for(nglocal, SPARTA_LAMBDA(int i) {
        l_temp[i] = l_array(i,index);
      });
    }
  }

  boltz = update->boltz;

  ParticleKokkos* particle_kk = ((ParticleKokkos*)particle);
  particle_kk->sync(Device,SPECIES_MASK);
  d_species = particle_kk->k_species.d_view;

  CollideVSSKokkos* collide_kk = ((CollideVSSKokkos*)collide);
  d_params_const = collide_kk->d_params_const;

  // compute mean free path for each grid cell
  // formula from Bird, eq 4.77

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeLambdaGrid_ComputePerGrid>(0,nglocal),*this);
  copymode = 0;

  // calculate per-cell Knudsen number

  if (knanyflag) {

    GridKokkos* grid_kk = ((GridKokkos*)grid);
    grid_kk->sync(Device,CELL_MASK);
    auto l_cells = grid_kk->k_cells.d_view;
    auto l_lambda_grid = d_lambda_grid;
    auto l_vector_grid = d_vector_grid;
    auto l_array_grid = d_array_grid;
    auto l_output_order = d_output_order;
    const int dimension = domain->dimension;
    auto l_knxflag = knxflag;
    auto l_knyflag = knyflag;
    auto l_knzflag = knzflag;
    auto l_knallflag = knallflag;
    auto l_noutputs = noutputs;

    Kokkos::parallel_for(nglocal, SPARTA_LAMBDA(int i) {
      const double lambda = l_lambda_grid(i);
      double sizex,sizey,sizez,sizeall;

      if (l_knxflag || l_knallflag)
        sizex = (l_cells[i].hi[0] - l_cells[i].lo[0]);

      if (l_knyflag || l_knallflag)
        sizey = (l_cells[i].hi[1] - l_cells[i].lo[1]);

      if (l_knzflag || (l_knallflag && dimension > 2))
        sizez = (l_cells[i].hi[2] - l_cells[i].lo[2]);

      if (l_knallflag) {
        sizeall = sizex + sizey;

        if (dimension == 2) sizeall *= 0.5;
        else {
          sizeall += sizez;
          sizeall /= 3.0;
        }
        if (l_noutputs == 1) l_vector_grid[i] = lambda / sizeall;
        else l_array_grid(i,l_output_order[KNALL]) = lambda / sizeall;
      }

      if (l_knxflag) {
        if (l_noutputs == 1) l_vector_grid[i] = lambda / sizex;
        else l_array_grid(i,l_output_order[KNX]) = lambda / sizex;
      }

      if (l_knyflag) {
        if (l_noutputs == 1) l_vector_grid[i] = lambda / sizey;
        l_array_grid(i,l_output_order[KNY]) = lambda / sizey;
      }

      if (l_knzflag) {
        if (l_noutputs == 1) l_vector_grid[i] = lambda / sizez;
        l_array_grid(i,l_output_order[KNZ]) = lambda / sizez;
      }
    });
  }

  if (noutputs == 1) {
    k_vector_grid.modify_device();
    k_vector_grid.sync_host();
  } else {
    k_array_grid.modify_device();
    k_array_grid.sync_host();
  }
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void ComputeLambdaGridKokkos::operator()(TagComputeLambdaGrid_ComputePerGrid, const int &i) const {
  double nrhosum,lambda,tau;
  nrhosum = lambda = tau = 0.0;
  for (int j = 0; j < nspecies; j++) {
    nrhosum += d_nrho(i,j);
    for (int k = 0; k < nspecies; k++) {
      const double dref = d_params_const(j,k).diam;
      const double tref = d_params_const(j,k).tref;
      const double omega = d_params_const(j,k).omega;
      const double mj = d_species[j].mass;
      const double mk = d_species[k].mass;
      const double mr = mj * mk / (mj + mk);

      if (tempwhich == NONE || d_temp[i] == 0.0) {
        if (lambdaflag)
          d_lambdainv(i,j) += (MY_PI * sqrt (1+mj/mk) * pow(dref,2.0) * d_nrho(i,k));
        if (tauflag)
          d_tauinv(i,j) += (2.0 * pow(dref,2.0) * d_nrho(i,k) * sqrt (2.0 * MY_PI * boltz * tref / mr));
      } else {
        if (lambdaflag)
          d_lambdainv(i,j) += (MY_PI * sqrt (1+mj/mk) * pow(dref,2.0) * d_nrho(i,k) * pow(tref/d_temp[i],omega-0.5));

        if (tauflag)
          d_tauinv(i,j) += (2.0 * pow(dref,2.0) * d_nrho(i,k) * sqrt (2.0 * MY_PI * boltz * tref / mr) * pow(d_temp[i]/tref,1.0-omega));
      }
    }
  }

  for (int j = 0; j < nspecies; j++) {
    if (lambdaflag && d_lambdainv(i,j) > 1e-30) lambda += d_nrho(i,j) / (nrhosum * d_lambdainv(i,j));
    if (tauflag && d_tauinv(i,j) > 1e-30) tau += d_nrho(i,j) / (nrhosum * d_tauinv(i,j));
  }

  // store per-grid lambda for possible later use in Knudsen numbers
  // lambdaflag may be set with no output of lambda

  if (lambdaflag) {
    if (lambda == 0.0) lambda = BIG;
    d_lambda_grid[i] = lambda;
    if (d_output_order[LAMBDA] >= 0) {
      if (noutputs == 1) d_vector_grid[i] = lambda;
      else d_array_grid(i,d_output_order[LAMBDA]) = lambda;
    }
  }

  if (tauflag) {
    if (tau == 0.0) tau = BIG;
    if (noutputs == 1) d_vector_grid[i] = tau;
    else d_array_grid(i,d_output_order[TAU]) = tau;
  }
}

/* ----------------------------------------------------------------------
   reallocate arrays if nglocal has changed
   called by init() and load balancer
------------------------------------------------------------------------- */

void ComputeLambdaGridKokkos::reallocate()
{
  if (grid->nlocal == nglocal) return;

  nglocal = grid->nlocal;

  memoryKK->destroy_kokkos(k_vector_grid,vector_grid);
  memoryKK->create_kokkos(k_vector_grid,vector_grid,nglocal,"lambda/grid:vector_grid");
  d_vector_grid = k_vector_grid.d_view;

  if (nrho_values > 1)
    d_array_grid1 = decltype(d_array_grid1)("lambda/grid:array_grid1",nglocal,nrho_values);

  if (noutputs > 1) {
    memoryKK->destroy_kokkos(k_array_grid,array_grid);
    memoryKK->create_kokkos(k_array_grid,array_grid,nglocal,noutputs,"lambda/grid:array_grid");
    d_array_grid = k_array_grid.d_view;
  }

  d_lambda_grid = decltype(d_lambda_grid)("lambda/grid:lambda_grid",nglocal);
  d_lambdainv = decltype(d_lambdainv)("lambda/grid:lambdainv",nglocal,ntotal);
  d_tauinv = decltype(d_tauinv)("lambda/grid:tauinv",nglocal,ntotal);
  d_nrho = decltype(d_nrho)("lambda/grid:nrho",nglocal,ntotal);

  memory->destroy(nrho);
  memory->create(nrho,nglocal,ntotal,"lambda/grid:nrho");

  if (tempwhich != NONE)
    d_temp = decltype(d_temp)("lambda/grid:temp",nglocal);
}

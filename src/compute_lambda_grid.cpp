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

#include "string.h"
#include "compute_lambda_grid.h"
#include "update.h"
#include "grid.h"
#include "domain.h"
#include "collide.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "math_const.h"
#include "memory.h"
#include "input.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{NONE,COMPUTE,FIX};
enum{LAMBDA,TAU,KNALL,KNX,KNY,KNZ};

#define INVOKED_PER_GRID 16
#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

ComputeLambdaGrid::ComputeLambdaGrid(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg < 5 || narg > 10)
    error->all(FLERR,"Illegal compute lambda/grid command");

  // parse three required input fields
  // customize a new keyword by adding to if statement

  id_temp = NULL;

  if (strncmp(arg[3],"c_",2) == 0 || strncmp(arg[3],"f_",2) == 0) {
    int n = strlen(arg[3]);
    id_temp = new char[n];
    strcpy(id_temp,&arg[3][2]);

    char *ptr = strchr(id_temp,'[');
    if (ptr) {
      if (id_temp[strlen(id_temp)-1] != ']')
        error->all(FLERR,"Invalid temp in compute lambda/grid command");
      tempindex = atoi(ptr+1);
      *ptr = '\0';
    } else tempindex = 0;

    if (strncmp(arg[3],"c_",2) == 0) tempwhich = COMPUTE;
    else tempwhich = FIX;

    if (tempwhich == COMPUTE) {
      int n = modify->find_compute(id_temp);
      if (n < 0)
        error->all(FLERR,"Could not find compute lambda/grid compute ID");
      if (modify->compute[n]->per_grid_flag == 0)
        error->all(FLERR,"Compute lambda/grid compute does not "
                   "compute per-grid info");
      if (tempindex == 0 && modify->compute[n]->size_per_grid_cols > 0)
        error->all(FLERR,
                   "Compute lambda/grid compute does not "
                   "compute per-grid vector");
      if (tempindex > 0 && modify->compute[n]->size_per_grid_cols == 0)
        error->all(FLERR,
                   "Compute lambda/grid compute does not "
                   "compute per-grid array");
      if (tempindex > 0 && tempindex > modify->compute[n]->size_per_grid_cols)
        error->all(FLERR,"Compute lambda compute vector is "
                   "accessed out-of-range");
    } else {
      int n = modify->find_fix(id_temp);
      if (n < 0) error->all(FLERR,"Could not find compute lambda/grid fix ID");
      if (modify->fix[n]->per_grid_flag == 0)
        error->all(FLERR,"Compute lambda/grid fix does not "
                   "compute per-grid info");
      if (tempindex == 0 && modify->fix[n]->size_per_grid_cols > 0)
        error->all(FLERR,"Compute lambda/grid fix does not "
                   "compute per-grid vector");
      if (tempindex > 0 && modify->fix[n]->size_per_grid_cols == 0)
        error->all(FLERR,"Compute lambda/grid fix does not "
                   "compute per-grid array");
      if (tempindex > 0 && tempindex > modify->fix[n]->size_per_grid_cols)
        error->all(FLERR,"Compute lambda/grid fix array is "
                   "accessed out-of-range");
    }
  } else if (strcmp(arg[3],"NULL") == 0) tempwhich = NONE;
  else error->all(FLERR,"Illegal compute lambda/grid command");

  lambdaflag = 0;
  tauflag = 0;
  knallflag = 0;
  knxflag = 0;
  knyflag = 0;
  knzflag = 0;
  knanyflag = 0;

  noutputs = narg - 4;

  int maxoutput = 6;
  output_order = new int[maxoutput];

  for (int i = 0; i < maxoutput; i++)
   output_order[i] = -1;

  int ioutput = 0;
  int iarg = 4;
  int dupflag = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"lambda") == 0) {
      if (output_order[LAMBDA] != -1) dupflag = 1;
      output_order[LAMBDA] = ioutput;
      lambdaflag = 1;
    } else if (strcmp(arg[iarg],"tau") == 0) {
      if (output_order[TAU] != -1) dupflag = 1;
      output_order[TAU] = ioutput;
      tauflag = 1;
    } else if (strcmp(arg[iarg],"knall") == 0) {
      if (output_order[KNALL] != -1) dupflag = 1;
      output_order[KNALL] = ioutput;
      knallflag = 1;
    } else if (strcmp(arg[iarg],"knx") == 0) {
      if (output_order[KNX] != -1) dupflag = 1;
      output_order[KNX] = ioutput;
      knxflag = 1;
    } else if (strcmp(arg[iarg],"kny") == 0) {
      if (output_order[KNY] != -1) dupflag = 1;
      output_order[KNY] = ioutput;
      knyflag = 1;
    } else if (strcmp(arg[iarg],"knz") == 0) {
      if (output_order[KNZ] != -1) dupflag = 1;
      output_order[KNZ] = ioutput;
      knzflag = 1;
    } else error->all(FLERR,"Illegal compute lambda/grid command");

    ioutput++;
    iarg++;
  }

  if (dupflag)
    error->all(FLERR,"Duplicated output in compute lambda/grid");

  if (knzflag && domain->dimension == 2)
    error->all(FLERR,"Cannot use compute lambda/grid knz for 2d simulation");

  // expand args if any have wildcard character "*"
  // this can reset nvalues

  int expand = 0;
  char **earg;
  nvalues = input->expand_args(1,&arg[2],1,earg);

  if (earg != &arg[2]) expand = 1;
  arg = earg;

  // parse values

  nrhowhich = new int[nvalues];
  nrhoindex = new int[nvalues];
  value2index = new int[nvalues];
  post_process = new int[nvalues];
  ids_nrho = new char*[nvalues];

  for (int i = 0; i < nvalues; i++) {
    if (arg[i][0] == 'c') nrhowhich[i] = COMPUTE;
    else if (arg[i][0] == 'f') nrhowhich[i] = FIX;
    else error->all(FLERR,"Illegal compute lambda/grid command");

    int n = strlen(arg[i]);
    char *suffix = new char[n];
    strcpy(suffix,&arg[i][2]);

    char *ptr = strchr(suffix,'[');
    if (ptr) {
      if (suffix[strlen(suffix)-1] != ']')
        error->all(FLERR,"Illegal compute lambda/grid command");
      nrhoindex[i] = atoi(ptr+1);
      *ptr = '\0';
    } else nrhoindex[i] = 0;

    n = strlen(suffix) + 1;
    ids_nrho[i] = new char[n];
    strcpy(ids_nrho[i],suffix);
    delete [] suffix;

    post_process[i] = 0;
    if (nrhowhich[i] == COMPUTE) {
      int icompute = modify->find_compute(ids_nrho[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for compute lambda/grid does not exist");
      post_process[i] =
        modify->compute[icompute]->post_process_grid_flag;
    }
  }

  // if wildcard expansion occurred, free earg memory from expand_args()
  // wait to do this until after file comment lines are printed

  if (expand) {
    for (int i = 0; i < nvalues; i++) delete [] earg[i];
    memory->sfree(earg);
  }

  // setup and error check

  for (int i = 0; i < nvalues; i++) {
    if (nrhowhich[i] == COMPUTE) {
      int icompute = modify->find_compute(ids_nrho[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for compute lambda/grid does not exist");
      if (modify->compute[icompute]->per_grid_flag == 0)
        error->all(FLERR,
                   "Compute lambda/grid compute does not calculate per-grid values");
      if (nrhoindex[i] == 0 &&
          modify->compute[icompute]->size_per_grid_cols != 0)
        error->all(FLERR,"Compute lambda/grid compute does not "
                   "calculate per-grid vector");
      if (nrhoindex[i] && modify->compute[icompute]->size_per_grid_cols == 0)
        error->all(FLERR,"Compute lambda/grid compute does not "
                   "calculate per-grid array");
      if (nrhoindex[i] &&
          nrhoindex[i] > modify->compute[icompute]->size_per_grid_cols)
        error->all(FLERR,"Compute lambda/grid compute array is accessed out-of-range");

    } else if (nrhowhich[i] == FIX) {
      int ifix = modify->find_fix(ids_nrho[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for compute lambda/grid does not exist");
      if (modify->fix[ifix]->per_grid_flag == 0)
        error->all(FLERR,"Compute lambda/grid fix does not calculate per-grid values");
      if (nrhoindex[i] == 0 && modify->fix[ifix]->size_per_grid_cols != 0)
        error->all(FLERR,
                   "Compute lambda/grid fix does not calculate per-grid vector");
      if (nrhoindex[i] && modify->fix[ifix]->size_per_grid_cols == 0)
        error->all(FLERR,
                   "Compute lambda/grid fix does not calculate per-grid array");
      if (nrhoindex[i] && nrhoindex[i] > modify->fix[ifix]->size_per_grid_cols)
        error->all(FLERR,"Compute lambda/grid fix array is accessed out-of-range");
    } else error->all(FLERR,"Illegal compute lambda/grid command");
  }

  // initialize data structures

  nparams = particle->nspecies;

  per_grid_flag = 1;
  if (noutputs > 1) size_per_grid_cols = noutputs;
  else size_per_grid_cols = 0;

  // knudsen number needs lambda

  knanyflag = (knallflag || knxflag || knyflag || knzflag);

  if (knanyflag && !lambdaflag) {
    lambdaflag = 1;
    output_order[LAMBDA] = 0;
  }

  nglocal = 0;
  vector_grid = NULL;
  array_grid = NULL;
  array_grid1 = NULL;
  nrho = NULL;
  temp = NULL;
  lambdainv = NULL;
  tauinv = NULL;

  // determine size of map/umap/uomap data structs and allocate them
  // tmax = max # of tally quantities for any value

  tmax = 1;
  for (int m = 0; m < nvalues; m++) {
    int n = -1;
    int j = nrhoindex[m];
    if (nrhowhich[m] != COMPUTE) continue;
    n = modify->find_compute(ids_nrho[m]);
    if (!modify->compute[n]->post_process_grid_flag) continue;
    double **array;
    int *cmap;
    int ncount = modify->compute[n]->query_tally_grid(j,array,cmap);
    tmax = MAX(tmax,ncount);
  }

  nmap = new int[nvalues];
  memory->create(map,nvalues,tmax,"compute lambda/grid:map");
  numap = new int[nvalues];
  memory->create(umap,nvalues,tmax,"compute lambda/grid:umap");
  memory->create(uomap,nvalues,tmax,"compute lambda/grid:uomap");

  // setup nmap/map and numap/umap/uomap data structs for all values
  // ntotal = total # of unique tally quantities = columns in tally array

  ntotal = 0;
  for (int m = 0; m < nvalues; m++) {
    int n = -1;
    int j = nrhoindex[m];

    int pflag = 0;
    if (nrhowhich[m] == COMPUTE) {
      n = modify->find_compute(ids_nrho[m]);
      if (modify->compute[n]->post_process_grid_flag) pflag = 1;
    }

    // if not a compute that post-processes,
    // add single new tally to nmap/map and numap/umap

    if (!pflag) {
      nmap[m] = 1;
      map[m][0] = ntotal;
      numap[m] = 1;
      umap[m][0] = ntotal;
      ntotal++;

    // else add all compute tallies to nmap/map
    // and only unique compute tallies to numap/umap/uomap

    } else {
      double **array;
      int *cmap;
      int ncount = modify->compute[n]->query_tally_grid(j,array,cmap);
      nmap[m] = numap[m] = 0;
      for (int i = 0; i < ncount; i++) {

        // set ucol = -1 if first time this compute quantity is tallied
        // else set to tally column that already tallies it

        int col = -1;
        for (int mm = 0; mm <= m; mm++) {
          if (nrhowhich[mm] != COMPUTE) continue;
          int nn = modify->find_compute(ids_nrho[mm]);
          if (!modify->compute[nn]->post_process_grid_flag) continue;
          if (n != nn) continue;  // not same compute
          for (int kk = 0; kk < numap[mm]; kk++)
            if (cmap[i] == uomap[mm][kk]) col = umap[mm][kk];
        }

        // if this quantity already tallied, just point nmap/map to it
        // else add to nmap/map and numap/umap/uomap

        if (col >= 0) {
          map[m][nmap[m]] = col;
          nmap[m]++;
        } else {
          map[m][nmap[m]] = ntotal;
          nmap[m]++;
          umap[m][numap[m]] = ntotal;
          uomap[m][numap[m]] = cmap[i];
          numap[m]++;
          ntotal++;
        }
      }
    }
  }

  if (nparams != ntotal)
      error->all(FLERR,"Number of species does not match size of compute vector or array");
  if (ntotal == 0)
    error->all(FLERR,"Cannot use compute lambda/grid command with no species defined");
}

/* ---------------------------------------------------------------------- */

ComputeLambdaGrid::~ComputeLambdaGrid()
{
  if (copymode) return;

  delete [] nrhowhich;
  delete [] nrhoindex;
  delete [] value2index;
  delete [] post_process;
  for (int i = 0; i < nvalues; i++) delete [] ids_nrho[i];
  delete [] ids_nrho;

  delete [] nmap;
  memory->destroy(map);
  delete [] numap;
  memory->destroy(umap);
  memory->destroy(uomap);

  if (nvalues == 1) memory->destroy(vector_grid);
  else memory->destroy(array_grid1);
  memory->destroy(array_grid);

  memory->destroy(nrho);
  memory->destroy(lambdainv);
  memory->destroy(tauinv);
  delete [] id_temp;
  memory->destroy(temp);
}

/* ---------------------------------------------------------------------- */

void ComputeLambdaGrid::init()
{
  reallocate();

  // initially read-in per-species params must match current species list

  if (nparams != particle->nspecies)
    error->all(FLERR,"VSS parameters do not match current species");
  if (ntotal != particle->nspecies)
    error->all(FLERR,"Compute array size does not match current species");

  // setup computes and fixes

  for (int m = 0; m < nvalues; m++) {
    if (nrhowhich[m] == COMPUTE) {
      int icompute = modify->find_compute(ids_nrho[m]);
      if (icompute < 0)
        error->all(FLERR,"Could not find compute lambda/grid compute ID");
      value2index[m] = icompute;

    } else if (nrhowhich[m] == FIX) {
      int ifix = modify->find_fix(ids_nrho[m]);
      if (ifix < 0)
        error->all(FLERR,"Could not find compute lambda/grid fix ID");
      value2index[m] = ifix;

    } else value2index[m] = -1;
  }

  if (tempwhich == COMPUTE) {
    int icompute = modify->find_compute(id_temp);
    if (icompute < 0)
      error->all(FLERR,"Could not find compute lambda/grid compute ID");
    ctemp = modify->compute[icompute];
  } else if (tempwhich == FIX) {
    int ifix = modify->find_fix(id_temp);
    if (ifix < 0)
      error->all(FLERR,"Could not find compute lambda/grid fix ID");
    ftemp = modify->fix[ifix];
  }

  // setup collision parameter

  if (!collide)
    error->all(FLERR,"Compute lambda/grid requires a "
               "collision style be defined");
}

/* ---------------------------------------------------------------------- */

void ComputeLambdaGrid::compute_per_grid()
{
  int i,j,k,n,itally;
  int ntally_col,kk;
  int *itmp;
  double **ctally;

  for (i = 0; i < nglocal; i++) {
    for (j = 0; j < ntotal; j++) {
        nrho[i][j] = 0.0;
        lambdainv[i][j] = 0.0;
        tauinv[i][j] = 0.0;
    }
  }

  invoked_per_grid = update->ntimestep;

  if (tempwhich == FIX && update->ntimestep % ftemp->per_grid_freq)
    error->all(FLERR,"Compute lambda/grid fix not computed at compatible time");

  // grab nrho and temp values from compute or fix
  // invoke nrho and temp computes as needed

  for (int m = 0; m < nvalues; m++) {
    n = value2index[m];
    j = nrhoindex[m];

    if (nrhowhich[m] == FIX && update->ntimestep % modify->fix[n]->per_grid_freq)
      error->all(FLERR,"Compute lambda/grid fix not computed at compatible time");

    if (nrhowhich[m] == COMPUTE) {
      Compute *compute = modify->compute[n];
      if (!(compute->invoked_flag & INVOKED_PER_GRID)) {
        compute->compute_per_grid();
        compute->invoked_flag |= INVOKED_PER_GRID;
      }

      // accumulate one or more compute values to umap columns of tally array
      // if compute does not post-process, access its vec/array grid directly
      // else access uomap columns in its ctally array

      if (post_process[m]) {
        ntally_col = numap[m];
        compute->query_tally_grid(j,ctally,itmp);
        for (i = 0; i < nglocal; i++) {
          for (itally = 0; itally < ntally_col; itally++) {
            k = umap[m][itally];
            kk = uomap[m][itally];
            nrho[i][k] = ctally[i][kk];
          }
        }

        k = umap[m][0];
        int jm1 = j - 1;
        if (nvalues == 1) {
            compute->post_process_grid(j,1,nrho,map[0],vector_grid,1);
            for (i = 0; i < nglocal; i++) nrho[i][k] = vector_grid[i];
        } else {
            compute->post_process_grid(j,1,nrho,map[m],&array_grid1[0][m],nvalues);
            for (i = 0; i < nglocal; i++) nrho[i][k] = array_grid1[i][jm1];
        }
      } else {
        k = umap[m][0];
        if (j == 0) {
          double *compute_vector = compute->vector_grid;
          for (i = 0; i < nglocal; i++)
            nrho[i][k] = compute_vector[i];
        } else {
          int jm1 = j - 1;
          double **compute_array = compute->array_grid;
          for (i = 0; i < nglocal; i++)
            nrho[i][k] = compute_array[i][jm1];
        }
      }

    // access fix fields, guaranteed to be ready

    } else if (nrhowhich[m] == FIX) {
      k = umap[m][0];
      if (j == 0) {
        double *fix_vector = modify->fix[n]->vector_grid;
        for (i = 0; i < nglocal; i++)
          nrho[i][k] = fix_vector[i];
      } else {
        int jm1 = j - 1;
        double **fix_array = modify->fix[n]->array_grid;
        for (i = 0; i < nglocal; i++) {
          nrho[i][k] = fix_array[i][jm1];
        }
      }
    }
  }

  if (tempwhich == COMPUTE) {
    if (!(ctemp->invoked_flag & INVOKED_PER_GRID)) {
      ctemp->compute_per_grid();
      ctemp->invoked_flag |= INVOKED_PER_GRID;
    }

    if (ctemp->post_process_grid_flag)
      ctemp->post_process_grid(tempindex,1,NULL,NULL,NULL,1);

    if (tempindex == 0 || ctemp->post_process_grid_flag)
      memcpy(temp,ctemp->vector_grid,nglocal*sizeof(double));
    else {
      int index = tempindex-1;
      double **array = ctemp->array_grid;
      for (int i = 0; i < nglocal; i++)
        temp[i] = array[i][index];
    }

  } else if (tempwhich == FIX) {
    if (tempindex == 0)
      memcpy(temp,ftemp->vector_grid,nglocal*sizeof(double));
    else {
      double **array = ftemp->array_grid;
      int index = tempindex-1;
      for (int i = 0; i < nglocal; i++)
        temp[i] = array[i][index];
    }
  }

  // compute mean free path for each grid cell
  // formula from Bird, eq 4.77

  double nrhosum,lambda,tau;

  for (int i = 0; i < nglocal; i++) {
    nrhosum = lambda = tau = 0.0;
    for (int j = 0; j < ntotal; j++) {
      nrhosum += nrho[i][j];
      for (int k = 0; k < ntotal; k++) {
        dref = collide->extract(j,k,"diam");
        tref = collide->extract(j,k,"tref");
        omega = collide->extract(j,k,"omega");
        mj = particle->species[j].mass;
        mk = particle->species[k].mass;
        mr = mj * mk / (mj + mk);
        if (tempwhich == NONE || temp[i] == 0.0) {
          if (lambdaflag)
            lambdainv[i][j] += (MY_PI * sqrt (1+mj/mk) * pow(dref,2.0) * nrho[i][k]);
          if (tauflag)
            tauinv[i][j] += (2.0 * pow(dref,2.0) * nrho[i][k] * sqrt (2.0 * MY_PI * update->boltz * tref / mr));
        } else {
          if (lambdaflag)
            lambdainv[i][j] += (MY_PI * sqrt (1+mj/mk) * pow(dref,2.0) * nrho[i][k] * pow(tref/temp[i],omega-0.5));

          if (tauflag)
            tauinv[i][j] += (2.0 * pow(dref,2.0) * nrho[i][k] * sqrt (2.0 * MY_PI * update->boltz * tref / mr) * pow(temp[i]/tref,1.0-omega));
        }
      }
    }

    for (int j = 0; j < ntotal; j++) {
      if (lambdaflag && lambdainv[i][j] > 1e-30) lambda += nrho[i][j] / (nrhosum * lambdainv[i][j]);
      if (tauflag && tauinv[i][j] > 1e-30) tau += nrho[i][j] / (nrhosum * tauinv[i][j]);
    }

    if (lambdaflag) {
      if (lambda == 0.0) lambda  = BIG;
      if (noutputs == 1 && !knanyflag) vector_grid[i] = lambda;
      else array_grid[i][output_order[LAMBDA]] = lambda;
    }

    if (tauflag) {
      if (tau == 0.0) tau = BIG;
      if (noutputs == 1) vector_grid[i] = tau;
      else array_grid[i][output_order[TAU]] = tau;
    }
  }

  // calculate per-cell Knudsen number

  if (!knanyflag) return;

  Grid::ChildCell *cells = grid->cells;
  int dimension = domain->dimension;
  double sizeall,sizex,sizey,sizez;

  for (int i = 0; i < nglocal; i++) {
    double lambda = array_grid[i][output_order[LAMBDA]];

    if (knxflag || knallflag)
      sizex = (cells[i].hi[0] - cells[i].lo[0]);

    if (knyflag || knallflag)
      sizey = (cells[i].hi[1] - cells[i].lo[1]);

    if (knzflag || (knallflag && dimension > 2))
      sizez = (cells[i].hi[2] - cells[i].lo[2]);

    if (knallflag) {
      sizeall = sizex + sizey;

      if (dimension == 2) sizeall *= 0.5;
      else {
        sizeall += sizez;
        sizeall /= 3.0;
      }
      if (noutputs == 1) vector_grid[i] = lambda / sizeall;
      else array_grid[i][output_order[KNALL]] = lambda / sizeall;
    }

    if (knxflag) {
      if (noutputs == 1) vector_grid[i] = lambda / sizex;
      else array_grid[i][output_order[KNX]] = lambda / sizex;
    }

    if (knyflag) {
      if (noutputs == 1) vector_grid[i] = lambda / sizey;
      array_grid[i][output_order[KNY]] = lambda / sizey;
    }

    if (knzflag) {
      if (noutputs == 1) vector_grid[i] = lambda / sizez;
      array_grid[i][output_order[KNZ]] = lambda / sizez;
    }
  }
}

/* ----------------------------------------------------------------------
   reallocate arrays if nglocal has changed
   called by init() and load balancer
------------------------------------------------------------------------- */

void ComputeLambdaGrid::reallocate()
{
  if (grid->nlocal == nglocal) return;

  nglocal = grid->nlocal;

  memory->destroy(vector_grid);
  memory->create(vector_grid,nglocal,"lambda/grid:vector_grid");
  for (int i = 0; i < nglocal; i++) vector_grid[i] = 0.0;

  if (nvalues > 1) {
    memory->destroy(array_grid1);
    memory->create(array_grid1,nglocal,nvalues,"lambda/grid:array_grid1");
    for (int i = 0; i < nglocal; i++)
      for (int m = 0; m < nvalues; m++) array_grid1[i][m] = 0.0;
  }

  if (noutputs > 1 || knanyflag) {
    memory->destroy(array_grid);
    memory->create(array_grid,nglocal,noutputs,"lambda/grid:array_grid");
    for (int i = 0; i < nglocal; i++)
      for (int m = 0; m < noutputs; m++) array_grid[i][m] = 0.0;
  }

  memory->destroy(lambdainv);
  memory->create(lambdainv,nglocal,ntotal,"lambda/grid:lambdainv");
  memory->destroy(tauinv);
  memory->create(tauinv,nglocal,ntotal,"lambda/grid:tauinv");

  memory->destroy(nrho);
  memory->create(nrho,nglocal,ntotal,"lambda/grid:nrho");

  if (tempwhich != NONE) {
    memory->destroy(temp);
    memory->create(temp,nglocal,"lambda/grid:temp");
  }
}

/* ----------------------------------------------------------------------
   memory usage of local grid-based array
------------------------------------------------------------------------- */

bigint ComputeLambdaGrid::memory_usage()
{
  bigint bytes;
  bytes = nglocal * sizeof(double);                            // vector_grid
  if (nvalues > 1)
    bytes = nglocal * nvalues * sizeof(double);                // array_grid1
  bytes += nglocal * noutputs * sizeof(double);                // array_grid
  bytes += 2 * nglocal * ntotal * sizeof(double);              // lambdainv + tauinv
  bytes += nglocal * ntotal * sizeof(double);                  // nrho
  if (tempwhich != NONE) bytes += nglocal * sizeof(double);    // temp
  return bytes;
}

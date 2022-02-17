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
#include "fix_field_forcing_grid.h"
#include "update.h"
#include "domain.h"
#include "grid.h"
#include "input.h"
#include "variable.h"
#include "random_knuth.h"
#include "random_mars.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixFieldForcingGrid::FixFieldForcingGrid(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  // require a 3d periodic box

  if (domain->dimension != 3) 
    error->all(FLERR,"Fix field/forcing/grid requires 3d simulation");

  int periodic[3];
  int pflag = domain->periodic(periodic);
  if (!domain->periodic(periodic))
    error->all(FLERR,"Fix field/forcing/grid requires fully periodic box");

  // parse args

  if (narg < 5) error->all(FLERR,"Illegal fix field/forcing/grid command");

  epsilon = atof(arg[2]);
  tau = atof(arg[3]);
  nkvec = atoi(arg[4]);

  if (nkvec <= 0) error->all(FLERR,"Fix field/forcing/grid requires Kspace vecs");
  if (narg-5 != 3*nkvec) 
    error->all(FLERR,"Illegal fix field/forcing/grid command");

  memory->create(kvecs,nkvec,3,"field/forcing/grid:kvecs");

  int iarg = 5;
  for (int i = 0; i < nkvec; i++) {
    kvecs[i][0] = atof(arg[iarg]);
    kvecs[i][1] = atof(arg[iarg+1]);
    kvecs[i][2] = atof(arg[iarg+2]);
    iarg += 3;
  }

  // fix settings

  per_grid_flag = 1;
  size_per_grid_cols = 3;
  per_grid_freq = 1;
  per_grid_field = 1;
  field_active[0] = field_active[1] = field_active[2] = 1;

  // random number generator, same on every proc
  
  random = new RanKnuth(update->ranmaster->uniform());

  // one-time K vector initialization

  memory->create(prefacN_r,nkvec,3,"field/forcing/grid:prefacN_r");
  memory->create(prefacN_i,nkvec,3,"field/forcing/grid:prefacN_i");
  memory->create(fr_prev,nkvec,3,"field/forcing/grid:fr_prev");
  memory->create(fi_prev,nkvec,3,"field/forcing/grid:fi_prev");

  // per-grid memory initialization

  maxgrid = 0;
  array_grid = NULL;
}

/* ---------------------------------------------------------------------- */

FixFieldForcingGrid::~FixFieldForcingGrid()
{  
  delete random;
  memory->destroy(kvecs);
  memory->destroy(prefacN_r);
  memory->destroy(prefacN_i);
  memory->destroy(fr_prev);
  memory->destroy(fi_prev);
  memory->destroy(array_grid);
}

/* ---------------------------------------------------------------------- */

int FixFieldForcingGrid::setmask()
{
  int mask = 0;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixFieldForcingGrid::init()
{
  // set initial grid values to zero in case dump is performed at step 0
  // NOTE: maybe this should only be done one-time

  if (grid->nlocal > maxgrid) {
    maxgrid = grid->maxlocal;
    memory->destroy(array_grid);
    memory->create(array_grid,maxgrid,3,"field/forcing:array_grid");
  }

  bigint nbytes = (bigint) grid->nlocal * 3;
  memset(&array_grid[0][0],0,nbytes*sizeof(double));

  // initialization of previous field values
  // NOTE: maybe this needs adjusting, maybe should only be done one-time

  laststep = update->ntimestep;

  for (int ik = 0; ik < nkvec; ik++) {
    fr_prev[ik][0] = fr_prev[ik][1] = fr_prev[ik][2] = 0.0;
    fi_prev[ik][0] = fi_prev[ik][1] = fi_prev[ik][2] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

void FixFieldForcingGrid::compute_field()
{
  // reallocate array_grid if necessary

  if (grid->nlocal > maxgrid) {
    maxgrid = grid->maxlocal;
    memory->destroy(array_grid);
    memory->create(array_grid,maxgrid,3,"field/forcing:array_grid");
  }

  // loop over kvecs to compute forcing coefficients prefacN_ri

  double dt = (update->ntimestep - laststep) * update->dt;
  laststep = update->ntimestep;

  for (int ik = 0; ik < nkvec; ik++) {
    double kx = kvecs[ik][0];
    double ky = kvecs[ik][1];
    double kz = kvecs[ik][2];
    double ksq = kx*kx + ky*ky + kz*kz;
    double rnxr = random->gaussian();
    double rnyr = random->gaussian();
    double rnzr = random->gaussian();
    double rnxi = random->gaussian();
    double rnyi = random->gaussian();
    double rnzi = random->gaussian();
      
    double fxr = (1.0-dt/tau) * fr_prev[ik][0] + sqrt(2.0*epsilon*dt)/tau * rnxr;
    double fyr = (1.0-dt/tau) * fr_prev[ik][1] + sqrt(2.0*epsilon*dt)/tau * rnyr;
    double fzr = (1.0-dt/tau) * fr_prev[ik][2] + sqrt(2.0*epsilon*dt)/tau * rnzr;
    
    double fxi = (1.0-dt/tau) * fi_prev[ik][0] + sqrt(2.0*epsilon*dt)/tau * rnxi;
    double fyi = (1.0-dt/tau) * fi_prev[ik][1] + sqrt(2.0*epsilon*dt)/tau * rnyi;
    double fzi = (1.0-dt/tau) * fi_prev[ik][2] + sqrt(2.0*epsilon*dt)/tau * rnzi;
    
    double frdotk = fxr*kx + fyr*ky + fzr*kz;
    double fidotk = fxi*kx + fyi*ky + fzi*kz;
    
    // solenoidal projection of real part

    prefacN_r[ik][0] = fxr - kx*frdotk/ksq;
    prefacN_r[ik][1] = fyr - ky*frdotk/ksq;
    prefacN_r[ik][2] = fzr - kz*frdotk/ksq;
    
    // solenoidal projection of imaginary part

    prefacN_i[ik][0] = fxi - kx * fidotk / ksq;
    prefacN_i[ik][1] = fyi - ky * fidotk / ksq;
    prefacN_i[ik][2] = fzi - kz * fidotk / ksq;

    // persist f real/image values from one field invocation to the next

    fr_prev[ik][0] = fxr;
    fr_prev[ik][1] = fyr;
    fr_prev[ik][2] = fzr;
    fi_prev[ik][0] = fxi;
    fi_prev[ik][1] = fyi;
    fi_prev[ik][2] = fzi;
  }

  // loop over grid cells
  // array_grid = current field at each grid cell center point

  double *prd = domain->prd;
  double *boxlo = domain->boxlo;

  Grid::ChildCell *cells = grid->cells;
  int nglocal = grid->nlocal;

  double xc,yc,zc;
  double kxc,kyc,kzc;

  for (int icell = 0; icell < nglocal; icell++) {

    // center point of grid cell

    xc = 0.5 * (cells[icell].lo[0] + cells[icell].hi[0]);
    yc = 0.5 * (cells[icell].lo[1] + cells[icell].hi[1]);
    zc = 0.5 * (cells[icell].lo[2] + cells[icell].hi[2]);

    // k xyz c = convert center pt to (-PI,PI) within each periodic dim

    kxc = (xc - boxlo[0]) / prd[0] * MY_2PI - MY_PI;
    kyc = (yc - boxlo[1]) / prd[1] * MY_2PI - MY_PI;
    kzc = (zc - boxlo[2]) / prd[2] * MY_2PI - MY_PI;

    // loop over KSpace vectors

    array_grid[icell][0] = 0.0;
    array_grid[icell][1] = 0.0;
    array_grid[icell][2] = 0.0;

    for (int ik = 0; ik < nkvec; ik++) {
      double kdotx = kvecs[ik][0]*kxc + kvecs[ik][1]*kyc + kvecs[ik][2]*kzc;
      double c = cos(kdotx);
      double s = sin(kdotx);
	  
      array_grid[icell][0] += 2.0 * (prefacN_r[ik][0]*c - prefacN_i[ik][0]*s);
      array_grid[icell][1] += 2.0 * (prefacN_r[ik][1]*c - prefacN_i[ik][1]*s);
      array_grid[icell][2] += 2.0 * (prefacN_r[ik][2]*c - prefacN_i[ik][2]*s);
    }
  }
}

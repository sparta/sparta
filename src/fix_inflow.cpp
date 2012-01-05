/* ----------------------------------------------------------------------
   DSMC - Sandia parallel DSMC code
   www.sandia.gov/~sjplimp/dsmc.html
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2011) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level DSMC directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_inflow.h"
#include "update.h"
#include "particle.h"
#include "mixture.h"
#include "domain.h"
#include "grid.h"
#include "comm.h"
#include "random_mars.h"
#include "random_park.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace DSMC_NS;
using namespace MathConst;

enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};     // same as Domain
enum{PERIODIC,OUTFLOW,SPECULAR};            // same as Domain

/* ---------------------------------------------------------------------- */

FixInflow::FixInflow(DSMC *dsmc, int narg, char **arg) :
  Fix(dsmc, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix inflow command");

  nevery = atoi(arg[2]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix inflow command");

  imix = particle->find_mixture(arg[3]);
  if (imix < 0) error->all(FLERR,"Fix inflow mixture ID does not exist");

  // extract list of faces

  faces[XLO] = faces[XHI] = faces[YLO] = faces[YHI] =
    faces[ZLO] = faces[ZHI] = 0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"all") == 0)
      faces[XLO] = faces[XHI] = faces[YLO] = faces[YHI] =
	faces[ZLO] = faces[ZHI] = 1;
    else if (strcmp(arg[iarg],"xlo") == 0) faces[XLO] = 1;
    else if (strcmp(arg[iarg],"xhi") == 0) faces[XHI] = 1;
    else if (strcmp(arg[iarg],"ylo") == 0) faces[YLO] = 1;
    else if (strcmp(arg[iarg],"yhi") == 0) faces[YHI] = 1;
    else if (strcmp(arg[iarg],"zlo") == 0) faces[ZLO] = 1;
    else if (strcmp(arg[iarg],"zhi") == 0) faces[ZHI] = 1;
    else break;
  }

  // error check

  if (domain->dimension == 2 && (faces[ZLO] || faces[ZHI])) 
    error->all(FLERR,"Cannot use fix inflow in z dimension for 2d simulation");

  // parse optional args

  np = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"n") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix inflow command");
      np = ATOBIGINT(arg[iarg+1]);
      if (np <= 0) error->all(FLERR,"Illegal fix inflow command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix inflow command");
  }

  // RNG

  int me = comm->me;
  random = new RanPark(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,me,100);

  // local storage

  cellface = NULL;
}

/* ---------------------------------------------------------------------- */

FixInflow::~FixInflow()
{
  delete random;
  memory->destroy(cellface);
}

/* ---------------------------------------------------------------------- */

int FixInflow::setmask()
{
  int mask = 0;
  mask |= START_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixInflow::init()
{
  // cannot inflow onto periodic boundary

  for (int i = 0; i < 6; i++)
    if (faces[i] && domain->bflag[i] == PERIODIC)
      error->all(FLERR,"Cannot use fix inflow on periodic boundary");

  // ncf = # of my cell/face pairs to insert onto

  int dimension = domain->dimension;

  Grid::OneCell *cells = grid->cells;
  int *mycells = grid->mycells;
  int nglocal = grid->nlocal;

  int icell;

  ncf = 0;
  for (int i = 0; i < nglocal; i++) {
    icell = mycells[i];
    if (cells[icell].neigh[XLO] < 0 && faces[XLO]) ncf++;
    if (cells[icell].neigh[XHI] < 0 && faces[XHI]) ncf++;
    if (cells[icell].neigh[YLO] < 0 && faces[XLO]) ncf++;
    if (cells[icell].neigh[YHI] < 0 && faces[XHI]) ncf++;
    if (dimension == 3) {
      if (cells[icell].neigh[ZLO] < 0 && faces[ZLO]) ncf++;
      if (cells[icell].neigh[ZHI] < 0 && faces[ZHI]) ncf++;
    }
  }

  // cellface = cell/face pair for all inserts I perform on my grid cells

  memory->destroy(cellface);
  memory->create(cellface,ncf,2,"inflow:cellface");

  ncf = 0;
  for (int i = 0; i < nglocal; i++) {
    icell = mycells[i];
    if (cells[icell].neigh[XLO] < 0 && faces[XLO]) {
      cellface[ncf][0] = icell; cellface[ncf++][1] = XLO;
    }
    if (cells[icell].neigh[XHI] < 0 && faces[XHI]) {
      cellface[ncf][0] = icell; cellface[ncf++][1] = XHI;
    }
    if (cells[icell].neigh[YLO] < 0 && faces[XLO]) {
      cellface[ncf][0] = icell; cellface[ncf++][1] = YLO;
    }
    if (cells[icell].neigh[YHI] < 0 && faces[XHI]) {
      cellface[ncf][0] = icell; cellface[ncf++][1] = YHI;
    }
    if (dimension == 3) {
      if (cells[icell].neigh[ZLO] < 0 && faces[ZLO]) {
	cellface[ncf][0] = icell; cellface[ncf++][1] = ZLO;
      }
      if (cells[icell].neigh[ZHI] < 0 && faces[ZHI]) {
	cellface[ncf][0] = icell; cellface[ncf++][1] = ZHI;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixInflow::start_of_step()
{
  if (update->ntimestep % nevery) return;

  int dimension = domain->dimension;

  Particle::OnePart *particles = particle->particles;
  Grid::OneCell *cells = grid->cells;

  int nspecies = particle->mixture[imix]->nspecies;
  int *species = particle->mixture[imix]->species;
  double *cummulative = particle->mixture[imix]->cummulative;
  double **vstream = particle->mixture[imix]->vstream;
  double *vscale = particle->mixture[imix]->vscale;

  int ilocal,icell,iface,npercell,ispecies;
  double x[3],v[3];
  double vol,ntarget,rn,vn,vr,theta1,theta2;
  double *lo,*hi;

  for (int i = 0; i < ncf; i++) {
    icell = cellface[i][0];
    iface = cellface[i][1];

    // local or global cell?
    // worry about normal dir dotted into vstream?

    //double area = 1.0;
    //double vol = area*dt*vstream;
    //ntarget = nme * volsum/volme - nprev;
    npercell = static_cast<int> (ntarget);
    if (random->uniform() < ntarget-npercell) npercell++;

    for (int m = 0; m < npercell; m++) {
      rn = random->uniform();
      ispecies = 0;
      while (cummulative[ispecies] < rn) ispecies++;

      x[0] = lo[0] + random->uniform() * (hi[0]-lo[0]);
      x[1] = lo[1] + random->uniform() * (hi[1]-lo[1]);
      x[2] = lo[2] + random->uniform() * (hi[2]-lo[2]);
      if (dimension == 2) x[2] = 0.0;

      vn = vscale[ispecies] * random->gaussian();
      vr = vscale[ispecies] * random->gaussian();
      theta1 = MY_2PI * random->uniform();
      theta2 = MY_2PI * random->uniform();
	
      v[0] = vstream[ispecies][0] + vn*cos(theta1);
      v[1] = vstream[ispecies][1] + vr*sin(theta2);
      v[2] = vstream[ispecies][2] + vr*cos(theta2);

      particle->add_particle(0,ispecies,icell,x,v);
    }
  }
}

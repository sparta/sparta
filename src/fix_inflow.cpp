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

#include "mpi.h"
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
enum{NO,YES};

/* ---------------------------------------------------------------------- */

FixInflow::FixInflow(DSMC *dsmc, int narg, char **arg) :
  Fix(dsmc, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix inflow command");

  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;

  imix = particle->find_mixture(arg[2]);
  if (imix < 0) error->all(FLERR,"Fix inflow mixture ID does not exist");

  // flag specified faces

  faces[XLO] = faces[XHI] = faces[YLO] = faces[YHI] =
    faces[ZLO] = faces[ZHI] = 0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"all") == 0) {
      if (domain->dimension == 3)
	faces[XLO] = faces[XHI] = faces[YLO] = faces[YHI] =
	  faces[ZLO] = faces[ZHI] = 1;
      else faces[XLO] = faces[XHI] = faces[YLO] = faces[YHI] = 1;
    } else if (strcmp(arg[iarg],"xlo") == 0) faces[XLO] = 1;
    else if (strcmp(arg[iarg],"xhi") == 0) faces[XHI] = 1;
    else if (strcmp(arg[iarg],"ylo") == 0) faces[YLO] = 1;
    else if (strcmp(arg[iarg],"yhi") == 0) faces[YHI] = 1;
    else if (strcmp(arg[iarg],"zlo") == 0) faces[ZLO] = 1;
    else if (strcmp(arg[iarg],"zhi") == 0) faces[ZHI] = 1;
    else break;
    iarg++;
  }

  // optional args

  np = 0;
  nevery = 1;
  perspecies = YES;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"n") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix inflow command");
      np = atoi(arg[iarg+1]);
      if (np <= 0) error->all(FLERR,"Illegal fix inflow command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"nevery") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix inflow command");
      nevery = atoi(arg[iarg+1]);
      if (nevery <= 0) error->all(FLERR,"Illegal fix inflow command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"perspecies") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix inflow command");
      if (strcmp(arg[iarg+1],"yes") == 0) perspecies = YES;
      else if (strcmp(arg[iarg+1],"no") == 0) perspecies = NO;
      else error->all(FLERR,"Illegal fix inflow command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix inflow command");
  }

  // error check

  if (domain->dimension == 2 && (faces[ZLO] || faces[ZHI])) 
    error->all(FLERR,"Cannot use fix inflow in z dimension for 2d simulation");
  if (np > 0 && perspecies == YES) 
    error->all(FLERR,"Cannot use fix inflow n > 0 with perspecies yes");

  // RNG

  int me = comm->me;
  random = new RanPark(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,me,100);

  // local storage

  cellface = NULL;
  ncf = 0;
  
  // counters

  nsingle = ntotal = 0;
}

/* ---------------------------------------------------------------------- */

FixInflow::~FixInflow()
{
  delete random;

  for (int i = 0; i < ncf; i++) delete [] cellface[i].ntargetsp;
  memory->sfree(cellface);
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
  // cannot inflow thru periodic boundary

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

  // cellface = per-face data struct for all inserts performed on my grid cells
  // indot = dot product of vstream with outward face normal
  // 2d vs 3d adjusts lo[2],hi[2] and area

  memory->sfree(cellface);
  cellface = (CellFace *) memory->smalloc(ncf*sizeof(CellFace),
					  "inflow:cellface");

  int nspecies = particle->mixture[imix]->nspecies;
  int *species = particle->mixture[imix]->species;

  for (int i = 0; i < ncf; i++) cellface[i].ntargetsp = new double[nspecies];

  double nrho = particle->mixture[imix]->nrho;
  double *vstream = particle->mixture[imix]->vstream;
  double fnum = update->fnum;
  double dt = update->dt;

  double area,indot;

  ncf = 0;
  for (int i = 0; i < nglocal; i++) {
    icell = mycells[i];
    if (cells[icell].neigh[XLO] < 0 && faces[XLO]) {
      cellface[ncf].icell = icell;
      cellface[ncf].ndim = 0;
      cellface[ncf].pdim1 = 1;
      cellface[ncf].pdim2 = 2;
      cellface[ncf].lo[0] = cells[icell].lo[0];
      cellface[ncf].hi[0] = cells[icell].lo[0];
      cellface[ncf].lo[1] = cells[icell].lo[1];
      cellface[ncf].hi[1] = cells[icell].hi[1];
      cellface[ncf].lo[2] = cells[icell].lo[2];
      cellface[ncf].hi[2] = cells[icell].hi[2];
      cellface[ncf].normal[0] = -1.0;
      cellface[ncf].normal[1] = 0.0;
      cellface[ncf].normal[2] = 0.0;
      
      area = (cells[icell].hi[1]-cells[icell].lo[1]) * 
	(cells[icell].hi[2]-cells[icell].lo[2]);

      if (dimension == 2) {
	cellface[ncf].lo[2] = 0.0;
	cellface[ncf].hi[2] = 0.0;
	area = (cells[icell].hi[1]-cells[icell].lo[1]);
      }

      cellface[ncf].ntarget = 0.0;
      indot = vstream[0]*cellface[ncf].normal[0] +
	vstream[1]*cellface[ncf].normal[1] +
	vstream[2]*cellface[ncf].normal[2];
      for (int isp = 0; isp < nspecies; isp++) {
	cellface[ncf].ntargetsp[isp] = mol_inflow(isp,indot);
	printf("SETUP %d %g\n",isp,cellface[ncf].ntargetsp[isp]);
	cellface[ncf].ntargetsp[isp] *= nrho*area*dt / fnum;
	printf("  SETUP %d %g\n",isp,cellface[ncf].ntargetsp[isp]);
	cellface[ncf].ntarget += cellface[ncf].ntargetsp[isp];
      }
      ncf++;
    }
    
    if (cells[icell].neigh[XHI] < 0 && faces[XHI]) {
      cellface[ncf].icell = icell;
      cellface[ncf].ndim = 0;
      cellface[ncf].pdim1 = 1;
      cellface[ncf].pdim2 = 2;
      cellface[ncf].lo[0] = cells[icell].hi[0];
      cellface[ncf].hi[0] = cells[icell].hi[0];
      cellface[ncf].lo[1] = cells[icell].lo[1];
      cellface[ncf].hi[1] = cells[icell].hi[1];
      cellface[ncf].lo[2] = cells[icell].lo[2];
      cellface[ncf].hi[2] = cells[icell].hi[2];
      cellface[ncf].normal[0] = 1.0;
      cellface[ncf].normal[1] = 0.0;
      cellface[ncf].normal[2] = 0.0;
      
      area = (cells[icell].hi[1]-cells[icell].lo[1]) * 
	(cells[icell].hi[2]-cells[icell].lo[2]);

      if (dimension == 2) {
	cellface[ncf].lo[2] = 0.0;
	cellface[ncf].hi[2] = 0.0;
	area = (cells[icell].hi[1]-cells[icell].lo[1]);
      }
      
      cellface[ncf].ntarget = 0.0;
      indot = vstream[0]*cellface[ncf].normal[0] +
	vstream[1]*cellface[ncf].normal[1] +
	vstream[2]*cellface[ncf].normal[2];
      for (int isp = 0; isp < nspecies; isp++) {
	cellface[ncf].ntargetsp[isp] = mol_inflow(isp,indot);
	cellface[ncf].ntargetsp[isp] *= nrho*area*dt / fnum;
	cellface[ncf].ntarget += cellface[ncf].ntargetsp[isp];
      }
      ncf++;
    }

    if (cells[icell].neigh[YLO] < 0 && faces[YLO]) {
      cellface[ncf].icell = icell;
      cellface[ncf].ndim = 1;
      cellface[ncf].pdim1 = 0;
      cellface[ncf].pdim2 = 2;
      cellface[ncf].lo[0] = cells[icell].lo[0];
      cellface[ncf].hi[0] = cells[icell].hi[0];
      cellface[ncf].lo[1] = cells[icell].lo[1];
      cellface[ncf].hi[1] = cells[icell].lo[1];
      cellface[ncf].lo[2] = cells[icell].lo[2];
      cellface[ncf].hi[2] = cells[icell].hi[2];
      cellface[ncf].normal[0] = 0.0;
      cellface[ncf].normal[1] = -1.0;
      cellface[ncf].normal[2] = 0.0;
      
      area = (cells[icell].hi[0]-cells[icell].lo[0]) * 
	(cells[icell].hi[2]-cells[icell].lo[2]);

      if (dimension == 2) {
	cellface[ncf].lo[2] = 0.0;
	cellface[ncf].hi[2] = 0.0;
	area = (cells[icell].hi[0]-cells[icell].lo[0]);
      }
      
      cellface[ncf].ntarget = 0.0;
      indot = vstream[0]*cellface[ncf].normal[0] +
	vstream[1]*cellface[ncf].normal[1] +
	vstream[2]*cellface[ncf].normal[2];
      for (int isp = 0; isp < nspecies; isp++) {
	cellface[ncf].ntargetsp[isp] = mol_inflow(isp,indot);
	cellface[ncf].ntargetsp[isp] *= nrho*area*dt / fnum;
	cellface[ncf].ntarget += cellface[ncf].ntargetsp[isp];
      }
      ncf++;
    }

    if (cells[icell].neigh[YHI] < 0 && faces[YHI]) {
      cellface[ncf].icell = icell;
      cellface[ncf].ndim = 1;
      cellface[ncf].pdim1 = 0;
      cellface[ncf].pdim2 = 2;
      cellface[ncf].lo[0] = cells[icell].lo[0];
      cellface[ncf].hi[0] = cells[icell].hi[0];
      cellface[ncf].lo[1] = cells[icell].hi[1];
      cellface[ncf].hi[1] = cells[icell].hi[1];
      cellface[ncf].lo[2] = cells[icell].lo[2];
      cellface[ncf].hi[2] = cells[icell].hi[2];
      cellface[ncf].normal[0] = 0.0;
      cellface[ncf].normal[1] = 1.0;
      cellface[ncf].normal[2] = 0.0;
      
      area = (cells[icell].hi[0]-cells[icell].lo[0]) * 
	(cells[icell].hi[2]-cells[icell].lo[2]);

      if (dimension == 2) {
	cellface[ncf].lo[2] = 0.0;
	cellface[ncf].hi[2] = 0.0;
	area = (cells[icell].hi[0]-cells[icell].lo[0]);
      }
      
      cellface[ncf].ntarget = 0.0;
      indot = vstream[0]*cellface[ncf].normal[0] +
	vstream[1]*cellface[ncf].normal[1] +
	vstream[2]*cellface[ncf].normal[2];
      for (int isp = 0; isp < nspecies; isp++) {
	cellface[ncf].ntargetsp[isp] = mol_inflow(isp,indot);
	cellface[ncf].ntargetsp[isp] *= nrho*area*dt / fnum;
	cellface[ncf].ntarget += cellface[ncf].ntargetsp[isp];
      }
      ncf++;
    }

    if (cells[icell].neigh[ZLO] < 0 && faces[ZLO]) {
      cellface[ncf].icell = icell;
      cellface[ncf].ndim = 2;
      cellface[ncf].pdim1 = 0;
      cellface[ncf].pdim2 = 1;
      cellface[ncf].lo[0] = cells[icell].lo[0];
      cellface[ncf].hi[0] = cells[icell].hi[0];
      cellface[ncf].lo[1] = cells[icell].lo[1];
      cellface[ncf].hi[1] = cells[icell].hi[1];
      cellface[ncf].lo[2] = cells[icell].lo[2];
      cellface[ncf].hi[2] = cells[icell].lo[2];
      cellface[ncf].normal[0] = 0.0;
      cellface[ncf].normal[1] = 0.0;
      cellface[ncf].normal[2] = -1.0;
      
      area = (cells[icell].hi[0]-cells[icell].lo[0]) * 
	(cells[icell].hi[1]-cells[icell].lo[1]);
      
      cellface[ncf].ntarget = 0.0;
      indot = vstream[0]*cellface[ncf].normal[0] +
	vstream[1]*cellface[ncf].normal[1] +
	vstream[2]*cellface[ncf].normal[2];
      for (int isp = 0; isp < nspecies; isp++) {
	cellface[ncf].ntargetsp[isp] = mol_inflow(isp,indot);
	cellface[ncf].ntargetsp[isp] *= nrho*area*dt / fnum;
	cellface[ncf].ntarget += cellface[ncf].ntargetsp[isp];
      }
      ncf++;
    }

    if (cells[icell].neigh[ZHI] < 0 && faces[ZHI]) {
      cellface[ncf].icell = icell;
      cellface[ncf].ndim = 2;
      cellface[ncf].pdim1 = 0;
      cellface[ncf].pdim2 = 1;
      cellface[ncf].lo[0] = cells[icell].lo[0];
      cellface[ncf].hi[0] = cells[icell].hi[0];
      cellface[ncf].lo[1] = cells[icell].lo[1];
      cellface[ncf].hi[1] = cells[icell].hi[1];
      cellface[ncf].lo[2] = cells[icell].hi[2];
      cellface[ncf].hi[2] = cells[icell].hi[2];
      cellface[ncf].normal[0] = 0.0;
      cellface[ncf].normal[1] = 0.0;
      cellface[ncf].normal[2] = 1.0;
      
      area = (cells[icell].hi[0]-cells[icell].lo[0]) * 
	(cells[icell].hi[1]-cells[icell].lo[1]);
      
      cellface[ncf].ntarget = 0.0;
      indot = vstream[0]*cellface[ncf].normal[0] +
	vstream[1]*cellface[ncf].normal[1] +
	vstream[2]*cellface[ncf].normal[2];
      for (int isp = 0; isp < nspecies; isp++) {
	cellface[ncf].ntargetsp[isp] = mol_inflow(isp,indot);
	cellface[ncf].ntargetsp[isp] *= nrho*area*dt / fnum;
	cellface[ncf].ntarget += cellface[ncf].ntargetsp[isp];
      }
      ncf++;
    }
  }
  
  // if Np > 0, npercell = # of insertions per cellface pair
  // set nthresh so as to achieve exactly Np insertions
  // cells > cells_with_no_extra need to insert 1 extra particle

  if (np > 0) {
    int all,nupto;
    MPI_Allreduce(&ncf,&all,1,MPI_INT,MPI_SUM,world);
    npercell = np / all;
    int cells_with_no_extra = all - (np % all);
    MPI_Scan(&ncf,&nupto,1,MPI_INT,MPI_SUM,world);
    if (cells_with_no_extra < nupto-ncf) nthresh = 0;
    else if (cells_with_no_extra >= nupto) nthresh = ncf;
    else nthresh = cells_with_no_extra - (nupto-ncf);
  }

  // cummulative counter

  ntotal = 0;
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
  double *vstream = particle->mixture[imix]->vstream;
  double *vscale = particle->mixture[imix]->vscale;

  int icell,ninsert,isp,ndim,pdim1,pdim2;
  double *lo,*hi,*normal;
  double x[3],v[3];
  double indot,rn,ntarget;
  double beta_un,normalized_distbn_fn,theta;

  // insert molecules by cell/face pair
  // ntarget/ninsert is either perspecies or for all species
  // for one molecule:
  //   x = random position on face
  //   v = randomized thermal velocity + vstream
  //       first stage: normal dimension (ndim)
  //       second stage: parallel dimensions (pdim1,pdim2)

  nsingle = 0;

  for (int i = 0; i < ncf; i++) {
    icell = cellface[i].icell;
    ndim = cellface[i].ndim;
    pdim1 = cellface[i].pdim1;
    pdim2 = cellface[i].pdim2;
    lo = cellface[i].lo;
    hi = cellface[i].hi;
    normal = cellface[i].normal;

    indot = vstream[0]*normal[0] + vstream[1]*normal[1] + vstream[2]*normal[2];

    if (perspecies == YES) {
      for (isp = 0; isp < nspecies; isp++) {
	ntarget = cellface[i].ntargetsp[isp];
	ninsert = static_cast<int> (ntarget);
	if (random->uniform() < ntarget-ninsert) ninsert++;

	//printf("AAA %d %d %g\n",icell,ninsert,ntarget);

	for (int m = 0; m < ninsert; m++) {
	  x[0] = lo[0] + random->uniform() * (hi[0]-lo[0]);
	  x[1] = lo[1] + random->uniform() * (hi[1]-lo[1]);
	  x[2] = lo[2] + random->uniform() * (hi[2]-lo[2]);

	  do {
	    do beta_un = (6.0*random->gaussian() - 3.0);
	    while (beta_un + indot < 0.0);
	    normalized_distbn_fn = 2.0 * (beta_un + indot) / 
	      (indot + sqrt(indot*indot + 2.0)) *
	      exp(0.5 + (0.5*indot)*(indot-sqrt(indot*indot + 2.0)) - 
		  beta_un*beta_un);
	  } while (normalized_distbn_fn < random->gaussian());
	  
	  v[ndim] = beta_un + vstream[0];
	  
	  theta = MY_PI * random->gaussian();
	  v[pdim1] = vscale[isp]*sin(theta) + vstream[1];
	  v[pdim2] = vscale[isp]*cos(theta) + vstream[2]; 
	  
	  particle->add_particle(0,isp,icell,x,v);
	}

	nsingle += ninsert;
      }

    } else {
      if (np == 0) {
	ntarget = cellface[i].ntarget;
	ninsert = static_cast<int> (ntarget);
	if (random->uniform() < ntarget-ninsert) ninsert++;
      } else {
	ninsert = npercell;
	if (i >= nthresh) ninsert++;
      }

      printf("AAA %d %d\n",icell,ninsert);

      for (int m = 0; m < ninsert; m++) {
	rn = random->uniform();
	isp = 0;
	while (cummulative[isp] < rn) isp++;

	x[0] = lo[0] + random->uniform() * (hi[0]-lo[0]);
	x[1] = lo[1] + random->uniform() * (hi[1]-lo[1]);
	x[2] = lo[2] + random->uniform() * (hi[2]-lo[2]);

	do {
	  do {
	    beta_un = (6.0*random->gaussian() - 3.0);
	    printf("CCC %g %g\n",beta_un,indot);
	  } while (beta_un + indot < 0.0);
	  normalized_distbn_fn = 2.0 * (beta_un + indot) / 
	    (indot + sqrt(indot*indot + 2.0)) *
	    exp(0.5 + (0.5*indot)*(indot-sqrt(indot*indot + 2.0)) - 
		beta_un*beta_un);
	  printf("BBB %g\n",normalized_distbn_fn);
	} while (normalized_distbn_fn < random->gaussian());
	
	v[ndim] = beta_un + vstream[0];
	
	theta = MY_PI * random->gaussian();
	v[pdim1] = vscale[isp]*sin(theta) + vstream[1];
	v[pdim2] = vscale[isp]*cos(theta) + vstream[2]; 

	particle->add_particle(0,isp,icell,x,v);
      }

      nsingle += ninsert;
    }
  }

  ntotal += nsingle;
}

/* ----------------------------------------------------------------------
   calculate flux of particles of species ISP entering a grid cell
   see Bird 1994, eq 4.22
------------------------------------------------------------------------- */

double FixInflow::mol_inflow(int isp, double indot)
{
  double *vscale = particle->mixture[imix]->vscale;
  double *fraction = particle->mixture[imix]->fraction;

  double inward_number_flux = 0.0;
  if (indot >= 0.0) {
    inward_number_flux = vscale[isp] * fraction[isp] *
      (exp(-indot*indot) + sqrt(MY_PI)*indot*(1.0 + erf(indot)));
    inward_number_flux += random->gaussian();
  }

  return inward_number_flux;
}

/* ----------------------------------------------------------------------
   return one-step or total count of particle insertions
------------------------------------------------------------------------- */

double FixInflow::compute_vector(int i)
{
  double one,all;
  
  if (i == 0) one = nsingle;
  else one = ntotal;
  MPI_Allreduce(&one,&all,1,MPI_DOUBLE,MPI_SUM,world);
  return all;
}

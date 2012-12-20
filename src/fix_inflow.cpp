/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2012) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "create_molecules.h"
#include "fix_inflow.h"
#include "update.h"
#include "particle.h"
#include "mixture.h"
#include "domain.h"
#include "grid.h"
#include "surf.h"
#include "comm.h"
#include "geometry.h"
#include "random_mars.h"
#include "random_park.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};         // same as Domain
enum{PERIODIC,OUTFLOW,REFLECT,SURFACE,AXISYM};  // same as Domain
enum{CELLUNKNOWN,CELLOUTSIDE,CELLINSIDE,CELLOVERLAP};   // same as Grid
enum{CORNERUNKNOWN,CORNEROUTSIDE,CORNERINSIDE,CORNEROVERLAP};  // same as Grid
enum{NO,YES};

/* ---------------------------------------------------------------------- */

FixInflow::FixInflow(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
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
  ncf = ncfmax = 0;
  
  // counters

  nsingle = ntotal = 0;
}

/* ---------------------------------------------------------------------- */

FixInflow::~FixInflow()
{
  delete random;

  for (int i = 0; i < ncfmax; i++) delete [] cellface[i].ntargetsp;
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
  int i,j,m,n,isp;

  // corners[i][j] = J corner points of face I of a grid cell
  // works for 2d quads and 3d hexes
  
  int corners[6][4] = {{0,2,4,6}, {1,3,5,7}, {0,1,4,5}, {2,3,6,7}, 
		       {0,1,2,3}, {4,5,6,7}};

  int nface_pts = 4;
  if (domain->dimension == 2) nface_pts = 2;

  // cannot inflow thru periodic boundary

  for (i = 0; i < 6; i++)
    if (faces[i] && domain->bflag[i] == PERIODIC)
      error->all(FLERR,"Cannot use fix inflow on periodic boundary");

  // allow[I][J] = 1 if my local cell I, face J allows insertions
  // only allow if face adjoins global boundary with inflow defined
  // if cell is CELLOUTSIDE, allow face
  // if cell is CELLINSIDE, disallow face
  // if cell is CELLOVERLAP:
  //   allow if any face corner point is CORNEROUTSIDE and none is CORNERINSIDE
  //   disallow if any pt of any cell line/tri touches face

  int dimension = domain->dimension;
  Surf::Point *pts = surf->pts;
  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  Grid::OneCell *cells = grid->cells;
  int **csurfs = grid->csurfs;
  int **cflags = grid->cflags;
  int *mycells = grid->mycells;
  int nglocal = grid->nlocal;

  int **allow;
  memory->create(allow,nglocal,6,"inflow:allow");

  int *flags;
  int icell,dim,extflag;
  double value;

  for (int m = 0; m < nglocal; m++) {
    icell = mycells[m];
    for (i = 0; i < 6; i++) {
      if (faces[i] && cells[icell].neigh[i] < 0) {
	if (cells[icell].type == CELLOUTSIDE) allow[m][i] = 1;
	else if (cells[icell].type == CELLINSIDE) allow[m][i] = 0;
	else if (cells[icell].type == CELLOVERLAP) {
	  allow[m][i] = 1;
	  flags = cflags[i];

	  extflag = 0;
	  for (j = 0; j < nface_pts; j++) {
	    if (flags[corners[i][j]] == CORNEROUTSIDE) extflag = 1;
	    else if (flags[corners[i][j]] == CORNERINSIDE) allow[m][i] = 0;
	  }
	  if (!extflag) allow[m][i] = 0;
	  
	  if (allow[m][i]) {
	    if (dimension == 2) {
	      for (j = 0; j < cells[icell].nsurf; j++) {
		n = csurfs[icell][j];
		if (Geometry::
		    line_quad_face_touch(pts[lines[n].p1].x,
					 pts[lines[n].p2].x,
					 i,cells[icell].lo,cells[icell].hi)) {
		  allow[m][i] = 0;
		  break;
		}
	      }
	    } else {
	      for (j = 0; j < cells[icell].nsurf; j++) {
		n = csurfs[icell][j];
		if (Geometry::
		    tri_hex_face_touch(pts[tris[n].p1].x,
				       pts[tris[n].p2].x,
				       pts[tris[n].p3].x,
				       i,cells[icell].lo,cells[icell].hi)) {
		  allow[m][i] = 0;
		  break;
		}
	      }
	    }
	  }
	}
      } else allow[m][i] = 0;
    }
  }

  // ncf = # of my cell/face pairs to insert onto
  // some may be eliminated later if no particles are actually inserted

  ncf = 0;
  for (m = 0; m < nglocal; m++) {
    icell = mycells[m];
    if (allow[m][XLO]) ncf++;
    if (allow[m][XHI]) ncf++;
    if (allow[m][YLO]) ncf++;
    if (allow[m][YHI]) ncf++;
    if (allow[m][ZLO]) ncf++;
    if (allow[m][ZHI]) ncf++;
  }

  // cellface = per-face data struct for all inserts performed on my grid cells
  // indot = dot product of vstream with outward face normal
  // skip the cellface if indot < 0.0, since no particles will be inserted
  // 2d vs 3d adjusts lo[2],hi[2] and area

  for (i = 0; i < ncfmax; i++) delete [] cellface[i].ntargetsp;
  memory->sfree(cellface);
  ncfmax = ncf;
  cellface = (CellFace *) memory->smalloc(ncfmax*sizeof(CellFace),
					  "inflow:cellface");

  int nspecies = particle->mixture[imix]->nspecies;
  int *species = particle->mixture[imix]->species;

  for (i = 0; i < ncfmax; i++)
    cellface[i].ntargetsp = new double[nspecies];

  double nrho = particle->mixture[imix]->nrho;
  double *vstream = particle->mixture[imix]->vstream;
  double fnum = update->fnum;
  double dt = update->dt;

  double area,indot;

  ncf = 0;
  for (m = 0; m < nglocal; m++) {
    icell = mycells[m];
    if (allow[m][XLO]) {
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
      if (indot < 0.0) continue;

      for (isp = 0; isp < nspecies; isp++) {
	cellface[ncf].ntargetsp[isp] = mol_inflow(isp,indot);
	cellface[ncf].ntargetsp[isp] *= nrho*area*dt / fnum;
	cellface[ncf].ntarget += cellface[ncf].ntargetsp[isp];
      }
      ncf++;
    }
    
    if (allow[m][XHI]) {
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
      if (indot < 0.0) continue;

      for (isp = 0; isp < nspecies; isp++) {
	cellface[ncf].ntargetsp[isp] = mol_inflow(isp,indot);
	cellface[ncf].ntargetsp[isp] *= nrho*area*dt / fnum;
	cellface[ncf].ntarget += cellface[ncf].ntargetsp[isp];
      }
      ncf++;
    }

    if (allow[m][YLO]) {
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
      if (indot < 0.0) continue;

      for (isp = 0; isp < nspecies; isp++) {
	cellface[ncf].ntargetsp[isp] = mol_inflow(isp,indot);
	cellface[ncf].ntargetsp[isp] *= nrho*area*dt / fnum;
	cellface[ncf].ntarget += cellface[ncf].ntargetsp[isp];
      }
      ncf++;
    }

    if (allow[m][YHI]) {
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
      if (indot < 0.0) continue;

      for (isp = 0; isp < nspecies; isp++) {
	cellface[ncf].ntargetsp[isp] = mol_inflow(isp,indot);
	cellface[ncf].ntargetsp[isp] *= nrho*area*dt / fnum;
	cellface[ncf].ntarget += cellface[ncf].ntargetsp[isp];
      }
      ncf++;
    }

    if (allow[m][ZLO]) {
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
      cellface[ncf].normal[2] = 1.0;
      
      area = (cells[icell].hi[0]-cells[icell].lo[0]) * 
	(cells[icell].hi[1]-cells[icell].lo[1]);
      
      cellface[ncf].ntarget = 0.0;
      indot = vstream[0]*cellface[ncf].normal[0] +
	vstream[1]*cellface[ncf].normal[1] +
	vstream[2]*cellface[ncf].normal[2];
      if (indot < 0.0) continue;

      for (isp = 0; isp < nspecies; isp++) {
	cellface[ncf].ntargetsp[isp] = mol_inflow(isp,indot);
	cellface[ncf].ntargetsp[isp] *= nrho*area*dt / fnum;
	cellface[ncf].ntarget += cellface[ncf].ntargetsp[isp];
      }
      ncf++;
    }

    if (allow[m][ZHI]) {
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
      cellface[ncf].normal[2] = -1.0;
      
      area = (cells[icell].hi[0]-cells[icell].lo[0]) * 
	(cells[icell].hi[1]-cells[icell].lo[1]);
      
      cellface[ncf].ntarget = 0.0;
      indot = vstream[0]*cellface[ncf].normal[0] +
	vstream[1]*cellface[ncf].normal[1] +
	vstream[2]*cellface[ncf].normal[2];
      if (indot < 0.0) continue;

      for (isp = 0; isp < nspecies; isp++) {
	cellface[ncf].ntargetsp[isp] = mol_inflow(isp,indot);
	cellface[ncf].ntargetsp[isp] *= nrho*area*dt / fnum;
	cellface[ncf].ntarget += cellface[ncf].ntargetsp[isp];
      }
      ncf++;
    }
  }

  // clean up

  memory->destroy(allow);

  // if Np > 0, npercell = # of insertions per active cellface pair
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

  int icell,ninsert,isp,ndim,pdim1,pdim2,ivib;
  double *lo,*hi,*normal;
  double x[3],v[3];
  double indot,rn,ntarget;
  double beta_un,normalized_distbn_fn,theta,erot;

  // insert molecules by cell/face pair
  // ntarget/ninsert is either perspecies or for all species
  // for one molecule:
  //   x = random position on face
  //   v = randomized thermal velocity + vstream
  //       first stage: normal dimension (ndim)
  //       second stage: parallel dimensions (pdim1,pdim2)

  // NOTE: if allow particle insertion on backflow boundaries
  //       then should worry about do while loops spinning endlessly
  //       due to difficulty of generating a valid particle to insert
  //       may especially happen if force Np insertions on backflow boundary

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
	  } while (normalized_distbn_fn < random->uniform());
	  
	  v[ndim] = beta_un + vstream[0];
	  
	  theta = MY_PI * random->gaussian();
	  v[pdim1] = vscale[isp]*sin(theta) + vstream[1];
	  v[pdim2] = vscale[isp]*cos(theta) + vstream[2]; 
          erot = particle->erot(isp,random);
          ivib = particle->evib(isp,random);
	  particle->add_particle(0,isp,icell,x,v,erot,ivib);
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

      //printf("AAA %d %d\n",icell,ninsert);

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
	    //printf("CCC %g %g\n",beta_un,indot);
	  } while (beta_un + indot < 0.0);
	  normalized_distbn_fn = 2.0 * (beta_un + indot) / 
	    (indot + sqrt(indot*indot + 2.0)) *
	    exp(0.5 + (0.5*indot)*(indot-sqrt(indot*indot + 2.0)) - 
		beta_un*beta_un);
	  //printf("BBB %g\n",normalized_distbn_fn);
	} while (normalized_distbn_fn < random->uniform());
	
	v[ndim] = beta_un*vscale[isp] + vstream[0];
	
	theta = MY_PI * random->uniform();
	v[pdim1] = vscale[isp]*sin(theta) + vstream[1];
	v[pdim2] = vscale[isp]*cos(theta) + vstream[2]; 
        erot = particle->erot(isp,random);
        ivib = particle->evib(isp,random);
	particle->add_particle(0,isp,icell,x,v,erot,ivib);
      }

      nsingle += ninsert;
    }
  }

  ntotal += nsingle;
}

/* ----------------------------------------------------------------------
   calculate flux of particles of species ISP entering a grid cell
   see Bird 1994, eq 4.22
   NOTE: could add option to insert particles on backflow boundaries
         when indot < 0.0
------------------------------------------------------------------------- */

double FixInflow::mol_inflow(int isp, double indot)
{
  double *vscale = particle->mixture[imix]->vscale;
  double *fraction = particle->mixture[imix]->fraction;

  if (indot < 0.0) error->one(FLERR,"Fix inflow used on outflow boundary");

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

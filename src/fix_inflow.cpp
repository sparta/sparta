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

#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_inflow.h"
#include "update.h"
#include "particle.h"
#include "mixture.h"
#include "domain.h"
#include "modify.h"
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
enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};           // same as Grid
enum{PKEEP,PINSERT,PDONE,PDISCARD,PENTRY,PEXIT,PSURF};   // several files
enum{NCHILD,NPARENT,NUNKNOWN,NPBCHILD,NPBPARENT,NPBUNKNOWN,NBOUND};  // Grid
enum{NO,YES};

#define DELTAFACE 10              // no smaller than 6
#define DELTAGRID 1000            // must be bigger than split cells per cell

/* ---------------------------------------------------------------------- */

FixInflow::FixInflow(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix inflow command");

  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  gridmigrate = 1;

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

  // error checks

  if (domain->dimension == 2 && (faces[ZLO] || faces[ZHI])) 
    error->all(FLERR,"Cannot use fix inflow in z dimension for 2d simulation");
  if (domain->axisymmetric && faces[YLO]) 
    error->all(FLERR,"Cannot use fix inflow on ylo face for "
               "axisymmetric model");
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
  c2f = NULL;
  nglocal = nglocalmax = 0;

  // counters

  nsingle = ntotal = 0;
}

/* ---------------------------------------------------------------------- */

FixInflow::~FixInflow()
{
  delete random;

  for (int i = 0; i < ncfmax; i++) delete [] cellface[i].ntargetsp;
  memory->sfree(cellface);
  memory->destroy(c2f);
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
  int i,j,m,n,isp,icell;
  double xface[3];

  int nspecies = particle->mixture[imix]->nspecies;
  double nrho = particle->mixture[imix]->nrho;
  double *vstream = particle->mixture[imix]->vstream;
  double *vscale = particle->mixture[imix]->vscale;
  double *fraction = particle->mixture[imix]->fraction;
  double fnum = update->fnum;
  double dt = update->dt;

  particle->exist = 1;

  // run-time error check

  if (domain->axisymmetric && faces[YHI] & vstream[1] != 0.0)
    error->all(FLERR,"Cannot use fix inflow on yhi for axisymmetric model "
               "if streaming velocity has a y-component");

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

  // warn if any inflow face does not have an inward normal
  //   in direction of streaming velocity

  double normal[3];
  int flag = 0;

  for (i = 0; i < 6; i++) {
    if (!faces[i]) continue;
    normal[0] = normal[1] = normal[2] = 0.0;
    if (i % 2 == 0) normal[i/2] = 1.0;
    else normal[i/2] = -1.0;
    double indot = vstream[0]*normal[0] + vstream[1]*normal[1] + 
      vstream[2]*normal[2];
    if (indot < 0.0) flag = 1;
  }

  if (flag && comm->me == 0)
    error->warning(FLERR,
                   "One or more fix inflow faces oppose streaming velocity");

  // c2f[I][J] = 1 if my local parent cell I, face J allows insertions
  // only allow if face adjoins global boundary with inflow defined
  // if cell is OUTSIDE, allow face
  // if cell is INSIDE, disallow face
  // if cell is OVERLAP:
  //   allow if any face corner point is OUTSIDE and none is INSIDE
  //   disallow if any pt of any cell line/tri touches face

  int dimension = domain->dimension;
  Surf::Point *pts = surf->pts;
  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  nglocal = grid->nlocal;

  // realloc c2f if necessary = pointers from cells/faces to cellface

  if (nglocal > nglocalmax) {
    memory->destroy(c2f);
    nglocalmax = nglocal;
    memory->create(c2f,nglocalmax,6,"inflow:c2f");
  }

  int *flags;
  int nmask,dim,extflag;
  double value;

  // set c2f to 1 if face is eligible for insertion, else 0

  for (icell = 0; icell < nglocal; icell++) {
    nmask = cells[icell].nmask;
    for (i = 0; i < 6; i++) {
      if (faces[i] && grid->neigh_decode(nmask,i) == NBOUND && 
          cells[icell].nsplit >= 1) {
	if (cinfo[icell].type == OUTSIDE) c2f[icell][i] = 1;
	else if (cinfo[icell].type == INSIDE) c2f[icell][i] = 0;
	else if (cinfo[icell].type == OVERLAP) {
	  c2f[icell][i] = 1;
	  flags = cinfo[icell].corner;

	  extflag = 0;
	  for (j = 0; j < nface_pts; j++) {
	    if (flags[corners[i][j]] == OUTSIDE) extflag = 1;
	    else if (flags[corners[i][j]] == INSIDE) c2f[icell][i] = 0;
	  }
	  if (!extflag) c2f[icell][i] = 0;

	  if (c2f[icell][i]) {
	    if (dimension == 2) {
	      for (j = 0; j < cells[icell].nsurf; j++) {
		n = cells[icell].csurfs[j];
		if (Geometry::
		    line_quad_face_touch(pts[lines[n].p1].x,
					 pts[lines[n].p2].x,
					 i,cells[icell].lo,cells[icell].hi)) {
		  c2f[icell][i] = 0;
		  break;
		}
	      }
	    } else {
	      for (j = 0; j < cells[icell].nsurf; j++) {
		n = cells[icell].csurfs[j];
		if (Geometry::
		    tri_hex_face_touch(pts[tris[n].p1].x,
				       pts[tris[n].p2].x,
				       pts[tris[n].p3].x,
				       i,cells[icell].lo,cells[icell].hi)) {
		  c2f[icell][i] = 0;
		  break;
		}
	      }
	    }
	  }
	}
      } else c2f[icell][i] = 0;
    }
  }

  // ncf = # of my cell/face pairs to insert onto
  // some may be eliminated later if no particles are actually inserted

  ncf = 0;
  for (icell = 0; icell < nglocal; icell++) {
    if (c2f[icell][XLO]) ncf++;
    if (c2f[icell][XHI]) ncf++;
    if (c2f[icell][YLO]) ncf++;
    if (c2f[icell][YHI]) ncf++;
    if (c2f[icell][ZLO]) ncf++;
    if (c2f[icell][ZHI]) ncf++;
  }

  // cellface = per-face data struct for all inserts performed on my grid cells
  // reallocate cellface since nspecies count of mixture may have changed
  // indot = dot product of vstream with inward face normal
  // skip cellface if indot < 0.0, to not allow any particles to be inserted
  // 2d vs 3d adjusts lo[2],hi[2] and area
  // convert c2f[I][J] = index into cellface of Jth face of cell I, -1 if none

  for (int i = 0; i < ncfmax; i++) delete [] cellface[i].ntargetsp;
  memory->sfree(cellface);
  ncfmax = ncf;
  cellface = (CellFace *) memory->smalloc(ncfmax*sizeof(CellFace),
                                          "inflow:cellface");
  memset(cellface,0,ncfmax*sizeof(CellFace));

  for (i = 0; i < ncfmax; i++)
    cellface[i].ntargetsp = new double[nspecies];

  double area,indot;

  ncf = 0;
  for (icell = 0; icell < nglocal; icell++) {

    if (c2f[icell][XLO]) {
      cellface[ncf].id = cells[icell].id;
      if (cells[icell].nsplit > 1) cellface[ncf].pcell = split(icell,XLO);
      else cellface[ncf].pcell = icell;
      cellface[ncf].icell = icell;
      cellface[ncf].iface = XLO;
      cellface[ncf].ndim = 0;
      cellface[ncf].pdim = 1;
      cellface[ncf].qdim = 2;
      cellface[ncf].lo[0] = cells[icell].lo[0];
      cellface[ncf].hi[0] = cells[icell].lo[0];
      cellface[ncf].lo[1] = cells[icell].lo[1];
      cellface[ncf].hi[1] = cells[icell].hi[1];
      cellface[ncf].lo[2] = cells[icell].lo[2];
      cellface[ncf].hi[2] = cells[icell].hi[2];
      cellface[ncf].normal[0] = 1.0;
      cellface[ncf].normal[1] = 0.0;
      cellface[ncf].normal[2] = 0.0;
      
      if (dimension == 3) 
        area = (cells[icell].hi[1]-cells[icell].lo[1]) * 
          (cells[icell].hi[2]-cells[icell].lo[2]);
      else if (domain->axisymmetric)
        area = (cells[icell].hi[1]*cells[icell].hi[1] -
                cells[icell].lo[1]*cells[icell].lo[1])*MY_PI;
      else {
	cellface[ncf].lo[2] = 0.0;
	cellface[ncf].hi[2] = 0.0;
	area = (cells[icell].hi[1]-cells[icell].lo[1]);
      }

      cellface[ncf].ntarget = 0.0;
      indot = vstream[0]*cellface[ncf].normal[0] +
	vstream[1]*cellface[ncf].normal[1] +
	vstream[2]*cellface[ncf].normal[2];
      if (indot >= 0.0) {
        for (isp = 0; isp < nspecies; isp++) {
          cellface[ncf].ntargetsp[isp] = 
            mol_inflow(indot,vscale[isp],fraction[isp]);
          cellface[ncf].ntargetsp[isp] *= nrho*area*dt / fnum;
          cellface[ncf].ntargetsp[isp] /= cinfo[icell].weight;
          cellface[ncf].ntarget += cellface[ncf].ntargetsp[isp];
        }
        c2f[icell][XLO] = ncf++;
      } else c2f[icell][XLO] = -1;
    } else c2f[icell][XLO] = -1;
    
    if (c2f[icell][XHI]) {
      cellface[ncf].id = cells[icell].id;
      if (cells[icell].nsplit > 1) cellface[ncf].pcell = split(icell,XHI);
      else cellface[ncf].pcell = icell;
      cellface[ncf].icell = icell;
      cellface[ncf].iface = XHI;
      cellface[ncf].ndim = 0;
      cellface[ncf].pdim = 1;
      cellface[ncf].qdim = 2;
      cellface[ncf].lo[0] = cells[icell].hi[0];
      cellface[ncf].hi[0] = cells[icell].hi[0];
      cellface[ncf].lo[1] = cells[icell].lo[1];
      cellface[ncf].hi[1] = cells[icell].hi[1];
      cellface[ncf].lo[2] = cells[icell].lo[2];
      cellface[ncf].hi[2] = cells[icell].hi[2];
      cellface[ncf].normal[0] = -1.0;
      cellface[ncf].normal[1] = 0.0;
      cellface[ncf].normal[2] = 0.0;
      
      if (dimension == 3) 
        area = (cells[icell].hi[1]-cells[icell].lo[1]) * 
          (cells[icell].hi[2]-cells[icell].lo[2]);
      else if (domain->axisymmetric)
        area = (cells[icell].hi[1]*cells[icell].hi[1] -
                cells[icell].lo[1]*cells[icell].lo[1])*MY_PI;
      else {
	cellface[ncf].lo[2] = 0.0;
	cellface[ncf].hi[2] = 0.0;
	area = (cells[icell].hi[1]-cells[icell].lo[1]);
      }
      
      cellface[ncf].ntarget = 0.0;
      indot = vstream[0]*cellface[ncf].normal[0] +
	vstream[1]*cellface[ncf].normal[1] +
	vstream[2]*cellface[ncf].normal[2];
      if (indot >= 0.0) {
        for (isp = 0; isp < nspecies; isp++) {
          cellface[ncf].ntargetsp[isp] =
            mol_inflow(indot,vscale[isp],fraction[isp]);
          cellface[ncf].ntargetsp[isp] *= nrho*area*dt / fnum;
          cellface[ncf].ntargetsp[isp] /= cinfo[icell].weight;
          cellface[ncf].ntarget += cellface[ncf].ntargetsp[isp];
        }
        c2f[icell][XHI] = ncf++;
      } else c2f[icell][XHI] = -1;
    } else c2f[icell][XHI] = -1;

    if (c2f[icell][YLO]) {
      cellface[ncf].id = cells[icell].id;
      if (cells[icell].nsplit > 1) cellface[ncf].pcell = split(icell,YLO);
      else cellface[ncf].pcell = icell;
      cellface[ncf].icell = icell;
      cellface[ncf].iface = YLO;
      cellface[ncf].ndim = 1;
      cellface[ncf].pdim = 0;
      cellface[ncf].qdim = 2;
      cellface[ncf].lo[0] = cells[icell].lo[0];
      cellface[ncf].hi[0] = cells[icell].hi[0];
      cellface[ncf].lo[1] = cells[icell].lo[1];
      cellface[ncf].hi[1] = cells[icell].lo[1];
      cellface[ncf].lo[2] = cells[icell].lo[2];
      cellface[ncf].hi[2] = cells[icell].hi[2];
      cellface[ncf].normal[0] = 0.0;
      cellface[ncf].normal[1] = 1.0;
      cellface[ncf].normal[2] = 0.0;
      
      // no axi-symmetry allowed

      if (dimension == 3) 
        area = (cells[icell].hi[0]-cells[icell].lo[0]) * 
          (cells[icell].hi[2]-cells[icell].lo[2]);
      else if (dimension == 2) {
	cellface[ncf].lo[2] = 0.0;
	cellface[ncf].hi[2] = 0.0;
	area = (cells[icell].hi[0]-cells[icell].lo[0]);
      }
      
      cellface[ncf].ntarget = 0.0;
      indot = vstream[0]*cellface[ncf].normal[0] +
	vstream[1]*cellface[ncf].normal[1] +
	vstream[2]*cellface[ncf].normal[2];
      if (indot >= 0.0) {
        for (isp = 0; isp < nspecies; isp++) {
          cellface[ncf].ntargetsp[isp] =
            mol_inflow(indot,vscale[isp],fraction[isp]);
          cellface[ncf].ntargetsp[isp] *= nrho*area*dt / fnum;
          cellface[ncf].ntargetsp[isp] /= cinfo[icell].weight;
          cellface[ncf].ntarget += cellface[ncf].ntargetsp[isp];
        }
        c2f[icell][YLO] = ncf++;
      } else c2f[icell][YLO] = -1;
    } else c2f[icell][YLO] = -1;

    if (c2f[icell][YHI]) {
      cellface[ncf].id = cells[icell].id;
      if (cells[icell].nsplit > 1) cellface[ncf].pcell = split(icell,YHI);
      else cellface[ncf].pcell = icell;
      cellface[ncf].icell = icell;
      cellface[ncf].iface = YHI;
      cellface[ncf].ndim = 1;
      cellface[ncf].pdim = 0;
      cellface[ncf].qdim = 2;
      cellface[ncf].lo[0] = cells[icell].lo[0];
      cellface[ncf].hi[0] = cells[icell].hi[0];
      cellface[ncf].lo[1] = cells[icell].hi[1];
      cellface[ncf].hi[1] = cells[icell].hi[1];
      cellface[ncf].lo[2] = cells[icell].lo[2];
      cellface[ncf].hi[2] = cells[icell].hi[2];
      cellface[ncf].normal[0] = 0.0;
      cellface[ncf].normal[1] = -1.0;
      cellface[ncf].normal[2] = 0.0;
      
      if (dimension == 3)
        area = (cells[icell].hi[0]-cells[icell].lo[0]) * 
          (cells[icell].hi[2]-cells[icell].lo[2]);
      else if (domain->axisymmetric)
         area = (cells[icell].hi[0]*cells[icell].hi[0] -
                 cells[icell].lo[0]*cells[icell].lo[0])*MY_PI;
      else {
	cellface[ncf].lo[2] = 0.0;
	cellface[ncf].hi[2] = 0.0;
	area = (cells[icell].hi[0]-cells[icell].lo[0]);
      }
      
      cellface[ncf].ntarget = 0.0;
      indot = vstream[0]*cellface[ncf].normal[0] +
	vstream[1]*cellface[ncf].normal[1] +
	vstream[2]*cellface[ncf].normal[2];
      if (indot >= 0.0) {
        for (isp = 0; isp < nspecies; isp++) {
          cellface[ncf].ntargetsp[isp] =
            mol_inflow(indot,vscale[isp],fraction[isp]);
          cellface[ncf].ntargetsp[isp] *= nrho*area*dt / fnum;
          cellface[ncf].ntargetsp[isp] /= cinfo[icell].weight;
          cellface[ncf].ntarget += cellface[ncf].ntargetsp[isp];
        }
        c2f[icell][YHI] = ncf++;
      } else c2f[icell][YHI] = -1;
    } else c2f[icell][YHI] = -1;

    if (c2f[icell][ZLO]) {
      cellface[ncf].id = cells[icell].id;
      if (cells[icell].nsplit > 1) cellface[ncf].pcell = split(icell,ZLO);
      else cellface[ncf].pcell = icell;
      cellface[ncf].icell = icell;
      cellface[ncf].iface = ZLO;
      cellface[ncf].ndim = 2;
      cellface[ncf].pdim = 0;
      cellface[ncf].qdim = 1;
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
      if (indot >= 0.0) {
        for (isp = 0; isp < nspecies; isp++) {
          cellface[ncf].ntargetsp[isp] =
            mol_inflow(indot,vscale[isp],fraction[isp]);
          cellface[ncf].ntargetsp[isp] *= nrho*area*dt / fnum;
          cellface[ncf].ntargetsp[isp] /= cinfo[icell].weight;
          cellface[ncf].ntarget += cellface[ncf].ntargetsp[isp];
        }
        c2f[icell][ZLO] = ncf++;
      } else c2f[icell][ZLO] = -1;
    } else c2f[icell][ZLO] = -1;

    if (c2f[icell][ZHI]) {
      cellface[ncf].id = cells[icell].id;
      if (cells[icell].nsplit > 1) cellface[ncf].pcell = split(icell,ZHI);
      else cellface[ncf].pcell = icell;
      cellface[ncf].icell = icell;
      cellface[ncf].iface = ZHI;
      cellface[ncf].ndim = 2;
      cellface[ncf].pdim = 0;
      cellface[ncf].qdim = 1;
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
      if (indot >= 0.0) {
        for (isp = 0; isp < nspecies; isp++) {
          cellface[ncf].ntargetsp[isp] =
            mol_inflow(indot,vscale[isp],fraction[isp]);
          cellface[ncf].ntargetsp[isp] *= nrho*area*dt / fnum;
          cellface[ncf].ntargetsp[isp] /= cinfo[icell].weight;
          cellface[ncf].ntarget += cellface[ncf].ntargetsp[isp];
        }
        c2f[icell][ZHI] = ncf++;
      } else c2f[icell][ZHI] = -1;
    } else c2f[icell][ZHI] = -1;
  }


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

/* ----------------------------------------------------------------------
   insert particles in grid cells with faces touching inflow boundaries
   see Bird ???, eq 12.5, for generation of inflow velocities
------------------------------------------------------------------------- */

void FixInflow::start_of_step()
{
  int pcell,ninsert,isp,ispecies,ndim,pdim,qdim,id;
  double *lo,*hi,*normal;
  double x[3],v[3];
  double indot,scosine,rn,ntarget,vr;
  double beta_un,normalized_distbn_fn,theta,erot,evib;
  Particle::OnePart *p;

  if (update->ntimestep % nevery) return;

  int dimension = domain->dimension;
  Particle::OnePart *particles = particle->particles;
  Grid::ChildCell *cells = grid->cells;
  double dt = update->dt;

  int nspecies = particle->mixture[imix]->nspecies;
  int *species = particle->mixture[imix]->species;
  double *cummulative = particle->mixture[imix]->cummulative;
  double *vstream = particle->mixture[imix]->vstream;
  double *vscale = particle->mixture[imix]->vscale;
  double temp_thermal = particle->mixture[imix]->temp_thermal;

  // insert particles by cell/face pair
  // ntarget/ninsert is either perspecies or for all species
  // for one particle:
  //   x = random position on face
  //   v = randomized thermal velocity + vstream
  //       first stage: normal dimension (ndim)
  //       second stage: parallel dimensions (pdim,qdim)

  // double while loop until randomized particle velocity meets 2 criteria
  // inner do-while loop:
  //   v = vstream-component + vthermal is into simulation box
  //   see Bird 1994, p 425
  // outer do-while loop:
  //   shift Maxwellian distribution by stream velocity component
  //   see Bird 1994, p 259, eq 12.5

  // NOTE:
  // currently not allowing particle insertion on backflow boundaries
  // enforced by indot >= 0.0 check in init()
  // could allow particle insertion on backflow boundaries
  //   when streaming velocity is small enough
  // need to insure two do-while loops below do not spin endlessly

  int nfix_add_particle = modify->n_add_particle;
  nsingle = 0;

  for (int i = 0; i < ncf; i++) {
    pcell = cellface[i].pcell;
    ndim = cellface[i].ndim;
    pdim = cellface[i].pdim;
    qdim = cellface[i].qdim;
    lo = cellface[i].lo;
    hi = cellface[i].hi;
    normal = cellface[i].normal;

    indot = vstream[0]*normal[0] + vstream[1]*normal[1] + vstream[2]*normal[2];

    if (perspecies == YES) {
      for (isp = 0; isp < nspecies; isp++) {
        ispecies = species[isp];
	ntarget = cellface[i].ntargetsp[isp]+random->uniform();
	ninsert = static_cast<int> (ntarget);
        scosine = indot / vscale[isp];

	for (int m = 0; m < ninsert; m++) {
	  x[0] = lo[0] + random->uniform() * (hi[0]-lo[0]);
	  x[1] = lo[1] + random->uniform() * (hi[1]-lo[1]);
	  if (dimension == 3) x[2] = lo[2] + random->uniform() * (hi[2]-lo[2]);
          else x[2] = 0.0;

	  do {
	    do beta_un = (6.0*random->gaussian() - 3.0);
	    while (beta_un + scosine < 0.0);
	    normalized_distbn_fn = 2.0 * (beta_un + scosine) / 
	      (scosine + sqrt(scosine*scosine + 2.0)) *
	      exp(0.5 + (0.5*scosine)*(scosine-sqrt(scosine*scosine + 2.0)) - 
		  beta_un*beta_un);
	  } while (normalized_distbn_fn < random->uniform());
	  
          v[ndim] = beta_un*vscale[isp]*normal[ndim] + vstream[ndim];

          theta = MY_2PI * random->gaussian();
          vr = vscale[isp] * sqrt(-log(random->uniform()));
          v[pdim] = vr * sin(theta) + vstream[pdim];
          v[qdim] = vr * cos(theta) + vstream[qdim];
          erot = particle->erot(ispecies,temp_thermal,random);
          evib = particle->evib(ispecies,temp_thermal,random);
          id = MAXSMALLINT*random->uniform();
	  particle->add_particle(id,ispecies,pcell,x,v,erot,evib);

          p = &particle->particles[particle->nlocal-1];
          p->flag = PINSERT;
          p->dtremain = dt * random->uniform();

          if (nfix_add_particle) 
            modify->add_particle(particle->nlocal-1,temp_thermal,vstream);
	}

	nsingle += ninsert;
      }

    } else {
      if (np == 0) {
	ntarget = cellface[i].ntarget+random->uniform();
	ninsert = static_cast<int> (ntarget);
      } else {
	ninsert = npercell;
	if (i >= nthresh) ninsert++;
      }

      for (int m = 0; m < ninsert; m++) {
	rn = random->uniform();
	isp = 0;
	while (cummulative[isp] < rn) isp++;
        ispecies = species[isp];
        scosine = indot / vscale[isp];

	x[0] = lo[0] + random->uniform() * (hi[0]-lo[0]);
	x[1] = lo[1] + random->uniform() * (hi[1]-lo[1]);
        if (dimension == 3) x[2] = lo[2] + random->uniform() * (hi[2]-lo[2]);
        else x[2] = 0.0;

	do {
	  do {
	    beta_un = (6.0*random->gaussian() - 3.0);
	  } while (beta_un + scosine < 0.0);
	  normalized_distbn_fn = 2.0 * (beta_un + scosine) / 
	    (scosine + sqrt(scosine*scosine + 2.0)) *
	    exp(0.5 + (0.5*scosine)*(scosine-sqrt(scosine*scosine + 2.0)) - 
		beta_un*beta_un);
	} while (normalized_distbn_fn < random->uniform());
	
        v[ndim] = beta_un*vscale[isp]*normal[ndim] + vstream[ndim];

        theta = MY_2PI * random->gaussian();
        vr = vscale[isp] * sqrt(-log(random->uniform()));
        v[pdim] = vr * sin(theta) + vstream[pdim];
        v[qdim] = vr * cos(theta) + vstream[qdim];
        erot = particle->erot(ispecies,temp_thermal,random);
        evib = particle->evib(ispecies,temp_thermal,random);
        id = MAXSMALLINT*random->uniform();
	particle->add_particle(id,ispecies,pcell,x,v,erot,evib);

        p = &particle->particles[particle->nlocal-1];
        p->flag = PINSERT;
        p->dtremain = dt * random->uniform();

        if (nfix_add_particle) 
          modify->add_particle(particle->nlocal-1,temp_thermal,vstream);
      }

      nsingle += ninsert;
    }
  }

  ntotal += nsingle;
}

/* ----------------------------------------------------------------------
   calculate flux of particles of a species with vscale/fraction
     entering a grid cell
   indot = vstream dotted into face normal, assumed to be >= 0.0
   scosine = s cos(theta) in Bird notation where vscale = 1/beta
   see Bird 1994, eq 4.22
------------------------------------------------------------------------- */

double FixInflow::mol_inflow(double indot, double vscale, double fraction)
{
  double scosine = indot / vscale;
  double inward_number_flux = vscale*fraction *
    (exp(-scosine*scosine) + sqrt(MY_PI)*scosine*(1.0 + erf(scosine))) / 
    (2*sqrt(MY_PI));
  return inward_number_flux;
}

/* ----------------------------------------------------------------------
   inserting into split cell parent ICELL on face FLAG
   determine which child split cell the face is part of
   face cannot be touched by surfs, so entire face is part of one split cell
   compute which via update->split() and return it
------------------------------------------------------------------------- */

int FixInflow::split(int icell, int flag)
{
  double x[3];

  Grid::ChildCell *cells = grid->cells;

  // x = center point on face

  x[0] = 0.5 * (cells[icell].lo[0] + cells[icell].hi[0]);
  x[1] = 0.5 * (cells[icell].lo[1] + cells[icell].hi[1]);
  x[2] = 0.5 * (cells[icell].lo[2] + cells[icell].hi[2]);
  if (flag == XLO) x[0] = cells[icell].lo[0];
  if (flag == XHI) x[0] = cells[icell].hi[0];
  if (flag == YLO) x[1] = cells[icell].lo[1];
  if (flag == YHI) x[1] = cells[icell].hi[1];
  if (flag == ZLO) x[2] = cells[icell].lo[2];
  if (flag == ZHI) x[2] = cells[icell].hi[2];
  if (domain->dimension == 2) x[2] = 0.0;

  int splitcell;
  if (domain->dimension == 2) splitcell = update->split2d(icell,x);
  else splitcell = update->split3d(icell,x);
  return splitcell;
}

/* ----------------------------------------------------------------------
   pack icell values for per-cell arrays into buf
   also pack cellface data for flagged faces
   return byte count of amount packed
   if memflag, only return count, do not fill buf
------------------------------------------------------------------------- */

int FixInflow::pack_grid_one(int icell, char *buf, int memflag)
{
  char *ptr = buf;

  int nspecies = particle->mixture[imix]->nspecies;

  if (memflag) memcpy(ptr,c2f[icell],6*sizeof(int));
  ptr += 6*sizeof(int);
  ptr = ROUNDUP(ptr);

  // loop over faces, pack cellface entry if it exists plus its vectors

  for (int iface = 0; iface < 6; iface++) {
    if (c2f[icell][iface] < 0) continue;
    int icf = c2f[icell][iface];
    if (memflag) memcpy(ptr,&cellface[icf],sizeof(CellFace));
    ptr += sizeof(CellFace);
    ptr = ROUNDUP(ptr);
    if (memflag) memcpy(ptr,cellface[icf].ntargetsp,nspecies*sizeof(double));
    ptr += nspecies*sizeof(double);
  }

  return ptr-buf;
}

/* ----------------------------------------------------------------------
   unpack icell values for per-cell arrays from buf
   also unpack cellface data for flagged faces
   return byte count of amount unpacked
------------------------------------------------------------------------- */

int FixInflow::unpack_grid_one(int icell, char *buf)
{
  char *ptr = buf;

  int nspecies = particle->mixture[imix]->nspecies;

  grow_percell(1);
  memcpy(c2f[icell],ptr,6*sizeof(int));
  ptr += 6*sizeof(int);
  ptr = ROUNDUP(ptr);
  nglocal++;

  // unpack c2f values
  // fill sub cells with -1

  int nsplit = grid->cells[icell].nsplit;
  if (nsplit > 1) {
    grow_percell(nsplit);
    for (int i = 0; i < nsplit; i++) {
      c2f[nglocal][0] = c2f[nglocal][1] = c2f[nglocal][2] = 
        c2f[nglocal][3] = c2f[nglocal][4] = c2f[nglocal][5] = -1;
      nglocal++;
    }
  }
  
  // unpack cellface entry for each face set in c2f
  // store ntargetsp to avoid overwriting allocated cellface.ntargetsp
  // reset c2f pointer into cellface
  // reset cellface.pcell and cellface.icell
  // pcell setting based on sub cells being immediately after the split cell

  for (int iface = 0; iface < 6; iface++) {
    if (c2f[icell][iface] < 0) continue;
    c2f[icell][iface] = ncf;
    grow_cellface(1);
    double *ntargetsp = cellface[ncf].ntargetsp;
    memcpy(&cellface[ncf],ptr,sizeof(CellFace));
    ptr += sizeof(CellFace);
    ptr = ROUNDUP(ptr);
    cellface[ncf].ntargetsp = ntargetsp;
    memcpy(cellface[ncf].ntargetsp,ptr,nspecies*sizeof(double));
    ptr += nspecies*sizeof(double);
    if (grid->cells[icell].nsplit == 1) cellface[ncf].pcell = icell;
    else cellface[ncf].pcell = split(icell,cellface[ncf].iface);
    cellface[ncf].icell = icell;
    ncf++;
  }

  return ptr-buf;
}

/* ----------------------------------------------------------------------
   compress per-cell arrays due to cells migrating to new procs
   criteria for keeping/discarding a cell is same as in Grid::compress()
   this keeps final ordering of per-cell arrays consistent with Grid class
------------------------------------------------------------------------- */

void FixInflow::compress_grid()
{
  int me = comm->me;
  Grid::ChildCell *cells = grid->cells;

  // compress cellface first, preserving/copying ntargetsp vector
  // reset c2f pointers to new cellface indices

  int nspecies = particle->mixture[imix]->nspecies;
  double *ntargetsp;

  int ncurrent = ncf;
  ncf = 0;
  for (int icf = 0; icf < ncurrent; icf++) {
    int icell = cellface[icf].icell;
    if (cells[icell].proc != me) continue;
    if (ncf != icf) {
      ntargetsp = cellface[ncf].ntargetsp;
      memcpy(&cellface[ncf],&cellface[icf],sizeof(CellFace));
      cellface[ncf].ntargetsp = ntargetsp;
      memcpy(ntargetsp,cellface[icf].ntargetsp,nspecies*sizeof(double));
    }
    c2f[icell][cellface[icf].iface] = ncf;
    ncf++;
  }

  // compress c2f second
  // keep an unsplit or split cell if staying on this proc
  // keep a sub cell if its split cell is staying on this proc
  // reset cellface.icell to new cell index
  // cellface.pcell reset in post_compress(),
  //   since don't know sub cell index at this point

  int nbytes = 6*sizeof(int);

  ncurrent = nglocal;
  nglocal = 0;
  for (int icell = 0; icell < ncurrent; icell++) {
    if (cells[icell].nsplit >= 1) {
      if (cells[icell].proc != me) continue;
    } else {
      int isplit = cells[icell].isplit;
      if (cells[grid->sinfo[isplit].icell].proc != me) continue;
    }

    if (nglocal != icell) memcpy(c2f[nglocal],c2f[icell],nbytes);

    if (c2f[nglocal][XLO] >= 0) cellface[c2f[nglocal][XLO]].icell = nglocal;
    if (c2f[nglocal][XHI] >= 0) cellface[c2f[nglocal][XHI]].icell = nglocal;
    if (c2f[nglocal][YLO] >= 0) cellface[c2f[nglocal][YLO]].icell = nglocal;
    if (c2f[nglocal][YHI] >= 0) cellface[c2f[nglocal][YHI]].icell = nglocal;
    if (c2f[nglocal][ZLO] >= 0) cellface[c2f[nglocal][ZLO]].icell = nglocal;
    if (c2f[nglocal][ZHI] >= 0) cellface[c2f[nglocal][ZHI]].icell = nglocal;

    nglocal++;
  }
}

/* ----------------------------------------------------------------------
   reset cellface.pcell for compressed cellface entries
   called from Grid::compress() after grid cells have been compressed
------------------------------------------------------------------------- */

void FixInflow::post_compress_grid()
{
  Grid::ChildCell *cells = grid->cells;

  for (int icf = 0; icf < ncf; icf++) {
    int icell = cellface[icf].icell;
    int nsplit = cells[icell].nsplit;
    if (nsplit == 1) cellface[icf].pcell = icell;
    else cellface[icf].pcell = split(icell,cellface[icf].iface);
  }
}

/* ----------------------------------------------------------------------
   insure c2f allocated long enough for N new cells
------------------------------------------------------------------------- */

void FixInflow::grow_percell(int n)
{
  if (nglocal+n < nglocalmax) return;
  nglocalmax += DELTAGRID;
  memory->grow(c2f,nglocalmax,6,"inflow:c2f");
}

/* ----------------------------------------------------------------------
   insure cellface allocated long enough for N new cell/face pairs
   also allocate new vector within each cellface
------------------------------------------------------------------------- */

void FixInflow::grow_cellface(int n)
{
  if (ncf+n < ncfmax) return;
  int oldmax = ncfmax;
  ncfmax += DELTAFACE;
  cellface = (CellFace *) memory->srealloc(cellface,ncfmax*sizeof(CellFace),
                                           "inflow:cellface");
  memset(&cellface[oldmax],0,(ncfmax-oldmax)*sizeof(CellFace));

  int nspecies = particle->mixture[imix]->nspecies;
  for (int i = oldmax; i < ncfmax; i++)
    cellface[i].ntargetsp = new double[nspecies];
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

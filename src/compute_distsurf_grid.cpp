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

#include "string.h"
#include "compute_distsurf_grid.h"
#include "update.h"
#include "grid.h"
#include "surf.h"
#include "domain.h"
#include "input.h"
#include "geometry.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

ComputeDistSurfGrid::
ComputeDistSurfGrid(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal compute distsurf/grid command");

  per_grid_flag = 1;
  size_per_grid_cols = 0;

  int igroup = surf->find_group(arg[2]);
  if (igroup < 0)
    error->all(FLERR,"Compute distsurf/grid command surface group "
               "does not exist");
  sgroupbit = surf->bitmask[igroup];

  // optional args

  sdir[0] = sdir[1] = sdir[2] = 0.0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"dir") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal distsurf/grid command");
      sdir[0] = input->numeric(FLERR,arg[iarg+1]);
      sdir[1] = input->numeric(FLERR,arg[iarg+2]);
      sdir[2] = input->numeric(FLERR,arg[iarg+3]);
      if (domain->dimension == 2 && sdir[2] != 0.0) 
        error->all(FLERR,"Illegal adapt command");
      iarg += 4;
    } else error->all(FLERR,"Illegal compute distsurf/grid command");
  }

  // setup per grid vector

  nglocal = 0;
  vector_grid = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeDistSurfGrid::~ComputeDistSurfGrid()
{
  memory->destroy(vector_grid);
}

/* ---------------------------------------------------------------------- */

void ComputeDistSurfGrid::init()
{
  reallocate();
}

/* ---------------------------------------------------------------------- */

void ComputeDistSurfGrid::compute_per_grid()
{
  int i,m,n,p1,p2,p3;
  int *csurfs,*csubs;
  double dist,mindist;
  double *lo,*hi;
  double cctr[3],cell2surf[3];

  invoked_per_grid = update->ntimestep;

  int dim = domain->dimension;
  Surf::Point *pts = surf->pts;
  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  int nline = surf->nline;
  int ntri = surf->ntri;

  // eflag = 0/1 flag = eligibility of each surf, based on normal vs sdir
  // slist = list of Nsurf eligible surf indices to use for computing distances

  int *eflag,*slist;
  int nsurf = 0;

  if (dim == 2) {
    memory->create(eflag,nline,"distsurf/grid:eflag");
    memory->create(slist,nline,"distsurf/grid:slist");
    for (i = 0; i < nline; i++) {
      eflag[i] = 0;
      if (!(lines[i].mask & sgroupbit)) continue;
      if (MathExtra::dot3(lines[i].norm,sdir) <= 0.0) {
        eflag[i] = 1;
        slist[nsurf++] = i;
      }
    }
  } else {
    memory->create(eflag,nline,"distsurf/grid:eflag");
    memory->create(slist,ntri,"distsurf/grid:slist");
    for (i = 0; i < ntri; i++) {
      eflag[i] = 0;
      if (!(tris[i].mask & sgroupbit)) continue;
      if (MathExtra::dot3(tris[i].norm,sdir) <= 0.0) {
        eflag[i] = 1;
        slist[nsurf++] = i;
      }
    }
  }

  // pre-compute center point of each eligible surf

  double invthird = 1.0/3.0;

  double **sctr;
  memory->create(sctr,nsurf,3,"distsurf/grid:sctr");
  for (i = 0; i < nsurf; i++) {
    m = slist[i];
    if (dim == 2) {
      if (!(lines[i].mask & sgroupbit)) continue;
      p1 = lines[m].p1;
      p2 = lines[m].p2;
      sctr[i][0] = 0.5 * (pts[p1].x[0] + pts[p2].x[0]);
      sctr[i][1] = 0.5 * (pts[p1].x[1] + pts[p2].x[1]);
      sctr[i][2] = 0.0;
    } else {
      if (!(tris[i].mask & sgroupbit)) continue;
      p1 = tris[m].p1;
      p2 = tris[m].p2;
      p2 = tris[m].p3;
      sctr[i][0] = invthird * (pts[p1].x[0] + pts[p2].x[0] + pts[p3].x[0]);
      sctr[i][1] = invthird * (pts[p1].x[1] + pts[p2].x[1] + pts[p3].x[1]);
      sctr[i][2] = invthird * (pts[p1].x[2] + pts[p2].x[2] + pts[p3].x[2]);
    }
  }

  // loop over my unsplit/split grid cells
  // if surfs in cell and any are eligible, dist = 0.0
  // else loop over eligible surfs in slist:
  //   if vector from cell center to surf center is against surf norm, skip it
  //   compute distance from cell to surf via Geometry method
  //   dist = minimum distance to any eligible surf
  // if assign dist to split cell, also assign dist to all its sub cells

  Grid::ChildCell *cells = grid->cells;
  Grid::SplitInfo *sinfo = grid->sinfo;

  for (int icell = 0; icell < nglocal; icell++) {
    if (cells[icell].nsplit < 1) continue;

    if (cells[icell].nsurf) {
      n = cells[icell].nsurf;
      csurfs = cells[icell].csurfs;
      for (i = 0; i < n; i++) {
        m = csurfs[i];
        if (eflag[m]) break;
      }
      if (i < n) {
        vector_grid[icell] = 0.0;
        continue;
      }
    }

    lo = cells[icell].lo;
    hi = cells[icell].hi;
    cctr[0] = 0.5 * (lo[0]+hi[0]);
    cctr[1] = 0.5 * (lo[1]+hi[1]);
    if (dim == 3) cctr[2] = 0.5 * (lo[2]+hi[2]);
    else cctr[2] = 0.0;

    mindist = BIG;
    for (i = 0; i < nsurf; i++) {
      m = slist[i];
      
      cell2surf[0] = sctr[i][0] - cctr[0];
      cell2surf[1] = sctr[i][1] - cctr[1];
      cell2surf[2] = sctr[i][2] - cctr[2];

      if (dim == 2) {
        if (MathExtra::dot3(cell2surf,lines[m].norm) > 0.0) continue;
        dist = Geometry::dist_line_quad(pts[lines[m].p1].x,pts[lines[m].p2].x,
                                        lo,hi);
      } else {
        if (MathExtra::dot3(cell2surf,tris[m].norm) > 0.0) continue;
        dist = Geometry::dist_tri_hex(pts[tris[m].p1].x,pts[tris[m].p2].x,
                                      pts[tris[m].p3].x,tris[m].norm,lo,hi);
      }

      mindist = MIN(mindist,dist);
    }
    
    vector_grid[icell] = mindist;

    // also set vector for sub-cells of split cell

    if (cells[icell].nsplit > 1) {
      n = cells[icell].nsplit;
      csubs = sinfo[cells[icell].isplit].csubs;
      for (i = 0; i < n; i++) {
        m = csubs[i];
        vector_grid[m] = mindist;
      }
    }
  }

  // clean up

  memory->destroy(eflag);
  memory->destroy(slist);
  memory->destroy(sctr);
}

/* ----------------------------------------------------------------------
   reallocate vector if nglocal has changed
   called by init() and load balancer
------------------------------------------------------------------------- */

void ComputeDistSurfGrid::reallocate()
{
  if (grid->nlocal == nglocal) return;

  nglocal = grid->nlocal;
  memory->destroy(vector_grid);
  memory->create(vector_grid,nglocal,"distsurf/grid:vector_grid");
}

/* ----------------------------------------------------------------------
   memory usage of local grid-based array
------------------------------------------------------------------------- */

bigint ComputeDistSurfGrid::memory_usage()
{
  bigint bytes;
  bytes = nglocal * sizeof(double);
  return bytes;
}

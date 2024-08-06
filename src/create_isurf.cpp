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

/* ----------------------------------------------------------------------
   Contributing author: Andrew Hong (Sandia)
------------------------------------------------------------------------- */

#include "comm.h"
#include "create_isurf.h"
#include "domain.h"
#include "error.h"
#include "fix_ablate.h"
#include "geometry.h"
#include "grid.h"
#include "input.h"
#include "math_extra.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "mpi.h"
#include "region.h"
#include "surf.h"

using namespace SPARTA_NS;
using namespace MathConst;

#define DELTAGRID 1024            // must be bigger than split cells per cell
#define DELTASEND 1024
#define EPSILON_GRID 1.0e-3

enum{CVAL,SVAL,IVAL,INVAL};
enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};         // same as Domain
enum{NCHILD,NPARENT,NUNKNOWN,NPBCHILD,NPBPARENT,NPBUNKNOWN,NBOUND};  // Update
enum{INOUT,VOXEL,AVE,INNER}; // how to mark corner values

/* ---------------------------------------------------------------------- */

CreateISurf::CreateISurf(SPARTA *sparta) : Pointers(sparta)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  dim = domain->dimension;

  if (dim == 2) {
    ncorner = 4;
    ninner = 4;
    nedge = 4;
  } else {
    ncorner = 8;
    ninner = 6;
    nedge = 12;
  }

  // for finding corner values

  cvalues = NULL;
  invalues = NULL;
  tmp_cvalues = NULL;
  tmp_invalues = NULL;
  svalues = NULL;
  mvalues = NULL;
  ivalues = NULL;

  // for communication

  ixyz = NULL;
  cghost = NULL;
  inghost = NULL;
  sghost = NULL;
  ighost = NULL;
  numsend = NULL;
  maxgrid = maxghost = 0;

  proclist = NULL;
  locallist = NULL;
  maxsend = 0;

  sbuf = NULL;
  maxsbuf = 0;
}

/* ---------------------------------------------------------------------- */

CreateISurf::~CreateISurf()
{
  memory->destroy(cvalues);
  memory->destroy(invalues);
  memory->destroy(tmp_cvalues);
  memory->destroy(tmp_invalues);
  memory->destroy(svalues);
  memory->destroy(mvalues);
  memory->destroy(ivalues);

  memory->destroy(ixyz);
  memory->destroy(cghost);
  memory->destroy(inghost);
  memory->destroy(sghost);
  memory->destroy(ighost);
  memory->destroy(numsend);

  memory->destroy(proclist);
  memory->destroy(locallist);

  memory->destroy(sbuf);
}

/* ---------------------------------------------------------------------- */

void CreateISurf::command(int narg, char **arg)
{
  if (!grid->exist)
    error->all(FLERR,"Cannot create_isurf before grid is defined");
  if (!surf->exist)
    error->all(FLERR,"Cannot create_isurf unless explicit surfaces exist");
  if (surf->implicit && surf->exist)
    error->all(FLERR,"Create_isurf Cannot have pre-existing implicit surfaces");
  if (!surf->distributed)
    error->all(FLERR,"Create_isurf requires distributed explicit surfaces");
  if (domain->axisymmetric)
    error->all(FLERR,"Cannot create_isurf for axisymmetric domains");

  if (narg < 4) error->all(FLERR,"Illegal create_isurf command");

  // grid group

  ggroup = grid->find_group(arg[0]);
  if (ggroup < 0) error->all(FLERR,"Create_isurf grid group ID does not exist");
  groupbit = grid->bitmask[ggroup];

  // ablate fix

  char *ablateID = arg[1];
  int ifix = modify->find_fix(ablateID);
  if (ifix < 0)
    error->all(FLERR,"Fix ID for create_isurf does not exist");
  if (strcmp(modify->fix[ifix]->style,"ablate") != 0)
    error->all(FLERR,"Fix for create_isurf is not a fix ablate");
  ablate = (FixAblate *) modify->fix[ifix];
  if (ggroup != ablate->igroup)
    error->all(FLERR,"Create_isurf group does not match fix ablate group");

  // threshold for corner values

  thresh = input->numeric(FLERR,arg[2]);
  if (thresh < 0 || thresh > 255)
    error->all(FLERR,"Create_isurf thresh must be bounded as (0,255)");
  int ithresh = static_cast<int> (thresh);
  if (ithresh == thresh)
    error->all(FLERR,"An integer value for create_isurf thresh is not allowed");

  // mode to determine corner values

  if (strcmp(arg[3],"inout") == 0) ctype = INOUT;
  else if (strcmp(arg[3],"voxel") == 0) ctype = VOXEL;
  else if (strcmp(arg[3],"ave") == 0) ctype = AVE;
  else if (strcmp(arg[3],"inner") == 0) ctype = INNER;
  else error->all(FLERR,"Create_isurf corner mode is invalid");

  // process optional command line args (nothing to process for now)

  //process_args(narg-4,&arg[4]);

  // check if grid group is a uniform grid

  grid->check_uniform_group(ggroup,nxyz,corner,xyzsize);
  nglocal = grid->nlocal;

  // set grid corner point values based on existing explicit surfs and thresh

  if (ctype == INNER) set_inner();
  else set_corners();

  // remove all explicit surfs

  remove_old();

  // update flags

  surf->implicit = 1;
  surf->exist = 1;

  // store corner point values in fix ablate instance
  // this call will also create implicit surfs
  // set pushflag = 0 so averaging option is not overridden

  tvalues = NULL;
  int pushflag = 0;
  char *sgroupID = NULL;
  if (ctype == INNER) {
    ablate->store_corners(nxyz[0],nxyz[1],nxyz[2],corner,xyzsize,
                          invalues,tvalues,thresh,sgroupID,pushflag);
  } else {
    ablate->store_corners(nxyz[0],nxyz[1],nxyz[2],corner,xyzsize,
                          cvalues,tvalues,thresh,sgroupID,pushflag);
  }

  if (ablate->nevery == 0) modify->delete_fix(ablateID);
}

/* ----------------------------------------------------------------------
   compute corner values from explicit surfs and thresh parameter
   cvalues used later in ablate to create the implicit surfaces
------------------------------------------------------------------------- */

void CreateISurf::set_corners()
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  // set cell spacing and indices

  memory->grow(ixyz,nglocal,3,"create_isurf:ixyz");

  for (int icell = 0; icell < nglocal; icell++) {
    ixyz[icell][0] = ixyz[icell][1] = ixyz[icell][2] = 0;
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    ixyz[icell][0] =
      static_cast<int> ((cells[icell].lo[0]-corner[0]) / xyzsize[0] + 0.5) + 1;
    ixyz[icell][1] =
      static_cast<int> ((cells[icell].lo[1]-corner[1]) / xyzsize[1] + 0.5) + 1;
    ixyz[icell][2] =
      static_cast<int> ((cells[icell].lo[2]-corner[2]) / xyzsize[2] + 0.5) + 1;
  }

  // first shift everything down by thresh
  // later shift back

  memory->create(cvalues,nglocal,ncorner,"createisurf:cvalues");
  memory->create(tmp_cvalues,nglocal,ncorner,"createisurf:tmp_cvalues");
  memory->create(mvalues,nglocal,ncorner,"createisurf:mvalues");
  memory->create(svalues,nglocal,ncorner,"createisurf:svalues");
  memory->create(ivalues,nglocal,ncorner,ninner,"createisurf:ivalues");

  // initialize

  for (int ic = 0; ic < nglocal; ic++) {
    for (int jc = 0; jc < ncorner; jc++) {
      tmp_cvalues[ic][jc] = -1.0;
      svalues[ic][jc] = -1;
      mvalues[ic][jc] = -1.0;
      for (int kc = 0; kc < ninner; kc++) ivalues[ic][jc][kc] = -1.0;
    }
  }


  // find intersections between edges and surfaces

  if (dim == 2) surface_edge2d();
  else surface_edge3d();

  // fill in side and corner values based on if grid cell is in or out

  set_inout();

  // sync side and intersection values between procs

  sync(SVAL);
  sync(IVAL);

  // find remaining side values based on neighbors

  int full;
  if (dim == 2) full = find_side_2d();
  else full = find_side_3d();

  // check all procs successfully filled side values

  int ofull;
  if(full) ofull = 0;
  else ofull = 1;
  int allofull;
  MPI_Allreduce(&ofull,&allofull,1,MPI_INT,MPI_SUM,world);
  if (allofull) {
    char str[128];
    sprintf(str,
            "Create_isurf could not determine whether some corner \
             values are inside or outside with respect to the surface");
    error->all(FLERR,str);
  }

  // from side (and intersection) values, determine corner point values

  set_cvalues();

  // sync corner point values between procs

  if (ctype == VOXEL) sync_voxels();
  else sync(CVAL);

  return;
}

/* ----------------------------------------------------------------------
   same as above but for inner values
------------------------------------------------------------------------- */

void CreateISurf::set_inner()
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  // set cell spacing and indices

  memory->grow(ixyz,nglocal,3,"create_isurf:ixyz");
  for (int icell = 0; icell < nglocal; icell++) {
    ixyz[icell][0] = ixyz[icell][1] = ixyz[icell][2] = 0;
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    ixyz[icell][0] =
      static_cast<int> ((cells[icell].lo[0]-corner[0]) / xyzsize[0] + 0.5) + 1;
    ixyz[icell][1] =
      static_cast<int> ((cells[icell].lo[1]-corner[1]) / xyzsize[1] + 0.5) + 1;
    ixyz[icell][2] =
      static_cast<int> ((cells[icell].lo[2]-corner[2]) / xyzsize[2] + 0.5) + 1;
  }

  // first shift everything down by thresh
  // later shift back

  memory->create(invalues,nglocal,ncorner,ninner,"createisurf:invalues");
  memory->create(tmp_invalues,nglocal,ncorner,ninner,"createisurf:tmp_invalues");
  memory->create(mvalues,nglocal,ncorner,"createisurf:mvalues");
  memory->create(svalues,nglocal,ncorner,"createisurf:svalues");
  memory->create(ivalues,nglocal,ncorner,ninner,"createisurf:ivalues");

  // initialize

  for (int ic = 0; ic < nglocal; ic++) {
    for (int jc = 0; jc < ncorner; jc++) {
      svalues[ic][jc] = -1;
      mvalues[ic][jc] = -1.0;
      for (int kc = 0; kc < ninner; kc++) ivalues[ic][jc][kc] = -1.0;
      for (int kc = 0; kc < ninner; kc++) tmp_invalues[ic][jc][kc] = -1.0;
    }
  }

  // find intersections between edges and surfaces

  if (dim == 2) surface_edge2d();
  else surface_edge3d();

  // fill in side and corner values based on if grid cell is in or out

  set_inout();

  // sync side and intersection values between procs

  sync(SVAL);
  sync(IVAL);

  // find remaining side values based on neighbors

  int full;
  if (dim == 2) full = find_side_2d();
  else full = find_side_3d();

  // check all procs successfully filled side values

  int ofull;
  if (full) ofull = 0;
  else ofull = 1;
  int allofull;
  MPI_Allreduce(&ofull,&allofull,1,MPI_INT,MPI_SUM,world);
  if (allofull) {
    char str[128];
    sprintf(str,
            "Create_isurf could not determine whether some corner \
             values are inside or outside with respect to the surface");
    error->all(FLERR,str);
  }

  // from side (and intersection) values, determine corner point values

  set_cvalues();

  // sync corner point values between procs

  sync(INVAL);

}

/* ----------------------------------------------------------------------
   finds intersections and sides for corners in cells with surfaces
   2d version
------------------------------------------------------------------------- */

void CreateISurf::surface_edge2d()
{
  int i,j,n1,n2;

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  double cl[3], ch[3]; // cell bounds
  double cx[8], cy[8], cz[8]; // store all corners
  double pi[3], pj[3]; // corner coordinate
  Surf::Line *line;
  Surf::Line *lines = surf->lines;
  surfint *csurfs;
  int isurf, nsurf, side, hitflag;
  double param, oparam;

  // initialize param to avoid error

  param = 0.0;

  // indices for corners making up the cell edges

  int ci[4], cj[4];
  ci[0] = 0; cj[0] = 1;
  ci[1] = 1; cj[1] = 3;
  ci[2] = 3; cj[2] = 2;
  ci[3] = 2; cj[3] = 0;

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    // if no surfs, continue

    nsurf = cells[icell].nsurf;
    if (!nsurf) continue;

    for (int d = 0; d < 3; d++) {
      cl[d] = cells[icell].lo[d];
      ch[d] = cells[icell].hi[d];
    }

    // store cell corners

    cx[0] = cx[2] = cl[0];
    cx[1] = cx[3] = ch[0];

    cy[0] = cy[1] = cl[1];
    cy[2] = cy[3] = ch[1];

    // zero out z-direction

    cz[0] = cz[1] = cz[2] = cz[3] = 0.0;

    // smallest cell length
    // should only have to do this once

    mind = MIN(ch[0]-cl[0], ch[1]-cl[1]);

    // iterate over all edges

    csurfs = cells[icell].csurfs;

    for (int ic = 0; ic < nedge; ic++) {
      i = ci[ic];
      pi[0] = cx[i];
      pi[1] = cy[i];
      pi[2] = cz[i];

      j = cj[ic];
      pj[0] = cx[j];
      pj[1] = cy[j];
      pj[2] = cz[j];

      // test all surfs+corners to see if any hit

      for (int m = 0; m < nsurf; m++) {
        isurf = csurfs[m];
        line = &lines[isurf];

        // always start with lower corner
        // side will always be either 0 or 1

        hitflag = corner_hit2d(pi, pj, line, param, side);

        // need to take care of values near machine precision

        if (param < EPSILON_GRID*mind) param = 0.0;
        if ((1.0-param) < EPSILON_GRID*mind) param = 1.0;
        oparam = 1.0-param;

        // once a hit is found

        if (hitflag) {
          if (ic == 0) {
            n1 = 1;
            n2 = 0;
          } else if (ic == 1) {
            n1 = 3;
            n2 = 2;
          } else if (ic == 2) {
            n1 = 0;
            n2 = 1;
          } else {
            n1 = 2;
            n2 = 3;
          }

          if (ivalues[icell][i][n1] > param || ivalues[icell][i][n1] < 0)
            ivalues[icell][i][n1] = param;
          if (ivalues[icell][j][n2] > oparam || ivalues[icell][j][n2] < 0)
            ivalues[icell][j][n2] = oparam;


          if ((mvalues[icell][i] < 0 || param <= mvalues[icell][i]) && 
               svalues[icell][i] != 2) {
            if (param == 0) svalues[icell][i] = 0;

            // conflicting sides from two surfaces meeting at corner

            else if (fabs(mvalues[icell][i]-param) < EPSILON_GRID
              && svalues[icell][i] != side) svalues[icell][i] = 2;
            else svalues[icell][i] = side;

            mvalues[icell][i] = param;
          }

          if ((mvalues[icell][j] < 0 || oparam <= mvalues[icell][j]) && 
               svalues[icell][j] != 2) {
            if (oparam == 0) svalues[icell][j] = 0;
            else if (fabs(mvalues[icell][j]-oparam) < EPSILON_GRID
              && svalues[icell][j] != !side) svalues[icell][j] = 2;
            else svalues[icell][j] = !side;

            mvalues[icell][j] = oparam;
          }
        } // end "if" hitflag
      } // end "for" for surfaces
    } // end "for" for corners in cell
  } // end "for" for grid cells

  return;
}

/* ----------------------------------------------------------------------
   finds intersections and sides for corners in cells with surfaces
   3d version
------------------------------------------------------------------------- */

void CreateISurf::surface_edge3d()
{
  int i,j,n1,n2;

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  double cl[3], ch[3]; // cell bounds
  double cx[8], cy[8], cz[8]; // store all corners
  double pi[3], pj[3]; // corner coordinate
  Surf::Tri *tri;
  Surf::Tri *tris = surf->tris;
  surfint *csurfs;
  int isurf, nsurf, hitflag, side;
  double param, oparam;

  // initialize param to avoid error

  param = 0.0;

  // indices for corners making up the cell edges

  int ci[12], cj[12];
  ci[0] = 0; cj[0] = 1;
  ci[1] = 2; cj[1] = 3;
  ci[2] = 4; cj[2] = 5;
  ci[3] = 6; cj[3] = 7;

  ci[4] = 0; cj[4] = 2;
  ci[5] = 1; cj[5] = 3;
  ci[6] = 4; cj[6] = 6;
  ci[7] = 5; cj[7] = 7;

  ci[8] = 0; cj[8] = 4;
  ci[9] = 1; cj[9] = 5;
  ci[10] = 2; cj[10] = 6;
  ci[11] = 3; cj[11] = 7;

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    // if no surfs, continue

    nsurf = cells[icell].nsurf;
    if (!nsurf) continue;

    for (int d = 0; d < 3; d++) {
      cl[d] = cells[icell].lo[d];
      ch[d] = cells[icell].hi[d];
    }

    cx[0] = cx[2] = cx[4] = cx[6] = cl[0];
    cx[1] = cx[3] = cx[5] = cx[7] = ch[0];

    cy[0] = cy[1] = cy[4] = cy[5] = cl[1];
    cy[2] = cy[3] = cy[6] = cy[7] = ch[1];

    cz[0] = cz[1] = cz[2] = cz[3] = cl[2];
    cz[4] = cz[5] = cz[6] = cz[7] = ch[2];

    // smallest cell length
    // should only have to do this once

    mind = MIN(ch[0]-cl[0], ch[1]-cl[1]);
    mind = MIN(mind, ch[2]-cl[2]);

    // determine corner values

    csurfs = cells[icell].csurfs;
    for (int ic = 0; ic < nedge; ic++) {
      i = ci[ic];
      pi[0] = cx[i];
      pi[1] = cy[i];
      pi[2] = cz[i];

      j = cj[ic];
      pj[0] = cx[j];
      pj[1] = cy[j];
      pj[2] = cz[j];

      // test all surfs+corners to see if any hit

      for (int m = 0; m < nsurf; m++) {
        isurf = csurfs[m];
        tri = &tris[isurf];
        hitflag = corner_hit3d(pi, pj, tri, param, side);

        // need to take care of values near machine precision

        if (param < EPSILON_GRID*mind) param = 0.0;
        if ((1.0-param) < EPSILON_GRID*mind) param = 1.0;
        oparam = 1.0-param;

        // once a hit is found

        if (hitflag) {
          if (ic<=3) {
            n1 = 1;
            n2 = 0;
          } else if (ic<=7) {
            n1 = 3;
            n2 = 2;
          } else {
            n1 = 5;
            n2 = 4;
          }

          if (ivalues[icell][i][n1] > param || ivalues[icell][i][n1] < 0)
            ivalues[icell][i][n1] = param;
          if (ivalues[icell][j][n2] > oparam || ivalues[icell][j][n2] < 0)
            ivalues[icell][j][n2] = oparam;

          if ((mvalues[icell][i] < 0 || param <= mvalues[icell][i]) &&
               svalues[icell][i] != 2) {
            if (param == 0) svalues[icell][i] = 0;

            // conflicting sides from two surfaces meeting at corner

            else if (fabs(mvalues[icell][i]-param) < EPSILON_GRID
              && svalues[icell][i] != side) svalues[icell][i] = 2;
            else svalues[icell][i] = side;
            mvalues[icell][i] = param;
          }

          if ((mvalues[icell][j] < 0 || oparam <= mvalues[icell][j]) && 
               svalues[icell][j] != 2) {
            if (oparam == 0) svalues[icell][j] = 0;
            else if (fabs(mvalues[icell][j]-oparam) < EPSILON_GRID
              && svalues[icell][j] != !side) svalues[icell][j] = 2;
            else svalues[icell][j] = !side;
            mvalues[icell][j] = oparam;
          }
        } // end "if" hitflag
      } // end "for" for surfaces
    } // end "for" for corners in cellfe
  } // end "for" for grid cells

  return;
}

/* ----------------------------------------------------------------------
   sync all copies of values associated with each corner to neighboring
   procs to ensure consistency between procs. Below, the behavior for 
   each option is shown:

   1) if SVAL, sets all svalues to max value (if one is inside, all are inside)
   2) if IVAL, sets all ivalues to min value (closest vertex to corner)
   3) if CVAL, sets all cvalues to max value (handled in previous routines)
   4) if INVAL, sets all invalues to max value (handled in previous routines)

   Temporary arrays are used for cvalues and ivalues to ensure so overwriting
   during the sync operation
------------------------------------------------------------------------- */

void CreateISurf::sync(int which)
{
  int i,j,ix,iy,iz,jx,jy,jz,ixfirst,iyfirst,izfirst,jcorner,jedge,jin;
  int icell,jcell,njcell;
  double dtotal[ninner], dtemp;

  comm_neigh_corners(which);

  // perform update of corner pts for all my owned grid cells
  //   using contributions from all cells that share the corner point
  // insure order of numeric operations will give exact same answer
  //   for all Ncorner duplicates of a corner point (stored by other cells)

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  for (icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    ix = ixyz[icell][0];
    iy = ixyz[icell][1];
    iz = ixyz[icell][2];

    // loop over corner points

    for (i = 0; i < ncorner; i++) {

      // ixyz first = offset from icell of lower left cell of 2x2x2 stencil
      //              that shares the Ith corner point

      ixfirst = (i % 2) - 1;
      iyfirst = (i/2 % 2) - 1;
      if (dim == 2) izfirst = 0;
      else izfirst = (i / 4) - 1;

      // loop over 2x2x2 stencil of cells that share the corner point
      // also works for 2d, since izfirst = 0

      for (j = 0; j < ninner; j++) dtotal[j] = -1.0;
      jcorner = ncorner;

      for (jz = izfirst; jz <= izfirst+1; jz++) {
        for (jy = iyfirst; jy <= iyfirst+1; jy++) {
          for (jx = ixfirst; jx <= ixfirst+1; jx++) {
            jcorner--;

            // check if neighbor cell is within bounds of ablate grid

            if (ix+jx < 1 || ix+jx > nxyz[0]) continue;
            if (iy+jy < 1 || iy+jy > nxyz[1]) continue;
            if (iz+jz < 1 || iz+jz > nxyz[2]) continue;

            // jcell = local index of (jx,jy,jz) neighbor cell of icell

            jcell = walk_to_neigh(icell,jx,jy,jz);

            // update total with one corner point of jcell
            // jcorner descends from ncorner

            if (jcell < nglocal) {
              if (which == SVAL) {
                if (svalues[jcell][jcorner] < 2)
                  dtotal[0] = MAX(dtotal[0],
                    static_cast<double>(svalues[jcell][jcorner]));
              } else if (which == IVAL) {
                for (jin = 0; jin < ninner; jin++) {
                  dtemp = ivalues[jcell][jcorner][jin];
                  if (dtemp >= 0) {
                    if (dtotal[jin] < 0) dtotal[jin] = dtemp;
                    else dtotal[jin] = MAX(dtotal[jin],dtemp);
                  }
                }
              } else if (which == CVAL) {
                dtotal[0] =
                  MAX(dtotal[0],tmp_cvalues[jcell][jcorner]);
              } else if (which == INVAL) {
                for (jin = 0; jin < ninner; jin++) {
                  dtemp = tmp_invalues[jcell][jcorner][jin];
                  if (dtemp >= 0) {
                    if (dtotal[jin] < 0) dtotal[jin] = dtemp;
                    else dtotal[jin] = MAX(dtotal[jin],dtemp);
                  }
                }
              }
            } else {
              if (which == SVAL) {
                if (sghost[jcell-nglocal][jcorner] < 2)
                  dtotal[0] =
                    MAX(dtotal[0],
                    static_cast<double>(sghost[jcell-nglocal][jcorner]));
              } else if (which == IVAL) {
                for (jin = 0; jin < ninner; jin++) {
                  dtemp = ighost[jcell-nglocal][jcorner][jin];
                  if (dtemp >= 0) {
                    if (dtotal[jin] < 0) dtotal[jin] = dtemp;
                    else dtotal[jin] = MAX(dtotal[jin],dtemp);
                  }
                }
              } else if (which == CVAL) {
                dtotal[0] =
                  MAX(dtotal[0],cghost[jcell-nglocal][jcorner]);
              } else if (which == INVAL) {
                for (jin = 0; jin < ninner; jin++) {
                  dtemp = inghost[jcell-nglocal][jcorner][jin];
                  if (dtemp >= 0) {
                    if (dtotal[jin] < 0) dtotal[jin] = dtemp;
                    else dtotal[jin] = MAX(dtotal[jin],dtemp);
                  }
                }
              }
            }

          } // end jx
        } // end jy
      } // end jz

      if (which == SVAL) svalues[icell][i] = static_cast<int>(dtotal[0]);
      else if (which == IVAL) {
        for (jin = 0; jin < ninner; jin++)
          ivalues[icell][i][jin] = dtotal[jin];
      } else if (which == CVAL) {
        cvalues[icell][i] = MAX(dtotal[0],0.0);
      } else if (which == INVAL) {
        for (jin = 0; jin < ninner; jin++)
          invalues[icell][i][jin] = MAX(dtotal[jin],0.0);
      }

    } // end corners
  } // end cells
  return;
}

/* ----------------------------------------------------------------------
   sync all copies of corner points values across procs for all owned grid cells
   algorithm by computing average
------------------------------------------------------------------------- */

void CreateISurf::sync_voxels()
{
  int i,j,ix,iy,iz,jx,jy,jz,ixfirst,iyfirst,izfirst,jcorner,jedge,jin;
  int icell,jcell,njcell;
  double vol_avg;

  comm_neigh_corners(CVAL);

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  for (icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    ix = ixyz[icell][0];
    iy = ixyz[icell][1];
    iz = ixyz[icell][2];

    // loop over corner points

    for (i = 0; i < ncorner; i++) {

      vol_avg = 0.0;

      ixfirst = (i % 2) - 1;
      iyfirst = (i/2 % 2) - 1;
      if (dim == 2) izfirst = 0;
      else izfirst = (i / 4) - 1;

      jcorner = ncorner;

      for (jz = izfirst; jz <= izfirst+1; jz++) {
        for (jy = iyfirst; jy <= iyfirst+1; jy++) {
          for (jx = ixfirst; jx <= ixfirst+1; jx++) {
            jcorner--;

            if (ix+jx < 1 || ix+jx > nxyz[0]) continue;
            if (iy+jy < 1 || iy+jy > nxyz[1]) continue;
            if (iz+jz < 1 || iz+jz > nxyz[2]) continue;

            jcell = walk_to_neigh(icell,jx,jy,jz);

            if (jcell < nglocal) vol_avg += tmp_cvalues[jcell][jcorner];
            else vol_avg += cghost[jcell-nglocal][jcorner];

          }
        }
      }

      // set corner value as average of cell volumes around stencil

      if (dim == 2) vol_avg *= 0.25;
      else vol_avg *= 0.125;

      cvalues[icell][i] = vol_avg;

    }
  }
  return;
}

/* ----------------------------------------------------------------------
   communicate my side, intersection, or corner values that are shared by
   neighbor cells each corner point is shared by N cells, less on borders
   done via irregular comm

   Modified versions of comm_neigh_corners from fix_ablate
------------------------------------------------------------------------- */

void CreateISurf::comm_neigh_corners(int which)
{
  int i,j,k,m,n,ix,iy,iz,ixfirst,iyfirst,izfirst,jx,jy,jz;
  int icell,ifirst,jcell,proc,ilocal;

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  memory->grow(numsend,nglocal,"createisurf:numsend");

  // make list of datums to send to neighbor procs
  // 8 or 26 cells surrounding icell need icell's cdelta info
  // but only if they are owned by a neighbor proc
  // insure icell is only sent once to same neighbor proc
  // also set proclist and locallist for each sent datum

  int nsend = 0;

  for (icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    ix = ixyz[icell][0];
    iy = ixyz[icell][1];
    iz = ixyz[icell][2];
    ifirst = nsend;

    // loop over 3x3x3 stencil of neighbor cells centered on icell

    for (jz = -1; jz <= 1; jz++) {
      for (jy = -1; jy <= 1; jy++) {
        for (jx = -1; jx <= 1; jx++) {

          // skip neigh = self

          if (jx == 0 && jy == 0 && jz == 0) continue;

          // check if neighbor cell is within bounds of ablate grid

          if (ix+jx < 1 || ix+jx > nxyz[0]) continue;
          if (iy+jy < 1 || iy+jy > nxyz[1]) continue;
          if (iz+jz < 1 || iz+jz > nxyz[2]) continue;

          // jcell = local index of (jx,jy,jz) neighbor cell of icell

          jcell = walk_to_neigh(icell,jx,jy,jz);

          // add a send list entry of icell to proc != me if haven't already

          proc = cells[jcell].proc;
          if (proc != me) {
            for (j = ifirst; j < nsend; j++)
              if (proc == proclist[j]) break;
            if (j == nsend) {
              if (nsend == maxsend) grow_send();
              proclist[nsend] = proc;
              locallist[nsend++] = cells[icell].id;
            }
          }
        }
      }
    }

    // # of neighbor procs to send icell to

    numsend[icell] = nsend - ifirst;
  }

  // realloc sbuf if necessary
  // ncomm = ilocal + Ncorner svalues (+ Ncorner*ninner)

  int ncomm;
  if (which == SVAL) ncomm = 1 + ncorner;
  else if (which == IVAL) ncomm = 1 + ncorner*ninner;
  else if (which == CVAL) ncomm = 1 + ncorner;
  else if (which == INVAL) ncomm = 1 + ncorner*ninner;

  if (nsend*ncomm > maxsbuf) {
    memory->destroy(sbuf);
    maxsbuf = nsend*ncomm;
    memory->create(sbuf,maxsbuf,"createisurf:sbuf");
  }

  // pack datums to send
  // datum = ilocal of neigh cell on other proc + Ncorner values

  nsend = 0;
  m = 0;

  for (icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    n = numsend[icell];
    for (i = 0; i < n; i++) {
      sbuf[m++] = ubuf(locallist[nsend]).d;
      if (which == SVAL) {
        for (j = 0; j < ncorner; j++)
          sbuf[m++] = static_cast<double> (svalues[icell][j]);
      } else if (which == IVAL) {
        for (j = 0; j < ncorner; j++)
          for (k = 0; k < ninner; k++)
            sbuf[m++] = ivalues[icell][j][k];
      } else if (which == CVAL) {
        for (j = 0; j < ncorner; j++)
          sbuf[m++] = tmp_cvalues[icell][j];
      } else if (which == INVAL) {
        for (j = 0; j < ncorner; j++)
          for (k = 0; k < ninner; k++)
            sbuf[m++] = tmp_invalues[icell][j][k];
      }
      nsend++;
    }
  }

  // perform irregular neighbor comm
  // Comm class manages rbuf memory

  double *rbuf;
  int nrecv = comm->irregular_uniform_neighs(nsend,proclist,(char *) sbuf,
                ncomm*sizeof(double),(char **) &rbuf);

  // realloc val_ghost if necessary

  if (grid->nghost > maxghost) {
    memory->destroy(inghost);
    memory->destroy(cghost);
    memory->destroy(sghost);
    memory->destroy(ighost);
    maxghost = grid->nghost;
    memory->create(inghost,maxghost,ncorner,ninner,"createisurf:inghost");
    memory->create(cghost,maxghost,ncorner,"createisurf:cghost");
    memory->create(sghost,maxghost,ncorner,"createisurf:sghost");
    memory->create(ighost,maxghost,ncorner,ninner,"createisurf:ighost");
  }

  // unpack received data into val_ghost = ghost cell corner points

  cellint cellID;
  Grid::MyHash *hash = grid->hash;

  m = 0;
  for (i = 0; i < nrecv; i++) {
    cellID = (cellint) ubuf(rbuf[m++]).u;
    ilocal = (*hash)[cellID];
    icell = ilocal - nglocal;

    if (which == SVAL) {
      for (j = 0; j < ncorner; j++)
        sghost[icell][j] = static_cast<int> (rbuf[m++]);
    } else if (which == IVAL) {
      for (j = 0; j < ncorner; j++)
        for (k = 0; k < ninner; k++)
          ighost[icell][j][k] = rbuf[m++];
    } else if (which == CVAL) {
      for (j = 0; j < ncorner; j++)
        cghost[icell][j] = rbuf[m++];
    } else if (which == INVAL) {
      for (j = 0; j < ncorner; j++)
        for (k = 0; k < ninner; k++)
          inghost[icell][j][k] = rbuf[m++];
    }
  }
}

/* ----------------------------------------------------------------------
   walk to neighbor of icell, offset by (jx,jy,jz)
   walk first by x, then by y, last by z
   return jcell = local index of neighbor cell
------------------------------------------------------------------------- */

int CreateISurf::walk_to_neigh(int icell, int jx, int jy, int jz)
{
  Grid::ChildCell *cells = grid->cells;

  int jcell = icell;

  if (jx < 0) {
    if (grid->neigh_decode(cells[jcell].nmask,XLO) != NCHILD)
      error->one(FLERR,"Walk to neighbor cell failed");
    jcell = cells[jcell].neigh[0];
  } else if (jx > 0) {
    if (grid->neigh_decode(cells[jcell].nmask,XHI) != NCHILD)
      error->one(FLERR,"Walk to neighbor cell failed");
    jcell = cells[jcell].neigh[1];
  }

  if (jy < 0) {
    if (grid->neigh_decode(cells[jcell].nmask,YLO) != NCHILD)
      error->one(FLERR,"Walk to neighbor cell failed");
    jcell = cells[jcell].neigh[2];
  } else if (jy > 0) {
    if (grid->neigh_decode(cells[jcell].nmask,YHI) != NCHILD)
      error->one(FLERR,"Walk to neighbor cell failed");
    jcell = cells[jcell].neigh[3];
  }

  if (jz < 0) {
    if (grid->neigh_decode(cells[jcell].nmask,ZLO) != NCHILD)
      error->one(FLERR,"Walk to neighbor cell failed");
    jcell = cells[jcell].neigh[4];
  } else if (jz > 0) {
    if (grid->neigh_decode(cells[jcell].nmask,ZHI) != NCHILD)
      error->one(FLERR,"Walk to neighbor cell failed");
    jcell = cells[jcell].neigh[5];
  }

  return jcell;
}

/* ----------------------------------------------------------------------
   reallocate send vectors
------------------------------------------------------------------------- */

void CreateISurf::grow_send()
{
  maxsend += DELTASEND;
  memory->grow(proclist,maxsend,"createisurf:proclist");
  memory->grow(locallist,maxsend,"createisurf:locallist");
}

/* ----------------------------------------------------------------------
   set corner values for cells fully in or out for both 2D and 3D
------------------------------------------------------------------------- */

void CreateISurf::set_inout()
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  int itype, sval;
  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    // itype = 1 - fully outside
    // itype = 2 - fully inside
    // itype = 3 - has surfaces
    // cannot just use types, if surface on corner
    // .. can be either 2 or 3 depending on which corner (ambiguity)

    itype = cinfo[icell].type;

    if (itype == 2) sval = 1; // fully inside so set all corner values to max
    else if (itype == 1) sval = 0; // fully outside so set all corners to min
    else continue;

    for (int m = 0; m < ncorner; m++)
      if (svalues[icell][m] < 0) svalues[icell][m] = sval;
  }
}


/* ----------------------------------------------------------------------
   resolve unknown side values (2D version)
------------------------------------------------------------------------- */

int CreateISurf::find_side_2d()
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  int filled, attempt;
  int jcorner, jcell;
  int ix, iy, ixfirst, iyfirst;
  int n1, n2, na1, na2;

  filled = attempt = 0;
  while (filled == 0) {
    filled = 1;
    for (int icell = 0; icell < nglocal; icell++) {
      if (!(cinfo[icell].mask & groupbit)) continue;
      if (cells[icell].nsplit <= 0) continue;

      ix = ixyz[icell][0];
      iy = ixyz[icell][1];

      for (int i = 0; i < ncorner; i++) {
        if (svalues[icell][i] == 0 || svalues[icell][i] == 1) continue;

        ixfirst = (i % 2) - 1;
        iyfirst = (i/2 % 2) - 1;

        // check around corner point

        jcorner = ncorner;

        for (int jy = iyfirst; jy <= iyfirst+1; jy++) {
          for (int jx = ixfirst; jx <= ixfirst+1; jx++) {
            jcorner--;

            // check if neighbor cell is within bounds of ablate grid

            if (ix+jx < 1 || ix+jx > nxyz[0]) continue;
            if (iy+jy < 1 || iy+jy > nxyz[1]) continue;

            // n are the corners next to corner i
            // na is n relative to corner i

            if (jcorner == 3) {
              n1 = 1; na1 = 2;
              n2 = 2; na2 = 0;
            } else if (jcorner == 2) {
              n1 = 0; na1 = 2;
              n2 = 3; na2 = 1;
            } else if (jcorner == 1) {
              n1 = 0; na1 = 0;
              n2 = 3; na2 = 3;
            } else {
              n1 = 1; na1 = 1;
              n2 = 2; na2 = 3;
            }

            // jcell = local index of (jx,jy,jz) neighbor cell of icell

            jcell = walk_to_neigh(icell,jx,jy,0);

            // compare with neighbor as reference

            if (jcell < nglocal) {

              // try first neighbor

              int itemp = svalues[jcell][n1];
              if (itemp == 0 || itemp == 1) {
                if (ivalues[jcell][jcorner][na1] <= 0)
                  svalues[icell][i] = itemp;
                else {
                  if (itemp == 0) svalues[icell][i] = 1;
                  else svalues[icell][i] = 0;
                }
                continue;
              }

              // try second neighbor

              itemp = svalues[jcell][n2];
              if (itemp == 0 || itemp == 1) {
                if (ivalues[jcell][jcorner][na2] <= 0)
                  svalues[icell][i] = itemp;
                else {
                  if (itemp == 0) svalues[icell][i] = 1;
                  else svalues[icell][i] = 0;
                }
                continue;
              }
            } else {

              // try first neighbor

              int itemp = sghost[jcell-nglocal][n1];
              if (itemp == 0 || itemp == 1) {
                if (ighost[jcell-nglocal][jcorner][na1] <= 0)
                  svalues[icell][i] = itemp;
                else {
                  if (itemp == 0) svalues[icell][i] = 1;
                  else svalues[icell][i] = 0;
                }
                continue;
              }

              // try second neighbor

              itemp = sghost[jcell-nglocal][n2];
              if (itemp == 0 || itemp == 1) {
                if (ighost[jcell-nglocal][jcorner][na2] <= 0)
                  svalues[icell][i] = itemp;
                else {
                  if (itemp == 0) svalues[icell][i] = 1;
                  else svalues[icell][i] = 0;
                }
                continue;
              }

            } // end jcell if nlocal
          } // end jx
        } // end jy
      } // end corners
    } // end icell

    for (int icell = 0; icell < nglocal; icell++) {
      if (!(cinfo[icell].mask & groupbit)) continue;
      if (cells[icell].nsplit <= 0) continue;
      for (int ic = 0; ic < ncorner; ic++)
        if (svalues[icell][ic] < 0) filled = 0;
    }

    // updating sides may inform on correct side for adjacent corners
    // should only iterate a few times

    attempt++;
    if (attempt > 20) return 0;

  } // end while

  return 1;
}

/* ----------------------------------------------------------------------
   resolve unknown side values (3D version)
------------------------------------------------------------------------- */

int CreateISurf::find_side_3d()
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  int filled, attempt;
  int jcorner, jcell;
  int ix, iy, iz, ixfirst, iyfirst, izfirst;
  int n1, n2, n3, na1, na2, na3;

  filled = attempt = 0;
  while (filled == 0) {
    filled = 1;
    for (int icell = 0; icell < nglocal; icell++) {
      if (!(cinfo[icell].mask & groupbit)) continue;
      if (cells[icell].nsplit <= 0) continue;

      ix = ixyz[icell][0];
      iy = ixyz[icell][1];
      iz = ixyz[icell][2];

      for (int i = 0; i < ncorner; i++) {
        if (svalues[icell][i] > -1 && svalues[icell][i] < 2) continue;

        ixfirst = (i % 2) - 1;
        iyfirst = (i/2 % 2) - 1;
        izfirst = (i / 4) - 1;

        // check around corner point

        jcorner = ncorner;

        for (int jz = izfirst; jz <= izfirst+1; jz++) {
          for (int jy = iyfirst; jy <= iyfirst+1; jy++) {
            for (int jx = ixfirst; jx <= ixfirst+1; jx++) {
              jcorner--;

              // check if neighbor cell is within bounds of ablate grid

              if (ix+jx < 1 || ix+jx > nxyz[0]) continue;
              if (iy+jy < 1 || iy+jy > nxyz[1]) continue;
              if (iz+jz < 1 || iz+jz > nxyz[2]) continue;

              // n are the corners next to corner i
              // na is n relative to corner i
              // these don't seem right

              if (jcorner == 7) { // 0
                n1 = 3; na1 = 4;
                n2 = 5; na2 = 2;
                n3 = 6; na3 = 0;
              } else if (jcorner == 6) {
                n1 = 2; na1 = 4;
                n2 = 4; na2 = 2;
                n3 = 7; na3 = 1;
              } else if (jcorner == 5) {
                n1 = 1; na1 = 4;
                n2 = 4; na2 = 0;
                n3 = 7; na3 = 3;
              } else if (jcorner == 4) {
                n1 = 0; na1 = 4;
                n2 = 5; na2 = 1;
                n3 = 6; na3 = 3;
              } else if (jcorner == 3) {
                n1 = 1; na1 = 2;
                n2 = 2; na2 = 0;
                n3 = 7; na3 = 5;
              } else if (jcorner == 2) {
                n1 = 0; na1 = 2;
                n2 = 3; na2 = 1;
                n3 = 6; na3 = 5;
              } else if (jcorner == 1) {
                n1 = 0; na1 = 0;
                n2 = 3; na2 = 3;
                n3 = 5; na3 = 5;
              } else {
                n1 = 1; na1 = 1;
                n2 = 2; na2 = 3;
                n3 = 4; na3 = 5;
              }

              // jcell = local index of (jx,jy,jz) neighbor cell of icell

              jcell = walk_to_neigh(icell,jx,jy,jz);

              // compare with neighbor as reference

              if (jcell < nglocal) {

                // try first neighbor

                int itemp = svalues[jcell][n1];
                if (itemp > -1 && itemp < 2) {
                  if (ivalues[jcell][jcorner][na1] <= 0)
                    svalues[icell][i] = itemp;
                  else {
                    if (itemp < 1) svalues[icell][i] = 1;
                    else svalues[icell][i] = 0;
                  }
                  continue;
                }

                // try second neighbor

                itemp = svalues[jcell][n2];
                if (itemp > -1 && itemp < 2) {
                  if (ivalues[jcell][jcorner][na2] <= 0)
                    svalues[icell][i] = itemp;
                  else {
                    if (itemp < 1) svalues[icell][i] = 1;
                    else svalues[icell][i] = 0;
                  }
                  continue;
                }

                // try third neighbor

                itemp = svalues[jcell][n3];
                if (itemp > -1 && itemp < 2) {
                  if (ivalues[jcell][jcorner][na3] <= 0)
                    svalues[icell][i] = itemp;
                  else {
                    if (itemp < 1) svalues[icell][i] = 1;
                    else svalues[icell][i] = 0;
                  }
                  continue;
                }
              } else {

                // try first neighbor

                int itemp = sghost[jcell-nglocal][n1];
                if (itemp > -1 && itemp < 2) {
                  if (ighost[jcell-nglocal][jcorner][na1] <= 0)
                    svalues[icell][i] = itemp;
                  else {
                    if (itemp < 1) svalues[icell][i] = 1;
                    else svalues[icell][i] = 0;
                  }
                  continue;
                }

                // try second neighbor

                itemp = sghost[jcell-nglocal][n2];
                if (itemp > -1 && itemp < 2) {
                  if (ighost[jcell-nglocal][jcorner][na2] <= 0)
                    svalues[icell][i] = itemp;
                  else {
                    if (itemp < 1) svalues[icell][i] = 1;
                    else svalues[icell][i] = 0;
                  }
                  continue;
                }

                // try third neighbor

                itemp = sghost[jcell-nglocal][n3];
                if (itemp > -1 && itemp < 2) {
                  if (ighost[jcell-nglocal][jcorner][na3] <= 0)
                    svalues[icell][i] = itemp;
                  else {
                    if (itemp < 1) svalues[icell][i] = 1;
                    else svalues[icell][i] = 0;
                  }
                  continue;
                }
              } // end jcell if nlocal
            } // end jx
          } // end jy
        } // end jz
      } // end corners
    } // end icell


    for (int icell = 0; icell < nglocal; icell++) {
      if (!(cinfo[icell].mask & groupbit)) continue;
      if (cells[icell].nsplit <= 0) continue;
      for (int i = 0; i < ncorner; i++)
        if (svalues[icell][i] < 0 || svalues[icell][i] > 1) filled = 0;
    }

    // updating sides may inform on correct side for adjacent corners
    // should only iterate a few times

    attempt++;
    if (attempt > 20) return 0;

  } // end while

  return 1;
}

/* ----------------------------------------------------------------------
   determine corner values for each grid cell. If the inout option is
   used, the corner value only depends on the side values. If the ave
   option is used, corner value is extrapolated from the intersection
   values between the surfaces and the cell edges
------------------------------------------------------------------------- */

void CreateISurf::set_cvalues()
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  cout = 0.0;
  cin = 255.0;

  // assign corner values based on choice

  if (ctype == INOUT) set_cvalues_inout();
  else if (ctype == VOXEL) set_cvalues_voxel();
  else if (ctype == AVE) set_cvalues_ave();
  else if (ctype == INNER) set_cvalues_inner();

  return;
}

/* ----------------------------------------------------------------------
   set corner point values depending on if corner is in or out of surface
------------------------------------------------------------------------- */

void CreateISurf::set_cvalues_inout()
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;
    for (int ic = 0; ic < ncorner; ic++) {
      if (svalues[icell][ic] == 0) tmp_cvalues[icell][ic] = cout;
      else if (svalues[icell][ic] == 1) tmp_cvalues[icell][ic] = cin;
      else error->one(FLERR,"Corner point not marked inside or outside");
    } // end corners
  } // end "for" for grid cells
}

/* ----------------------------------------------------------------------
   find corner point values using averaging
------------------------------------------------------------------------- */

void CreateISurf::set_cvalues_voxel()
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  double dx,dy,dz;
  double cvol, sfrac;

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    cvol = cinfo[icell].volume;
    dx = cells[icell].hi[0] - cells[icell].lo[0];
    dy = cells[icell].hi[1] - cells[icell].lo[1];
    dz = cells[icell].hi[2] - cells[icell].lo[2];
    sfrac = (dx*dy*dz - cvol) / (dx*dy*dz);

    if (sfrac < 0.0 || sfrac > 1.0) error->one(FLERR,"bad volume");

    for (int ic = 0; ic < ncorner; ic++) {
      if (svalues[icell][ic] == 0) tmp_cvalues[icell][ic] = cout;
      else if (svalues[icell][ic] == 1) tmp_cvalues[icell][ic] = sfrac*cin;
      //else error->one(FLERR,"bad svalues");
    }
  } // end grid cells

  return;
}

/* ----------------------------------------------------------------------
   find corner point values using averaging
------------------------------------------------------------------------- */

void CreateISurf::set_cvalues_ave()
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  int nval;
  double ivalsum;
  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;
    for (int ic = 0; ic < ncorner; ic++) {

      ivalsum = 0.0;
      if (svalues[icell][ic] == 0) tmp_cvalues[icell][ic] = cout;
      else {
        nval = 0;
        for (int iin = 0; iin < ninner; iin++) {
          if (ivalues[icell][ic][iin] >= 0) {
            ivalsum += ivalues[icell][ic][iin];
            nval++;
          }
        }

        // no intersections

        if (nval == 0) tmp_cvalues[icell][ic] = cin;
        else {
          ivalsum /= nval;
          if (ivalsum > 1.0) error->one(FLERR,"Calculated vertex location outside cell");
          tmp_cvalues[icell][ic] = MAX(param2cval(ivalsum,0.0),0.0);
        }

      } // end svalues
    } // end corners
  } // end grid cells
}

/* ----------------------------------------------------------------------
   set inner values of corner points
------------------------------------------------------------------------- */

void CreateISurf::set_cvalues_inner()
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  // first set corner point values of fully inside and outside cells

  int refsval, allsame;
  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;
    refsval = svalues[icell][0];
    allsame = 1;
    for (int ic = 1; ic < ncorner; ic++)
      if (svalues[icell][ic] != refsval)
        allsame = 0;

    if (allsame) {
      for (int ic = 0; ic < ncorner; ic++) {
        for (int iin = 0; iin < ninner; iin++) {
          if (refsval==1)
            tmp_invalues[icell][ic][iin] = cin;
          else
            tmp_invalues[icell][ic][iin] = cout;
        }
      }
    }
  } // end cells

  // now handle the overlap cells

  double ival, cval;

  // value of inside corner point next to surface

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;
    for (int ic = 0; ic < ncorner; ic++) {

      for (int k = 0; k < ninner; k++) {

        // bound the intersection values

        ival = ivalues[icell][ic][k];
        
        if (ival <= 0) { // no intersection this edge
          if (svalues[icell][ic] == 0) cval = cout;
          else cval = cin;
        } else if (svalues[icell][ic] == 1) {
          cval = param2cval(ival,0.0);
        } else {
          cval = param2cval(ival,255.0);
        }

        tmp_invalues[icell][ic][k] = cval;

      }
    } // end ncorner
  } // end cells

  return;
}

/* ----------------------------------------------------------------------
   determines if surface is inline with cell edge
   side = 0,1 -> [out,in]
   2d version
------------------------------------------------------------------------- */

int CreateISurf::corner_hit2d(double *p1, double *p2,
                              Surf::Line *line, double &param, int &side)
{
  // try intersect first

  int h, tside;
  double tparam;
  double d3dum[3];
  h = Geometry::line_line_intersect(p1,p2,line->p1,line->p2,line->norm,
                                    d3dum,tparam,tside);
  if (h) {
    if (tside == 1 || tside == 2 || tside == 5) {
      side = 1;
      param = tparam;
      return true;
    } else {
      side = 0;
      param = tparam;
      return true;
    }
  }

  // perturbed points

  double p1p[3];
  double p2p[3];

  double dx[8], dy[8];
  double dp = mind/1000.0;
  dx[0] = dx[1] = dx[6] = dp;
  dx[2] = dx[5] = 0;
  dx[3] = dx[4] = dx[7] = -dp;

  dy[0] = dy[2] = dy[7] = dp;
  dy[1] = dy[4] = 0;
  dy[3] = dy[5] = dy[6] = -dp;

  for (int i = 0; i < 8; i++) {
    p1p[2] = p1[2];
    p2p[2] = p2[2];

    p1p[0] = p1[0] + dx[i];
    p1p[1] = p1[1] + dy[i];

    p2p[0] = p2[0] + dx[i];
    p2p[1] = p2[1] + dy[i];

    h = Geometry::line_line_intersect(p1p,p2p,line->p1,line->p2,
                                      line->norm,d3dum,tparam,tside);
    if (h) {
      if (tside == 1 || tside == 2 || tside == 5) {

        //side = 2;

        side = 1;
        if (tparam<0.5) param = 0.0;
        else param = 1.0;
        return true;
      } else {

        //side = 2;

        side = 0;
        if (tparam<0.5) param = 0.0;
        else param = 1.0;
        return true;
      }
    }
  }

  // true miss

  return false;
}

/* ----------------------------------------------------------------------
   determines if surface is inline with cell face
   side = 0,1 -> [out,in]
   3d version
------------------------------------------------------------------------- */

int CreateISurf::corner_hit3d(double *p1, double *p2,
                              Surf::Tri* tri, double &param, int &side)
{
  // try intersect first

  int h, tside;
  double tparam;
  double d3dum[3];
  h = Geometry::line_tri_intersect(p1,p2,tri->p1,tri->p2,tri->p3,
      tri->norm,d3dum,tparam,tside);

  if (h) {
    if (tside == 1 || tside == 2 || tside == 5) {
      side = 1;
      param = tparam;
      return true;
    } else {
      side = 0;
      param = tparam;
      return true;
    }
  }

  // if miss, maybe surface very close to corner/edge
  // perturbed points

  double p1p[3];
  double p2p[3];

  double dx[26], dy[26], dz[26];
  double dp = mind/1000.0;
  dx[0] = dx[2] = dx[3] = dx[4] =
    dx[15] = dx[16] = dx[18] = dx[23] = dx[24] = dp;
  dx[7] = dx[9] = dx[10] = dx[11] =
    dx[14] = dx[17] = dx[19] = dx[22] = dx[25] = -dp;
  dx[1] = dx[5] = dx[6] = dx[8] =
    dx[12] = dx[13] = dx[20] = dx[21] = 0;

  dy[0] = dy[1] = dy[3] = dy[5] =
    dy[14] = dy[16] = dy[19] = dy[21] = dy[25] = dp;
  dy[7] = dy[8] = dy[10] = dy[12] =
    dy[15] = dy[17] = dy[18] = dy[20] = dy[24] = -dp;
  dy[2] = dy[4] = dy[6] = dy[9] =
    dy[11] = dy[13] = dy[22] = dy[23] = 0;

  dz[0] = dz[1] = dz[2] = dz[6] =
    dz[14] = dz[15] = dz[17] = dz[20] = dz[22] = dp;
  dz[7] = dz[8] = dz[9] = dz[13] =
    dz[16] = dz[18] = dz[19] = dz[21] = dz[23] = -dp;
  dz[3] = dz[4] = dz[5] = dz[10] =
    dz[11] = dz[12] = dz[24] = dz[25] = 0;

  for (int i = 0; i < 26; i++) {
    p1p[0] = p1[0] + dx[i];
    p1p[1] = p1[1] + dy[i];
    p1p[2] = p1[2] + dz[i];

    p2p[0] = p2[0] + dx[i];
    p2p[1] = p2[1] + dy[i];
    p2p[2] = p2[2] + dz[i];

    h = Geometry::line_tri_intersect(p1p,p2p,tri->p1,tri->p2,tri->p3,
        tri->norm,d3dum,tparam,tside);
    if (h) {
      if (tside == 1 || tside == 2 || tside == 5) {
        side = 1;
        if (tparam<0.5) param = 0.0;
        else param = 1.0;
        return true;
      } else {
        side = 0;
        if (tparam<0.5) param = 0.0;
        else param = 1.0;
        return true;
      }
    }
  }

  // true miss

  return false;
}

/* ----------------------------------------------------------------------
   find inside corner value from outside corner value and intersection
   value
------------------------------------------------------------------------- */

double CreateISurf::param2cval(double param, double v1)
{
  // param is proportional to cell length so
  // ... lo = 0; hi = 1
  // trying to find v0
  // param = (thresh  - v0) / (v1 - v0)

  double v0 = (thresh - v1*param) / (1.0 - param);

  // bound by limits

  v0 = MIN(v0,cin);
  v0 = MAX(v0,cout);

  return v0;
}

/* ----------------------------------------------------------------------
   removes all explicit surfaces
   copied from RemoveSurf
------------------------------------------------------------------------- */

void CreateISurf::remove_old()
{
  if (particle->exist) particle->sort();

  cuvalues = NULL;
  int *index_custom = new int[surf->ncustom];
  int ncustom = 0;

  if (surf->ncustom) {
    for (int i = 0; i < surf->ncustom; i++) {
      if (!surf->ename[i]) continue;
      index_custom[ncustom++] = i;
    }
  }

  llines = NULL;
  ltris = NULL;
  int nsurf = surf->nown;

  int nbytes;
  if (dim == 2) nbytes = sizeof(Surf::Line);
  else nbytes = sizeof(Surf::Tri);

  if (dim == 2) {
    llines = (Surf::Line *) memory->smalloc(nsurf*nbytes,"createisurf:lines");
    memcpy(llines,surf->mylines,nsurf*nbytes);
  } else {
    ltris = (Surf::Tri *) memory->smalloc(nsurf*nbytes,"createisurf:ltris");
    memcpy(ltris,surf->mytris,nsurf*nbytes);
  }

  surf->add_surfs(1,0,llines,ltris,ncustom,index_custom,cuvalues);

  memory->sfree(llines);
  memory->sfree(ltris);
  memory->destroy(cuvalues);
  delete [] index_custom;

  // check that remaining surfs are still watertight

  if (dim == 2) surf->check_watertight_2d();
  else surf->check_watertight_3d();

  // reset grid due to changing surfs
  // assign surfs to grid cells

  surf->setup_owned();
  grid->unset_neighbors();
  grid->remove_ghosts();

  // reassign split cell particles to parent split cell

  if (particle->exist && grid->nsplitlocal) {
    Grid::ChildCell *cells = grid->cells;
    for (int icell = 0; icell < nglocal; icell++)
      if (cells[icell].nsplit > 1)
        grid->combine_split_cell_particles(icell,1);
  }

  grid->clear_surf();

  //if (surf->exist) grid->surf2grid(1); // surf shouldn't exist

  // reassign particles in split cells to sub cell owner

  if (particle->exist && grid->nsplitlocal) {
    Grid::ChildCell *cells = grid->cells;
    for (int icell = 0; icell < nglocal; icell++)
      if (cells[icell].nsplit > 1)
        grid->assign_split_cell_particles(icell);
  }

  // re-setup owned and ghost cell info

  grid->setup_owned();
  grid->acquire_ghosts();
  grid->reset_neighbors();
  comm->reset_neighbors();

  grid->set_inout();
  grid->type_check();
}

/* ----------------------------------------------------------------------
   process command line args
------------------------------------------------------------------------- */

void CreateISurf::process_args(int narg, char **arg)
{
  return;
}

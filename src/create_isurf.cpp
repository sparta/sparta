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

#define EPSILON_GRID 1.0e-3

/* ---------------------------------------------------------------------- */

CreateISurf::CreateISurf(SPARTA *sparta) : Pointers(sparta)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
}

/* ---------------------------------------------------------------------- */

CreateISurf::~CreateISurf()
{

}

/* ---------------------------------------------------------------------- */

void CreateISurf::command(int narg, char **arg)
{
  if (!grid->exist)
    error->all(FLERR,"Cannot create_isurf before grid is defined");
  if(!surf->exist)
    error->all(FLERR,"Must read in surface first with read_surf");
  if (surf->implicit && surf->exist)
    error->all(FLERR,"Cannot have pre-existing implicit surfaces");
  if (!surf->distributed)
    error->all(FLERR,"Explicit surface must be distributed");

  int dim = domain->dimension;

  if (narg != 4) error->all(FLERR,"Illegal create_isurf command");

  // grid group
  ggroup = grid->find_group(arg[0]);
  if (ggroup < 0) error->all(FLERR,"Read_isurf grid group ID does not exist");

  // ablate fix
  char *ablateID = arg[1];
  int ifix = modify->find_fix(ablateID);
  if (ifix < 0)
    error->all(FLERR,"Fix ID for read_surf does not exist");
  if (strcmp(modify->fix[ifix]->style,"ablate") != 0)
    error->all(FLERR,"Fix for read_surf is not a fix ablate");
  ablate = (FixAblate *) modify->fix[ifix];
  if (ggroup != ablate->igroup)
    error->all(FLERR,"Read_surf group does not match fix ablate group");

  // threshold for corner value
  thresh = input->numeric(FLERR,arg[2]);
  if (thresh < 0 || thresh > 255)
    error->all(FLERR,"Create_isurf thresh must be bounded as (0,255)");

  // mode to determine corner values
  if (strcmp(arg[3],"inout") == 0) aveFlag = 0;
  else if (strcmp(arg[3],"ave") == 0) aveFlag = 1;
  else error->all(FLERR,"Unknown surface corner mode called");

  //if(aveFlag && comm->nprocs > 1)
  //  error->all(FLERR,"Create_isurf averaging not possible in parallel");

  // nxyz already takes into account subcells
  // find corner values for all grid cells initially
  // only store those within ggroup when calling ablate->store
  // 0 check uniform grid for entire domain
  grid->check_uniform_group(0,nxyz,corner,xyzsize);

  // corner points needs one more in each dim
  Nxyz = (nxyz[0]+1)*(nxyz[1]+1);
  if(dim == 3) Nxyz *= (nxyz[2]+1);

  // find all corner values
  set_corners();
  memory->destroy(cvalues);
  MPI_Barrier(world);

  // remove old explicit surfaces
  remove_old();

  surf->implicit = 1;
  surf->exist = 1;

  tvalues = NULL; // TODO: Add per-surface type
  int cpushflag = 0; // don't push
  char *sgroupID = arg[0];

  ablate->store_corners(nxyz[0],nxyz[1],nxyz[2],corner,xyzsize,
                  icvalues,tvalues,thresh,sgroupID,cpushflag);

  if (ablate->nevery == 0) modify->delete_fix(ablateID);
  MPI_Barrier(world);
}

/* ----------------------------------------------------------------------
   Determines corner values given an explicit surface. Cvalues then used
   later in ablate to create the implicit surfaces
------------------------------------------------------------------------- */

void CreateISurf::set_corners()
{
  // first shift everything down by thresh
  // later shift back
  cout = 0.0;
  cin = 255.0;

  // cvalues is the grid (no overlaps)
  // icvalues is setting each corner for each cell based on cvalues
  // .. icvalues is cvalues in read_isurf and fix_ablate
  memory->create(cvalues,Nxyz,"createisurf:cvalues");
  for (int i = 0; i < Nxyz; i++) cvalues[i] = -1.0;

  // mvalues stores minimum param
  // also used to determine side
  memory->create(mvalues,Nxyz,"createisurf:mvalues");
  // svalues stores side of minimum param
  memory->create(svalues,Nxyz,"createisurf:svalues");
  // array to keep track of param between corners

  // find corner values for all cells with surfaces
  int npairs; // number of points around corner point

  if(domain->dimension==2) {
    memory->create(icvalues,grid->nlocal,4,"createisurf:icvalues");
    for (int i = 0; i < grid->nlocal; i++) {
      for (int j = 0; j < 4; j++)
        icvalues[i][j] = 0.0;
    }

    // each corner has 4 neighbors
    // -x, +x, -y, +y
    npairs = 4;
    memory->create(ivalues,Nxyz,npairs,"createisurf:ivalues");
    
    for(int ic = 0; ic < Nxyz; ic++) {
      svalues[ic] = -1;
      mvalues[ic] = -1.0;
      for(int jc = 0; jc < npairs; jc++) ivalues[ic][jc] = -1.0;
    }

    // find intersections between edges and surfaces
    surface_edge2d();
  } else {
    memory->create(icvalues,grid->nlocal,8,"createisurf:icvalues");
    for (int i = 0; i < grid->nlocal; i++) {
      for (int j = 0; j < 8; j++)
        icvalues[i][j] = 0.0;
    }

    // each corner has 6 neighbors
    // -x, +x, -y, +y, -z, +z
    npairs = 6;
    memory->create(ivalues,Nxyz,npairs,"createisurf:ivalues");
    
    for(int ic = 0; ic < Nxyz; ic++) {
      svalues[ic] = -1;
      mvalues[ic] = -1.0;
      for(int jc = 0; jc < npairs; jc++) ivalues[ic][jc] = -1.0;
    }
    
    for(int ic = 0; ic < Nxyz; ic++) {
      svalues[ic] = -1;
      mvalues[ic] = -1.0;
      for(int jc = 0; jc < npairs; jc++) ivalues[ic][jc] = -1.0;
    }

    surface_edge3d();
  }

  // fill in corner values based on if grid cell is in or out
  set_inout();
  MPI_Barrier(world);

  // sync svalues
  // might need to sync ivalues as well
  if (me == 0)
    if (screen) fprintf(screen,"Syncing intermediate values ...\n");

  int lsval, gsval;
  double lival, gival;
  int lsgn, gsgn;
  for(int ic = 0; ic < Nxyz; ic++) {
    // sync side values
    if(svalues[ic] < 2) lsval = svalues[ic];
    else lsval = -1;
    MPI_Allreduce(&lsval,&gsval,1,MPI_INT,MPI_MAX,world);
    if(gsval > svalues[ic]) svalues[ic] = gsval;

    // sync intersection values
    for(int jc = 0; jc < npairs; jc++) {
      lival = ivalues[ic][jc];
      if(lival > 0) lsgn = 1;
      else {
        lsgn = 0;
        lival = 0.0;
      }

      MPI_Allreduce(&lsgn,&gsgn,1,MPI_INT,MPI_SUM,world);
      MPI_Allreduce(&lival,&gival,1,MPI_DOUBLE,MPI_SUM,world);
      gsgn = MAX(gsgn,1);
      ivalues[ic][jc] = gival/static_cast<double>(gsgn);
      if(ivalues[ic][jc]>1.0 && me==0) {
        printf("ival: %4.3e; nval: %i\n", ivalues[ic][jc], gsgn);
        error->one(FLERR,"ival too big");
      }
    }
  }
  MPI_Barrier(world);

  // find remaining corner values based on neighbors
  if (me == 0)
    if (screen) fprintf(screen,"Cleanup corners ...\n");
  cleanup();

  // free up memory since these are no longer needed
  memory->destroy(ivalues);
  memory->destroy(mvalues);
  memory->destroy(svalues);

  // initialize corner point matrix for fix-ablate
  if(domain->dimension==2) {
    memory->create(icvalues,grid->nlocal,4,"createisurf:icvalues");
    for (int i = 0; i < grid->nlocal; i++) {
      for (int j = 0; j < 4; j++)
        icvalues[i][j] = 0.0;
    }
  } else {
    memory->create(icvalues,grid->nlocal,8,"createisurf:icvalues");
    for (int i = 0; i < grid->nlocal; i++) {
      for (int j = 0; j < 8; j++)
        icvalues[i][j] = 0.0;
    }
  }

  // DEBUG: check no negative values
  for(int ic = 0; ic < Nxyz; ic++)
    if(cvalues[ic] < 0) printf("negative cvalue\n");

  // sync cvalues
  // may not be necessary
  if (me == 0)
    if (screen) fprintf(screen,"Syncing cvalues ...\n");
  /*double lcval, maxcval;
  for(int ic = 0; ic < Nxyz; ic++) {
    lcval = cvalues[ic];
    MPI_Allreduce(&lcval,&maxcval,1,MPI_DOUBLE,MPI_MAX,world);
    cvalues[ic] = maxcval;
  }*/

  double lcval, gcval;
  for(int ic = 0; ic < Nxyz; ic++) {
    lcval = cvalues[ic];
    if(lcval >= thresh) lsgn = 1;
    else {
      lsgn = 0;
      lcval = 0.0;
    }
    MPI_Allreduce(&lsgn,&gsgn,1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&lcval,&gcval,1,MPI_DOUBLE,MPI_SUM,world);
    gsgn = MAX(gsgn,1);
    cvalues[ic] = gcval/gsgn;
    if(cvalues[ic]>255.0 && me==0) {
      printf("cval: %4.3e; nval: %i\n", gcval, gsgn);
      error->one(FLERR,"cval too big");
    }
    //if(me==0 && cvalues[ic] > thresh)
    //  printf("ic: %i; cval: %4.2e\n", ic, cvalues[ic]);
  }

  MPI_Barrier(world);

  // store into icvalues for generating implicit surface
  Grid::ChildCell *cells = grid->cells;
  int xyzcell, ic[3];
  double lc[3];
  for(int icell = 0; icell < grid->nlocal; icell++) {

    lc[0] = cells[icell].lo[0];
    lc[1] = cells[icell].lo[1];
    lc[2] = cells[icell].lo[2];

    xyzcell = get_cxyz(ic,lc);
    icvalues[icell][0] = cvalues[xyzcell];
    xyzcell = get_corner(ic[0]+1, ic[1], ic[2]);
    icvalues[icell][1] = cvalues[xyzcell];
    xyzcell = get_corner(ic[0], ic[1]+1, ic[2]);
    icvalues[icell][2] = cvalues[xyzcell];
    xyzcell = get_corner(ic[0]+1, ic[1]+1, ic[2]);
    icvalues[icell][3] = cvalues[xyzcell];

    if(icvalues[icell][0] < 0 ||
       icvalues[icell][1] < 0 ||
       icvalues[icell][2] < 0 ||
       icvalues[icell][3] < 0) printf("negative corner value 2d\n");

    if(domain->dimension==3) {
      xyzcell = get_corner(ic[0], ic[1], ic[2]+1);
      icvalues[icell][4] = cvalues[xyzcell];
      xyzcell = get_corner(ic[0]+1, ic[1], ic[2]+1);
      icvalues[icell][5] = cvalues[xyzcell];
      xyzcell = get_corner(ic[0], ic[1]+1, ic[2]+1);
      icvalues[icell][6] = cvalues[xyzcell];
      xyzcell = get_corner(ic[0]+1, ic[1]+1, ic[2]+1);
      icvalues[icell][7] = cvalues[xyzcell];

      if(icvalues[icell][4] < 0 ||
         icvalues[icell][5] < 0 ||
         icvalues[icell][6] < 0 ||
         icvalues[icell][7] < 0) printf("negative corner value 3d\n");
    }
  }
  return;
}

/* ----------------------------------------------------------------------
   Finds intersections and sides for corners in cells with surfaces
------------------------------------------------------------------------- */

void CreateISurf::surface_edge2d()
{
  int i,j,c1,c2,n1,n2;

  Grid::ChildCell *cells = grid->cells;
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

  // number of edges in a cell
  int nedge = 4;
  // indices for corners making up the cell edges
  int ci[4], cj[4];
  ci[0] = 0; cj[0] = 1;
  ci[1] = 1; cj[1] = 3;
  ci[2] = 3; cj[2] = 2;
  ci[3] = 2; cj[3] = 0;

  for(int icell = 0; icell < grid->nlocal; icell++) {

    // if no surfs, continue
    nsurf = cells[icell].nsurf;
    if(!nsurf) continue;

    for(int d = 0; d < 3; d++) {
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

    // determine corner values
    csurfs = cells[icell].csurfs;
    for(int ic = 0; ic < nedge; ic++) {
      i = ci[ic];
      pi[0] = cx[i];
      pi[1] = cy[i];
      pi[2] = cz[i];

      j = cj[ic];
      pj[0] = cx[j];
      pj[1] = cy[j];
      pj[2] = cz[j];

      // get index for p1 and p2
      // c1 is always lower from above
      c1 = get_corner(pi[0], pi[1], pi[2]);
      c2 = get_corner(pj[0], pj[1], pj[2]);

      // test all surfs+corners to see if any hit
      for (int m = 0; m < nsurf; m++) {
        isurf = csurfs[m];
        line = &lines[isurf];

        // always start with lower corner
        // side will always be either 0 or 1
        hitflag = corner_hit2d(pi, pj, line, param, side);

        // need to take care of values near machine precision
        if(param < EPSILON_GRID*mind) param = 0.0;
        if((1.0-param) < EPSILON_GRID*mind) param = 1.0;
        oparam = 1.0-param;

        // once a hit is found
        if(hitflag) {
          if(ic==0) {
            n1 = 1;
            n2 = 0;
          } else if(ic==1) {
            n1 = 3;
            n2 = 2;
          } else if(ic==2) {
            n1 = 0;
            n2 = 1;
          } else {
            n1 = 2;
            n2 = 3;
          }

          if(ivalues[c1][n1] > param || ivalues[c1][n1] < 0) ivalues[c1][n1] = param;
          if(ivalues[c2][n2] > oparam || ivalues[c2][n2] < 0) ivalues[c2][n2] = oparam;

          if(mvalues[c1] < 0 || param <= mvalues[c1]) {
            if(param == 0) svalues[c1] = 0;
            else if(svalues[c1] == 2) 0; // do nothing
            // conflicting sides from two surfaces meeting at corner
            else if(fabs(mvalues[c1]-param) < EPSILON_GRID && svalues[c1] != side) svalues[c1] = 2;
            else svalues[c1] = side;

            mvalues[c1] = param;
          }

          if(mvalues[c2] < 0 || oparam <= mvalues[c2]) {
            if(oparam == 0) svalues[c2] = 0;
            else if(svalues[c2] == 2) 0; // do nothing
            else if(fabs(mvalues[c2]-oparam) < EPSILON_GRID && svalues[c2] != !side) svalues[c2] = 2;
            else svalues[c2] = !side;

            mvalues[c2] = oparam;
          }
        } // end "if" hitflag
      }// end "for" for surfaces
    }// end "for" for corners in cell
  }// end "for" for grid cells

  return;
}

/* ----------------------------------------------------------------------
   Same as above but for 3d
------------------------------------------------------------------------- */

void CreateISurf::surface_edge3d()
{
  int i,j,c1,c2,n1,n2;

  Grid::ChildCell *cells = grid->cells;
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

  // number of edges in a cell
  int nedge = 12;
  // indices for corners making up the cell edges
  int ci[12], cj[12];
  ci[0] = 0; cj[0] = 1;
  ci[1] = 1; cj[1] = 3;
  ci[2] = 3; cj[2] = 2;
  ci[3] = 2; cj[3] = 0;
  ci[4] = 0; cj[4] = 4;
  ci[5] = 1; cj[5] = 5;
  ci[6] = 3; cj[6] = 7;
  ci[7] = 2; cj[7] = 6;
  ci[8] = 4; cj[8] = 5;
  ci[9] = 5; cj[9] = 7;
  ci[10] = 7; cj[10] = 6;
  ci[11] = 6; cj[11] = 4;

  for(int icell = 0; icell < grid->nlocal; icell++) {

    // if no surfs, continue
    nsurf = cells[icell].nsurf;
    if(!nsurf) continue;

    for(int d = 0; d < 3; d++) {
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
    for(int ic = 0; ic < nedge; ic++) {
      i = ci[ic];
      pi[0] = cx[i];
      pi[1] = cy[i];
      pi[2] = cz[i];

      j = cj[ic];
      pj[0] = cx[j];
      pj[1] = cy[j];
      pj[2] = cz[j];

      // get corners
      c1 = get_corner(pi[0], pi[1], pi[2]);
      c2 = get_corner(pj[0], pj[1], pj[2]);

      // test all surfs+corners to see if any hit
      for (int m = 0; m < nsurf; m++) {
        isurf = csurfs[m];
        tri = &tris[isurf];
        hitflag = corner_hit3d(pi, pj, tri, param, side);

        // need to take care of values near machine precision
        if(param < EPSILON_GRID*mind) param = 0.0;
        if((1.0-param) < EPSILON_GRID*mind) param = 1.0;
        oparam = 1.0-param;

        // once a hit is found
        if(hitflag) {

          if(ic==0 || ic==8) {
            n1 = 1;
            n2 = 0;
          } else if(ic==1 || ic==9) {
            n1 = 3;
            n2 = 2;
          } else if(ic==2 || ic==10) {
            n1 = 0;
            n2 = 1;
          } else if(ic==3 || ic==11) {
            n1 = 2;
            n2 = 3;
          } else {
            n1 = 5;
            n2 = 4;
          }

          if(ivalues[c1][n1] > param || ivalues[c1][n1] < 0) ivalues[c1][n1] = param;
          if(ivalues[c2][n2] > oparam || ivalues[c2][n2] < 0) ivalues[c2][n2] = oparam;

          if(mvalues[c1] < 0 || param <= mvalues[c1]) {
            if(param == 0) svalues[c1] = 0;
            else if(svalues[c1] == 2) 0; // do nothing
            // conflicting sides from two surfaces meeting at corner
            else if(fabs(mvalues[c1]-param) < EPSILON_GRID && svalues[c1] != side) svalues[c1] = 2;
            else svalues[c1] = side;

            mvalues[c1] = param;
          }

          if(mvalues[c2] < 0 || oparam <= mvalues[c2]) {
            if(oparam == 0) svalues[c2] = 0;
            else if(svalues[c2] == 2) 0; // do nothing
            else if(fabs(mvalues[c2]-oparam) < EPSILON_GRID && svalues[c2] != !side) svalues[c2] = 2;
            else svalues[c2] = !side;

            mvalues[c2] = oparam;
          }

        } // end "if" hitflag
      }// end "for" for surfaces
    }// end "for" for corners in cellfe
  }// end "for" for grid cells

  return;
}

/* ----------------------------------------------------------------------
  Set corner values for cells fully in or out for both 2D and 3D
------------------------------------------------------------------------- */

void CreateISurf::set_inout()
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  double cl[3], ch[3]; // cell bounds
  int itype, sval, xyzcell, cxyz[3];
  for(int icell = 0; icell < grid->nlocal; icell++) {

    cl[0] = cells[icell].lo[0];
    cl[1] = cells[icell].lo[1];
    cl[2] = cells[icell].lo[2];
    xyzcell = get_cxyz(cxyz,cl);

    // itype = 1 - fully outside
    // itype = 2 - fully inside
    // itype = 3 - has surfaces
    // cannot just use types, if surface on corner
    // .. can be either 2 or 3 depending on which corner (ambiguity)
    itype = cinfo[icell].type;

    // fully inside so set all corner values to max
    if(itype==2) {
      sval = 1;
    // fully outside so set all corners to min
    } else if (itype==1) {
      sval = 0;
    } else {
      continue;
    }

    // set corners if not already set
    if(svalues[xyzcell] < 0) svalues[xyzcell] = sval;
    xyzcell = get_corner(cxyz[0]+1, cxyz[1], cxyz[2]);
    if(svalues[xyzcell] < 0) svalues[xyzcell] = sval;
    xyzcell = get_corner(cxyz[0], cxyz[1]+1, cxyz[2]);
    if(svalues[xyzcell] < 0) svalues[xyzcell] = sval;
    xyzcell = get_corner(cxyz[0]+1, cxyz[1]+1, cxyz[2]);
    if(svalues[xyzcell] < 0) svalues[xyzcell] = sval;

    if(domain->dimension==3) {
      xyzcell = get_corner(cxyz[0], cxyz[1], cxyz[2]+1);
      if(svalues[xyzcell] < 0) svalues[xyzcell] = sval;
      xyzcell = get_corner(cxyz[0]+1, cxyz[1], cxyz[2]+1);
      if(svalues[xyzcell] < 0) svalues[xyzcell] = sval;
      xyzcell = get_corner(cxyz[0], cxyz[1]+1, cxyz[2]+1);
      if(svalues[xyzcell] < 0) svalues[xyzcell] = sval;
      xyzcell = get_corner(cxyz[0]+1, cxyz[1]+1, cxyz[2]+1);
      if(svalues[xyzcell] < 0) svalues[xyzcell] = sval;
    }
  }
  return;
}

/* ----------------------------------------------------------------------
  Resolve unknown side values and fill in remaining corner values
------------------------------------------------------------------------- */

void CreateISurf::cleanup()
{
  int j;
  int filled, attempt;
  filled = attempt = 0;
  while(filled==0) {
    filled=1;
    for(int i = 0; i < Nxyz; i++) {
      if(svalues[i] == 0 || svalues[i] == 1) continue;

      // try x-neighbor
      j = i - 1;
      if(j >= 0 && (svalues[j] == 0 || svalues[j] == 1)) {
        if(ivalues[i][0] <= 0) svalues[i] = svalues[j];
        else {
          if(svalues[j] == 0) svalues[i] = 1;
          else svalues[i] = 0;
        }
        continue;
      } 

      j = i + 1;
      if(j < Nxyz && (svalues[j] == 0 || svalues[j] == 1)) {
        if(ivalues[i][1] <= 0) svalues[i] = svalues[j];
        else {
          if(svalues[j] == 0) svalues[i] = 1;
          else svalues[i] = 0;
        }
        continue;
      }

      // try y-neighbor
      j = i - (nxyz[0]+1);
      if(j >= 0 && (svalues[j] == 0 || svalues[j] == 1)) {
        if(ivalues[i][2] <= 0) svalues[i] = svalues[j];
        else {
          if(svalues[j] == 0) svalues[i] = 1;
          else svalues[i] = 0;
        }
        continue;
      }

      j = i + (nxyz[0]+1);
      if(j < Nxyz && (svalues[j] == 0 || svalues[j] == 1)) {
        if(ivalues[i][3] <= 0) svalues[i] = svalues[j];
        else {
          if(svalues[j] == 0) svalues[i] = 1;
          else svalues[i] = 0;
        }
        continue;
      }

      if(domain->dimension=3) {
        // try z-neighbor
        j = i - (nxyz[0]+1)*(nxyz[1]+1);
        if(j >= 0 && (svalues[j] == 0 || svalues[j] == 1)) {
          if(ivalues[i][4] <= 0) svalues[i] = svalues[j];
          else {
            if(svalues[j] == 0) svalues[i] = 1;
            else svalues[i] = 0;
          }
          continue;
        }

        j = i + (nxyz[0]+1)*(nxyz[1]+1);
        if(j < Nxyz && (svalues[j] == 0 || svalues[j] == 1)) {
          if(ivalues[i][5] <= 0) svalues[i] = svalues[j];
          else {
            if(svalues[j] == 0) svalues[i] = 1;
            else svalues[i] = 0;
          }
          continue;
        }
      }

      /*
      // try x-neighbor
      j = i - 1;
      if(j >= 0 && ivalues[i][0] < 0) {
        svalues[i] = svalues[j];
        continue;
      }

      j = i + 1;
      if(j < Nxyz && ivalues[i][1] < 0) {
        svalues[i] = svalues[j];
        continue;
      }

      // try y-neighbor
      j = i - (nxyz[0]+1);
      if(j >= 0 && ivalues[i][2] < 0) {
        svalues[i] = svalues[j];
        continue;
      }

      j = i + (nxyz[0]+1);
      if(j < Nxyz && ivalues[i][3] < 0) {
        svalues[i] = svalues[j];
        continue;
      }

      if(domain->dimension=3) {
        // try z-neighbor
        j = i - (nxyz[0]+1)*(nxyz[1]+1);
        if(j >= 0 && ivalues[i][4] < 0) {
          svalues[i] = svalues[j];
          continue;
        }

        j = i + (nxyz[0]+1)*(nxyz[1]+1);
        if(j < Nxyz && ivalues[i][5] < 0) {
          svalues[i] = svalues[j];
          continue;
        }
      }
      */

    }

    for(int i = 0; i < Nxyz; i++)
      if(svalues[i] < 0) filled = 0;

    attempt++;
    if(attempt>20) error->one(FLERR,"Cannot find reference");

  }

  for(int i = 0; i < Nxyz; i++)
    if(svalues[i] < 0) error->one(FLERR,"bad sval");

  int npairs;
  if(domain->dimension==2) npairs = 4;
  else npairs = 6;

  // initially set inside as cin and outside as cout
  if(aveFlag) {
    int nval;
    double ivalsum;
    for(int i = 0; i < Nxyz; i++) {
      ivalsum = 0.0;
      if(svalues[i] == 0) cvalues[i] = cout;
      else {
        nval = 0;
        for(int j = 0; j < npairs; j++) {
          if(ivalues[i][j] >= 0) {
            ivalsum += ivalues[i][j];
            nval++;
          }
        }

        // no intersections
        if(nval == 0) cvalues[i] = cin;
        else {
          ivalsum /= nval;
          if(ivalsum > 1.0) error->one(FLERR,"over 1");
          // add small buffer to avoid very small volumes/areas
          /*if(domain->dimension == 3) {
            ivalsum = MAX(ivalsum, 0.045);
            ivalsum = MIN(ivalsum, 0.955);
          } else {
            ivalsum = MAX(ivalsum, 0.01);
            ivalsum = MIN(ivalsum, 0.99);
          }*/
          cvalues[i] = param2in(ivalsum,0.0);
        }
      }
    }// end "for" for grid cells

  } else {
    for(int i = 0; i < Nxyz; i++) {
      if(svalues[i] == 0) cvalues[i] = cout;
      else if(svalues[i] == 1) cvalues[i] = cin;
      else error->one(FLERR,"bad svalues");
    }// end "for" for grid cells
  }
}

/* ----------------------------------------------------------------------
   Determines if surface is inline with cell edge.

   Side = 0,1 -> [out,in]
------------------------------------------------------------------------- */

int CreateISurf::corner_hit2d(double *p1, double *p2,
    Surf::Line *line, double &param, int &side)
{
  // try intersect first
  int h, tside;
  double tparam;
  double d3dum[3];
  h = Geometry::
    line_line_intersect(p1,p2,line->p1,line->p2,line->norm,d3dum,tparam,tside);
  if(h) {
    if(tside == 1 || tside == 2 || tside == 5) {
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

  for(int i = 0; i < 8; i++) {
    p1p[2] = p1[2];
    p2p[2] = p2[2];

    p1p[0] = p1[0] + dx[i];
    p1p[1] = p1[1] + dy[i];

    p2p[0] = p2[0] + dx[i];
    p2p[1] = p2[1] + dy[i];

    h = Geometry::
      line_line_intersect(p1p,p2p,line->p1,line->p2,line->norm,d3dum,tparam,tside);
    if(h) {
      if(tside == 1 || tside == 2 || tside == 5) {
        //side = 2;
        side = 1;
        if(tparam<0.5) param = 0.0;
        else param = 1.0;
        return true;
      } else {
        //side = 2;
        side = 0;
        if(tparam<0.5) param = 0.0;
        else param = 1.0;
        return true;
      }
    }
  }

  // true miss
  return false;
}

/* ----------------------------------------------------------------------
   Same as corner_hit2d but for 3d
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

  if(h) {
    if(tside == 1 || tside == 2 || tside == 5) {
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

  for(int i = 0; i < 26; i++) {
    p1p[0] = p1[0] + dx[i];
    p1p[1] = p1[1] + dy[i];
    p1p[2] = p1[2] + dz[i];

    p2p[0] = p2[0] + dx[i];
    p2p[1] = p2[1] + dy[i];
    p2p[2] = p2[2] + dz[i];
 
    h = Geometry::line_tri_intersect(p1p,p2p,tri->p1,tri->p2,tri->p3,
        tri->norm,d3dum,tparam,tside);
    if(h) {
      if(tside == 1 || tside == 2 || tside == 5) {
        side = 1;
        if(tparam<0.5) param = 0.0;
        else param = 1.0;
        return true;
      } else {
        side = 0;
        if(tparam<0.5) param = 0.0;
        else param = 1.0;
        return true;
      }
    }
  }

  // true miss
  return false;
}


/* ----------------------------------------------------------------------
   Find inside corner value from corner value
------------------------------------------------------------------------- */
double CreateISurf::param2in(double param, double v1)
{
  double v0;
  // param is proportional to cell length so 
  // ... lo = 0; hi = 1
  // trying to find v0
  // param = (thresh  - v0) / (v1 - v0)
  if(param == 1.0) return 255.0;
  v0 = (thresh - v1*param) / (1.0 - param);

  // bound by limits
  //v0 = MAX(v0,thresh);
  v0 = MIN(v0,255.0);
  return v0;
}

/* ----------------------------------------------------------------------
   Determines corner values given an explicit surface. Cvalues then used
   later in ablate to create the implicit surfaces
------------------------------------------------------------------------- */

int CreateISurf::get_cxyz(int *ic, double *lc)
{

  // find lower bounds
  double lo[3];
  lo[0] = domain->boxlo[0];
  lo[1] = domain->boxlo[1];
  lo[2] = domain->boxlo[2];

  // shift by lower bounds
  double lclo[3];
  for(int d = 0; d < 3; d++) {
    lclo[d] = lc[d]-lo[d];
    ic[d] = static_cast<int> (lclo[d] / xyzsize[d] + 0.5);
  }

  int icell = get_corner(ic[0], ic[1], ic[2]);

  return icell;
}

/* ----------------------------------------------------------------------
   Get cell from coordinates
------------------------------------------------------------------------- */

int CreateISurf::get_cell(int icx, int icy, int icz)
{
  int icell;
  if(domain->dimension == 2) icell = icx + icy*nxyz[0];
  else icell = icx + icy*nxyz[0] + icz*nxyz[0]*nxyz[1];

  if(icell >= nxyz[0]*nxyz[1]*nxyz[2] || icell < 0) error->one(FLERR,"bad cell from int");

  return icell;
}

/* ----------------------------------------------------------------------
   Get corner values from coordinates
------------------------------------------------------------------------- */

int CreateISurf::get_corner(int icx, int icy, int icz)
{
  int icell;
  if(domain->dimension == 2) icell = icx + icy*(nxyz[0]+1);
  else icell = icx + icy*(nxyz[0]+1) + icz*(nxyz[0]+1)*(nxyz[1]+1);

  if(icell >= Nxyz || icell < 0) {
    printf("icell: %i\n", icell);
    error->one(FLERR,"bad corner from int");
  }

  return icell;
}

/* ----------------------------------------------------------------------
   Get corner values from coordinates
------------------------------------------------------------------------- */

int CreateISurf::get_corner(double dcx, double dcy, double dcz)
{
  // find lower bounds
  double lo[3];
  lo[0] = domain->boxlo[0];
  lo[1] = domain->boxlo[1];
  lo[2] = domain->boxlo[2];

  // shift by lower bounds
  double lclo[3];
  int ic[3];
  double lc[3];
  lc[0] = dcx;
  lc[1] = dcy; 
  lc[2] = dcz;
  for(int d = 0; d < 3; d++) {
    lclo[d] = lc[d]-lo[d];
    ic[d] = static_cast<int> (lclo[d] / xyzsize[d] + 0.5);
  }

  int icell;
  if(domain->dimension == 2) icell = ic[0] + ic[1]*(nxyz[0]+1);
  else icell = ic[0] + ic[1]*(nxyz[0]+1) + ic[2]*(nxyz[0]+1)*(nxyz[1]+1);

  if(icell >= Nxyz || icell < 0) {
    printf("dc: %4.3e, %4.3e, %4.3e\n", dcx, dcy, dcz);
    printf("icell: %i\n", icell);
    error->one(FLERR,"bad corner from double");
  }

  return icell;
}

/* ----------------------------------------------------------------------
   Removes old explicit surfaces
------------------------------------------------------------------------- */
void CreateISurf::remove_old()
{
  // copied from remove_surf.cpp
  if (me == 0)
    if (screen) fprintf(screen,"Removing explicit surfs ...\n");

  if (particle->exist) particle->sort();
  MPI_Barrier(world);

  // finish implementing custom values
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
  if (domain->dimension == 2) nbytes = sizeof(Surf::Line);
  else nbytes = sizeof(Surf::Tri);

  if (domain->dimension == 2) {
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

  if (domain->dimension == 2) surf->check_watertight_2d();
  else surf->check_watertight_3d();

  MPI_Barrier(world);

  // reset grid due to changing surfs
  // assign surfs to grid cells

  surf->setup_owned();
  grid->unset_neighbors();
  grid->remove_ghosts();

  // reassign split cell particles to parent split cell

  if (particle->exist && grid->nsplitlocal) {
    Grid::ChildCell *cells = grid->cells;
    int nglocal = grid->nlocal;
    for (int icell = 0; icell < nglocal; icell++)
      if (cells[icell].nsplit > 1)
        grid->combine_split_cell_particles(icell,1);
  }

  grid->clear_surf();
  //if (surf->exist) grid->surf2grid(1); // surf shouldn't exist

  // reassign particles in split cells to sub cell owner

  if (particle->exist && grid->nsplitlocal) {
    Grid::ChildCell *cells = grid->cells;
    int nglocal = grid->nlocal;
    for (int icell = 0; icell < nglocal; icell++)
      if (cells[icell].nsplit > 1)
        grid->assign_split_cell_particles(icell);
  }

  MPI_Barrier(world);
  double time4 = MPI_Wtime();

  // re-setup owned and ghost cell info

  grid->setup_owned();
  grid->acquire_ghosts();
  grid->reset_neighbors();
  comm->reset_neighbors();

  grid->set_inout();
  grid->type_check();

  MPI_Barrier(world);

  if (me == 0)
    if (screen) fprintf(screen,"Finished deleting old explicit surfaces\n");
}

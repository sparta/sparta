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

#include "marching_squares_dev.h"
#include "grid.h"
#include "surf.h"
#include "error.h"

using namespace SPARTA_NS;

enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};           // several files

/* ---------------------------------------------------------------------- */

MarchingSquaresDev::MarchingSquaresDev(SPARTA *sparta, int ggroup_caller,
                                 double thresh_caller) :
  Pointers(sparta)
{
  ggroup = ggroup_caller;
  thresh = thresh_caller;
}

/* ----------------------------------------------------------------------
   create 2d implicit surfs from grid point values
   follows https://en.wikipedia.org/wiki/Marching_squares
   see 2 sections: Basic algorithm and Disambiguation of saddle points
     treating open circles as flow volume, solid circles as material
     NOTE: Wiki page numbers points counter-clockwise
           SPARTA numbers them in x, then in y
           so bit2 and bit3 are swapped below
           this gives case #s here consistent with Wiki page
   process each grid cell independently
   4 corner points open/solid -> 2^4 = 16 cases
   cases infer 0,1,2 line segments in each grid cell
   order 2 points in each line segment to give normal into flow volume
   treat two saddle point cases (my 9,6) (Wiki 5,10)
     based on ave value at cell center
------------------------------------------------------------------------- */

void MarchingSquaresDev::invoke(double ***cvalues, int *svalues)
{
  int i,ipt,isurf,nsurf,which;
  double v00,v01,v10,v11;
  double i00,i01,i10,i11; // intersection values
  int bit0,bit1,bit2,bit3;
  double ave;
  double *lo,*hi;
  surfint surfID;
  surfint *ptr;

  double pt[4][3];
  pt[0][2] = pt[1][2] = pt[2][2] = pt[3][2] = 0.0;

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  MyPage<surfint> *csurfs = grid->csurfs;
  int nglocal = grid->nlocal;
  int groupbit = grid->bitmask[ggroup];

  bigint maxsurfID = 0;
  if (sizeof(surfint) == 4) maxsurfID = MAXSMALLINT;
  if (sizeof(surfint) == 8) maxsurfID = MAXBIGINT;

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;
    lo = cells[icell].lo;
    hi = cells[icell].hi;

    printf("icell - ms: %i\n", icell);

    // cvalues are ordered lower-left, lower-right, upper-left, upper-right
    // Vyx encodes this as 0/1 in each dim

    v00 = cvalues[icell][0][0];
    v01 = cvalues[icell][1][0];
    v10 = cvalues[icell][2][0];
    v11 = cvalues[icell][3][0];

    // set corner point to maxval of adjacent values
    for (i = 1; i < 4; i++) {
      if(v00 < thresh) v00 = MAX(v00, cvalues[icell][0][i]);
      if(v01 < thresh) v01 = MAX(v01, cvalues[icell][1][i]);
      if(v10 < thresh) v10 = MAX(v10, cvalues[icell][2][i]);
      if(v11 < thresh) v11 = MAX(v11, cvalues[icell][3][i]);
    }

    // intersections on the cell edges
    // notice that any edge will only have one intersection

    i00 = i01 = i10 = i11 = -1.0;

    if(cvalues[icell][1][0] >= thresh || v00 >= thresh)
      i00 = interpolate(v00,cvalues[icell][1][0],lo[0],hi[0]);

    if(cvalues[icell][3][2] >= thresh || v01 >= thresh)
      i01 = interpolate(v01,cvalues[icell][3][2],lo[1],hi[1]);

    if(cvalues[icell][2][1] >= thresh || v11 >= thresh)
      i10 = interpolate(v11,cvalues[icell][2][1],lo[0],hi[0]);

    if(cvalues[icell][0][3] >= thresh || v10 >= thresh)
      i11 = interpolate(v10,cvalues[icell][0][3],lo[1],hi[1]);

    // make last 2 bits consistent with Wiki page (see NOTE above)

    bit0 = v00 <= thresh ? 0 : 1;
    bit1 = v01 <= thresh ? 0 : 1;
    bit2 = v11 <= thresh ? 0 : 1;
    bit3 = v10 <= thresh ? 0 : 1;

    // TODO: Cases 5 and 10 may not be correct since more than one surface
    // Need to use minimum

    which = (bit3 << 3) + (bit2 << 2) + (bit1 << 1) + bit0;

    switch (which) {

    case 0:
      nsurf = 0;
      break;

    case 1:
      nsurf = 1;
      pt[0][0] = lo[0];
      //pt[0][1] = interpolate(v00,v10,lo[1],hi[1]);
      //pt[1][0] = interpolate(v00,v01,lo[0],hi[0]);
      pt[0][1] = i11;
      pt[1][0] = i00;
      pt[1][1] = lo[1];
      break;

    case 2:
      nsurf = 1;
      //pt[0][0] = interpolate(v00,v01,lo[0],hi[0]);
      pt[0][0] = i00;
      pt[0][1] = lo[1];
      pt[1][0] = hi[0];
      //pt[1][1] = interpolate(v01,v11,lo[1],hi[1]);
      pt[1][1] = i01;
      break;

    case 3:
      nsurf = 1;
      pt[0][0] = lo[0];
      //pt[0][1] = interpolate(v00,v10,lo[1],hi[1]);
      pt[0][1] = i11;
      pt[1][0] = hi[0];
      //pt[1][1] = interpolate(v01,v11,lo[1],hi[1]);
      pt[1][1] = i01;
      break;

    case 4:
      nsurf = 1;
      pt[0][0] = hi[0];
      //pt[0][1] = interpolate(v01,v11,lo[1],hi[1]);
      pt[0][1] = i01;
      //pt[1][0] = interpolate(v10,v11,lo[0],hi[0]);
      pt[1][0] = i10;
      pt[1][1] = hi[1];
      break;

    case 5:
      nsurf = 2;
      ave = 0.25 * (v00 + v01 + v10 + v11);
      if (ave > thresh) {
        pt[0][0] = lo[0];
        //pt[0][1] = interpolate(v00,v10,lo[1],hi[1]);
        //pt[1][0] = interpolate(v10,v11,lo[0],hi[0]);
        pt[0][1] = i11;
        pt[1][0] = i10;
        pt[1][1] = hi[1];
        pt[2][0] = hi[0];
        pt[2][1] = i11;
        pt[3][0] = i10;
        pt[3][1] = lo[1];
      } else {
        pt[0][0] = lo[0];
        //pt[0][1] = interpolate(v00,v10,lo[1],hi[1]);
        //pt[1][0] = interpolate(v00,v01,lo[0],hi[0]);
        pt[0][1] = i11;
        pt[1][0] = i00;
        pt[1][1] = lo[1];
        pt[2][0] = hi[0];
        //pt[2][1] = interpolate(v01,v11,lo[1],hi[1]);
        //pt[3][0] = interpolate(v10,v11,lo[0],hi[0]);
        pt[2][1] = i01;
        pt[3][0] = i10;
        pt[3][1] = hi[1];
      }
      break;

    case 6:
      nsurf = 1;
      //pt[0][0] = interpolate(v00,v01,lo[0],hi[0]);
      pt[0][0] = i00;
      pt[0][1] = lo[1];
      //pt[1][0] = interpolate(v10,v11,lo[0],hi[0]);
      pt[1][0] = i10;
      pt[1][1] = hi[1];
      break;

    case 7:
      nsurf = 1;
      pt[0][0] = lo[0];
      //pt[0][1] = interpolate(v00,v10,lo[1],hi[1]);
      pt[0][1] = i11;
      //pt[1][0] = interpolate(v10,v11,lo[0],hi[0]);
      pt[1][0] = i10;
      pt[1][1] = hi[1];
      break;

    case 8:
      nsurf = 1;
      //pt[0][0] = interpolate(v10,v11,lo[0],hi[0]);
      pt[0][0] = i10;
      pt[0][1] = hi[1];
      pt[1][0] = lo[0];
      //pt[1][1] = interpolate(v00,v10,lo[1],hi[1]);
      pt[1][1] = i11;
      break;

    case 9:
      nsurf = 1;
      //pt[0][0] = interpolate(v10,v11,lo[0],hi[0]);
      pt[0][0] = i10;
      pt[0][1] = hi[1];
      //pt[1][0] = interpolate(v00,v01,lo[0],hi[0]);
      pt[1][0] = i00;
      pt[1][1] = lo[1];
      break;

    case 10:
      nsurf = 2;
      ave = 0.25 * (v00 + v01 + v10 + v11);
      if (ave > thresh) {
        //pt[0][0] = interpolate(v00,v01,lo[0],hi[0]);
        pt[0][0] = i00;
        pt[0][1] = lo[1];
        pt[1][0] = lo[0];
        //pt[1][1] = interpolate(v00,v10,lo[1],hi[1]);
        pt[1][1] = i11;
        //pt[2][0] = interpolate(v10,v11,lo[0],hi[0]);
        pt[2][0] = i10;
        pt[2][1] = hi[1];
        pt[3][0] = hi[0];
        //pt[3][1] = interpolate(v01,v11,lo[1],hi[1]);
        pt[3][1] = i01;
      } else {
        //pt[0][0] = interpolate(v10,v11,lo[0],hi[0]);
        pt[0][0] = i10;
        pt[0][1] = hi[1];
        pt[1][0] = lo[0];
        //pt[1][1] = interpolate(v00,v10,lo[1],hi[1]);
        pt[1][1] = i11;
        //pt[2][0] = interpolate(v00,v01,lo[0],hi[0]);
        pt[2][0] = i00;
        pt[2][1] = lo[1];
        pt[3][0] = hi[0];
        //pt[3][1] = interpolate(v01,v11,lo[1],hi[1]);
        pt[3][1] = i01;
      }
      break;

    case 11:
      nsurf = 1;
      //pt[0][0] = interpolate(v10,v11,lo[0],hi[0]);
      pt[0][0] = i10;
      pt[0][1] = hi[1];
      pt[1][0] = hi[0];
      //pt[1][1] = interpolate(v01,v11,lo[1],hi[1]);
      pt[1][1] = i01;
      break;

    case 12:
      nsurf = 1;
      pt[0][0] = hi[0];
      //pt[0][1] = interpolate(v01,v11,lo[1],hi[1]);
      pt[0][1] = i01;
      pt[1][0] = lo[0];
      //pt[1][1] = interpolate(v00,v10,lo[1],hi[1]);
      pt[1][1] = i11;
      break;

    case 13:
      nsurf = 1;
      pt[0][0] = hi[0];
      //pt[0][1] = interpolate(v01,v11,lo[1],hi[1]);
      pt[0][1] = i01;
      //pt[1][0] = interpolate(v00,v01,lo[0],hi[0]);
      pt[1][0] = i00;
      pt[1][1] = lo[1];
      break;

    case 14:
      nsurf = 1;
      //pt[0][0] = interpolate(v00,v01,lo[0],hi[0]);
      pt[0][0] = i00;
      pt[0][1] = lo[1];
      pt[1][0] = lo[0];
      //pt[1][1] = interpolate(v00,v10,lo[1],hi[1]);
      pt[1][1] = i11;
      break;

    case 15:
      nsurf = 0;
      break;
    }

    /*printf("v: %4.3e,%4.3e,%4.3e,%4.3e\n",
      v00,v01,v10,v11);
    printf("which: %i; i: %4.3e,%4.3e,%4.3e,%4.3e\n",
      which, i00,i01,i10,i11);
    printf("cvalues[0]: %4.3e,%4.3e,%4.3e,%4.3e\n",
      cvalues[icell][0][0], cvalues[icell][0][1],
      cvalues[icell][0][2], cvalues[icell][0][3]);
    printf("cvalues[1]: %4.3e,%4.3e,%4.3e,%4.3e\n",
      cvalues[icell][1][0], cvalues[icell][1][1],
      cvalues[icell][1][2], cvalues[icell][1][3]);
    printf("cvalues[2]: %4.3e,%4.3e,%4.3e,%4.3e\n",
      cvalues[icell][2][0], cvalues[icell][2][1],
      cvalues[icell][2][2], cvalues[icell][2][3]);
    printf("cvalues[3]: %4.3e,%4.3e,%4.3e,%4.3e\n",
      cvalues[icell][3][0], cvalues[icell][3][1],
      cvalues[icell][3][2], cvalues[icell][3][3]);*/

    printf("pt-0: %4.3e, %4.3e\n", pt[0][0], pt[0][1]);
    printf("pt-1: %4.3e, %4.3e\n", pt[1][0], pt[1][1]);

    // DEBUG: Check for bad end points
    if(pt[0][0] < lo[0] || pt[0][1] < lo[1] ||
       pt[1][0] < lo[0] || pt[1][1] < lo[1])
        error->one(FLERR,"Outside low cell");

    if(pt[0][0] > hi[0] || pt[0][1] > hi[1] ||
       pt[1][0] > hi[0] || pt[1][1] > hi[1])
        error->one(FLERR,"Outside hi cell");

    if(which == 5 || which == 10) {
      if(pt[2][0] < lo[0] || pt[2][1] < lo[1] ||
         pt[3][0] < lo[0] || pt[3][1] < lo[1])
      error->one(FLERR,"Outside low cell");

      if(pt[2][0] > hi[0] || pt[2][1] > hi[1] ||
         pt[3][0] > hi[0] || pt[3][1] > hi[1])
          error->one(FLERR,"Outside hi cell");
    }

    // populate Grid and Surf data structs
    // points will be duplicated, not unique
    // surf ID = cell ID for all surfs in cell
    // check if uint cell ID overflows int surf ID

    if (nsurf) {
      if (cells[icell].id > maxsurfID)
        error->one(FLERR,"Grid cell ID overflows implicit surf ID");
      surfID = cells[icell].id;
    }

    ptr = csurfs->get(nsurf);

    ipt = 0;
    for (i = 0; i < nsurf; i++) {
      if (svalues) surf->add_line(surfID,svalues[icell],pt[ipt],pt[ipt+1]);
      else surf->add_line(surfID,1,pt[ipt],pt[ipt+1]);
      ipt += 2;
      isurf = surf->nlocal - 1;
      ptr[i] = isurf;
    }

    cells[icell].nsurf = nsurf;
    if (nsurf) {
      cells[icell].csurfs = ptr;
      cinfo[icell].type = OVERLAP;
    }
  }
}

/* ----------------------------------------------------------------------
   interpolate function used by both marching squares and cubes
   lo/hi = coordinates of end points of edge of square
   v0/v1 = values at lo/hi end points
   value = interpolated coordinate for thresh value
------------------------------------------------------------------------- */

double MarchingSquaresDev::interpolate(double v0, double v1, double lo, double hi)
{
  double value = lo + (hi-lo)*(thresh-v0)/(v1-v0);
  value = MAX(value,lo);
  value = MIN(value,hi);
  return value;
}

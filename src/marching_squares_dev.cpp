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

void MarchingSquaresDev::invoke(double ***cvalues, double ***ivalues, int *svalues)
{
  int i,ipt,isurf,nsurf,which;
  double v00,v01,v10,v11;
  double i0, i1, i2, i3;
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

    // cvalues are ordered lower-left, lower-right, upper-left, upper-right
    // Vyx encodes this as 0/1 in each dim

    v00 = v01 = v10 = v11 = 0.0;

    // set corner point to average of adjacent values

    for (i = 0; i < 4; i++) {
      v00 += cvalues[icell][0][i];
      v01 += cvalues[icell][1][i];
      v10 += cvalues[icell][2][i];
      v11 += cvalues[icell][3][i];
    }

    v00 /= 4;
    v01 /= 4;
    v10 /= 4;
    v11 /= 4;

    if(v00 > thresh) v00 = 255.0;
    else v00 = 0.0;
    if(v01 > thresh) v01 = 255.0;
    else v01 = 0.0;
    if(v10 > thresh) v10 = 255.0;
    else v10 = 0.0;
    if(v11 > thresh) v11 = 255.0;
    else v11 = 0.0;

    // intersection of surfaces on cell edges

    i0 = ivalues[icell][0][1];
    i1 = ivalues[icell][1][3];
    i2 = ivalues[icell][2][1];
    i3 = ivalues[icell][0][3];

    if (i0<0) i0 = lo[0];
    else i0 = lo[0] + (hi[0]-lo[0])*i0;
    if (i1<0) i1 = lo[1];
    else i1 = lo[1] + (hi[1]-lo[1])*i1;
    if (i2<0) i2 = lo[0];
    else i2 = lo[0] + (hi[0]-lo[0])*i2;
    if (i3<0) i3 = lo[1];
    else i3 = lo[1] + (hi[1]-lo[1])*i3;

    //printf("i: %4.3e, %4.3e, %4.3e, %4.3e\n", i0, i1, i2, i3);
    //printf("l: %4.3e, %4.3e, %4.3e, %4.3e\n", lo[0],lo[1], lo[0], lo[1]);
    //printf("u: %4.3e, %4.3e, %4.3e, %4.3e\n", hi[0],hi[1], hi[0], hi[1]);

    if(i0 < lo[0] || i0 > hi[0]) error->one(FLERR,"out");
    if(i1 < lo[1] || i1 > hi[1]) error->one(FLERR,"out");
    if(i2 < lo[0] || i2 > hi[0]) error->one(FLERR,"out");
    if(i3 < lo[1] || i3 > hi[1]) error->one(FLERR,"out");

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
      pt[0][1] = i3;
      pt[1][0] = i0;
      pt[1][1] = lo[1];
      break;

    case 2:
      nsurf = 1;
      pt[0][0] = i0;
      pt[0][1] = lo[1];
      pt[1][0] = hi[0];
      pt[1][1] = i1;
      break;

    case 3:
      nsurf = 1;
      pt[0][0] = lo[0];
      pt[0][1] = i3;
      pt[1][0] = hi[0];
      pt[1][1] = i1;
      break;

    case 4:
      nsurf = 1;
      pt[0][0] = hi[0];
      pt[0][1] = i1;
      pt[1][0] = i2;
      pt[1][1] = hi[1];
      break;

    case 5:
      nsurf = 2;
      ave = 0.25 * (v00 + v01 + v10 + v11);
      if (ave > thresh) {
       pt[0][0] = lo[0];
        pt[0][1] = i3;
        pt[1][0] = i2;
        pt[1][1] = hi[1];
        pt[2][0] = hi[0];
        pt[2][1] = i1;
        pt[3][0] = i0;
        pt[3][1] = lo[1];
      } else {
        pt[0][0] = lo[0];
        pt[0][1] = i3;
        pt[1][0] = i0;
        pt[1][1] = lo[1];
        pt[2][0] = hi[0];
        pt[2][1] = i1;
        pt[3][0] = i2;
        pt[3][1] = hi[1];
      }
      break;

    case 6:
      nsurf = 1;
      pt[0][0] = i0;
      pt[0][1] = lo[1];
      pt[1][0] = i2;
      pt[1][1] = hi[1];
      break;

    case 7:
      nsurf = 1;
      pt[0][0] = lo[0];
      pt[0][1] = i3;
      pt[1][0] = i2;
      pt[1][1] = hi[1];
      break;

    case 8:
      nsurf = 1;
      pt[0][0] = i2;
      pt[0][1] = hi[1];
      pt[1][0] = lo[0];
      pt[1][1] = i3;
      break;

    case 9:
      nsurf = 1;
      pt[0][0] = i2;
      pt[0][1] = hi[1];
      pt[1][0] = i0;
      pt[1][1] = lo[1];
      break;

    case 10:
      nsurf = 2;
      ave = 0.25 * (v00 + v01 + v10 + v11);
      if (ave > thresh) {
        pt[0][0] = i0;
        pt[0][1] = lo[1];
        pt[1][0] = lo[0];
        pt[1][1] = i3;
        pt[2][0] = i2;
        pt[2][1] = hi[1];
        pt[3][0] = hi[0];
        pt[3][1] = i1;
      } else {
        pt[0][0] = i2;
        pt[0][1] = hi[1];
        pt[1][0] = lo[0];
        pt[1][1] = i3;
        pt[2][0] = i0;
        pt[2][1] = lo[1];
        pt[3][0] = hi[0];
        pt[3][1] = i1;
      }
      break;

    case 11:
      nsurf = 1;
      pt[0][0] = i2;
      pt[0][1] = hi[1];
      pt[1][0] = hi[0];
      pt[1][1] = i1;
      break;

    case 12:
      nsurf = 1;
      pt[0][0] = hi[0];
      pt[0][1] = i1;
      pt[1][0] = lo[0];
      pt[1][1] = i3;
      break;

    case 13:
      nsurf = 1;
      pt[0][0] = hi[0];
      pt[0][1] = i1;
      pt[1][0] = i0;
      pt[1][1] = lo[1];
      break;

    case 14:
      nsurf = 1;
      pt[0][0] = i0;
      pt[0][1] = lo[1];
      pt[1][0] = lo[0];
      pt[1][1] = i3;
      break;

    case 15:
      nsurf = 0;
      break;
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

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

#include "math.h"
#include "math_extra.h"
#include "string.h"
#include "marching_cubes.h"
#include "grid.h"
#include "surf.h"
#include "irregular.h"
#include "lookup_table.h"
#include "geometry.h"
#include "my_page.h"
#include "memory.h"
#include "error.h"

// DEBUG
#include "update.h"

using namespace SPARTA_NS;

// prototype for non-class function

int compare_indices(const void *, const void *);

enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};           // several files
enum{NCHILD,NPARENT,NUNKNOWN,NPBCHILD,NPBPARENT,NPBUNKNOWN,NBOUND};  // Grid

#define DELTA 128
#define BIG 1.0e20
#define EPSILON 1.0e-16

/* ----------------------------------------------------------------------
   Same as above but uses inner values. Also if there are ambiguities,
   the corner values corresponding to the intersections are first found
   then the ambiguity tests are performed
------------------------------------------------------------------------- */

void MarchingCubes::invoke(double ***cvalues, int *svalues, int **mcflags)
{
  int i,j,ipt,isurf,nsurf,icase,which;
  surfint surfID;
  surfint *ptr;

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

    // nsurf = # of tris in cell
    // cvalues[8] = 8 corner point values, each is 0 to 255 inclusive
    // thresh = value between 0 and 255 to threshhold on
    // lo[3] = lower left corner pt of grid cell
    // hi[3] = upper right corner pt of grid cell
    // pt = list of 3*nsurf points that are the corner pts of each tri

    // cvalues are ordered
    // bottom-lower-left, bottom-lower-right,
    // bottom-upper-left, bottom-upper-right
    // top-lower-left, top-lower-right, top-upper-left, top-upper-right
    // Vzyx encodes this as 0/1 in each dim

    // temporarily store all inner values

    for (i = 0; i < 8; i++)
      for (j = 0; j < 6; j++)
        inval[i][j] = cvalues[icell][i][j];

    // use averages for now

    v000 = v001 = v010 = v011 = 0.0;
    v100 = v101 = v110 = v111 = 0.0;

    // ordering in cvalues different from loop up table
    // manually change for consistency

    for (j = 0; j < 6; j++) {
      v000 += inval[0][j];
      v001 += inval[1][j];
      v010 += inval[2][j];
      v011 += inval[3][j];
      v100 += inval[4][j];
      v101 += inval[5][j];
      v110 += inval[6][j];
      v111 += inval[7][j];
    }

    v000 /= 6.0;
    v001 /= 6.0;
    v010 /= 6.0;
    v011 /= 6.0;
    v100 /= 6.0;
    v101 /= 6.0;
    v110 /= 6.0;
    v111 /= 6.0;

    v000iso = v000 - thresh;
    v001iso = v001 - thresh;
    v010iso = v010 - thresh;
    v011iso = v011 - thresh;
    v100iso = v100 - thresh;
    v101iso = v101 - thresh;
    v110iso = v110 - thresh;
    v111iso = v111 - thresh;

    // intersection of surfaces on all cell edges

    i0  = interpolate(inval[0][1],inval[1][0],lo[0],hi[0]);
    i1  = interpolate(inval[1][3],inval[3][2],lo[1],hi[1]);
    i2  = interpolate(inval[2][1],inval[3][0],lo[0],hi[0]);
    i3  = interpolate(inval[0][3],inval[2][2],lo[1],hi[1]);

    i4  = interpolate(inval[4][1],inval[5][0],lo[0],hi[0]);
    i5  = interpolate(inval[5][3],inval[7][2],lo[1],hi[1]);
    i6  = interpolate(inval[6][1],inval[7][0],lo[0],hi[0]);
    i7  = interpolate(inval[4][3],inval[6][2],lo[1],hi[1]);

    i8  = interpolate(inval[0][5],inval[4][4],lo[2],hi[2]);
    i9  = interpolate(inval[1][5],inval[5][4],lo[2],hi[2]);
    i10 = interpolate(inval[3][5],inval[7][4],lo[2],hi[2]);
    i11 = interpolate(inval[2][5],inval[6][4],lo[2],hi[2]);

    // make bits 2, 3, 6 and 7 consistent with Lewiner paper (see NOTE above)

    bit0 = v000 <= thresh ? 0 : 1;
    bit1 = v001 <= thresh ? 0 : 1;
    bit2 = v011 <= thresh ? 0 : 1;
    bit3 = v010 <= thresh ? 0 : 1;
    bit4 = v100 <= thresh ? 0 : 1;
    bit5 = v101 <= thresh ? 0 : 1;
    bit6 = v111 <= thresh ? 0 : 1;
    bit7 = v110 <= thresh ? 0 : 1;

    which = (bit7 << 7) + (bit6 << 6) + (bit5 << 5) + (bit4 << 4) +
      (bit3 << 3) + (bit2 << 2) + (bit1 << 1) + bit0;

    // icase = case of the active cube in [0..15]

    icase = cases[which][0];
    config = cases[which][1];
    subconfig = 0;

    switch (icase) {
    case  0:
      nsurf = 0;
      break;

    case  1:
      nsurf = add_triangle_inner(tiling1[config], 1);
      break;

    case  2:
      nsurf = add_triangle_inner(tiling2[config], 2);
      break;

    case  3:
      if (test_face(test3[config]))
        nsurf = add_triangle_inner(tiling3_2[config], 4); // 3.2
      else
        nsurf = add_triangle_inner(tiling3_1[config], 2); // 3.1
      break;

    case  4:
      if (modified_test_interior(test4[config],icase))
        nsurf = add_triangle_inner(tiling4_1[config], 2); // 4.1.1
      else
        nsurf = add_triangle_inner(tiling4_2[config], 6); // 4.1.2
      break;

    case  5:
      nsurf = add_triangle_inner(tiling5[config], 3);
      break;

    case  6:
      if (test_face(test6[config][0]))
        nsurf = add_triangle_inner(tiling6_2[config], 5); // 6.2
      else {
        if (modified_test_interior(test6[config][1],icase))
          nsurf = add_triangle_inner(tiling6_1_1[config], 3); // 6.1.1
        else {
          nsurf = add_triangle_inner(tiling6_1_2[config], 9); // 6.1.2
        }
      }
      break;

    case  7:
      if (test_face(test7[config][0])) subconfig +=  1;
      if (test_face(test7[config][1])) subconfig +=  2;
      if (test_face(test7[config][2])) subconfig +=  4;
      switch (subconfig) {
      case 0:
        nsurf = add_triangle_inner(tiling7_1[config], 3); break;
      case 1:
        nsurf = add_triangle_inner(tiling7_2[config][0], 5); break;
      case 2:
        nsurf = add_triangle_inner(tiling7_2[config][1], 5); break;
      case 3:
        nsurf = add_triangle_inner(tiling7_3[config][0], 9); break;
      case 4:
        nsurf = add_triangle_inner(tiling7_2[config][2], 5); break;
      case 5:
        nsurf = add_triangle_inner(tiling7_3[config][1], 9); break;
      case 6:
        nsurf = add_triangle_inner(tiling7_3[config][2], 9); break;
      case 7:
        if (test_interior(test7[config][3],icase))
          nsurf = add_triangle_inner(tiling7_4_2[config], 9);
        else
          nsurf = add_triangle_inner(tiling7_4_1[config], 5);
        break;
      };
      break;

    case  8:
      nsurf = add_triangle_inner(tiling8[config], 2);
      break;

    case  9:
      nsurf = add_triangle_inner(tiling9[config], 4);
      break;

    case 10:
      if (test_face(test10[config][0])) {
        if (test_face(test10[config][1]))
          nsurf = add_triangle_inner(tiling10_1_1_[config], 4); // 10.1.1
        else {
          nsurf = add_triangle_inner(tiling10_2[config], 8); // 10.2
        }
      } else {
        if (test_face(test10[config][1])) {
          nsurf = add_triangle_inner(tiling10_2_[config], 8); // 10.2
        } else {
          if (test_interior(test10[config][2],icase))
            nsurf = add_triangle_inner(tiling10_1_1[config], 4); // 10.1.1
          else
            nsurf = add_triangle_inner(tiling10_1_2[config], 8); // 10.1.2
        }
      }
      break;

    case 11:
      nsurf = add_triangle_inner(tiling11[config], 4);
      break;

    case 12:
      if (test_face(test12[config][0])) {
        if (test_face(test12[config][1]))
          nsurf = add_triangle_inner(tiling12_1_1_[config], 4); // 12.1.1
        else {
          nsurf = add_triangle_inner(tiling12_2[config], 8); // 12.2
        }
      } else {
        if (test_face(test12[config][1])) {
          nsurf = add_triangle_inner(tiling12_2_[config], 8); // 12.2
        } else {
          if (test_interior(test12[config][2],icase))
            nsurf = add_triangle_inner(tiling12_1_1[config], 4); // 12.1.1
          else
            nsurf = add_triangle_inner(tiling12_1_2[config], 8); // 12.1.2
        }
      }
      break;

    case 13:
      if (test_face(test13[config][0])) subconfig +=  1;
      if (test_face(test13[config][1])) subconfig +=  2;
      if (test_face(test13[config][2])) subconfig +=  4;
      if (test_face(test13[config][3])) subconfig +=  8;
      if (test_face(test13[config][4])) subconfig += 16;
      if (test_face(test13[config][5])) subconfig += 32;

      switch (subconfig13[subconfig]) {
      case 0:/* 13.1 */
        nsurf = add_triangle_inner(tiling13_1[config], 4); break;

      case 1:/* 13.2 */
        nsurf = add_triangle_inner(tiling13_2[config][0], 6); break;
      case 2:/* 13.2 */
        nsurf = add_triangle_inner(tiling13_2[config][1], 6); break;
      case 3:/* 13.2 */
        nsurf = add_triangle_inner(tiling13_2[config][2], 6); break;
      case 4:/* 13.2 */
        nsurf = add_triangle_inner(tiling13_2[config][3], 6); break;
      case 5:/* 13.2 */
        nsurf = add_triangle_inner(tiling13_2[config][4], 6); break;
      case 6:/* 13.2 */
        nsurf = add_triangle_inner(tiling13_2[config][5], 6); break;

      case 7:/* 13.3 */
        nsurf = add_triangle_inner(tiling13_3[config][0], 10); break;
      case 8:/* 13.3 */
        nsurf = add_triangle_inner(tiling13_3[config][1], 10); break;
      case 9:/* 13.3 */
        nsurf = add_triangle_inner(tiling13_3[config][2], 10); break;
      case 10:/* 13.3 */
        nsurf = add_triangle_inner(tiling13_3[config][3], 10); break;
      case 11:/* 13.3 */
        nsurf = add_triangle_inner(tiling13_3[config][4], 10); break;
      case 12:/* 13.3 */
        nsurf = add_triangle_inner(tiling13_3[config][5], 10); break;
      case 13:/* 13.3 */
        nsurf = add_triangle_inner(tiling13_3[config][6], 10); break;
      case 14:/* 13.3 */
        nsurf = add_triangle_inner(tiling13_3[config][7], 10); break;
      case 15:/* 13.3 */
        nsurf = add_triangle_inner(tiling13_3[config][8], 10); break;
      case 16:/* 13.3 */
        nsurf = add_triangle_inner(tiling13_3[config][9], 10); break;
      case 17:/* 13.3 */
        nsurf = add_triangle_inner(tiling13_3[config][10], 10); break;
      case 18:/* 13.3 */
        nsurf = add_triangle_inner(tiling13_3[config][11], 10); break;

      case 19:/* 13.4 */
        nsurf = add_triangle_inner(tiling13_4[config][0], 12); break;
      case 20:/* 13.4 */
        nsurf = add_triangle_inner(tiling13_4[config][1], 12); break;
      case 21:/* 13.4 */
        nsurf = add_triangle_inner(tiling13_4[config][2], 12); break;
      case 22:/* 13.4 */
        nsurf = add_triangle_inner(tiling13_4[config][3], 12); break;

      case 23:/* 13.5 */
        subconfig = 0;
        if (interior_test_case13())
          nsurf = add_triangle_inner(tiling13_5_1[config][0], 6);
        else
          nsurf = add_triangle_inner(tiling13_5_2[config][0], 10);
        break;

      case 24:/* 13.5 */
        subconfig = 1;
        if (interior_test_case13())
          nsurf = add_triangle_inner(tiling13_5_1[config][1], 6);
        else
          nsurf = add_triangle_inner(tiling13_5_2[config][1], 10);
        break;

      case 25:/* 13.5 */
        subconfig = 2;
        if (interior_test_case13())
          nsurf = add_triangle_inner(tiling13_5_1[config][2], 6);
        else
          nsurf = add_triangle_inner(tiling13_5_2[config][2], 10);
        break;

      case 26:/* 13.5 */
        subconfig = 3;
        if (interior_test_case13())
          nsurf = add_triangle_inner(tiling13_5_1[config][3], 6);
        else
          nsurf = add_triangle_inner(tiling13_5_2[config][3], 10);
        break;

      case 27:/* 13.3 */
        nsurf = add_triangle_inner(tiling13_3_[config][0], 10); break;
      case 28:/* 13.3 */
        nsurf = add_triangle_inner(tiling13_3_[config][1], 10); break;
      case 29:/* 13.3 */
        nsurf = add_triangle_inner(tiling13_3_[config][2], 10); break;
      case 30:/* 13.3 */
        nsurf = add_triangle_inner(tiling13_3_[config][3], 10); break;
      case 31:/* 13.3 */
        nsurf = add_triangle_inner(tiling13_3_[config][4], 10); break;
      case 32:/* 13.3 */
        nsurf = add_triangle_inner(tiling13_3_[config][5], 10); break;
      case 33:/* 13.3 */
        nsurf = add_triangle_inner(tiling13_3_[config][6], 10); break;
      case 34:/* 13.3 */
        nsurf = add_triangle_inner(tiling13_3_[config][7], 10); break;
      case 35:/* 13.3 */
        nsurf = add_triangle_inner(tiling13_3_[config][8], 10); break;
      case 36:/* 13.3 */
        nsurf = add_triangle_inner(tiling13_3_[config][9], 10); break;
      case 37:/* 13.3 */
        nsurf = add_triangle_inner(tiling13_3_[config][10], 10); break;
      case 38:/* 13.3 */
        nsurf = add_triangle_inner(tiling13_3_[config][11], 10); break;

      case 39:/* 13.2 */
        nsurf = add_triangle_inner(tiling13_2_[config][0], 6); break;
      case 40:/* 13.2 */
        nsurf = add_triangle_inner(tiling13_2_[config][1], 6); break;
      case 41:/* 13.2 */
        nsurf = add_triangle_inner(tiling13_2_[config][2], 6); break;
      case 42:/* 13.2 */
        nsurf = add_triangle_inner(tiling13_2_[config][3], 6); break;
      case 43:/* 13.2 */
        nsurf = add_triangle_inner(tiling13_2_[config][4], 6); break;
      case 44:/* 13.2 */
        nsurf = add_triangle_inner(tiling13_2_[config][5], 6); break;

      case 45:/* 13.1 */
        nsurf = add_triangle_inner(tiling13_1_[config], 4); break;

      default:
        print_cube();
        error->one(FLERR,"Marching cubes - impossible case 13");
      }
      break;

    case 14:
      nsurf = add_triangle_inner(tiling14[config], 4);
      break;
    };

    // store 4 MC labels for FixAblate caller

    mcflags[icell][0] = icase;
    mcflags[icell][1] = config;
    mcflags[icell][2] = subconfig;
    mcflags[icell][3] = nsurf;

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
      if (svalues) surf->add_tri(surfID,svalues[icell],
                                 pt[ipt+2],pt[ipt+1],pt[ipt]);
      else surf->add_tri(surfID,1,pt[ipt+2],pt[ipt+1],pt[ipt]);
      ipt += 3;
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
   adding inner triangles
------------------------------------------------------------------------- */

int MarchingCubes::add_triangle_inner(int *trig, int n)
{
  for(int t = 0; t < 3*n; t++) {
    switch (trig[t]) {
    case 0:
      pt[t][0] = i0;
      pt[t][1] = lo[1];
      pt[t][2] = lo[2];
      break;
    case 1:
      pt[t][0] = hi[0];
      pt[t][1] = i1;
      pt[t][2] = lo[2];
      break;
    case 2:
      pt[t][0] = i2;
      pt[t][1] = hi[1];
      pt[t][2] = lo[2];
      break;
    case 3:
      pt[t][0] = lo[0];
      pt[t][1] = i3;
      pt[t][2] = lo[2];
      break;
    case 4:
      pt[t][0] = i4;
      pt[t][1] = lo[1];
      pt[t][2] = hi[2];
      break;
    case 5:
      pt[t][0] = hi[0];
      pt[t][1] = i5;
      pt[t][2] = hi[2];
      break;
    case 6:
      pt[t][0] = i6;
      pt[t][1] = hi[1];
      pt[t][2] = hi[2];
      break;
    case 7:
      pt[t][0] = lo[0];
      pt[t][1] = i7;
      pt[t][2] = hi[2];
      break;
    case 8:
      pt[t][0] = lo[0];
      pt[t][1] = lo[1];
      pt[t][2] = i8;
      break;
    case 9:
      pt[t][0] = hi[0];
      pt[t][1] = lo[1];
      pt[t][2] = i9;
      break;
    case 10:
      pt[t][0] = hi[0];
      pt[t][1] = hi[1];
      pt[t][2] = i10;
      break;
    case 11:
      pt[t][0] = lo[0];
      pt[t][1] = hi[1];
      pt[t][2] = i11;
      break;
    case 12: {
      int u = 0;
      pt[t][0] = pt[t][1] = pt[t][2] = 0.0;
      if (bit0 ^ bit1) {
        ++u;
        pt[t][0] += i0;
        pt[t][1] += lo[1];
        pt[t][2] += lo[2];
      }
      if (bit1 ^ bit2) {
        ++u;
        pt[t][0] += hi[0];
        pt[t][1] += i1;
        pt[t][2] += lo[2];
      }
      if (bit2 ^ bit3) {
        ++u;
        pt[t][0] += i2;
        pt[t][1] += hi[1];
        pt[t][2] += lo[2];
      }
      if (bit3 ^ bit0) {
        ++u;
        pt[t][0] += lo[0];
        pt[t][1] += i3;
        pt[t][2] += lo[2];
      }
      if (bit4 ^ bit5) {
        ++u;
        pt[t][0] += i4;
        pt[t][1] += lo[1];
        pt[t][2] += hi[2];
      }
      if (bit5 ^ bit6) {
        ++u;
        pt[t][0] += hi[0];
        pt[t][1] += i5;
        pt[t][2] += hi[2];
      }
      if (bit6 ^ bit7) {
        ++u;
        pt[t][0] += i6;
        pt[t][1] += hi[1];
        pt[t][2] += hi[2];
      }
      if (bit7 ^ bit4) {
        ++u;
        pt[t][0] += lo[0];
        pt[t][1] += i7;
        pt[t][2] += hi[2];
      }
      if (bit0 ^ bit4) {
        ++u;
        pt[t][0] += lo[0];
        pt[t][1] += lo[1];
        pt[t][2] += i8;
      }
      if (bit1 ^ bit5) {
        ++u;
        pt[t][0] += hi[0];
        pt[t][1] += lo[1];
        pt[t][2] += i9;
      }
      if (bit2 ^ bit6) {
        ++u;
        pt[t][0] += hi[0];
        pt[t][1] += hi[1];
        pt[t][2] += i10;
      }
      if (bit3 ^ bit7) {
        ++u;
        pt[t][0] += lo[0];
        pt[t][1] += hi[1];
        pt[t][2] += i11;
      }

      pt[t][0] /= static_cast<double> (u);
      pt[t][1] /= static_cast<double> (u);
      pt[t][2] /= static_cast<double> (u);
      break;
    }

    default:
      break;
    }
  }

  return n;
}

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

#include "spatype.h"
#include "stdlib.h"
#include "string.h"
#include "fix_ablate.h"
#include "update.h"
#include "geometry.h"
#include "grid.h"
#include "domain.h"
#include "decrement_lookup_table.h"
#include "comm.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "output.h"
#include "input.h"
#include "variable.h"
#include "dump.h"
#include "marching_squares.h"
#include "marching_cubes.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{CVALUE,CDELTA,NVERT};

#define EPSILON 1.0e-4            // this is on a scale of 0 to 255

/* ----------------------------------------------------------------------
   Part 1 of 2 for multi-point decrement. Determine total amount to
   decrement in each interface corner point (an outside corner point
   is connected to an inside one and vice versa). Also determines
   number of vertices around each corner in each cell
      - cell decrement evenly divided between each interface corner
      - at this point, negative corner points are ok! Handled in part 2
      - outside corner points (all neighbors are also outside) are not
        touched
------------------------------------------------------------------------- */

void FixAblate::decrement_multid_outside()
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  int i,k;
  double total,perout,Ninterface,Nout;

  int i_cneigh;
  int *neighbors;

  int nsurf;
  double ninter;

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    for (i = 0; i < ncorner; i++) {
      cdelta[icell][i] = 0.0;
      nvert[icell][i] = 0.0;
    }

    nsurf = cells[icell].nsurf;
    if (!nsurf) continue; // if no surfs, no interface points

    // find which corners in the cell are inside, outside, and interface
    // output how many interface points there are

    if (dim == 2) mark_corners_2d(icell);
    else mark_corners_3d(icell);

    // perout is how much to decrement at each interface point

    Ninterface = find_ninter();
    total = celldelta[icell];
    perout = total/Ninterface;

    // iterate to find the number of vertices around each corner
    // also assign perout to the interface points

    for (i = 0; i < ncorner; i++) {

      // outside points have zero vertices so can skip

      if (refcorners[i] == -1) continue;

      // manually check for vertices
      // 0 -> 1,2,4
      // 1 -> 3,5
      // 2 -> 3,6
      // 3 -> 7
      // 4 -> 5,6
      // 5 -> 7
      // 6 -> 7

      if (refcorners[i] == 0) { // this corner is interface
        if (i == 0) {
          if (refcorners[1] == 1) nvert[icell][i] += 1.0;
          if (refcorners[2] == 1) nvert[icell][i] += 1.0;
        } else if (i == 1) {
          if (refcorners[3] == 1) nvert[icell][i] += 1.0;
        } else if (i==2) {
          if (refcorners[3] == 1) nvert[icell][i] += 1.0;
        }

        if (dim == 3) {
          if (i == 0) {
            if (refcorners[4] == 1) nvert[icell][i] += 1.0;
          } else if (i == 1) {
            if (refcorners[5] == 1) nvert[icell][i] += 1.0;
          } else if (i==2) {
            if (refcorners[6] == 1) nvert[icell][i] += 1.0;
          } else if (i==3) {
            if (refcorners[7] == 1) nvert[icell][i] += 1.0;
          } else if (i==4) {
            if (refcorners[5] == 1) nvert[icell][i] += 1.0;
            if (refcorners[6] == 1) nvert[icell][i] += 1.0;
          } else if (i==5) {
            if (refcorners[7] == 1) nvert[icell][i] += 1.0;
          } else if (i==6) {
            if (refcorners[7] == 1) nvert[icell][i] += 1.0;
          }
        }
      } else { // this corner is inside
        if (i == 0) {
          if (refcorners[1] == 0) nvert[icell][1] += 1.0;
          if (refcorners[2] == 0) nvert[icell][2] += 1.0;
        } else if (i == 1) {
          if (refcorners[3] == 0) nvert[icell][3] += 1.0;
        } else if (i==2) {
          if (refcorners[3] == 0) nvert[icell][3] += 1.0;
        }

        if (dim == 3) {
          if (i == 0) {
            if (refcorners[4] == 0) nvert[icell][4] += 1.0;
          } else if (i == 1) {
            if (refcorners[5] == 0) nvert[icell][5] += 1.0;
          } else if (i==2) {
            if (refcorners[6] == 0) nvert[icell][6] += 1.0;
          } else if (i==3) {
            if (refcorners[7] == 0) nvert[icell][7] += 1.0;
          } else if (i==4) {
            if (refcorners[5] == 0) nvert[icell][5] += 1.0;
            if (refcorners[6] == 0) nvert[icell][6] += 1.0;
          } else if (i==5) {
            if (refcorners[7] == 0) nvert[icell][7] += 1.0;
          } else if (i==6) {
            if (refcorners[7] == 0) nvert[icell][7] += 1.0;
          }
        }
      }

      if (refcorners[i] != 1) continue;

      // find number of interface points

      neighbors = corner_neighbor[i];

      ninter = 0;
      for (k = 0; k < dim; k++)
        if (refcorners[neighbors[k]] == 0) ninter++;

      if (ninter == 0) continue;

      Nout = 0;
      for (k = 0; k < dim; k++)
        if (refcorners[neighbors[k]] == 0) Nout++;

      if (Nout == 0) error->one(FLERR,"No outside neghbors");

      for (k = 0; k < dim; k++) {
        i_cneigh = neighbors[k];

        if (refcorners[i_cneigh] == 0)
          cdelta[icell][i_cneigh] += perout/Nout;
      }

    } // end corners
  } // end cells
}

/* ----------------------------------------------------------------------
   Sync and update interface corner values. Some values may be negative,
    this is OK. Refer to sync() for more details about logic.
------------------------------------------------------------------------- */

void FixAblate::sync_multid_outside()
{
  int i,ix,iy,iz,jx,jy,jz,ixfirst,iyfirst,izfirst,jcorner;
  int icell,jcell;
  double total;

  comm_neigh_corners(CDELTA);

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  for (icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    ix = ixyz[icell][0];
    iy = ixyz[icell][1];
    iz = ixyz[icell][2];

    for (i = 0; i < ncorner; i++) {

      ixfirst = (i % 2) - 1;
      iyfirst = (i/2 % 2) - 1;
      if (dim == 2) izfirst = 0;
      else izfirst = (i / 4) - 1;

      total = 0.0;
      jcorner = ncorner;

      for (jz = izfirst; jz <= izfirst+1; jz++) {
        for (jy = iyfirst; jy <= iyfirst+1; jy++) {
          for (jx = ixfirst; jx <= ixfirst+1; jx++) {
            jcorner--;

            if (ix+jx < 1 || ix+jx > nx) continue;
            if (iy+jy < 1 || iy+jy > ny) continue;
            if (iz+jz < 1 || iz+jz > nz) continue;

            jcell = walk_to_neigh(icell,jx,jy,jz);

            if (jcell < nglocal) total += cdelta[jcell][jcorner];
            else total += cdelta_ghost[jcell-nglocal][jcorner];

          } // end jx
        } // end jy
      } // end jz

      cvalues[icell][i] -= total;

    } // end corners
  } // end cells
}

/* ----------------------------------------------------------------------
   Part 2 of 2 for multi-point decrement. Based on the value of the
   neighboring interface values, update the inside corner point
      - inside corner points must always be > 0.0
------------------------------------------------------------------------- */

void FixAblate::decrement_multid_inside()
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  int i,j,ix,iy,iz,jx,jy,jz,ixfirst,iyfirst,izfirst,jcorner;
  int icell,jcell;
  int *neighbors,i_cneigh;
  double total_remain,nvertices;

  // find total number of vertices around each corner point
  // required to pass the correct amount from each interfact to inside point

  comm_neigh_corners(NVERT);

  for (icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    ix = ixyz[icell][0];
    iy = ixyz[icell][1];
    iz = ixyz[icell][2];

    for (i = 0; i < ncorner; i++) cdelta[icell][i] = 0.0;

    // finds which corner points are inside, interface, and outside
    // uses a lookup table. For 3D, this is decrement_lookup_table.h

    if (dim == 2) mark_corners_2d(icell);
    else mark_corners_3d(icell);

    for (i = 0; i < ncorner; i++) {

      /*---------------------------------------------------------------*/
      // sync operation to find total number of vertices around each corner

      ixfirst = (i % 2) - 1;
      iyfirst = (i/2 % 2) - 1;
      if (dim == 2) izfirst = 0;
      else izfirst = (i / 4) - 1;

      nvertices = 0.0;
      jcorner = ncorner;

      for (jz = izfirst; jz <= izfirst+1; jz++) {
        for (jy = iyfirst; jy <= iyfirst+1; jy++) {
          for (jx = ixfirst; jx <= ixfirst+1; jx++) {
            jcorner--;

            if (ix+jx < 1 || ix+jx > nx) continue;
            if (iy+jy < 1 || iy+jy > ny) continue;
            if (iz+jz < 1 || iz+jz > nz) continue;

            jcell = walk_to_neigh(icell,jx,jy,jz);

            if (jcell < nglocal) nvertices += nvert[jcell][jcorner];
            else nvertices += nvert_ghost[jcell-nglocal][jcorner];

          } // end jx
        } // end jy
      } // end jz

      /*---------------------------------------------------------------*/

      // evenly distribute the underflow from the interface points to the
      // adjacent inside corner points

      if (cvalues[icell][i] < 0) {

        // each cell edge is touched by two (four) cells in 2D (3D)
        // scale appropriately to avoid multi-counting

        if (dim == 2) nvertices *= 0.5;
        else nvertices *= 0.25;

        // determines which corners are the neighbors of "i"

        neighbors = corner_neighbor[i];

        // decrement neighboring inside points by the underflow value

        for (j = 0; j < dim; j++) {
          i_cneigh = neighbors[j];

          if(refcorners[i_cneigh] == 1) {
            total_remain = fabs(cvalues[icell][i])/ nvertices;
            cdelta[icell][i_cneigh] += total_remain;
          }
        } // end dim

        // interface points has passed all decrement so zero it out now

        cvalues[icell][i] = 0.0;

      } // end if for negative cvalues
    } // end corner
  } // end cells
}

/* ----------------------------------------------------------------------
   same as sync_multi_outside() but not the inside corner points
   are updated
------------------------------------------------------------------------- */

void FixAblate::sync_multid_inside()
{
  int i,ix,iy,iz,jx,jy,jz,ixfirst,iyfirst,izfirst,jcorner;
  int icell,jcell;
  double total;

  comm_neigh_corners(CDELTA);

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  for (icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    ix = ixyz[icell][0];
    iy = ixyz[icell][1];
    iz = ixyz[icell][2];

    for (i = 0; i < ncorner; i++) {

      ixfirst = (i % 2) - 1;
      iyfirst = (i/2 % 2) - 1;
      if (dim == 2) izfirst = 0;
      else izfirst = (i / 4) - 1;

      total = 0.0;
      jcorner = ncorner;

      for (jz = izfirst; jz <= izfirst+1; jz++) {
        for (jy = iyfirst; jy <= iyfirst+1; jy++) {
          for (jx = ixfirst; jx <= ixfirst+1; jx++) {
            jcorner--;

            if (ix+jx < 1 || ix+jx > nx) continue;
            if (iy+jy < 1 || iy+jy > ny) continue;
            if (iz+jz < 1 || iz+jz > nz) continue;

            jcell = walk_to_neigh(icell,jx,jy,jz);

            if (jcell < nglocal) {
              if (cdelta[jcell][jcorner] > 0)
                total += cdelta[jcell][jcorner];
            } else {
              if (cdelta_ghost[jcell-nglocal][jcorner] > 0)
                total += cdelta_ghost[jcell-nglocal][jcorner];
            }

          } // end jx
        } // end jy
      } // end jz

      // update the inside corner points based on the negative values
      // passed from interface points

      if (total > 0.0) {

        // the interface corner points were sync'd previously. The total
        // decrement will be counted by the number of cells which share
        // a common edge. This is 2 in 2D and 4 in 3D

        if (dim == 2) total *= 0.5;
        else total *= 0.25;

        cvalues[icell][i] -= total;
      }

      // if decrement is too large, zero out the inside corner point
      // should not underflow significantly. if it does, simulation parameters
      // need to be changed because ablation time scale is much larger than
      // the time step

      if (cvalues[icell][i] < 0) cvalues[icell][i] = 0.0;

    } // end corners
  } // end cells
}

/* ----------------------------------------------------------------------
   version of epsilon_adjust for inner indices
------------------------------------------------------------------------- */

void FixAblate::epsilon_adjust_multiv()
{
  int allin,mixflag;

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    for (int i = 0; i < ncorner; i++) {

      // check all in or out
      if (mvalues[icell][i][0] > thresh) allin = 1;
      else allin = 0;

      mixflag = 0;
      for (int j = 0; j < nmultiv; j++) {
        if (mvalues[icell][i][j] <= thresh && allin) mixflag = 1;
        if (mvalues[icell][i][j] > thresh && !allin) mixflag = 1;
      }

      // if mixflag = 1, inner indices in disagreement in terms of side
      // set to all out (inside can become out but not vice versa)

      if (mixflag) {

        for (int j = 0; j < nmultiv; j++)
          mvalues[icell][i][j] = thresh;

      // all out
      } else if (!allin) {
        for (int j = 0; j < nmultiv; j++)
          if (mvalues[icell][i][j] > thresh)
            mvalues[icell][i][j] = thresh;
      }

    } // end corner
  } // end cells
}

/* ----------------------------------------------------------------------
   ensure each corner point value is not too close to threshold
   this avoids creating tiny or zero-size surface elements
   corner_inside_min and corner_outside_max are set in store_corners()
     via epsilon method or isosurface stuffing method
------------------------------------------------------------------------- */

void FixAblate::decrement_multiv()
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  int i,j,imin;
  double minvalue,total;
  double iavg;

  // total = full amount to decrement from cell
  // cdelta[icell] = amount to decrement from each corner point of icell

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    for (i = 0; i < ncorner; i++)
      for (j = 0; j < nmultiv; j++)
        mdelta[icell][i][j] = 0.0;

    total = celldelta[icell];

    while (total > 0.0) {

      imin = -1;
      minvalue = 256.0;
      for (i = 0; i < ncorner; i++) {

        iavg = 0;
        for (j = 0; j < nmultiv; j++) iavg += mvalues[icell][i][j];
        iavg /= nmultiv;

        if (iavg > 0.0 && iavg < minvalue && mdelta[icell][i][0] == 0.0) {
          imin = i;
          minvalue = iavg;
        }

      } // end corner

      if (imin == -1) break;

      if (total < minvalue) {
        for (j = 0; j < nmultiv; j++)
          mdelta[icell][imin][j] += total;
        total = 0.0;
      } else {
        for (j = 0; j < nmultiv; j++)
          mdelta[icell][imin][j] = minvalue;
        total -= minvalue;
      }
    }

  }
}

/* ----------------------------------------------------------------------
   version of sync() for inner indices
------------------------------------------------------------------------- */

void FixAblate::sync_multiv()
{
  int i,j,ix,iy,iz,jx,jy,jz,ixfirst,iyfirst,izfirst,jcorner;
  int icell,jcell;
  double total[nmultiv];

  comm_neigh_corners(CDELTA);

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

      /*-----------------------------------------------------------*/

      // ixyz first = offset from icell of lower left cell of 2x2x2 stencil
      //              that shares the Ith corner point

      ixfirst = (i % 2) - 1;
      iyfirst = (i/2 % 2) - 1;
      if (dim == 2) izfirst = 0;
      else izfirst = (i / 4) - 1;

      // loop over 2x2x2 stencil of cells that share the corner point
      // also works for 2d, since izfirst = 0

      for (j = 0; j < nmultiv; j++) total[j] = 0.0;

      jcorner = ncorner;

      for (jz = izfirst; jz <= izfirst+1; jz++) {
        for (jy = iyfirst; jy <= iyfirst+1; jy++) {
          for (jx = ixfirst; jx <= ixfirst+1; jx++) {
            jcorner--;

            // check if neighbor cell is within bounds of ablate grid

            if (ix+jx < 1 || ix+jx > nx) continue;
            if (iy+jy < 1 || iy+jy > ny) continue;
            if (iz+jz < 1 || iz+jz > nz) continue;

            // jcell = local index of (jx,jy,jz) neighbor cell of icell

            jcell = walk_to_neigh(icell,jx,jy,jz);

            // update total with one corner point of jcell
            // jcorner descends from ncorner

            for (j = 0; j < nmultiv; j++) {
              if (jcell < nglocal) total[j] += mdelta[jcell][jcorner][j];
              else total[j] += mdelta_ghost[jcell-nglocal][jcorner][j];
            }
          }
        }
      }

      /*-----------------------------------------------------------*/

      // now decrement corners

      for (j = 0; j < nmultiv; j++) {
        mvalues[icell][i][j] -= total[j];
        if (mvalues[icell][i][j] < 0.0) mvalues[icell][i][j] = 0.0;
      }

    } // end corners
  } // end cells
}

/* ----------------------------------------------------------------------
   find how much to decrement outside corner points
------------------------------------------------------------------------- */

void FixAblate::decrement_multiv_multid_outside()
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  int i,j,k,icell;
  double total,perout,Ninterface,Nout;

  int i_in,oin,i_cneigh;
  int *ineighbors,*neighbors;

  int nsurf;
  int ninter; // inside and connected to interface point; total interface

  // find total to decrement from each corner
  // total = full amount to decrement from cell
  // cdelta[icell] = amount to decrement from each corner point of icell

  for (icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    for (i = 0; i < ncorner; i++)
      for (j = 0; j < nmultiv; j++) mdelta[icell][i][j] = 0.0;

    nsurf = cells[icell].nsurf;
    if (!nsurf) continue; // if no surfs, no interface points

    if (dim == 2) mark_corners_2d(icell);
    else mark_corners_3d(icell);

    // count number of inside points and total amount to decrement

    Ninterface = find_ninter();
    total = celldelta[icell];
    perout = total/Ninterface;

    for (i = 0; i < ncorner; i++) {

      if (refcorners[i] != 1) continue;

      neighbors = corner_neighbor[i];
      ineighbors = inner_neighbor[i];

      // find inside interface points by checking number
      // ... of adjacent outside interface points
      ninter = 0;
      for (k = 0; k < dim; k++)
        if (refcorners[neighbors[k]] == 0) ninter++;

      if (ninter == 0) continue;

      // zero out norm component if no outside interface point there

      Nout = 0;
      for (k = 0; k < dim; k++)
        if (refcorners[neighbors[k]] == 0) Nout++;

      // scale values so their sum is one (L1 norm over L2 norm)

      if (Nout == 0) error->one(FLERR,"No outside neghbors");

      // update inner indices of interface point connected to inside point

      for (k = 0; k < dim; k++) {
        i_in = ineighbors[k];
        i_cneigh = neighbors[k];

        if (i_in == 0) oin = 1;
        else if (i_in == 1) oin = 0;
        else if (i_in == 2) oin = 3;
        else if (i_in == 3) oin = 2;
        else if (i_in == 4) oin = 5;
        else if (i_in == 5) oin = 4;

        if (refcorners[i_cneigh] == 0)
          mdelta[icell][i_cneigh][oin] += perout/Nout;
      } // end dim

    } // end corner
  } // end cells
}

/* ----------------------------------------------------------------------
   version of sync_multi_outside() for inner indices
------------------------------------------------------------------------- */

void FixAblate::sync_multiv_multid_outside()
{
  int i,j,ix,iy,iz,jx,jy,jz,ixfirst,iyfirst,izfirst,jcorner;
  int icell,jcell;
  double total[6];

  comm_neigh_corners(CDELTA);

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  for (icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    ix = ixyz[icell][0];
    iy = ixyz[icell][1];
    iz = ixyz[icell][2];

    for (i = 0; i < ncorner; i++) {

      ixfirst = (i % 2) - 1;
      iyfirst = (i/2 % 2) - 1;
      if (dim == 2) izfirst = 0;
      else izfirst = (i / 4) - 1;

      for (j = 0; j < 6; j++) total[j] = 0.0;

      jcorner = ncorner;

      for (jz = izfirst; jz <= izfirst+1; jz++) {
        for (jy = iyfirst; jy <= iyfirst+1; jy++) {
          for (jx = ixfirst; jx <= ixfirst+1; jx++) {
            jcorner--;

            if (ix+jx < 1 || ix+jx > nx) continue;
            if (iy+jy < 1 || iy+jy > ny) continue;
            if (iz+jz < 1 || iz+jz > nz) continue;

            jcell = walk_to_neigh(icell,jx,jy,jz);

            for (j = 0; j < nmultiv; j++) {
              if (jcell < nglocal) total[j] += mdelta[jcell][jcorner][j];
              else total[j] += mdelta_ghost[jcell-nglocal][jcorner][j];
            }

          } // end jx
        } // end jy
      } // end jz

      for (j = 0; j < nmultiv; j++)
        mvalues[icell][i][j] -= total[j];

    } // end corners
  } // end cells
}

/* ----------------------------------------------------------------------
   version of decrement_multi_inside() for inner indices
------------------------------------------------------------------------- */

void FixAblate::decrement_multiv_multid_inside()
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  int i,j,icell;
  int i_in,o_in,i_cneigh;
  int *ineighbors,*neighbors;
  double total_remain;

  for (icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    for (i = 0; i < ncorner; i++)
      for (j = 0; j < nmultiv; j++) mdelta[icell][i][j] = 0.0;

    if (dim == 2) mark_corners_2d(icell);
    else mark_corners_3d(icell);

    for (i = 0; i < ncorner; i++) {

      // only check inside corners

      if (mvalues[icell][i][0] < 0) continue;

      // corner and inner neighbors

      neighbors = corner_neighbor[i];
      ineighbors = inner_neighbor[i];

      for (j = 0; j < dim; j++) {
        i_in = ineighbors[j];
        i_cneigh = neighbors[j];

        if (i_in == 0) o_in = 1;
        else if (i_in == 1) o_in = 0;
        else if (i_in == 2) o_in = 3;
        else if (i_in == 3) o_in = 2;
        else if (i_in == 4) o_in = 5;
        else if (i_in == 5) o_in = 4;
        else error->one(FLERR,"Bad inner index");

        total_remain = mvalues[icell][i_cneigh][o_in];
        if (total_remain < 0)
          mdelta[icell][i][i_in] += fabs(total_remain);
      } // end dim

    } // end corner

    // zero out negative values

    for (i = 0; i < ncorner; i++)
      for (j = 0; j < nmultiv; j++)
        if (mvalues[icell][i][j] < 0.0) mvalues[icell][i][j] = 0.0;

  } // end cells
}

/* ----------------------------------------------------------------------
   version of sync_multi_inside() for inner indices
------------------------------------------------------------------------- */

void FixAblate::sync_multiv_multid_inside()
{
  int i,j,ix,iy,iz,jx,jy,jz,ixfirst,iyfirst,izfirst,jcorner;
  int icell,jcell;
  double total[nmultiv];

  comm_neigh_corners(CDELTA);

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  for (icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    ix = ixyz[icell][0];
    iy = ixyz[icell][1];
    iz = ixyz[icell][2];

    for (i = 0; i < ncorner; i++) {

      ixfirst = (i % 2) - 1;
      iyfirst = (i/2 % 2) - 1;
      if (dim == 2) izfirst = 0;
      else izfirst = (i / 4) - 1;

      for (j = 0; j < nmultiv; j++) total[j] = 0.0;

      jcorner = ncorner;

      for (jz = izfirst; jz <= izfirst+1; jz++) {
        for (jy = iyfirst; jy <= iyfirst+1; jy++) {
          for (jx = ixfirst; jx <= ixfirst+1; jx++) {
            jcorner--;

            if (ix+jx < 1 || ix+jx > nx) continue;
            if (iy+jy < 1 || iy+jy > ny) continue;
            if (iz+jz < 1 || iz+jz > nz) continue;

            jcell = walk_to_neigh(icell,jx,jy,jz);

            if (jcell < nglocal) {

              for (j = 0; j < nmultiv; j++)
                if (mdelta[jcell][jcorner][j] > 0)
                  total[j] += mdelta[jcell][jcorner][j];

            } else {

              for (j = 0; j < nmultiv; j++)
                if (mdelta_ghost[jcell-nglocal][jcorner][j] > 0)
                  total[j] += mdelta_ghost[jcell-nglocal][jcorner][j];

            }

          } // jz
        } // jy
      } // jx

      for (j = 0; j < nmultiv; j++) {
        if (total[j] > 0.0) {
          if (dim == 2) total[j] *= 0.5;
          else total[j] *= 0.25;
          mvalues[icell][i][j] -= total[j];
        }
      }

      for (j = 0; j < nmultiv; j++)
        if (mvalues[icell][i][j] < 0.0) mvalues[icell][i][j] = 0.0;

    } // end corners
  } // end cells
}

/* ----------------------------------------------------------------------
   determines how many interface points and marks each corner point as
   either inside, outside, or interface.

   1 - inside
   -1 - outside
   0 - interface

   All neighbors of inside points have values >= thresh. Similarly, all
   neighbors of outside poitns have values < thresh. An interface point
   is has a value < thresh but at least one of its neighbors has a value
   >= thresh. In other words, there is a vertex located between an interface
   point and an inside point

   refcorners stores corner point type and follows SPARTA corner ordering
------------------------------------------------------------------------- */
void FixAblate::mark_corners_2d(int icell)
{
  for (int i = 0; i < ncorner; i++)
    refcorners[i] = 0;

  // mark inside corners first

  int Nin = 0;

  if (multi_val_flag) {
    if (mvalues[icell][0][0] > thresh) {
      refcorners[0] = 1;
      Nin++;
    }

    if (mvalues[icell][1][0] > thresh) {
      refcorners[1] = 1;
      Nin++;
    }

    if (mvalues[icell][2][0] > thresh) {
      refcorners[2] = 1;
      Nin++;
    }

    if (mvalues[icell][3][0] > thresh) {
      refcorners[3] = 1;
      Nin++;
    }

  } else {
    if (cvalues[icell][0] > thresh) {
      refcorners[0] = 1;
      Nin++;
    }

    if (cvalues[icell][1] > thresh) {
      refcorners[1] = 1;
      Nin++;
    }

    if (cvalues[icell][2] > thresh) {
      refcorners[2] = 1;
      Nin++;
    }

    if (cvalues[icell][3] > thresh) {
      refcorners[3] = 1;
      Nin++;
    }
  }

  // built-in lookup table for marking outside and interface points

  if (Nin == 0) {
    refcorners[0] = refcorners[1] = refcorners[2] = refcorners[3] = -1;
  } else if (Nin == 1) {
    if (refcorners[0] == 1) refcorners[3] = -1;
    else if (refcorners[1] == 1) refcorners[2] = -1;
    else if (refcorners[2] == 1) refcorners[1] = -1;
    else refcorners[0] = -1;
  }

}

/* ----------------------------------------------------------------------
   same as mark_corners_2d() for 3d
   - bit stores whether corner value is below or above thresh

   IMPORTANT: decrement_lookup_table follows "standard" ordering reported
   in the literature and is different from SPARTA. Refer to diagrams in
   decrement_lookup_table
------------------------------------------------------------------------- */

void FixAblate::mark_corners_3d(int icell)
{

  // find which corners are inside or outside the surface

  int i,bit[8],which;
  if (multi_val_flag) {
    for (i = 0; i < ncorner; i++)
      bit[i] = mvalues[icell][i][0] <= thresh ? 0 : 1;
  } else {
    for (i = 0; i < ncorner; i++)
      bit[i] = cvalues[icell][i] <= thresh ? 0 : 1;
  }

  // find configuration case (whic)
  // note: that because there are two different orderings, some indices
  // are switched

  which = (bit[6] << 7) + (bit[7] << 6) + (bit[5] << 5) + (bit[4] << 4) +
    (bit[2] << 3) + (bit[3] << 2) + (bit[1] << 1) + bit[0];

  // find how many inside and interface points based on configuration
  // note: Nin + Ninterface does not always equal ncorner

  int Nin = Ninside[which]; // inside (includes inside interface)
  int Nointerface = Noutside[which]; // outside interface

  int *incorners, *outcorners;
  incorners = inside[which]; // returns list of inside corner points
  outcorners = outside[which]; // returns list of oustide cornr points

  // mark the outside, then inside, then interface points

  for (i = 0; i < ncorner; i++)
    refcorners[i] = -1;

  for (i = 0; i < Nin; i++) {
    if (incorners[i] == 2) refcorners[3] = 1;
    else if (incorners[i] == 3) refcorners[2] = 1;
    else if (incorners[i] == 6) refcorners[7] = 1;
    else if (incorners[i] == 7) refcorners[6] = 1;
    else refcorners[incorners[i]] = 1;
  }

  for (i = 0; i < Nointerface; i++) {
    if (outcorners[i] == 2) refcorners[3] = 0;
    else if (outcorners[i] == 3) refcorners[2] = 0;
    else if (outcorners[i] == 6) refcorners[7] = 0;
    else if (outcorners[i] == 7) refcorners[6] = 0;
    else refcorners[outcorners[i]] = 0;
  }
}

/* ----------------------------------------------------------------------
   finds number of inside interface points
------------------------------------------------------------------------- */

int FixAblate::find_ninter()
{
  int total = 0;
  int add;
  int *neighbors, i_cneigh;
  for (int i = 0; i < ncorner; i++) {
    if (refcorners[i] == 1) {
      neighbors = corner_neighbor[i];
      add = 0;
      for (int j = 0; j < dim; j++) {
        i_cneigh = neighbors[j];
        if (refcorners[i_cneigh] == 0) add = 1;
      }
      if (add) total++;
    }
  }
  return total;
}

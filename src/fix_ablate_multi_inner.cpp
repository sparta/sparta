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

  int i,j,k,m,Ninterface;
  double total,perout;

  int i_inner,oinner,i_neigbor_corner,i_neighbor_corner;
  int *ineighbors,*neighbors;

  Surf::Line *line;
  Surf::Line *lines = surf->lines;
  Surf::Tri *tri;
  Surf::Tri *tris = surf->tris;

  surfint *csurfs;
  double norm[3],pt[3],totalnorm;
  double dist,smindist;
  int nsurf,isurf;
  int nin; // inside and connected to interface point
  double totalperout, scale, ninter;

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    for (i = 0; i < ncorner; i++) {
      cdelta[icell][i] = 0.0;
      nvert[icell][i] = 0.0;
    }

    nsurf = cells[icell].nsurf;
    if (!nsurf) continue; // if no surfs, no interface points
    csurfs = cells[icell].csurfs;

    // find which corners in the cell are inside, outside, and interface
    // output how many interface points there are

    if (dim == 2) mark_corners_2d(icell);
    else mark_corners_3d(icell);

    // if no interface points, then points are either all outside or all inside
    // no mechanism for ablation if there is no surfaces so can skip

    //if (Ninterface == 0) continue;

    // perout is how much to decrement at each interface point

    nin = find_nin();
    total = celldelta[icell];
    perout = total/nin;//Ninterface;

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

      // find location of corner point

      if (i == 0) {
        pt[0] = cells[icell].lo[0];
        pt[1] = cells[icell].lo[1];
        pt[2] = cells[icell].lo[2];
      } else if (i == 1) {
        pt[0] = cells[icell].hi[0];
        pt[1] = cells[icell].lo[1];
        pt[2] = cells[icell].lo[2];
      } else if (i == 2) {
        pt[0] = cells[icell].lo[0];
        pt[1] = cells[icell].hi[1];
        pt[2] = cells[icell].lo[2];
      } else if (i == 3) {
        pt[0] = cells[icell].hi[0];
        pt[1] = cells[icell].hi[1];
        pt[2] = cells[icell].lo[2];
      } else if (i == 4) {
        pt[0] = cells[icell].lo[0];
        pt[1] = cells[icell].lo[1];
        pt[2] = cells[icell].hi[2];
      } else if (i == 5) {
        pt[0] = cells[icell].hi[0];
        pt[1] = cells[icell].lo[1];
        pt[2] = cells[icell].hi[2];
      } else if (i == 6) {
        pt[0] = cells[icell].lo[0];
        pt[1] = cells[icell].hi[1];
        pt[2] = cells[icell].hi[2];
      } else {
        pt[0] = cells[icell].hi[0];
        pt[1] = cells[icell].hi[1];
        pt[2] = cells[icell].hi[2];
      }

      // find closest surface and record the norm

      smindist = -1;
      norm[0] = norm[1] = norm[2] = 0.0;
      for (m = 0; m < nsurf; m++) {
        isurf = csurfs[m];
        
        if (dim == 2) {
          line = &lines[isurf];
          dist = Geometry::distsq_point_line(pt, line->p1, line->p2);
          dist = fabs(dist);
          if (smindist < 0 || dist < smindist) {
            smindist = dist;
            norm[0] = line->norm[0];
            norm[1] = line->norm[1];
            norm[2] = line->norm[2];
          }
        } else {
          tri = &tris[isurf];
          dist = Geometry::distsq_point_tri(pt, tri->p1, tri->p2, tri->p3, tri->norm);
          dist = fabs(dist);
          if (smindist < 0 || dist < smindist) {
            smindist = dist;
            norm[0] = tri->norm[0];
            norm[1] = tri->norm[1];
            norm[2] = tri->norm[2];
          }
        }
      } // end surfs

      // scale values of norm so their sum is one (L1 norm over L2 norm)

      if (dim == 2) {
        totalnorm = fabs(norm[0]) + fabs(norm[1]);
        norm[0] /= totalnorm;
        norm[1] /= totalnorm;
        norm[2] = 0.0;
      } else {
        totalnorm = fabs(norm[0]) + fabs(norm[1]) + fabs(norm[2]);
        norm[0] /= totalnorm;
        norm[1] /= totalnorm;
        norm[2] /= totalnorm;
      }

      // update interface point connected to inside point

      neighbors = corner_neighbor[i];
      ineighbors = inner_neighbor[i];

      totalperout = 0.0;
      ninter = 0;
      for (k = 0; k < dim; k++) {
        i_inner = ineighbors[k];
        i_neighbor_corner = neighbors[k];

        if (refcorners[i_neighbor_corner] == 0) {
          if (i_inner < 2) {
            cdelta[icell][i_neighbor_corner] += perout*fabs(norm[0]);
            totalperout += perout*fabs(norm[0]);
            ninter++;
          } else if (i_inner < 4) {
            cdelta[icell][i_neighbor_corner] += perout*fabs(norm[1]);
            totalperout += perout*fabs(norm[1]);
            ninter++;
          } else {
            cdelta[icell][i_neighbor_corner] += perout*fabs(norm[2]);
            totalperout += perout*fabs(norm[2]);
            ninter++;
          }
        }
      } // end dim

      // scale decrement

      totalperout = (perout-totalperout)/ninter;
      for (k = 0; k < dim; k++) {
        i_neighbor_corner = neighbors[k];
        if (refcorners[i_neighbor_corner] == 0) {
          cdelta[icell][i_neighbor_corner] += totalperout;
          cdelta[icell][i_neighbor_corner] =
            MAX(cdelta[icell][i_neighbor_corner],0.0);
        }
      } // end dim
      

      //if (refcorners[i] != 0) continue;

      //cdelta[icell][i] += perout;

    } // end corners
  } // end cells
}

/* ----------------------------------------------------------------------
   sync and update interface corner values. some values may be negative.
   this is ok. refer to sync() for more details about logic
------------------------------------------------------------------------- */

void FixAblate::sync_multid_outside()
{
  int i,j,ix,iy,iz,jx,jy,jz,ixfirst,iyfirst,izfirst,jcorner;
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
  int Nneigh,*neighbors,i_neighbor_corner;
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
          i_neighbor_corner = neighbors[j];

          if(refcorners[i_neighbor_corner] == 1) {
            total_remain = fabs(cvalues[icell][i])/ nvertices;
            cdelta[icell][i_neighbor_corner] += total_remain;
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
  int i,j,ix,iy,iz,jx,jy,jz,ixfirst,iyfirst,izfirst,jcorner;
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
   store grid corner point and type values in cvalues and tvalues
   then create implicit surfaces
   called by ReadIsurf when corner point grid is read in
------------------------------------------------------------------------- */

void FixAblate::store_corners(int nx_caller, int ny_caller, int nz_caller,
                              double *cornerlo_caller, double *xyzsize_caller,
                              double ***ivalues_caller, int *tvalues_caller,
                              double thresh_caller, char *sgroupID, int pushflag)
{
  storeflag = 1;
  multi_val_flag = 1;

  nx = nx_caller;
  ny = ny_caller;
  nz = nz_caller;
  cornerlo[0] = cornerlo_caller[0];
  cornerlo[1] = cornerlo_caller[1];
  cornerlo[2] = cornerlo_caller[2];
  xyzsize[0] = xyzsize_caller[0];
  xyzsize[1] = xyzsize_caller[1];
  xyzsize[2] = xyzsize_caller[2];
  thresh = thresh_caller;

  tvalues_flag = 0;
  if (tvalues_caller) tvalues_flag = 1;

  if (sgroupID) {
    int sgroup = surf->find_group(sgroupID);
    if (sgroup < 0) sgroup = surf->add_group(sgroupID);
    sgroupbit = surf->bitmask[sgroup];
  } else sgroupbit = 0;

  // allocate per-grid cell data storage

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  nglocal = grid->nlocal;

  grow_percell(0);

  // copy caller values into local values of FixAblate

  for (int icell = 0; icell < nglocal; icell++) {
    for (int m = 0; m < ncorner; m++)
      for (int n = 0; n < nmultiv; n++)
        ivalues[icell][m][n] = ivalues_caller[icell][m][n];
    if (tvalues_flag) tvalues[icell] = tvalues_caller[icell];
  }

  // set all values to either min or max value

  if (minmaxflag) {
    for (int icell = 0; icell < nglocal; icell++) {
      for (int m = 0; m < ncorner; m++) {
        for (int n = 0; n < nmultiv; n++) {
          if (ivalues[icell][m][n] < thresh) ivalues[icell][m][n] = 0.0;
          else ivalues[icell][m][n] = 255.0;
        }
      }
    }
  }

  // set ix,iy,iz indices from 1 to Nxyz for each of my owned grid cells
  // same logic as ReadIsurf::create_hash()

  for (int i = 0; i < nglocal; i++)
    ixyz[i][0] = ixyz[i][1] = ixyz[i][2] = 0;

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    ixyz[icell][0] =
      static_cast<int> ((cells[icell].lo[0]-cornerlo[0]) / xyzsize[0] + 0.5) + 1;
    ixyz[icell][1] =
      static_cast<int> ((cells[icell].lo[1]-cornerlo[1]) / xyzsize[1] + 0.5) + 1;
    ixyz[icell][2] =
      static_cast<int> ((cells[icell].lo[2]-cornerlo[2]) / xyzsize[2] + 0.5) + 1;
  }

  if (mindist == 0.0) {
    corner_inside_min = thresh + EPSILON;
    corner_outside_max = thresh - EPSILON;
  } else {
    corner_inside_min = (thresh - 0.0*mindist) / (1.0-mindist);
    corner_outside_max = (thresh - 255.0*mindist) / (1.0-mindist);
  }

  corner_inside_min = MIN(corner_inside_min,255.0);
  corner_outside_max = MAX(corner_outside_max,0.0);

  epsilon_adjust_multiv();

  // create marching squares/cubes classes, now that have group & threshold

  if (dim == 2) ms = new MarchingSquares(sparta,igroup,thresh);
  else mc = new MarchingCubes(sparta,igroup,thresh);

  // create implicit surfaces

  create_surfs(1);
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
      if (ivalues[icell][i][0] > corner_inside_min) allin = 1;
      else allin = 0;

      mixflag = 0;
      for (int j = 0; j < nmultiv; j++) {
        if (ivalues[icell][i][j] < corner_inside_min && allin) mixflag = 1;
        if (ivalues[icell][i][j] > corner_inside_min && !allin) mixflag = 1;
      }

      // if mixflag = 1, inner indices in disagreement in terms of side
      // set to all out (inside can become out but not vice versa)

      if (mixflag) {

        for (int j = 0; j < nmultiv; j++)
          ivalues[icell][i][j] = corner_outside_max;

      // all out
      } else if (!allin) {
        for (int j = 0; j < nmultiv; j++)
          if (ivalues[icell][i][j] > corner_outside_max)
            ivalues[icell][i][j] = corner_outside_max;
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

  int i,j,i_inner,imin;
  double minvalue,total;
  double iavg;

  // total = full amount to decrement from cell
  // cdelta[icell] = amount to decrement from each corner point of icell

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    for (i = 0; i < ncorner; i++)
      for (j = 0; j < nmultiv; j++)
        idelta[icell][i][j] = 0.0;

    total = celldelta[icell];

    while (total > 0.0) {

      imin = -1;
      minvalue = 256.0;
      for (i = 0; i < ncorner; i++) {

        iavg = 0;
        for (j = 0; j < nmultiv; j++) iavg += ivalues[icell][i][j];
        iavg /= nmultiv;

        if (iavg > 0.0 && iavg < minvalue && idelta[icell][i][0] == 0.0) {
          imin = i;
          minvalue = iavg;
        }

      } // end corner

      if (imin == -1) break;

      if (total < minvalue) {
        for (j = 0; j < nmultiv; j++)
          idelta[icell][imin][j] += total;
        total = 0.0;
      } else {
        for (j = 0; j < nmultiv; j++)
          idelta[icell][imin][j] = minvalue;
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
  double total[nmultiv],inner_total;

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
              if (jcell < nglocal) total[j] += idelta[jcell][jcorner][j];
              else total[j] += idelta_ghost[jcell-nglocal][jcorner][j];
            }
          }
        }
      }

      /*-----------------------------------------------------------*/

      // now decrement corners

      for (j = 0; j < nmultiv; j++) {
        ivalues[icell][i][j] -= total[j];
        if (ivalues[icell][i][j] < 0.0) ivalues[icell][i][j] = 0.0;
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

  int i,j,k,m,icell;
  int i_inner,oinner,i_neigbor_corner,i_neighbor_corner;
  int *ineighbors,*neighbors;
  double total,perout;

  Surf::Line *line;
  Surf::Line *lines = surf->lines;
  Surf::Tri *tri;
  Surf::Tri *tris = surf->tris;

  surfint *csurfs;
  double norm[3],pt[3],totalnorm;
  double dist,smindist;
  int nsurf,isurf;
  int nin, ninter; // inside and connected to interface point; total interface
  double totalperout;

  // find total to decrement from each corner
  // total = full amount to decrement from cell
  // cdelta[icell] = amount to decrement from each corner point of icell

  for (icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    for (i = 0; i < ncorner; i++)
      for (j = 0; j < nmultiv; j++) idelta[icell][i][j] = 0.0;

    nsurf = cells[icell].nsurf;
    if (!nsurf) continue; // if no surfs, no interface points
    csurfs = cells[icell].csurfs;

    if (dim == 2) mark_corners_2d(icell);
    else mark_corners_3d(icell);

    //if (Ninterface == 0) continue;

    // count number of inside points and total amount to decrement

    nin = find_nin();
    total = celldelta[icell];
    perout = total/nin;

    for (i = 0; i < ncorner; i++) {

      if (refcorners[i] != 1) continue;

      // find location of corner point

      if (i == 0) {
        pt[0] = cells[icell].lo[0];
        pt[1] = cells[icell].lo[1];
        pt[2] = cells[icell].lo[2];
      } else if (i == 1) {
        pt[0] = cells[icell].hi[0];
        pt[1] = cells[icell].lo[1];
        pt[2] = cells[icell].lo[2];
      } else if (i == 2) {
        pt[0] = cells[icell].lo[0];
        pt[1] = cells[icell].hi[1];
        pt[2] = cells[icell].lo[2];
      } else if (i == 3) {
        pt[0] = cells[icell].hi[0];
        pt[1] = cells[icell].hi[1];
        pt[2] = cells[icell].lo[2];
      } else if (i == 4) {
        pt[0] = cells[icell].lo[0];
        pt[1] = cells[icell].lo[1];
        pt[2] = cells[icell].hi[2];
      } else if (i == 5) {
        pt[0] = cells[icell].hi[0];
        pt[1] = cells[icell].lo[1];
        pt[2] = cells[icell].hi[2];
      } else if (i == 6) {
        pt[0] = cells[icell].lo[0];
        pt[1] = cells[icell].hi[1];
        pt[2] = cells[icell].hi[2];
      } else {
        pt[0] = cells[icell].hi[0];
        pt[1] = cells[icell].hi[1];
        pt[2] = cells[icell].hi[2];
      }

      // find closest surface and record the norm

      smindist = -1;
      norm[0] = norm[1] = norm[2] = 0.0;
      for (m = 0; m < nsurf; m++) {
        isurf = csurfs[m];
        
        if (dim == 2) {
          line = &lines[isurf];
          dist = Geometry::distsq_point_line(pt, line->p1, line->p2);
          dist = fabs(dist);
          if (smindist < 0 || dist < smindist) {
            smindist = dist;
            norm[0] = line->norm[0];
            norm[1] = line->norm[1];
            norm[2] = line->norm[2];
          }
        } else {
          tri = &tris[isurf];
          dist = Geometry::distsq_point_tri(pt, tri->p1, tri->p2, tri->p3, tri->norm);
          dist = fabs(dist);
          if (smindist < 0 || dist < smindist) {
            smindist = dist;
            norm[0] = tri->norm[0];
            norm[1] = tri->norm[1];
            norm[2] = tri->norm[2];
          }
        }
      } // end surfs
      
      // scale values of norm so their sum is one (L1 norm over L2 norm)

      if (dim == 2) {
        totalnorm = fabs(norm[0]) + fabs(norm[1]);
        norm[0] /= totalnorm;
        norm[1] /= totalnorm;
        norm[2] = 0.0;
      } else {
        totalnorm = fabs(norm[0]) + fabs(norm[1]) + fabs(norm[2]);
        norm[0] /= totalnorm;
        norm[1] /= totalnorm;
        norm[2] /= totalnorm;
      }

      // update inner indices of interface point connected to inside point

      neighbors = corner_neighbor[i];
      ineighbors = inner_neighbor[i];

      totalperout = 0.0;
      ninter = 0;
      for (k = 0; k < dim; k++) {
        i_inner = ineighbors[k];
        i_neighbor_corner = neighbors[k];

        if (i_inner == 0) oinner = 1;
        else if (i_inner == 1) oinner = 0;
        else if (i_inner == 2) oinner = 3;
        else if (i_inner == 3) oinner = 2;
        else if (i_inner == 4) oinner = 5;
        else if (i_inner == 5) oinner = 4;

        if (refcorners[i_neighbor_corner] == 0) {
          if (i_inner < 2) {
            idelta[icell][i_neighbor_corner][oinner] += perout*fabs(norm[0]);
            totalperout += perout*fabs(norm[0]);
            ninter++;
          } else if (i_inner < 4) {
            idelta[icell][i_neighbor_corner][oinner] += perout*fabs(norm[1]);
            totalperout += perout*fabs(norm[1]);
            ninter++;
          } else {
            idelta[icell][i_neighbor_corner][oinner] += perout*fabs(norm[2]);
            totalperout += perout*fabs(norm[2]);
            ninter++;
          }
        }
      } // end dim

      // scale decrement

      totalperout = (perout-totalperout)/ninter;
      for (k = 0; k < dim; k++) {
        i_inner = ineighbors[k];
        i_neighbor_corner = neighbors[k];

        if (i_inner == 0) oinner = 1;
        else if (i_inner == 1) oinner = 0;
        else if (i_inner == 2) oinner = 3;
        else if (i_inner == 3) oinner = 2;
        else if (i_inner == 4) oinner = 5;
        else if (i_inner == 5) oinner = 4;

        if (refcorners[i_neighbor_corner] == 0) {
          idelta[icell][i_neighbor_corner][oinner] += totalperout;
          idelta[icell][i_neighbor_corner][oinner] =
            MAX(idelta[icell][i_neighbor_corner][oinner],0.0);
        }
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
              if (jcell < nglocal) total[j] += idelta[jcell][jcorner][j];
              else total[j] += idelta_ghost[jcell-nglocal][jcorner][j];
            }

          } // end jx
        } // end jy
      } // end jz

      for (j = 0; j < nmultiv; j++)
        ivalues[icell][i][j] -= total[j];

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
  int i_inner,oinner,i_neighbor_corner;
  int *ineighbors,*neighbors;
  double total_remain;

  for (icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    for (i = 0; i < ncorner; i++)
      for (j = 0; j < nmultiv; j++) idelta[icell][i][j] = 0.0;

    if (dim == 2) mark_corners_2d(icell);
    else mark_corners_3d(icell);

    for (i = 0; i < ncorner; i++) {

      // only check inside corners

      if (ivalues[icell][i][0] < 0) continue;

      // corner and inner neighbors

      neighbors = corner_neighbor[i];
      ineighbors = inner_neighbor[i];

      for (j = 0; j < dim; j++) {
        i_inner = ineighbors[j];
        i_neighbor_corner = neighbors[j];

        if (i_inner == 0) oinner = 1;
        else if (i_inner == 1) oinner = 0;
        else if (i_inner == 2) oinner = 3;
        else if (i_inner == 3) oinner = 2;
        else if (i_inner == 4) oinner = 5;
        else if (i_inner == 5) oinner = 4;
        else error->one(FLERR,"Bad inner index");

        total_remain = ivalues[icell][i_neighbor_corner][oinner];
        if (total_remain < 0)
          idelta[icell][i][i_inner] += fabs(total_remain);
      } // end dim

    } // end corner

    // zero out negative values

    for (i = 0; i < ncorner; i++)
      for (j = 0; j < nmultiv; j++)
        if (ivalues[icell][i][j] < 0.0) ivalues[icell][i][j] = 0.0;

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
                if (idelta[jcell][jcorner][j] > 0)
                  total[j] += idelta[jcell][jcorner][j];

            } else {

              for (j = 0; j < nmultiv; j++)
                if (idelta_ghost[jcell-nglocal][jcorner][j] > 0)
                  total[j] += idelta_ghost[jcell-nglocal][jcorner][j];

            }

          } // jz
        } // jy
      } // jx

      for (j = 0; j < nmultiv; j++) {
        if (total[j] > 0.0) {
          if (dim == 2) total[j] *= 0.5;
          else total[j] *= 0.25;
          ivalues[icell][i][j] -= total[j];
        }
      }

      for (j = 0; j < nmultiv; j++)
        if (ivalues[icell][i][j] < 0.0) ivalues[icell][i][j] = 0.0;

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
int FixAblate::mark_corners_2d(int icell)
{
  for (int i = 0; i < ncorner; i++)
    refcorners[i] = 0;

  // mark inside corners first

  int Nin = 0;

  if (multi_val_flag) {
    if (ivalues[icell][0][0] > thresh) {
      refcorners[0] = 1;
      Nin++;
    }

    if (ivalues[icell][1][0] > thresh) {
      refcorners[1] = 1;
      Nin++;
    }

    if (ivalues[icell][2][0] > thresh) {
      refcorners[2] = 1;
      Nin++;
    }

    if (ivalues[icell][3][0] > thresh) {
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

  int Ninterface;
  if (Nin == 0) {
    Ninterface = 0;
    refcorners[0] = refcorners[1] = refcorners[2] = refcorners[3] = -1;
  } else if (Nin == 1) {
    Ninterface = 2;
    if (refcorners[0] == 1) refcorners[3] = -1;
    else if (refcorners[1] == 1) refcorners[2] = -1;
    else if (refcorners[2] == 1) refcorners[1] = -1;
    else refcorners[0] = -1;
  }
  else if (Nin == 2) Ninterface = 2;
  else if (Nin == 3) Ninterface = 1;
  else Ninterface = 0;

  return Ninterface;
}


/* ----------------------------------------------------------------------
   same as mark_corners_2d() for 3d
   - bit stores whether corner value is below or above thresh

   IMPORTANT: decrement_lookup_table follows "standard" ordering reported
   in the literature and is different from SPARTA. Refer to diagrams in
   decrement_lookup_table
------------------------------------------------------------------------- */

int FixAblate::mark_corners_3d(int icell)
{

  // find which corners are inside or outside the surface

  int i,bit[8],which;
  if (multi_val_flag) {
    for (i = 0; i < ncorner; i++)
      bit[i] = ivalues[icell][i][0] <= thresh ? 0 : 1;
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

  int Nin = Ninside[which];
  int Ninterface = Noutside[which];

  int *incorners, *outcorners;
  incorners = inside[which]; // returns list of inside corner points
  outcorners = outside[which]; // returns list of oustide cornr points

  // mark the outside, then inside, then interface points

  int icorner;
  for (i = 0; i < ncorner; i++)
    refcorners[i] = -1;

  for (i = 0; i < Nin; i++) {
    if (incorners[i] == 2) refcorners[3] = 1;
    else if (incorners[i] == 3) refcorners[2] = 1;
    else if (incorners[i] == 6) refcorners[7] = 1;
    else if (incorners[i] == 7) refcorners[6] = 1;
    else refcorners[incorners[i]] = 1;
  }

  for (i = 0; i < Ninterface; i++) {
    if (outcorners[i] == 2) refcorners[3] = 0;
    else if (outcorners[i] == 3) refcorners[2] = 0;
    else if (outcorners[i] == 6) refcorners[7] = 0;
    else if (outcorners[i] == 7) refcorners[6] = 0;
    else refcorners[outcorners[i]] = 0;
  }

  return Ninterface;
}

/* ----------------------------------------------------------------------
   finds number of inside points touch an interface point
------------------------------------------------------------------------- */

int FixAblate::find_nin()
{
  int total = 0;
  int add; 
  int *neighbors, i_neighbor_corner;
  for (int i = 0; i < ncorner; i++) {
    if (refcorners[i] == 1) {
      neighbors = corner_neighbor[i];
      add = 0;
      for (int j = 0; j < dim; j++) {
        i_neighbor_corner = neighbors[j];
        if (refcorners[i_neighbor_corner] == 0) add = 1;
      }
      if (add) total++;
    }
  }
  return total;
}

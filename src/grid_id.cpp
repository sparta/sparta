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

#include "string.h"
#include "grid.h"
#include "error.h"

using namespace SPARTA_NS;

enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};         // same as Domain

// operations with grid cell IDs
// IMPORTANT: must be careful to use cellint (32 or 64 bit) in these cases
//   any cell ID, child or parent
//   any mask for bitwise operations
//   for ichild = 1 to Nxyz for subcell index within a parent cell
//   to multiply 2 or 3 nx,ny,nz indices

/* ----------------------------------------------------------------------
   compute indices of child cell within parent cell that contains point X
   lo/hi = corner points of parent cell
   nx,ny,nz = number of child cells within parent
   return ix,iy,iz for which child cell the point is in
   definition of inside is >= lo boundary and < hi boundary
   ix,iy,iz range from 0 to Nxyz-1 inclusive
------------------------------------------------------------------------- */

void Grid::id_point_child(double *x, double *lo, double *hi,
                          int nx, int ny, int nz, int &ix, int &iy, int &iz)
{
  // ix,iy,iz = child cell indices within parent lo/hi cell
  // inverse of master equation in id_child_lohi() for cell boundaries
  // for point on or eps from cell boundary, can produce round-off error

  ix = static_cast<int> ((x[0]-lo[0]) * nx/(hi[0]-lo[0]));
  iy = static_cast<int> ((x[1]-lo[1]) * ny/(hi[1]-lo[1]));
  iz = static_cast<int> ((x[2]-lo[2]) * nz/(hi[2]-lo[2]));

  // insure indices match grid cell boundaries in case of round-off error
  // via master equation id_child_lohi() that defines cell boundaries

  double edge;

  edge = lo[0] + ix*(hi[0]-lo[0])/nx;
  if (x[0] < edge) ix--;
  edge = lo[0] + (ix+1)*(hi[0]-lo[0])/nx;
  if (x[0] >= edge) ix++;

  edge = lo[1] + iy*(hi[1]-lo[1])/ny;
  if (x[1] < edge) iy--;
  edge = lo[1] + (iy+1)*(hi[1]-lo[1])/ny;
  if (x[1] >= edge) iy++;

  edge = lo[2] + iz*(hi[2]-lo[2])/nz;
  if (x[2] < edge) iz--;
  edge = lo[2] + (iz+1)*(hi[2]-lo[2])/nz;
  if (x[2] >= edge) iz++;

  // insure indices range from 0 to Nxyz-1 inclusive

  ix = MAX(ix,0);
  ix = MIN(ix,nx-1);
  iy = MAX(iy,0);
  iy = MIN(iy,ny-1);
  iz = MAX(iz,0);
  iz = MIN(iz,nz-1);
}

/* ----------------------------------------------------------------------
   calculate parentID of a childID at level
   return parentID
------------------------------------------------------------------------- */

cellint Grid::id_parent_of_child(cellint childID, int level)
{
  // mask = all 1s for parent bits of ID

  int parentbits = plevels[level-1].nbits;
  cellint mask = (1L << parentbits) - 1;
  cellint parentID = childID & mask;

  return parentID;
}

/* ----------------------------------------------------------------------
   find child cell within parentID which contains point X
   level = level of parent cell
   oplo/ophi = original parent cell corner pts
   pt X must be inside or on any boundary of parent cell
   recurse from parent downward until find a child cell or reach maxlevel
   if find child cell this proc stores (owned or ghost), return its local index
   else return -1 for unknown
------------------------------------------------------------------------- */

int Grid::id_find_child(cellint parentID, int plevel,
                        double *oplo, double *ophi, double *x)
{
  int ix,iy,iz,nx,ny,nz;
  double plo[3],phi[3],clo[3],chi[3];
  cellint childID,ichild;

  cellint id = parentID;
  int level = plevel;
  double *lo = oplo;
  double *hi = ophi;

  while (level < maxlevel) {
    nx = plevels[level].nx;
    ny = plevels[level].ny;
    nz = plevels[level].nz;

    id_point_child(x,lo,hi,nx,ny,nz,ix,iy,iz);
    ichild = (cellint) iz*nx*ny + (cellint) iy*nx + ix + 1;
    childID = (ichild << plevels[level].nbits) | id;

    if (hash->find(childID) != hash->end()) return (*hash)[childID];

    id = childID;
    id_child_lohi(level,lo,hi,ichild,clo,chi);
    plo[0] = clo[0]; plo[1] = clo[1]; plo[2] = clo[2];
    phi[0] = chi[0]; phi[1] = chi[1]; phi[2] = chi[2];
    lo = plo; hi = phi;
    level++;
  }

  return -1;
}

/* ----------------------------------------------------------------------
   compute cell ID of the child cell in a full-box uniform grid at level
   xyz grid = indices (0 to N-1) of the child cell within full-box grid
   return the cell ID
------------------------------------------------------------------------- */

cellint Grid::id_uniform_level(int level, int xgrid, int ygrid, int zgrid)
{
  int ix,iy,iz,nx,ny,nz;
  int parentbits;
  cellint ichild;

  cellint childID = 0;

  for (int ilevel = level-1; ilevel >= 0; ilevel--) {
    nx = plevels[ilevel].nx;
    ny = plevels[ilevel].ny;
    nz = plevels[ilevel].nz;

    ix = xgrid % nx;
    iy = ygrid % ny;
    iz = zgrid % nz;

    ichild = (cellint) iz*ny*nx + (cellint) iy*nx + ix + 1;
    parentbits = plevels[ilevel].nbits;
    childID |= ichild << parentbits;

    xgrid /= nx;
    ygrid /= ny;
    zgrid /= nz;
  }

  return childID;
}

/* ----------------------------------------------------------------------
   find child cell within a uniform grid at level which contains point X
     definition of "contains" depends on lohi flag
   level = level of uniform grid og child cells
   lohi = 0 for lower boundary, 1 for upper boundary
     for 0: if pt is on lower boundary (in any dim) of child cell,
            return index of adjacent cell (index is decremented)
     for 1: if pt is on upper boundary (in any dim) of child cell,
            return index of adjacent cell (index is incremented)
     incremented/decremented indices cannot be outside simulation box
   boxlo/boxhi = root level simulation box
   pt X must be inside or on any boundary of simulation box
   recurse downward to level
   return xyz grid = indices (0 to N-1) of child cell
     in full-box uniform grid at level
------------------------------------------------------------------------- */

void Grid::id_find_child_uniform_level(int level, int lohi,
                                       double *boxlo, double *boxhi, double *x,
                                       int &xgrid, int &ygrid, int &zgrid)
{
  int ix,iy,iz,nx,ny,nz;
  double plo[3],phi[3],clo[3],chi[3];
  cellint childID,ichild;

  cellint id = 0;
  int plevel = 0;
  double *lo = boxlo;
  double *hi = boxhi;
  xgrid = ygrid = zgrid = 0;

  while (plevel < level) {
    nx = plevels[plevel].nx;
    ny = plevels[plevel].ny;
    nz = plevels[plevel].nz;

    id_point_child(x,lo,hi,nx,ny,nz,ix,iy,iz);

    xgrid = nx*xgrid + ix;
    ygrid = ny*ygrid + iy;
    zgrid = nz*zgrid + iz;

    ichild = (cellint) iz*nx*ny + (cellint) iy*nx + ix + 1;
    childID = (ichild << plevels[plevel].nbits) | id;

    id = childID;
    id_child_lohi(plevel,lo,hi,ichild,clo,chi);
    plo[0] = clo[0]; plo[1] = clo[1]; plo[2] = clo[2];
    phi[0] = chi[0]; phi[1] = chi[1]; phi[2] = chi[2];
    lo = plo; hi = phi;
    plevel++;
  }

  // x >= lower edge and < upper edge of lo/hi cell with indices ix,iy,iz
  // if lohi = 0 and pt is on lower edge, decrement index
  // still require 0 <= index <= N-1

  if (lohi == 0) {
    if (x[0] == lo[0] && ix != 0) xgrid--;
    if (x[1] == lo[1] && iy != 0) ygrid--;
    if (x[2] == lo[2] && iz != 0) zgrid--;
  }
}

/* ----------------------------------------------------------------------
   calculate ID of neighbor cell of childID across its face
   neighbor cell is at same level
   assume neigbor cell has same parent
   return neighID if that is the case (doesn't check for periodic BC)
   return 0 if not same parent
------------------------------------------------------------------------- */

cellint Grid::id_neigh_same_parent(cellint id, int level, int face)
{
  int ixyz[3];
  int *nxyz;

  // mask = all 1s for parent bits of ID

  int parentbits = plevels[level-1].nbits;
  cellint mask = (1L << parentbits) - 1;
  cellint parent = id & mask;

  // mask = all 1s for child bits of ID

  int childbits = plevels[level-1].newbits;
  mask = (1L << childbits) - 1;
  cellint ichild = (id >> plevels[level-1].nbits) & mask;

  nxyz = &plevels[level-1].nx;

  ichild--;
  ixyz[0] = ichild % nxyz[0];
  ixyz[1] = (ichild / nxyz[0]) % nxyz[1];
  ixyz[2] = ichild / ((cellint) nxyz[0]*nxyz[1]);

  int dim = face / 2;
  if (face % 2) ixyz[dim]++;
  else ixyz[dim]--;
  if (ixyz[dim] < 0) return 0;
  if (ixyz[dim] == nxyz[dim]) return 0;


  ichild = (cellint) ixyz[2]*nxyz[1]*nxyz[0] +
    (cellint) ixyz[1]*nxyz[0] + ixyz[0] + 1;
  cellint neighID = (ichild << parentbits) | parent;

  return neighID;
}

/* ----------------------------------------------------------------------
   calculate ID of neighbor cell of cell ID across its face
   neighbor cell is at same level
   no assumption that neighbor cell has same parent
   assume periodic BC
   return neighID
------------------------------------------------------------------------- */

cellint Grid::id_neigh_same_level(cellint id, int level, int face)
{
  int ix,iy,iz,childbits,parentbits;
  cellint ichild,mask;
  int ixyz[3],fullxyz[3];
  int *nxyz;

  // ixyz[3] = indices of cell ID in fully resolved fullxyz grid at level
  // each index = 0 to N-1

  ixyz[0] = ixyz[1] = ixyz[2] = 0;
  fullxyz[0] = fullxyz[1] = fullxyz[2] = 1;

  for (int ilevel = 0; ilevel < level; ilevel++) {
    nxyz = &plevels[ilevel].nx;

    fullxyz[0] *= nxyz[0];
    fullxyz[1] *= nxyz[1];
    fullxyz[2] *= nxyz[2];

    ixyz[0] *= nxyz[0];
    ixyz[1] *= nxyz[1];
    ixyz[2] *= nxyz[2];

    // mask = all 1s for child bits at ilevel

    childbits = plevels[ilevel].newbits;
    mask = (1L << childbits) - 1;
    ichild = (id >> plevels[ilevel].nbits) & mask;

    ichild--;
    ixyz[0] += ichild % nxyz[0];
    ixyz[1] += (ichild / nxyz[0]) % nxyz[1];
    ixyz[2] += ichild / ((cellint) nxyz[0]*nxyz[1]);
  }

  // alter ixyz[3] to be indices of neighbor cell
  // apply periodic boundary conditions based on fullxyz[3]
  // each altered index = 0 to N-1

  int dim = face / 2;
  if (face % 2) ixyz[dim]++;
  else ixyz[dim]--;
  if (ixyz[dim] < 0) ixyz[dim] = fullxyz[dim] - 1;
  if (ixyz[dim] == fullxyz[dim]) ixyz[dim] = 0;

  // neighID = neighbor cell ID generated from ixyz[3]

  cellint neighID = 0;

  for (int ilevel = level-1; ilevel >= 0; ilevel--) {
    nxyz = &plevels[ilevel].nx;

    ix = ixyz[0] % nxyz[0];
    iy = ixyz[1] % nxyz[1];
    iz = ixyz[2] % nxyz[2];

    ichild = (cellint) iz*nxyz[1]*nxyz[0] + (cellint) iy*nxyz[0] + ix + 1;
    parentbits = plevels[ilevel].nbits;
    neighID |= ichild << parentbits;

    ixyz[0] /= nxyz[0];
    ixyz[1] /= nxyz[1];
    ixyz[2] /= nxyz[2];
  }

  return neighID;
}

/* ----------------------------------------------------------------------
   calculate a child ID of parent cell at level
   choose child in center of face of parent
   return child ID
------------------------------------------------------------------------- */

cellint Grid::id_refine(cellint parentID, int plevel, int face)
{
  int *nxyz = &plevels[plevel].nx;

  int ixyz[3];
  ixyz[0] = nxyz[0] / 2;
  ixyz[1] = nxyz[1] / 2;
  ixyz[2] = nxyz[2] / 2;

  if (face % 2) ixyz[face/2] = nxyz[face/2] - 1;
  else ixyz[face/2] = 0;

  cellint ichild = (cellint) ixyz[2]*nxyz[1]*nxyz[0] +
    (cellint) ixyz[1]*nxyz[0] + ixyz[0] + 1;

  int parentbits = plevels[plevel].nbits;
  cellint childID = parentID | (ichild << parentbits);
  return childID;
}

/* ----------------------------------------------------------------------
   calculate parent ID of child cell at level
   return parent ID
------------------------------------------------------------------------- */

cellint Grid::id_coarsen(cellint childID, int level)
{
  // mask = all 1s for parent bits of child

  int parentbits = plevels[level-1].nbits;
  cellint mask = (1L << parentbits) - 1;
  cellint parentID = childID & mask;
  return parentID;
}

/* ----------------------------------------------------------------------
   calculate which ichild the childID of parentID is
   plevel = grid level of parent, not child
   return ichild = 1 to Nxyz
   return a cellint in case Nxyz is huge, e.g. at level 1
------------------------------------------------------------------------- */

cellint Grid::id_ichild(cellint childID, cellint parentID, int plevel)
{
  int parentbits = plevels[plevel].nbits;
  cellint ichild = childID >> parentbits;
  return ichild;
}

/* ----------------------------------------------------------------------
   return level the cell ID is at from 1 to Nlevels
   calculate by recursing from root until id = 0
   if id still not 0 at maxlevel, return -1
------------------------------------------------------------------------- */

int Grid::id_level(cellint id)
{
  int newbits;
  cellint mask,ichild;

  int level = 0;
  while (level < maxlevel) {
    newbits = plevels[level].newbits;
    mask = (1L << newbits) - 1;
    ichild = id & mask;
    id = id >> newbits;
    if (!id) break;
    level++;
  }

  return level+1;
}

/* ----------------------------------------------------------------------
   compute lo/hi extent of a specific child cell within a parent cell
   plevel = level of parent
   plo/phi = parent cell corner points
   ichild ranges from 1 to Nx*Ny*Nz within parent cell
   return clo/chi corner points, caller must allocate them
------------------------------------------------------------------------- */

void Grid::id_child_lohi(int plevel, double *plo, double *phi,
                         cellint ichild, double *clo, double *chi)
{
  int nx = plevels[plevel].nx;
  int ny = plevels[plevel].ny;
  int nz = plevels[plevel].nz;

  ichild--;
  int ix = ichild % nx;
  int iy = (ichild/nx) % ny;
  int iz = ichild / ((cellint) nx*ny);

  clo[0] = plo[0] + ix*(phi[0]-plo[0])/nx;
  clo[1] = plo[1] + iy*(phi[1]-plo[1])/ny;
  clo[2] = plo[2] + iz*(phi[2]-plo[2])/nz;

  chi[0] = plo[0] + (ix+1)*(phi[0]-plo[0])/nx;
  chi[1] = plo[1] + (iy+1)*(phi[1]-plo[1])/ny;
  chi[2] = plo[2] + (iz+1)*(phi[2]-plo[2])/nz;

  if (ix == nx-1) chi[0] = phi[0];
  if (iy == ny-1) chi[1] = phi[1];
  if (iz == nz-1) chi[2] = phi[2];
}

/* ----------------------------------------------------------------------
   compute lo/hi extent of cell ID at level
   boxlo/boxhi = global simulation box (level 0)
   return lo/hi corner points, caller must allocate them
------------------------------------------------------------------------- */

void Grid::id_lohi(cellint id, int level, double *boxlo, double *boxhi,
                   double *lo, double *hi)
{
  int childbits;
  cellint ichild,mask;
  double plo[3],phi[3];

  plo[0] = boxlo[0]; plo[1] = boxlo[1]; plo[2] = boxlo[2];
  phi[0] = boxhi[0]; phi[1] = boxhi[1]; phi[2] = boxhi[2];

  for (int ilevel = 0; ilevel < level; ilevel++) {

    // mask = all 1s for child bits at ilevel

    childbits = plevels[ilevel].newbits;
    mask = (1L << childbits) - 1;
    ichild = (id >> plevels[ilevel].nbits) & mask;

    id_child_lohi(ilevel,plo,phi,ichild,lo,hi);

    plo[0] = lo[0]; plo[1] = lo[1]; plo[2] = lo[2];
    phi[0] = hi[0]; phi[1] = hi[1]; phi[2] = hi[2];
  }
}

/* ----------------------------------------------------------------------
   calculate # of bits needed to store values from 1 to Nx*Ny*Nz
------------------------------------------------------------------------- */

int Grid::id_bits(int nx, int ny, int nz)
{
  bigint n = (bigint) nx*ny*nz;
  bigint nstore = 1;
  int nbits = 1;

  while (nstore < n) {
    nstore = 2*nstore + 1;
    nbits++;
  }

  return nbits;
}

/* ----------------------------------------------------------------------
   convert cell ID from number to string with levels separated by dashes
   walk grid hierarchy from root to ID, generating one substr at each level
   NOTE: should append letter for sub cells, but not enough info to do it
         don't know the parent cell (which might be on a different proc)
         one way to do it would be for dump_grid to have a pack_idstr
         which sets some hi-bits to encode the sub cell #
         but this would use up some extra bits in the ID
------------------------------------------------------------------------- */

void Grid::id_num2str(cellint id, char *str)
{
  int newbits;
  cellint mask,ichild;

  // sprintf child bits one level at a time until id = 0

  int offset = 0;
  int level = 0;

  while (1) {
    newbits = plevels[level].newbits;
    mask = (1L << newbits) - 1;
    ichild = id & mask;
    sprintf(&str[offset],"%d",ichild);
    offset = strlen(str);
    id = id >> newbits;
    if (!id) return;
    str[offset++] = '-';
    level++;
  }
}

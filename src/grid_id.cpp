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
#include "grid.h"
#include "error.h"

using namespace SPARTA_NS;

enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};         // same as Domain

// operations with grid cell IDs

/* ----------------------------------------------------------------------
   find child cell within parentID which contains pt X
   level = level of parent cell
   oplo/ophi = original parent cell corner pts
   pt X can be inside or on any boundary of parent cell
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
    ix = static_cast<int> ((x[0]-lo[0]) * nx/(hi[0]-lo[0]));
    iy = static_cast<int> ((x[1]-lo[1]) * ny/(hi[1]-lo[1]));
    iz = static_cast<int> ((x[2]-lo[2]) * nz/(hi[2]-lo[2]));
    if (ix == nx) ix--;
    if (iy == ny) iy--;
    if (iz == nz) iz--;

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
   find parent of a child or parent ID
   loop from root thru parents until match the input ID
   return local index of parent cell
   also return ichild = 1 to Nx*Ny*Nz for index of child within parent
   ichild is cellint in case Nx*Ny*Nz exceeds 32-bit int
------------------------------------------------------------------------- */

int Grid::id_find_parent(cellint id, cellint &ichild)
{
  // TMP
  /*
  int level,nbits,newbits,index;
  cellint idparent,idnew,mask;
  ParentCell *p;

  int iparent = 0;
  while (1) {
    p = &pcells[iparent];
    idparent = p->id;
    level = p->level;
    nbits = level_nbits[level];
    newbits = level_newbits[level];

    // ichild = the newbits above nbits in id

    mask = (1L << newbits) - 1;
    ichild = (id >> nbits) & mask;
    idnew = idparent | (ichild << nbits);
    if (idnew == id) break;
    if (hash->find(idnew) == hash->end())
      error->one(FLERR,"Grid::id_find_new() not in hash");
    index = (*hash)[idnew];
    iparent = -index-1;
  }

  return iparent;
  */
  
  return 0;
}

/* ----------------------------------------------------------------------
   convert cell ID from string to number and return it
   idstr = "0" is special case, return 0
   otherwise, walk hierarchy from root to child so can shift bits correctly
   return 0 if error:
     parent cell does not exist
     field exceeds bit count of parent
------------------------------------------------------------------------- */

cellint Grid::id_str2num(char *idstr)
{
  // TMP
  /*
  int iparent,level;
  cellint ichild,nxyz;
  ParentCell *p;
  
  if (idstr[0] == '0') return 0;

  cellint id = 0;
  char *word = idstr;
  char *ptr = strchr(word,'-');

  while (1) {
    if (hash->find(id) == hash->end()) return -1;
    iparent = (*hash)[id];
    p = &pcells[iparent];
    if (ptr) *ptr = '\0';
    ichild = ATOCELLINT(word);
    level = p->level;
    nxyz = (cellint) level_xyz[level][0] * 
      level_xyz[level][0] * level_xyz[level][0];
    if (ichild == 0 || ichild > nxyz) {
      if (ptr) *ptr = '-';
      return 0;
    }
    id |= ichild << level_nbits[level];
    if (!ptr) break;
    *ptr = '-';
    word = ptr+1;
    ptr = strchr(word,'-');
  }

  return id;
  */

  return 0;
}

/* ----------------------------------------------------------------------
   convert cell ID from number to string
   walk hierarchy from root to child so can shift bits correctly
   NOTE: should append letter for sub cells, but not enough info to do it
         don't know the parent cell (which might be on a different proc)
         one way to do it would be for dump_grid to have a pack_idstr
         which sets some hi-bits to encode the sub cell #
         but this would use up some extra bits in the ID
------------------------------------------------------------------------- */

void Grid::id_num2str(cellint id, char *str)
{
  // TMP
  /*
  int index,level;
  cellint mask,idlevel;

  cellint idparent = 0;
  int iparent = 0;
  int offset = 0;

  while (1) {
    level = pcells[iparent].level;
    mask = (1L << level_newbits[level]) - 1;
    idlevel = id & mask;
    sprintf(&str[offset],"%d",idlevel);
    offset = strlen(str);
    id = id >> level_newbits[level];
    if (!id) return;
    str[offset++] = '-';
    idparent = idparent | (idlevel << level_nbits[level]);
    index = (*hash)[idparent];
    iparent = -index-1;
  }
  */
}

/* ----------------------------------------------------------------------
   split cell ID string into parent ID and child ID strings
   if ID = "0", return parent as empty string
   if ID has only a child field, return "0" as parent
------------------------------------------------------------------------- */

void Grid::id_pc_split(char *idstr, char *parent, char *child)
{
  if (idstr[0] == '0') {
    parent[0] = '\0';
    strcpy(child,idstr);
    return;
  }

  char *ptr = strrchr(idstr,'-');

  if (!ptr) {
    strcpy(parent,"0");
    strcpy(child,idstr);
    return;
  }

  *ptr = '\0';
  strcpy(parent,idstr);
  strcpy(child,ptr+1);
  *ptr = '-';
  return;
}

/* ----------------------------------------------------------------------
   return index of child cell that is in icorner of iparent cell
   recurse down thru grid hierarchy until find child cell
   hash contains all owned/ghost child cells and parent cells
   if I don't store child cell as owned or ghost, return -1 for unknown
------------------------------------------------------------------------- */

int Grid::id_child_from_parent_corner(int iparent, int icorner)
{
  // TMP
  /*
  ParentCell *p = &pcells[iparent];

  int level = p->level;
  int nx = level_xyz[level][0];
  int ny = level_xyz[level][1];
  int nz = level_xyz[level][2];

  int ix = (icorner % 2) ? nx-1 : 0;
  int iy = ((icorner/2) % 2) ? ny-1 : 0;
  int iz = (icorner / 4) ? nz-1 : 0;
  
  cellint ichild = (cellint) iz*nx*ny + (cellint) iy*nx + ix + 1;
  cellint idchild = p->id | (ichild << level_nbits[level]);

  if (hash->find(idchild) == hash->end()) return -1;
  int index = (*hash)[idchild];
  if (index > 0) return index-1;
  return id_child_from_parent_corner(-index-1,icorner);
  */

  return 0;
}

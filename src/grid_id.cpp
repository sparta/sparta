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

// operations with grid cell IDs

/* ----------------------------------------------------------------------
   find child cell within iparent cell assumed to be contain pt x
   recurse down thru parents until reach a child cell
   pt can be on any boundary of parent cell
   if I don't store child cell as owned or ghost, return -1 for unknown
   else return local index of child cell
------------------------------------------------------------------------- */

int Grid::id_find_child(int iparent, double *x)
{
  while (1) {
    ParentCell *p = &pcells[iparent];
    double *lo = p->lo;
    double *hi = p->hi;
    int nx = p->nx;
    int ny = p->ny;
    int nz = p->nz;
    int ix = static_cast<int> ((x[0]-lo[0]) * nx/(hi[0]-lo[0]));
    int iy = static_cast<int> ((x[1]-lo[1]) * ny/(hi[1]-lo[1]));
    int iz = static_cast<int> ((x[2]-lo[2]) * nz/(hi[2]-lo[2]));
    if (ix == nx) ix--;
    if (iy == ny) iy--;
    if (iz == nz) iz--;

    cellint ichild = (cellint) iz*nx*ny + (cellint) iy*nx + ix + 1;
    cellint idchild = p->id | (ichild << p->nbits);

    if (hash->find(idchild) == hash->end()) return -1;
    int index = (*hash)[idchild];
    if (index > 0) return index-1;
    iparent = -index-1;
  }
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
  int nbits,newbits,index;
  cellint idparent,idnew,mask;
  ParentCell *p;

  int iparent = 0;
  while (1) {
    p = &pcells[iparent];
    idparent = p->id;
    nbits = p->nbits;
    newbits = p->newbits;

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
}

/* ----------------------------------------------------------------------
   convert cell ID from string to number and return it
   idstr = "0" is special case, return 0
   otherwise, walk hierarchy from root to child so can shift bits correctly
   return -1 if error:
     parent cell does not exist
     field exceeds bit count of parent
------------------------------------------------------------------------- */

cellint Grid::id_str2num(char *idstr)
{
  if (idstr[0] == '0') return 0;

  cellint id = 0;
  char *word = idstr;
  char *ptr = strchr(word,'-');

  while (1) {
    if (hash->find(id) == hash->end()) return -1;
    int iparent = (*hash)[id];
    ParentCell *p = &pcells[iparent];
    if (ptr) *ptr = '\0';
    cellint ichild = ATOCELLINT(word);
    if (ichild == 0 || ichild > ((cellint) p->nx*p->ny*p->nz)) {
      if (ptr) *ptr = '-';
      return -1;
    }
    id |= ichild << p->nbits;
    if (!ptr) break;
    *ptr = '-';
    word = ptr+1;
    ptr = strchr(word,'-');
  }

  return id;
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
  cellint mask,idlevel;

  cellint idparent = 0;
  int iparent = 0;
  int offset = 0;

  while (1) {
    mask = (1L << pcells[iparent].newbits) - 1;
    idlevel = id & mask;
    sprintf(&str[offset],"%d",idlevel);
    offset = strlen(str);
    id = id >> pcells[iparent].newbits;
    if (!id) return;
    str[offset++] = '-';
    idparent = idparent | (idlevel << pcells[iparent].nbits);
    int index = (*hash)[idparent];
    iparent = -index-1;
  }
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
   compute lo/hi extent of a specific child cell within a parent cell
   iparent = index of parent cell, has Nx by Ny by Nz children
   ichild ranges from 1 to Nx*Ny*Nz
------------------------------------------------------------------------- */

void Grid::id_child_lohi(int iparent, cellint ichild, double *lo, double *hi)
{
  ParentCell *p = &pcells[iparent];
  ichild--;

  int nx = p->nx;
  int ny = p->ny;
  int nz = p->nz;

  int ix = ichild % nx;
  int iy = (ichild/nx) % ny;
  int iz = ichild / ((bigint) nx*ny);

  double *plo = p->lo;
  double *phi = p->hi;

  lo[0] = plo[0] + ix*(phi[0]-plo[0])/nx;
  lo[1] = plo[1] + iy*(phi[1]-plo[1])/ny;
  lo[2] = plo[2] + iz*(phi[2]-plo[2])/nz;

  hi[0] = plo[0] + (ix+1)*(phi[0]-plo[0])/nx;
  hi[1] = plo[1] + (iy+1)*(phi[1]-plo[1])/ny;
  hi[2] = plo[2] + (iz+1)*(phi[2]-plo[2])/nz;

  if (ix == nx-1) hi[0] = phi[0];
  if (iy == ny-1) hi[1] = phi[1];
  if (iz == nz-1) hi[2] = phi[2];
}

/* ----------------------------------------------------------------------
   determine # of bits needed to store values from 1 to Nx*Ny*Nz
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
   recurse to find smallest child or parent cell that contains pt x and 
     fully covers face in dim of olo/ohi cell
   only called with x that is inside simulation box
   icell = parent cell index
   return ID of child or parent cell
------------------------------------------------------------------------- */

cellint Grid::id_find_face(double *x, int icell, int dim, 
                           double *olo, double *ohi)
{
  ParentCell *p = &pcells[icell];
  double *lo = p->lo;
  double *hi = p->hi;

  // go down a level to find (ix,iy,iz) of new cell that contains pt x

  int nx = p->nx;
  int ny = p->ny;
  int nz = p->nz;
  int ix = static_cast<int> ((x[0]-lo[0]) * nx/(hi[0]-lo[0]));
  int iy = static_cast<int> ((x[1]-lo[1]) * ny/(hi[1]-lo[1]));
  int iz = static_cast<int> ((x[2]-lo[2]) * nz/(hi[2]-lo[2]));
  if (ix == nx) ix--;
  if (iy == ny) iy--;
  if (iz == nz) iz--;

  // calculate lo/hi of new cell
  // exact same math as in id_child_lohi()

  double newlo[3],newhi[3];

  newlo[0] = lo[0] + ix*(hi[0]-lo[0])/nx;
  newlo[1] = lo[1] + iy*(hi[1]-lo[1])/ny;
  newlo[2] = lo[2] + iz*(hi[2]-lo[2])/nz;

  newhi[0] = lo[0] + (ix+1)*(hi[0]-lo[0])/nx;
  newhi[1] = lo[1] + (iy+1)*(hi[1]-lo[1])/ny;
  newhi[2] = lo[2] + (iz+1)*(hi[2]-lo[2])/nz;

  if (ix == nx-1) newhi[0] = hi[0];
  if (iy == ny-1) newhi[1] = hi[1];
  if (iz == nz-1) newhi[2] = hi[2];

  // if new cell does not fully overlap olo/ohi face, return parent ID

  if (dim != 0 && (newlo[0] > olo[0] || newhi[0] < ohi[0])) return p->id;
  if (dim != 1 && (newlo[1] > olo[1] || newhi[1] < ohi[1])) return p->id;
  if (dim != 2 && (newlo[2] > olo[2] || newhi[2] < ohi[2])) return p->id;

  // id = ID of new cell
  // if I don't store new ID, it's a child ID, return it
  // if I do store new ID, determine if parent or child cell
  // if child, return it
  // if parent, recurse

  cellint ichild = (cellint) iz*nx*ny + (cellint) iy*nx + ix + 1;
  cellint id = p->id | (ichild << p->nbits);
  if (hash->find(id) == hash->end()) return id;
  icell = (*hash)[id];
  if (icell > 0) return id;
  return id_find_face(x,-icell-1,dim,olo,ohi);
}

/* ----------------------------------------------------------------------
   return index of child cell that is in icorner of iparent cell
   recurse down thru grid hierarchy until find child cell
   hash contains all owned/ghost child cells and parent cells
   if I don't store child cell as owned or ghost, return -1 for unknown
------------------------------------------------------------------------- */

int Grid::id_child_from_parent_corner(int iparent, int icorner)
{
  ParentCell *p = &pcells[iparent];

  int ix = (icorner % 2) ? p->nx-1 : 0;
  int iy = ((icorner/2) % 2) ? p->ny-1 : 0;
  int iz = (icorner / 4) ? p->nz-1 : 0;
  
  cellint ichild = (cellint) iz*p->nx*p->ny + (cellint) iy*p->nx + ix + 1;
  cellint idchild = p->id | (ichild << p->nbits);

  if (hash->find(idchild) == hash->end()) return -1;
  int index = (*hash)[idchild];
  if (index > 0) return index-1;
  return id_child_from_parent_corner(-index-1,icorner);
}

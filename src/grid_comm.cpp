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
#include "particle.h"
#include "collide.h"
#include "modify.h"

// DEBUG
#include "comm.h"

using namespace SPARTA_NS;

// grid cell communication

/* ----------------------------------------------------------------------
   pack single icell into buf
   only called for owned unsplit and split cells, not for sub cells
   include its cinfo, surfs, split info, sub cells if necessary
   ownflag = 1/0 = owned or ghost cell
     for owned cell, also pack cinfo and auxiliary collision/fix info
     if nsplit < 0, is an empty ghost cell, pack only cells
   molflag = 0/1 = no/yes to also pack particles
   memflag = 0/1 = no/yes to actually pack into buf, 0 = just length
   return length of packing in bytes
------------------------------------------------------------------------- */

int Grid::pack_one(int icell, char *buf, 
                   int ownflag, int molflag, int memflag)
{
  char *ptr = buf;

  // pack child cell and its csurfs and cinfo
  // just pack child cell if nsurf < 0 since is EMPTY ghost

  if (memflag) memcpy(ptr,&cells[icell],sizeof(ChildCell));
  ptr += sizeof(ChildCell);
  ptr = ROUNDUP(ptr);

  if (cells[icell].nsurf < 0) return ptr - buf;

  if (cells[icell].nsurf) {
    int nsurf = cells[icell].nsurf;
    if (memflag) memcpy(ptr,cells[icell].csurfs,nsurf*sizeof(int));
    ptr += nsurf*sizeof(int);
    ptr = ROUNDUP(ptr);
  }

  if (ownflag) {
    if (memflag) memcpy(ptr,&cinfo[icell],sizeof(ChildInfo));
    ptr += sizeof(ChildInfo);
    ptr = ROUNDUP(ptr);
  }

  // if split cell, pack sinfo and sinfo.csplits and sinfo.csubs
  // but not sub cells themselves

  if (cells[icell].nsplit > 1) {
    int isplit = cells[icell].isplit;
    if (memflag) memcpy(ptr,&sinfo[isplit],sizeof(SplitInfo));
    ptr += sizeof(SplitInfo);
    ptr = ROUNDUP(ptr);

    int nsurf = cells[icell].nsurf;
    if (memflag) memcpy(ptr,sinfo[isplit].csplits,nsurf*sizeof(int));
    ptr += nsurf*sizeof(int);
    ptr = ROUNDUP(ptr);

    int nsplit = cells[icell].nsplit;
    if (memflag) memcpy(ptr,sinfo[isplit].csubs,nsplit*sizeof(int));
    ptr += nsplit*sizeof(int);
    ptr = ROUNDUP(ptr);
  }

  // pack collision and fix info for migrating cell

  if (ownflag) {
    if (collide) {
      ptr += collide->pack_grid_one(icell,ptr,memflag);
      ptr = ROUNDUP(ptr);
    }
    if (modify->n_pergrid) {
      ptr += modify->pack_grid_one(icell,ptr,memflag);
      ptr = ROUNDUP(ptr);
    }
  }

  // pack particles for unsplit cell or split cell

  if (!molflag) return ptr - buf;

  ncustom = particle->ncustom;
  nbytes_particle = sizeof(Particle::OnePart);
  nbytes_custom = particle->sizeof_custom();
  nbytes_total = nbytes_particle + nbytes_custom;

  ptr += pack_particles(icell,ptr,memflag);

  // pack particles of sub cells

  if (cells[icell].nsplit > 1) {
    int isplit = cells[icell].isplit;
    int nsplit = cells[icell].nsplit;
    for (int i = 0; i < nsplit; i++) {
      int m = sinfo[isplit].csubs[i];
      ptr += pack_particles(m,ptr,memflag);
    }
  }

  return ptr - buf;
}

/* ----------------------------------------------------------------------
   static callback function for ring communication
   uppack received ghost cells into cell lists
------------------------------------------------------------------------- */

void Grid::unpack_ghosts(int nsize, char *buf)
{
  int n = 0;
  while (n < nsize)
    n += gptr->unpack_one(&buf[n],0,0);
}

/* ----------------------------------------------------------------------
   unpack single icell from buf
   include its cinfo, surfs, split info, sub cells if necessary
   ownflag = 1/0 = owned or ghost cell
     for owned cell, also unpack cinfo and auxiliary collision/fix info
   molflag = 0/1 = no/yes to also unpack particles
   return length of unpacking in bytes
------------------------------------------------------------------------- */

int Grid::unpack_one(char *buf, int ownflag, int molflag)
{
  char *ptr = buf;

  // unpack child cell and its csurfs and cinfo, add to ghost cells

  int icell;
  if (ownflag) icell = nlocal;
  else icell = nlocal + nghost;
  grow_cells(1,ownflag);
  if (ownflag) nlocal++;
  else nghost++;

  memcpy(&cells[icell],ptr,sizeof(ChildCell));
  ptr += sizeof(ChildCell);
  ptr = ROUNDUP(ptr);

  if (cells[icell].nsurf < 0) return ptr - buf;

  if (ownflag) {
    cells[icell].proc = me;
    cells[icell].ilocal = icell;
  }

  if (cells[icell].nsurf) {
    int nsurf = cells[icell].nsurf;
    cells[icell].csurfs = csurfs->vget();
    memcpy(cells[icell].csurfs,ptr,nsurf*sizeof(int));
    csurfs->vgot(nsurf);
    ptr += nsurf*sizeof(int);
    ptr = ROUNDUP(ptr);
  }

  if (ownflag) {
    memcpy(&cinfo[icell],ptr,sizeof(ChildInfo));
    ptr += sizeof(ChildInfo);
    ptr = ROUNDUP(ptr);
  }

  // if split cell, unpack sinfo and sinfo.csplits and sinfo.csubs
  // create Nsplit sub cells
  // use sinfo.csubs to set cells.ilocal for new sub cells
  // create new csub for new sub cell indices

  if (cells[icell].nsplit > 1) {
    int isplit;
    if (ownflag) isplit = nsplitlocal;
    else isplit = nsplitlocal + nsplitghost;
    cells[icell].isplit = isplit;
    add_split_cell(ownflag);
    memcpy(&sinfo[isplit],ptr,sizeof(SplitInfo));
    ptr += sizeof(SplitInfo);
    ptr = ROUNDUP(ptr);

    sinfo[isplit].icell = icell;
    int nsurf = cells[icell].nsurf;
    sinfo[isplit].csplits = csplits->vget();
    memcpy(sinfo[isplit].csplits,ptr,nsurf*sizeof(int));
    csplits->vgot(nsurf);
    ptr += nsurf*sizeof(int);
    ptr = ROUNDUP(ptr);

    int nsplit = cells[icell].nsplit;
    sinfo[isplit].csubs = csubs->vget();
    memcpy(sinfo[isplit].csubs,ptr,nsplit*sizeof(int));
    csubs->vgot(nsplit);
    ptr += nsplit*sizeof(int);
    ptr = ROUNDUP(ptr);

    int isub;
    for (int i = 0; i < nsplit; i++) {
      if (ownflag) isub = nlocal;
      else isub = nlocal + nghost;
      add_sub_cell(icell,ownflag);
      cells[isub].ilocal = sinfo[isplit].csubs[i];
      cells[isub].nsplit = -i;
      sinfo[isplit].csubs[i] = isub;
    }

  } else {
    if (ownflag) nunsplitlocal++;
    else nunsplitghost++;
  }

  // unpack collision and fix info for new grid cell

  if (ownflag) {
    if (collide) {
      ptr += collide->unpack_grid_one(icell,ptr);
      ptr = ROUNDUP(ptr);
    }
    if (modify->n_pergrid) {
      ptr += modify->unpack_grid_one(icell,ptr);
      ptr = ROUNDUP(ptr);
    }
  }
  
  // unpack particles, for unsplit cell or split cell

  if (!molflag) return ptr - buf;

  ptr += unpack_particles(ptr,icell);

  // unpack particles of sub cells

  if (cells[icell].nsplit > 1) {
    int isplit = cells[icell].isplit;
    int nsplit = cells[icell].nsplit;
    for (int i = 0; i < nsplit; i++) {
      int m = sinfo[isplit].csubs[i];
      ptr += unpack_particles(ptr,m);
    }
  }

  return ptr - buf;
}

/* ----------------------------------------------------------------------
   pack particles of one cell = icell into buf
   memflag = 0/1 = no/yes to perform actual packing into buf, 0 = just length
   flag packed particles with icell = -1, so can delete via particle->compress()
   return length of packing in bytes
------------------------------------------------------------------------- */

int Grid::pack_particles(int icell, char *buf, int memflag)
{
  char *ptr = buf;

  int np = cinfo[icell].count;
  if (memflag) *((int *) ptr) = np;
  ptr += sizeof(int);
  ptr = ROUNDUP(ptr);
  if (!np) return ptr-buf;

  if (memflag) {
    Particle::OnePart *particles = particle->particles;
    int *next = particle->next;
    int ip = cinfo[icell].first;

    while (ip >= 0) {
      memcpy(ptr,&particles[ip],nbytes_particle);
      ptr += nbytes_particle;
      if (ncustom) {
        particle->pack_custom(ip,ptr);
        ptr += nbytes_custom;
      }
      particles[ip].icell = -1;
      ip = next[ip];
    }
  } else ptr += np * nbytes_total;

  ptr = ROUNDUP(ptr);
  return ptr-buf;
}

/* ----------------------------------------------------------------------
   unpack particles from buf into one cell = icell
   point unpacked particles to icell
   return length of unpacking in bytes
------------------------------------------------------------------------- */

int Grid::unpack_particles(char *buf, int icell)
{
  char *ptr = buf;

  int np = *((int *) ptr);

  ptr += sizeof(int);
  ptr = ROUNDUP(ptr);
  if (!np) return ptr-buf;

  particle->grow(np);

  Particle::OnePart *particles = particle->particles;
  int nplocal = particle->nlocal;

  if (ncustom) {
    int n = nplocal;
    for (int i = 0; i < np; i++) {
      memcpy(&particles[n],ptr,nbytes_particle);
      ptr += nbytes_particle;
      particle->unpack_custom(ptr,n);
      ptr += nbytes_custom;
      n++;
    }
  } else {
    memcpy(&particles[nplocal],ptr,np*nbytes_particle);
    ptr += np*nbytes_particle;
  }

  ptr = ROUNDUP(ptr);

  int npnew = nplocal + np;
  for (int i = nplocal; i < npnew; i++) particles[i].icell = icell;
  particle->nlocal = npnew;

  return ptr-buf;
}

/* ----------------------------------------------------------------------
   commpress owned grid cells due to some cells migrating to other procs
   operates on cells, cinfo, sinfo
   integer lists (csurfs,csplits,csubs) are re-created from scratch
------------------------------------------------------------------------- */

void Grid::compress()
{
  // must compress per-cell arrays in collide and fixes before
  // cells data structure changes

  if (collide) collide->compress_grid();
  if (modify->n_pergrid) modify->compress_grid(0);

  // copy of integer lists
  // create new lists

  MyPage<int> *csurfs_old = csurfs;
  MyPage<int> *csplits_old = csplits;
  MyPage<int> *csubs_old = csubs;

  csurfs = NULL; csplits = NULL; csubs = NULL;
  allocate_surf_arrays();

  // compress cells and cinfo and sinfo arrays
  // 4 cases:
  // discard unsplit or sub cell:
  //   continue
  // keep unsplit or sub cell:
  //   memcpy of cells and cinfo, setup csurfs
  //   reset cells.ilocal and sinfo.csubs
  //   increment nlocal, nunsplitlocal, nsublocal
  // discard split cell:
  //   flag all sub cells with other proc
  //   continue
  // keep split cell:
  //   memcpy of cells and cinfo and sinfo, setup csurfs/csplits/csubs
  //   reset sinfo.icell
  //   reset cells.isplit, ditto for sub cells
  //   increment nlocal, nsplitlocal
  // relies on sub cells appearing in cells array after their split cell

  int ncurrent = nlocal;
  nlocal = nunsplitlocal = nsplitlocal = nsublocal = 0;

  for (int icell = 0; icell < ncurrent; icell++) {

    // unsplit and sub cells

    if (cells[icell].nsplit <= 1) {
      if (cells[icell].proc != me) continue;

      if (icell != nlocal) {
        memcpy(&cells[nlocal],&cells[icell],sizeof(ChildCell));
        memcpy(&cinfo[nlocal],&cinfo[icell],sizeof(ChildInfo));
      }

      cells[nlocal].ilocal = nlocal;

      // for unsplit cell, new copy of all csurfs indices
      // for sub cell, just copy csurfs ptr from its split cell

      if (cells[nlocal].nsurf) {
        if (cells[nlocal].nsplit == 1) {
          int *oldcsurfs = cells[icell].csurfs;
          cells[nlocal].csurfs = csurfs->vget();
          memcpy(cells[nlocal].csurfs,oldcsurfs,
                 cells[nlocal].nsurf*sizeof(int));
          csurfs->vgot(cells[nlocal].nsurf);
        } else
          cells[nlocal].csurfs = 
            cells[sinfo[cells[nlocal].isplit].icell].csurfs;
      }

      if (cells[nlocal].nsplit <= 0) {
        sinfo[cells[nlocal].isplit].csubs[-cells[nlocal].nsplit] = nlocal;
        nsublocal++;
      } else nunsplitlocal++;

      nlocal++;

    // split cells

    } else {
      if (cells[icell].proc != me) {
        int isplit = cells[icell].isplit;
        int nsplit = cells[icell].nsplit;
        for (int i = 0; i < nsplit; i++) {
          int m = sinfo[isplit].csubs[i];
          cells[m].proc = cells[icell].proc;
        }
        continue;
      }

      if (icell != nlocal) {
        memcpy(&cells[nlocal],&cells[icell],sizeof(ChildCell));
        memcpy(&cinfo[nlocal],&cinfo[icell],sizeof(ChildInfo));
      }

      cells[nlocal].ilocal = nlocal;
      
      // new copy of all csurfs indices

      int *oldcsurfs = cells[icell].csurfs;
      cells[nlocal].csurfs = csurfs->vget();
      memcpy(cells[nlocal].csurfs,oldcsurfs,
             cells[nlocal].nsurf*sizeof(int));
      csurfs->vgot(cells[nlocal].nsurf);

      // compress sinfo

      int isplit = cells[nlocal].isplit;
      if (isplit != nsplitlocal)
        memcpy(&sinfo[nsplitlocal],&sinfo[isplit],sizeof(SplitInfo));
      cells[nlocal].isplit = nsplitlocal;
      sinfo[nsplitlocal].icell = nlocal;

      // new copy of all csplits indices

      int *oldcsplits = sinfo[isplit].csplits;
      sinfo[nsplitlocal].csplits = csplits->vget();
      memcpy(sinfo[nsplitlocal].csplits,oldcsplits,
             cells[nlocal].nsurf*sizeof(int));
      csplits->vgot(cells[nlocal].nsurf);

      // new csubs list of length nsplit
      // values unset for now, set one-by-one when sub cells are compressed

      int *oldcsubs = sinfo[isplit].csubs;
      sinfo[nsplitlocal].csubs = csubs->vget();
      csubs->vgot(cells[nlocal].nsplit);

      // point each sub cell in old sinfo csubs at new compressed sinfo index

      int nsplit = cells[nlocal].nsplit;
      for (int i = 0; i < nsplit; i++)
        cells[oldcsubs[i]].isplit = nsplitlocal;

      nsplitlocal++;
      nlocal++;
    }
  }

  // delete old integer lists

  delete csurfs_old;
  delete csplits_old;
  delete csubs_old;

  // repoint particles in all grid cells to new icell index
  // particles are still sorted and have not yet been compressed
  // so count/first values in cinfo are still valid

  Particle::OnePart *particles = particle->particles;
  int *next = particle->next;

  int ip;
  for (int icell = 0; icell < nlocal; icell++) {
    ip = cinfo[icell].first;
    while (ip >= 0) {
      particles[ip].icell = icell;
      ip = next[ip];
    }
  }

  // some fixes have post-compress operations to perform

  if (modify->n_pergrid) modify->compress_grid(1);
}

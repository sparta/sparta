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
#include "particle.h"
#include "domain.h"
#include "surf.h"
#include "collide.h"
#include "modify.h"
#include "adapt_grid.h"

using namespace SPARTA_NS;

// grid cell communication

/* ----------------------------------------------------------------------
   pack single icell into buf
   only called for owned unsplit and split cells, not for sub cells
   include its cinfo, surfs, split info, sub cells if necessary
   ownflag = 1/0 = receiver stores cell as owned or ghost
     for owned cell, also pack cinfo and auxiliary collision/fix info
     if nsurf < 0, is an empty ghost cell, pack only the cell
   partflag = 0/1 = no/yes to also pack particles
   surfflag = 0/1 = no/yes to also pack surface info
     surfflag = no always called with and ownflag = 0 and partflag = 0
   memflag = 0/1 = no/yes to actually pack into buf, 0 = just length
   return length of packing in bytes
   called from Grid::acquire_ghosts(), Comm::migrate_cells()
------------------------------------------------------------------------- */

int Grid::pack_one(int icell, char *buf,
                   int ownflag, int partflag, int surfflag, int memflag)
{
  char *ptr = buf;

  // pack child cell data struct, csurf ptr will be reset when unpacked

  if (memflag) memcpy(ptr,&cells[icell],sizeof(ChildCell));
  ptr += sizeof(ChildCell);
  ptr = ROUNDUP(ptr);

  // pack any custom grid data

  if (ncustom) {
    pack_custom(icell,buf);
    ptr += nbytes_custom;
    ptr = ROUNDUP(ptr);
  }

  // no surfs or any other info
  // ditto for sending empty ghost

  if (!surfflag) return ptr - buf;
  if (cells[icell].nsurf < 0) return ptr - buf;

  // if nsurfs, pack different info for explicit vs implicit surfs
  // explicit all: just list of csurf indices
  // explicit distributed and implicit: list of lines or triangles

  int nsurf = cells[icell].nsurf;
  if (nsurf) {
    if (!surf->implicit && !surf->distributed) {
      if (memflag) memcpy(ptr,cells[icell].csurfs,nsurf*sizeof(surfint));
      ptr += nsurf*sizeof(surfint);
      ptr = ROUNDUP(ptr);
    } else if (domain->dimension == 2) {
      Surf::Line *lines = surf->lines;
      int sizesurf = sizeof(Surf::Line);
      surfint *csurfs = cells[icell].csurfs;
      for (int m = 0; m < nsurf; m++) {
        int isurf = csurfs[m];
        if (memflag) memcpy(ptr,&lines[isurf],sizesurf);
        ptr += sizesurf;
        ptr = ROUNDUP(ptr);
      }
    } else {
      Surf::Tri *tris = surf->tris;
      int sizesurf = sizeof(Surf::Tri);
      surfint *csurfs = cells[icell].csurfs;
      for (int m = 0; m < nsurf; m++) {
        int isurf = csurfs[m];
        if (memflag) memcpy(ptr,&tris[isurf],sizesurf);
        ptr += sizesurf;
        ptr = ROUNDUP(ptr);
      }
    }
  }

  if (ownflag) {
    if (memflag) memcpy(ptr,&cinfo[icell],sizeof(ChildInfo));
    ptr += sizeof(ChildInfo);
    ptr = ROUNDUP(ptr);
  }

  // if split cell, pack sinfo and sinfo.csplits and sinfo.csubs
  // if ownflag, also pack volumes from sub cells themselves

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

    if (ownflag) {
      int *csubs = sinfo[isplit].csubs;
      double *dptr = (double *) ptr;
      if (memflag)
        for (int i = 0; i < nsplit; i++) dptr[i] = cinfo[csubs[i]].volume;
      ptr += nsplit*sizeof(double);
    }
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

  if (!partflag) return ptr - buf;

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
   unpack received ghost cells into cell lists
------------------------------------------------------------------------- */

void Grid::unpack_ghosts(int nsize, char *buf, void *ptr)
{
  Grid *gptr = (Grid *) ptr;
  int surfflag = gptr->unpack_ghosts_surfflag;

  int n = 0;
  while (n < nsize)
    n += gptr->unpack_one(&buf[n],0,0,surfflag);
}

/* ----------------------------------------------------------------------
   unpack single icell from buf
   include its cinfo, surfs, split info, sub cells if necessary
   ownflag = 1/0 = receiver stores cell as owned or ghost
     for owned cell, also unpack cinfo and auxiliary collision/fix info
   partflag = 0/1 = no/yes to also unpack particles
   surfflag = 0/1 = no/yes to also unpack surface info
     surfflag = 0 always called with ownflag = 0 and partflag = 0
   return length of unpacking in bytes
------------------------------------------------------------------------- */

int Grid::unpack_one(char *buf,
                     int ownflag, int partflag, int surfflag, int sortflag)
{
  char *ptr = buf;

  // unpack child cell as owned or ghost

  int icell;
  if (ownflag) icell = nlocal;
  else icell = nlocal + nghost;
  grow_cells(1,ownflag);
  if (ownflag) nlocal++;
  else nghost++;

  memcpy(&cells[icell],ptr,sizeof(ChildCell));
  ptr += sizeof(ChildCell);
  ptr = ROUNDUP(ptr);

  if (ownflag) {
    cells[icell].proc = me;
    cells[icell].ilocal = icell;
  }

  // pack any custom grid data

  if (ncustom) {
    unpack_custom(buf,icell);
    ptr += nbytes_custom;
    ptr = ROUNDUP(ptr);
  }

  // no surfs or any other info
  // ditto for EMPTY ghost with nsurf < 0
  // reset other fields for ghost cell (csurfs, nsplit, isplit)

  if (!surfflag || cells[icell].nsurf < 0) {
    cells[icell].csurfs = NULL;
    cells[icell].nsplit = 1;
    cells[icell].isplit = -1;
    return ptr - buf;
  }

  // if nsurfs, unpack different info for explicit vs implicit
  // explicit all: list of csurf indices
  // implicit: add entire list of lines or triangles
  // explicit distributed:
  //    list of lines or triangles
  //    check hash each to see if can skip b/c already have it
  //    add new surfs to hash

  int nsurf = cells[icell].nsurf;
  if (nsurf) {
    cells[icell].csurfs = csurfs->vget();

    // explicit all surfs

    if (!surf->implicit && !surf->distributed) {
      memcpy(cells[icell].csurfs,ptr,nsurf*sizeof(surfint));
      ptr += nsurf*sizeof(surfint);
      ptr = ROUNDUP(ptr);

    // implicit surfs

    } else if (surf->implicit) {
      if (domain->dimension == 2) {
        int sizesurf = sizeof(Surf::Line);
        surfint *csurfs = cells[icell].csurfs;
        for (int m = 0; m < nsurf; m++) {
          Surf::Line *line = (Surf::Line *) ptr;
          surf->add_line_copy(ownflag,line);
          if (ownflag) csurfs[m] = surf->nlocal-1;
          else csurfs[m] = surf->nlocal+surf->nghost-1;
          ptr += sizesurf;
          ptr = ROUNDUP(ptr);
        }
      } else {
        int sizesurf = sizeof(Surf::Tri);
        surfint *csurfs = cells[icell].csurfs;
        for (int m = 0; m < nsurf; m++) {
          Surf::Tri *tri = (Surf::Tri *) ptr;
          surf->add_tri_copy(ownflag,tri);
          if (ownflag) csurfs[m] = surf->nlocal-1;
          else csurfs[m] = surf->nlocal+surf->nghost-1;
          ptr += sizesurf;
          ptr = ROUNDUP(ptr);
        }
      }

    // explicit distributed surfs

    } else {
      Surf::MySurfHash *shash = surf->hash;

      if (domain->dimension == 2) {
        int sizesurf = sizeof(Surf::Line);
        surfint *csurfs = cells[icell].csurfs;
        for (int m = 0; m < nsurf; m++) {
          Surf::Line *line = (Surf::Line *) ptr;
          if (shash->find(line->id) == shash->end()) {
            surf->add_line_copy(ownflag,line);
            if (ownflag) csurfs[m] = surf->nlocal-1;
            else csurfs[m] = surf->nlocal+surf->nghost-1;
            (*shash)[line->id] = csurfs[m];
          } else csurfs[m] = (*shash)[line->id];
          ptr += sizesurf;
          ptr = ROUNDUP(ptr);
        }
      } else {
        int sizesurf = sizeof(Surf::Tri);
        surfint *csurfs = cells[icell].csurfs;
        for (int m = 0; m < nsurf; m++) {
          Surf::Tri *tri = (Surf::Tri *) ptr;
          if (shash->find(tri->id) == shash->end()) {
            surf->add_tri_copy(ownflag,tri);
            if (ownflag) csurfs[m] = surf->nlocal-1;
            else csurfs[m] = surf->nlocal+surf->nghost-1;
            (*shash)[tri->id] = csurfs[m];
          } else csurfs[m] = (*shash)[tri->id];
          ptr += sizesurf;
          ptr = ROUNDUP(ptr);
        }
      }
    }

    csurfs->vgot(nsurf);
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
  // if ownflag, also unpack volumes from sub cells themselves

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

    double *dptr;
    if (ownflag) {
      dptr = (double *) ptr;
      ptr += nsplit*sizeof(double);
    }

    int isub;
    for (int i = 0; i < nsplit; i++) {
      if (ownflag) isub = nlocal;
      else isub = nlocal + nghost;
      add_sub_cell(icell,ownflag);
      cells[isub].ilocal = sinfo[isplit].csubs[i];
      cells[isub].nsplit = -i;
      if (ownflag) cinfo[isub].volume = dptr[i];
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

  if (!partflag) return ptr - buf;

  ptr += unpack_particles(ptr,icell,sortflag);

  // unpack particles of sub cells

  if (cells[icell].nsplit > 1) {
    int isplit = cells[icell].isplit;
    int nsplit = cells[icell].nsplit;
    for (int i = 0; i < nsplit; i++) {
      int m = sinfo[isplit].csubs[i];
      ptr += unpack_particles(ptr,m,sortflag);
    }
  }

  return ptr - buf;
}

/* ----------------------------------------------------------------------
   pack single icell into buf with grid adaptation info
   memflag = 0/1 = no/yes to actually pack into buf, 0 = just length
   flag all packed particles for deletion
   return length of packing in bytes
------------------------------------------------------------------------- */

int Grid::pack_one_adapt(char *inbuf, char *buf, int memflag)
{
  int isub;

  AdaptGrid::SendAdapt *s = (AdaptGrid::SendAdapt *) inbuf;
  char *ptr = buf;

  int owner = s->owner;
  int icell = s->icell;
  int nsurf = s->nsurf;
  int np = s->np;

  // pack SendAdapt data struct

  if (memflag) memcpy(ptr,s,sizeof(AdaptGrid::SendAdapt));
  ptr += sizeof(AdaptGrid::SendAdapt);
  ptr = ROUNDUP(ptr);

  // if parent cell owner also owns the child cell being packed, then done
  // otherwise pack surfs and particles for comm to another proc

  if (owner == me) return ptr-buf;

  // pack explicit surf info
  // for non-distributed, just pack indices, since every proc stores all surfs
  // for distributed, pack the surfs themselves

  if (nsurf) {
    int isurf;
    if (!surf->distributed) {
      if (memflag) memcpy(ptr,cells[icell].csurfs,nsurf*sizeof(surfint));
      ptr += nsurf*sizeof(surfint);
    } else if (domain->dimension == 2) {
      Surf::Line *lines = surf->lines;
      int sizesurf = sizeof(Surf::Line);
      surfint *csurfs = cells[icell].csurfs;
      for (int m = 0; m < nsurf; m++) {
        isurf = csurfs[m];
        if (memflag) memcpy(ptr,&lines[isurf],sizesurf);
        ptr += sizesurf;
      }
    } else {
      Surf::Tri *tris = surf->tris;
      int sizesurf = sizeof(Surf::Tri);
      surfint *csurfs = cells[icell].csurfs;
      for (int m = 0; m < nsurf; m++) {
        isurf = csurfs[m];
        if (memflag) memcpy(ptr,&tris[isurf],sizesurf);
        ptr += sizesurf;
      }
    }
    ptr = ROUNDUP(ptr);
  }

  // done if no particles

  if (np == 0) return ptr-buf;

  // pack particles

  if (memflag) {

    // pack particles for unsplit cell
    // flag each particle for deletion

    Particle::OnePart *particles = particle->particles;
    int *next = particle->next;

    if (cells[icell].nsplit == 1) {
      int ip = cinfo[icell].first;
      while (ip >= 0) {
        memcpy(ptr,&particles[ip],nbytes_particle);
        ptr += nbytes_particle;
        if (ncustom_particle) {
          particle->pack_custom(ip,ptr);
          ptr += nbytes_particle_custom;
        }
        particles[ip].icell = -1;
        ip = next[ip];
      }

    // pack particles for all sub cells
    // flag each particle for deletion

    } else {
      int isplit = cells[icell].isplit;
      int nsplit = cells[icell].nsplit;
      for (int i = 0; i < nsplit; i++) {
        isub = sinfo[isplit].csubs[i];

        int ip = cinfo[isub].first;
        while (ip >= 0) {
          memcpy(ptr,&particles[ip],nbytes_particle);
          ptr += nbytes_particle;
          if (ncustom_particle) {
            particle->pack_custom(ip,ptr);
            ptr += nbytes_particle_custom;
          }
          particles[ip].icell = -1;
          ip = next[ip];
        }
      }
    }

  } else ptr += np * nbytes_particle_total;

  ptr = ROUNDUP(ptr);
  return ptr-buf;
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
      if (ncustom_particle) {
        particle->pack_custom(ip,ptr);
        ptr += nbytes_particle_custom;
      }
      particles[ip].icell = -1;
      ip = next[ip];
    }
  } else ptr += np * nbytes_particle_total;

  ptr = ROUNDUP(ptr);
  return ptr-buf;
}

/* ----------------------------------------------------------------------
   unpack particles from buf into one cell = icell
   point unpacked particles to icell
   if sortflag: keep particles sorted via setting
     particle->cinfo.first/count and particle->next values for added particles
   return length of unpacking in bytes
------------------------------------------------------------------------- */

int Grid::unpack_particles(char *buf, int icell, int sortflag)
{
  char *ptr = buf;

  int np = *((int *) ptr);
  ptr += sizeof(int);
  ptr = ROUNDUP(ptr);
  if (!np) return ptr-buf;

  particle->grow(np);

  Particle::OnePart *particles = particle->particles;
  int nplocal = particle->nlocal;

  if (ncustom_particle) {
    int n = nplocal;
    for (int i = 0; i < np; i++) {
      memcpy(&particles[n],ptr,nbytes_particle);
      ptr += nbytes_particle;
      particle->unpack_custom(ptr,n);
      ptr += nbytes_particle_custom;
      n++;
    }
  } else {
    memcpy(&particles[nplocal],ptr,np*nbytes_particle);
    ptr += np * nbytes_particle_total;
  }

  ptr = ROUNDUP(ptr);

  int npnew = nplocal + np;
  for (int i = nplocal; i < npnew; i++) particles[i].icell = icell;
  particle->nlocal = npnew;

  if (sortflag) {
    particle->grow_next();

    cinfo[icell].first = nplocal;
    cinfo[icell].count = np;

    for (int i = nplocal; i < npnew-1; i++)
      particle->next[i] = i+1;
    particle->next[npnew-1] = -1;
  }

  return ptr-buf;
}

/* ----------------------------------------------------------------------
   unpack Np particles from buf
   called from AdaptGrid for one cell's particles sent to another proc
------------------------------------------------------------------------- */

void Grid::unpack_particles_adapt(int np, char *buf)
{
  particle->grow(np);

  Particle::OnePart *particles = particle->particles;
  int nplocal = particle->nlocal;

  if (ncustom_particle) {
    char *ptr = buf;
    for (int i = 0; i < np; i++) {
      memcpy(&particles[nplocal],ptr,nbytes_particle);
      ptr += nbytes_particle;
      particle->unpack_custom(ptr,nplocal);
      ptr += nbytes_particle_custom;
      nplocal++;
    }
  } else {
    memcpy(&particles[nplocal],buf,np*nbytes_particle);
    nplocal += np;
  }

  particle->nlocal = nplocal;

  // insure particle->next list grows to accommodate new particles
  // since particle sorting is maintained during this stage of adaptation

  particle->grow_next();
}

/* ----------------------------------------------------------------------
   commpress owned grid cells due to some cells migrating to other procs
     or due to deletion from AdaptGrid
   removed cells are marked with proc != me, could be another proc or -1
   operates on cells, cinfo, sinfo
   integer lists (csurfs,csplits,csubs) are re-created from scratch
     copy old info from kept cells into new data structs
   also resets particle cells
     particles MUST be sorted and uncompressed
     when done, particles are still sorted in new compressed grid cells
------------------------------------------------------------------------- */

void Grid::compress()
{
  // copy of integer lists
  // create new lists

  MyPage<surfint> *csurfs_old = csurfs;
  MyPage<int> *csplits_old = csplits;
  MyPage<int> *csubs_old = csubs;

  csurfs = NULL; csplits = NULL; csubs = NULL;
  allocate_surf_arrays();

  // compress cells and cinfo and sinfo arrays
  // 4 cases:
  // discard unsplit or sub cell:
  //   b/c assigned to another procs
  //   continue
  // keep unsplit or sub cell:
  //   memcpy of cells and cinfo, setup csurfs
  //   reset cells.ilocal and sinfo.csubs
  //   increment nlocal, nunsplitlocal, nsublocal
  // discard split cell:
  //   b/c assigned to another proc
  //   flag all its sub cells with other proc (come later in list)
  //   continue
  // keep split cell:
  //   memcpy of cells and cinfo and sinfo, setup csurfs/csplits/csubs
  //   reset sinfo.icell
  //   reset cells.isplit, ditto for sub cells
  //   increment nlocal, nsplitlocal
  // relies on sub cells appearing in cells array after their split cell
  //   no need for sub cells of one split cell to be contiguous before or after

  int ncurrent = nlocal;
  nlocal = nunsplitlocal = nsplitlocal = nsublocal = 0;

  for (int icell = 0; icell < ncurrent; icell++) {

    // unsplit and sub cells

    if (cells[icell].nsplit <= 1) {
      if (cells[icell].proc != me) continue;

      // copy cell from nlocal to icell
      // collide and fixes also need to do the same

      if (icell != nlocal) {
        memcpy(&cells[nlocal],&cells[icell],sizeof(ChildCell));
        memcpy(&cinfo[nlocal],&cinfo[icell],sizeof(ChildInfo));
        if (collide) collide->copy_grid_one(icell,nlocal);
        if (modify->n_pergrid) modify->copy_grid_one(icell,nlocal);
      }

      cells[nlocal].ilocal = nlocal;

      // for unsplit cell, new copy of all csurfs indices
      // for sub cell, just copy csurfs ptr from its split cell

      if (cells[nlocal].nsurf) {
        if (cells[nlocal].nsplit == 1) {
          surfint *oldcsurfs = cells[icell].csurfs;
          cells[nlocal].csurfs = csurfs->vget();
          memcpy(cells[nlocal].csurfs,oldcsurfs,
                 cells[nlocal].nsurf*sizeof(surfint));
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

      // discard and mark sub cells (appear later in list) with other proc

      if (cells[icell].proc != me) {
        int isplit = cells[icell].isplit;
        int nsplit = cells[icell].nsplit;
        for (int i = 0; i < nsplit; i++) {
          int m = sinfo[isplit].csubs[i];
          cells[m].proc = cells[icell].proc;
        }
        continue;
      }

      // copy cell from nlocal to icell
      // collide and fixes also need to do the same

      if (icell != nlocal) {
        memcpy(&cells[nlocal],&cells[icell],sizeof(ChildCell));
        memcpy(&cinfo[nlocal],&cinfo[icell],sizeof(ChildInfo));
        if (collide) collide->copy_grid_one(icell,nlocal);
        if (modify->n_pergrid) modify->copy_grid_one(icell,nlocal);
      }

      cells[nlocal].ilocal = nlocal;

      // new copy of all csurfs indices

      surfint *oldcsurfs = cells[icell].csurfs;
      cells[nlocal].csurfs = csurfs->vget();
      memcpy(cells[nlocal].csurfs,oldcsurfs,
             cells[nlocal].nsurf*sizeof(surfint));
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

  // reset final grid cell count in collide and fixes

  if (collide) collide->reset_grid_count(nlocal);
  if (modify->n_pergrid) modify->reset_grid_count(nlocal);

  // delete old integer lists

  delete csurfs_old;
  delete csplits_old;
  delete csubs_old;

  // repoint particles in all remaining grid cells to new icell indices
  // assumes particles are sorted and have not yet been compressed,
  //   so count/first values in compressed cinfo data struct are still valid
  // when done, particles are still sorted

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
}

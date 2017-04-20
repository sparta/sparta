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

#include "grid.h"
#include "particle.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "adapt_grid.h"
#include "cut3d.h"
#include "cut2d.h"
#include "error.h"

using namespace SPARTA_NS;

enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};   // several files

/* ----------------------------------------------------------------------
   refine child icell into Nx by Ny by Nz new child cells
   caller already created new pcells[iparent] from child icell
   old child icell still exists, will be deleted when grid compressed later
------------------------------------------------------------------------- */

void Grid::refine_cell(int icell, int iparent,
                       int nx, int ny, int nz, int *childlist,
                       Cut2d *cut2d, Cut3d *cut3d)
{
  int i,j,m,ix,iy,iz,offset,ichild;
  int dim,ncorner,mark,ip,ipnew;
  cellint id;
  double lo[3],hi[3];

  dim = domain->dimension;
  ncorner = 8;
  if (dim == 2) ncorner = 4;

  // remove icell from grid hash

  hash->erase(cells[icell].id);

  // loop over creation of new child cells

  m = 0;
  for (iz = 0; iz < nz; iz++)
    for (iy = 0; iy < ny; iy++)
      for (ix = 0; ix < nx; ix++) {
        childlist[m] = nlocal;
        m++;
        id = pcells[iparent].id | ((cellint) m << pcells[iparent].nbits);
        id_child_lohi(iparent,m,lo,hi);
        add_child_cell(id,iparent,lo,hi);
        weight_one(nlocal-1);

        // add each new child cell to grid hash
        
        (*hash)[id] = nlocal;

        // update any per grid fixes for the new child cell

        if (modify->n_pergrid) modify->add_grid_one(nlocal-1,0);

        // if surfs in parent cell, intersect them with child cell
        // add_child_cell marked child type as OUTSIDE
        //   correct only if no surfs in parent
        //   if surfs in parent, mark all child cells as UNKNOWN
        // if surfs in child, surf2grid_one will mark type = OVERLAP
        //   and will mark corner points of new child cell

        if (cells[icell].nsurf) {
          ichild = nlocal - 1;
          cinfo[ichild].type = UNKNOWN;
          surf2grid_one(0,ichild,icell,-1,cut3d,cut2d);
          
          // update any per grid fixes for the newly created sub cells
          
          if (modify->n_pergrid)
            for (i = ichild+1; i < nlocal; i++)
              modify->add_grid_one(i,0);
        }
      }

  // use parent corner pts to mark type of new child cells with no surfs
  //   as INSIDE/OUTSIDE, no need to mark corner pts of an INSIDE/OUTSIDE cell
  // child cells with surfs were marked above, including their corner pts
  // if mark as OUTSIDE, no need to set volume since add_cell() set it

  if (cells[icell].nsurf) {
    for (i = 0; i < ncorner; i++) {
      if (i == 0) offset = 0;
      else if (i == 1) offset = nx-1;
      else if (i == 2) offset = ny*nx - nx;
      else if (i == 3) offset = nx*ny - 1;
      else if (i == 4) offset = (nz-1)*ny*nx;
      else if (i == 5) offset = (nz-1)*ny*nx + nx-1;
      else if (i == 6) offset = (nz-1)*ny*nx + ny*nx - nx;
      else if (i == 7) offset = (nz-1)*ny*nx + nx*ny - 1;
      m = childlist[offset];

      if (cinfo[m].type == OVERLAP) continue;

      // parent corner pt could be UNKNOWN if iterating
      //   and have not yet re-marked all new grid cells via set_inout()
      // otherwise mark cell type as INSIDE or OUTSIDE
      //   to match parent corner, since new child has no surfs

      mark = cinfo[icell].corner[i];
      if (mark == UNKNOWN) continue;
      cinfo[m].type = mark;
    }
  }

  // if old child was a split cell,
  // add all its sub cell particles to old split cell
  
  if (cells[icell].nsplit > 1)
    combine_split_cell_particles(icell,0);
  
  // assign particles in single old child cell to many new child cells
  // build linked list of particles within each new child cell
  //   via particle next and cinfo count and first
  //   cinfo->mask enables this by storing last particle in cell (so far)
  //   then reset cinfo->mask = 1 (all) for newly created child cells
  // must allow for each new child cell to be a split cell
  // NOTE: this is tricky code, document it better
  //       is inverse of coarsen_cell() logic
  // NOTE: may not need to build new linkeds list if will re-sort anyway
  
  Particle::OnePart *particles = particle->particles;
  int *next = particle->next;
  Particle::OnePart *p;

  ip = cinfo[icell].first;
  while (ip >= 0) {
    p = &particles[ip];
    
    // ichild = new unsplit or split cell the particle is now in
    
    ichild = id_find_child(iparent,p->x);
    if (ichild < 0) {
      printf("BAD CHILD %d: %d " CELLINT_FORMAT " %d " CELLINT_FORMAT ": %d\n",
             me,icell,cells[icell].id,iparent,pcells[iparent].id,ichild);
      error->one(FLERR,"Adapt particle remap indexed bad child cell");
    }
    
    // if ichild is split cell:
    // use split2d/3d to find which sub cell particle is in
    // reset ichild = index of sub cell
    
    if (cells[ichild].nsplit > 1) {
      if (dim == 3) ichild = update->split3d(ichild,p->x);
      else ichild = update->split2d(ichild,p->x);
    }
    
    // re-assign particle to ichild
    
    particles[ip].icell = ichild;
    cinfo[ichild].count++;
    if (cinfo[ichild].first < 0) cinfo[ichild].first = ip;
    else next[cinfo[ichild].mask] = ip;
    cinfo[ichild].mask = ip;
    ipnew = next[ip];
    next[ip] = -1;
    ip = ipnew;
  }
  
  // erase particles in old child cell that has become a parent
  
  cinfo[icell].count = 0;
  cinfo[icell].first = -1;
}

/* ----------------------------------------------------------------------
   coarsen parent cell iparent with its nchild child cells
   create a new child cell
------------------------------------------------------------------------- */

void Grid::coarsen_cell(int iparent, int nchild, 
                        int *proc, int *index, int *recv,
                        AdaptGrid *ag,
                        Cut2d *cut2d, Cut3d *cut3d)
{
  int i,j,k,m,icell,ichild,inew,nsurf,ip,ipnew,ns;
  int *ptr,*cs;

  AdaptGrid::SendAdapt **sa_header = ag->sa_header;
  int **sa_csurfs = ag->sa_csurfs;
  char **sa_particles = ag->sa_particles;

  int dim = domain->dimension;
  int maxsurfpercell = grid->maxsurfpercell;

  // add parent as new child cell at end of my cells
  // add new child cell to grid hash
  
  add_child_cell(pcells[iparent].id,pcells[iparent].iparent,
                 pcells[iparent].lo,pcells[iparent].hi);
  weight_one(nlocal-1);
        
  (*hash)[pcells[iparent].id] = nlocal;
  inew = nlocal - 1;

  // update any per grid fixes for the new child cell
  // add new child cell to newcells for later processing
  
  if (modify->n_pergrid) modify->add_grid_one(nlocal-1,0);

  // if any children have surfs, add union to new child cell
  // NOTE: any way to avoid N^2 loop for unique surfs
  
  nsurf = 0;
  ptr = csurfs->vget();

  for (m = 0; m < nchild; m++) {
    if (proc[m] == me) {
      icell = index[m];
      if (cells[icell].nsurf == 0) continue;
      ns = cells[icell].nsurf;
      cs = cells[icell].csurfs;
    } else {
      ns = sa_header[recv[m]]->nsurf;
      cs = sa_csurfs[recv[m]];
    }
    
    for (j = 0; j < ns; j++) {
      for (k = 0; k < nsurf; k++)
        if (ptr[k] == cs[j]) break;
      if (k < nsurf) continue;
      if (nsurf == maxsurfpercell)
        error->one(FLERR,"Too many surfs in coarsened cell");
      ptr[nsurf++] = cs[j];
    }
  }

  // new child cell has surfs
  // surf2grid_one() cuts and splits the new cell
  // corners of new child cell were set by add_cell() to UNKNOWN
  // if has surfs, type = OVERLAP and surf2grid_one() will set corners

  if (nsurf) {
    cinfo[inew].type = OVERLAP;
    cells[inew].nsurf = nsurf;
    cells[inew].csurfs = ptr;
    csurfs->vgot(nsurf);
    
    surf2grid_one(1,inew,-1,nsurf,cut3d,cut2d);
    //cells = cells;
    //cinfo = cinfo;
    //sinfo = sinfo;
    
    // update any per grid fixes for newly created sub cells
      
    if (modify->n_pergrid)
      for (i = inew+1; i < nlocal; i++)
        modify->add_grid_one(i,0);

  // if new cell has no surfs,
  // set type = INSIDE if one of its children was INSIDE

  } else {
    if (proc[0] == me) {
      if (cinfo[index[0]].type == INSIDE) cinfo[inew].type = INSIDE;
    } else {
      if (sa_header[recv[0]]->type == INSIDE) cinfo[inew].type = INSIDE;
    }
  }

  // if any of old children I own is a split cell,
  //   add all its sub cell particles to split cell parent
  // not needed for child cells received from other procs
  
  for (m = 0; m < nchild; m++) {
    if (proc[m] != me) continue;
    ichild = index[m];
    if (cells[ichild].nsplit == 1) continue;
    combine_split_cell_particles(ichild,0);
  }

  // assign particles in many old child cells to single new child cell
  // each of children is owned by me or part of recvlist
  // build linked list of particles within new child cell
  //   via particle next and cinfo count and first
  //   cinfo->mask enables this by storing last particle in cell (so far)
  //   then reset cinfo->mask = 1 (all) for newly created child cells
  // must allow for new child cell to be a split cell
  // NOTE: this is tricky code, document it better
  //       is inverse of refine_cell() logic
  // NOTE: may not need to build new linked list if will re-sort anyway

  Particle::OnePart *particles = particle->particles;
  int *next = particle->next;
  Particle::OnePart *p;

  for (m = 0; m < nchild; m++) {
    
    // old child cell that I own
    
    if (proc[m] == me) {
      ichild = index[m];
      ip = cinfo[ichild].first;
      
      while (ip >= 0) {
        p = &particles[ip];
        
        // if new child is not split: assign particle to icell
        // else:
        //   use split2d/3d to find which sub cell particle is in
        //   assign particle to sub cell
        
        if (cells[inew].nsplit == 1) icell = inew;
        else {
          if (dim == 3) icell = update->split3d(inew,p->x);
          else icell = update->split2d(inew,p->x);
        }
        
        // re-assign particle to icell
        
        particles[ip].icell = icell;
        cinfo[icell].count++;
        if (cinfo[icell].first < 0) cinfo[icell].first = ip;
        else next[cinfo[icell].mask] = ip;
        cinfo[icell].mask = ip;
        ipnew = next[ip];
        next[ip] = -1;
        ip = ipnew;
      }
      
      // erase particles in old child cell that is being deleted
      
      cinfo[ichild].count = 0;
      cinfo[ichild].first = -1;
      
    // old child cell sent from another proc
      
    } else {
      int np = sa_header[recv[m]]->np;
      unpack_particles_adapt(np,sa_particles[recv[m]]);
      particles = particle->particles;
      next = particle->next;
      int nplocal = particle->nlocal;
      
      for (ip = nplocal-np; ip < nplocal; ip++) {
        
        // if new child is not split: assign particle to icell
        // else:
        //   use split2d/3d to find which sub cell particle is in
        //   assign particle to sub cell
        
        if (cells[inew].nsplit == 1) icell = inew;
        else {
          if (dim == 3) icell = update->split3d(inew,particles[ip].x);
          else icell = update->split2d(inew,particles[ip].x);
        }
        
        // assign particle to icell
        
        particles[ip].icell = icell;
        cinfo[icell].count++;
        if (cinfo[icell].first < 0) cinfo[icell].first = ip;
        else next[cinfo[icell].mask] = ip;
        cinfo[icell].mask = ip;
        next[ip] = -1;
      }
    }
  }
}

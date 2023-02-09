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
   icell still exists, will be deleted when grid is compressed later
   childlist is workspace of length Nx by Ny by Nz, passed by caller
------------------------------------------------------------------------- */

void Grid::refine_cell(int icell, int *childlist, Cut2d *cut2d, Cut3d *cut3d)
{
  int i,m,ix,iy,iz,offset,ichild;
  int dim,ncorner,mark,ip,ipnew;
  cellint childID;
  double lo[3],hi[3];
  double *plo,*phi;

  dim = domain->dimension;
  ncorner = 8;
  if (dim == 2) ncorner = 4;

  cellint parentID = cells[icell].id;
  int plevel = cells[icell].level;

  int pbits = plevels[plevel].nbits;
  int nx = plevels[plevel].nx;
  int ny = plevels[plevel].ny;
  int nz = plevels[plevel].nz;

  // loop over creation of new child cells
  // add new children to hash, so id_find_child() for particles will work
  // set plo/phi inside loop b/c cells can be realloced by add_child_cell()

  m = 0;
  for (iz = 0; iz < nz; iz++)
    for (iy = 0; iy < ny; iy++)
      for (ix = 0; ix < nx; ix++) {
        childlist[m] = nlocal;
        m++;
        childID = ((cellint) m << pbits) | parentID;
        plo = cells[icell].lo;
        phi = cells[icell].hi;
        id_child_lohi(plevel,plo,phi,m,lo,hi);
        add_child_cell(childID,plevel+1,lo,hi);
        weight_one(nlocal-1);
        (*hash)[childID] = nlocal-1;

        // update any per grid fixes for the new child cell

        if (modify->n_pergrid) modify->add_grid_one();

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
              modify->add_grid_one();
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

  // if parent was a split cell,
  // add all its sub cell particles to parent split cell

  if (cells[icell].nsplit > 1)
    combine_split_cell_particles(icell,0);

  // assign particles in parent cell to many new child cells
  // build linked list of particles within each new child cell
  //   via particle next and cinfo count and first and mask
  //   cinfo->mask enables this by storing last particle in cell (so far)
  //   caller will reset cinfo->mask for newly created child cells
  // must allow for each new child cell to be a split cell
  // NOTE: may not need to build new linked lists if will re-sort anyway

  Particle::OnePart *particles = particle->particles;
  int *next = particle->next;
  Particle::OnePart *p;

  plo = cells[icell].lo;
  phi = cells[icell].hi;

  ip = cinfo[icell].first;
  while (ip >= 0) {
    p = &particles[ip];

    // ichild = index of new unsplit or split child cell the particle is now in

    ichild = id_find_child(parentID,plevel,plo,phi,p->x);

    if (ichild < 0) {
      printf("BAD CHILD %d: %d " CELLINT_FORMAT ": %g %g %g\n",
             me,icell,parentID,p->x[0],p->x[1],p->x[2]);
      error->one(FLERR,"Adapt refine particle could not be mapped to child cell");
    }

    // if ichild is split cell:
    //   use split2d/3d to find which sub cell particle is in
    // reset ichild = index of sub cell

    if (cells[ichild].nsplit > 1) {
      if (dim == 3) ichild = update->split3d(ichild,p->x);
      else ichild = update->split2d(ichild,p->x);
    }

    // re-assign particle to ichild

    particles[ip].icell = ichild;

    // create new linked lists for particles in each child cell
    // necessary b/c grid cells will be compressed after refinement

    cinfo[ichild].count++;
    if (cinfo[ichild].first < 0) cinfo[ichild].first = ip;
    else next[cinfo[ichild].mask] = ip;
    cinfo[ichild].mask = ip;
    ipnew = next[ip];
    next[ip] = -1;
    ip = ipnew;
  }
}

/* ----------------------------------------------------------------------
   create a new coarnsed child cell with parentID
   use info from its nchild child cells
------------------------------------------------------------------------- */

void Grid::coarsen_cell(cellint parentID, int plevel, double *plo, double *phi,
                        int nchild, int *index,
                        int *nsurf_child, int *npart_child,
                        void **surf_child, char **part_child,
                        Cut2d *cut2d, Cut3d *cut3d)
{
  int i,j,k,m,icell,jcell,ns,ip,ipnew;

  int dim = domain->dimension;
  int maxsurfpercell = grid->maxsurfpercell;

  // add parent as new child cell at end of my cells
  // add new child to hash

  add_child_cell(parentID,plevel,plo,phi);
  weight_one(nlocal-1);
  (*hash)[parentID] = nlocal-1;
  int newcell = nlocal - 1;

  // update any per grid fixes for the new child cell

  if (modify->n_pergrid) modify->add_grid_one();

  // if any children have surfs, add union of surfs to new child cell
  // shash enables union by storing surfIDs added to new cell
  // process explicit surf info in surf_child vectors

  int nsurf = 0;
  surfint *ptr = csurfs->vget();

  // for non-distributed surfs:
  // surf_child = vectors of indices into global list owned by all procs

  if (!surf->distributed) {
    MySurfHash shash;
    Surf::Line *lines = surf->lines;
    Surf::Tri *tris = surf->tris;
    surfint *cs;
    surfint surfID;
    int ns;

    for (m = 0; m < nchild; m++) {
      if (index[m] >= 0) {
        icell = index[m];
        ns = cells[icell].nsurf;
        cs = cells[icell].csurfs;
      } else {
        ns = nsurf_child[m];
        cs = (surfint *) surf_child[m];
      }

      for (i = 0; i < ns; i++) {
        if (dim == 2) surfID = lines[cs[i]].id;
        else surfID = tris[cs[i]].id;
        if (shash.find(surfID) == shash.end()) {
          shash[surfID] = 0;
          if (nsurf == maxsurfpercell)
            error->one(FLERR,"Too many surfs in coarsened cell");
          ptr[nsurf++] = cs[i];
        }
      }
    }

  // for distributed surfs:
  // if this proc owns child:
  //   surf_child = vector of indices into nlocal list owned by this proc
  // if this proc does not own child:
  //   surf_child = list of lines or tris
  //   add surf to this proc's nlocal list if not already owned

  } else {
    MySurfHash shash;
    Surf::Line *lines;
    Surf::Tri *tris;
    Surf::MySurfHash *surfhash = surf->hash;
    surfint *cs;
    surfint surfID;
    int ilocal;
    int ns;

    for (m = 0; m < nchild; m++) {

      // locally owned cell, so indices

      if (index[m] >= 0) {
        icell = index[m];
        ns = cells[icell].nsurf;
        cs = cells[icell].csurfs;
        if (dim == 2) lines = surf->lines;
        else tris = surf->tris;

        for (i = 0; i < ns; i++) {
          if (dim == 2) surfID = lines[cs[i]].id;
          else surfID = tris[cs[i]].id;
          if (shash.find(surfID) == shash.end()) {
            shash[surfID] = 0;
            if (nsurf == maxsurfpercell)
              error->one(FLERR,"Too many surfs in coarsened cell");
            ptr[nsurf++] = cs[i];
          }
        }

      // communicated cell, so lines or tris

      } else {

        ns = nsurf_child[m];
        if (dim == 2) lines = (Surf::Line *) surf_child[m];
        else tris = (Surf::Tri *) surf_child[m];

        for (i = 0; i < ns; i++) {
          if (dim == 2) surfID = lines[i].id;
          else surfID = tris[i].id;
          if (shash.find(surfID) == shash.end()) {
            shash[surfID] = 0;
            if (surfhash->find(surfID) != surfhash->end())
              ilocal = (*surfhash)[surfID];
            else {
              if (dim == 2) surf->add_line_copy(1,&lines[i]);
              else surf->add_tri_copy(1,&tris[i]);
              ilocal = surf->nlocal-1;
            }
            if (nsurf == maxsurfpercell)
              error->one(FLERR,"Too many surfs in coarsened cell");
            ptr[nsurf++] = ilocal;
          }
        }
      }
    }
  }

  // new child cell has surfs
  // surf2grid_one() cuts and splits the new cell
  // corners of new child cell were set by add_cell() to UNKNOWN
  // if has surfs, type = OVERLAP and surf2grid_one() will set corners

  if (nsurf) {
    cinfo[newcell].type = OVERLAP;
    cells[newcell].nsurf = nsurf;
    cells[newcell].csurfs = ptr;
    csurfs->vgot(nsurf);

    surf2grid_one(1,newcell,-1,nsurf,cut3d,cut2d);

    // update any per grid fixes for newly created sub cells

    if (modify->n_pergrid)
      for (i = newcell+1; i < nlocal; i++)
        modify->add_grid_one();
  }

  // assign particles in old child cells to single new child cell
  // new child cell may be a split cell

  Particle::OnePart *particles = particle->particles;
  int *next = particle->next;
  Particle::OnePart *p;

  for (m = 0; m < nchild; m++) {

    // old child cell that this proc owns
    // if it is a split cell:
    //   add all its sub cell particles to split cell

    if (index[m] >= 0) {
      jcell = index[m];
      if (cells[jcell].nsplit > 1) combine_split_cell_particles(jcell,0);
      ip = cinfo[jcell].first;

      while (ip >= 0) {
        p = &particles[ip];

        // if new child is not split: assign particle to icell
        // else:
        //   use split2d/3d to find which sub cell particle is in
        //   assign particle to sub cell

        if (cells[newcell].nsplit == 1) icell = newcell;
        else {
          if (dim == 3) icell = update->split3d(newcell,p->x);
          else icell = update->split2d(newcell,p->x);
        }

        // re-assign particle to icell

        particles[ip].icell = icell;

        // create new linked list for all particles in icell (newcell)
        // necessary b/c grid cells will be compressed after coarsening

        cinfo[icell].count++;
        if (cinfo[icell].first < 0) cinfo[icell].first = ip;
        else next[cinfo[icell].mask] = ip;
        cinfo[icell].mask = ip;
        ipnew = next[ip];
        next[ip] = -1;
        ip = ipnew;
      }

    // old child cell that another proc communicated to this proc

    } else {
      int np = npart_child[m];
      unpack_particles_adapt(np,part_child[m]);
      particles = particle->particles;
      next = particle->next;
      int nplocal = particle->nlocal;

      for (ip = nplocal-np; ip < nplocal; ip++) {

        // if new child is not split: assign particle to icell
        // else:
        //   use split2d/3d to find which sub cell particle is in
        //   assign particle to sub cell

        if (cells[newcell].nsplit == 1) icell = newcell;
        else {
          if (dim == 3) icell = update->split3d(newcell,particles[ip].x);
          else icell = update->split2d(newcell,particles[ip].x);
        }

        // assign particle to icell

        particles[ip].icell = icell;

        // create new linked list for all particles in newcell
        // necessary b/c grid cells will be compressed after coarsening

        cinfo[icell].count++;
        if (cinfo[icell].first < 0) cinfo[icell].first = ip;
        else next[cinfo[icell].mask] = ip;
        cinfo[icell].mask = ip;
        next[ip] = -1;
      }
    }
  }
}

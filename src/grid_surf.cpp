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
#include "domain.h"
#include "update.h"
#include "comm.h"
#include "particle.h"
#include "surf.h"
#include "cut2d.h"
#include "cut3d.h"
#include "irregular.h"
#include "math_const.h"
#include "hashlittle.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

#define BIG 1.0e20
#define MAXSPLITPERCELL 10

enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};   // several files

// operations for surfaces in grid cells

/* ----------------------------------------------------------------------
   allocate page data structs to hold variable-length surf and cell info
------------------------------------------------------------------------- */

void Grid::allocate_surf_arrays()
{
  delete csurfs;
  delete csplits;
  delete csubs;

  csurfs = new MyPage<int>(maxsurfpercell,MAX(100*maxsurfpercell,1024));
  csplits = new MyPage<int>(maxsurfpercell,MAX(100*maxsurfpercell,1024));
  csubs = new MyPage<int>(MAXSPLITPERCELL,128);
}

void Grid::allocate_cell_arrays()
{
  delete cpsurf;
  cpsurf = new MyPage<int>(maxcellpersurf,MAX(100*maxcellpersurf,1024));
}

/* ----------------------------------------------------------------------
   request N-length vector from csubs
   called by ReadRestart
------------------------------------------------------------------------- */

int *Grid::csubs_request(int n)
{
  int *ptr = csubs->vget();
  cpsurf->vgot(n);
  return ptr;
}

/* ----------------------------------------------------------------------
   enumerate intersections of isurf with any child grid cell
   do this by recursively walking down grid tree
------------------------------------------------------------------------- */

int Grid::find_overlaps(int isurf, cellint *list)
{
  int ncell = 0;

  if (domain->dimension == 2) {
    Surf::Line *line = &surf->lines[isurf];

    double sbox[2][2];
    sbox[0][0] = MIN(line->p1[0],line->p2[0]);
    sbox[0][1] = MIN(line->p1[1],line->p2[1]);
    sbox[1][0] = MAX(line->p1[0],line->p2[0]);
    sbox[1][1] = MAX(line->p1[1],line->p2[1]);
    
    recurse2d(isurf,sbox,0,ncell,list);

  } else {
    Surf::Tri *tri = &surf->tris[isurf];

    double sbox[2][3];
    sbox[0][0] = MIN(tri->p1[0],tri->p2[0]);
    sbox[0][0] = MIN(tri->p3[0],sbox[0][0]);
    sbox[0][1] = MIN(tri->p1[1],tri->p2[1]);
    sbox[0][1] = MIN(tri->p3[1],sbox[0][1]);
    sbox[0][2] = MIN(tri->p1[2],tri->p2[2]);
    sbox[0][2] = MIN(tri->p3[2],sbox[0][2]);

    sbox[1][0] = MAX(tri->p1[0],tri->p2[0]);
    sbox[1][0] = MAX(tri->p3[0],sbox[1][0]);
    sbox[1][1] = MAX(tri->p1[1],tri->p2[1]); 
    sbox[1][1] = MAX(tri->p3[1],sbox[1][1]);
    sbox[1][2] = MAX(tri->p1[2],tri->p2[2]);
    sbox[1][2] = MAX(tri->p3[2],sbox[1][2]);
    
    recurse3d(isurf,sbox,0,ncell,list);
  }

  return ncell;
}

/* ----------------------------------------------------------------------
   enumerate intersections of isurf with any child grid cell in iparent cell
------------------------------------------------------------------------- */

void Grid::recurse2d(int iline, double sbox[][2], int iparent, 
                     int &n, int *list)
{
  int ix,iy,ichild,newparent,index,parentflag,overlap;
  cellint idchild;
  double boundary;
  double newsbox[2][2];
  double clo[3],chi[3];

  double *p1 = surf->lines[iline].p1;
  double *p2 = surf->lines[iline].p2;

  ParentCell *p = &pcells[iparent];
  double *plo = p->lo;
  double *phi = p->hi;
  int nx = p->nx;
  int ny = p->ny;

  // ilo,ihi jlo,jhi = indices for range of grid cells overlapped by sbox
  // overlap means point is inside grid cell or touches boundary
  // same equation as in Grid::id_find_child()

  int ilo = static_cast<int> ((sbox[0][0]-plo[0]) * nx/(phi[0]-plo[0]));
  int ihi = static_cast<int> ((sbox[1][0]-plo[0]) * nx/(phi[0]-plo[0]));
  int jlo = static_cast<int> ((sbox[0][1]-plo[1]) * ny/(phi[1]-plo[1]));
  int jhi = static_cast<int> ((sbox[1][1]-plo[1]) * ny/(phi[1]-plo[1]));
  
  // decrement indices by 1 if touching lower boundary
  // should be no need to increment if touching upper boundary
  //   since ihi,jhi will then already index cell beyond boundary
  // same equation as in Grid::id_child_lohi()

  boundary = plo[0] + ilo*(phi[0]-plo[0])/nx;
  if (boundary == sbox[0][0]) ilo--;
  boundary = plo[1] + jlo*(phi[1]-plo[1])/ny;
  if (boundary == sbox[0][1]) jlo--;

  // insure each index is between 0 and N-1 inclusive

  ilo = MAX(ilo,0);
  ilo = MIN(ilo,nx-1);
  ihi = MAX(ihi,0);
  ihi = MIN(ihi,nx-1);
  jlo = MAX(jlo,0);
  jlo = MIN(jlo,ny-1);
  jhi = MAX(jhi,0);
  jhi = MIN(jhi,ny-1);

  // loop over range of grid cells between ij lo/hi inclusive
  // if cell is a child cell, compute overlap via surf2gridne()
  // else it is a parent cell, so recurse
  // set newsbox to intersection of sbox with new parent cell

  for (iy = jlo; iy <= jhi; iy++) {
    for (ix = ilo; ix <= ihi; ix++) {
      ichild = iy*nx + ix + 1;
      idchild = p->id | (ichild << p->nbits);
      grid->id_child_lohi(iparent,ichild,clo,chi);

      if (hash->find(idchild) == hash->end()) parentflag = 0;
      else if ((*hash)[idchild] >= 0) parentflag = 0;
      else parentflag = 1;
      
      if (parentflag) {
        index = (*hash)[idchild];
        newparent = -index-1;
        newsbox[0][0] = MAX(sbox[0][0],clo[0]);
        newsbox[0][1] = MAX(sbox[0][1],clo[1]);
        newsbox[1][0] = MIN(sbox[1][0],chi[0]);
        newsbox[1][1] = MIN(sbox[1][1],chi[1]);
        recurse2d(iline,newsbox,newparent,n,list);
      } else { 
        overlap = cut2d->surf2grid_one(p1,p2,clo,chi);
        if (!overlap) continue;
        if (n < maxcellpersurf) list[n] = idchild;
        n++;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   enumerate intersections of isurf with any child grid cell in iparent cell
------------------------------------------------------------------------- */

void Grid::recurse3d(int itri, double sbox[][3], int iparent, 
                     int &n, int *list)
{
  int ix,iy,iz,ichild,newparent,index,parentflag,overlap;
  cellint idchild;
  double boundary;
  double newsbox[2][3];
  double clo[3],chi[3];

  double *p1 = surf->tris[itri].p1;
  double *p2 = surf->tris[itri].p2;
  double *p3 = surf->tris[itri].p3;

  ParentCell *p = &pcells[iparent];
  double *plo = p->lo;
  double *phi = p->hi;
  int nx = p->nx;
  int ny = p->ny;
  int nz = p->nz;

  // ilo,ihi jlo,jhi = indices for range of grid cells overlapped by sbox
  // overlap means point is inside grid cell or touches boundary
  // same equation as in Grid::id_find_child()

  int ilo = static_cast<int> ((sbox[0][0]-plo[0]) * nx/(phi[0]-plo[0]));
  int ihi = static_cast<int> ((sbox[1][0]-plo[0]) * nx/(phi[0]-plo[0]));
  int jlo = static_cast<int> ((sbox[0][1]-plo[1]) * ny/(phi[1]-plo[1]));
  int jhi = static_cast<int> ((sbox[1][1]-plo[1]) * ny/(phi[1]-plo[1]));
  int klo = static_cast<int> ((sbox[0][2]-plo[2]) * nz/(phi[2]-plo[2]));
  int khi = static_cast<int> ((sbox[1][2]-plo[2]) * nz/(phi[2]-plo[2]));
  
  // decrement indices by 1 if touching lower boundary
  // should be no need to increment if touching upper boundary
  //   since ihi,jhi will then already index cell beyond boundary
  // same equation as in Grid::id_child_lohi()

  boundary = plo[0] + ilo*(phi[0]-plo[0])/nx;
  if (boundary == sbox[0][0]) ilo--;
  boundary = plo[1] + jlo*(phi[1]-plo[1])/ny;
  if (boundary == sbox[0][1]) jlo--;
  boundary = plo[2] + klo*(phi[2]-plo[2])/nz;
  if (boundary == sbox[0][2]) klo--;

  // insure each index is between 0 and N-1 inclusive

  ilo = MAX(ilo,0);
  ilo = MIN(ilo,nx-1);
  ihi = MAX(ihi,0);
  ihi = MIN(ihi,nx-1);
  jlo = MAX(jlo,0);
  jlo = MIN(jlo,ny-1);
  jhi = MAX(jhi,0);
  jhi = MIN(jhi,ny-1);
  klo = MAX(klo,0);
  klo = MIN(klo,nz-1);
  khi = MAX(khi,0);
  khi = MIN(khi,nz-1);

  // loop over range of grid cells between ij lo/hi inclusive
  // if cell is a child cell, compute overlap via surf2grid_one()
  // else it is a parent cell, so recurse
  // set newsbox to intersection of sbox with new parent cell

  for (iz = klo; iz <= khi; iz++) {
    for (iy = jlo; iy <= jhi; iy++) {
      for (ix = ilo; ix <= ihi; ix++) {
        ichild = iz*nx*ny + iy*nx + ix + 1;
        idchild = p->id | (ichild << p->nbits);
        grid->id_child_lohi(iparent,ichild,clo,chi);

        if (hash->find(idchild) == hash->end()) parentflag = 0;
        else if ((*hash)[idchild] >= 0) parentflag = 0;
        else parentflag = 1;
      
        if (parentflag) {
          index = (*hash)[idchild];
          newparent = -index-1;
          newsbox[0][0] = MAX(sbox[0][0],clo[0]);
          newsbox[0][1] = MAX(sbox[0][1],clo[1]);
          newsbox[0][2] = MAX(sbox[0][2],clo[2]);
          newsbox[1][0] = MIN(sbox[1][0],chi[0]);
          newsbox[1][1] = MIN(sbox[1][1],chi[1]);
          newsbox[1][2] = MIN(sbox[1][2],chi[2]);
          recurse3d(itri,newsbox,newparent,n,list);
        } else { 
          overlap = cut3d->surf2grid_one(p1,p2,p3,clo,chi);
          if (!overlap) continue;
          if (n < maxcellpersurf) list[n] = idchild;
          n++;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   NEW VERSION
   map surf elements into owned grid cells
   if subflag = 1, create new owned split and sub cells as needed
     called from ReadSurf, RemoveSurf, MoveSurf
   if subflag = 0, split/sub cells already exist
     called from ReadRestart
   in cells: set nsurf, csurfs, nsplit, isplit
   in cinfo: for cells with surfs, set type, corner, volume 
   initialize sinfo as needed
------------------------------------------------------------------------- */

void Grid::surf2grid2(int subflag, int outflag)
{
  int i,j,m,icell,isurf;

  int dim = domain->dimension;

  if (dim == 3) cut3d = new Cut3d(sparta);
  else cut2d = new Cut2d(sparta,domain->axisymmetric);

  // insure parent/child IDs are in hash

  if (!hashfilled) rehash();

  // compute overlap of surfs I own with all grid cells
  // info stored in ncell, celllist
  // ncell = # of cells that overlap with surf bbox

  double t1 = MPI_Wtime();

  int nprocs = comm->nprocs;
  int ntotal,ncell;

  if (dim == 2) ntotal = surf->nline;
  else ntotal = surf->ntri;

  int nsurf = ntotal / nprocs;
  if (me < nsurf % nprocs) nsurf++;

  int *cellcount;
  memory->create(cellcount,nsurf,"surf2grid2:cellcount");
  cellint **celllist = 
    (cellint **) memory->smalloc(nsurf*sizeof(cellint *),"surf2grid2:celllist");

  cellint *ptr;

  int badcount = 0;
  m = 0;
  for (isurf = comm->me; isurf < ntotal; isurf += nprocs) {
    ptr = cpsurf->vget();
    ncell = find_overlaps(isurf,ptr);

    if (ncell > maxcellpersurf) {
      cpsurf->vgot(0);
      cellcount[m] = ncell;
      badcount++;
      continue;
    }

    cellcount[m] = ncell;
    celllist[m] = ptr;
    cpsurf->vgot(ncell);
    m++;
  }

  int mytotal = 0;
  for (i = 0; i < nsurf; i++)
    mytotal += cellcount[i];

  int allcount;
  MPI_Allreduce(&mytotal,&allcount,1,MPI_INT,MPI_SUM,world);
  if (me == 0) printf("Count of surf/cell intersects = %d\n",allcount);

  MPI_Allreduce(&badcount,&allcount,1,MPI_INT,MPI_SUM,world);
  if (allcount) {
    if (me == 0) printf("All bad count = %d\n",allcount);
    error->all(FLERR,"Cells per surf exceeded limit");
  }

  double t2 = MPI_Wtime();
  printf("TIME first %g\n",t2-t1);

  // -----------------------------------------------------
  // rendezvous to convert list of cells per surf to list of surfs per cell
  // -----------------------------------------------------

  // ncount = # of my datums to send
  // include nlocal datums with owner of each grid cell

  int ncount = nlocal;
  for (i = 0; i < nsurf; i++)
    ncount += cellcount[i];

  int *proclist;
  memory->create(proclist,ncount,"surf2grid2:proclist");
  InRvous *inbuf = (InRvous *) memory->smalloc((bigint) ncount*sizeof(InRvous),
                                               "surf2grid:inbuf");
    
  // setup input buf to rendezvous comm
  // input datums = pairs of surfIDs and cellIDs
  // owning proc for each datum = random hash of cellID
  // one datum for each owned cell: datum = owning proc, cellID
  // one datum for each surf/cell pair: datum = cellID, surfID

  m = 0;
  for (i = 0; i < nlocal; i++) {
    proclist[m] = hashlittle(&cells[i].id,sizeof(cellint),0) % nprocs;
    inbuf[m].me = me;
    inbuf[m].cellID = cells[i].id;
    inbuf[m].surfID = 0;
    m++;
  }

  surfint surfID;

  for (i = 0; i < nsurf; i++) {
    if (dim == 2) surfID = surf->lines[me+i*nprocs].id;
    else surfID = surf->tris[me+i*nprocs].id;
    for (j = 0; j < cellcount[i]; j++) {
      proclist[m] = hashlittle(&celllist[i][j],sizeof(cellint),0) % nprocs;
      inbuf[m].me = -1;
      inbuf[m].cellID = celllist[i][j];
      inbuf[m].surfID = surfID;
      m++;
    }
  }

  // perform rendezvous operation
  // each proc owns random subset of cells
  // receives all info to form and return their surf lists

  char *buf;
  int nreturn = comm->rendezvous(ncount,proclist,(char *) inbuf,sizeof(InRvous),
                                 rendezvous_surflist,
                                 buf,sizeof(OutRvous),(void *) this);
  OutRvous *outbuf = (OutRvous *) buf;

  memory->destroy(proclist);
  memory->sfree(inbuf);

  // set cells nsurf and surfs for all my cells based on output info
  // output datums = pairs of cellIDs and surfIDs
  // process in 2 passes to first count surfs/cell and allocate csurfs

  // NOTE: is this part necessary?
  // is nsurf,csurfs already 0,NULL?

  for (icell = 0; icell < nlocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    cells[icell].nsurf = 0;
    cells[icell].csurfs = NULL;
  }

  for (m = 0; m < nreturn; m++) {
    icell = (*hash)[outbuf[m].cellID];
    cells[icell].nsurf++;
  }

  for (i = 0; i < nlocal; i++) {
    if (cells[icell].nsplit <= 0) continue;
    if (cells[i].nsurf)
      cells[icell].csurfs = csurfs->get(cells[i].nsurf);
    cells[icell].nsurf = 0;
  }

  // NOTE: assigning a surfint to a local int

  for (m = 0; m < nreturn; m++) {
    icell = (*hash)[outbuf[m].cellID];
    cells[icell].csurfs[cells[icell].nsurf++] = outbuf[m].surfID;
  }

  memory->sfree(outbuf);

  //memory->destroy(cellcount);
  //memory->sfree(celllist);
  //delete cpsurf;
  //cpsurf = NULL;

  double t3 = MPI_Wtime();
  printf("TIME second %g\n",t3-t2);

  surf2grid_half(subflag,outflag);
}

/* ----------------------------------------------------------------------
   on rendezvous proc:
   create list of surfs for each grid cell
   send (cellID,list) back to owning proc of grid cell
   do this in 2 passes to count sizes for each grid cell
   throw error if too many surfs/cell
---------------------------------------------------------------------- */

int Grid::rendezvous_surflist(int n, char *inbuf,
                              int *&proclist, char *&outbuf,
                              void *ptr)
{
  int i,m;

  Grid *gptr = (Grid *) ptr;
  Memory *memory = gptr->memory;
  Error *error = gptr->error;

  // initialize hash
  // ncount = number of atoms assigned to me
  // key = atom ID
  // value = index into Ncount-length data structure

  InRvous *in = (InRvous *) inbuf;
  MyHash *hash;
  
  int ncount = 0;
  for (i = 0; i < n; i++)
    if (in[i].me >= 0)
      (*hash)[in[i].cellID] = ncount++;

  // procowner = caller proc that owns each cell
  // surfcount = count of surfs per cell

  int *procowner,*surfcount;
  memory->create(procowner,ncount,"special:procowner");
  memory->create(surfcount,ncount,"special:surfcount");
  for (m = 0; m < ncount; m++) surfcount[m] = 0;

  for (i = 0; i < n; i++) { 
    m = hash->find(in[i].cellID)->second;
    if (in[i].me >= 0) procowner[m] = in[i].me;
    else surfcount[m]++;
  }

  // error check if any cell exceeds max surf count

  int maxsurfpercell = gptr->maxsurfpercell;
  for (m = 0; m < ncount; m++)
    if (surfcount[m] > maxsurfpercell)
      error->one(FLERR,"Too many surfs in one cell");

  // pass list of OutRvous datums back to comm->rendezvous

  int nout = 0;
  for (m = 0; m < ncount; m++) nout += surfcount[m];

  memory->create(proclist,nout,"special:proclist");
  OutRvous *out = (OutRvous *)
    memory->smalloc((bigint) nout*sizeof(OutRvous),"special:out");

  nout = 0;
  for (i = 0; i < n; i++) {
    if (in[i].me >= 0) continue;
    m = hash->find(in[i].cellID)->second;
    proclist[nout] = procowner[m];
    out[nout].cellID = in[i].cellID;
    out[nout].surfID = in[i].surfID;
    nout++;
  }

  outbuf = (char *) out;

  // clean up
  // Comm::rendezvous will delete proclist and out (outbuf)

  memory->destroy(procowner);
  memory->destroy(surfcount);

  return nout;
}

/* ----------------------------------------------------------------------
   map surf elements into owned grid cells
   if subflag = 1, create new owned split and sub cells as needed
     called from ReadSurf, RemoveSurf, MoveSurf
   if subflag = 0, split/sub cells already exist
     called from ReadRestart
   in cells: set nsurf, csurfs, nsplit, isplit
   in cinfo: for cells with surfs, set type, corner, volume 
   initialize sinfo as needed
------------------------------------------------------------------------- */

void Grid::surf2grid(int subflag, int outflag)
{
  int nsurf;
  int *ptr;
  double *lo,*hi;

  int dim = domain->dimension;

  double *slo = surf->bblo;
  double *shi = surf->bbhi;

  if (dim == 3) cut3d = new Cut3d(sparta);
  else cut2d = new Cut2d(sparta,domain->axisymmetric);

  // compute overlap of surfs with each cell I own
  // info stored in nsurf,csurfs
  // skip if nsplit <= 0 b/c split cells could exist if restarting

  double t1 = MPI_Wtime();

  for (int icell = 0; icell < nlocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;

    lo = cells[icell].lo;
    hi = cells[icell].hi;
    if (!box_overlap(lo,hi,slo,shi)) continue;

    ptr = csurfs->vget();

    if (dim == 3)
      nsurf = cut3d->surf2grid(cells[icell].id,cells[icell].lo,cells[icell].hi,
                               ptr,maxsurfpercell);
    else
      nsurf = cut2d->surf2grid(cells[icell].id,cells[icell].lo,cells[icell].hi,
                               ptr,maxsurfpercell);

    if (nsurf < 0) error->one(FLERR,"Too many surfs in one cell");
    if (nsurf) {
      cinfo[icell].type = OVERLAP;
      cells[icell].nsurf = nsurf;
      cells[icell].csurfs = ptr;
      csurfs->vgot(nsurf);
    }
  }

  double t2 = MPI_Wtime();
  printf("TIME %g\n",t2-t1);

  surf2grid_half(subflag,outflag);
}

/* ----------------------------------------------------------------------
   2nd half of surf2grid
   map surf elements into owned grid cells
   if subflag = 1, create new owned split and sub cells as needed
     called from ReadSurf, RemoveSurf, MoveSurf
   if subflag = 0, split/sub cells already exist
     called from ReadRestart
   in cells: set nsurf, csurfs, nsplit, isplit
   in cinfo: for cells with surfs, set type, corner, volume 
   initialize sinfo as needed
------------------------------------------------------------------------- */

void Grid::surf2grid_half(int subflag, int outflag)
{
  int i,isub,nsurf,nsplitone,xsub;
  int *surfmap,*ptr;
  double *lo,*hi,*vols;
  double xsplit[3];
  ChildCell *c;
  SplitInfo *s;

  int dim = domain->dimension;

  //double t2 = MPI_Wtime();
  //printf("TIME %g\n",t2-t1);

  if (outflag) surf2grid_stats();

  // compute cut volume and possible split of each grid cell by surfs
  // decrement nunsplitlocal if convert an unsplit cell to split cell
  // if nsplitone > 1, create new split cell sinfo and sub-cells
  // skip if nsplit <= 0 b/c split cells could exist if restarting

  int ncurrent = nlocal;
  for (int icell = 0; icell < ncurrent; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    if (cinfo[icell].type != OVERLAP) continue;

    surfmap = csplits->vget();
    c = &cells[icell];

    if (dim == 3)
      nsplitone = cut3d->split(c->id,c->lo,c->hi,c->nsurf,c->csurfs,
                               vols,surfmap,cinfo[icell].corner,xsub,xsplit);
    else
      nsplitone = cut2d->split(c->id,c->lo,c->hi,c->nsurf,c->csurfs,
                               vols,surfmap,cinfo[icell].corner,xsub,xsplit);

    if (nsplitone == 1) {
      cinfo[icell].volume = vols[0];
    
    } else if (subflag) {
      cells[icell].nsplit = nsplitone;
      nunsplitlocal--;
      
      cells[icell].isplit = nsplitlocal;
      add_split_cell(1);
      s = &sinfo[nsplitlocal-1];
      s->icell = icell;
      s->csplits = surfmap;
      s->xsub = xsub;
      s->xsplit[0] = xsplit[0];
      s->xsplit[1] = xsplit[1];
      if (dim == 3) s->xsplit[2] = xsplit[2];
      else s->xsplit[2] = 0.0;

      ptr = s->csubs = csubs->vget();

      for (i = 0; i < nsplitone; i++) {
        isub = nlocal;
        add_sub_cell(icell,1);
        cells[isub].nsplit = -i;
        cinfo[isub].volume = vols[i];
        ptr[i] = isub;
      }
      
      csplits->vgot(cells[icell].nsurf);
      csubs->vgot(nsplitone);

    } else {
      if (cells[icell].nsplit != nsplitone) {
        printf("BAD %d %d: %d %d\n",icell,cells[icell].id,
               nsplitone,cells[icell].nsplit);
        error->one(FLERR,
                   "Inconsistent surface to grid mapping in read_restart");
      }

      s = &sinfo[cells[icell].isplit];
      s->csplits = surfmap;
      s->xsub = xsub;
      s->xsplit[0] = xsplit[0];
      s->xsplit[1] = xsplit[1];
      if (dim == 3) s->xsplit[2] = xsplit[2];
      else s->xsplit[2] = 0.0;

      ptr = s->csubs;
      for (i = 0; i < nsplitone; i++) {
        isub = ptr[i];
        cells[isub].nsurf = cells[icell].nsurf;
        cells[isub].csurfs = cells[icell].csurfs;
        cinfo[isub].volume = vols[i];
      }
      
      csplits->vgot(cells[icell].nsurf);
    }
  }

  //double t3 = MPI_Wtime();
  //printf("TIME %g\n",t3-t2);

  // stats on pushed cells and unmarked corner points in OVERLAP cells

  if (outflag) {
    int npushmax;
    int *npushcell;
    if (dim == 3) {
      npushmax = cut3d->npushmax;
      npushcell = cut3d->npushcell;
    } else {
      npushmax = cut2d->npushmax;
      npushcell = cut2d->npushcell;
    }
    int *npushall = new int[npushmax+1];
    MPI_Allreduce(npushcell,npushall,npushmax+1,MPI_INT,MPI_SUM,world);
    if (comm->me == 0) {
      if (screen) {
        fprintf(screen,"  ");
        for (int i = 1; i <= npushmax; i++)
          fprintf(screen,"%d ",npushall[i]);
        fprintf(screen,"= number of pushed cells\n");
      }
      if (logfile) {
        fprintf(logfile,"  ");
        for (int i = 1; i <= npushmax; i++)
          fprintf(logfile,"%d ",npushall[i]);
        fprintf(logfile,"= number of pushed cells\n");
      }
    }
    delete [] npushall;

    int noverlap = 0;
    int ncorner = 0;
    for (int icell = 0; icell < nlocal; icell++) {
      if (cells[icell].nsplit <= 0) continue;
      if (cinfo[icell].type == OVERLAP) {
        noverlap++;
        if (cinfo[icell].corner[0] == UNKNOWN) ncorner++;
      }
    }
    int ncornerall,noverlapall;
    MPI_Allreduce(&ncorner,&ncornerall,1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&noverlap,&noverlapall,1,MPI_INT,MPI_SUM,world);
    if (comm->me == 0) {
      if (screen) fprintf(screen,"  %d %d = cells overlapping surfs, "
                          "overlap cells with unmarked corner pts\n",
                          noverlapall,ncornerall);
      if (logfile) fprintf(logfile,"  %d %d = cells overlapping surfs, "
                           "overlap cells with unmarked corner pts\n",
                           noverlapall,ncornerall);
    }
  }

  //double t4 = MPI_Wtime();
  //printf("TIME %g\n",t4-t3);

  // clean up

  if (dim == 3) delete cut3d;
  else delete cut2d;
}

/* ----------------------------------------------------------------------
   map surf elements into a single grid cell = icell
   flag = 0 for grid refinement, 1 for grid coarsening
   in cells: set nsurf, csurfs, nsplit, isplit
   in cinfo: set type, corner, volume
   initialize sinfo as needed
------------------------------------------------------------------------- */

void Grid::surf2grid_one(int flag, int icell, int iparent, int nsurf_caller,
                         Cut3d *cut3d, Cut2d *cut2d)
{
  int nsurf,isub,xsub,nsplitone;
  int *ptr;
  double xsplit[3];
  double *vols;

  int dim = domain->dimension;

  // identify surfs in new cell only for grid refinement

  if (flag == 0) {
    ptr = csurfs->vget();
    if (dim == 3)
      nsurf = cut3d->surf2grid_list(cells[icell].id,
                                    cells[icell].lo,cells[icell].hi,
                                    cells[iparent].nsurf,cells[iparent].csurfs,
                                    ptr,maxsurfpercell);
    else
      nsurf = cut2d->surf2grid_list(cells[icell].id,
                                    cells[icell].lo,cells[icell].hi,
                                    cells[iparent].nsurf,cells[iparent].csurfs,
                                    ptr,maxsurfpercell);
    
    if (nsurf == 0) return;
    if (nsurf < 0) error->one(FLERR,"Too many surfs in one refined cell");

    cinfo[icell].type = OVERLAP;
    cells[icell].nsurf = nsurf;
    cells[icell].csurfs = ptr;
    csurfs->vgot(nsurf);

  } else nsurf = nsurf_caller;
  
  // split check done for both refinement and coarsening
  
  int *surfmap = csplits->vget();
  ChildCell *c = &cells[icell];
  
  if (dim == 3)
    nsplitone = cut3d->split(c->id,c->lo,c->hi,c->nsurf,c->csurfs,
                             vols,surfmap,cinfo[icell].corner,
                             xsub,xsplit);
  else
    nsplitone = cut2d->split(c->id,c->lo,c->hi,c->nsurf,c->csurfs,
                             vols,surfmap,cinfo[icell].corner,
                             xsub,xsplit);
  
  if (nsplitone == 1) {
    if (cinfo[icell].corner[0] != UNKNOWN)
      cinfo[icell].volume = vols[0];
    
  } else {
    c->nsplit = nsplitone;
    nunsplitlocal--;
    
    c->isplit = grid->nsplitlocal;
    add_split_cell(1);
    SplitInfo *s = &sinfo[nsplitlocal-1];
    s->icell = icell;
    s->csplits = surfmap;
    s->xsub = xsub;
    s->xsplit[0] = xsplit[0];
    s->xsplit[1] = xsplit[1];
    if (dim == 3) s->xsplit[2] = xsplit[2];
    else s->xsplit[2] = 0.0;
    
    ptr = s->csubs = csubs->vget();
    
    for (int i = 0; i < nsplitone; i++) {
      isub = nlocal;
      add_sub_cell(icell,1);
      cells[isub].nsplit = -i;
      cinfo[isub].volume = vols[i];
      ptr[i] = isub;
    }
    
    csplits->vgot(nsurf);
    csubs->vgot(nsplitone);
  }
}

/* ----------------------------------------------------------------------
   remove all surf info from owned grid cells and reset cell volumes
   also remove sub cells by compressing grid cells list
   set cell type and corner flags to UNKNOWN or OUTSIDE
   called before reassigning surfs to grid cells
   changes cells data structure since sub cells are removed
   if particles exist, reassign them to new cells
------------------------------------------------------------------------- */

void Grid::clear_surf()
{
  int dimension = domain->dimension;
  int ncorner = 8;
  if (dimension == 2) ncorner = 4;
  double *lo,*hi;

  // if surfs no longer exist, set cell type to OUTSIDE, else UNKNOWN
  // set corner points of every cell to UNKNOWN

  int celltype = UNKNOWN;
  if (!surf->exist) celltype = OUTSIDE;

  int nlocal_prev = nlocal;

  int icell = 0;
  while (icell < nlocal) {
    if (cells[icell].nsplit <= 0) {
      memcpy(&cells[icell],&cells[nlocal-1],sizeof(ChildCell));
      memcpy(&cinfo[icell],&cinfo[nlocal-1],sizeof(ChildInfo));
      nlocal--;
    } else {
      cells[icell].ilocal = icell;
      cells[icell].nsurf = 0;
      cells[icell].csurfs = NULL;
      cells[icell].nsplit = 1;
      cells[icell].isplit = -1;
      cinfo[icell].type = celltype;
      for (int m = 0; m < ncorner; m++) cinfo[icell].corner[m] = UNKNOWN;
      lo = cells[icell].lo;
      hi = cells[icell].hi;
      if (dimension == 3) 
        cinfo[icell].volume = (hi[0]-lo[0]) * (hi[1]-lo[1]) * (hi[2]-lo[2]);
      else if (domain->axisymmetric)
        cinfo[icell].volume = MY_PI * (hi[1]*hi[1]-lo[1]*lo[1]) * (hi[0]-lo[0]);
      else
        cinfo[icell].volume = (hi[0]-lo[0]) * (hi[1]-lo[1]);
      icell++;
    }
  }

  // if particles exist and local cell count changed
  // repoint particles to new icell indices
  // assumes particles are sorted,
  //   so count/first values in compressed cinfo data struct are still valid
  // when done, particles are still sorted

  if (particle->exist && nlocal < nlocal_prev) {
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

  // reset csurfs and csplits and csubs so can refill

  csurfs->reset();
  csplits->reset();
  csubs->reset();

  // reset all cell counters
  
  nunsplitlocal = nlocal;
  nsplitlocal = nsublocal = 0;
}

/* ----------------------------------------------------------------------
   remove all surf info from owned grid cells and reset cell volumes
   do NOT remove sub cells by compressing grid cells list
   called from read_restart before reassigning surfs to grid cells
   sub cells already exist from restart file
------------------------------------------------------------------------- */

void Grid::clear_surf_restart()
{
  // reset current grid cells as if no surfs existed
  // just skip sub cells
  // set values in cells/cinfo as if no surfaces, including volume

  int dimension = domain->dimension;
  int ncorner = 8;
  if (dimension == 2) ncorner = 4;
  double *lo,*hi;

  for (int icell = 0; icell < nlocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    cinfo[icell].type = UNKNOWN;
    for (int m = 0; m < ncorner; m++) cinfo[icell].corner[m] = UNKNOWN;
    lo = cells[icell].lo;
    hi = cells[icell].hi;
    if (dimension == 3) 
      cinfo[icell].volume = (hi[0]-lo[0]) * (hi[1]-lo[1]) * (hi[2]-lo[2]);
    else if (domain->axisymmetric)
      cinfo[icell].volume = MY_PI * (hi[1]*hi[1]-lo[1]*lo[1]) * (hi[0]-lo[0]);
    else
      cinfo[icell].volume = (hi[0]-lo[0]) * (hi[1]-lo[1]);
  }
}

/* ----------------------------------------------------------------------
   combine all particles in sub cells of a split icell to be in split cell
   assumes particles are sorted, returns them sorted in icell
   if relabel = 1, also change icell value for each particle, else do not
------------------------------------------------------------------------- */

void Grid::combine_split_cell_particles(int icell, int relabel)
{
  int ip,iplast,jcell;

  int nsplit = cells[icell].nsplit;
  int *mycsubs = sinfo[cells[icell].isplit].csubs;
  int count = 0;
  int first = -1;

  int *next = particle->next;

  for (int i = 0; i < nsplit; i++) {
    jcell = mycsubs[i];
    count += cinfo[jcell].count;
    if (cinfo[jcell].first < 0) continue;
    
    if (first < 0) first = cinfo[jcell].first;
    else next[iplast] = cinfo[jcell].first;
    
    ip = cinfo[jcell].first;
    while (ip >= 0) {
      iplast = ip;
      ip = next[ip];
    }
  }

  cinfo[icell].count = count;
  cinfo[icell].first = first;

  // repoint each particle now in parent split cell to the split cell

  if (relabel) {
    Particle::OnePart *particles = particle->particles;
    ip = first;
    while (ip >= 0) {
      particles[ip].icell = icell;
      ip = next[ip];
    }
  }
}

/* ----------------------------------------------------------------------
   assign all particles in split icell to appropriate sub cells
   assumes particles are sorted, are NOT sorted by sub cell when done
   also change particle icell label
------------------------------------------------------------------------- */

void Grid::assign_split_cell_particles(int icell)
{
  int ip,jcell;

  int dim = domain->dimension;
  Particle::OnePart *particles = particle->particles;
  int *next = particle->next;

  ip = cinfo[icell].first;
  while (ip >= 0) {
    if (dim == 3) jcell = update->split3d(icell,particles[ip].x);
    else jcell = update->split2d(icell,particles[ip].x);
    particles[ip].icell = jcell;
    ip = next[ip];
  }

  cinfo[icell].count = 0;
  cinfo[icell].first = -1;
}

/* ----------------------------------------------------------------------
   new surf-based methods
------------------------------------------------------------------------- */

void Grid::surf2grid_stats()
{
  double cmax,len;
  int dimension = domain->dimension;

  int scount = 0;
  int stotal = 0;
  int smax = 0;
  double sratio = BIG;
  for (int icell = 0; icell < nlocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    if (cells[icell].nsurf) scount++;
    stotal += cells[icell].nsurf;
    smax = MAX(smax,cells[icell].nsurf);
    
    cmax = MAX(cells[icell].hi[0] - cells[icell].lo[0],
	       cells[icell].hi[1] - cells[icell].lo[1]);
    if (dimension == 3) 
      cmax = MAX(cmax,cells[icell].hi[2] - cells[icell].lo[2]);
    
    if (dimension == 2) {
      for (int i = 0; i < cells[icell].nsurf; i++) {
	len = surf->line_size(cells[icell].csurfs[i]);
	sratio = MIN(sratio,len/cmax);
      }
    } else if (dimension == 3) {
      for (int i = 0; i < cells[icell].nsurf; i++) {
	surf->tri_size(cells[icell].csurfs[i],len);
	sratio = MIN(sratio,len/cmax);
      }
    }
  }
  
  int scountall,stotalall,smaxall;
  double sratioall;
  MPI_Allreduce(&scount,&scountall,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&stotal,&stotalall,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&smax,&smaxall,1,MPI_INT,MPI_MAX,world);
  MPI_Allreduce(&sratio,&sratioall,1,MPI_DOUBLE,MPI_MIN,world);
  
  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"  %d = cells with surfs\n",scountall);
      fprintf(screen,"  %d = total surfs in all grid cells\n",stotalall);
      fprintf(screen,"  %d = max surfs in one grid cell\n",smaxall);
      fprintf(screen,"  %g = min surf-size/cell-size ratio\n",sratioall);
    }
    if (logfile) {
      fprintf(logfile,"  %d = cells with surfs\n",scountall);
      fprintf(logfile,"  %d = total surfs in all grid cells\n",stotalall);
      fprintf(logfile,"  %d = max surfs in one grid cell\n",smaxall);
      fprintf(logfile,"  %g = min surf-size/cell-size ratio\n",sratioall);
    }
  }
}

/* ---------------------------------------------------------------------- */

void Grid::flow_stats()
{
  int i;

  int outside = 0;
  int inside = 0;
  int overlap = 0;
  int maxsplitone = 0;
  double cellvolume = 0.0;

  for (int icell = 0; icell < nlocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    if (cinfo[icell].type == OUTSIDE) outside++;
    else if (cinfo[icell].type == INSIDE) inside++;
    else if (cinfo[icell].type == OVERLAP) overlap++;
    maxsplitone = MAX(maxsplitone,cells[icell].nsplit);
  }

  // sum volume for unsplit and sub cells
  // skip split cells and INSIDE cells

  for (int icell = 0; icell < nlocal; icell++) {
    if (cells[icell].nsplit > 1) continue;
    if (cinfo[icell].type != INSIDE) cellvolume += cinfo[icell].volume;
  }

  int outall,inall,overall,maxsplitall;
  double cellvolumeall;
  MPI_Allreduce(&outside,&outall,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&inside,&inall,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&overlap,&overall,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&maxsplitone,&maxsplitall,1,MPI_INT,MPI_MAX,world);
  MPI_Allreduce(&cellvolume,&cellvolumeall,1,MPI_DOUBLE,MPI_SUM,world);

  double flowvolume = flow_volume();

  int *tally = new int[maxsplitall];
  int *tallyall = new int[maxsplitall];
  for (i = 0; i < maxsplitall; i++) tally[i] = 0;

  for (int icell = 0; icell < nlocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    if (cinfo[icell].type == OVERLAP) tally[cells[icell].nsplit-1]++;
  }

  MPI_Allreduce(tally,tallyall,maxsplitall,MPI_INT,MPI_SUM,world);

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"  %d %d %d = cells outside/inside/overlapping surfs\n",
	      outall,inall,overall);
      fprintf(screen," ");
      for (i = 0; i < maxsplitall; i++) fprintf(screen," %d",tallyall[i]);
      fprintf(screen," = surf cells with 1,2,etc splits\n");
      fprintf(screen,"  %.15g %.15g = cell-wise and global flow volume\n",
              cellvolumeall,flowvolume);
    }
    if (logfile) {
      fprintf(logfile,"  %d %d %d = cells outside/inside/overlapping surfs\n",
	      outall,inall,overall);
      fprintf(logfile," ");
      for (i = 0; i < maxsplitall; i++) fprintf(logfile," %d",tallyall[i]);
      fprintf(logfile," = surf cells with 1,2,etc splits\n");
      fprintf(logfile,"  %g %g = cell-wise and global flow volume\n",
              cellvolumeall,flowvolume);
    }
  }

  delete [] tally;
  delete [] tallyall;
}

/* ----------------------------------------------------------------------
   compute flow volume for entire box, using global list of surfs
   volume for one surf is projection to lower z face
   NOTE: this does not work if any surfs are clipped to zlo or zhi faces in 3d
         this does not work if any surfs are clipped to ylo or yhi faces in 3d
         need to add contribution due to closing surfs on those faces
         fairly easy to add in 2d, not so easy in 3d
------------------------------------------------------------------------- */

double Grid::flow_volume()
{
  double zarea;
  double *p1,*p2,*p3;

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  double volume = 0.0;

  if (domain->dimension == 3) {
    for (int i = 0; i < surf->ntri; i++) {
      p1 = tris[i].p1;
      p2 = tris[i].p2;
      p3 = tris[i].p3;
      zarea = 0.5 * ((p2[0]-p1[0])*(p3[1]-p1[1]) - (p2[1]-p1[1])*(p3[0]-p1[0]));
      volume -= zarea * ((p1[2]+p2[2]+p3[2])/3.0 - boxlo[2]);
    }
    if (volume <= 0.0) 
      volume += (boxhi[0]-boxlo[0]) * (boxhi[1]-boxlo[1]) * 
        (boxhi[2]-boxlo[2]); 
 
  // axisymmetric "volume" of line segment = volume of truncated cone
  // PI/3 (y1^2 + y1y2 + y2^2) (x2-x1)

  } else if (domain->axisymmetric) {
    for (int i = 0; i < surf->nline; i++) {
      p1 = lines[i].p1;
      p2 = lines[i].p2;
      volume -= 
        MY_PI3 * (p1[1]*p1[1] + p1[1]*p2[1] + p2[1]*p2[1]) * (p2[0]-p1[0]);
    }
    if (volume <= 0.0) 
      volume += MY_PI * boxhi[1]*boxhi[1] * (boxhi[0]-boxlo[0]);

  } else {
    for (int i = 0; i < surf->nline; i++) {
      p1 = lines[i].p1;
      p2 = lines[i].p2;
      volume -= (0.5*(p1[1]+p2[1]) - boxlo[1]) * (p2[0]-p1[0]);
    }
    if (volume <= 0.0) volume += (boxhi[0]-boxlo[0]) * (boxhi[1]-boxlo[1]); 
  }
  
  return volume;
}

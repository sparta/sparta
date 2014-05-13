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
#include "geometry.h"
#include "domain.h"
#include "comm.h"
#include "irregular.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

// NOTE: set DELTA to large value when done debugging

#define DELTA 10000
#define BIG 1.0e20

// default value, can be overridden by global command

#define MAXSURFPERCELL 100

enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};         // same as Domain
enum{PERIODIC,OUTFLOW,REFLECT,SURFACE,AXISYM};  // same as Domain

// cell is entirely outside/inside surfs or has any overlap with surfs
// corner pt is outside/inside surfs or is on a surf

enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};   // several files
enum{NCHILD,NPARENT,NUNKNOWN,NPBCHILD,NPBPARENT,NPBUNKNOWN,NBOUND};  // Update
enum{NOWEIGHT,VOLWEIGHT,RADWEIGHT};

// allocate space for static class variable

Grid *Grid::gptr;

// corners[i][j] = J corner points of face I of a grid cell
// works for 2d quads and 3d hexes

int corners[6][4] = {{0,2,4,6}, {1,3,5,7}, {0,1,4,5}, {2,3,6,7}, 
                     {0,1,2,3}, {4,5,6,7}};

/* ---------------------------------------------------------------------- */

Grid::Grid(SPARTA *sparta) : Pointers(sparta)
{
  exist = exist_ghost = clumped = 0;
  MPI_Comm_rank(world,&me);

  ncell = nunsplit = nsplit = nsub = 0;

  nlocal = nghost = maxlocal = maxcell = 0;
  nsplitlocal = nsplitghost = maxsplit = 0;
  nsublocal = nsubghost = 0;
  nparent = maxparent = 0;

  cells = NULL;
  cinfo = NULL;
  sinfo = NULL;
  pcells = NULL;

  maxbits = 8*sizeof(cellint)-1;

  maxsurfpercell = MAXSURFPERCELL;
  csurfs = NULL; csplits = NULL; csubs = NULL;
  allocate_surf_arrays();

  neighshift[XLO] = 0;
  neighshift[XHI] = 3;
  neighshift[YLO] = 6;
  neighshift[YHI] = 9;
  neighshift[ZLO] = 12;
  neighshift[ZHI] = 15;

  neighmask[XLO] = 7 << neighshift[XLO];
  neighmask[XHI] = 7 << neighshift[XHI];
  neighmask[YLO] = 7 << neighshift[YLO];
  neighmask[YHI] = 7 << neighshift[YHI];
  neighmask[ZLO] = 7 << neighshift[ZLO];
  neighmask[ZHI] = 7 << neighshift[ZHI];

  cutoff = -1.0;
  cellweightflag = NOWEIGHT;

  // allocate hash for cell IDs

#ifdef SPARTA_MAP
  hash = new std::map<cellint,int>();
#else
  hash = new std::tr1::unordered_map<cellint,int>();
#endif

  hashfilled = 0;
}

/* ---------------------------------------------------------------------- */

Grid::~Grid()
{
  memory->sfree(cells);
  memory->sfree(cinfo);
  memory->sfree(sinfo);
  memory->sfree(pcells);

  delete csurfs;
  delete csplits;
  delete csubs;
  delete hash;
}

/* ----------------------------------------------------------------------
   add a single owned child cell to cells and cinfo
   assume no surfs at this point, so cell and corners are OUTSIDE
   neighs and nmask will be set later
------------------------------------------------------------------------- */

void Grid::add_child_cell(cellint id, int iparent, double *lo, double *hi)
{
  grow_cells(1,1);

  int ncorner;
  if (domain->dimension == 3) ncorner = 8;
  else ncorner = 4;

  ChildCell *c = &cells[nlocal];
  c->id = id;
  c->iparent = iparent;
  c->proc = me;
  c->ilocal = nlocal;
  c->lo[0] = lo[0];
  c->lo[1] = lo[1];
  c->lo[2] = lo[2];
  c->hi[0] = hi[0];
  c->hi[1] = hi[1];
  c->hi[2] = hi[2];
  c->nsurf = 0;
  c->csurfs = NULL;
  c->nsplit = 1;
  c->isplit = -1;

  ChildInfo *ci = &cinfo[nlocal];
  ci->count = 0;
  ci->first = -1;
  ci->type = OUTSIDE;
  for (int i = 0; i < ncorner; i++) ci->corner[i] = OUTSIDE;
  ci->weight = 1.0;

  if (domain->dimension == 3) 
    ci->volume = (hi[0]-lo[0]) * (hi[1]-lo[1]) * (hi[2]-lo[2]);
  else if (domain->axisymmetric)
    ci->volume = MY_PI * (hi[1]*hi[1]-lo[1]*lo[1]) * (hi[0]-lo[0]);
  else
    ci->volume = (hi[0]-lo[0]) * (hi[1]-lo[1]);

  // increment both since are adding an unsplit cell

  nunsplitlocal++;
  nlocal++;
}

/* ----------------------------------------------------------------------
   add a single parent cell to pcells
   iparent = index of parent of this parent cell
------------------------------------------------------------------------- */

void Grid::add_parent_cell(cellint id, int iparent,
                           int nx, int ny, int nz, double *lo, double *hi)
{
  grow_pcells(1);

  ParentCell *p = &pcells[nparent];
  p->id = id;
  if (iparent >= 0) {
    p->level = pcells[iparent].level + 1;
    p->nbits = pcells[iparent].nbits + pcells[iparent].newbits;
  } else p->level = p->nbits = 0;
  p->newbits = id_bits(nx,ny,nz);
  p->iparent = iparent;

  if (p->nbits + p->newbits > maxbits) 
    error->all(FLERR,"Cell ID has too many bits");

  p->nx = nx;
  p->ny = ny;
  p->nz = nz;
  p->lo[0] = lo[0]; p->lo[1] = lo[1]; p->lo[2] = lo[2]; 
  p->hi[0] = hi[0]; p->hi[1] = hi[1]; p->hi[2] = hi[2]; 

  nparent++;
}

/* ----------------------------------------------------------------------
   add a single split cell to sinfo
   ownflag = 1/0 if split cell is owned or ghost
   sinfo fields will be set by caller
------------------------------------------------------------------------- */

void Grid::add_split_cell(int ownflag)
{
  grow_sinfo(1);
  if (ownflag) nsplitlocal++;
  else nsplitghost++;
}

/* ----------------------------------------------------------------------
   add a sub cell to cells and cinfo (if owned)
   ownflag = 1/0 if sub cell is owned or ghost
   icell = index of split cell that contains sub cell
   copy cells and cinfo for split cell to new sub cell
   sub cell fields that should be different will be set by caller
------------------------------------------------------------------------- */

void Grid::add_sub_cell(int icell, int ownflag)
{
  grow_cells(1,1);

  int inew;
  if (ownflag) inew = nlocal;
  else inew = nlocal + nghost;

  memcpy(&cells[inew],&cells[icell],sizeof(ChildCell));
  if (ownflag) memcpy(&cinfo[inew],&cinfo[icell],sizeof(ChildInfo));

  if (ownflag) {
    nsublocal++;
    nlocal++;
  } else {
    nsubghost++;
    nghost++;
  }
}

/* ----------------------------------------------------------------------
   remove ghosts and any allocated data they have
   NOTE: currently, cells just have ptrs to vectors that will be cleared
         but in future, may need to de-register the vectors
------------------------------------------------------------------------- */

void Grid::remove_ghosts()
{
  exist_ghost = 0;
  nghost = nunsplitghost = nsplitghost = nsubghost = 0;
}

/* ----------------------------------------------------------------------
   set grid stats and cpart list of owned cells with particles
------------------------------------------------------------------------- */

void Grid::setup_owned()
{
  // global counts for ncell, nunsplit, nsplit, nsub

  bigint one = nlocal - nsublocal;
  MPI_Allreduce(&one,&ncell,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  one = nunsplitlocal = nlocal - nsplitlocal - nsublocal;
  MPI_Allreduce(&one,&nunsplit,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  MPI_Allreduce(&nsplitlocal,&nsplit,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&nsublocal,&nsub,1,MPI_INT,MPI_SUM,world);
  one = nunsplitlocal = nlocal - nsplitlocal - nsublocal;

  // set cell_epsilon to 1/2 the smallest dimension of any grid cell

  int dimension = domain->dimension;

  double eps = BIG;
  for (int i = 0; i < nlocal; i++) {
    if (cells[i].nsplit <= 0) continue;
    eps = MIN(eps,cells[i].hi[0]-cells[i].lo[0]);
    eps = MIN(eps,cells[i].hi[1]-cells[i].lo[1]);
    if (dimension == 3) eps = MIN(eps,cells[i].hi[2]-cells[i].lo[2]);
  }

  MPI_Allreduce(&eps,&cell_epsilon,1,MPI_DOUBLE,MPI_MIN,world);
  cell_epsilon *= 0.5;
}

/* ----------------------------------------------------------------------
   acquire ghost cells from local cells of other procs
   method used depends on ghost cutoff
   no-op if grid is not clumped and want to acquire only nearby ghosts
------------------------------------------------------------------------- */

void Grid::acquire_ghosts()
{
  if (cutoff < 0.0) acquire_ghosts_all();
  else if (clumped) acquire_ghosts_near();
}

/* ----------------------------------------------------------------------
   acquire ghost cells from local cells of other procs
   use ring comm to get copy of all other cells in global system
------------------------------------------------------------------------- */

void Grid::acquire_ghosts_all()
{
  exist_ghost = 1;
  nempty = 0;

  // compute total # of ghosts so can pre-allocate cells array

  int nghost_new;
  MPI_Allreduce(&nlocal,&nghost_new,1,MPI_INT,MPI_SUM,world);
  nghost_new -= nlocal;
  grow_cells(nghost_new,0);

  // create buf for holding all of my cells, not including sub cells

  int sendsize = 0;
  for (int icell = 0; icell < nlocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    sendsize += pack_one(icell,NULL,0,0,0);
  }

  char *sbuf;
  memory->create(sbuf,sendsize,"grid:sbuf");
  memset(sbuf,0,sendsize);

  // pack each unsplit or split cell
  // subcells will be packed by split cell

  sendsize = 0;
  for (int icell = 0; icell < nlocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    sendsize += pack_one(icell,&sbuf[sendsize],0,0,1);
  }

  // circulate buf of my grid cells around to all procs
  // unpack augments my ghost cells with info from other procs

  gptr = this;
  comm->ring(sendsize,sizeof(char),sbuf,1,unpack_ghosts,NULL,0);

  memory->destroy(sbuf);
}

/* ----------------------------------------------------------------------
   acquire ghost cells from local cells of other procs
   use irregular comm to only get copy of nearby cells
   within extended bounding box = bounding box of owned cells + cutoff
------------------------------------------------------------------------- */

void Grid::acquire_ghosts_near()
{
  exist_ghost = 1;

  // bb lo/hi = bounding box for my owned cells

  int i;
  double bblo[3],bbhi[3];
  double *lo,*hi;

  for (i = 0; i < 3; i++) {
    bblo[i] = BIG;
    bbhi[i] = -BIG;
  }
  
  for (int icell = 0; icell < nlocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    lo = cells[icell].lo;
    hi = cells[icell].hi;
    for (i = 0; i < 3; i++) {
      bblo[i] = MIN(bblo[i],lo[i]);
      bbhi[i] = MAX(bbhi[i],hi[i]);
    }
  }

  // ebb lo/hi = bbox + grid cutoff
  // trim to simulation box in non-periodic dims
  // if bblo/hi is at periodic boundary and cutoff is 0.0,
  //   add cell_epsilon to insure ghosts across periodic boundary acquired,
  //   else may be UNKNOWN to owned cell

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  double *prd = domain->prd;
  int *bflag = domain->bflag;

  double ebblo[3],ebbhi[3];
  for (i = 0; i < 3; i++) {
    ebblo[i] = bblo[i] - cutoff;
    ebbhi[i] = bbhi[i] + cutoff;
    if (bflag[2*i] != PERIODIC) ebblo[i] = MAX(ebblo[i],boxlo[i]);
    if (bflag[2*i] != PERIODIC) ebbhi[i] = MIN(ebbhi[i],boxhi[i]);
    if (bflag[2*i] == PERIODIC && bblo[i] == boxlo[i] && cutoff == 0.0)
      ebblo[i] -= cell_epsilon;
    if (bflag[2*i] == PERIODIC && bbhi[i] == boxhi[i] && cutoff == 0.0)
      ebbhi[i] += cell_epsilon;
  }

  // box = ebbox split across periodic BC
  // 27 is max number of periodic images in 3d

  Box box[27];
  int nbox = box_periodic(ebblo,ebbhi,box);

  // boxall = collection of boxes from all procs

  int me = comm->me;
  int nprocs = comm->nprocs;

  int nboxall;
  MPI_Allreduce(&nbox,&nboxall,1,MPI_INT,MPI_SUM,world);

  int *recvcounts,*displs;
  memory->create(recvcounts,nprocs,"grid:recvcounts");
  memory->create(displs,nprocs,"grid:displs");

  int nsend = nbox*sizeof(Box);
  MPI_Allgather(&nsend,1,MPI_INT,recvcounts,1,MPI_INT,world);
  displs[0] = 0;
  for (i = 1; i < nprocs; i++) displs[i] = displs[i-1] + recvcounts[i-1];

  Box *boxall = new Box[nboxall];
  MPI_Allgatherv(box,nbox*sizeof(Box),MPI_CHAR,
                 boxall,recvcounts,displs,MPI_CHAR,world);

  memory->destroy(recvcounts);
  memory->destroy(displs);

  // nlist = # of boxes that overlap with my bbox, skipping self boxes
  // list = indices into boxall of overlaps
  // overlap = true overlap or just touching

  int nlist = 0;
  int *list;
  memory->create(list,nboxall,"grid:list");

  for (i = 0; i < nboxall; i++) {
    if (boxall[i].proc == me) continue;
    if (box_overlap(bblo,bbhi,boxall[i].lo,boxall[i].hi)) list[nlist++] = i;
  }

  // loop over my owned cells, not including sub cells
  // each may overlap with multiple boxes in list
  // on 1st pass, just tally memory to send copies of my cells
  // use lastproc to insure a cell only overlaps once per other proc
  // if oflag = 2 = my cell just touches box,
  // so flag grid cell as EMPTY ghost by setting nsurf = -1

  int j,oflag,lastproc,nsurf_hold;

  nsend = 0;
  int sendsize = 0;

  for (int icell = 0; icell < nlocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    lo = cells[icell].lo;
    hi = cells[icell].hi;
    lastproc = -1;
    for (i = 0; i < nlist; i++) {
      j = list[i];
      oflag = box_overlap(lo,hi,boxall[j].lo,boxall[j].hi);
      if (!oflag) continue;
      if (boxall[j].proc == lastproc) continue;
      lastproc = boxall[j].proc;

      if (oflag == 2) {
        nsurf_hold = cells[icell].nsurf;
        cells[icell].nsurf = -1;
      }
      sendsize += pack_one(icell,NULL,0,0,0);
      if (oflag == 2) cells[icell].nsurf = nsurf_hold;
      nsend++;
    }
  }

  // create send buf and auxiliary irregular comm vectors

  char *sbuf;
  memory->create(sbuf,sendsize,"grid:sbuf");
  memset(sbuf,0,sendsize);

  int *proclist,*sizelist;
  memory->create(proclist,nsend,"grid:proclist");
  memory->create(sizelist,nsend,"grid:sizelist");

  // on 2nd pass over local cells, fill the send buf
  // use lastproc to insure a cell only overlaps once per other proc
  // if oflag = 2 = my cell just touches box,
  // so flag grid cell as EMPTY ghost by setting nsurf = -1

  nsend = 0;
  sendsize = 0;
  for (int icell = 0; icell < nlocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    lo = cells[icell].lo;
    hi = cells[icell].hi;
    lastproc = -1;
    for (i = 0; i < nlist; i++) {
      j = list[i];
      oflag = box_overlap(lo,hi,boxall[j].lo,boxall[j].hi);
      if (!oflag) continue;
      if (boxall[j].proc == lastproc) continue;
      lastproc = boxall[j].proc;

      if (oflag == 2) {
        nsurf_hold = cells[icell].nsurf;
        cells[icell].nsurf = -1;
      }
      sizelist[nsend] = pack_one(icell,&sbuf[sendsize],0,0,1);
      if (oflag == 2) cells[icell].nsurf = nsurf_hold;
      proclist[nsend] = lastproc;
      sendsize += sizelist[nsend];
      nsend++;
    }
  }

  // clean up

  memory->destroy(list);
  delete [] boxall;

  // perform irregular communication of list on ghost cells

  Irregular *irregular = new Irregular(sparta);
  int recvsize;
  int nrecv = irregular->create_data_variable(nsend,proclist,sizelist,
                                              recvsize,comm->commsortflag);

  char *rbuf;
  memory->create(rbuf,recvsize,"grid:rbuf");
  memset(rbuf,0,recvsize);

  irregular->exchange_variable(sbuf,sizelist,rbuf);
  delete irregular;

  // unpack received grid cells as ghost cells

  int offset = 0;
  for (i = 0; i < nrecv; i++)
    offset += grid->unpack_one(&rbuf[offset],0,0);

  // more clean up

  memory->destroy(proclist);
  memory->destroy(sizelist);
  memory->destroy(sbuf);
  memory->destroy(rbuf);

  // set nempty = # of EMPTY ghost cells I store

  nempty = 0;
  for (int icell = nlocal; icell < nlocal+nghost; icell++)
    if (cells[icell].nsurf < 0) nempty++;
}

/* ----------------------------------------------------------------------
   check for overlap of 2 orthongal boxes, alo/ahi and blo/bhi
   if no overlap, return 0
   if non-grazing overlap, return 1
   if grazing overlap, return 2
   works for 2d and 3d, since in 2d both boxes have zhi > zlo
------------------------------------------------------------------------- */

void Grid::box_intersect(double *alo, double *ahi, double *blo, double *bhi,
                         double *lo, double *hi)
{
  lo[0] = MAX(alo[0],blo[0]);
  hi[0] = MIN(ahi[0],bhi[0]);
  lo[1] = MAX(alo[1],blo[1]);
  hi[1] = MIN(ahi[1],bhi[1]);
  lo[2] = MAX(alo[2],blo[2]);
  hi[2] = MIN(ahi[2],bhi[2]);
}

/* ----------------------------------------------------------------------
   check for overlap of 2 orthongal boxes, alo/ahi and blo/bhi
   if no overlap, return 0
   if non-grazing overlap, return 1
   if grazing overlap, return 2
   works for 2d and 3d, since in 2d both boxes have zhi > zlo
------------------------------------------------------------------------- */

int Grid::box_overlap(double *alo, double *ahi, double *blo, double *bhi)
{
  double lo[3],hi[3];
  box_intersect(alo,ahi,blo,bhi,lo,hi);

  if (lo[0] > hi[0]) return 0;
  if (lo[1] > hi[1]) return 0;
  if (lo[2] > hi[2]) return 0;

  if (lo[0] == hi[0]) return 2;
  if (lo[1] == hi[1]) return 2;
  if (lo[2] == hi[2]) return 2;

  return 1;
}

/* ----------------------------------------------------------------------
   split lo/hi box into multilple boxes straddling periodic boundaries
   return # of split boxes and list of new boxes in box
------------------------------------------------------------------------- */

int Grid::box_periodic(double *lo, double *hi, Box *box)
{
  int ilo,ihi,jlo,jhi,klo,khi;
  ilo = ihi = jlo = jhi = klo = khi = 0;

  int *bflag = domain->bflag;
  if (bflag[0] == PERIODIC) {
    ilo = -1; ihi = 1;
  }
  if (bflag[2] == PERIODIC) {
    jlo = -1; jhi = 1;
  }
  if (bflag[4] == PERIODIC && domain->dimension == 3) {
    klo = -1; khi = 1;
  }

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  double *prd = domain->prd;

  double plo[3],phi[3],olo[3],ohi[3];

  int n = 0;

  int i,j,k;
  for (k = klo; k <= khi; k++)
    for (j = jlo; j <= jhi; j++)
      for (i = ilo; i <= ihi; i++) {
        plo[0] = lo[0] + i*prd[0];
        phi[0] = hi[0] + i*prd[0];
        plo[1] = lo[1] + j*prd[1];
        phi[1] = hi[1] + j*prd[1];
        plo[2] = lo[2] + k*prd[2];
        phi[2] = hi[2] + k*prd[2];
        box_intersect(plo,phi,boxlo,boxhi,olo,ohi);
        if (olo[0] >= ohi[0] || olo[1] >= ohi[1] || olo[2] >= ohi[2]) continue;

        box[n].proc = me;
        box[n].lo[0] = olo[0];
        box[n].lo[1] = olo[1];
        box[n].lo[2] = olo[2];
        box[n].hi[0] = ohi[0];
        box[n].hi[1] = ohi[1];
        box[n].hi[2] = ohi[2];
        n++;
      }

  return n;
}

/* ----------------------------------------------------------------------
   find and set neighbor indices for all owned and ghost cells
   when done, hash table is valid for all owned and ghost cells
   no-op if ghosts don't exist
------------------------------------------------------------------------- */

void Grid::find_neighbors()
{
  int icell,index,nmask,boundary,periodic;
  cellint *neigh;
  cellint id;
  double *lo,*hi;
  double out[3];

  if (!exist_ghost) return;

  int dim = domain->dimension;
  int *bflag = domain->bflag;
  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  // hash all child and parent cell IDs
  // key = ID, value = index+1 for child cells, value = -(index+1) for parents
  // skip sub cells, since want neighbors to be the split cell

  hash->clear();
  hashfilled = 1;

  for (icell = 0; icell < nlocal+nghost; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    (*hash)[cells[icell].id] = icell+1;
  }
  for (icell = 0; icell < nparent; icell++)
    (*hash)[pcells[icell].id] = -(icell+1);
  
  // set neigh flags and nmask for each owned and ghost child cell
  // sub cells have same lo/hi as split cell, so their neigh info is the same

  for (icell = 0; icell < nlocal+nghost; icell++) {
    lo = cells[icell].lo;
    hi = cells[icell].hi;
    neigh = cells[icell].neigh;
    nmask = 0;

    // generate a point cell_epsilon away from face midpoint, respecting PBC
    // id_find_face() walks from root cell to find the parent or child cell
    //   furthest down heirarchy containing pt and entire lo/hi face of icell

    // XLO

    out[1] = 0.5 * (lo[1] + hi[1]);
    out[2] = 0.5 * (lo[2] + hi[2]);

    if (lo[0] == boxlo[0]) boundary = 1;
    else boundary = 0;
    if (bflag[XLO] == PERIODIC) periodic = 1;
    else periodic = 0;

    if (boundary && !periodic) {
      neigh[XLO] = 0;
      nmask = neigh_encode(NBOUND,nmask,XLO);
    } else {
      if (boundary) out[0] = boxhi[0] - cell_epsilon;
      else out[0] = lo[0] - cell_epsilon;

      id = id_find_face(out,0,0,lo,hi);

      if (hash->find(id) == hash->end()) {
        neigh[XLO] = id;
        if (boundary) nmask = neigh_encode(NPBUNKNOWN,nmask,XLO);
        else nmask = neigh_encode(NUNKNOWN,nmask,XLO);
      } else {
        index = (*hash)[id];
        if (index > 0) {
          neigh[XLO] = index-1;
          if (boundary) nmask = neigh_encode(NPBCHILD,nmask,XLO);
          else nmask = neigh_encode(NCHILD,nmask,XLO);
        } else {
          neigh[XLO] = -index-1;
          if (boundary) nmask = neigh_encode(NPBPARENT,nmask,XLO);
          else nmask = neigh_encode(NPARENT,nmask,XLO);
        }
      }
    }

    // XHI

    if (hi[0] == boxhi[0]) boundary = 1;
    else boundary = 0;
    if (bflag[XHI] == PERIODIC) periodic = 1;
    else periodic = 0;

    if (boundary && !periodic) {
      neigh[XHI] = 0;
      nmask = neigh_encode(NBOUND,nmask,XHI);
    } else {
      if (boundary) out[0] = boxlo[0] + cell_epsilon;
      else out[0] = hi[0] + cell_epsilon;

      id = id_find_face(out,0,0,lo,hi);

      if (hash->find(id) == hash->end()) {
        neigh[XHI] = id;
        if (boundary) nmask = neigh_encode(NPBUNKNOWN,nmask,XHI);
        else nmask = neigh_encode(NUNKNOWN,nmask,XHI);
      } else {
        index = (*hash)[id];
        if (index > 0) {
          neigh[XHI] = index-1;
          if (boundary) nmask = neigh_encode(NPBCHILD,nmask,XHI);
          else nmask = neigh_encode(NCHILD,nmask,XHI);
        } else {
          neigh[XHI] = -index-1;
          if (boundary) nmask = neigh_encode(NPBPARENT,nmask,XHI);
          else nmask = neigh_encode(NPARENT,nmask,XHI);
        }
      }
    }

    // YLO

    out[0] = 0.5 * (lo[0] + hi[0]);
    out[2] = 0.5 * (lo[2] + hi[2]);

    if (lo[1] == boxlo[1]) boundary = 1;
    else boundary = 0;
    if (bflag[YLO] == PERIODIC) periodic = 1;
    else periodic = 0;

    if (boundary && !periodic) {
      neigh[YLO] = 0;
      nmask = neigh_encode(NBOUND,nmask,YLO);
    } else {
      if (boundary) out[1] = boxhi[1] - cell_epsilon;
      else out[1] = lo[1] - cell_epsilon;

      id = id_find_face(out,0,1,lo,hi);

      if (hash->find(id) == hash->end()) {
        neigh[YLO] = id;
        if (boundary) nmask = neigh_encode(NPBUNKNOWN,nmask,YLO);
        else nmask = neigh_encode(NUNKNOWN,nmask,YLO);
      } else {
        index = (*hash)[id];
        if (index > 0) {
          neigh[YLO] = index-1;
          if (boundary) nmask = neigh_encode(NPBCHILD,nmask,YLO);
          else nmask = neigh_encode(NCHILD,nmask,YLO);
        } else {
          neigh[YLO] = -index-1;
          if (boundary) nmask = neigh_encode(NPBPARENT,nmask,YLO);
          else nmask = neigh_encode(NPARENT,nmask,YLO);
        }
      }
    }

    // YHI

    if (hi[1] == boxhi[1]) boundary = 1;
    else boundary = 0;
    if (bflag[YHI] == PERIODIC) periodic = 1;
    else periodic = 0;

    if (boundary && !periodic) {
      neigh[YHI] = 0;
      nmask = neigh_encode(NBOUND,nmask,YHI);
    } else {
      if (boundary) out[1] = boxlo[1] + cell_epsilon;
      else out[1] = hi[1] + cell_epsilon;

      id = id_find_face(out,0,1,lo,hi);

      if (hash->find(id) == hash->end()) {
        neigh[YHI] = id;
        if (boundary) nmask = neigh_encode(NPBUNKNOWN,nmask,YHI);
        else nmask = neigh_encode(NUNKNOWN,nmask,YHI);
      } else {
        index = (*hash)[id];
        if (index > 0) {
          neigh[YHI] = index-1;
          if (boundary) nmask = neigh_encode(NPBCHILD,nmask,YHI);
          else nmask = neigh_encode(NCHILD,nmask,YHI);
        } else {
          neigh[YHI] = -index-1;
          if (boundary) nmask = neigh_encode(NPBPARENT,nmask,YHI);
          else nmask = neigh_encode(NPARENT,nmask,YHI);
        }
      }
    }

    // ZLO
    // treat boundary as non-periodic if 2d, so is flagged as NBOUND

    out[0] = 0.5 * (lo[0] + hi[0]);
    out[1] = 0.5 * (lo[1] + hi[1]);

    if (lo[2] == boxlo[2]) boundary = 1;
    else boundary = 0;
    if (bflag[ZLO] == PERIODIC && dim == 3) periodic = 1;
    else periodic = 0;

    if (boundary && !periodic) {
      neigh[ZLO] = 0;
      nmask = neigh_encode(NBOUND,nmask,ZLO);
    } else {
      if (boundary) out[2] = boxhi[2] - cell_epsilon;
      else out[2] = lo[2] - cell_epsilon;

      id = id_find_face(out,0,2,lo,hi);

      if (hash->find(id) == hash->end()) {
        neigh[ZLO] = id;
        if (boundary) nmask = neigh_encode(NPBUNKNOWN,nmask,ZLO);
        else nmask = neigh_encode(NUNKNOWN,nmask,ZLO);
      } else {
        index = (*hash)[id];
        if (index > 0) {
          neigh[ZLO] = index-1;
          if (boundary) nmask = neigh_encode(NPBCHILD,nmask,ZLO);
          else nmask = neigh_encode(NCHILD,nmask,ZLO);
        } else {
          neigh[ZLO] = -index-1;
          if (boundary) nmask = neigh_encode(NPBPARENT,nmask,ZLO);
          else nmask = neigh_encode(NPARENT,nmask,ZLO);
        }
      }
    }

    // ZHI
    // treat boundary as non-periodic if 2d, so is flagged as NBOUND

    if (hi[2] == boxhi[2]) boundary = 1;
    else boundary = 0;
    if (bflag[ZHI] == PERIODIC && dim == 3) periodic = 1;
    else periodic = 0;

    if (boundary && !periodic) {
      neigh[ZHI] = 0;
      nmask = neigh_encode(NBOUND,nmask,ZHI);
    } else {
      if (boundary) out[2] = boxlo[2] + cell_epsilon;
      else out[2] = hi[2] + cell_epsilon;

      id = id_find_face(out,0,2,lo,hi);

      if (hash->find(id) == hash->end()) {
        neigh[ZHI] = id;
        if (boundary) nmask = neigh_encode(NPBUNKNOWN,nmask,ZHI);
        else nmask = neigh_encode(NUNKNOWN,nmask,ZHI);
      } else {
        index = (*hash)[id];
        if (index > 0) {
          neigh[ZHI] = index-1;
          if (boundary) nmask = neigh_encode(NPBCHILD,nmask,ZHI);
          else nmask = neigh_encode(NCHILD,nmask,ZHI);
        } else {
          neigh[ZHI] = -index-1;
          if (boundary) nmask = neigh_encode(NPBPARENT,nmask,ZHI);
          else nmask = neigh_encode(NPARENT,nmask,ZHI);
        }
      }
    }

    cells[icell].nmask = nmask;
  }

  // insure no UNKNOWN neighbors for owned cell
  // else cannot move particle to new proc to continue move

  int n1,n2,n3,n4,n5,n6;

  int flag = 0;
  for (icell = 0; icell < nlocal; icell++) {
    nmask = cells[icell].nmask;
    n1 = neigh_decode(nmask,XLO);
    n2 = neigh_decode(nmask,XHI);
    n3 = neigh_decode(nmask,YLO);
    n4 = neigh_decode(nmask,YHI);
    n5 = neigh_decode(nmask,ZLO);
    n6 = neigh_decode(nmask,ZHI);
    if (n1 == NUNKNOWN || n2 == NUNKNOWN || n3 == NUNKNOWN || 
        n4 == NUNKNOWN || n5 == NUNKNOWN || n6 == NUNKNOWN ) flag++;
    if (n1 == NPBUNKNOWN || n2 == NPBUNKNOWN || n3 == NPBUNKNOWN || 
        n4 == NPBUNKNOWN || n5 == NPBUNKNOWN || n6 == NPBUNKNOWN ) flag++;
  }

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);

  if (flagall) {
    char str[128];
    sprintf(str,"Owned cells with unknown neighbors = %d",flagall);
    error->all(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   change neigh[] in owned cells that are local indices to global cell IDs
   nmask is not changed
   called before cells data structure changes
     e.g. due to cell migration or new surfs inducing split cells
   can later use reset_neighbor() instead of find_neighbor() to
     change neigh[] back to local indices
   no-op if ghosts don't exist
------------------------------------------------------------------------- */

void Grid::unset_neighbors()
{
  if (!exist_ghost) return;

  int dimension = domain->dimension;

  // no change in neigh[] needed if nflag = NUNKNOWN, NPBUNKNOWN, or NBOUND

  int i,index,nmask,nflag;
  cellint *neigh;

  for (int icell = 0; icell < nlocal; icell++) {
    neigh = cells[icell].neigh;
    nmask = cells[icell].nmask;
    
    for (i = 0; i < 6; i++) {
      index = neigh[i];
      nflag = neigh_decode(nmask,i);
      if (nflag == NCHILD || nflag == NPBCHILD) 
        neigh[i] = cells[index].id;
      else if (nflag == NPARENT || nflag == NPBPARENT) 
        neigh[i] = pcells[index].id;
    }
  }
}

/* ----------------------------------------------------------------------
   reset neigh[] and nmask for owned and ghost cells after call to unset()
   neigh[] currently stores cell IDs
   convert them to local indices using hash table
   when done, hash table is valid for all owned and ghost cells
   no-op if ghosts don't exist
------------------------------------------------------------------------- */

void Grid::reset_neighbors()
{
  if (!exist_ghost) return;

  // hash all child and parent cell IDs
  // key = ID, value = index+1 for child cells, value = -(index+1) for parents
  // skip sub cells, since want neighbors to be the split cell

  hash->clear();
  hashfilled = 1;

  for (int icell = 0; icell < nlocal+nghost; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    (*hash)[cells[icell].id] = icell+1;
  }
  for (int icell = 0; icell < nparent; icell++)
    (*hash)[pcells[icell].id] = -(icell+1);

  // set neigh[] and nmask of each owned and ghost child cell
  // hash lookup can reset nmask to CHILD or UNKNOWN
  // no change in neigh[] or nmask needed if nflag = NBOUND

  int i,nmask,nflag;
  cellint *neigh;

  for (int icell = 0; icell < nlocal+nghost; icell++) {
    neigh = cells[icell].neigh;
    nmask = cells[icell].nmask;

    for (i = 0; i < 6; i++) {
      nflag = neigh_decode(nmask,i);
      if (nflag == NCHILD || nflag == NPBCHILD) {
        if (hash->find(neigh[i]) == hash->end()) {
          if (nflag == NCHILD) nmask = neigh_encode(NUNKNOWN,nmask,i);
          else nmask = neigh_encode(NPBUNKNOWN,nmask,i);
        } else neigh[i] = (*hash)[neigh[i]] - 1; 

      } else if (nflag == NPARENT || nflag == NPBPARENT) {
        neigh[i] = -(*hash)[neigh[i]] - 1; 

      } else if (nflag == NUNKNOWN || nflag == NPBUNKNOWN) {
        if (hash->find(neigh[i]) != hash->end()) {
          neigh[i] = (*hash)[neigh[i]] - 1; 
          if (nflag == NUNKNOWN) nmask = neigh_encode(NCHILD,nmask,i);
          else nmask = neigh_encode(NPBCHILD,nmask,i);
        }
      }
    }

    cells[icell].nmask = nmask;
  }
}

/* ----------------------------------------------------------------------
   set type and corner flags of all owned cells
------------------------------------------------------------------------- */

void Grid::set_inout()
{
  //set_inout_old();
  set_inout_new();
}

/* ----------------------------------------------------------------------
   set type and corner flags of all owned cells
------------------------------------------------------------------------- */

void Grid::set_inout_old()
{
  int i,m,n,icell,icorner,iface,ineigh,jcorner,jcell,ivalue,offset,iconnect;
  int nflag,ncorner,nface,nface_pts;
  int flag,mark,progress;
  double *ilo,*ihi,*jlo,*jhi;

  int faceflip[6] = {XHI,XLO,YHI,YLO,ZHI,ZLO};

  int me = comm->me;
  int dimension = domain->dimension;
  if (dimension == 3) {
    ncorner = 8;
    nface = 6;
    nface_pts = 4;
  } else {
    ncorner = 4;
    nface = 4;
    nface_pts = 2;
  }

  // enumerate connections for each corner point of each owned cell

  int **firstconnect,**nvalues;
  memory->create(firstconnect,nlocal,ncorner,"grid:firstconnect");
  memory->create(nvalues,nlocal,ncorner,"grid:nvalues");
  
  // 3d = 3 possible connections per corner per cell 
  // 2d = 2 possible connections per corner per cell 

  int maxconnect = dimension*ncorner*nlocal;
  Connect *connect = 
    (Connect *) memory->smalloc(maxconnect*sizeof(Connect),"grid:connect");
  int nconnect = 0;

  for (icell = 0; icell < nlocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    ilo = cells[icell].lo;
    ihi = cells[icell].hi;

    for (icorner = 0; icorner < ncorner; icorner++) {
      firstconnect[icell][icorner] = nconnect;

      for (iface = 0; iface < nface; iface++) {
        for (m = 0; m < nface_pts; m++) {
          if (corners[iface][m] != icorner) continue;

          // neigh cell must be child or parent that I know
          // if a parent cell find the child cell owning the matching jcorner pt

          nflag = neigh_decode(cells[icell].nmask,iface);
          
          if (nflag == NCHILD || nflag == NPBCHILD) {
            jcell = cells[icell].neigh[iface];
            jlo = cells[jcell].lo;
            jhi = cells[jcell].hi;
            jcorner = corner_compare(icorner,ilo,ihi,faceflip[iface],jlo,jhi);
            if (jcorner < 0) continue;
          } else if (nflag == NPARENT || nflag == NPBPARENT) {
            jcell = cells[icell].neigh[iface];
            jlo = pcells[jcell].lo;
            jhi = pcells[jcell].hi;
            jcorner = corner_compare(icorner,ilo,ihi,faceflip[iface],jlo,jhi);
            if (jcorner < 0) continue;
            jcell = id_child_from_parent_corner(jcell,jcorner);
            if (jcell < 0) continue;
          } else continue;

          connect[nconnect].icorner = jcorner;
          connect[nconnect].icell = cells[jcell].ilocal;
          connect[nconnect].proc = cells[jcell].proc;
          nconnect++;
        }
      }

      nvalues[icell][icorner] = nconnect - firstconnect[icell][icorner];
    }
  }

  // create irregular communicator for exchanging off-proc corner connections

  int nsend = 0;
  for (i = 0; i < nconnect; i++)
    if (connect[i].proc != me) nsend++;

  int *proclist;
  memory->create(proclist,nsend,"grid:proclist");

  nsend = 0;
  for (i = 0; i < nconnect; i++)
    if (connect[i].proc != me) proclist[nsend++] = connect[i].proc;

  Irregular *irregular = new Irregular(sparta);
  int nrecv = irregular->create_data_uniform(nsend,proclist,comm->commsortflag);
  Connect *sbuf,*rbuf;
  sbuf = (Connect *) memory->smalloc(nsend*sizeof(Connect),"grid:sbuf");
  rbuf = (Connect *) memory->smalloc(nrecv*sizeof(Connect),"grid:rbuf");

  // loop until no more progress marking corner and cell flags I own

  while (1) {
    progress = 0;

    // mark connected corner points with my corner point value
    // repeat until loop produces no changes

    while (1) {
      mark = 0;

      for (icell = 0; icell < nlocal; icell++) {
        if (cells[icell].nsplit <= 0) continue;
        for (icorner = 0; icorner < ncorner; icorner++) {
          ivalue = cinfo[icell].corner[icorner];
          
          offset = firstconnect[icell][icorner];
          n = nvalues[icell][icorner];
          for (m = 0; m < n; m++) {
            iconnect = m + offset;
            if (connect[iconnect].proc != me) {
              connect[iconnect].newvalue = ivalue;
              continue;
            } else {
              flag = mark_corner(ivalue,connect[iconnect].icell,
                                 connect[iconnect].icorner);
              if (flag < 0) {
                int jcell = connect[iconnect].icell;
                int jcorner = connect[iconnect].icorner;
                char str[128];
                sprintf(str,"Grid in/out self-mark error %d for "
                        "icell %d, icorner %d, connect %d %d, other cell %d, "
                        "other corner %d, values %d %d\n",
                        flag,icell,icorner,m,iconnect,jcell,jcorner,
                        ivalue,cinfo[jcell].corner[jcorner]);
                error->one(FLERR,str);
              }
              if (flag) mark = 1;
            }
          }
        }
      }
      
      if (mark) progress = 1;
      if (!mark) break;
    }

    // communicate off-proc connections

    nsend = 0;
    for (i = 0; i < nconnect; i++)
      if (connect[i].proc != me)
        memcpy(&sbuf[nsend++],&connect[i],sizeof(Connect));

    // perform irregular comm

    irregular->exchange_uniform((char *) sbuf,sizeof(Connect),(char *) rbuf);

    // use received connection to update cornerflags of my cells
    // set progress if updates change any of my cells

    for (m = 0; m < nrecv; m++) {
      flag = mark_corner(rbuf[m].newvalue,rbuf[m].icell,rbuf[m].icorner);
      if (flag < 0) {
        char str[128];
        sprintf(str,"Grid in/out other-mark error %d\n",flag);
        error->one(FLERR,str);
      }
      if (flag) progress = 1;
    }
    
    // if no progress made by any processor, done

    int progress_any;
    MPI_Allreduce(&progress,&progress_any,1,MPI_INT,MPI_SUM,world);
    if (!progress_any) break;
  }

  // set type and cflags for all sub cells from split cell it belongs to

  int j,splitcell;

  for (icell = 0; icell < nlocal; icell++) {
    if (cells[icell].nsplit > 0) continue;
    splitcell = sinfo[cells[icell].isplit].icell;
    cinfo[icell].type = cinfo[splitcell].type;
    for (j = 0; j < ncorner; j++)
      cinfo[icell].corner[j] = cinfo[splitcell].corner[j];
  }

  // clean up

  delete irregular;
  memory->destroy(proclist);
  memory->sfree(sbuf);
  memory->sfree(rbuf);

  memory->destroy(firstconnect);
  memory->destroy(nvalues);
  memory->sfree(connect);
}

/* ----------------------------------------------------------------------
   set type in/out for all owned cells
   mark cell and corner flags of the 4 neighbors of each cell
   UNKNOWN cell can be marked
   OVERLAP cell with UNKNOWN cornerflags can be marked
     all the cornerflags will be UNKNOWN, can be set to all INSIDE/OUTSIDE
     if set to OUTSIDE, also set flow area = entire cell
   NOTE: are periodic vs non-periodic BC handled correctly?
------------------------------------------------------------------------- */

void Grid::set_inout_new()
{
  int i,j,m,icell,jcell,pcell,itype,jtype,marktype;
  int nflag,iface,icorner,ctype,ic;
  int iset,nset,nsetnew;
  int nsend,maxsend,nrecv,maxrecv;
  int *set,*setnew,*jcorner;
  int *cflags;
  int *proclist;
  double xcorner[3];
  Connect2 *sbuf,*rbuf;

  int faceflip[6] = {XHI,XLO,YHI,YLO,ZHI,ZLO};

  int me = comm->me;
  int dimension = domain->dimension;
  int ncorner,nface,nface_pts;
  if (dimension == 3) {
    ncorner = 8;
    nface = 6;
    nface_pts = 4;
  } else {
    ncorner = 4;
    nface = 4;
    nface_pts = 2;
  }

  // create irregular communicator for exchanging off-processor cell types

  Irregular *irregular = new Irregular(sparta);

  // create set1 and set2 lists so can swap between them

  int *set1,*set2;
  memory->create(set1,nlocal*nface,"grid:set1");
  memory->create(set2,nlocal*nface,"grid:set2");

  // initial set list = overlapped cells

  nset = nsetnew = 0;
  set = set1;
  setnew = set2;

  for (icell = 0; icell < nlocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    if (cinfo[icell].type == OVERLAP) set[nset++] = icell;
  }

  // loop until no more cell types are set

  maxsend = maxrecv = 0;
  proclist = NULL;
  sbuf = rbuf = NULL;

  while (1) {
    nsend = 0;

    // process local set list
    // for unmarked neighbor cells I own, mark them and add to new set list
    // for neighbor cells I don't own, add to comm list

    while (nset) {
      nsetnew = 0;
      for (iset = 0; iset < nset; iset++) {
        icell = set[iset];
        itype = cinfo[icell].type;
        cflags = cinfo[icell].corner;

        for (iface = 0; iface < nface; iface++) {
          if (itype == OVERLAP) {
            ctype = cflags[corners[iface][0]];
            for (m = 1; m < nface_pts; m++)
              if (cflags[corners[iface][m]] != ctype) break;
            if (m < nface_pts) continue;
            if (ctype == OUTSIDE) marktype = OUTSIDE;
            else if (ctype == INSIDE) marktype = INSIDE;
            else continue;
          } else marktype = itype;

          jcell = cells[icell].neigh[iface];
          nflag = neigh_decode(cells[icell].nmask,iface);

          if (nflag == NCHILD || nflag == NPBCHILD) {
            if (cells[jcell].proc == me) {
              jtype = cinfo[jcell].type;
              if (jtype == UNKNOWN) {
                cinfo[jcell].type = marktype;
                setnew[nsetnew++] = jcell;
              } else if (jtype == OVERLAP) {
                jcorner = cinfo[jcell].corner;
                if (jcorner[0] != UNKNOWN) continue;
                for (icorner = 0; icorner < ncorner; icorner++)
                  jcorner[icorner] = marktype;
                if (marktype == OUTSIDE) {
                  double *lo = cells[jcell].lo;
                  double *hi = cells[jcell].hi;
                  if (dimension == 3) 
                    cinfo[jcell].volume = 
                      (hi[0]-lo[0]) * (hi[1]-lo[1]) * (hi[2]-lo[2]);
                  else if (domain->axisymmetric)
                    cinfo[jcell].volume = 
                      MY_PI * (hi[1]*hi[1]-lo[1]*lo[1]) * (hi[0]-lo[0]);
                  else
                    cinfo[jcell].volume = (hi[0]-lo[0]) * (hi[1]-lo[1]);
                }
                setnew[nsetnew++] = jcell;
              } else if (marktype != jtype) {
                error->one(FLERR,"Cell type mis-match when marking on self");
              }
            } else {
              if (nsend == maxsend) {
                maxsend += DELTA;
                memory->grow(proclist,maxsend,"grid:proclist");
                sbuf = (Connect2 *) 
                  memory->srealloc(sbuf,maxsend*sizeof(Connect2),"grid:sbuf");
              }
              proclist[nsend] = cells[jcell].proc;
              sbuf[nsend].newvalue = marktype;
              sbuf[nsend].icell = cells[jcell].ilocal;
              nsend++;
            }

          } else if (nflag == NPARENT || nflag == NPBPARENT) {
            pcell = jcell;
            for (m = 0; m < nface_pts; m++) {
              ic = corners[iface][m];
              if (ic % 2) xcorner[0] = cells[icell].hi[0];
              else xcorner[0] = cells[icell].lo[0];
              if ((ic/2) % 2) xcorner[1] = cells[icell].hi[1];
              else xcorner[1] = cells[icell].lo[1];
              if (ic/4) xcorner[2] = cells[icell].hi[2];
              else xcorner[2] = cells[icell].lo[2];

              if (nflag == NPBPARENT) 
                domain->uncollide(faceflip[iface],xcorner);

              jcell = id_find_child(pcell,xcorner);
              if (jcell < 0) error->one(FLERR,"Parent cell child missing");

              if (cells[jcell].proc == me) {
                jtype = cinfo[jcell].type;
                if (jtype == UNKNOWN) {
                  cinfo[jcell].type = marktype;
                  setnew[nsetnew++] = jcell;
                } else if (jtype == OVERLAP) {
                  jcorner = cinfo[jcell].corner;
                  if (jcorner[0] != UNKNOWN) continue;
                  for (icorner = 0; icorner < ncorner; icorner++)
                    jcorner[icorner] = marktype;
                  if (marktype == OUTSIDE) {
                    double *lo = cells[jcell].lo;
                    double *hi = cells[jcell].hi;
                    if (dimension == 3) 
                      cinfo[jcell].volume = 
                        (hi[0]-lo[0]) * (hi[1]-lo[1]) * (hi[2]-lo[2]);
                    else if (domain->axisymmetric)
                      cinfo[jcell].volume = 
                        MY_PI * (hi[1]*hi[1]-lo[1]*lo[1]) * (hi[0]-lo[0]);
                    else
                      cinfo[jcell].volume = (hi[0]-lo[0]) * (hi[1]-lo[1]);
                  }
                  setnew[nsetnew++] = jcell;
                } else if (marktype != jtype) {
                  error->one(FLERR,"Cell type mis-match when marking on self");
                }
              } else {
                if (nsend == maxsend) {
                  maxsend += DELTA;
                  memory->grow(proclist,maxsend,"grid:proclist");
                  sbuf = (Connect2 *) 
                    memory->srealloc(sbuf,maxsend*sizeof(Connect2),"grid:sbuf");
                }
                proclist[nsend] = cells[jcell].proc;
                sbuf[nsend].newvalue = marktype;
                sbuf[nsend].icell = cells[jcell].ilocal;
                nsend++;
              }
            }
          } else continue;
        }
      }

      // swap set lists

      nset = nsetnew;
      if (set == set1) {
        set = set2;
        setnew = set1;
      } else {
        set = set1;
        setnew = set2;
      }
    }
    
    // if no proc has info to communicate, then done iterating

    int anysend;
    MPI_Allreduce(&nsend,&anysend,1,MPI_INT,MPI_MAX,world);
    if (!anysend) break;

    // perform irregular comm of each proc's comm list
    // realloc rbuf as needed

    nrecv = irregular->create_data_uniform(nsend,proclist,comm->commsortflag);
    if (nrecv > maxrecv) {
      memory->sfree(rbuf);
      maxrecv = nrecv;
      rbuf = (Connect2 *) memory->smalloc(maxrecv*sizeof(Connect2),"grid:rbuf");
    }

    irregular->exchange_uniform((char *) sbuf,sizeof(Connect2),(char *) rbuf);

    // set cell types based on received info

    for (i = 0; i < nrecv; i++) {
      marktype = rbuf[i].newvalue;
      jcell = rbuf[i].icell;
      jtype = cinfo[jcell].type;
      if (jtype == UNKNOWN) {
        cinfo[jcell].type = marktype;
        set[nset++] = jcell;
      } else if (jtype == OVERLAP) {
        continue;
      } else if (marktype != jtype) {
        error->one(FLERR,"Cell type mis-match when marking on neigh proc");
      }
    }  
  }

  // set type and cflags for all sub cells from split cell it belongs to

  int splitcell;

  for (icell = 0; icell < nlocal; icell++) {
    if (cells[icell].nsplit > 0) continue;
    splitcell = sinfo[cells[icell].isplit].icell;
    cinfo[icell].type = cinfo[splitcell].type;
    for (j = 0; j < ncorner; j++)
      cinfo[icell].corner[j] = cinfo[splitcell].corner[j];
  }

  // clean up

  delete irregular;
  memory->destroy(set1);
  memory->destroy(set2);
  memory->destroy(proclist);
  memory->sfree(sbuf);
  memory->sfree(rbuf);
}

/* ----------------------------------------------------------------------
   look for matching corner pt on jface
------------------------------------------------------------------------- */

int Grid::corner_compare(int icorner, double *ilo, double *ihi, 
                         int jface, double *jlo, double *jhi)
{
  double x[3],y[3];

  int nface_pts = 4;
  if (domain->dimension == 2) nface_pts = 2;

  corner2point(icorner,ilo,ihi,x);
  for (int i = 0; i < nface_pts; i++) {
    corner2point(corners[jface][i],jlo,jhi,y);
    if (x[0] == y[0] && x[1] == y[1] && x[2] == y[2]) return corners[jface][i];
  }

  return -1;
}

/* ----------------------------------------------------------------------
   look for matching corner pt on jface
------------------------------------------------------------------------- */

void Grid::corner2point(int icorner, double *lo, double *hi, double *pt)
{
  if (icorner % 2 == 0) pt[0] = lo[0];
  else pt[0] = hi[0];
  if ((icorner/2) % 2 == 0) pt[1] = lo[1];
  else pt[1] = hi[1];
  if (icorner < 4) pt[2] = lo[2];
  else pt[2] = hi[2];
}

/* ----------------------------------------------------------------------
   mark icorner of icell with value
   possibly update cell type and other corner flags
   return 0 if nothing changed
   return 1 if anything changed
   return < 0 if error condition
   error -1 = corner already has a different value
   error -2 = some other corner already has a different value
   error -3 = cell was already marked INSIDE or OUTSIDE
              which means all its corners should have already been marked
------------------------------------------------------------------------- */

int Grid::mark_corner(int value, int icell, int icorner)
{
  if (value == UNKNOWN) return 0;

  int *corner = cinfo[icell].corner;
  if (corner[icorner] == value) return 0;
  if (corner[icorner] != UNKNOWN) return -1;
  corner[icorner] = value;

  // mark no further if icell has surfs

  if (cells[icell].nsurf) return 1;

  // value = INSIDE or OUTSIDE
  // mark all corners and cell type with value

  int ncorner = 8;
  if (domain->dimension == 2) ncorner = 4;

  for (int i = 0; i < ncorner; i++) {
    if (corner[i] == OVERLAP) continue;
    if (corner[i] == UNKNOWN) corner[i] = value;
    else if (corner[i] != value) return -2;
    else corner[i] = value;
  }

  // check and set cell type

  if (cinfo[icell].type == OVERLAP) return 1;
  if (cinfo[icell].type != UNKNOWN) return -3;
  if (value == INSIDE) cinfo[icell].type = INSIDE;
  else cinfo[icell].type = OUTSIDE;
  return 1;
}

/* ----------------------------------------------------------------------
   set maxlevel and check if the hierarchical grid is uniform
   uniform means all child cells have parent at same level
   thus finest level grid is effectively uniform across entire domain
   NOTE: also need to enforce that Nx,Ny,Nz 
         of all parents at same level are same?
------------------------------------------------------------------------- */

void Grid::check_uniform()
{
  // maxlevel = max level of any child cell in grid

  maxlevel = 0;
  for (int i = 0; i < nparent; i++)
    maxlevel = MAX(maxlevel,pcells[i].level);
  maxlevel++;

  // grid is uniform only if parents of all child cells are at same level

  int plevel = -1;
  for (int i = 0; i < nlocal; i++)
    plevel = MAX(plevel,pcells[cells[i].iparent].level);

  int all;
  MPI_Allreduce(&plevel,&all,1,MPI_INT,MPI_MAX,world);

  uniform = 1;
  for (int i = 0; i < nlocal; i++)
    if (pcells[cells[i].iparent].level != all) uniform = 0;

  MPI_Allreduce(&uniform,&all,1,MPI_INT,MPI_MIN,world);
  if (!all) uniform = 0;

  if (uniform) {
    int *lflag = new int[maxlevel];
    for (int i = 0; i < maxlevel; i++) lflag[i] = 0;

    int level;
    unx = uny = unz = 1;

    for (int i = 0; i < nparent; i++) {
      level = pcells[i].level;
      if (lflag[level]) continue;
      lflag[level] = 1;
      unx *= pcells[i].nx;
      uny *= pcells[i].ny;
      unz *= pcells[i].nz;
    }
    delete [] lflag;
  }
}

/* ----------------------------------------------------------------------
   require all cell types be set
   check corner pts of cells with surfs
   require corner pts on global box boundary be set
   warn if interior corner pts are not set
------------------------------------------------------------------------- */

void Grid::type_check()
{
  int i,m,icell;

  // check cell types

  int unknown = 0;
  for (int icell = 0; icell < nlocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    if (cinfo[icell].type == UNKNOWN) unknown++;
  }
  int unknownall;
  MPI_Allreduce(&unknown,&unknownall,1,MPI_INT,MPI_SUM,world);

  if (unknownall) {
    char str[128];
    sprintf(str,"Grid cells marked as unknown = %d",unknownall);
    //error->all(FLERR,str);
    if (me == 0) error->warning(FLERR,str);
  }

  // check corner flags of cells that are OVERLAP
  // warn if any interior corner flags are not set
  // error if any corner flags on global boundaries are unset

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  int dimension = domain->dimension;

  int ncorner = 4;
  if (dimension == 3) ncorner = 8;

  double x[3];
  int inside = 0;
  int outside = 0;

  for (int icell = 0; icell < nlocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    if (cinfo[icell].type != OVERLAP) continue;
    for (i = 0; i < ncorner; i++) {
      if (cinfo[icell].corner[i] != UNKNOWN) continue;
      if (i % 2 == 0) x[0] = cells[icell].lo[0];
      else x[0] = cells[icell].hi[0];
      if ((i/2) % 2 == 0) x[1] = cells[icell].lo[1];
      else x[1] = cells[icell].hi[1];
      if (dimension == 3) {
        if (i/4 == 0) x[2] = cells[icell].lo[2];
        else x[2] = cells[icell].hi[2];
      } else x[2] = 0.0;

      if (Geometry::point_on_hex(x,boxlo,boxhi)) outside++;
      else inside++;
    }
  }

  int insideall;
  MPI_Allreduce(&inside,&insideall,1,MPI_INT,MPI_SUM,world);
  if (insideall) {
    char str[128];
    sprintf(str,"Grid cell interior corner points marked as unknown = %d",
	    insideall);
    if (comm->me == 0) error->warning(FLERR,str);
  }

  int outsideall;
  MPI_Allreduce(&outside,&outsideall,1,MPI_INT,MPI_SUM,world);
  if (outsideall) {
    char str[128];
    sprintf(str,"Grid cell corner points on boundary marked as unknown = %d",
	    outsideall);
    error->all(FLERR,str);
  }

  flow_stats();
}

/* ----------------------------------------------------------------------
   set static cellwise fnum weights based on uncut volumes
   for volume, use volume of cell, whether axisymmetric or not
   for radius, use radius of cell centroid from axisymmetric axis
   called from input script
------------------------------------------------------------------------- */

void Grid::weight(int narg, char **arg)
{
  int i;
  double *lo,*hi;

  if (!exist) error->all(FLERR,"Cannot weight cells before grid is defined");
  if (narg != 1) error->all(FLERR,"Illegal weight command");

  int dimension = domain->dimension;
  int axisymmetric = domain->axisymmetric;

  if (strcmp(arg[0],"none") == 0) {
    cellweightflag = NOWEIGHT;
    for (i = 0; i < nlocal; i++) cinfo[i].weight = 1.0;

  } else if (strcmp(arg[0],"volume") == 0) {
    cellweightflag = VOLWEIGHT;

    for (int i = 0; i < nlocal; i++) {
      lo = cells[i].lo;
      hi = cells[i].hi;
      if (dimension == 3) 
        cinfo[i].weight = (hi[0]-lo[0]) * (hi[1]-lo[1]) * (hi[2]-lo[2]);
      else if (axisymmetric)
        cinfo[i].weight = MY_PI * (hi[1]*hi[1]-lo[1]*lo[1]) * (hi[0]-lo[0]);
      else
        cinfo[i].weight = (hi[0]-lo[0]) * (hi[1]-lo[1]);
    }

  } else if (strcmp(arg[0],"radius") == 0) {
    if (!axisymmetric) 
      error->all(FLERR,"Cannot use weight cell radius unless axisymmetric");
    cellweightflag = RADWEIGHT;

    for (int i = 0; i < nlocal; i++) {
      lo = cells[i].lo;
      hi = cells[i].hi;
      cinfo[i].weight = 0.5*(hi[1]+lo[1]) * (hi[0]-lo[0]);
    }

  } else error->all(FLERR,"Illegal weight command");
}

///////////////////////////////////////////////////////////////////////////
// grow cell list data structures
///////////////////////////////////////////////////////////////////////////

/* ----------------------------------------------------------------------
   insure cells and cinfo can hold N and M new cells respectively
------------------------------------------------------------------------- */

void Grid::grow_cells(int n, int m)
{
  if (nlocal+nghost+n >= maxcell) {
    int oldmax = maxcell;
    while (maxcell < nlocal+nghost+n) maxcell += DELTA;
    cells = (ChildCell *)
      memory->srealloc(cells,maxcell*sizeof(ChildCell),"grid:cells");
    memset(&cells[oldmax],0,(maxcell-oldmax)*sizeof(ChildCell));
  }

  if (nlocal+m >= maxlocal) {
    int oldmax = maxlocal;
    while (maxlocal < nlocal+m) maxlocal += DELTA;
    cinfo = (ChildInfo *)
      memory->srealloc(cinfo,maxlocal*sizeof(ChildInfo),"grid:cinfo");
    memset(&cinfo[oldmax],0,(maxlocal-oldmax)*sizeof(ChildInfo));
  }
}

/* ----------------------------------------------------------------------
   insure pcells can hold N new parent cells
------------------------------------------------------------------------- */

void Grid::grow_pcells(int n)
{
  if (nparent+n >= maxparent) {
    while (maxparent < nparent+n) maxparent += DELTA;
    pcells = (ParentCell *)
      memory->srealloc(pcells,maxparent*sizeof(ParentCell),"grid:pcells");
  }
}

/* ----------------------------------------------------------------------
   insure sinfo can hold N new split cells
------------------------------------------------------------------------- */

void Grid::grow_sinfo(int n)
{
  if (nsplitlocal+nsplitghost+n >= maxsplit) {
    while (maxsplit < nsplitlocal+nsplitghost+n) maxsplit += DELTA;
    sinfo = (SplitInfo *)
      memory->srealloc(sinfo,maxsplit*sizeof(SplitInfo),"grid:sinfo");
  }
}

/* ---------------------------------------------------------------------- */

bigint Grid::memory_usage()
{
  bigint bytes = maxcell * sizeof(ChildCell);
  bytes += maxlocal * sizeof(ChildInfo);
  bytes += maxsplit * sizeof(SplitInfo);
  bytes += nparent * sizeof(ParentCell);
  bytes += csurfs->size();
  bytes += csplits->size();

  return bytes;
}

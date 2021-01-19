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
#include "region.h"
#include "surf.h"
#include "comm.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "output.h"
#include "dump.h"
#include "irregular.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

// DEBUG
#include "update.h"

using namespace SPARTA_NS;
using namespace MathConst;

#define DELTA 8192
#define DELTAPARENT 1024
#define LARGE 256000
#define BIG 1.0e20
#define MAXGROUP 32
#define MAXLEVEL 32

// default values, can be overridden by global command

#define MAXSURFPERCELL  100
#define MAXSPLITPERCELL 10

enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};         // same as Domain
enum{PERIODIC,OUTFLOW,REFLECT,SURFACE,AXISYM};  // same as Domain
enum{REGION_ALL,REGION_ONE,REGION_CENTER};      // same as Surf
enum{PERAUTO,PERCELL,PERSURF};                  // several files

// cell type = OUTSIDE/INSIDE/OVERLAP if entirely outside/inside surfs
//   or has any overlap with surfs including grazing or touching
// corner point = OUTSIDE/INSIDE (not OVERLAP) if outside/inside
//   if exactly on surf, is still marked OUTSIDE or INSIDE by cut2d/cut3d
//   corner pts are only set if cell type = OVERLAP

enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};           // several files
enum{NCHILD,NPARENT,NUNKNOWN,NPBCHILD,NPBPARENT,NPBUNKNOWN,NBOUND};  // Update
enum{NOWEIGHT,VOLWEIGHT,RADWEIGHT};

// corners[i][j] = J corner points of face I of a grid cell
// works for 2d quads and 3d hexes

int corners[6][4] = {{0,2,4,6}, {1,3,5,7}, {0,1,4,5}, {2,3,6,7},
                     {0,1,2,3}, {4,5,6,7}};

/* ---------------------------------------------------------------------- */

Grid::Grid(SPARTA *sparta) : Pointers(sparta)
{
  exist = exist_ghost = clumped = 0;
  MPI_Comm_rank(world,&me);

  gnames = (char **) memory->smalloc(MAXGROUP*sizeof(char *),"grid:gnames");
  bitmask = (int *) memory->smalloc(MAXGROUP*sizeof(int),"grid:bitmask");
  inversemask = (int *) memory->smalloc(MAXGROUP*sizeof(int),
                                        "grid:inversemask");

  for (int i = 0; i < MAXGROUP; i++) bitmask[i] = 1 << i;
  for (int i = 0; i < MAXGROUP; i++) inversemask[i] = bitmask[i] ^ ~0;

  ngroup = 1;
  int n = strlen("all") + 1;
  gnames[0] = new char[n];
  strcpy(gnames[0],"all");

  ncell = nunsplit = nsplit = nsub = 0;

  nlocal = nghost = maxlocal = maxcell = 0;
  nsplitlocal = nsplitghost = maxsplit = 0;
  nsublocal = nsubghost = 0;
  nparent = maxparent = 0;
  maxlevel = 0;
  plevel_limit = MAXLEVEL;

  cells = NULL;
  cinfo = NULL;
  sinfo = NULL;
  pcells = NULL;

  plevels = new ParentLevel[MAXLEVEL];
  memset(plevels,0,MAXLEVEL*sizeof(ParentLevel));

  surfgrid_algorithm = PERAUTO;
  maxsurfpercell = MAXSURFPERCELL;
  maxsplitpercell = MAXSPLITPERCELL;
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

  hash = new MyHash();
  hashfilled = 0;

  copy = copymode = 0;
}

/* ---------------------------------------------------------------------- */

Grid::~Grid()
{
  if (copy || copymode) return;

  for (int i = 0; i < ngroup; i++) delete [] gnames[i];
  memory->sfree(gnames);
  memory->sfree(bitmask);
  memory->sfree(inversemask);

  memory->sfree(cells);
  memory->sfree(cinfo);
  memory->sfree(sinfo);
  memory->sfree(pcells);

  delete [] plevels;

  delete csurfs;
  delete csplits;
  delete csubs;
  delete hash;
}

/* ----------------------------------------------------------------------
   remove existing grid
   called by AdaptGrid when it reads a new grid from file
------------------------------------------------------------------------- */

void Grid::remove()
{
  memory->sfree(cells);
  memory->sfree(cinfo);
  memory->sfree(sinfo);
  memory->sfree(pcells);

  delete csurfs;
  delete csplits;
  delete csubs;

  exist_ghost = clumped = 0;
  ncell = nunsplit = nsplit = nsub = 0;
  nlocal = nghost = maxlocal = maxcell = 0;
  nsplitlocal = nsplitghost = maxsplit = 0;
  nsublocal = nsubghost = 0;
  maxlevel = 0;

  hash->clear();
  hashfilled = 0;

  cells = NULL;
  cinfo = NULL;
  sinfo = NULL;

  csurfs = NULL; csplits = NULL; csubs = NULL;
  allocate_surf_arrays();

  // NOTE: what about cutoff and cellweightflag
}

/* ----------------------------------------------------------------------
   store copy of Particle class settings
------------------------------------------------------------------------- */

void Grid::init()
{
  ncustom = particle->ncustom;
  nbytes_particle = sizeof(Particle::OnePart);
  nbytes_custom = particle->sizeof_custom();
  nbytes_total = nbytes_particle + nbytes_custom;
}

/* ----------------------------------------------------------------------
   add a single owned child cell to cells and cinfo
   assume no surfs at this point, so cell is OUTSIDE, corners are ignored
   neighs and nmask will be set later
------------------------------------------------------------------------- */

void Grid::add_child_cell(cellint id, int level, double *lo, double *hi)
{
  grow_cells(1,1);

  int ncorner;
  if (domain->dimension == 3) ncorner = 8;
  else ncorner = 4;

  ChildCell *c = &cells[nlocal];
  c->id = id;
  c->proc = me;
  c->ilocal = nlocal;
  c->level = level;
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
  ci->mask = 1;
  ci->type = OUTSIDE;
  for (int i = 0; i < ncorner; i++) ci->corner[i] = UNKNOWN;
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
   called during a run when per-processor list of grid cells may have changed
   trigger fixes, computes, dumps to change their allocated per-grid data
------------------------------------------------------------------------- */

void Grid::notify_changed()
{
  if (modify->n_pergrid) modify->grid_changed();

  Compute **compute = modify->compute;
  for (int i = 0; i < modify->ncompute; i++) {
    if (compute[i]->per_grid_flag) compute[i]->reallocate();
    if (compute[i]->per_surf_flag) compute[i]->reallocate();
  }

  for (int i = 0; i < output->ndump; i++)
    output->dump[i]->reset_grid_count();
}

/* ----------------------------------------------------------------------
   return minlevel = minlevel of any grid cell
------------------------------------------------------------------------- */

int Grid::set_minlevel()
{
  int mylevel = MAXLEVEL;
  for (int i = 0; i < nlocal; i++)
    mylevel = MIN(mylevel,cells[i].level);

  int minlevel;
  MPI_Allreduce(&mylevel,&minlevel,1,MPI_INT,MPI_MIN,world);
  return minlevel;
}

/* ----------------------------------------------------------------------
   set maxlevel = maxlevel of any grid cell
------------------------------------------------------------------------- */

void Grid::set_maxlevel()
{
  int mylevel = 0;
  for (int i = 0; i < nlocal; i++)
    mylevel = MAX(mylevel,cells[i].level);

  MPI_Allreduce(&mylevel,&maxlevel,1,MPI_INT,MPI_MAX,world);
}

/* ----------------------------------------------------------------------
   set grid stats for owned cells
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
   remove ghost grid cells and any allocated data they have
     currently, cells just have ptrs into pages that are deallocated separately
   also remove ghost surfaces, either explicit or implicit
------------------------------------------------------------------------- */

void Grid::remove_ghosts()
{
  hashfilled = 0;
  exist_ghost = 0;
  nghost = nunsplitghost = nsplitghost = nsubghost = 0;
  surf->remove_ghosts();
}

/* ----------------------------------------------------------------------
   acquire ghost cells from local cells of other procs
   if surfs are distributed, also acquire ghost cell surfs
     explicit distributed surfs require use of hash
   method used depends on ghost cutoff
   no-op if grid is not clumped and want to acquire only nearby ghosts
------------------------------------------------------------------------- */

void Grid::acquire_ghosts(int surfflag)
{
  if (surf->distributed && !surf->implicit) surf->rehash();

  if (cutoff < 0.0) acquire_ghosts_all(surfflag);
  else if (clumped) acquire_ghosts_near(surfflag);
  else if (comm->me == 0)
    error->warning(FLERR,"Could not acquire nearby ghost cells b/c "
                   "grid partition is not clumped");

  if (surf->distributed && !surf->implicit) {
    surf->hash->clear();
    surf->hashfilled = 0;
  }
}

/* ----------------------------------------------------------------------
   acquire ghost cells from local cells of other procs
   use ring comm to get copy of all other cells in global system
------------------------------------------------------------------------- */

void Grid::acquire_ghosts_all(int surfflag)
{
  exist_ghost = 1;
  nempty = 0;

  // compute total # of ghosts so can pre-allocate cells array
  // issue a memory warning if grid cell count >= LARGE and
  //   user has not specified a grid cutoff

  int nghost_new;
  MPI_Allreduce(&nlocal,&nghost_new,1,MPI_INT,MPI_SUM,world);

  if (nghost_new >= LARGE && comm->nprocs > 1 && comm->me == 0)
    error->warning(FLERR,"Per-processor grid cell memory will be large "
                   "because global gridcut < 0.0");

  nghost_new -= nlocal;
  grow_cells(nghost_new,0);

  // create buf for holding all of my cells, not including sub cells

  bigint bsendsize = 0;
  for (int icell = 0; icell < nlocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    bsendsize += pack_one(icell,NULL,0,0,surfflag,0);
  }

  if (bsendsize > MAXSMALLINT)
    error->one(FLERR,"Acquire ghosts all send buffer exceeds 2 GB");
  int sendsize = bsendsize;

  char *sbuf;
  memory->create(sbuf,sendsize,"grid:sbuf");
  memset(sbuf,0,sendsize);

  // pack each unsplit or split cell
  // subcells will be packed by split cell

  sendsize = 0;
  for (int icell = 0; icell < nlocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    sendsize += pack_one(icell,&sbuf[sendsize],0,0,surfflag,1);
  }

  // circulate buf of my grid cells around to all procs
  // unpack augments my ghost cells with info from other procs

  unpack_ghosts_surfflag = surfflag;
  comm->ring(sendsize,sizeof(char),sbuf,1,unpack_ghosts,NULL,0,(void *) this);

  memory->destroy(sbuf);
}

/* ----------------------------------------------------------------------
   acquire ghost cells from local cells of other procs
   use irregular comm to only get copy of nearby cells
   within extended bounding box = bounding box of owned cells + cutoff
------------------------------------------------------------------------- */

void Grid::acquire_ghosts_near(int surfflag)
{
  if (update->have_mem_limit())
    return acquire_ghosts_near_less_memory(surfflag);

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
  MPI_Allgatherv(box,nsend,MPI_CHAR,boxall,recvcounts,displs,MPI_CHAR,world);

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
  bigint bsendsize = 0;

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
      bsendsize += pack_one(icell,NULL,0,0,surfflag,0);
      if (oflag == 2) cells[icell].nsurf = nsurf_hold;
      nsend++;
    }
  }

  if (bsendsize > MAXSMALLINT)
    error->one(FLERR,"Acquire ghosts near send buffer exceeds 2 GB");
  int sendsize = bsendsize;

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
      sizelist[nsend] = pack_one(icell,&sbuf[sendsize],0,0,surfflag,1);
      if (oflag == 2) cells[icell].nsurf = nsurf_hold;
      proclist[nsend] = lastproc;
      sendsize += sizelist[nsend];
      nsend++;
    }
  }

  // clean up

  memory->destroy(list);
  delete [] boxall;

  // perform irregular communication of list of ghost cells

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
    offset += grid->unpack_one(&rbuf[offset],0,0,surfflag);

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
   acquire ghost cells from local cells of other procs
   use irregular comm to only get copy of nearby cells
   within extended bounding box = bounding box of owned cells + cutoff
------------------------------------------------------------------------- */

void Grid::acquire_ghosts_near_less_memory(int surfflag)
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
  MPI_Allgatherv(box,nsend,MPI_CHAR,boxall,recvcounts,displs,MPI_CHAR,world);

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

  int icell_start = 0;
  int icell_end = nlocal;
  int not_done = 1;

  while (not_done) {
    nsend = 0;
    bigint bsendsize = 0;

    for (int icell = icell_start; icell < nlocal; icell++) {
      icell_end = icell+1;
      if (cells[icell].nsplit <= 0) continue;
      lo = cells[icell].lo;
      hi = cells[icell].hi;
      lastproc = -1;
      int break_flag = 0;
      for (i = 0; i < nlist; i++) {
        j = list[i];
        oflag = box_overlap(lo,hi,boxall[j].lo,boxall[j].hi);
        if (!oflag) continue;
        if (boxall[j].proc == lastproc) continue;
        lastproc = boxall[j].proc;

        int n = pack_one(icell,NULL,0,0,surfflag,0);
        if (n > 0 && bsendsize > 0 && bsendsize+n > update->global_mem_limit) {
          icell_end -= 1;
          break_flag = 1;
          break;
        }
        bsendsize += n;
        nsend++;
      }
      if (break_flag) break;
    }

    if (bsendsize > MAXSMALLINT)
      error->one(FLERR,"Acquire ghosts near send buffer exceeds 2 GB");
    int sendsize = bsendsize;

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
    for (int icell = icell_start; icell < icell_end; icell++) {
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
        sizelist[nsend] = pack_one(icell,&sbuf[sendsize],0,0,surfflag,1);
        if (oflag == 2) cells[icell].nsurf = nsurf_hold;
        proclist[nsend] = lastproc;
        sendsize += sizelist[nsend];
        nsend++;
      }
    }

    // perform irregular communication of list of ghost cells

    Irregular *irregular = new Irregular(sparta);
    int recvsize;
    int nrecv = irregular->create_data_variable(nsend,proclist,sizelist,
                                                recvsize,comm->commsortflag);

    char *rbuf;
    memory->create(rbuf,recvsize,"grid:rbuf");
    memset(rbuf,0,recvsize);

    irregular->exchange_variable(sbuf,sizelist,rbuf);
    delete irregular;
    irregular = NULL;

    // unpack received grid cells as ghost cells

    int offset = 0;
    for (i = 0; i < nrecv; i++)
      offset += grid->unpack_one(&rbuf[offset],0,0,surfflag);

    // more clean up

    memory->destroy(proclist);
    memory->destroy(sizelist);
    memory->destroy(sbuf);
    memory->destroy(rbuf);

    icell_start = icell_end;
    int not_done_local = icell_start < nlocal;
    MPI_Allreduce(&not_done_local,&not_done,1,MPI_INT,MPI_SUM,world);

  } // end while loop

  // clean up

  memory->destroy(list);
  delete [] boxall;

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
   hash all my child IDs, owned + ghost
------------------------------------------------------------------------- */

void Grid::rehash()
{
  // hash all owned/ghost child cell IDs
  // key = ID, value = index for child cells
  // skip sub cells

  hash->clear();

  for (int icell = 0; icell < nlocal+nghost; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    (*hash)[cells[icell].id] = icell;
  }

  hashfilled = 1;
}

/* ----------------------------------------------------------------------
   find and set neighbor indices for all owned and ghost cells
   when done, hash table is valid for all owned and ghost cells
   no-op if ghosts don't exist
------------------------------------------------------------------------- */

void Grid::find_neighbors()
{
  int icell,idim,iface,ilevel,level,nmask,boundary,periodic;
  int found,face_touching,unknownflag;
  cellint id,neighID,refineID,coarsenID;
  cellint *neigh;
  double *lo,*hi;

  if (!exist_ghost) return;

  int dimension = domain->dimension;
  int *bflag = domain->bflag;
  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  // insure cell IDs (owned + ghost) are hashed

  rehash();

  // clear parent cells data structure since will (re)build it here

  nparent = 0;

  // set neigh flags and nmask for each owned and ghost child cell
  // skip sub cells, loop below copies their info from original split cell

  int nunknown = 0;

  for (icell = 0; icell < nlocal+nghost; icell++) {
    if (cells[icell].nsplit <= 0) continue;

    id = cells[icell].id;
    level = cells[icell].level;
    lo = cells[icell].lo;
    hi = cells[icell].hi;
    neigh = cells[icell].neigh;
    nmask = 0;
    unknownflag = 0;

    for (iface = 0; iface < 6; iface++) {
      idim = iface/2;

      // set boundary and periodic flags for this face
      // treat 2d Z boundaries as non-periodic

      if (iface % 2 == 0 && lo[idim] == boxlo[idim]) boundary = 1;
      else if (iface % 2 == 1 && hi[idim] == boxhi[idim]) boundary = 1;
      else boundary = 0;
      if (bflag[iface] == PERIODIC) periodic = 1;
      else periodic = 0;
      if (dimension == 2 && (iface == ZLO || iface == ZHI))
	periodic = 0;

      // face = non-periodic boundary, neighbor is BOUND

      if (boundary && !periodic) {
	neigh[iface] = 0;
	nmask = neigh_encode(NBOUND,nmask,iface);
	continue;
      }

      // neighID = ID of neighbor cell at same level as icell

      neighID = id_neigh_same_parent(id,level,iface);
      if (neighID == 0) neighID = id_neigh_same_level(id,level,iface);

      // if in hash, neighbor is CHILD

      if (hash->find(neighID) != hash->end()) {
	neigh[iface] = (*hash)[neighID];
	if (!boundary) nmask = neigh_encode(NCHILD,nmask,iface);
	else nmask = neigh_encode(NPBCHILD,nmask,iface);
	continue;
      }

      // refine from neighID until reach maxlevel
      // look for a child cell on touching face I do own (or ghost)
      // if find one, neighbor is PARENT, add its ID and lo/hi to local pcells

      face_touching = iface % 2 ? iface-1 : iface+1;
      refineID = neighID;
      ilevel = level;
      found = 0;

      while (ilevel < maxlevel) {
	refineID = id_refine(refineID,ilevel,face_touching);
	if (hash->find(refineID) != hash->end()) {
	  neigh[iface] = nparent;
	  if (!boundary) nmask = neigh_encode(NPARENT,nmask,iface);
	  else nmask = neigh_encode(NPBPARENT,nmask,iface);
	
	  if (nparent == maxparent) grow_pcells();
	  pcells[nparent].id = neighID;
	  id_lohi(neighID,level,boxlo,boxhi,pcells[nparent].lo,pcells[nparent].hi);
	  nparent++;
	
	  found = 1;
	  break;
	} else ilevel++;
      }
      if (found) continue;

      // coarsen from neighID until reach top level
      // if find one, neighbor is CHILD

      coarsenID = neighID;
      ilevel = level;
      found = 0;

      while (ilevel > 1) {
	coarsenID = id_coarsen(coarsenID,ilevel);
	if (hash->find(coarsenID) != hash->end()) {
	  neigh[iface] = (*hash)[coarsenID];
	  if (!boundary) nmask = neigh_encode(NCHILD,nmask,iface);
	  else nmask = neigh_encode(NPBCHILD,nmask,iface);
	  found = 1;
	  break;
	} else ilevel--;
      }
      if (found) continue;

      // found nothing, so UNKNOWN neighbor
      // should never happend for an owned cell (error check below)
      // can happen for a ghost cell

      neigh[iface] = 0;
      if (!boundary) nmask = neigh_encode(NUNKNOWN,nmask,iface);
      else nmask = neigh_encode(NPBUNKNOWN,nmask,iface);

      if (icell < nlocal) unknownflag = 1;
    }

    cells[icell].nmask = nmask;
    if (unknownflag) nunknown++;
  }

  // for sub cells, copy neighbor info from original split cell

  int m,splitcell;

  for (icell = 0; icell < nlocal+nghost; icell++) {
    if (cells[icell].nsplit >= 1) continue;
    splitcell = sinfo[cells[icell].isplit].icell;
    for (m = 0; m < 6; m++)
      cells[icell].neigh[m] = cells[splitcell].neigh[m];
    cells[icell].nmask = cells[splitcell].nmask;
  }

  // error if any UNKNOWN neighbors for an owned cell
  // cannot move particle to new proc to continue move

  int nall;
  MPI_Allreduce(&nunknown,&nall,1,MPI_INT,MPI_SUM,world);

  if (nall) {
    char str[128];
    sprintf(str,"Owned cells with unknown neighbors = %d",nall);
    error->all(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   change neigh[] in owned cells that are local indices to global cell IDs
   nmask is not changed
   called before cells data structure changes
     e.g. due to cell migration or new surfs inducing split cells
   can later use reset_neighbors() instead of find_neighbors() to
     change neigh[] back to local indices
   no-op if ghosts don't exist
------------------------------------------------------------------------- */

void Grid::unset_neighbors()
{
  if (!exist_ghost) return;

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
        neigh[i] = pcells[neigh[i]].id;
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

  // insure all cell IDs (owned + ghost) are hashed

  rehash();

  // clear parent cells data structure since will (re)build it here

  nparent = 0;

  // set neigh[] and nmask of each owned and ghost child cell
  // hash lookup can reset nmask to CHILD or UNKNOWN
  // no change in neigh[] or nmask needed if nflag = NBOUND

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  int i,level,nmask,nflag;
  cellint *neigh;

  for (int icell = 0; icell < nlocal+nghost; icell++) {
    level = cells[icell].level;
    neigh = cells[icell].neigh;
    nmask = cells[icell].nmask;

    for (i = 0; i < 6; i++) {
      nflag = neigh_decode(nmask,i);
      if (nflag == NCHILD || nflag == NPBCHILD) {
        if (hash->find(neigh[i]) == hash->end()) {
          if (nflag == NCHILD) nmask = neigh_encode(NUNKNOWN,nmask,i);
          else nmask = neigh_encode(NPBUNKNOWN,nmask,i);
        } else neigh[i] = (*hash)[neigh[i]];

      } else if (nflag == NPARENT || nflag == NPBPARENT) {
	if (nparent == maxparent) grow_pcells();
	pcells[nparent].id = neigh[i];
	id_lohi(neigh[i],level,boxlo,boxhi,pcells[nparent].lo,pcells[nparent].hi);
        neigh[i] = nparent;
	nparent++;

      } else if (nflag == NUNKNOWN || nflag == NPBUNKNOWN) {
        if (hash->find(neigh[i]) != hash->end()) {
          neigh[i] = (*hash)[neigh[i]];
          if (nflag == NUNKNOWN) nmask = neigh_encode(NCHILD,nmask,i);
          else nmask = neigh_encode(NPBCHILD,nmask,i);
        }
      }
    }

    cells[icell].nmask = nmask;
  }
}

/* ----------------------------------------------------------------------
   set type = INSIDE/OUTSIDE for all owned cells if not already set
   input:
     all cells marked as OVERLAP or UNKNOWN
     most OVERLAP cells have their corner points marked as INSIDE/OUTSIDE
       a corner point is never marked as OVERLAP
     some OVERLAP cells may have all corner points UNKNOWN
       due to Cut2d/Cut3d split() failing to mark them
       b/c surf(s) only touch cell at point(s) on cell surface
         and not enough info to determine if cell otherwise inside vs outside
       for these cells, volume = 0.0
   output:
     all UNKNOWN cells marked as INSIDE/OUTSIDE
     mark corner pts of any OVERLAP cells with UNKNOWN corner pts
       to be all INSIDE or all OUTSIDE, set cell volume accordingly
   output status will be checked when type_check() is invoked
     errors or warnings will be triggered
   NOTE: if UNKNOWN corner pts of OVERLAP cells are not set,
     may just be flagged as warning in type_check(),
     but it means a non-zero volume for a cell that
       is effectively OUTSIDE will not be set correctly
------------------------------------------------------------------------- */

void Grid::set_inout()
{
  int i,j,m,icell,jcell,itype,jtype,marktype;
  int nflag,iface,icorner,ctype,ic;
  int iset,nset,nsetnew;
  int nsend,maxsend,nrecv,maxrecv;
  int *set,*setnew,*jcorner;
  int *cflags;
  int *proclist;
  double xcorner[3];
  ParentCell *pcell;
  Connect *sbuf,*rbuf;

  if (!exist_ghost)
    error->all(FLERR,"Cannot mark grid cells as inside/outside surfs because "
               "ghost cells do not exist");

  // if all surfs are transparent, just mark all cells as OUTSIDE

  if (surf->all_transparent()) {
    for (icell = 0; icell < nlocal; icell++) cinfo[icell].type = OUTSIDE;
    return;
  }

  // set dimensional dependent quantities

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

  // initial set list = overlapped cells with corner values which are set

  nset = nsetnew = 0;
  set = set1;
  setnew = set2;

  for (icell = 0; icell < nlocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    if (cinfo[icell].type == OVERLAP && cinfo[icell].corner[0] != UNKNOWN)
      set[nset++] = icell;
  }

  // loop until no more cell types are set

  maxsend = maxrecv = 0;
  proclist = NULL;
  sbuf = rbuf = NULL;

  while (1) {
    nsend = 0;

    // process local set list
    // for unmarked neighbor cells I own, mark them and add to new set list
    // for neighbor cells I would mark but don't own, add to comm list

    while (nset) {
      nsetnew = 0;
      for (iset = 0; iset < nset; iset++) {
        icell = set[iset];
        itype = cinfo[icell].type;
        cflags = cinfo[icell].corner;

        // loop over my cell's faces

        for (iface = 0; iface < nface; iface++) {

          // if I am an OUTSIDE/INSIDE cell: marktype = my itype
          // if I am an OVERLAP cell:
          //   can only mark neighbor if all corner pts on face are same
          //   marktype = value of those corner pts = OUTSIDE or INSIDE

          if (itype != OVERLAP) marktype = itype;
          else {
            ctype = cflags[corners[iface][0]];
            for (m = 1; m < nface_pts; m++)
              if (cflags[corners[iface][m]] != ctype) break;
            if (m < nface_pts) continue;
            if (ctype == OUTSIDE) marktype = OUTSIDE;
            else if (ctype == INSIDE) marktype = INSIDE;
            else continue;
          }

          jcell = cells[icell].neigh[iface];
          nflag = neigh_decode(cells[icell].nmask,iface);

          if (nflag == NCHILD || nflag == NPBCHILD) {

            // this proc owns neighbor cell
            // perform this marking logic (here and below):
            // (1) if jtype = UNKNOWN:
            //     set jtype = marktype = OUTSIDE/INSIDE
            //     add jcell to new set
            // (2) else if jtype = OVERLAP:
            //     skip if jcell's corners are already marked
            //     set jcell corner pts = marktype = OUTSIDE/INSIDE
            //     also set its volume to full cell or zero
            //     do not add jcell to new set, since is OVERLAP cell
            //       maybe could add, but not sure need to,
            //       and could be marked wrong ?
            // (3) else jcell = OUTSIDE/INSIDE:
            //     if jcell mark is different than marktype,
            //     error b/c markings are inconsistent

            if (cells[jcell].proc == me) {
              jtype = cinfo[jcell].type;

              if (jtype == UNKNOWN) {
                cinfo[jcell].type = marktype;
                setnew[nsetnew++] = jcell;
              } else if (jtype == OVERLAP) {
                // don't think this restriction needed (Mar 2017)
                //if (itype == OVERLAP) continue;
                jcorner = cinfo[jcell].corner;
                if (jcorner[0] != UNKNOWN) continue;
                for (icorner = 0; icorner < ncorner; icorner++)
                  jcorner[icorner] = marktype;
                if (marktype == INSIDE) cinfo[jcell].volume = 0.0;
                else if (marktype == OUTSIDE) {
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
                // do not add to setnew, see comment above
                //setnew[nsetnew++] = jcell;
              } else {
                if (jtype != marktype) {
                  printf("CELL1 proc %d icell %d id " CELLINT_FORMAT " iface %d "
                         "jcell %d id " CELLINT_FORMAT " marktype %d jtype %d\n",
                         me,icell,cells[icell].id,iface,jcell,cells[jcell].id,
                         marktype,jtype);
                  error->one(FLERR,"Cell type mis-match when marking on self");
                }
              }

            // this proc does not own neighbor cell
            // pack sbuf with info to send

            } else {
              if (nsend == maxsend) {
                maxsend += DELTA;
                memory->grow(proclist,maxsend,"grid:proclist");
                sbuf = (Connect *)
                  memory->srealloc(sbuf,maxsend*sizeof(Connect),"grid:sbuf");
              }
              proclist[nsend] = cells[jcell].proc;
              sbuf[nsend].itype = itype;
              sbuf[nsend].marktype = marktype;
              sbuf[nsend].jcell = cells[jcell].ilocal;
              nsend++;
            }

          } else if (nflag == NPARENT || nflag == NPBPARENT) {

            // neighbor is a parent cell, not a child cell
            // for each corner points of icell face:
            //   find jcell = child cell owner of the face corner pt
            //   if I own the child cell, mark it in same manner as above

	    pcell = &pcells[jcell];
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

              jcell = id_find_child(pcell->id,cells[icell].level,
				    pcell->lo,pcell->hi,xcorner);
              if (jcell < 0) error->one(FLERR,"Parent cell child missing");

              // this proc owns neighbor cell
              // perform same marking logic as above

              if (cells[jcell].proc == me) {
                jtype = cinfo[jcell].type;

                if (jtype == UNKNOWN) {
                  cinfo[jcell].type = marktype;
                  setnew[nsetnew++] = jcell;
                } else if (jtype == OVERLAP) {
                  // don't think this restriction needed (Mar 2017)
                  //if (itype == OVERLAP) continue;
                  jcorner = cinfo[jcell].corner;
                  if (jcorner[0] != UNKNOWN) continue;
                  for (icorner = 0; icorner < ncorner; icorner++)
                    jcorner[icorner] = marktype;
                  if (marktype == INSIDE) cinfo[jcell].volume = 0.0;
                  else if (marktype == OUTSIDE) {
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
                  // do not add to setnew, see comment above
                  //setnew[nsetnew++] = jcell;
                } else {
                  if (jtype != marktype) {
                    printf("CELL2 proc %d icell %d id " CELLINT_FORMAT
                           " iface %d "
                           "jcell %d id " CELLINT_FORMAT " marktype %d "
                           "jtype %d\n",
                           me,icell,cells[icell].id,iface,jcell,cells[jcell].id,
                           marktype,jtype);
                    error->one(FLERR,
                               "Cell type mis-match when marking on self");
                  }
                }

              // this proc does not own neighbor cell
              // pack sbuf with info to send

              } else {
                if (nsend == maxsend) {
                  maxsend += DELTA;
                  memory->grow(proclist,maxsend,"grid:proclist");
                  sbuf = (Connect *)
                    memory->srealloc(sbuf,maxsend*sizeof(Connect),"grid:sbuf");
                }
                proclist[nsend] = cells[jcell].proc;
                sbuf[nsend].itype = itype;
                sbuf[nsend].marktype = marktype;
                sbuf[nsend].jcell = cells[jcell].ilocal;
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
      rbuf = (Connect *) memory->smalloc(maxrecv*sizeof(Connect),"grid:rbuf");
    }

    irregular->exchange_uniform((char *) sbuf,sizeof(Connect),(char *) rbuf);

    // this proc received info to mark neighbor cell it owns
    // perform same marking logic as above

    for (i = 0; i < nrecv; i++) {
      itype = rbuf[i].itype;
      marktype = rbuf[i].marktype;
      jcell = rbuf[i].jcell;
      jtype = cinfo[jcell].type;

      if (jtype == UNKNOWN) {
        cinfo[jcell].type = marktype;
        set[nset++] = jcell;
      } else if (jtype == OVERLAP) {
        // don't think this restriction needed (Mar 2017)
        //if (itype == OVERLAP) continue;
        jcorner = cinfo[jcell].corner;
        if (jcorner[0] != UNKNOWN) continue;
        for (icorner = 0; icorner < ncorner; icorner++)
          jcorner[icorner] = marktype;
        if (marktype == INSIDE) cinfo[jcell].volume = 0.0;
        else if (marktype == OUTSIDE) {
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
        // do not add to setnew, see comment above
        //set[nset++] = jcell;
      } else {
        if (marktype != jtype) {
          printf("CELL3 me %d jcell %d id " CELLINT_FORMAT
                 " marktype %d jtype %d\n",
                 me,jcell,cells[jcell].id,marktype,jtype);
          error->one(FLERR,"Cell type mis-match when marking on neigh proc");
        }
      }
    }
  }

  // NOTE: at this point could make a final attempt to mark
  //   any remaining UNKNOWN corner pts of an overlap cell
  //   to avoid warnings and errors in type_check()
  // when doing grid adaptation, a new cell could possibly still be unmarked,
  //   because its already marked neighbors which are non-OVERLAP cells
  //     will not try to mark it
  //   this is in contrast to first-time marking, when sweep entire grid
  // the final attempt logic could do this:
  //   instead of having marked cells mark their unmarked neighbors
  //   have unmarked cells look at their neighbors to acquire markings
  //   not necessary when doing initial full-grid sweep,
  //     but may be necessary when doing incremental adaptation

  // all done with marking
  // set type and cflags for all sub cells from split cell it belongs to

  int splitcell;

  for (icell = 0; icell < nlocal; icell++) {
    if (cells[icell].nsplit > 0) continue;
    splitcell = sinfo[cells[icell].isplit].icell;
    cinfo[icell].type = cinfo[splitcell].type;
    for (j = 0; j < ncorner; j++)
      cinfo[icell].corner[j] = cinfo[splitcell].corner[j];
  }

  // set volume of cells that are now INSIDE to 0.0
  // this allows error check in Collide and FixGridCheck for particles
  //   in zero-volume cells

  for (icell = 0; icell < nlocal; icell++)
    if (cinfo[icell].type == INSIDE) cinfo[icell].volume = 0.0;

  /*
  printf("POST INOUT %d: %d\n",comm->me,grid->nlocal);
  for (int i = 0; i < grid->nlocal; i++) {
    Grid::ChildCell *g = &grid->cells[i];
    if (g->id == 52)
    printf("ICELL %d: %d id %d pid %d lo %g %g "
           "hi %g %g type %d corners %d %d %d %d vol %g\n",
           comm->me,i,g->id,pcells[g->iparent].id,
           g->lo[0],
           g->lo[1],
           g->hi[0],
           g->hi[1],
           grid->cinfo[i].type,
           grid->cinfo[i].corner[0],
           grid->cinfo[i].corner[1],
           grid->cinfo[i].corner[2],
           grid->cinfo[i].corner[3],grid->cinfo[i].volume);
  }
  */

  // clean up

  delete irregular;
  memory->destroy(set1);
  memory->destroy(set2);
  memory->destroy(proclist);
  memory->sfree(sbuf);
  memory->sfree(rbuf);
}

/* ----------------------------------------------------------------------
   check if the hierarchical grid is uniform
   uniform means all child cells are at same level
   thus finest level grid is effectively uniform across entire domain
   if uniform, set unx,uny,unz
------------------------------------------------------------------------- */

void Grid::check_uniform()
{
  // uniform = 1 only if all child cells are at maxlevel

  int myuniform = 1;
  for (int i = 0; i < nlocal; i++)
    if (cells[i].level != maxlevel) myuniform = 0;

  MPI_Allreduce(&myuniform,&uniform,1,MPI_INT,MPI_MIN,world);

  // if grid is uniform, compute grid extent in each dimension

  if (uniform) {
    unx = uny = unz = 1;
    for (int i = 0; i < maxlevel; i++) {
      unx *= plevels[i].nx;
      uny *= plevels[i].ny;
      unz *= plevels[i].nz;
    }
  }
}

/* ----------------------------------------------------------------------
   require all cell types be set
   check if all corner pts of OVERLAP cells are set
     error if corner pts on global box boundary are not set
     warning if corner pts interior to simulation box are not set
   check that no OUTSIDE cell has zero volume
   flag = 1 to output flow stats (default)
------------------------------------------------------------------------- */

void Grid::type_check(int outflag)
{
  int i;

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
    error->all(FLERR,str);
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

      if (Geometry::point_on_hex(x,boxlo,boxhi)) {
        printf("BAD CORNER icell %d id %d type %d "
               "icorner %d x %g %g %g cflags %d %d %d %d\n",
               icell,cells[icell].id,cinfo[icell].type,i,x[0],x[1],x[2],
               cinfo[icell].corner[0],
               cinfo[icell].corner[1],
               cinfo[icell].corner[2],
               cinfo[icell].corner[3]);
        outside++;
      }
      else inside++;
    }
  }

  int insideall;
  MPI_Allreduce(&inside,&insideall,1,MPI_INT,MPI_SUM,world);
  if (insideall) {
    char str[128];
    sprintf(str,"Grid cell interior corner points marked as unknown "
            "(volume will be wrong if cell is effectively outside) = %d",
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

  int volzero = 0;
  for (int icell = 0; icell < nlocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    if (cinfo[icell].type == OUTSIDE && cinfo[icell].volume == 0.0) volzero++;
  }
  int volzeroall;
  MPI_Allreduce(&volzero,&volzeroall,1,MPI_INT,MPI_SUM,world);
  if (outsideall) {
    char str[128];
    sprintf(str,"Grid cells marked outside, but with zero volume = %d",
	    volzeroall);
    error->all(FLERR,str);
  }

  if (outflag) flow_stats();
}

/* ----------------------------------------------------------------------
   set static cellwise fnum weights based on uncut volumes
   for volume, use volume of cell, whether axisymmetric or not
   for radius, use radius of cell centroid from axisymmetric axis
   weight() called from input script and read_restart (with narg = -1)
   weight_one() is called only for adapted cells from grid_adapt
------------------------------------------------------------------------- */

void Grid::weight(int narg, char **arg)
{
  int i;

  if (!exist) error->all(FLERR,"Cannot weight cells before grid is defined");
  if (narg > 0 && narg != 1) error->all(FLERR,"Illegal weight command");

  // if called from read_restart with narg = -1, cellweightflag is already set

  if (narg == 1) {
    if (strcmp(arg[0],"none") == 0) cellweightflag = NOWEIGHT;
    else if (strcmp(arg[0],"volume") == 0) cellweightflag = VOLWEIGHT;
    else if (strcmp(arg[0],"radius") == 0) cellweightflag = RADWEIGHT;
    else error->all(FLERR,"Illegal weight command");
  }

  if (cellweightflag == RADWEIGHT && !domain->axisymmetric)
    error->all(FLERR,"Cannot use weight cell radius unless axisymmetric");

  // set per-cell weights

  for (i = 0; i < nlocal; i++) weight_one(i);
}

void Grid::weight_one(int icell)
{
  double *lo,*hi;

  int dimension = domain->dimension;
  int axisymmetric = domain->axisymmetric;

  if (cellweightflag == NOWEIGHT) {
    cinfo[icell].weight = 1.0;
  } else if (cellweightflag == VOLWEIGHT) {
    lo = cells[icell].lo;
    hi = cells[icell].hi;
    if (dimension == 3)
      cinfo[icell].weight = (hi[0]-lo[0]) * (hi[1]-lo[1]) * (hi[2]-lo[2]);
    else if (axisymmetric)
      cinfo[icell].weight = MY_PI * (hi[1]*hi[1]-lo[1]*lo[1]) * (hi[0]-lo[0]);
    else
      cinfo[icell].weight = (hi[0]-lo[0]) * (hi[1]-lo[1]);
  } else if (cellweightflag == RADWEIGHT) {
    lo = cells[icell].lo;
    hi = cells[icell].hi;
    cinfo[icell].weight = 0.5*(hi[1]+lo[1]) * (hi[0]-lo[0]);
  }
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
   grow pcells
   NOTE: need this function for Kokkos
------------------------------------------------------------------------- */

void Grid::grow_pcells()
{
  maxparent += DELTAPARENT;
  pcells = (ParentCell *)
     memory->srealloc(pcells,maxparent*sizeof(ParentCell),"grid:pcells");
}

/* ----------------------------------------------------------------------
   insure sinfo can hold N new split cells
------------------------------------------------------------------------- */

void Grid::grow_sinfo(int n)
{
  if (nsplitlocal+nsplitghost+n >= maxsplit) {
    int oldmax = maxsplit;
    while (maxsplit < nsplitlocal+nsplitghost+n) maxsplit += DELTA;
    sinfo = (SplitInfo *)
      memory->srealloc(sinfo,maxsplit*sizeof(SplitInfo),"grid:sinfo");
    memset(&sinfo[oldmax],0,(maxsplit-oldmax)*sizeof(SplitInfo));
  }
}

/* ----------------------------------------------------------------------
   group grid command called via input script
------------------------------------------------------------------------- */

void Grid::group(int narg, char **arg)
{
  int i,flag;
  double x[3];

  if (narg < 3) error->all(FLERR,"Illegal group command");

  int dimension = domain->dimension;

  int igroup = find_group(arg[0]);
  if (igroup < 0) igroup = add_group(arg[0]);
  int bit = bitmask[igroup];

  // style = region
  // add grid cell to group if in region

  if (strcmp(arg[2],"region") == 0) {
    if (narg != 5) error->all(FLERR,"Illegal group command");
    int iregion = domain->find_region(arg[3]);
    if (iregion == -1) error->all(FLERR,"Group region ID does not exist");
    Region *region = domain->regions[iregion];

    int rstyle;
    if (strcmp(arg[4],"all") == 0) rstyle = REGION_ALL;
    else if (strcmp(arg[4],"one") == 0) rstyle = REGION_ONE;
    else if (strcmp(arg[4],"center") == 0) rstyle = REGION_CENTER;
    else error->all(FLERR,"Illegal group command");

    if (dimension == 2) x[2] = 0.0;

    if (rstyle == REGION_ALL) {
      for (i = 0; i < nlocal; i++) {
        flag = 1;

        if (dimension == 3) x[2] = cells[i].lo[2];
        x[0] = cells[i].lo[0];
        x[1] = cells[i].lo[1];
        if (!region->match(x)) flag = 0;
        x[0] = cells[i].hi[0];
        x[1] = cells[i].lo[1];
        if (!region->match(x)) flag = 0;
        x[0] = cells[i].lo[0];
        x[1] = cells[i].hi[1];
        if (!region->match(x)) flag = 0;
        x[0] = cells[i].hi[0];
        x[1] = cells[i].hi[1];
        if (!region->match(x)) flag = 0;

        if (dimension == 3) x[2] = cells[i].hi[2];
        if (dimension == 3) {
          x[0] = cells[i].lo[0];
          x[1] = cells[i].lo[1];
          if (!region->match(x)) flag = 0;
          x[0] = cells[i].hi[0];
          x[1] = cells[i].lo[1];
          if (!region->match(x)) flag = 0;
          x[0] = cells[i].lo[0];
          x[1] = cells[i].hi[1];
          if (!region->match(x)) flag = 0;
          x[0] = cells[i].hi[0];
          x[1] = cells[i].hi[1];
          if (!region->match(x)) flag = 0;
        }

        if (flag) cinfo[i].mask |= bit;
      }

    } else if (rstyle == REGION_ONE) {
      for (i = 0; i < nlocal; i++) {
        flag = 0;

        if (dimension == 3) x[2] = cells[i].lo[2];
        x[0] = cells[i].lo[0];
        x[1] = cells[i].lo[1];
        if (region->match(x)) flag = 1;

        x[0] = cells[i].hi[0];
        x[1] = cells[i].lo[1];
        if (region->match(x)) flag = 1;
        x[0] = cells[i].lo[0];
        x[1] = cells[i].hi[1];
        if (region->match(x)) flag = 1;
        x[0] = cells[i].hi[0];
        x[1] = cells[i].hi[1];
        if (region->match(x)) flag = 1;

        if (dimension == 3) x[2] = cells[i].hi[2];
        if (dimension == 3) {
          x[0] = cells[i].lo[0];
          x[1] = cells[i].lo[1];
          if (region->match(x)) flag = 1;
          x[0] = cells[i].hi[0];
          x[1] = cells[i].lo[1];
          if (region->match(x)) flag = 1;
          x[0] = cells[i].lo[0];
          x[1] = cells[i].hi[1];
          if (region->match(x)) flag = 1;
          x[0] = cells[i].hi[0];
          x[1] = cells[i].hi[1];
          if (region->match(x)) flag = 1;
        }

        if (flag) cinfo[i].mask |= bit;
      }

    } else if (rstyle == REGION_CENTER) {
      for (i = 0; i < nlocal; i++) {
        x[0] = 0.5 * (cells[i].lo[0] + cells[i].hi[0]);
        x[1] = 0.5 * (cells[i].lo[1] + cells[i].hi[1]);
        if (dimension == 3) x[2] = 0.5 * (cells[i].lo[2] + cells[i].hi[2]);
        if (region->match(x)) cinfo[i].mask |= bit;
      }
    }

  // style = subtract

  } else if (strcmp(arg[2],"subtract") == 0) {
    if (narg < 5) error->all(FLERR,"Illegal group command");

    int length = narg-3;
    int *list = new int[length];

    int jgroup;
    for (int iarg = 3; iarg < narg; iarg++) {
      jgroup = find_group(arg[iarg]);
      if (jgroup == -1) error->all(FLERR,"Group ID does not exist");
      list[iarg-3] = jgroup;
    }

    // add to group if in 1st group in list

    int otherbit = bitmask[list[0]];

    for (i = 0; i < nlocal; i++)
      if (cinfo[i].mask & otherbit) cinfo[i].mask |= bit;

    // remove grid cells if they are in any of the other groups
    // AND with inverse mask removes the grid cell from group

    int inverse = inversemask[igroup];

    for (int ilist = 1; ilist < length; ilist++) {
      otherbit = bitmask[list[ilist]];
      for (i = 0; i < nlocal; i++)
        if (cinfo[i].mask & otherbit) cinfo[i].mask &= inverse;
    }

    delete [] list;

  // style = union

  } else if (strcmp(arg[2],"union") == 0) {
    if (narg < 4) error->all(FLERR,"Illegal group command");

    int length = narg-3;
    int *list = new int[length];

    int jgroup;
    for (int iarg = 3; iarg < narg; iarg++) {
      jgroup = find_group(arg[iarg]);
      if (jgroup == -1) error->all(FLERR,"Group ID does not exist");
      list[iarg-3] = jgroup;
    }

    // add to group if in any other group in list

    int otherbit;

    for (int ilist = 0; ilist < length; ilist++) {
      otherbit = bitmask[list[ilist]];
      for (i = 0; i < nlocal; i++)
        if (cinfo[i].mask & otherbit) cinfo[i].mask |= bit;
    }

    delete [] list;

  // style = intersect

  } else if (strcmp(arg[2],"intersect") == 0) {
    if (narg < 5) error->all(FLERR,"Illegal group command");

    int length = narg-3;
    int *list = new int[length];

    int jgroup;
    for (int iarg = 3; iarg < narg; iarg++) {
      jgroup = find_group(arg[iarg]);
      if (jgroup == -1) error->all(FLERR,"Group ID does not exist");
      list[iarg-3] = jgroup;
    }

    // add to group if in all groups in list

    int otherbit,ok,ilist;

    for (i = 0; i < nlocal; i++) {
      ok = 1;
      for (ilist = 0; ilist < length; ilist++) {
        otherbit = bitmask[list[ilist]];
        if ((cinfo[i].mask & otherbit) == 0) ok = 0;
      }
      if (ok) cinfo[i].mask |= bit;
    }

    delete [] list;

  // style = clear
  // remove all grid cells from group

  } else if (strcmp(arg[2],"clear") == 0) {
    if (igroup == 0) error->all(FLERR,"Cannot clear group all");
    int inversebits = inversemask[igroup];

    for (i = 0; i < nlocal; i++) cinfo[i].mask &= inversebits;
  }

  // print stats for changed group

  bigint n = 0;
  for (i = 0; i < nlocal; i++)
    if (cinfo[i].mask & bit) n++;

  bigint nall;
  MPI_Allreduce(&n,&nall,1,MPI_SPARTA_BIGINT,MPI_SUM,world);

  if (comm->me == 0) {
    if (screen)
      fprintf(screen,BIGINT_FORMAT " grid cells in group %s\n",
              nall,gnames[igroup]);
    if (logfile)
      fprintf(logfile,BIGINT_FORMAT " grid cells in group %s\n",
              nall,gnames[igroup]);
  }
}

/* ----------------------------------------------------------------------
   add a new grid group ID, assumed to be unique
------------------------------------------------------------------------- */

int Grid::add_group(const char *id)
{
  if (ngroup == MAXGROUP)
    error->all(FLERR,"Cannot have more than 32 surface groups");

  int n = strlen(id) + 1;
  gnames[ngroup] = new char[n];
  strcpy(gnames[ngroup],id);
  ngroup++;
  return ngroup-1;
}

/* ----------------------------------------------------------------------
   find a grid group ID
   return index of group or -1 if not found
------------------------------------------------------------------------- */

int Grid::find_group(const char *id)
{
  int igroup;
  for (igroup = 0; igroup < ngroup; igroup++)
    if (strcmp(id,gnames[igroup]) == 0) break;
  if (igroup == ngroup) return -1;
  return igroup;
}

/* ----------------------------------------------------------------------
   check if a grid igroup is a uniform grid
   no cells with surfs
   all cells are at same level
   group forms a contiguous 3d block of cells
   return count of my cells in the group
   set nxyz = extent of 3d block in each dim
   set corner = lower left corner of 3d block
   set xyzsize = size of one grid cell (all are same size)
------------------------------------------------------------------------- */

int Grid::check_uniform_group(int igroup, int *nxyz,
                              double *corner, double *xyzsize)
{
  double lo[3],hi[3],onesize[3];

  int sflag = 0;
  int minlev = maxlevel;
  int maxlev = 0;

  int count = 0;
  lo[0] = domain->boxhi[0];
  lo[1] = domain->boxhi[1];
  lo[2] = domain->boxhi[2];
  hi[0] = domain->boxlo[0];
  hi[1] = domain->boxlo[1];
  hi[2] = domain->boxlo[2];

  int groupbit = bitmask[igroup];

  for (int icell = 0; icell < nlocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsurf) sflag++;
    minlev = MIN(minlev,cells[icell].level);
    maxlev = MAX(maxlev,cells[icell].level);
    lo[0] = MIN(lo[0],cells[icell].lo[0]);
    lo[1] = MIN(lo[1],cells[icell].lo[1]);
    lo[2] = MIN(lo[2],cells[icell].lo[2]);
    hi[0] = MAX(hi[0],cells[icell].hi[0]);
    hi[1] = MAX(hi[1],cells[icell].hi[1]);
    hi[2] = MAX(hi[2],cells[icell].hi[2]);
    count++;
  }

  // check that no cells already have surfs

  int allsflag;
  MPI_Allreduce(&sflag,&allsflag,1,MPI_INT,MPI_SUM,world);
  if (allsflag) {
    char str[128];
    sprintf(str,
            "Read_isurfs adding surfs to %d cells which already have surfs",
            allsflag);
    error->all(FLERR,str);
  }

  // check that all cells are at same level

  int allminlev,allmaxlev;
  MPI_Allreduce(&minlev,&allminlev,1,MPI_INT,MPI_MIN,world);
  MPI_Allreduce(&maxlev,&allmaxlev,1,MPI_INT,MPI_MAX,world);
  if (allminlev != allmaxlev)
    error->all(FLERR,"Read_isurfs grid group is not all at uniform level");

  // check that cell count matches a contiguous block of cells
  // xyzsize = size of one cell at allmaxlev

  xyzsize[0] = domain->xprd;
  xyzsize[1] = domain->yprd;
  xyzsize[2] = domain->zprd;

  for (int ilevel = 0; ilevel < allmaxlev; ilevel++) {
    xyzsize[0] /= plevels[ilevel].nx;
    xyzsize[1] /= plevels[ilevel].ny;
    xyzsize[2] /= plevels[ilevel].nz;
  }

  double alllo[3],allhi[3];
  MPI_Allreduce(lo,&alllo,3,MPI_DOUBLE,MPI_MIN,world);
  MPI_Allreduce(hi,&allhi,3,MPI_DOUBLE,MPI_MAX,world);

  corner[0] = alllo[0];
  corner[1] = alllo[1];
  corner[2] = alllo[2];

  nxyz[0] = static_cast<int> ((allhi[0]-alllo[0])/xyzsize[0] + 0.5);
  nxyz[1] = static_cast<int> ((allhi[1]-alllo[1])/xyzsize[1] + 0.5);
  nxyz[2] = static_cast<int> ((allhi[2]-alllo[2])/xyzsize[2] + 0.5);

  bigint allbcount;
  bigint bcount = count;
  MPI_Allreduce(&bcount,&allbcount,1,MPI_SPARTA_BIGINT,MPI_SUM,world);

  if ((bigint) nxyz[0] * nxyz[1]*nxyz[2] != allbcount)
    error->all(FLERR,"Read_isurfs grid group is not a contiguous brick");

  return count;
}

/* ----------------------------------------------------------------------
   proc 0 writes grid info and group info to restart file
------------------------------------------------------------------------- */

void Grid::write_restart(FILE *fp)
{
  fwrite(&maxlevel,sizeof(int),1,fp);
  fwrite(plevels,sizeof(ParentLevel),maxlevel,fp);

  fwrite(&ngroup,sizeof(int),1,fp);

  int n;
  for (int i = 0; i < ngroup; i++) {
    n = strlen(gnames[i]) + 1;
    fwrite(&n,sizeof(int),1,fp);
    fwrite(gnames[i],sizeof(char),n,fp);
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads grid info and group info from restart file
   bcast to other procs
------------------------------------------------------------------------- */

void Grid::read_restart(FILE *fp)
{
  // read_restart may have reset maxsurfpercell in its header() method
  // if so, need to reallocate surf arrays to correct max length

  if (maxsurfpercell != MAXSURFPERCELL) allocate_surf_arrays();

  // read level info

  if (me == 0) {
    fread(&maxlevel,sizeof(int),1,fp);
    fread(plevels,sizeof(ParentLevel),maxlevel,fp);
  }
  MPI_Bcast(&maxlevel,1,MPI_INT,0,world);
  MPI_Bcast(plevels,maxlevel*sizeof(ParentLevel),MPI_CHAR,0,world);

  // if any exist, clear existing group names, before reading new ones

  for (int i = 0; i < ngroup; i++) delete [] gnames[i];

  if (me == 0) fread(&ngroup,sizeof(int),1,fp);
  MPI_Bcast(&ngroup,1,MPI_INT,0,world);

  int n;
  for (int i = 0; i < ngroup; i++) {
    if (me == 0) fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    gnames[i] = new char[n];
    if (me == 0) fread(gnames[i],sizeof(char),n,fp);
    MPI_Bcast(gnames[i],n,MPI_CHAR,0,world);
  }
}

/* ----------------------------------------------------------------------
   return size of child grid restart info for this proc
   using count of all owned cells
  // NOTE: worry about N overflowing int, and in IROUNDUP ???
------------------------------------------------------------------------- */

int Grid::size_restart()
{
  int n = 2*sizeof(int);
  n = IROUNDUP(n);
  n += nlocal * sizeof(cellint);
  n = IROUNDUP(n);
  n += nlocal * sizeof(int);
  n = IROUNDUP(n);
  n += nlocal * sizeof(int);
  n = IROUNDUP(n);
  return n;
}

/* ----------------------------------------------------------------------
   return size of child grid restart info
   using nlocal_restart count of all owned cells
------------------------------------------------------------------------- */

int Grid::size_restart(int nlocal_restart)
{
  int n = 2*sizeof(int);
  n = IROUNDUP(n);
  n += nlocal_restart * sizeof(cellint);
  n = IROUNDUP(n);
  n += nlocal_restart * sizeof(int);
  n = IROUNDUP(n);
  n += nlocal_restart * sizeof(int);
  n = IROUNDUP(n);
  return n;
}

/* ----------------------------------------------------------------------
   pack my child grid info into buf
   nlocal, clumped as scalars
   ID, level, nsplit as vectors for all owned cells
   // NOTE: worry about N overflowing int, and in IROUNDUP ???
------------------------------------------------------------------------- */

int Grid::pack_restart(char *buf)
{
  int n;

  int *ibuf = (int *) buf;
  ibuf[0] = nlocal;
  ibuf[1] = clumped;
  n = 2*sizeof(int);
  n = IROUNDUP(n);

  cellint *cbuf = (cellint *) &buf[n];
  for (int i = 0; i < nlocal; i++)
    cbuf[i] = cells[i].id;
  n += nlocal * sizeof(cellint);
  n = IROUNDUP(n);

  ibuf = (int *) &buf[n];
  for (int i = 0; i < nlocal; i++)
    ibuf[i] = cells[i].level;
  n += nlocal * sizeof(int);
  n = IROUNDUP(n);

  ibuf = (int *) &buf[n];
  for (int i = 0; i < nlocal; i++)
    ibuf[i] = cells[i].nsplit;
  n += nlocal * sizeof(int);
  n = IROUNDUP(n);

  return n;
}

/* ----------------------------------------------------------------------
   unpack child grid info into restart storage
   nlocal_restart, clumped as scalars
   id_restart, nsplit_restart as vectors
   allocate vectors here, will be deallocated by ReadRestart
------------------------------------------------------------------------- */

int Grid::unpack_restart(char *buf)
{
  int n;

  int *ibuf = (int *) buf;
  nlocal_restart = ibuf[0];
  clumped = ibuf[1];
  n = 2*sizeof(int);
  n = IROUNDUP(n);

  memory->create(id_restart,nlocal_restart,"grid:id_restart");
  memory->create(level_restart,nlocal_restart,"grid:nlevel_restart");
  memory->create(nsplit_restart,nlocal_restart,"grid:nsplit_restart");

  cellint *cbuf = (cellint *) &buf[n];
  for (int i = 0; i < nlocal_restart; i++)
    id_restart[i] = cbuf[i];
  n += nlocal_restart * sizeof(cellint);
  n = IROUNDUP(n);

  ibuf = (int *) &buf[n];
  for (int i = 0; i < nlocal_restart; i++)
    level_restart[i] = ibuf[i];
  n += nlocal_restart * sizeof(int);
  n = IROUNDUP(n);

  ibuf = (int *) &buf[n];
  for (int i = 0; i < nlocal_restart; i++)
    nsplit_restart[i] = ibuf[i];
  n += nlocal_restart * sizeof(int);
  n = IROUNDUP(n);

  return n;
}

/* ---------------------------------------------------------------------- */

bigint Grid::memory_usage()
{
  bigint bytes = maxcell * sizeof(ChildCell);
  bytes += maxlocal * sizeof(ChildInfo);
  bytes += maxsplit * sizeof(SplitInfo);
  bytes += csurfs->size();
  bytes += csplits->size();

  return bytes;
}

/* ---------------------------------------------------------------------- */

void Grid::debug()
{
  for (int i = 0; i < nlocal; i++) {
    printf("GRID %d " CELLINT_FORMAT ": \n",i,cells[i].id);
    printf("  neigh " CELLINT_FORMAT " " CELLINT_FORMAT " "
           CELLINT_FORMAT " " CELLINT_FORMAT " "
           CELLINT_FORMAT " " CELLINT_FORMAT "\n",
           cells[i].neigh[0],cells[i].neigh[1],cells[i].neigh[2],
           cells[i].neigh[3],cells[i].neigh[4],cells[i].neigh[5]);
    printf("  lohi %g %g %g: %g %g %g\n",
           cells[i].lo[0],cells[i].lo[1],cells[i].lo[2],
           cells[i].hi[0],cells[i].hi[1],cells[i].hi[2]);
    printf("  nsurf %d:",cells[i].nsurf);
    for (int j = 0; j < cells[i].nsurf; j++)
      printf(" %d",cells[i].csurfs[j]);
    printf("\n");
    printf("  nsplit %d isplit %d\n",cells[i].nsplit,cells[i].isplit);
    printf("  type %d corner %d %d %d %d %d %d %d %d\n",
           cinfo[i].type,
           cinfo[i].corner[0],cinfo[i].corner[1],cinfo[i].corner[2],
           cinfo[i].corner[3],cinfo[i].corner[4],cinfo[i].corner[5],
           cinfo[i].corner[6],cinfo[i].corner[7]);
    printf("  volume %g\n",cinfo[i].volume);
  }
}

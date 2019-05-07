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
#include "my_page.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

// prototype for non-class function

int compare_surfIDs(const void *, const void *);

#define BIG 1.0e20
#define CHUNK 16

enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};   // several files
enum{PERAUTO,PERCELL,PERSURF};          // several files

// operations for surfaces in grid cells

/* ----------------------------------------------------------------------
   map surf elements to grid cells for explicit surfs, distributed or not
   via one of two algorithms
   cell_alg = original, loop over my cells, check all surfs within bbox
   surf_alg = Jan19, loop over N/P surfs, find small set of cells each overlaps,
              perform rendezvous comm to convert cells per surf to surfs per cell
   for distributed surfs, have to use surf_alg
   PERAUTO option chooses based on total nsurfs vs nprocs
   called from Readsurf, MoveSurf, RemoveSurf, ReadRestart, and FixMoveSurf
------------------------------------------------------------------------- */

void Grid::surf2grid(int subflag, int outflag)
{
  if (surf->distributed)
    surf2grid_surf_algorithm(subflag,outflag);
  else if (surfgrid_algorithm == PERAUTO) {
    if (comm->nprocs > surf->nsurf) surf2grid_cell_algorithm(outflag);
    else surf2grid_surf_algorithm(subflag,outflag);
  } else if (surfgrid_algorithm == PERCELL) {
    surf2grid_cell_algorithm(outflag);
  } else if (surfgrid_algorithm == PERSURF) {
    surf2grid_surf_algorithm(subflag,outflag);
  }
  
  // now have nsurf,csurfs list of local surfs that overlap each cell
  // compute cut volume and split info for each cell

  surf2grid_split(subflag,outflag);
}

/* ----------------------------------------------------------------------
   find surfs that overlap owned grid cells
   algorithm: for each of my cells, check all surfs
   in cells: set nsurf, csurfs
   in cinfo: set type=OVERLAP for cells with surfs
------------------------------------------------------------------------- */

void Grid::surf2grid_cell_algorithm(int outflag)
{
  int nsurf;
  double t1,t2;
  surfint *ptr;
  double *lo,*hi;

  int dim = domain->dimension;

  if (outflag) {
    MPI_Barrier(world);
    t1 = MPI_Wtime();
  }

  surf->setup_bbox();

  double *slo = surf->bblo;
  double *shi = surf->bbhi;

  if (dim == 3) cut3d = new Cut3d(sparta);
  else cut2d = new Cut2d(sparta,domain->axisymmetric);

  // compute overlap of surfs with each cell I own
  // info stored in nsurf,csurfs
  // skip if nsplit <= 0 b/c split cells could exist if restarting

  int max = 0;
  
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

    if (nsurf > maxsurfpercell) {
      max = MAX(max,nsurf);
      csurfs->vgot(0);
    } else if (nsurf) {
      csurfs->vgot(nsurf);
      cells[icell].nsurf = nsurf;
      cells[icell].csurfs = ptr;
      cinfo[icell].type = OVERLAP;
    }
  }

  // error if surf count exceeds maxsurfpercell in any cell

  int maxall;
  MPI_Allreduce(&max,&maxall,1,MPI_INT,MPI_MAX,world);
  if (maxall) {
    if (me == 0) printf("Max surfs in any cell = %d\n",maxall);
    error->all(FLERR,"Too many surfs in one cell - set global surfmax");
  }

  // timing info
  
  if (outflag) {
    MPI_Barrier(world);
    t2 = MPI_Wtime();
    tmap = t2-t1;
    trvous1 = trvous2 = 0.0;
  }

  if (outflag) surf2grid_stats();
}

/* ----------------------------------------------------------------------
   find surfs that overlap owned grid cells
   algorithm: for each of my N/P surfs, find small set of overlapping cells
     use rendezvous algorithm to convert cells per surf to surfs per cell
     for distributed, use rvous alg to acquire nlocal surfs from nown owners
   in cells: set nsurf, csurfs
   in cinfo: set type=OVERLAP for cells with surfs
------------------------------------------------------------------------- */

void Grid::surf2grid_surf_algorithm(int subflag, int outflag)
{
  int i,j,m,icell,isurf;
  double t1,t2,t3,t4;

  int dim = domain->dimension;
  int distributed = surf->distributed;

  if (dim == 3) cut3d = new Cut3d(sparta);
  else cut2d = new Cut2d(sparta,domain->axisymmetric);

  if (outflag) {
    MPI_Barrier(world);
    t1 = MPI_Wtime();
  }

  // if no hash or subflag, reset hash for parent/child IDs
  // needed b/c callers called clear_surf before surf2grid
  //   to wipe out split cells and compress local cell list

  if (!hashfilled || subflag) rehash();

  // stage 1 ----------------------------------------------------
  // nsurf = my N/P surfs for finding overlaps with any cells
  // if distributed, my surfs = nown mylines/mytris
  // else, my surfs = every 1/Pth from lines/tris
  // cellcount = # of cells my surfs overlap with
  // celllist = list of cell IDs my surfs overlap with

  Surf::Line *surf_lines;
  Surf::Tri *surf_tris;
  int nsurf,istart,istop,idelta;
  int nprocs = comm->nprocs;

  if (distributed) {
    surf_lines = surf->mylines;
    surf_tris = surf->mytris;
    nsurf = surf->nown;
    istart = 0;
    istop = nsurf;
    idelta = 1;
  } else {
    surf_lines = surf->lines;
    surf_tris = surf->tris;
    int ntotal = surf->nsurf;
    nsurf = ntotal / nprocs;
    if (me < ntotal % nprocs) nsurf++;
    istart = comm->me;
    istop = ntotal;
    idelta = nprocs;
  }

  int *cellcount;
  memory->create(cellcount,nsurf,"surf2grid2:cellcount");
  cellint **celllist = 
    (cellint **) memory->smalloc(nsurf*sizeof(cellint *),"surf2grid2:celllist");

  // for my surfs: populate cellcount and celllist via find_overlaps()
  // done recursively using bounding box of surf and grid structure
  // overlap test is performed for each cell that overlaps surf bbox
  //   via call to cut2d/3d->surf2grid_one()

  int ncell;
  cellint *ptr;

  int max = 0;
  m = 0;
  for (isurf = istart; isurf < istop; isurf += idelta) {
    ptr = cpsurf->vget();
    ncell = find_overlaps(isurf,ptr);

    if (ncell > maxcellpersurf) {
      max = MAX(max,ncell);
      cpsurf->vgot(0);
    } else {
      cpsurf->vgot(ncell);
      cellcount[m] = ncell;
      celllist[m] = ptr;
      m++;
    }
  }

  // error if cell count exceeds maxcellpersurf for any surf

  int maxall;
  MPI_Allreduce(&max,&maxall,1,MPI_INT,MPI_MAX,world);
  if (maxall) {
    if (me == 0) printf("Max cells overlapping any surf = %d\n",maxall);
    error->all(FLERR,"Too many cells overlap one surf - set global cellmax");
  }

  if (outflag) {
    MPI_Barrier(world);
    t2 = MPI_Wtime();
  }

  // stage 2 ----------------------------------------------------
  // rendezvous to convert list of cells per surf to list of surfs per cell
  // ncount = # of my datums to send to rendevous procs
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
    inbuf[m].proc = me;
    inbuf[m].cellID = cells[i].id;
    inbuf[m].surfID = 0;
    m++;
  }

  surfint surfID;

  for (i = 0; i < nsurf; i++) {
    if (distributed) {
      if (dim == 2) surfID = surf_lines[i].id;
      else surfID = surf_tris[i].id;
    } else {
      if (dim == 2) surfID = surf_lines[me+i*nprocs].id;
      else surfID = surf_tris[me+i*nprocs].id;
    }

    for (j = 0; j < cellcount[i]; j++) {
      proclist[m] = hashlittle(&celllist[i][j],sizeof(cellint),0) % nprocs;
      inbuf[m].proc = -1;
      inbuf[m].cellID = celllist[i][j];
      inbuf[m].surfID = surfID;
      m++;
    }
  }

  memory->destroy(cellcount);
  memory->sfree(celllist);
  cpsurf->reset();

  // perform rendezvous operation
  // each proc owns random subset of cells
  // receives all info to form and return their surf lists

  char *buf;
  int nreturn = comm->rendezvous(1,ncount,(char *) inbuf,sizeof(InRvous),
                                 0,proclist,rendezvous_surflist,
                                 0,buf,sizeof(OutRvous),(void *) this);
  OutRvous *outbuf = (OutRvous *) buf;

  memory->destroy(proclist);
  memory->sfree(inbuf);

  // set cells nsurf for all my cells based on output info
  // output datums = pairs of cellIDs and surfIDs
  // first pass counts surfs/cell and allocates csurfs
  // skip if nsplit <= 0 b/c split cells could already exist if restarting

  for (icell = 0; icell < nlocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    cells[icell].nsurf = 0;
    cells[icell].csurfs = NULL;
  }

  for (m = 0; m < nreturn; m++) {
    icell = (*hash)[outbuf[m].cellID] - 1;
    cells[icell].nsurf++;
  }

  max = 0;
  for (icell = 0; icell < nlocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    nsurf = cells[icell].nsurf;
    max = MAX(max,nsurf);
    if (nsurf) cells[icell].csurfs = csurfs->get(nsurf);
    cells[icell].nsurf = 0;
  }

  // error if surf count exceeds maxsurfpercell in any cell

  MPI_Allreduce(&max,&maxall,1,MPI_INT,MPI_MAX,world);

  if (maxall > maxsurfpercell) {
    if (me == 0) printf("Max surfs in any cell = %d\n",maxall);
    error->all(FLERR,"Too many surfs in one cell - set global surfmax");
  }

  // non-distributed surfs:
  // store surfIDs as local indices in each cell's csurfs list
  // for performance, sort each cell's csurfs list, same order as cell alg

  if (!distributed) {
    for (m = 0; m < nreturn; m++) {
      icell = (*hash)[outbuf[m].cellID] - 1;
      cells[icell].csurfs[cells[icell].nsurf++] = outbuf[m].surfID - 1;
    }

    for (icell = 0; icell < nlocal; icell++)
      if (cells[icell].nsurf) {
        cinfo[icell].type = OVERLAP;
        qsort(cells[icell].csurfs,cells[icell].nsurf,
              sizeof(surfint),compare_surfIDs);
      }

    memory->sfree(outbuf);
  }

  if (outflag) {
    MPI_Barrier(world);
    t3 = MPI_Wtime();
    tmap = t2-t1;
    trvous1 = t3-t2;
  }

  // stage 3 ----------------------------------------------------
  // distributed surfs
  // rendezvous operation to obtain nlocal surfs for each proc
  // each grid cell requests a surf from proc that owns surf in mylines/mytris
  //   use hash to only do this once per surf
  // receive the surf and store in nlocal lines/tris

  if (distributed) {

    // hash for surf IDs that are in my owned cells

    MySurfHash shash;
    MyIterator it;

    // loop over output of previous rvous operation
    // ninput = # of unique surfs I need for my owned grid cells
    // store IDs of those surfs in hash

    ncount = 0;
    for (m = 0; m < nreturn; m++) {
      if (shash.find(outbuf[m].surfID) == shash.end()) {
        shash[outbuf[m].surfID] = 0;
        ncount++;
      }
    }

    // allocate memory for rvous input

    int *proclist2;
    memory->create(proclist2,ncount,"surf2grid:proclist2");
    InRvous2 *inbuf2 = 
      (InRvous2 *) memory->smalloc((bigint) ncount*sizeof(InRvous2),
                                   "surf2grid:inbuf2");

    // create rvous inputs
    // proclist2 = owner of each surf

    ncount = 0;
    for (it = shash.begin(); it != shash.end(); ++it) {
      surfID = it->first;
      proclist2[ncount] = (surfID-1) % nprocs;
      inbuf2[ncount].proc = me;
      inbuf2[ncount].surfID = surfID;
      ncount++;
    }

    // perform rendezvous operation
    // each proc owns subset of surfs
    // receives all surf requests to return surf to each proc who needs it
    
    char *outbuf2;
    int outbytes;
    if (dim == 2) outbytes = sizeof(OutRvous2line);
    else outbytes = sizeof(OutRvous2tri);

    int nreturn2 = comm->rendezvous(1,ncount,(char *) inbuf2,sizeof(InRvous2),
                                    0,proclist2,rendezvous_surfrequest,
                                    0,outbuf2,outbytes,(void *) this);
    
    memory->destroy(proclist2);
    memory->sfree(inbuf2);

    // copy entire rendezvous output buf into realloced Surf lines/tris

    surf->nlocal = surf->nghost = 0;
    surf->nmax = surf->nlocal = nreturn2;
    surf->grow();

    if (dim == 2) memcpy(surf->lines,outbuf2,nreturn2*sizeof(Surf::Line));
    else memcpy(surf->tris,outbuf2,nreturn2*sizeof(Surf::Tri));

    memory->sfree(outbuf2);

    // reset Surf hash to point to surf list in lines/tris

    if (dim == 2) {
      Surf::Line *lines = surf->lines;
      for (i = 0; i < nreturn2; i++) {
        surfID = lines[i].id;
        shash[surfID] = i;
      }
    } else {
      Surf::Tri *tris = surf->tris;
      for (i = 0; i < nreturn2; i++) {
        surfID = tris[i].id;
        shash[surfID] = i;
      }
    }

    // loop over first rendezvous outbuf
    // store surfIDs as local indices in each cell's csurfs list
    // for performance, sort each cell's csurfs list, same order as cell alg
    
    for (m = 0; m < nreturn; m++) {
      icell = (*hash)[outbuf[m].cellID] - 1;
      isurf = shash[outbuf[m].surfID];
      cells[icell].csurfs[cells[icell].nsurf++] = isurf;
    }

    for (icell = 0; icell < nlocal; icell++)
      if (cells[icell].nsurf) {
        cinfo[icell].type = OVERLAP;
        qsort(cells[icell].csurfs,cells[icell].nsurf,
              sizeof(surfint),compare_surfIDs);
      }

    memory->sfree(outbuf);
  }

  if (outflag) {
    MPI_Barrier(world);
    t4 = MPI_Wtime();
    trvous2 = t4-t3;
  }

  if (outflag) surf2grid_stats();
}

/* ----------------------------------------------------------------------
   compute split cells for implicit surfs
   surfs per cell already created
   called from ReadISurf
------------------------------------------------------------------------- */

void Grid::surf2grid_implicit(int subflag, int outflag)
{
  int dim = domain->dimension;
  if (dim == 3) cut3d = new Cut3d(sparta);
  else cut2d = new Cut2d(sparta,domain->axisymmetric);

  tmap = trvous1 = trvous2 = 0.0;

  if (outflag) surf2grid_stats();
  surf2grid_split(subflag,outflag);
}

/* ----------------------------------------------------------------------
   compute cut volume of each cell and any split cell info
   nsurf and csurfs list for each grid cell have already been computed
   if subflag = 1, create new owned split and sub cells as needed
     called from ReadSurf, RemoveSurf, MoveSurf
   if subflag = 0, split/sub cells already exist
     called from ReadRestart
   in cells: set nsplit, isplit
   in cinfo: set corner, volume 
   initialize sinfo as needed
------------------------------------------------------------------------- */

void Grid::surf2grid_split(int subflag, int outflag)
{
  int i,isub,nsurf,nsplitone,xsub;
  int *surfmap,*ptr;
  double t1,t2;
  double *lo,*hi,*vols;
  double xsplit[3];
  ChildCell *c;
  SplitInfo *s;

  int dim = domain->dimension;

  if (outflag) {
    MPI_Barrier(world);
    t1 = MPI_Wtime();
  }

  // compute cut volume and possible split of each grid cell by surfs
  // decrement nunsplitlocal if convert an unsplit cell to split cell
  // if nsplitone > 1, create new split cell sinfo and sub-cells
  // skip if nsplit <= 0 b/c split cells could exist if restarting

  int max = 0;
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
      if (nsplitone > maxsplitpercell) {
	max = MAX(max,nsplitone);
	csplits->vgot(0);
	
      } else {
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
	
	csubs->vgot(nsplitone);
	csplits->vgot(cells[icell].nsurf);
      }
      
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

  // error if split count exceeds maxsplitpercell for any cell

  int maxall;
  MPI_Allreduce(&max,&maxall,1,MPI_INT,MPI_MAX,world);
  if (maxall) {
    if (me == 0) printf("Max split cells in any cell = %d\n",maxall);
    error->all(FLERR,"Too many split cells in a single cell - "
               "set global splitmax");
  }

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

  // clean up

  if (dim == 3) delete cut3d;
  else delete cut2d;

  if (outflag) {
    MPI_Barrier(world);
    t2 = MPI_Wtime();
    tsplit = t2-t1;
  }
}

/* ----------------------------------------------------------------------
   map surf elements into a single grid cell = icell
   flag = 0 for grid refinement, 1 for grid coarsening
   in cells: set nsurf, csurfs, nsplit, isplit
   in cinfo: set type, corner, volume
   initialize sinfo as needed
   called from AdaptGrid
------------------------------------------------------------------------- */

void Grid::surf2grid_one(int flag, int icell, int iparent, int nsurf_caller,
                         Cut3d *cut3d, Cut2d *cut2d)
{
  int nsurf,isub,xsub,nsplitone;
  int *iptr;
  surfint *sptr;
  double xsplit[3];
  double *vols;

  int dim = domain->dimension;

  // identify surfs in new cell only for grid refinement

  if (flag == 0) {
    sptr = csurfs->vget();
    if (dim == 3)
      nsurf = cut3d->surf2grid_list(cells[icell].id,
                                    cells[icell].lo,cells[icell].hi,
                                    cells[iparent].nsurf,cells[iparent].csurfs,
                                    sptr,maxsurfpercell);
    else
      nsurf = cut2d->surf2grid_list(cells[icell].id,
                                    cells[icell].lo,cells[icell].hi,
                                    cells[iparent].nsurf,cells[iparent].csurfs,
                                    sptr,maxsurfpercell);

    if (nsurf == 0) return;
    if (nsurf > maxsurfpercell) {
      printf("Surfs in one refined cell = %d\n",nsurf);
      error->one(FLERR,"Too many surfs in one refined cell - set global surfmax");
    }
    
    cinfo[icell].type = OVERLAP;
    cells[icell].nsurf = nsurf;
    cells[icell].csurfs = sptr;
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
    
    iptr = s->csubs = csubs->vget();
    
    for (int i = 0; i < nsplitone; i++) {
      isub = nlocal;
      add_sub_cell(icell,1);
      cells[isub].nsplit = -i;
      cinfo[isub].volume = vols[i];
      iptr[i] = isub;
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
   re-allocate page data structs to hold variable-length surf and cell info
   used for mapping surfaces to grid cells
------------------------------------------------------------------------- */

void Grid::allocate_surf_arrays()
{
  delete csurfs;
  delete csplits;
  delete csubs;

  csurfs = new MyPage<surfint>(maxsurfpercell,MAX(100*maxsurfpercell,1024));
  csplits = new MyPage<int>(maxsurfpercell,MAX(100*maxsurfpercell,1024));
  csubs = new MyPage<int>(maxsplitpercell,MAX(100*maxsplitpercell,128));
}

void Grid::allocate_cell_arrays()
{
  delete cpsurf;
  cpsurf = new MyPage<cellint>(maxcellpersurf,MAX(100*maxcellpersurf,1024));
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

// ----------------------------------------------------------------------
// private class methods
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   rendezvous decomposition computation
   create list of surfs for each grid cell
   recv (surfID,cellID) pairs from owner of each surf
     also (me,cellID) from owner of each grid cell
   send (cellID,surfID) pairs back to owning proc of each grid cell
   do this in 2 passes to count sizes for each grid cell
---------------------------------------------------------------------- */

int Grid::rendezvous_surflist(int n, char *inbuf, int &flag,
                              int *&proclist, char *&outbuf,
                              void *ptr)
{
  int i,m;

  Grid *gptr = (Grid *) ptr;
  Memory *memory = gptr->memory;

  // initialize hash
  // ncount = number of atoms assigned to me
  // key = atom ID
  // value = index into Ncount-length data structure

  InRvous *in = (InRvous *) inbuf;
  MyHash rhash;
  
  int ncount = 0;
  for (i = 0; i < n; i++)
    if (in[i].proc >= 0)
      rhash[in[i].cellID] = ncount++;

  // procowner = caller proc that owns each cell
  // surfcount = count of surfs per cell

  int *procowner,*surfcount;
  memory->create(procowner,ncount,"gridsurf:procowner");
  memory->create(surfcount,ncount,"gridsurf:surfcount");
  for (m = 0; m < ncount; m++) surfcount[m] = 0;

  for (i = 0; i < n; i++) { 
    m = rhash.find(in[i].cellID)->second;
    if (in[i].proc >= 0) procowner[m] = in[i].proc;
    else surfcount[m]++;
  }

  // pass list of OutRvous datums back to comm->rendezvous

  int nout = 0;
  for (m = 0; m < ncount; m++) nout += surfcount[m];

  memory->create(proclist,nout,"gridsurf:proclist");
  OutRvous *out = (OutRvous *)
    memory->smalloc((bigint) nout*sizeof(OutRvous),"gridsurf:out");

  nout = 0;
  for (i = 0; i < n; i++) {
    if (in[i].proc >= 0) continue;
    m = rhash.find(in[i].cellID)->second;
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

  // flag = 2: new outbuf

  flag = 2;
  return nout;
}

/* ----------------------------------------------------------------------
   rendezvous decomposition computation
   return surf info for each requested surf element
   recv (proc,surfID) pairs from requestor of each surf
   send memcpy of surf info back to the requesting proc
---------------------------------------------------------------------- */

int Grid::rendezvous_surfrequest(int n, char *inbuf, int &flag,
                                 int *&proclist, char *&outbuf,
                                 void *ptr)
{
  int i,m;
  
  Grid *gptr = (Grid *) ptr;
  Surf::Line *mylines = gptr->surf->mylines;
  Surf::Tri *mytris = gptr->surf->mytris;
  Memory *memory = gptr->memory;
  int dim = gptr->domain->dimension;
  int nprocs = gptr->comm->nprocs;

  InRvous2 *in = (InRvous2 *) inbuf;

  // proclist = list of procs who made requests
  // outbuf = list of lines or tris to pass back to comm->rendezvous

  memory->create(proclist,n,"gridsurf:proclist");

  int surfbytes;
  if (dim == 2) surfbytes = sizeof(OutRvous2line);
  else surfbytes = sizeof(OutRvous2tri);
  outbuf = (char *) memory->smalloc((bigint) n*surfbytes,"gridsurf:outbuf");

  if (dim == 2) {
    for (i = 0; i < n; i++) {
      proclist[i] = in[i].proc;
      m = (in[i].surfID-1) / nprocs;
      memcpy(&outbuf[(bigint) i*surfbytes],&mylines[m],surfbytes);
    }
  } else {
    for (i = 0; i < n; i++) {
      proclist[i] = in[i].proc;
      m = (in[i].surfID-1) / nprocs;
      memcpy(&outbuf[(bigint) i*surfbytes],&mytris[m],surfbytes);
    }
  }

  // Comm::rendezvous will delete proclist and outbuf
  // flag = 2: new outbuf

  flag = 2;
  return n;
}

/* ----------------------------------------------------------------------
   enumerate intersections of isurf with any child grid cell
   recurse 2d/3d do this by recursively walking down the grid tree
   called by surf2grid_surf_algorithm()
------------------------------------------------------------------------- */

int Grid::find_overlaps(int isurf, cellint *list)
{
  double slo[3],shi[3];

  int ncell = 0;

  if (domain->dimension == 2) {
    Surf::Line *line;
    if (surf->distributed) line = &surf->mylines[isurf];
    else line = &surf->lines[isurf];

    slo[0] = MIN(line->p1[0],line->p2[0]);
    slo[1] = MIN(line->p1[1],line->p2[1]);
    shi[0] = MAX(line->p1[0],line->p2[0]);
    shi[1] = MAX(line->p1[1],line->p2[1]);
    
    recurse2d(isurf,slo,shi,0,ncell,list);

  } else {
    Surf::Tri *tri;
    if (surf->distributed) tri = &surf->mytris[isurf];
    else tri = &surf->tris[isurf];

    double sbox[2][3];
    slo[0] = MIN(tri->p1[0],tri->p2[0]);
    slo[0] = MIN(tri->p3[0],slo[0]);
    slo[1] = MIN(tri->p1[1],tri->p2[1]);
    slo[1] = MIN(tri->p3[1],slo[1]);
    slo[2] = MIN(tri->p1[2],tri->p2[2]);
    slo[2] = MIN(tri->p3[2],slo[2]);

    shi[0] = MAX(tri->p1[0],tri->p2[0]);
    shi[0] = MAX(tri->p3[0],shi[0]);
    shi[1] = MAX(tri->p1[1],tri->p2[1]); 
    shi[1] = MAX(tri->p3[1],shi[1]);
    shi[2] = MAX(tri->p1[2],tri->p2[2]);
    shi[2] = MAX(tri->p3[2],shi[2]);
    
    recurse3d(isurf,slo,shi,0,ncell,list);
  }

  return ncell;
}

/* ----------------------------------------------------------------------
   enumerate intersections of isurf with any child grid cell in iparent cell
   called by surf2grid_surf_algorithm() via find_overlaps()
   iline = which surf element
   slo/shi = bounding box around line segment
   iparent = current parent cell
   n, list = growing list of cell IDs this line overlaps with
   cut2d->surf2grid_one() is used to determine actual overlap
------------------------------------------------------------------------- */

void Grid::recurse2d(int iline, double *slo, double *shi, int iparent, 
                     int &n, cellint *list)
{
  int ix,iy,newparent,index,parentflag,overlap;
  cellint ichild,idchild;
  double celledge;
  double newslo[2],newshi[2];
  double clo[3],chi[3];

  double *p1,*p2;
  if (surf->distributed) {
    p1 = surf->mylines[iline].p1;
    p2 = surf->mylines[iline].p2;
  } else {
    p1 = surf->lines[iline].p1;
    p2 = surf->lines[iline].p2;
  }

  ParentCell *p = &pcells[iparent];
  double *plo = p->lo;
  double *phi = p->hi;
  int nx = p->nx;
  int ny = p->ny;

  // ilo,ihi jlo,jhi = indices for range of grid cells overlapped by surf bbox
  // overlap means point is inside grid cell or touches boundary
  // same equation as in Grid::id_find_child()

  int ilo = static_cast<int> ((slo[0]-plo[0]) * nx/(phi[0]-plo[0]));
  int ihi = static_cast<int> ((shi[0]-plo[0]) * nx/(phi[0]-plo[0]));
  int jlo = static_cast<int> ((slo[1]-plo[1]) * ny/(phi[1]-plo[1]));
  int jhi = static_cast<int> ((shi[1]-plo[1]) * ny/(phi[1]-plo[1]));

  // augment indices if sbox touches or slightly overlaps cell edges
  // same equation as in Grid::id_child_lohi()

  celledge = plo[0] + ilo*(phi[0]-plo[0])/nx;
  if (slo[0] <= celledge) ilo--;
  celledge = plo[0] + (ihi+1)*(phi[0]-plo[0])/nx;
  if (shi[0] >= celledge) ihi++;

  celledge = plo[1] + jlo*(phi[1]-plo[1])/ny;
  if (slo[1] <= celledge) jlo--;
  celledge = plo[1] + (jhi+1)*(phi[1]-plo[1])/ny;
  if (shi[1] >= celledge) jhi++;

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
  // if cell is a child cell, compute overlap via surf2grid_one()
  // else it is a parent cell, so recurse
  // set newslo/newshi to intersection of slo/shi with new parent cell

  for (iy = jlo; iy <= jhi; iy++) {
    for (ix = ilo; ix <= ihi; ix++) {
      ichild = (cellint) iy*nx + ix + 1;
      idchild = p->id | (ichild << p->nbits);
      grid->id_child_lohi(iparent,ichild,clo,chi);

      if (hash->find(idchild) == hash->end()) parentflag = 0;
      else if ((*hash)[idchild] >= 0) parentflag = 0;
      else parentflag = 1;
      
      if (parentflag) {
        index = (*hash)[idchild];
        newparent = -index-1;
        newslo[0] = MAX(slo[0],clo[0]);
        newslo[1] = MAX(slo[1],clo[1]);
        newshi[0] = MIN(shi[0],chi[0]);
        newshi[1] = MIN(shi[1],chi[1]);
        recurse2d(iline,newslo,newshi,newparent,n,list);
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
   called by surf2grid_surf_algorithm() via find_overlaps()
   itri = which surf element
   slo/shi = bounding box around triangle
   iparent = current parent cell
   n, list = growing list of cell IDs this tri overlaps with
   cut2d->surf2grid_one() is used to determine actual overlap
------------------------------------------------------------------------- */

void Grid::recurse3d(int itri, double *slo, double *shi, int iparent, 
                     int &n, cellint *list)
{
  int ix,iy,iz,newparent,index,parentflag,overlap;
  cellint ichild,idchild;
  double celledge;
  double newslo[3],newshi[3];
  double clo[3],chi[3];

  double *p1,*p2,*p3;
  if (surf->distributed) {
    p1 = surf->mytris[itri].p1;
    p2 = surf->mytris[itri].p2;
    p3 = surf->mytris[itri].p3;
  } else {
    p1 = surf->tris[itri].p1;
    p2 = surf->tris[itri].p2;
    p3 = surf->tris[itri].p3;
  }

  ParentCell *p = &pcells[iparent];
  double *plo = p->lo;
  double *phi = p->hi;
  int nx = p->nx;
  int ny = p->ny;
  int nz = p->nz;

  // ilo,ihi jlo,jhi = indices for range of grid cells overlapped by surf bbox
  // overlap means point is inside grid cell or touches boundary
  // same equation as in Grid::id_find_child()

  int ilo = static_cast<int> ((slo[0]-plo[0]) * nx/(phi[0]-plo[0]));
  int ihi = static_cast<int> ((shi[0]-plo[0]) * nx/(phi[0]-plo[0]));
  int jlo = static_cast<int> ((slo[1]-plo[1]) * ny/(phi[1]-plo[1]));
  int jhi = static_cast<int> ((shi[1]-plo[1]) * ny/(phi[1]-plo[1]));
  int klo = static_cast<int> ((slo[2]-plo[2]) * nz/(phi[2]-plo[2]));
  int khi = static_cast<int> ((shi[2]-plo[2]) * nz/(phi[2]-plo[2]));
  
  // augment indices if sbox touches or slightly overlaps cell edges
  // same equation as in Grid::id_child_lohi()

  celledge = plo[0] + ilo*(phi[0]-plo[0])/nx;
  if (slo[0] <= celledge) ilo--;
  celledge = plo[0] + (ihi+1)*(phi[0]-plo[0])/nx;
  if (shi[0] >= celledge) ihi++;

  celledge = plo[1] + jlo*(phi[1]-plo[1])/ny;
  if (slo[1] <= celledge) jlo--;
  celledge = plo[1] + (jhi+1)*(phi[1]-plo[1])/ny;
  if (shi[1] >= celledge) jhi++;

  celledge = plo[2] + klo*(phi[2]-plo[2])/nz;
  if (slo[2] <= celledge) klo--;
  celledge = plo[2] + (khi+1)*(phi[2]-plo[2])/nz;
  if (shi[2] >= celledge) khi++;

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
  // set newslo/newshi to intersection of slo/shi with new parent cell

  for (iz = klo; iz <= khi; iz++) {
    for (iy = jlo; iy <= jhi; iy++) {
      for (ix = ilo; ix <= ihi; ix++) {
        ichild = (cellint) iz*nx*ny + (cellint) iy*nx + ix + 1;
        idchild = p->id | (ichild << p->nbits);
        grid->id_child_lohi(iparent,ichild,clo,chi);

        if (hash->find(idchild) == hash->end()) parentflag = 0;
        else if ((*hash)[idchild] >= 0) parentflag = 0;
        else parentflag = 1;
      
        if (parentflag) {
          index = (*hash)[idchild];
          newparent = -index-1;
          newslo[0] = MAX(slo[0],clo[0]);
          newslo[1] = MAX(slo[1],clo[1]);
          newslo[2] = MAX(slo[2],clo[2]);
          newshi[0] = MIN(shi[0],chi[0]);
          newshi[1] = MIN(shi[1],chi[1]);
          newshi[2] = MIN(shi[2],chi[2]);
          recurse3d(itri,newslo,newshi,newparent,n,list);
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
   comparison function invoked by qsort() called by surf2grid_surf_algorithm()
   used to sort the csurfs list of a single cell
   this is not a class method
------------------------------------------------------------------------- */

int compare_surfIDs(const void *iptr, const void *jptr)
{
  int i = *((surfint *) iptr);
  int j = *((surfint *) jptr);
  if (i < j) return -1;
  if (i > j) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   output stats on overlap of surfs with grid cells
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

        // NOTE: this line is bad for distributed surfs
        //       b/c calling line_size that will access lines
        //       which is not yet initialized
        // maybe should do rendezvous already to get them?

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

/* ----------------------------------------------------------------------
   output stats on aggregate flow volume and cells in/out of the flow
------------------------------------------------------------------------- */

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
   compute flow volume for entire box, using list of surfs
   volume for one surf is projection to lower z face
   NOTE: this does not work if any surfs are clipped to zlo or zhi faces in 3d
         this does not work if any surfs are clipped to ylo or yhi faces in 3d
         need to add contribution due to closing surfs on those faces
         fairly easy to add in 2d, not so easy in 3d
------------------------------------------------------------------------- */

double Grid::flow_volume()
{
  double zarea,volall;
  double *p1,*p2,*p3;

  int n;
  Surf::Line *lines;
  Surf::Tri *tris;

  if (domain->dimension == 3) {
    if (!surf->implicit && surf->distributed) {
      tris = surf->mytris;
      n = surf->nown;
    } else {
      tris = surf->tris;
      n = surf->nlocal;
    }
  } else {
    if (!surf->implicit && surf->distributed) {
      lines = surf->mylines;
      n = surf->nown;
    } else {
      lines = surf->lines;
      n = surf->nlocal;
    }
  }

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  double volume = 0.0;

  if (domain->dimension == 3) {
    for (int i = 0; i < n; i++) {
      p1 = tris[i].p1;
      p2 = tris[i].p2;
      p3 = tris[i].p3;
      zarea = 0.5 * ((p2[0]-p1[0])*(p3[1]-p1[1]) - (p2[1]-p1[1])*(p3[0]-p1[0]));
      volume -= zarea * ((p1[2]+p2[2]+p3[2])/3.0 - boxlo[2]);
    }

    if (surf->distributed)
      MPI_Allreduce(&volume,&volall,1,MPI_DOUBLE,MPI_SUM,world);
    else volall = volume;

    if (volall <= 0.0) 
      volall += (boxhi[0]-boxlo[0]) * (boxhi[1]-boxlo[1]) * 
        (boxhi[2]-boxlo[2]); 
 
  // axisymmetric "volume" of line segment = volume of truncated cone
  // PI/3 (y1^2 + y1y2 + y2^2) (x2-x1)

  } else if (domain->axisymmetric) {
    for (int i = 0; i < n; i++) {
      p1 = lines[i].p1;
      p2 = lines[i].p2;
      volume -= 
        MY_PI3 * (p1[1]*p1[1] + p1[1]*p2[1] + p2[1]*p2[1]) * (p2[0]-p1[0]);
    }

    if (surf->distributed)
      MPI_Allreduce(&volume,&volall,1,MPI_DOUBLE,MPI_SUM,world);
    else volall = volume;

    if (volall <= 0.0) 
      volall += MY_PI * boxhi[1]*boxhi[1] * (boxhi[0]-boxlo[0]);

  } else {
    for (int i = 0; i < n; i++) {
      p1 = lines[i].p1;
      p2 = lines[i].p2;
      volume -= (0.5*(p1[1]+p2[1]) - boxlo[1]) * (p2[0]-p1[0]);
    }

    if (surf->distributed)
      MPI_Allreduce(&volume,&volall,1,MPI_DOUBLE,MPI_SUM,world);
    else volall = volume;

    if (volall <= 0.0) volall += (boxhi[0]-boxlo[0]) * (boxhi[1]-boxlo[1]); 
  }

  return volall;
}

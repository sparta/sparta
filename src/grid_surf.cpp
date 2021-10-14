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
#include "modify.h"
#include "collide.h"
#include "surf.h"
#include "cut2d.h"
#include "cut3d.h"
#include "irregular.h"
#include "geometry.h"
#include "math_const.h"
#include "hashlittle.h"
#include "my_page.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

// prototype for non-class function

int compare_surfIDs(const void *, const void *);

#define BIG 1.0e20
#define CHUNK 16
#define EPSSURF 1.0e-4
#define DELTA_SEND 16384

enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};         // several files
enum{PERAUTO,PERCELL,PERSURF};                // several files
enum{SOUTSIDE,SINSIDE,ONSURF2OUT,ONSURF2IN};  // several files (changed 2 words)

// operations for surfaces in grid cells

/* ----------------------------------------------------------------------
   map surf elements to grid cells for explicit surfs (distributed or not)
   via one of two algorithms
   cell_alg = original, loop over my cells, check all surfs within bbox
   surf_alg = Jan19, loop over N/P surfs, find small set of cells each overlaps,
              perform rendezvous comm to convert cells per surf to surfs per cell
   new_alg = Nov20, no more parent cells, rendezvous alg at each level of grid
   new2_alg = Mar21, modified new alg, use recursive tree to speed overlap finding
     this one is now surf_alg
   for distributed surfs, have to use surf alg
   PERAUTO option chooses based on total nsurfs vs nprocs
   see info on subflag, outflag options with surf2grid_split()
   called from Readsurf, MoveSurf, RemoveSurf, ReadRestart, and FixMoveSurf
------------------------------------------------------------------------- */

void Grid::surf2grid(int subflag, int outflag)
{
  if (surf->distributed) {
    surf2grid_surf_algorithm(outflag);
  } else if (surfgrid_algorithm == PERAUTO) {
    if (comm->nprocs > surf->nsurf) surf2grid_cell_algorithm(outflag);
    else {
      surf2grid_surf_algorithm(outflag);
    }
  } else if (surfgrid_algorithm == PERCELL) {
    surf2grid_cell_algorithm(outflag);
  } else if (surfgrid_algorithm == PERSURF) {
    surf2grid_surf_algorithm(outflag);
  }

  // now have nsurf,csurfs list of local surfs that overlap each cell
  // compute cut volume and split info for each cell

  surf2grid_split(subflag,outflag);
}

/* ----------------------------------------------------------------------
   compute split cells for implicit surfs
   surfs per cell already created
   called from ReadISurf and FixAblate
------------------------------------------------------------------------- */

void Grid::surf2grid_implicit(int subflag, int outflag)
{
  int dim = domain->dimension;
  if (dim == 3) cut3d = new Cut3d(sparta);
  else cut2d = new Cut2d(sparta,domain->axisymmetric);

  tmap = tcomm1 = tcomm2 = tcomm3 = 0.0;

  if (outflag) surf2grid_stats();
  surf2grid_split(subflag,outflag);
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

    // add nsplitone sub cells
    // collide and fixes also need to add cells

    for (int i = 0; i < nsplitone; i++) {
      isub = nlocal;
      add_sub_cell(icell,1);
      if (collide) collide->add_grid_one();
      if (modify->n_pergrid) modify->add_grid_one();
      cells[isub].nsplit = -i;
      cinfo[isub].volume = vols[i];
      iptr[i] = isub;
    }

    csplits->vgot(nsurf);
    csubs->vgot(nsplitone);
  }
}

/* ----------------------------------------------------------------------
   find surfs that overlap owned grid cells, only for non-distributed surfs
   algorithm: for each of my cells, check all surfs
   in cells: set nsurf, csurfs
   in cinfo: set type=OVERLAP for cells with surfs
------------------------------------------------------------------------- */

void Grid::surf2grid_cell_algorithm(int outflag)
{
  int i,nsurf,nontrans;
  double t1,t2;
  surfint *ptr;
  double *lo,*hi;

  int dim = domain->dimension;

  if (outflag) {
    MPI_Barrier(world);
    t1 = MPI_Wtime();
  }

  surf->bbox_all();

  double *slo = surf->bblo;
  double *shi = surf->bbhi;

  if (dim == 3) cut3d = new Cut3d(sparta);
  else cut2d = new Cut2d(sparta,domain->axisymmetric);

  // compute overlap of all surfs with each cell I own
  // info stored in nsurf,csurfs
  // skip if nsplit <= 0 b/c split cells could exist if restarting

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;

  int max = 0;

  for (int icell = 0; icell < nlocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;

    // skip grid cell if outside bounding box of all surfs

    lo = cells[icell].lo;
    hi = cells[icell].hi;
    if (!box_overlap(lo,hi,slo,shi)) continue;

    // cut2d/3d surf2grid finds intersection of all surfs with a single grid cell

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

      // only mark cell as OVERLAP if has a non-transparent surf element

      nontrans = 0;

      if (dim == 2) {

        for (i = 0; i < nsurf; i++) {
          if (!lines[ptr[i]].transparent) nontrans = 1;
          break;
        }
      } else {
        for (i = 0; i < nsurf; i++) {
          if (!tris[ptr[i]].transparent) nontrans = 1;
          break;
        }
      }

      if (nontrans) cinfo[icell].type = OVERLAP;
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
    tcomm1 = tcomm2 = tcomm3 = tcomm4 = 0.0;
  }

  if (outflag) surf2grid_stats();
}

/* ----------------------------------------------------------------------
   find surfs that overlap owned grid cells, for non-distributed or distributed
   algorithm:
     each proc responsible for subset of surfs
     loop over levels of hierarchical grid
     conceptual uniform grid at that level overlayed on bounding box for all surfs
     partition uniform grid via RCB into sub-grids per proc
     irregular comm of 2 sets of data to the RCB procs
       (a) surfs I own with a bbox that overlaps any procs RCB sub-grid
       (b) child cells I own at level to owning proc in RCB decomp
     each proc can then identify surf/grid intersections for its RCB grid cells
       done recursively by dropping each surf down conceptual tree of parent cells
       until reach child cells that exist at level
     irregular comm the surf/grid intersection pairs back to procs that own grid cells
   in cells: set nsurf, csurfs
   in cinfo: set type=OVERLAP for cells with surfs
------------------------------------------------------------------------- */

void Grid::surf2grid_surf_algorithm(int outflag)
{
  int i,j,n,ix,iy,iz,icell,isurf;
  cellint childID,parentID;
  double t1,t2,t3,t4,t5;
  double ctr[3];
  Irregular *irregular;

  double *boxlo,*boxhi;             // corner points of entire simulation box
  double *allsurflo,*allsurfhi;     // bounding box for all surfs in system
  int unilo[3],unihi[3];            // indices of corner cells of uniform grid
                                    //   at level which encompasses allsurf_lo/hi
  int myunilo[3],myunihi[3];        // my RCB portion of unilo/hi, inclusive
  double slo[3],shi[3];             // bounding box around one surf
  int sunilo[3],sunihi[3];          // indices of corner cells of uniform grid
                                    //   at level which encompasses a single surf
  int glo[3],ghi[3];                // corner indices of a grid box
  double bblo[3],bbhi[3];           // corners of a bounding box
  double rcblo[3],rcbhi[3];         // corners of my RCB box
  GridTree *gtree;                  // tree of RCB cuts for partitioning uniform grid

  // data structs for communication

  struct Send2 {
    cellint childID;
    int proc,icell;
  };

  struct Send3 {
    surfint surfID;
    int icell;
  };

  struct RCBlohi {
    double lo[3],hi[3];             // corner pts for RCB child cells
  };

  int me = comm->me;
  int nprocs = comm->nprocs;
  int dim = domain->dimension;
  int distributed = surf->distributed;

  boxlo = domain->boxlo;
  boxhi = domain->boxhi;

  surf->bbox_all();                 // bounding box for all surfs
  allsurflo = surf->bblo;
  allsurfhi = surf->bbhi;

  if (dim == 3) cut3d = new Cut3d(sparta);
  else cut2d = new Cut2d(sparta,domain->axisymmetric);

  int *plist;
  memory->create(plist,nprocs,"surf2grid:plist");
  gtree = (GridTree *) memory->smalloc(nprocs*sizeof(GridTree),"surf2grid:gtree");

  // data structs for 3 rendezvous comms

  int *proclist1 = NULL;
  int *proclist2 = NULL;
  int *proclist3 = NULL;
  char *sbuf1 = NULL;
  Send2 *sbuf2 = NULL;
  Send3 *sbuf3 = NULL;
  int maxsend1 = 0;
  int maxsend2 = 0;
  int maxsend3 = 0;
  int **pairs = NULL;
  int maxpair = 0;

  // which set of Lines or Tris to process, distributed or not

  Surf::Line *lines;
  Surf::Tri *tris;
  int nsurf,istart,istop,idelta,nbytes_surf;

  if (distributed) {
    lines = surf->mylines;
    tris = surf->mytris;
    nsurf = surf->nown;
    istart = 0;
    istop = nsurf;
    idelta = 1;
  } else {
    lines = surf->lines;
    tris = surf->tris;
    int ntotal = surf->nsurf;
    nsurf = ntotal / nprocs;
    if (me < ntotal % nprocs) nsurf++;
    istart = comm->me;
    istop = ntotal;
    idelta = nprocs;
  }

  if (dim == 2) nbytes_surf = sizeof(Surf::Line);
  else nbytes_surf = sizeof(Surf::Tri);

  // loop over levels of grid
  // at each iteration, operate only on child cells that exist at that level

  tmap = tcomm1 = tcomm2 = tcomm3 = 0.0;

  int minlevel = set_minlevel();

  for (int level = minlevel; level <= maxlevel; level++) {

    if (outflag) {
      MPI_Barrier(world);
      t1 = MPI_Wtime();
    }

    // compute extent of uniform grid at level which overlaps surf bbox
    // unilo/hi = inclusive range of grid box of overlapping grid box

    id_find_child_uniform_level(level,0,boxlo,boxhi,allsurflo,
				unilo[0],unilo[1],unilo[2]);
    id_find_child_uniform_level(level,1,boxlo,boxhi,allsurfhi,
				unihi[0],unihi[1],unihi[2]);

    // compute a recursive decomp (RCB) of the uniform grid box
    // gtree = tree of RCB cuts, cuts are along grid planes
    // myunilo/hi = inclusive range of my portion of grid box
    // rcblo/hi = corner points of my RCB box

    partition_grid(0,nprocs-1,unilo[0],unihi[0],unilo[1],unihi[1],
		   unilo[2],unihi[2],gtree);
    myunilo[0] = unilo[0]; myunihi[0] = unihi[0];
    myunilo[1] = unilo[1]; myunihi[1] = unihi[1];
    myunilo[2] = unilo[2]; myunihi[2] = unihi[2];
    mybox(me,0,nprocs-1,myunilo[0],myunihi[0],myunilo[1],
	  myunihi[1],myunilo[2],myunihi[2],gtree);

    childID = id_uniform_level(level,myunilo[0],myunilo[1],myunilo[2]);
    id_lohi(childID,level,boxlo,boxhi,rcblo,bbhi);
    childID = id_uniform_level(level,myunihi[0],myunihi[1],myunihi[2]);
    id_lohi(childID,level,boxlo,boxhi,bblo,rcbhi);

    // first irregular comm
    // loop over my surfs:
    //   compute single surf bbox as a brick of uniform grid cells at this level
    //   drop bbox down RCB tree to identify set of RCB procs the surf overlaps
    // send copy of surf geometry to RCB procs
    // nrecv1 = # of surfs I have copy of in RCB decomp
    // NOTE: this comm might be faster in Rvous mode?

    int nsend = 0;

    for (isurf = istart; isurf < istop; isurf += idelta) {
      if (dim == 2) surf->bbox_one(&lines[isurf],slo,shi);
      else surf->bbox_one(&tris[isurf],slo,shi);
      id_find_child_uniform_level(level,0,boxlo,boxhi,slo,
				  sunilo[0],sunilo[1],sunilo[2]);
      id_find_child_uniform_level(level,1,boxlo,boxhi,shi,
				  sunihi[0],sunihi[1],sunihi[2]);

      // drop trimmed surf box on RCB tree
      // return list of procs whose RCB subbox it overlaps

      int np = 0;
      box_drop(sunilo,sunihi,0,nprocs-1,gtree,np,plist);
      if (!np) continue;

      for (i = 0; i < np; i++) {
	if (nsend == maxsend1) {
	  maxsend1 += DELTA_SEND;
	  memory->grow(proclist1,maxsend1,"surf2grid:proclist1");
	  if (dim == 2)
	    sbuf1 = (char *) memory->srealloc(sbuf1,maxsend1*sizeof(Surf::Line),
					      "surf2grid:sbuf1");
	  else
	    sbuf1 = (char *) memory->srealloc(sbuf1,maxsend1*sizeof(Surf::Tri),
					      "surf2grid:sbuf1");
	}
	proclist1[nsend] = plist[i];
	if (dim == 2) memcpy(&sbuf1[nsend*nbytes_surf],&lines[isurf],nbytes_surf);
	else memcpy(&sbuf1[nsend*nbytes_surf],&tris[isurf],nbytes_surf);
	nsend++;
      }
    }

    irregular = new Irregular(sparta);
    int nrecv1 = irregular->create_data_uniform(nsend,proclist1,1);
    char *rbuf1 = (char *) memory->smalloc(nrecv1*nbytes_surf,"surf2grid:rbuf");
    irregular->exchange_uniform(sbuf1,nbytes_surf,rbuf1);
    delete irregular;

    if (outflag) {
      MPI_Barrier(world);
      t2 = MPI_Wtime();
      tcomm1 += t2-t1;
    }

    // second irregular comm
    // identify which RCB proc owns each of my child cells at this level
    // send childID and my proc ID to RCB procs
    // nrecv2 = # of child cells I have copy of in RCB decomp

    int cx,cy,cz;
    double ctr[3];

    nsend = 0;

    for (icell = 0; icell < nlocal; icell++) {
      if (cells[icell].level != level) continue;
      if (cells[icell].nsplit <= 0) continue;

      ctr[0] = 0.5 * (cells[icell].lo[0] + cells[icell].hi[0]);
      ctr[1] = 0.5 * (cells[icell].lo[1] + cells[icell].hi[1]);
      ctr[2] = 0.5 * (cells[icell].lo[2] + cells[icell].hi[2]);
      id_find_child_uniform_level(level,0,boxlo,boxhi,ctr,cx,cy,cz);

      // glo/hi = single cell grid box

      glo[0] = cx; ghi[0] = cx;
      glo[1] = cy; ghi[1] = cy;
      glo[2] = cz; ghi[2] = cz;

      int np = 0;
      box_drop(glo,ghi,0,nprocs-1,gtree,np,plist);
      if (np != 1) error->one(FLERR,"Box drop of grid cell failed");

      if (nsend == maxsend2) {
	maxsend2 += DELTA_SEND;
	memory->grow(proclist2,maxsend2,"surf2grid:proclist2");
	sbuf2 = (Send2 *) memory->srealloc(sbuf2,maxsend2*sizeof(Send2),
					  "surf2grid:sbuf2");
      }

      proclist2[nsend] = plist[0];
      sbuf2[nsend].childID = cells[icell].id;
      sbuf2[nsend].proc = me;
      sbuf2[nsend].icell = icell;
      nsend++;
    }

    irregular = new Irregular(sparta);
    int nrecv2 = irregular->create_data_uniform(nsend,proclist2,1);
    Send2 *rbuf2 = (Send2 *) memory->smalloc(nrecv2*sizeof(Send2),"surf2grid:rbuf2");
    irregular->exchange_uniform((char *) sbuf2,sizeof(Send2),(char *) rbuf2);
    delete irregular;

    if (outflag) {
      MPI_Barrier(world);
      t3 = MPI_Wtime();
      tcomm2 += t3-t2;
    }

    // chash = hash with cell IDs I own in RCB decomp
    //   key = childID, value = index in my RCB list of child cells
    // phash = hash with all parent cell IDs of child cells
    //   key = parentID, value not used
    //   b/c RCB grid box is compact, size of phash should be small
    // rcblohi = lo/hi extents of each child cells

    MyHash *chash = new MyHash();
    MyHash *phash = new MyHash();
    RCBlohi *rcblohi =
      (RCBlohi *) memory->smalloc(nrecv2*sizeof(RCBlohi),"surf2grid:rcblohi");

    for (i = 0; i < nrecv2; i++) {
      childID = rbuf2[i].childID;
      (*chash)[childID] = i;
      id_lohi(childID,level,boxlo,boxhi,rcblohi[i].lo,rcblohi[i].hi);

      for (int ilevel = level; ilevel > 0; ilevel--) {
	parentID = id_parent_of_child(childID,ilevel);
	if (phash->find(parentID) != phash->end()) break;
	(*phash)[parentID] = 0;
	childID = parentID;
      }
    }

    // in RCB decomp, compute intersections between:
    //   my RCB child cells (only those that exist) and
    //   set of RCB surfs that overlap my RCB grid box
    // append results one by one to pairs = surf/grid intersections
    // loop over surfs:
    //   check if surf actually intersects with RCB box, else skip it
    //     could just be the bbox of surf overlaps with RCB box
    //   bblo/hi = intersection of surf bbox with RCB box
    //   recurse2d/3d starts with bblo/hi within parentID = 0 (sim box)
    //     will find every child cell in chash that surf intersects with
    //     checks for actual intersection are via cut2d/cut3d
    //   build list of pairs, one pair = surf index, cell index
    //     both are indices into received RCB surf/cells data

    Surf::Line *rcblines;
    Surf::Tri *rcbtris;
    if (dim == 2) rcblines = (Surf::Line *) rbuf1;
    else rcbtris = (Surf::Tri *) rbuf1;

    int npair = 0;
    int overlap;

    if (dim == 2) {
      for (i = 0; i < nrecv1; i++) {

	// skip surf if it does not intersect my RCB box

	overlap = cut2d->surf2grid_one(rcblines[i].p1,rcblines[i].p2,rcblo,rcbhi);
	if (!overlap) continue;

	// slo/hi = bbox around one surf
	
	surf->bbox_one(&rcblines[i],slo,shi);

	// bblo/hi = overlap of surf bbox with my RCB box
	
	bblo[0] = MAX(slo[0],rcblo[0]);
	bblo[1] = MAX(slo[1],rcblo[1]);
	bbhi[0] = MIN(shi[0],rcbhi[0]);
	bbhi[1] = MIN(shi[1],rcbhi[1]);
	bblo[2] = 0.0;
	bbhi[2] = 0.0;

	// find all my RCB child cells this surf intersects

	recurse2d(0,0,boxlo,boxhi,i,&rcblines[i],bblo,bbhi,
		  npair,maxpair,pairs,chash,phash);
      }

    } else {
      for (i = 0; i < nrecv1; i++) {

	// skip surf if it does not intersect my RCB box

	overlap = cut3d->surf2grid_one(rcbtris[i].p1,rcbtris[i].p2,rcbtris[i].p3,
				       rcblo,rcbhi);
	if (!overlap) continue;

	// slo/hi = bbox around one surf
	
	surf->bbox_one(&rcbtris[i],slo,shi);

	// bblo/hi = overlap of surf bbox with my RCB box
	
	bblo[0] = MAX(slo[0],rcblo[0]);
	bblo[1] = MAX(slo[1],rcblo[1]);
	bblo[2] = MAX(slo[2],rcblo[2]);
	bbhi[0] = MIN(shi[0],rcbhi[0]);
	bbhi[1] = MIN(shi[1],rcbhi[1]);
	bbhi[2] = MIN(shi[2],rcbhi[2]);

	// find all my RCB child cells this surf intersects

	recurse3d(0,0,boxlo,boxhi,i,&rcbtris[i],bblo,bbhi,
		  npair,maxpair,pairs,chash,phash);
      }
    }

    if (outflag) {
      MPI_Barrier(world);
      t4 = MPI_Wtime();
      tmap += t4-t3;
    }

    // third irregular comm
    // send each surf/grid intersection pair back to proc that owns grid cell

    int surfindex,cellindex;

    nsend = 0;

    for (i = 0; i < npair; i++) {
      if (nsend == maxsend3) {
	maxsend3 += DELTA_SEND;
	memory->grow(proclist3,maxsend3,"surf2grigd:proclist3");
	sbuf3 = (Send3 *) memory->srealloc(sbuf3,maxsend3*sizeof(Send3),
					  "surf2grid:sbuf3");
      }

      surfindex = pairs[i][0];
      cellindex = pairs[i][1];
      proclist3[i] = rbuf2[cellindex].proc;
      if (dim == 2) sbuf3[i].surfID = rcblines[surfindex].id;
      else sbuf3[i].surfID = rcbtris[surfindex].id;
      sbuf3[i].icell = rbuf2[cellindex].icell;
      nsend++;
    }

    irregular = new Irregular(sparta);
    int nrecv3 = irregular->create_data_uniform(nsend,proclist3,1);
    Send3 *rbuf3 = (Send3 *) memory->smalloc(nrecv3*sizeof(Send3),
					     "surf2grid:rbuf3");
    irregular->exchange_uniform((char *) sbuf3,sizeof(Send3),(char *) rbuf3);
    delete irregular;

    // process received cell/surf pairs back in simulation decomposition
    // set nsurf and csurfs for each cell (only child cells at this level)
    // 1st pass: count surfs in each cell, then allocate csurfs in each cell
    // 2nd pass: fill each cell's csurf list

    for (i = 0; i < nrecv3; i++) {
      icell = rbuf3[i].icell;
      cells[icell].nsurf++;
    }

    // skip sub cells since may exist in a restart

    for (icell = 0; icell < nlocal; icell++) {
      if (cells[icell].level != level) continue;
      if (cells[icell].nsplit <= 0) continue;
      nsurf = cells[icell].nsurf;
      if (nsurf) {
	if (nsurf > maxsurfpercell)
	  error->one(FLERR,"Too many surfs in one cell - set global surfmax");
	cells[icell].csurfs = csurfs->get(nsurf);
	cells[icell].nsurf = 0;
      }
    }

    for (i = 0; i < nrecv3; i++) {
      icell = rbuf3[i].icell;
      nsurf = cells[icell].nsurf;
      cells[icell].csurfs[nsurf] = rbuf3[i].surfID;
      cells[icell].nsurf++;
    }

    if (outflag) {
      MPI_Barrier(world);
      t5 = MPI_Wtime();
      tcomm3 += t5-t4;
    }

    // clean up for this level iteration

    memory->sfree(rbuf1);
    memory->sfree(rbuf2);
    memory->sfree(rbuf3);
    memory->sfree(rcblohi);
    delete chash;
    delete phash;
  }

  if (outflag) {
    MPI_Barrier(world);
    t1 = MPI_Wtime();
  }

  // clean up after all iterations

  memory->destroy(proclist1);
  memory->destroy(proclist2);
  memory->destroy(proclist3);
  memory->sfree(sbuf1);
  memory->sfree(sbuf2);
  memory->sfree(sbuf3);
  memory->destroy(pairs);
  memory->destroy(plist);
  memory->sfree(gtree);

  // non-distributed surfs:
  // each cell's csurf list currently stores surf IDs
  // convert them indices into global list stored by each proc
  // shash used to store IDs of entire global list

  if (!distributed) {
    lines = surf->lines;
    tris = surf->tris;
    int nslocal = surf->nlocal;

    MySurfHash shash;
    surfint *list;

    if (dim == 2) {
      for (i = 0; i < nslocal; i++)
	shash[lines[i].id] = i;
    } else {
      for (i = 0; i < nslocal; i++)
	shash[tris[i].id] = i;
    }

    for (icell = 0; icell < nlocal; icell++) {
      if (!cells[icell].nsurf) continue;
      if (cells[icell].nsplit <= 0) continue;

      list = cells[icell].csurfs;
      n = cells[icell].nsurf;

      for (i = 0; i < n; i++)
	list[i] = shash[list[i]];
    }
  }

  // distributed surfs:
  // rendezvous operation to obtain nlocal surfs for each proc
  //   these are the surfs that intersect child cells this proc owns
  // each grid cell requests a surf from proc that owns surf in mylines/mytris
  //   use shash to only do this once per surf
  // receive the surf and store in nlocal lines/tris

  if (distributed) {

    // ncount = # of unique surfs I need for my owned grid cells
    // store IDs of those surfs in shash

    MySurfHash shash;
    MyIterator it;
    surfint *list;
    int ncount = 0;

    for (icell = 0; icell < nlocal; icell++) {
      if (!cells[icell].nsurf) continue;
      if (cells[icell].nsplit <= 0) continue;

      list = cells[icell].csurfs;
      n = cells[icell].nsurf;

      for (i = 0; i < n; i++)
	if (shash.find(list[i]) == shash.end()) {
	  shash[list[i]] = 0;
	  ncount++;
	}
    }

    // allocate memory for rvous input

    int *proclist;
    memory->create(proclist,ncount,"surf2grid:proclist2");
    InRvous *inbuf =
      (InRvous *) memory->smalloc((bigint) ncount*sizeof(InRvous),
				  "surf2grid:inbuf");

    // create rvous inputs
    // proclist = owner of each surf

    surfint surfID;

    ncount = 0;
    for (it = shash.begin(); it != shash.end(); ++it) {
      surfID = it->first;
      proclist[ncount] = (surfID-1) % nprocs;
      inbuf[ncount].proc = me;
      inbuf[ncount].surfID = surfID;
      ncount++;
    }

    // perform rendezvous operation
    // each proc owns subset of surfs
    // receives all surf requests to return surf to each proc who needs it

    char *outbuf;
    int outbytes;
    if (dim == 2) outbytes = sizeof(OutRvousLine);
    else outbytes = sizeof(OutRvousTri);

    int nreturn = comm->rendezvous(1,ncount,(char *) inbuf,sizeof(InRvous),
				   0,proclist,rendezvous_surfrequest,
				   0,outbuf,outbytes,(void *) this);

    memory->destroy(proclist);
    memory->sfree(inbuf);

    // copy entire rendezvous output buf into realloced Surf lines/tris

    surf->nlocal = surf->nghost = 0;
    int nmax_old = surf->nmax;
    surf->nmax = surf->nlocal = nreturn;
    surf->grow(nmax_old);

    if (dim == 2) memcpy(surf->lines,outbuf,nreturn*sizeof(Surf::Line));
    else memcpy(surf->tris,outbuf,nreturn*sizeof(Surf::Tri));

    memory->sfree(outbuf);

    // reset Surf hash to point to surf list in lines/tris

    Surf::Line *lines = surf->lines;
    Surf::Tri *tris = surf->tris;

    if (dim == 2) {
      for (i = 0; i < nreturn; i++) {
        surfID = lines[i].id;
        shash[surfID] = i;
      }
    } else {
      for (i = 0; i < nreturn; i++) {
        surfID = tris[i].id;
        shash[surfID] = i;
      }
    }

    // reset csurfs list for each of my owned cells
    // from storing surfID to storing local index of that surfID

    for (icell = 0; icell < nlocal; icell++) {
      if (!cells[icell].nsurf) continue;
      if (cells[icell].nsplit <= 0) continue;

      list = cells[icell].csurfs;
      n = cells[icell].nsurf;

      for (i = 0; i < n; i++)
	list[i] = shash[list[i]];
    }
  }

  // for performance, sort each cell's csurfs list, same order as cell alg
  // mark cells with surfs as OVERLAP, only if has a non-transparent surf

  lines = surf->lines;
  tris = surf->tris;

  surfint *list;
  int nontrans;

  for (icell = 0; icell < nlocal; icell++) {
    if (!cells[icell].nsurf) continue;
    if (cells[icell].nsplit <= 0) continue;

    qsort(cells[icell].csurfs,cells[icell].nsurf,
	  sizeof(surfint),compare_surfIDs);

    list = cells[icell].csurfs;
    n = cells[icell].nsurf;
    nontrans = 0;

    if (dim == 2) {
      for (i = 0; i < n; i++) {
	if (!lines[list[i]].transparent) nontrans = 1;
	break;
      }
    } else {
      for (i = 0; i < n; i++) {
	if (!tris[list[i]].transparent) nontrans = 1;
	break;
      }
    }

    if (nontrans) cinfo[icell].type = OVERLAP;
  }

  if (outflag) {
    MPI_Barrier(world);
    t2 = MPI_Wtime();
    tcomm4 = t2-t1;
  }
}

/* ----------------------------------------------------------------------
   compute cut volume of each cell and any split cell info
   nsurf and csurfs list for each grid cell have already been computed
   if subflag = 1, create new owned split and sub cells as needed
     called from ReadSurf, RemoveSurf, MoveSurf, FixAblate
   if subflag = 0, split/sub cells already exist
     called from ReadRestart
   outflag = 1 for timing and statistics info
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

        // add nsplitone sub cells
        // collide and fixes also need to add cells

	for (i = 0; i < nsplitone; i++) {
	  isub = nlocal;
	  add_sub_cell(icell,1);
          if (collide) collide->add_grid_one();
          if (modify->n_pergrid) modify->add_grid_one();
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

  // stats on unmarked corner points in OVERLAP cells

  if (outflag) {
    int noverlap = 0;
    int ncorner = 0;
    for (int icell = 0; icell < nlocal; icell++) {
      if (cells[icell].nsplit <= 0) continue;
      if (cinfo[icell].type == OVERLAP) {
        noverlap++;
        if (cinfo[icell].corner[0] == UNKNOWN) ncorner++;
      }
    }

    bigint bncorner = ncorner;
    bigint bnoverlap = noverlap;
    bigint ncornerall,noverlapall;
    MPI_Allreduce(&bncorner,&ncornerall,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
    MPI_Allreduce(&bnoverlap,&noverlapall,1,MPI_SPARTA_BIGINT,MPI_SUM,world);

    if (comm->me == 0) {
      if (screen) fprintf(screen,"  " BIGINT_FORMAT " " BIGINT_FORMAT
                          " = cells overlapping surfs, "
                          "overlap cells with unmarked corner pts\n",
                          noverlapall,ncornerall);
      if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " " BIGINT_FORMAT
                           " = cells overlapping surfs, "
                           "overlap cells with unmarked corner pts\n",
                           noverlapall,ncornerall);
    }
  }

  // print info on unusual surf split cases

  if (dim == 3) {
    int ntiny = cut3d->ntiny;
    int alltiny;
    MPI_Allreduce(&ntiny,&alltiny,1,MPI_INT,MPI_SUM,world);
    if (alltiny && comm->me == 0) {
      if (screen) fprintf(screen,"  %d tiny edges removed\n",alltiny);
      if (logfile) fprintf(logfile,"  %d tiny edges removed\n",alltiny);
    }

    int nshrink = cut3d->nshrink;
    int allshrink;
    MPI_Allreduce(&nshrink,&allshrink,1,MPI_INT,MPI_SUM,world);
    if (allshrink && comm->me == 0) {
      if (screen) 
        fprintf(screen,"  %d cells shrunk to enable splitting\n",allshrink);
      if (logfile)
        fprintf(logfile,"  %d cells shrunk to enable splitting\n",allshrink);
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
   enumerate all child cells in chash which a single line intersects with
   done recursively, 1st call from surf2grid_surf_algorithm() uses parentID = root
   phash stores IDs of all parent cells for child cells in chash
   bblo/hi = portion of bounding box for surf that is wholly within parentID
   parentID = parent cell
   level = level of parent cell
   plo/phi = corner points of parent cell
   surfindex, line = surf element and its local index in caller
   npair, maxpair, pair = growing list of I,J gridcell/surf overlap pairs
   chash, phash = hashes of cell IDs for child and parent cells in this proc's RCB box
   cut2d->surf2grid_one() is used to determine actual overlap
------------------------------------------------------------------------- */

void Grid::recurse2d(cellint parentID, int level, double *plo, double *phi,
		     int surfindex, Surf::Line *line, double *bblo, double *bbhi,
                     int &npair, int &maxpair, int **&pairs,
		     MyHash *chash, MyHash *phash)
{
  int ix,iy,cflag,pflag,overlap;
  cellint ichild,childID;
  double celledge;
  double clo[3],chi[3];
  double newlo[3],newhi[3];

  double *p1 = line->p1;
  double *p2 = line->p2;

  int nx = plevels[level].nx;
  int ny = plevels[level].ny;
  int nbits = plevels[level].nbits;

  // ij lohi = indices for range of child cells overlapped by surf bbox
  // overlap = surf bbox include any interior of grid cell or touches its boundary
  // id_point_child() returns cell indices for >= lo and < hi
  // so check if lower range should be decremented
  
  int ilo,ihi,jlo,jhi,klo,khi;
  id_point_child(bblo,plo,phi,nx,ny,1,ilo,jlo,klo);
  id_point_child(bbhi,plo,phi,nx,ny,1,ihi,jhi,khi);

  celledge = plo[0] + ilo*(phi[0]-plo[0])/nx;
  if (bblo[0] <= celledge && ilo > 0) ilo--;
  celledge = plo[1] + jlo*(phi[1]-plo[1])/ny;
  if (bblo[1] <= celledge && jlo > 0) jlo--;
  
  // loop over range of grid cells between ij lohi inclusive
  // if cell is neither a child or parent cell in chash/phash, skip it
  // if line does not intersect cell, skip it
  // if cell is a child, add intersection to pairs list
  // if cell is a parent:
  //   recurse using new lohi for intersection of surf bbox with new parent cell

  for (iy = jlo; iy <= jhi; iy++) {
    for (ix = ilo; ix <= ihi; ix++) {
      ichild = (cellint) iy*nx + ix + 1;
      childID = parentID | (ichild << nbits);

      if (chash->find(childID) == chash->end()) cflag = 0;
      else cflag = 1;
      if (phash->find(childID) == phash->end()) pflag = 0;
      else pflag = 1;
      if (!cflag && !pflag) continue;

      grid->id_child_lohi(level,plo,phi,ichild,clo,chi);
      overlap = cut2d->surf2grid_one(p1,p2,clo,chi);
      if (!overlap) continue;

      if (cflag) {
	if (npair == maxpair) {
	  maxpair += DELTA_SEND;
	  memory->grow(pairs,maxpair,2,"surf2grid:pairs");
	}
	pairs[npair][0] = surfindex;
	pairs[npair][1] = (*chash)[childID];
	npair++;
	continue;
      }

      if (pflag) {
	newlo[0] = MAX(bblo[0],clo[0]);
	newlo[1] = MAX(bblo[1],clo[1]);
	newhi[0] = MIN(bbhi[0],chi[0]);
	newhi[1] = MIN(bbhi[1],chi[1]);
	recurse2d(childID,level+1,clo,chi,surfindex,line,newlo,newhi,
		  npair,maxpair,pairs,chash,phash);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   enumerate all child cells in chash which a single tri intersects with
   done recursively, 1st call from surf2grid_surf_algorithm() uses parentID = root
     phash stores IDs of all parent cells for child cells in chash
   bblo/hi = portion of bounding box for surf that is wholly within parentID
   parentID = parent cell
   level = level of parent cell
   plo/phi = corner points of parent cell
   surfindex, tri = surf element and its local index in caller
   npair, maxpair, pair = growing list of I,J gridcell/surf overlap pairs
   chash, phash = hashes of cell IDs for child and parent cells in this proc's RCB box
   cut3d->surf2grid_one() is used to determine actual overlap
------------------------------------------------------------------------- */

void Grid::recurse3d(cellint parentID, int level, double *plo, double *phi,
		     int surfindex, Surf::Tri *tri, double *bblo, double *bbhi,
                     int &npair, int &maxpair, int **&pairs,
		     MyHash *chash, MyHash *phash)
{
  int ix,iy,iz,cflag,pflag,overlap;
  cellint ichild,childID;
  double celledge;
  double clo[3],chi[3];
  double newlo[3],newhi[3];

  double *p1 = tri->p1;
  double *p2 = tri->p2;
  double *p3 = tri->p3;

  int nx = plevels[level].nx;
  int ny = plevels[level].ny;
  int nz = plevels[level].nz;
  int nbits = plevels[level].nbits;

  // ijk lohi = indices for range of child cells overlapped by surf bbox
  // overlap = surf bbox include any interior of grid cell or touches its boundary
  // id_point_child() returns cell indices for >= lo and < hi
  // so check if lower range should be decremented
  
  int ilo,ihi,jlo,jhi,klo,khi;
  id_point_child(bblo,plo,phi,nx,ny,nz,ilo,jlo,klo);
  id_point_child(bbhi,plo,phi,nx,ny,nz,ihi,jhi,khi);

  celledge = plo[0] + ilo*(phi[0]-plo[0])/nx;
  if (bblo[0] <= celledge && ilo > 0) ilo--;
  celledge = plo[1] + jlo*(phi[1]-plo[1])/ny;
  if (bblo[1] <= celledge && jlo > 0) jlo--;
  celledge = plo[2] + klo*(phi[2]-plo[2])/nz;
  if (bblo[2] <= celledge && klo > 0) klo--;

  // loop over range of grid cells between ij lohi inclusive
  // if cell is neither a child or parent cell in chash/phash, skip it
  // if tri does not intersect cell, skip it
  // if cell is a child, add intersectino to pairs list
  // if pairs is a parent:
  //   recurse using new lohi for intersection of surf bbox with new parent cell

  for (iz = klo; iz <= khi; iz++) {
    for (iy = jlo; iy <= jhi; iy++) {
      for (ix = ilo; ix <= ihi; ix++) {
        ichild = (cellint) iz*nx*ny + (cellint) iy*nx + ix + 1;
	childID = parentID | (ichild << nbits);

	if (chash->find(childID) == chash->end()) cflag = 0;
	else cflag = 1;
	if (phash->find(childID) == phash->end()) pflag = 0;
	else pflag = 1;
	if (!cflag && !pflag) continue;

	grid->id_child_lohi(level,plo,phi,ichild,clo,chi);
	overlap = cut3d->surf2grid_one(p1,p2,p3,clo,chi);
	if (!overlap) continue;

	if (cflag) {
	  if (npair == maxpair) {
	    maxpair += DELTA_SEND;
	    memory->grow(pairs,maxpair,2,"surf2grid:pairs");
	  }
	  pairs[npair][0] = surfindex;
	  pairs[npair][1] = (*chash)[childID];
	  npair++;
	  continue;
	}

	if (pflag) {
	  newlo[0] = MAX(bblo[0],clo[0]);
	  newlo[1] = MAX(bblo[1],clo[1]);
	  newlo[2] = MAX(bblo[2],clo[2]);
	  newhi[0] = MIN(bbhi[0],chi[0]);
	  newhi[1] = MIN(bbhi[1],chi[1]);
	  newhi[2] = MIN(bbhi[2],chi[2]);
	  recurse3d(childID,level+1,clo,chi,surfindex,tri,newlo,newhi,
		    npair,maxpair,pairs,chash,phash);
	}
      }
    }
  }
}

/* ----------------------------------------------------------------------
   recursive method to partition a block of uniform grid cells
   uses RCB to create one sub-block per processor
   xyz lo/hi = extent of initial grid block
   proc lower/upper = 0 and Nprocs-1 initially
   output: gtree = tree of RCB cuts
------------------------------------------------------------------------- */

void Grid::partition_grid(int proclower, int procupper,
			  int xlo, int xhi, int ylo, int yhi, int zlo, int zhi,
			  GridTree *gtree)
{
  // end recursion when partition is a single proc

  if (proclower == procupper) return;

  int procmid = proclower + (procupper-proclower) / 2 + 1;
  int nplower = procmid-proclower;
  int npupper = procupper-procmid + 1;

  int xrange = xhi-xlo + 1;
  int yrange = yhi-ylo + 1;
  int zrange = zhi-zlo + 1;

  if (xrange >= yrange && xrange >= zrange) {
    int mid = xlo + static_cast<int> ((0.5*nplower/npupper) * xrange);
    gtree[procmid].dim = 0;
    gtree[procmid].cut = mid;
    partition_grid(proclower,procmid-1,xlo,mid-1,ylo,yhi,zlo,zhi,gtree);
    partition_grid(procmid,procupper,mid,xhi,ylo,yhi,zlo,zhi,gtree);
  } else if (yrange >= zrange) {
    int mid = ylo + static_cast<int> ((0.5*nplower/npupper) * yrange);
    gtree[procmid].dim = 1;
    gtree[procmid].cut = mid;
    partition_grid(proclower,procmid-1,xlo,xhi,ylo,mid-1,zlo,zhi,gtree);
    partition_grid(procmid,procupper,xlo,xhi,mid,yhi,zlo,zhi,gtree);
  } else {
    int mid = zlo + static_cast<int> ((0.5*nplower/npupper) * zrange);
    gtree[procmid].dim = 2;
    gtree[procmid].cut = mid;
    partition_grid(proclower,procmid-1,xlo,xhi,ylo,yhi,zlo,mid-1,gtree);
    partition_grid(procmid,procupper,xlo,xhi,ylo,yhi,mid,zhi,gtree);
  }
}

/* ----------------------------------------------------------------------
   recursive method to identify my sub-block uniform grid cell block
   xyz lohi = extent of initial grid block
   proc lower/upper = 0 and Nprocs-1 initially
   traverses RCB gtree of cuts to zoom in on this processor
   output: xyz lohi will be overwritten with this proc's sub-block extent
------------------------------------------------------------------------- */

void Grid::mybox(int me, int proclower, int procupper,
		 int &xlo, int &xhi, int &ylo, int &yhi, int &zlo, int &zhi,
		 GridTree *gtree)
{
  // end recursion when partition is a single proc

  if (proclower == procupper) return;

  int procmid = proclower + (procupper-proclower) / 2 + 1;

  if (me < procmid) {
    if (gtree[procmid].dim == 0) xhi = gtree[procmid].cut-1;
    else if (gtree[procmid].dim == 1) yhi = gtree[procmid].cut-1;
    else zhi = gtree[procmid].cut-1;
    mybox(me,proclower,procmid-1,xlo,xhi,ylo,yhi,zlo,zhi,gtree);
  } else {
    if (gtree[procmid].dim == 0) xlo = gtree[procmid].cut;
    else if (gtree[procmid].dim == 1) ylo = gtree[procmid].cut;
    else zlo = gtree[procmid].cut;
    mybox(me,procmid,procupper,xlo,xhi,ylo,yhi,zlo,zhi,gtree);
  }
}

/* ----------------------------------------------------------------------
   recursive method to drop a grid box down RCB tree to identify procs it overlaps
   lo/hi = indices (0 to N-1) of uniform grid block for range of box
   proc lower/upper = 0 and Nprocs-1 initially
   traverses RCB gtree to drop box on one or both sides of each cut
   output: noverlap = # of proc sub-boxes the box overlaps
   output: overlap = list of proc IDs the box overlaps
   overlap vector is allocated by caller
------------------------------------------------------------------------- */

void Grid::box_drop(int *lo, int *hi, int proclower, int procupper,
		    GridTree *gtree, int &noverlap, int *overlap)
{
  // end recursion when partition is a single proc
  // add proc to overlap list

  if (proclower == procupper) {
    overlap[noverlap++] = proclower;
    return;
  }

  // drop box on each side of cut it extends beyond
  // use of < and >= criteria are important
  // procmid = 1st processor in upper half of partition
  //         = location in tree that stores this cut
  // dim = 0,1,2 dimension of cut
  // cut = position of cut

  int procmid = proclower + (procupper - proclower) / 2 + 1;
  int dim = gtree[procmid].dim;
  int cut = gtree[procmid].cut;

  if (lo[dim] < cut)
    box_drop(lo,hi,proclower,procmid-1,gtree,noverlap,overlap);
  if (hi[dim] >= cut)
    box_drop(lo,hi,procmid,procupper,gtree,noverlap,overlap);
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

  InRvous *in = (InRvous *) inbuf;

  // proclist = list of procs who made requests
  // outbuf = list of lines or tris to pass back to comm->rendezvous

  memory->create(proclist,n,"gridsurf:proclist");

  int surfbytes;
  if (dim == 2) surfbytes = sizeof(OutRvousLine);
  else surfbytes = sizeof(OutRvousTri);
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

  hashfilled = 0;

  // if surfs no longer exist, set cell type to OUTSIDE, else UNKNOWN
  // set corner points of every cell to UNKNOWN

  int celltype = UNKNOWN;
  if (!surf->exist) celltype = OUTSIDE;

  // compress cell list
  // collide and fixes also need to do the same

  int nlocal_prev = nlocal;

  int icell = 0;
  while (icell < nlocal) {
    if (cells[icell].nsplit <= 0) {
      if (icell != nlocal-1) {
        memcpy(&cells[icell],&cells[nlocal-1],sizeof(ChildCell));
        memcpy(&cinfo[icell],&cinfo[nlocal-1],sizeof(ChildInfo));
        if (collide && collide->ngroups) collide->copy_grid_one(nlocal-1,icell);
        if (modify->n_pergrid) modify->copy_grid_one(nlocal-1,icell);
      }
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

  // reset final grid cell count in collide and fixes

  if (collide) collide->reset_grid_count(nlocal);
  if (modify->n_pergrid) modify->reset_grid_count(nlocal);

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
   assign all particles in a split icell to appropriate sub cells
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
   check if particle at x is outside any surfs in icell
   icell can be split or unsplit cell, not a sub cell
   if outside, return 1, else return 0
   // NOTE: make this only for implicit surfs
------------------------------------------------------------------------- */

int Grid::outside_surfs(int icell, double *x,
                        Cut3d *cut3d_caller, Cut2d *cut2d_caller)
{
  if (cells[icell].nsurf == 0) {
    if (cinfo[icell].type == INSIDE) return 0;
    else return 1;
  }

  // set xnew to midpt of first line or center pt of first triangle
  // for implicit surfs this is guaranteed to be a pt in or on icell
  // then displace it by EPSSURF in the line/tri norm direction
  // reason for this:
  //   want to insure an inside particle is flagged
  //   requires a ray from inside particle x to xnew intersects a surf
  //   if no intersection, logic below assumes particle is outside
  //   if xnew is midpt of tri, then an inside particle may have no intersection
  //     due to round-off

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  surfint *csurfs = cells[icell].csurfs;

  int dim = domain->dimension;
  double xnew[3],edge[3];
  double edgelen,minedge,displace;

  int isurf = csurfs[0];

  if (dim == 2) {
    xnew[0] = 0.5 * (lines[isurf].p1[0] + lines[isurf].p2[0]);
    xnew[1] = 0.5 * (lines[isurf].p1[1] + lines[isurf].p2[1]);
    xnew[2] = 0.0;

    MathExtra::sub3(lines[isurf].p1,lines[isurf].p2,edge);
    minedge = MathExtra::len3(edge);

  } else {
    double onethird = 1.0/3.0;
    xnew[0] = onethird *
      (tris[isurf].p1[0] + tris[isurf].p2[0] + tris[isurf].p3[0]);
    xnew[1] = onethird *
      (tris[isurf].p1[1] + tris[isurf].p2[1] + tris[isurf].p3[1]);
    xnew[2] = onethird *
      (tris[isurf].p1[2] + tris[isurf].p2[2] + tris[isurf].p3[2]);

    MathExtra::sub3(tris[isurf].p1,tris[isurf].p2,edge);
    edgelen = MathExtra::len3(edge);
    minedge = edgelen;
    MathExtra::sub3(tris[isurf].p2,tris[isurf].p3,edge);
    edgelen = MathExtra::len3(edge);
    minedge = MIN(minedge,edgelen);
    MathExtra::sub3(tris[isurf].p3,tris[isurf].p1,edge);
    minedge = MIN(minedge,edgelen);
  }

  displace = EPSSURF * minedge;
  xnew[0] += displace*tris[isurf].norm[0];
  xnew[1] += displace*tris[isurf].norm[1];
  xnew[2] += displace*tris[isurf].norm[2];

  // loop over surfs, ray-trace from x to xnew, see which surf is hit first
  // if no surf is hit (roundoff), assume particle is outside

  int m,cflag,hitflag,side,minside;
  double param,minparam;
  double xc[3];
  Surf::Line *line;
  Surf::Tri *tri;

  int nsurf = cells[icell].nsurf;

  cflag = 0;
  minparam = 2.0;
  for (m = 0; m < nsurf; m++) {
    isurf = csurfs[m];
    if (dim == 3) {
      tri = &tris[isurf];
      hitflag = Geometry::
        line_tri_intersect(x,xnew,tri->p1,tri->p2,tri->p3,
                           tri->norm,xc,param,side);
    } else {
      line = &lines[isurf];
      hitflag = Geometry::
        line_line_intersect(x,xnew,line->p1,line->p2,line->norm,xc,param,side);
    }
    if (hitflag && param < minparam) {
      cflag = 1;
      minparam = param;
      minside = side;
    }
  }

  if (!cflag) return 1;
  if (minside == SOUTSIDE || minside == ONSURF2OUT) return 1;
  return 0;
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

/* ----------------------------------------------------------------------
   request N-length vector from csubs
   called by ReadRestart
------------------------------------------------------------------------- */

int *Grid::csubs_request(int n)
{
  int *ptr = csubs->vget();
  csubs->vgot(n);
  return ptr;
}

// ----------------------------------------------------------------------
// private class methods
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   comparison function invoked by qsort() called by surf2grid algorithm
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

  bigint bscount = scount;
  bigint bstotal = stotal;
  bigint scountall,stotalall;
  int smaxall;
  double sratioall;
  MPI_Allreduce(&bscount,&scountall,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  MPI_Allreduce(&bstotal,&stotalall,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  MPI_Allreduce(&smax,&smaxall,1,MPI_INT,MPI_MAX,world);
  MPI_Allreduce(&sratio,&sratioall,1,MPI_DOUBLE,MPI_MIN,world);

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"  " BIGINT_FORMAT " = cells with surfs\n",scountall);
      fprintf(screen,"  " BIGINT_FORMAT
              " = total surfs in all grid cells\n",stotalall);
      fprintf(screen,"  %d = max surfs in one grid cell\n",smaxall);
      fprintf(screen,"  %g = min surf-size/cell-size ratio\n",sratioall);
    }
    if (logfile) {
      fprintf(logfile,"  " BIGINT_FORMAT " = cells with surfs\n",scountall);
      fprintf(logfile,"  " BIGINT_FORMAT
              " = total surfs in all grid cells\n",stotalall);
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
   volume for one surf is projection to lower z face (3d) or y face (2d)
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

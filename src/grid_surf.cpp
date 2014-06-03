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
#include "comm.h"
#include "surf.h"
#include "cut2d.h"
#include "cut3d.h"
#include "math_const.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

#define BIG 1.0e20
#define MAXSPLITPERCELL 10

enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};   // several files

// operations for surfaces in grid cells

/* ----------------------------------------------------------------------
   allocate page data structs to hold variable-length surf info
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

/* ----------------------------------------------------------------------
   map surf elements into owned grid cells
   if subflag = 1, create new owned split and sub cells as needed
     called from ReadSurf
   if subflag = 0, split/sub cells already exist
     called from ReadRestart
   in cells: set nsurf, csurfs, nsplit, isplit
   in cinfo: set type, corner, volume
   initialize sinfo as needed
------------------------------------------------------------------------- */

void Grid::surf2grid(int subflag)
{
  int i,isub,nsurf,nsplit,xsub;
  int *surfmap,*ptr;
  int corner[8];
  double *lo,*hi,*vols;
  double xsplit[3];
  ChildCell *c;
  SplitInfo *s;
  Cut2d *cut2d;
  Cut3d *cut3d;

  int dim = domain->dimension;

  double *slo = surf->bblo;
  double *shi = surf->bbhi;

  if (dim == 3) cut3d = new Cut3d(sparta);
  else cut2d = new Cut2d(sparta,domain->axisymmetric);

  // compute overlap of surfs with each cell I own
  // info stored in nsurf,csurfs

  //double t1 = MPI_Wtime();

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

  //double t2 = MPI_Wtime();
  //printf("TIME %g\n",t2-t1);
  
  surf2grid_stats();

  // compute cut volume and possible split of each grid cell by surfs
  // decrement nunsplitlocal if convert an unsplit cell to split cell
  // if nsplit > 1, create new split cell sinfo and sub-cells

  int ncurrent = nlocal;
  for (int icell = 0; icell < ncurrent; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    if (cells[icell].nsurf == 0) continue;

    surfmap = csplits->vget();
    c = &cells[icell];

    if (dim == 3)
      nsplit = cut3d->split(c->id,c->lo,c->hi,c->nsurf,c->csurfs,
                            vols,surfmap,cinfo[icell].corner,xsub,xsplit);
    else
      nsplit = cut2d->split(c->id,c->lo,c->hi,c->nsurf,c->csurfs,
                            vols,surfmap,cinfo[icell].corner,xsub,xsplit);

    if (nsplit == 1) cinfo[icell].volume = vols[0];

    else if (subflag) {
      cells[icell].nsplit = nsplit;
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

      for (i = 0; i < nsplit; i++) {
        isub = nlocal;
        add_sub_cell(icell,1);
        cells[isub].nsplit = -i;
        cinfo[isub].volume = vols[i];
        ptr[i] = isub;
      }
      
      csplits->vgot(cells[icell].nsurf);
      csubs->vgot(nsplit);

    } else {
      if (cells[icell].nsplit != nsplit) error->one(FLERR,"");
      
      s = &sinfo[cells[icell].isplit];
      s->csplits = surfmap;
      s->xsub = xsub;
      s->xsplit[0] = xsplit[0];
      s->xsplit[1] = xsplit[1];
      if (dim == 3) s->xsplit[2] = xsplit[2];
      else s->xsplit[2] = 0.0;

      ptr = s->csubs;
      for (i = 0; i < nsplit; i++) {
        isub = ptr[i];
        cinfo[isub].volume = vols[i];
      }
      
      csplits->vgot(cells[icell].nsurf);
    }
  }

  //double t3 = MPI_Wtime();
  //printf("TIME %g\n",t3-t2);

  if (dim == 3) delete cut3d;
  else delete cut2d;
}

/* ----------------------------------------------------------------------
   remove all surf info from owned grid cells
   called before reassigning surfs to grid cells
   changes cells data structure since sub cells are removed
------------------------------------------------------------------------- */

void Grid::clear_surf()
{
  // reset current grid cells as if no surfs existed
  // discard sub cells
  // compact cells and cinfo arrays
  // set values in cells/cinfo as if no surfaces

  int dimension = domain->dimension;
  int ncorner = 8;
  if (dimension == 2) ncorner = 4;
  double *lo,*hi;

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
      icell++;
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
   remove all surf info from owned grid cells
   called before reassigning surfs to grid cells
   changes cells data structure since sub cells are removed
------------------------------------------------------------------------- */

void Grid::clear_surf_restart()
{
  // reset current grid cells as if no surfs existed
  // discard sub cells
  // compact cells and cinfo arrays
  // set values in cells/cinfo as if no surfaces

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
   new surf-based methods
------------------------------------------------------------------------- */

void Grid::surf2grid_stats()
{
  int i,m;
  double cmax,len,area;
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
	area = surf->tri_size(cells[icell].csurfs[i],len);
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
  int i,icell;

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

/* ---------------------------------------------------------------------- */

double Grid::flow_volume()
{
  double zarea;
  double *p1,*p2,*p3;

  Surf::Point *pts = surf->pts;
  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  double volume = 0.0;

  if (domain->dimension == 3) {
    for (int i = 0; i < surf->ntri; i++) {
      p1 = pts[tris[i].p1].x;
      p2 = pts[tris[i].p2].x;
      p3 = pts[tris[i].p3].x;
      zarea = 0.5 * ((p2[0]-p1[0])*(p3[1]-p1[1]) - (p2[1]-p1[1])*(p3[0]-p1[0]));
      volume -= zarea * ((p1[2]+p2[2]+p3[2])/3.0 - boxlo[2]);
    }
    if (volume <= 0.0) 
      volume += (boxhi[0]-boxlo[0]) * (boxhi[1]-boxlo[1]) * 
        (boxhi[2]-boxlo[2]); 
 
  } else if (domain->axisymmetric) {
    for (int i = 0; i < surf->nline; i++) {
      p1 = pts[lines[i].p1].x;
      p2 = pts[lines[i].p2].x;
      volume -= 
        MY_PI3 * (p1[1]*p1[1] + p1[1]*p2[1] + p2[1]*p2[1]) * (p2[0]-p1[0]);
    }
    if (volume <= 0.0) 
      volume += MY_PI * boxhi[1]*boxhi[1] * (boxhi[0]-boxlo[0]);

  } else {
    for (int i = 0; i < surf->nline; i++) {
      p1 = pts[lines[i].p1].x;
      p2 = pts[lines[i].p2].x;
      volume -= (0.5*(p1[1]+p2[1]) - boxlo[1]) * (p2[0]-p1[0]);
    }
    if (volume <= 0.0) volume += (boxhi[0]-boxlo[0]) * (boxhi[1]-boxlo[1]); 
  }
  
  return volume;
}

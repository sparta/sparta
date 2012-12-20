/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2012) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "cut2d.h"
#include "domain.h"
#include "grid.h"
#include "surf.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{CELLUNKNOWN,CELLOUTSIDE,CELLINSIDE,CELLOVERLAP};   // same as Grid
enum{CORNERUNKNOWN,CORNEROUTSIDE,CORNERINSIDE,CORNEROVERLAP};  // samde as Grid

enum{POINTINSIDE,POINTBORDER,POINTOUTSIDE};
enum{LOOPINSIDE,LOOPBORDER};

enum{INSIDE,OUTSIDE};
enum{ASCEND,DESCEND};

#define ENTRY -1
#define EXIT 1
#define EPSILON 1.0e-6
#define TOLERANCE 1.0e-10

#define VERBOSE 0

/* ---------------------------------------------------------------------- */

Cut2d::Cut2d(SPARTA *sparta) : Pointers(sparta) {}

/* ---------------------------------------------------------------------- */

Cut2d::~Cut2d() {}

/* ----------------------------------------------------------------------
   compute intersections of surfs with grid cells
   for now, each proc computes for all cells
   done via 2 loops, one to count intersections, one to populate csurfs
   sets nsurf,csurfs for every grid cell with indices into global surf list
   also allocates csplits for every grid cell
------------------------------------------------------------------------- */

void Cut2d::surf2grid()
{
  int i,j,k,m,icell;
  int ilo,ihi,jlo,jhi,klo,khi;
  double lo[3],hi[3];
  double *x1,*x2,*x3;

  Grid::OneCell *cells = grid->cells;
  int ncell = grid->ncell;
  int nx = grid->nx;
  int ny = grid->ny;
  double xdeltainv = grid->xdeltainv;
  double ydeltainv = grid->ydeltainv;

  // epsilon = EPSILON fraction of largest box length

  double bmax = MAX(domain->xprd,domain->yprd);
  double epsilon = EPSILON * bmax;

  // count[M] = # of surfs overlapping grid cell M

  int *count;
  memory->create(count,ncell,"grid:count");
  for (m = 0; m < ncell; m++) count[m] = 0;

  // tally count by double loop over surfs and grid cells within surf bbox
  // lo/hi = bounding box around surf
  // ijk lo/hi = grid index bounding box around surf
  // add epsilon to insure surf is counted in any cell it touches
  // icell = index of a grid cell within bounding box
  // NOTE: this logic is specific to regular Nx by Ny by Nz grid

  Surf::Point *pts = surf->pts;
  Surf::Line *lines = surf->lines;
  double xlo = domain->boxlo[0];
  double ylo = domain->boxlo[0];

  int nsurf = surf->nline;

  for (m = 0; m < nsurf; m++) {
    x1 = pts[lines[m].p1].x;
    x2 = pts[lines[m].p2].x;

    lo[0] = MIN(x1[0],x2[0]);
    hi[0] = MAX(x1[0],x2[0]);

    lo[1] = MIN(x1[1],x2[1]);
    hi[1] = MAX(x1[1],x2[1]);

    ilo = MAX(0,static_cast<int> ((lo[0]-xlo-epsilon)*xdeltainv));
    ihi = MIN(nx-1,static_cast<int> ((hi[0]-xlo+epsilon)*xdeltainv));
    jlo = MAX(0,static_cast<int> ((lo[1]-ylo-epsilon)*ydeltainv));
    jhi = MIN(ny-1,static_cast<int> ((hi[1]-ylo+epsilon)*ydeltainv));

    for (i = ilo; i <= ihi; i++)
      for (j = jlo; j <= jhi; j++) {
        icell = j*nx + i;
        if (cliptest(x1,x2,cells[icell].lo,cells[icell].hi))
          count[icell]++;
      }
  }

  // allocate ragged csurfs and csplits array
  // csurfs[I][J] = index of Jth surf in global cell I
  // csplits will be setup in split()

  memory->destroy(grid->csurfs);
  memory->destroy(grid->csplits);
  memory->create_ragged(grid->csurfs,ncell,count,"grid:csurfs");
  memory->create_ragged(grid->csplits,ncell,count,"grid:csplits");

  // NOTE: this logic is specific to regular Nx by Ny by Nz grid
  // populate csurfs with same double loop

  int **csurfs = grid->csurfs;
  for (m = 0; m < ncell; m++) count[m] = 0;
  
  for (m = 0; m < nsurf; m++) {
    x1 = pts[lines[m].p1].x;
    x2 = pts[lines[m].p2].x;

    lo[0] = MIN(x1[0],x2[0]);
    hi[0] = MAX(x1[0],x2[0]);

    lo[1] = MIN(x1[1],x2[1]);
    hi[1] = MAX(x1[1],x2[1]);

    ilo = MAX(0,static_cast<int> ((lo[0]-xlo-epsilon)*xdeltainv));
    ihi = MIN(nx-1,static_cast<int> ((hi[0]-xlo+epsilon)*xdeltainv));
    jlo = MAX(0,static_cast<int> ((lo[1]-ylo-epsilon)*ydeltainv));
    jhi = MIN(ny-1,static_cast<int> ((hi[1]-ylo+epsilon)*ydeltainv));

    for (i = ilo; i <= ihi; i++)
      for (j = jlo; j <= jhi; j++) {
        icell = j*nx + i;
        if (cliptest(x1,x2,cells[icell].lo,cells[icell].hi))
          csurfs[icell][count[icell]++] = m;
      }
  }

  // set per-cell surf count

  for (icell = 0; icell < ncell; icell++)
    cells[icell].nsurf = count[icell];

  // clean up
    
  memory->destroy(count);
}

/* ----------------------------------------------------------------------
   calculate split area(s) of a single cell containing lines
   for each grid cell: set type, cflags, nsplit, volume
   for now, each proc works on all cells
------------------------------------------------------------------------- */

void Cut2d::split()
{
  int i,j,m;
  int cornerflag[4];

  Grid::OneCell *cells = grid->cells;
  int **csurfs = grid->csurfs;
  int **csplits = grid->csplits;
  int **cflags = grid->cflags;
  int ncell = grid->ncell;
  int me = comm->me;

  for (i = 0; i < ncell; i++) {
    nsurf = cells[i].nsurf;
    if (nsurf) {
      if (VERBOSE) {
        printf("\nCELL %d %g %g %g %g\n",cells[i].id,
               cells[i].lo[0],cells[i].hi[0],cells[i].lo[1],cells[i].hi[1]);
      }
      // TEMP for VERBOSE output
      icell = i;
      line2pl(nsurf,csurfs[i]);
      weiler_intersect(cells[i].lo,cells[i].hi,cornerflag);
      weiler_walk(cells[i].lo,cells[i].hi,cornerflag);
      loop2pg(cells[i].lo,cells[i].hi,cornerflag);
      if (pg.n > 1) surf2pg(nsurf,csurfs[i],csplits[i]);

      if (cells[i].proc == me) {
        cells[i].type = CELLOVERLAP;
        m = cells[i].local;
        cflags[m][0] = cornerflag[0];
        cflags[m][1] = cornerflag[1];
        cflags[m][2] = cornerflag[2];
        cflags[m][3] = cornerflag[3];
      }

      if (VERBOSE) {
        printf("SPLIT %d %d: %d %d %d %d\n",cells[i].id,cells[i].type,
               cornerflag[0],cornerflag[1],cornerflag[2],cornerflag[3]);
        printf("  AREAS:");
        for (j = 0; j < pg.n; j++)
          printf(" %g",pg[j].area);
        printf("\n");
        if (pg.n > 1) {
          printf("  SURFMAP:");
          for (j = 0; j < nsurf; j++)
            printf(" %d",csplits[i][j]);
          printf("\n");
        }
      }

      // here is where need to create split cells with volumes
      // needs to be done for all cells, not just owned cells

      cells[i].nsplit = pg.n;
      if (pg.n == 1) cells[i].volume = pg[0].area;
      else {
        cells[i].volume = 0.0;
        for (j = 0; j < pg.n; j++)
          cells[i].volume += pg[j].area;
      }

    } else {
      if (cells[i].proc == me) {
        cells[i].type = CELLUNKNOWN;
        m = cells[i].local;
        cflags[m][0] = cflags[m][1] = cflags[m][2] = cflags[m][3] = 
          CORNERUNKNOWN;
      }
      cells[i].nsplit = 1;
      cells[i].volume = (cells[i].hi[0]-cells[i].lo[0]) * 
        (cells[i].hi[1]-cells[i].lo[1]);
    }
  }

  if (VERBOSE) printf("\n");
}

/* ----------------------------------------------------------------------
   build PL = list of polylines
   from list of N lines that intersect one cell
------------------------------------------------------------------------- */

void Cut2d::line2pl(int n, int *list)
{
  int i,m,pt;
  PLone *datum;
  MyList<PLone*> *one;

  Surf::Line *lines = surf->lines;

  ppool.reset();

  used.grow(n);
  startpts.grow(n);
  endpts.grow(n);
  for (i = 0; i < n; i++) {
    used[i] = 0;
    startpts[i] = lines[list[i]].p1;
    endpts[i] = lines[list[i]].p2;
  }

  int unused = 1;
  int usedzero = 0;

  int npl = 0;
  while (unused) {
    i = usedzero;
    used[i] = 1;
    pl.grow(npl+1);
    one = &pl[npl];
    one->reset();
    datum = ppool.get();
    datum->iline = list[i];
    one->append(datum);

    while (1) {
      pt = lines[one->last->iline].p2;
      for (i = 0; i < n; i++)
        if (pt == startpts[i]) break;
      if (i == n) break;
      if (used[i]) break;
      used[i] = 1;
      datum = ppool.get();
      datum->iline = list[i];
      one->append(datum);
    }
    while (1) {
      pt = lines[one->first->iline].p1;
      for (i = 0; i < n; i++)
        if (pt == endpts[i]) break;
      if (i == n) break;
      if (used[i]) break;
      used[i] = 1;
      datum = ppool.get();
      datum->iline = list[i];
      one->prepend(datum);
    }
    npl++;


    for (usedzero = 0; usedzero < n; usedzero++)
      if (used[usedzero] == 0) break;
    if (usedzero == n) unused = 0;
  }

  pl.n = npl;

  if (VERBOSE) {
    printf("LIST OF LINES %d:",n);
    for (i = 0; i < n; i++) printf(" %d",list[i]);
    printf("\n");
    printf("PL %d\n",npl);
    for (i = 0; i < npl; i++) {
      one = &pl[i];
      printf("  ONE %d:",one->n);
      for (PLone *p = one->first; p; p = p->next)
        printf(" %d",p->iline);
      printf("\n");
    }
  }
}

/* ----------------------------------------------------------------------
   Weiler-Atherton algorithm for intersecting polylines with grid cell
   build and return CPTS and OPTS
   CPTS = 4 cell corner pts, interleaved with intersection pts
   OPTS = one set of object pts for each polyline
          pts in polyline, interleaved with intersection pts
   also return cornerflags
------------------------------------------------------------------------- */

void Cut2d::weiler_intersect(double *lo, double *hi, int *cornerflag)
{
  int i,ipt,iline,ic,flag,bflag;
  int nopt;
  int nentry,nexit,lastflag;
  double a[2],b[2];
  double *p,*q,*norm;
  Opt *opt;
  Cpt *cpt;
  MyList<PLone*> *polyline;
  PLone *poly;

  Surf::Line *lines = surf->lines;
  Surf::Point *pts = surf->pts;

  cornerflag[0] = cornerflag[1] = cornerflag[2] = cornerflag[3] = 
    CORNERUNKNOWN;

  cpool.reset();
  cpts.reset();
  opts.grow(pl.n);

  cpt = cpool.get();
  cpt->flag = 0;
  cpt->dot = 0.0;
  cpt->ipl = -1;
  cpt->oindex = -1;
  cpt->x[0] = lo[0]; cpt->x[1] = lo[1];
  cpts.append(cpt);
  cindex[0] = cpt;

  cpt = cpool.get();
  cpt->flag = 0;
  cpt->dot = 0.0;
  cpt->ipl = -1;
  cpt->oindex = -1;
  cpt->x[0] = hi[0]; cpt->x[1] = lo[1];
  cpts.append(cpt);
  cindex[1] = cpt;

  cpt = cpool.get();
  cpt->flag = 0;
  cpt->dot = 0.0;
  cpt->ipl = -1;
  cpt->oindex = -1;
  cpt->x[0] = hi[0]; cpt->x[1] = hi[1];
  cpts.append(cpt);
  cindex[2] = cpt;

  cpt = cpool.get();
  cpt->flag = 0;
  cpt->dot = 0.0;
  cpt->ipl = -1;
  cpt->oindex = -1;
  cpt->x[0] = lo[0]; cpt->x[1] = hi[1];
  cpts.append(cpt);
  cindex[3] = cpt;

  cpt = cpool.get();
  cpt->flag = 0;
  cpt->dot = 0.0;
  cpt->ipl = -1;
  cpt->oindex = -1;
  cpt->x[0] = lo[0]; cpt->x[1] = lo[1];
  cpts.append(cpt);
  cindex[4] = cpt;

  int anyinterior = 0;

  for (int ipl = 0; ipl < pl.n; ipl++) {
    polyline = &pl[ipl];
    MyVec<Opt> &one = opts[ipl];
    one.grow(3*polyline->n + 1);
    nopt = 0;
    
    iline = polyline->first->iline;
    ipt = lines[iline].p1;
    opt = &one[nopt++];
    opt->flag = 0;
    opt->x[0] = pts[ipt].x[0];
    opt->x[1] = pts[ipt].x[1];
    opt->cindex = iline;

    flag = ptflag(pts[ipt].x,lo,hi);
    if (flag == POINTINSIDE) {
      lastflag = INSIDE;
      anyinterior = 1;
    } else lastflag = OUTSIDE;
    
    for (poly = polyline->first; poly; poly = poly->next) {
      iline = poly->iline;
      p = pts[lines[iline].p1].x;
      q = pts[lines[iline].p2].x;
      clip(p,q,lo,hi,a,b);
      if (VERBOSE) printf("CLIP %d %g %g %g %g\n",ipl,a[0],a[1],b[0],b[1]);

      ic = corner(a,lo,hi);
      if (ic >= 0) cornerflag[ic] = CORNEROVERLAP;
      ic = corner(b,lo,hi);
      if (ic >= 0) cornerflag[ic] = CORNEROVERLAP;

      bflag = ptflag(b,lo,hi);
      
      // starting from INSIDE, add exit pt if hit border
      
      if (lastflag == INSIDE) {
        if (bflag == POINTBORDER) {
          anyinterior = 1;
          norm = lines[iline].norm;
          interleave(b,norm,lo,hi,EXIT,ipl,nopt);
          opt = &one[nopt++];
          opt->flag = EXIT;
          opt->x[0] = b[0];
          opt->x[1] = b[1];
          opt->cindex = iline;
          lastflag = OUTSIDE;
        }

      // starting from OUTSIDE, add entry pt on entry/exit
          
      } else {
        if (bflag == POINTINSIDE) {
          anyinterior = 1;
          norm = lines[iline].norm;
          interleave(a,norm,lo,hi,ENTRY,ipl,nopt);
          opt = &one[nopt++];
          opt->flag = ENTRY;
          opt->x[0] = a[0];
          opt->x[1] = a[1];
          opt->cindex = iline;
          lastflag = INSIDE;
        } else {
          if (!sameborder(a,b,lo,hi)) {
            anyinterior = 1;
            norm = lines[iline].norm;
            interleave(a,norm,lo,hi,ENTRY,ipl,nopt);
            opt = &one[nopt++];
            opt->flag = ENTRY;
            opt->x[0] = a[0];
            opt->x[1] = a[1];
            opt->cindex = iline;
            interleave(b,norm,lo,hi,EXIT,ipl,nopt);
            opt = &one[nopt++];
            opt->flag = EXIT;
            opt->x[0] = b[0];
            opt->x[1] = b[1];
            opt->cindex = iline;
          }
        }
      }

      opt = &one[nopt++];
      opt->flag = 0;
      opt->x[0] = q[0];
      opt->x[1] = q[1];
      opt->cindex = iline;
    }
    
    one.n = nopt;

    // check for matching entry/exit pts in the list

    nentry = nexit = 0;
    
    for (i = 0; i < nopt; i++) {
      if (one[i].flag == ENTRY) nentry++;
      if (one[i].flag == EXIT) nexit++;
    }
    if (nentry != nexit) error->one(FLERR,"OPTS entry/exit mis-match");
  }
  
  opts.n = pl.n;

  if (VERBOSE) {
    printf("CORNER FLAGS %d: %d %d %d %d\n",grid->cells[icell].id,
           cornerflag[0],cornerflag[1],cornerflag[2],cornerflag[3]);
    printf("CPTS %d\n",cpts.n);
    int nn = 0;
    for (Cpt *c = cpts.first; c; c = c->next) {
      printf("  %d: %d %g %g %g %d %d\n",nn,c->flag,c->x[0],c->x[1],
             c->dot,c->ipl,c->oindex);
      nn++;
    }
    printf("OPTS %d\n",opts.n);
    for (i = 0; i < opts.n; i++) {
      printf("  ONE %d\n",opts[i].n);
      int nn = 0;
      for (int j = 0; j < opts[i].n; j++) {
        printf("  %d: %d %g %g %d\n",nn,opts[i][j].flag,
               opts[i][j].x[0],opts[i][j].x[1],opts[i][j].cindex);
        nn++;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   interleave an intersection pt into CPTS
   flag = ENTRY/EXIT
   ipl,oindex = which polyline and index within polyline of the new pt
------------------------------------------------------------------------- */

void Cut2d::interleave(double *pt, double *norm, double *lo, double *hi,
                       int flag, int ipl, int oindex)
{
  int dim,order;
  double cnorm[2],cvalue;
  Cpt *start,*stop;
  Cpt *cpt;

  if (pt[1] == lo[1] && pt[0] < hi[0]) {
    start = cindex[0];
    stop = cindex[1];
    cnorm[0] = 0.0;
    cnorm[1] = 1.0;
    dim = 0;
    order = ASCEND;
  } else if (pt[0] == hi[0] && pt[1] < hi[1]) {
    start = cindex[1];
    stop = cindex[2];
    cnorm[0] = -1.0;
    cnorm[1] = 0.0;
    dim = 1;
    order = ASCEND;
  } else if (pt[1] == hi[1] && pt[0] > lo[0]) {
    start = cindex[2];
    stop = cindex[3];
    cnorm[0] = 0.0;
    cnorm[1] = -1.0;
    dim = 0;
    order = DESCEND;
  } else {
    start = cindex[3];
    stop = cindex[4];
    cnorm[0] = 1.0;
    cnorm[1] = 0.0;
    dim = 1;
    order = DESCEND;
  }

  double dot = norm[0]*cnorm[0] + norm[1]*cnorm[1];
  if (flag == EXIT) dot = -dot;
  double value = pt[dim];

  cpt = start;
  for (cpt = start; cpt; cpt = cpt->next) {
    cvalue = cpt->x[dim];
    // values are different, order by ASCEND or DESCEND
    if (order == ASCEND) {
      if (value < cvalue-TOLERANCE) break;
      if (value > cvalue+TOLERANCE) continue;
    } else {
      if (value > cvalue-TOLERANCE) break;
      if (value < cvalue+TOLERANCE) continue;
    }
    // values are equal (within TOLERANCE), corner pt comes first
    if (cpt->flag == 0) continue;
    // values are equal, order by dot product (-1 to 1)
    // dot = line-norm dot cell-edge-norm for entry pt, negative for exit pt
    if (dot < cpt->dot-TOLERANCE) break;
    if (dot > cpt->dot+TOLERANCE) continue;
    // 2 dot products are equal (within TOLERANCE)
    // entry pt comes before exit pt
    if (flag == ENTRY) break;
    // end of loop
    if (cpt == stop) break;
  }

  Cpt *cptnew = cpool.get();
  cptnew->flag = flag;
  cptnew->x[0] = pt[0];
  cptnew->x[1] = pt[1];
  cptnew->dot = dot;
  cptnew->ipl = ipl;
  cptnew->oindex = oindex;
  cpts.insert(cptnew,cpt->prev,cpt);
}
  
/* ----------------------------------------------------------------------
   walk the dual Weiler-Atherton data structures to find loops
   startpts = list of start pts for loop walks in all OPTS
     either an entry pt, or 1st pt of OPT that has no entry pts
   loop = LOOPINSIDE/LOOPBORDER, area, list of line indices in loop
------------------------------------------------------------------------- */

void Cut2d::weiler_walk(double *lo, double *hi, int *cornerflag)
{
  int i,ifirst,iopt,ioptfirst;
  int nlist,done,noexit,ic;
  double area;
  Entrypt *entry,*last,*ept;
  Opt *opt;
  double *p0,*p1;
  Cpt *cpt,*cptprev;

  epool.reset();
  entrypts.reset();
  for (int iopt = 0; iopt < opts.n; iopt++) {
    last = entrypts.last;
    MyVec<Opt> &one = opts[iopt];
    for (i = 0; i < one.n; i++)
      if (one[i].flag == ENTRY) {
        entry = epool.get();
        entry->iopt = iopt;
        entry->index = i;
        entrypts.append(entry);
      }
    if (entrypts.last == last && ptflag(one[0].x,lo,hi) == POINTINSIDE) {
      entry = epool.get();
      entry->iopt = iopt;
      entry->index = 0;
      entrypts.append(entry);
    }
  }

  int nloop = 0;

  while (entrypts.first) {
    iopt = ioptfirst = entrypts.first->iopt;
    i = ifirst = entrypts.first->index;

    // toggle walk between OPTS and CPTS until return to start pt
    // circular walk within OPTS until reach an exit pt or start pt
    // circular walk within CPTS until reach an entry pt
    // OPTS or CPTS edge contributes to area
    // only OPTS edge is added to list
    // if walk corner pt in CPTS:
    //   mark it CORNEROUTSIDE, if not already CORNEROVERLAP

    nloop++;
    loops.grow(nloop);
    Loop &loop = loops[nloop-1];
    MyVec<int> &lines = loop.lines;
    lines.grow(nsurf);
    nlist = 0;

    done = 0;
    noexit = 1;
    area = 0.0;

    while (1) {
      for (ept = entrypts.first; ept; ept = ept->next)
        if (ept->iopt == iopt && ept->index == i) break;
      if (!ept) {
        //print "BAD CELL:",cell[0],cell[1:]
        //print "BAD CPTS",cpts
        //print "BAD OPTS",opts
        error->one(FLERR,"WA walk error");
      }
      entrypts.remove(ept);

      while (1) {
        i++;
        if (i == opts[iopt].n) i = 0;
        if (i == ifirst && iopt == ioptfirst) {
          done = 1;
          break;
        }
        opt = &opts[iopt][i];
        
        p0 = opts[iopt][i-1].x;
        p1 = opt->x;
        if (p0[0] < p1[0])
          area -= (0.5*(p0[1]+p1[1]) - lo[1]) * (p1[0]-p0[0]);
        else
          area += (0.5*(p0[1]+p1[1]) - lo[1]) * (p0[0]-p1[0]);
        if (!nlist || lines[nlist-1] != opt->cindex)
          lines[nlist++] = opt->cindex;

        if (opt->flag == EXIT) {
          noexit = 0;
          break;
        }
      }

      if (done) break;

      for (cpt = cpts.first; cpt; cpt = cpt->next)
        if (iopt == cpt->ipl && i == cpt->oindex) break;

      while (1) {
        cptprev = cpt;
        cpt = cpt->next;
        if (!cpt) cpt = cpts.first;

        if (cpt->flag == 0) {
          ic = corner(cpt->x,lo,hi);
          if (cornerflag[ic] != CORNEROVERLAP) cornerflag[ic] = CORNEROUTSIDE;
        }
        p0 = cptprev->x;
        p1 = cpt->x;
        if (p0[0] < p1[0])
          area -= (0.5*(p0[1]+p1[1]) - lo[1]) * (p1[0]-p0[0]);
        else
          area += (0.5*(p0[1]+p1[1]) - lo[1]) * (p0[0]-p1[0]);

        if (cpt->flag == ENTRY) {
          iopt = cpt->ipl;
          i = cpt->oindex;
          break;
        }
      }
        
      if (i == ifirst && iopt == ioptfirst) break;
    }

    lines.n = nlist;
    if (noexit) loop.flag = LOOPINSIDE;
    else loop.flag = LOOPBORDER;
    loop.area = area;
  }
  
  loops.n = nloop;

  if (VERBOSE) {
    printf("LOOPS %d\n",loops.n);
    for (i = 0; i < loops.n; i++) {
      printf("  %d: %d %g:",i,loops[0].flag,loops[i].area);
      for (int j = 0; j < loops[i].lines.n; j++)
        printf(" %d",loops[i].lines[j]);
      printf("\n");
    }
  }
}

/* ----------------------------------------------------------------------
   convert loops into PG = split cells by combining loops if needed
   mark CORNERUNKNOWN corner flags if possible
   treat zero-area loop as a negative area (unique to 2d)
------------------------------------------------------------------------- */

void Cut2d::loop2pg(double *lo, double *hi, int *cornerflag)
{
  int i,j,unknown,nlines;

  // if any loops are CELLBORDER, remaining CORNERUNKNOWN pts are CORNERINSIDE
  // else if there are interior loops:
  //   if all areas are negative, remaining CORNERUNKNOWN pts are CORNEROUTSIDE
  //   if all areas are positive, remaining CORNERUNKNOWN pts are CORNERINSIDE

  int border = 0;
  int interior = 0;
  int positive = 0;
  int negative = 0;

  for (i = 0; i < loops.n; i++) {
    if (loops[i].flag == LOOPBORDER) border++;
    else if (loops[i].flag == LOOPINSIDE) interior++;
    if (loops[i].area > 0.0) positive++;
    else if (loops[i].area <= 0.0) negative++;
  }
  
  if (border) {
    unknown = 0;
    for (i = 0; i < 4; i++)
      if (cornerflag[i] == CORNERUNKNOWN) unknown = 1;
    if (unknown)
      for (i = 0; i < 4; i++)
        if (cornerflag[i] == CORNERUNKNOWN) cornerflag[i] = CORNERINSIDE;
  } else if (interior) {
    if (positive == 0)
      for (i = 0; i < 4; i++)
        if (cornerflag[i] == CORNERUNKNOWN) cornerflag[i] = CORNEROUTSIDE;
    if (negative == 0)
      for (i = 0; i < 4; i++)
        if (cornerflag[i] == CORNERUNKNOWN) cornerflag[i] = CORNERINSIDE;
  }  

  // if no positive areas:
  //   1 split cell = cell + all loops
  //   area = cell area + sum of all loop areas (negative)
  // if 1 positive area:
  //   1 split cell = positive loop + all negative loops
  //   area = positive area + sum of all loop areas (negative)
  // if multiple positive areas:
  //   each positive-area loop is a PG
  //   allow no negative areas (hard to compute which loop they are inside)
  
  if (positive <= 1) {
    pg.grow(1);
    MyVec<int> &lines = pg[0].lines;
    lines.grow(nsurf);
    nlines = 0;
    pg[0].area = 0.0;
    if (positive == 0) 
      pg[0].area += (hi[0]-lo[0]) * (hi[1]-lo[1]);
    for (i = 0; i < loops.n; i++) {
      pg[0].area += loops[i].area;
      for (j = 0; j < loops[i].lines.n; j++)
        lines[nlines++] = loops[i].lines[j];
    }
    if (pg[0].area <= 0.0)
      error->one(FLERR,"Single area is negative, inverse donut");
    lines.n = nlines;
    pg.n = 1;

  } else {
    if (negative) {
      //areas = [one[1] for one in loop];
      //print "Mixed areas",cell[0],areas;
      error->one(FLERR,"More than one positive area with a negative area");
    }
    pg.grow(loops.n);
    for (i = 0; i < loops.n; i++) {
      pg[i].area = loops[i].area;
      MyVec<int> &lines = pg[i].lines;
      lines.grow(loops[i].lines.n);
      for (j = 0; j < loops[i].lines.n; j++)
        lines[j] = loops[i].lines[j];
      lines.n = loops[i].lines.n;
    }
    pg.n = loops.n;
  }

  if (VERBOSE) {
    printf("PG %d\n",pg.n);
    for (i = 0; i < pg.n; i++) {
      printf("  %d: %g:",i,pg[i].area);
      for (j = 0; j < pg[i].lines.n; j++)
        printf(" %d",pg[i].lines[j]);
      printf("\n");
    }
  }
}

/* ----------------------------------------------------------------------
   assign each line index in list to one of the split cells in PG
   set smap[i] = which PG the Ith line index is assigned to
   smap[i] = -1 if the line segment did not end up in a PG
------------------------------------------------------------------------- */
        
void Cut2d::surf2pg(int n, int *list, int *smap)
{
  int i,j,k,m;
  
  for (i = 0; i < n; i++) smap[i] = -1;
  for (i = 0; i < pg.n; i++) {
    MyVec<int> &lines = pg[i].lines;
    for (j = 0; j < lines.n; j++) {
      m = lines[j];
      for (k = 0; k < n; k++)
        if (m == list[k]) break;
      if (k == n) error->one(FLERR,"Invalid line in PG for mapping to slist");
      smap[k] = i;
    }
  }

  // NOTE: could do this more efficiently if stored list index
  // in PG, rather than original line index
}

/* ----------------------------------------------------------------------
   clip test of line segment PQ against cell with corners LO,HI
   return 1 if intersects, 0 if not
------------------------------------------------------------------------- */

int Cut2d::cliptest(double *p, double *q, double *lo, double *hi)
{
  if (p[0] >= lo[0] && p[0] <= hi[0] &&
      p[1] >= lo[1] && p[1] <= hi[1]) return 1;
  if (q[0] >= lo[0] && q[0] <= hi[0] &&
      q[1] >= lo[1] && q[1] <= hi[1]) return 1;

  double x,y;
  double a[2],b[2];
  a[0] = p[0]; a[1] = p[1];
  b[0] = q[0]; b[1] = q[1];

  if (a[0] < lo[0] && b[0] < lo[0]) return 0;
  if (a[0] < lo[0] || b[0] < lo[0]) {
    y = a[1] + (lo[0]-a[0])/(b[0]-a[0])*(b[1]-a[1]);
    if (a[0] < lo[0]) {
      a[0] = lo[0]; a[1] = y;
    } else {
      b[0] = lo[0]; b[1] = y;
    }
  }
  if (a[0] > hi[0] && b[0] > hi[0]) return 0;
  if (a[0] > hi[0] || b[0] > hi[0]) {
    y = a[1] + (hi[0]-a[0])/(b[0]-a[0])*(b[1]-a[1]);
    if (a[0] > hi[0]) {
      a[0] = hi[0]; a[1] = y;
    } else {
      b[0] = hi[0]; b[1] = y;
    }
  }
  if (a[1] < lo[1] && b[1] < lo[1]) return 0;
  if (a[1] < lo[1] || b[1] < lo[1]) {
    x = a[0] + (lo[1]-a[1])/(b[1]-a[1])*(b[0]-a[0]);
    if (a[1] < lo[1]) {
      a[0] = x; a[1] = lo[1];
    } else {
      b[0] = x; b[1] = lo[1];
    }
  }
  if (a[1] > hi[1] && b[1] > hi[1]) return 0;
  if (a[1] > hi[1] || b[1] > hi[1]) {
    x = a[0] + (hi[1]-a[1])/(b[1]-a[1])*(b[0]-a[0]);
    if (a[1] > hi[1]) {
      a[0] = x; a[1] = hi[1];
    } else {
      b[0] = x; b[1] = hi[1];
    }
  }
     
  return 1;
}

/* ----------------------------------------------------------------------
   clip test line segment PQ against cell with corners LO,HI
   PQ is known to intersect cell
   return AB = clipped segment
------------------------------------------------------------------------- */

void Cut2d::clip(double *p, double *q, double *lo, double *hi,
                double *a, double *b)
{
  double x,y;

  a[0] = p[0]; a[1] = p[1];
  b[0] = q[0]; b[1] = q[1];

  if (p[0] >= lo[0] && p[0] <= hi[0] &&
      p[1] >= lo[1] && p[1] <= hi[1] &&
      q[0] >= lo[0] && q[0] <= hi[0] &&
      q[1] >= lo[1] && q[1] <= hi[1]) return;

  if (a[0] < lo[0] || b[0] < lo[0]) {
    y = a[1] + (lo[0]-a[0])/(b[0]-a[0])*(b[1]-a[1]);
    if (a[0] < lo[0]) {
      a[0] = lo[0]; a[1] = y;
    } else {
      b[0] = lo[0]; b[1] = y;
    }
  }
  if (a[0] > hi[0] || b[0] > hi[0]) {
    y = a[1] + (hi[0]-a[0])/(b[0]-a[0])*(b[1]-a[1]);
    if (a[0] > hi[0]) {
      a[0] = hi[0]; a[1] = y;
    } else {
      b[0] = hi[0]; b[1] = y;
    }
  }
  if (a[1] < lo[1] || b[1] < lo[1]) {
    x = a[0] + (lo[1]-a[1])/(b[1]-a[1])*(b[0]-a[0]);
    if (a[1] < lo[1]) {
      a[0] = x; a[1] = lo[1];
    } else {
      b[0] = x; b[1] = lo[1];
    }
  }
  if (a[1] > hi[1] || b[1] > hi[1]) {
    x = a[0] + (hi[1]-a[1])/(b[1]-a[1])*(b[0]-a[0]);
    if (a[1] > hi[1]) {
      a[0] = x; a[1] = hi[1];
    } else {
      b[0] = x; b[1] = hi[1];
    }
  }
}

/* ----------------------------------------------------------------------
   return 1 if pt1 and pt2 are on same border of cell, 0 if not
   assume pt1 and pt2 are both border pts
------------------------------------------------------------------------- */
  
int Cut2d::sameborder(double *pt1, double *pt2, double *lo, double *hi)
{
  if (pt1[0] == lo[0] and pt2[0] == lo[0]) return 1;
  if (pt1[0] == hi[0] and pt2[0] == hi[0]) return 1;
  if (pt1[1] == lo[1] and pt2[1] == lo[1]) return 1;
  if (pt1[1] == hi[1] and pt2[1] == hi[1]) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   check if pt X,Y is a corner pt of cell
   return 0,1,2,3 = LL,LR,UL,UR if it is, -1 if not
------------------------------------------------------------------------- */

int Cut2d::corner(double *pt, double *lo, double *hi)
{
  if (pt[1] == lo[1]) {
    if (pt[0] == lo[0]) return 0;
    else if (pt[0] == hi[0]) return 1;
  } else if (pt[1] == hi[1]) {
    if (pt[0] == lo[0]) return 2;
    else if (pt[0] == hi[0]) return 3;
  }
  return -1;
}

/* ----------------------------------------------------------------------
   check if pt is inside or outside or on cell border
   return POINTOUTSIDE,POINTINSIDE,POINTBORDER
------------------------------------------------------------------------- */

int Cut2d::ptflag(double *pt, double *lo, double *hi) 
{
  if (pt[0] < lo[0] || pt[0] > hi[0] || pt[1] < lo[1] || pt[1] > hi[1])
    return POINTOUTSIDE;
  if (pt[0] > lo[0] && pt[0] < hi[0] && pt[1] > lo[1] && pt[1] < hi[1])
    return POINTINSIDE;
  return POINTBORDER;
}

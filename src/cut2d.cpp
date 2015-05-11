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

#include "math.h"
#include "string.h"
#include "cut2d.h"
#include "surf.h"
#include "math_const.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};     // several files
enum{EXTERIOR,INTERIOR,BORDER,INTBORD};
enum{ENTRY,EXIT,TWO,CORNER};              // same as Cut3d

#define EPSCELL 1.0e-10    // tolerance for pushing surf pts to cell surface

// cell ID for 2d or 3d cell

//#define VERBOSE
#define VERBOSE_ID 109084

/* ---------------------------------------------------------------------- */

Cut2d::Cut2d(SPARTA *sparta, int caller_axisymmetric) : Pointers(sparta)
{
  axisymmetric = caller_axisymmetric;

  pushflag = surf->pushflag;
  npush = 0;
}

/* ----------------------------------------------------------------------
   compute intersections of a grid cell with all surfs
   csurfs = indices into global surf list
   return nsurf = # of surfs
   return -1 if nsurf > max
------------------------------------------------------------------------- */

int Cut2d::surf2grid(cellint id_caller, double *lo_caller, double *hi_caller, 
                     int *surfs_caller, int max)
{
  id = id_caller;
  lo = lo_caller;
  hi = hi_caller;
  surfs = surfs_caller;

  Surf::Point *pts = surf->pts;
  Surf::Line *lines = surf->lines;
  int nline = surf->nline;

  double *x1,*x2;

  nsurf = 0;
  for (int m = 0; m < nline; m++) {
    x1 = pts[lines[m].p1].x;
    x2 = pts[lines[m].p2].x;

    if (MAX(x1[0],x2[0]) < lo[0]) continue;
    if (MIN(x1[0],x2[0]) > hi[0]) continue;
    if (MAX(x1[1],x2[1]) < lo[1]) continue;
    if (MIN(x1[1],x2[1]) > hi[1]) continue;

    if (cliptest(x1,x2)) {
      if (nsurf == max) return -1;
      surfs[nsurf++] = m;
    }
  }

  return nsurf;
}

/* ----------------------------------------------------------------------
   clip test of line segment PQ against cell with corners LO,HI
   return 1 if intersects, 0 if not
------------------------------------------------------------------------- */

int Cut2d::cliptest(double *p, double *q)
{
  double x,y;

  if (p[0] >= lo[0] && p[0] <= hi[0] &&
      p[1] >= lo[1] && p[1] <= hi[1]) return 1;
  if (q[0] >= lo[0] && q[0] <= hi[0] &&
      q[1] >= lo[1] && q[1] <= hi[1]) return 1;

  double a[2],b[2];
  a[0] = p[0]; a[1] = p[1];
  b[0] = q[0]; b[1] = q[1];

  // if requested, push line pts from just outside cell surface to surface
  // must be done now so to insure a line slightly outside cell 
  //   does not actually intersect cell when pushed in build_clines(),
  //   else might not be added to list of lines intersecting cell

  if (pushflag % 2) {
    push_to_cell(a);
    push_to_cell(b);
  }

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
   clip line segment PQ against cell with corners CLO,CHI
   line may or may not intersect cell (due to rounding)
   return # of clipped points, can be 0,1,2
   return clipped points in cpath as series of x,y pairs
   called externally, depends on no class variables
   no push of line points is done, since external caller uses true surface
   duplicate points in cpath are deleted
------------------------------------------------------------------------- */

int Cut2d::clip_external(double *p, double *q, double *clo, double *chi, 
                         double *cpath)
{
  double x,y;

  // PQ is interior to cell

  if (p[0] >= clo[0] && p[0] <= chi[0] &&
      p[1] >= clo[1] && p[1] <= chi[1] &&
      q[0] >= clo[0] && q[0] <= chi[0] &&
      q[1] >= clo[1] && q[1] <= chi[1]) {
    cpath[0] = p[0];
    cpath[1] = p[1];
    cpath[2] = q[0];
    cpath[3] = q[1];
    return 2;
  }

  double a[2],b[2];
  a[0] = p[0]; a[1] = p[1];
  b[0] = q[0]; b[1] = q[1];

  if (a[0] < clo[0] && b[0] < clo[0]) return 0;
  if (a[0] < clo[0] || b[0] < clo[0]) {
    y = a[1] + (clo[0]-a[0])/(b[0]-a[0])*(b[1]-a[1]);
    if (a[0] < clo[0]) {
      a[0] = clo[0]; a[1] = y;
    } else {
      b[0] = clo[0]; b[1] = y;
    }
  }
  if (a[0] > chi[0] && b[0] > chi[0]) return 0;
  if (a[0] > chi[0] || b[0] > chi[0]) {
    y = a[1] + (chi[0]-a[0])/(b[0]-a[0])*(b[1]-a[1]);
    if (a[0] > chi[0]) {
      a[0] = chi[0]; a[1] = y;
    } else {
      b[0] = chi[0]; b[1] = y;
    }
  }
  if (a[1] < clo[1] && b[1] < clo[1]) return 0;
  if (a[1] < clo[1] || b[1] < clo[1]) {
    x = a[0] + (clo[1]-a[1])/(b[1]-a[1])*(b[0]-a[0]);
    if (a[1] < clo[1]) {
      a[0] = x; a[1] = clo[1];
    } else {
      b[0] = x; b[1] = clo[1];
    }
  }
  if (a[1] > chi[1] && b[1] > chi[1]) return 0;
  if (a[1] > chi[1] || b[1] > chi[1]) {
    x = a[0] + (chi[1]-a[1])/(b[1]-a[1])*(b[0]-a[0]);
    if (a[1] > chi[1]) {
      a[0] = x; a[1] = chi[1];
    } else {
      b[0] = x; b[1] = chi[1];
    }
  }

  cpath[0] = a[0];
  cpath[1] = a[1];
  cpath[2] = b[0];
  cpath[3] = b[1];

  if (a[0] == b[0] && a[1] == b[1]) return 1;
  return 2;
}

/* ----------------------------------------------------------------------
   calculate cut area of a grid cell that contains nsurf lines
   also calculate if cell is split into distinct sub-areas by lines
   return nsplit = # of splits, 1 for no split
   return areas = ptr to vector of areas = one area per split
     if nsplit = 1, cut area
     if nsplit > 1, one area per split cell
   if nsplit > 1, also return:
     surfmap = which sub-cell (0 to Nsurfs-1) each surf is in
             = -1 if not in any sub-cell (i.e. discarded from clines)
     xsub = which sub-cell (0 to Nsplit-1) xsplit is in
     xsplit = coords of a point in the split cell
------------------------------------------------------------------------- */

int Cut2d::split(cellint id_caller, double *lo_caller, double *hi_caller, 
                 int nsurf_caller, int *surfs_caller,
                 double *&areas_caller, int *surfmap, 
                 int *corners, int &xsub, double *xsplit)
{
  id = id_caller;
  lo = lo_caller;
  hi = hi_caller;
  nsurf = nsurf_caller;
  surfs = surfs_caller;

  int grazeflag = build_clines();

  if (clines.n == 0) {
    areas.grow(1);
    areas[0] = 0.0;
    if (grazeflag) corners[0] = corners[1] = corners[2] = corners[3] = INSIDE;
    areas_caller = &areas[0];
    return 1;
  }

  weiler_build();
  weiler_loops();
  loop2pg();
  
  int nsplit = pgs.n;
  if (nsplit > 1) {
    create_surfmap(surfmap);
    xsub = split_point(surfmap,xsplit);
  }
  
  // set corners = OUTSIDE if corner pt is in list of points in PGs
  // else set corners = INSIDE

  corners[0] = corners[1] = corners[2] = corners[3] = INSIDE;

  int iloop,nloop,mloop,ipt,npt,mpt;

  int npg = pgs.n;
  for (int ipg = 0; ipg < npg; ipg++) {
    nloop = pgs[ipg].n;
    mloop = pgs[ipg].first;
    for (iloop = 0; iloop < nloop; iloop++) {
      npt = loops[mloop].n;
      mpt = loops[mloop].first;
      for (ipt = 0; ipt < npt; ipt++) {
        if (points[mpt].corner >= 0)
          corners[points[mpt].corner] = OUTSIDE;
        mpt = points[mpt].next;
      }
      mloop = loops[mloop].next;
    }
  }

  // store areas in vector so can return ptr to it

  areas.grow(nsplit);
  for (int i = 0; i < nsplit; i++) areas[i] = pgs[i].area;
  areas_caller = &areas[0];

  return nsplit;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void Cut2d::split_face(int id_caller, int iface, double *onelo, double *onehi)
{
  id = id_caller;
  lo = onelo;
  hi = onehi;

#ifdef VERBOSE
  if (id == VERBOSE_ID) {
    printf("Clines for cell %d, face %d\n",id,iface);
    print_clines();
  }
#endif

  weiler_build();

#ifdef VERBOSE
  if (id == VERBOSE_ID) {
    printf("Points for cell %d, face %d\n",id,iface);
    print_points();
  }
#endif

  weiler_loops();
  loop2pg();

#ifdef VERBOSE
  if (id == VERBOSE_ID) {
    printf("Loops for cell %d, face %d\n",id,iface);
    print_loops();
  }
#endif
}

/* ----------------------------------------------------------------------
   create clines = list of lines clipped to cell
------------------------------------------------------------------------- */

int Cut2d::build_clines()
{
  int m;
  double p1[2],p2[2];
  double *x,*y,*norm;
  Surf::Line *line;
  Cline *cline;

  Surf::Point *pts = surf->pts;
  Surf::Line *lines = surf->lines;

  clines.grow(nsurf);
  clines.n = 0;

  int grazeflag = 0;
  int n = 0;

  for (int i = 0; i < nsurf; i++) {
    m = surfs[i];
    line = &lines[m];
    memcpy(p1,pts[line->p1].x,2*sizeof(double));
    memcpy(p2,pts[line->p2].x,2*sizeof(double));

    if (pushflag) {
      npush += push_to_cell(p1);
      npush += push_to_cell(p2);
    }

    cline = &clines[n];
    cline->line = i;

    // clip PQ to cell and store as XY in Cline

    x = cline->x;
    y = cline->y;
    clip(p1,p2,x,y);

    // discard clipped line if only one point

    if (x[0] == y[0] && x[1] == y[1]) continue;
    
    // discard clipped line if lies along one cell edge
    //   and normal is not into cell
    // leave grazeflag incremented in this case

    grazeflag++;
    if (ptflag(x) == BORDER && ptflag(y) == BORDER) {
      int edge = sameedge(x,y);
      if (edge) {
        norm = line->norm;
        if (edge == 1 and norm[0] < 0.0) continue;
        if (edge == 2 and norm[0] > 0.0) continue;
        if (edge == 3 and norm[1] < 0.0) continue;
        if (edge == 4 and norm[1] > 0.0) continue;
        grazeflag--;
      }
    }

    n++;
  }

  clines.n = n;
  return grazeflag;
}

/* ----------------------------------------------------------------------
   create Weiler/Atherton data structure within points
   add each CLINE, with next walking from ENTRY pt to EXIT pt
   check that CLINES has valid set of points
     no pt should appear no more than once as first or more than once as last
     every singlet pt should be on cell border
   add 4 corner pts to points
     new CORNER pt if doesn't exist
     else already an ENTRY/EXIT pt (not a TWO pt)
   create counter-clockwise linked list of cell perimeter pts
     use cprev,cnext,side,value fields of points
   set next field for points in linked list to complete loops
     do not change next for ENTRY pts
------------------------------------------------------------------------- */

void Cut2d::weiler_build()
{
  int i,j;
  int firstpt,lastpt,nextpt;
  double *pt;

  // insure space in points for 2 pts per cline + 4 corner pts

  int nlines = clines.n;
  points.grow(2*nlines + 4);
  points.n = 0;

  // add each cline end pt to points
  // set x,type,next,line for each pt
  // NOTE: could use hash to find existing pts in O(1) time
  //   see Brian Adams code (howto/c++) for doing this on a vec of doubles

  int npt = 0;

  for (i = 0; i < nlines; i++) {

    // 1st point in cline

    pt = clines[i].x;
    for (j = 0; j < npt; j++)
      if (pt[0] == points[j].x[0] && pt[1] == points[j].x[1]) break;

    // pt already exists

    if (j < npt) {
      if (points[j].type == ENTRY || points[j].type == TWO) 
        error->one(FLERR,"Point appears first in more than one CLINE");
      points[j].type = TWO;
      points[j].line = clines[i].line;
      firstpt = j;
    }

    // new pt

    else {
      points[npt].x[0] = pt[0];
      points[npt].x[1] = pt[1];
      points[npt].type = ENTRY;
      points[npt].line = clines[i].line;
      firstpt = npt;
      npt++;
    }

    // 2nd point in cline

    pt = clines[i].y;
    for (j = 0; j < npt; j++)
      if (pt[0] == points[j].x[0] && pt[1] == points[j].x[1]) break;

    // pt already exists

    if (j < npt) {
      if (points[j].type == EXIT || points[j].type == TWO)
        error->one(FLERR,"Point appears last in more than one CLINE");
      points[j].type = TWO;
      points[firstpt].next = j;
    }

    // new pt

    else {
      points[npt].x[0] = pt[0];
      points[npt].x[1] = pt[1];
      points[npt].type = EXIT;
      points[firstpt].next = npt;
      npt++;
    }
  }

  // error check that every singlet point is on cell border

  for (i = 0; i < npt; i++)
    if (points[i].type != TWO && ptflag(points[i].x) != BORDER)
      error->one(FLERR,"Singlet CLINES point not on cell border");

  // add 4 cell CORNER pts to points
  // only if corner pt is not already an ENTRY or EXIT pt
  // if a TWO pt, still add corner pt as CORNER
  // each corner pt is at beginning of side it is on
  // NOTE: could use hash to find existing pts in O(1) time
  // corner flag = 0,1,2,3 for LL,LR,UL,UR, same as in Grid::ChildInfo,
  //   but ordering of corner pts in linked list is LL,LR,UR,UL
  // side = 0,1,2,3 for lower,right,upper,left = traversal order in linked list
  // value = x-coord for lower/upper sides, y-coord for left,right sides

  double cpt[2];
  int ipt1,ipt2,ipt3,ipt4;

  for (i = 0; i < npt; i++) points[i].corner = -1;

  cpt[0] = lo[0]; cpt[1] = lo[1];
  for (j = 0; j < npt; j++)
    if (cpt[0] == points[j].x[0] && cpt[1] == points[j].x[1]) break;
  if (j == npt || points[j].type == TWO) {
    points[npt].x[0] = cpt[0];
    points[npt].x[1] = cpt[1];
    points[npt].type = CORNER;
    ipt1 = npt++;
  } else ipt1 = j;
  points[ipt1].corner = 0;
  points[ipt1].side = 0;
  points[ipt1].value = lo[0];

  cpt[0] = hi[0]; cpt[1] = lo[1];
  for (j = 0; j < npt; j++)
    if (cpt[0] == points[j].x[0] && cpt[1] == points[j].x[1]) break;
  if (j == npt || points[j].type == TWO) {
    points[npt].x[0] = cpt[0];
    points[npt].x[1] = cpt[1];
    points[npt].type = CORNER;
    ipt2 = npt++;
  } else ipt2 = j;
  points[ipt2].corner = 1;
  points[ipt2].side = 1;
  points[ipt2].value = lo[1];

  cpt[0] = hi[0]; cpt[1] = hi[1];
  for (j = 0; j < npt; j++)
    if (cpt[0] == points[j].x[0] && cpt[1] == points[j].x[1]) break;
  if (j == npt || points[j].type == TWO) {
    points[npt].x[0] = cpt[0];
    points[npt].x[1] = cpt[1];
    points[npt].type = CORNER;
    ipt3 = npt++;
  } else ipt3 = j;
  points[ipt3].corner = 3;
  points[ipt3].side = 2;
  points[ipt3].value = hi[0];

  cpt[0] = lo[0]; cpt[1] = hi[1];
  for (j = 0; j < npt; j++)
    if (cpt[0] == points[j].x[0] && cpt[1] == points[j].x[1]) break;
  if (j == npt || points[j].type == TWO) {
    points[npt].x[0] = cpt[0];
    points[npt].x[1] = cpt[1];
    points[npt].type = CORNER;
    ipt4 = npt++;
  } else ipt4 = j;
  points[ipt4].corner = 2;
  points[ipt4].side = 3;
  points[ipt4].value = hi[1];

  // n = final # of points

  points.n = npt;

  // create counter-clockwise linked list of 4 pts around cell perimeter

  firstpt = ipt1;
  lastpt = ipt4;

  points[ipt1].cprev = -1;
  points[ipt1].cnext = ipt2;
  points[ipt2].cprev = ipt1;
  points[ipt2].cnext = ipt3;
  points[ipt3].cprev = ipt2;
  points[ipt3].cnext = ipt4;
  points[ipt4].cprev = ipt3;
  points[ipt4].cnext = -1;

  // add all non-corner ENTRY/EXIT pts to counter-clockwise linked list
  // side = 0,1,2,3 for lower,right,upper,left sides
  // value = coord of pt along the side it is on

  int ipt,iprev,side;
  double value;

  for (i = 0; i < npt; i++) {
    if (points[i].type == TWO || points[i].type == CORNER) continue;
    if (points[i].corner >= 0) continue;

    side = whichside(points[i].x);
    if (side % 2) value = points[i].x[1];
    else value = points[i].x[0];

    // interleave Ith point into linked list between firstpt and lastpt
    // insertion location is between iprev and ipt
    // special logic if inserting at end of list

    ipt = firstpt;
    while (ipt >= 0) {
      if (side < points[ipt].side) break;
      if (side == points[ipt].side) {
        if (side < 2) {
          if (value < points[ipt].value) break;
        } else {
          if (value > points[ipt].value) break;
        }
      }
      iprev = ipt;
      ipt = points[ipt].cnext;
    }

    points[i].side = side;
    points[i].value = value;
    points[i].cprev = iprev;
    points[i].cnext = ipt;

    points[iprev].cnext = i;
    if (ipt >= 0) points[ipt].cprev = i;
    else lastpt = i;
  }

  // set next field for cell perimeter points in linked list
  // this completes loops for next fields already set between ENTRY/EXIT pts
  // do not reset next for ENTRY pts
  // after loop, explicitly connect lastpt to firstpt

  ipt = firstpt;
  while (ipt >= 0) {
    nextpt = points[ipt].cnext;
    if (points[ipt].type != ENTRY) points[ipt].next = nextpt;
    ipt = nextpt;
  }
  if (points[lastpt].type != ENTRY) points[lastpt].next = firstpt;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void Cut2d::weiler_loops()
{
  // used = 0/1 flag for whether a point is already part of a loop

  int n = points.n;
  used.grow(n);
  for (int i = 0; i < n; i++) used[i] = 0;
  used.n = n;

  // iterate over all pts
  // start a loop at any unused pt
  // walk loop via next field:
  //   mark pts as used
  //   compute area along the way
  // stop when reach initial pt:
  //   closed loop, add to loops data structure
  //   loop of just 4 corner pts is a valid loop
  // stop when reach used pt:
  //   discard loop, just traversed non-loop corner pts

  int ipt,iflag,cflag,ncount,firstpt,nextpt;
  double area;
  double *x,*y;

  int nloop = 0;

  for (int i = 0; i < n; i++) {
    if (used[i]) continue;
    area = 0.0;
    iflag = cflag = 1;
    ncount = 0;

    ipt = firstpt = i;
    x = points[ipt].x;

    while (!used[ipt]) {
      used[ipt] = 1;
      ncount++;
      if (points[ipt].type != TWO) iflag = 0;
      if (points[ipt].type != CORNER) cflag = 0;
      nextpt = points[ipt].next;
      y = points[nextpt].x;
      if (axisymmetric)
        area -= MY_PI3 * (x[1]*x[1] + x[1]*y[1] + y[1]*y[1]) * (y[0]-x[0]);
      else area -= (0.5*(x[1]+y[1]) - lo[1]) * (y[0]-x[0]);
      x = y;
      ipt = nextpt;
      if (ipt == firstpt) break;
    }
    if (ipt != firstpt) continue;
    
    loops.grow(nloop+1);
    loops[nloop].area = area;
    loops[nloop].active = 1;
    if (iflag) loops[nloop].flag = INTERIOR;
    else if (cflag) loops[nloop].flag = BORDER;
    else loops[nloop].flag = INTBORD;
    loops[nloop].n = ncount;
    loops[nloop].first = firstpt;
    nloop++;
  }

  loops.n = nloop;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void Cut2d::loop2pg()
{
  int positive = 0;
  int negative = 0;

  int nloop = loops.n;
  for (int i = 0; i < nloop; i++)
    if (loops[i].area > 0.0) positive++;
    else negative++;
  if (positive == 0) error->one(FLERR,"No positive areas in cell");
  if (positive > 1 && negative) 
    error->one(FLERR,"More than one positive area with a negative area");

  pgs.grow(positive);

  // if multiple positive loops, mark BORDER loop as inactive if exists
  // don't want entire cell border to be a loop in this case

  if (positive > 1) {
    for (int i = 0; i < nloop; i++)
      if (loops[i].flag == BORDER) {
        loops[i].active = 0;
        positive--;
      }
  }

  // positive = 1 means 1 PG with area = sum of all pos/neg loops
  // positive > 1 means each loop is a PG

  if (positive == 1) {
    double area = 0.0;
    int prev = -1;
    int count = 0;
    int first;

    for (int i = 0; i < nloop; i++) {
      if (!loops[i].active) continue;
      area += loops[i].area;
      count++;
      if (prev < 0) first = i;
      else loops[prev].next = i;
      prev = i;
    }
    loops[prev].next = -1;

    if (area < 0.0) error->one(FLERR,"Single area is negative, inverse donut");

    pgs[0].area = area;
    pgs[0].n = count;
    pgs[0].first = first;

  } else {
    int m = 0;
    for (int i = 0; i < nloop; i++) {
      if (!loops[i].active) continue;
      pgs[m].area = loops[i].area;
      pgs[m].n = 1;
      pgs[m].first = i;
      m++;
      loops[i].next = -1;
    }
  }

  pgs.n = positive;
}

/* ----------------------------------------------------------------------
   assign each line index in list to one of the split cells in PG
   return surfmap[i] = which PG the Ith line index is assigned to
   set surfmap[i] = -1 if the line did not end up in a PG
     could have been discarded from CLINES
     due to touching cell or lying along a cell edge
------------------------------------------------------------------------- */

void Cut2d::create_surfmap(int *surfmap)
{
  for (int i = 0; i < nsurf; i++) surfmap[i] = -1;

  int iloop,nloop,mloop,ipt,npt,mpt;

  int npg = pgs.n;
  for (int ipg = 0; ipg < npg; ipg++) {
    nloop = pgs[ipg].n;
    mloop = pgs[ipg].first;
    for (iloop = 0; iloop < nloop; iloop++) {
      npt = loops[mloop].n;
      mpt = loops[mloop].first;
      for (ipt = 0; ipt < npt; ipt++) {
        if (points[mpt].type == TWO || points[mpt].type == ENTRY)
          surfmap[points[mpt].line] = ipg;
        mpt = points[mpt].next;
      }
      mloop = loops[mloop].next;
    }
  }
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

int Cut2d::split_point(int *surfmap, double *xsplit)
{
  int iline;
  double *x1,*x2;
  double a[2],b[2];

  Surf::Point *pts = surf->pts;
  Surf::Line *lines = surf->lines;

  // if end pt of any line with non-negative surfmap is in/on cell, return

  for (int i = 0; i < nsurf; i++) {
    if (surfmap[i] < 0) continue;
    iline = surfs[i];
    x1 = pts[lines[iline].p1].x;
    x2 = pts[lines[iline].p2].x;
    if (ptflag(x1) != EXTERIOR) {
      xsplit[0] = x1[0]; xsplit[1] = x1[1];
      return surfmap[i];
    }
    if (ptflag(x2) != EXTERIOR) {
      xsplit[0] = x2[0]; xsplit[1] = x2[1];
      return surfmap[i];
    }
  }

  // clip 1st line with non-negative surfmap to cell, and return clip point

  for (int i = 0; i < nsurf; i++) {
    if (surfmap[i] < 0) continue;
    iline = surfs[i];
    x1 = pts[lines[iline].p1].x;
    x2 = pts[lines[iline].p2].x;
    clip(x1,x2,a,b);
    xsplit[0] = a[0]; xsplit[1] = a[1];
    return surfmap[i];
  }

  error->one(FLERR,"Could not find split point in split cell");
  return -1;
}

/* ----------------------------------------------------------------------
   clip test line segment PQ against cell with corners LO,HI
   PQ is known to intersect cell
   return AB = clipped segment
------------------------------------------------------------------------- */

void Cut2d::clip(double *p, double *q, double *a, double *b)
{
  double x,y;

  a[0] = p[0]; a[1] = p[1];
  b[0] = q[0]; b[1] = q[1];

  if (p[0] >= lo[0] && p[0] <= hi[0] &&
      p[1] >= lo[1] && p[1] <= hi[1] &&
      q[0] >= lo[0] && q[0] <= hi[0] &&
      q[1] >= lo[1] && q[1] <= hi[1]) return;

  // if requested, push line pts from just outside cell surface to surface
  // must be done now so to insure a line slightly outside cell 
  //   does not actually intersect cell when pushed in build_clines(),
  //   else might not be added to list of lines intersecting cell

  if (pushflag % 2) {
    push_to_cell(a);
    push_to_cell(b);
  }

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
   check if pt is inside or outside or on cell border
   return EXTERIOR,BORDER,INTERIOR
------------------------------------------------------------------------- */

int Cut2d::ptflag(double *pt)
{
  double x = pt[0];
  double y = pt[1];
  if (x < lo[0] || x > hi[0] || y < lo[1] || y > hi[1]) return EXTERIOR;
  if (x > lo[0] && x < hi[0] && y > lo[1] && y < hi[1]) return INTERIOR;
  return BORDER;
}

/* ----------------------------------------------------------------------
   check if pt is within EPSCELL of cell surface
   if so, push to cell surface
   for pushflag = 1, can be inside or outside by EPSCELL
   for pushflag = 2, only inside
   for pushflag = 3, only outside
------------------------------------------------------------------------- */

int Cut2d::push_to_cell(double *pt)
{
  double x = pt[0];
  double y = pt[1];
  double epsx = EPSCELL * (hi[0]-lo[0]);
  double epsy = EPSCELL * (hi[1]-lo[1]);

  if (pushflag % 2) {
    if (x < lo[0]-epsx || x > hi[0]+epsx || 
        y < lo[1]-epsy || y > hi[1]+epsy) return 0;
  } else {
    if (x < lo[0] || x > hi[0] || y < lo[1] || y > hi[1]) return 0;
  }

  if (pushflag == 1) {
    if (fabs(x-lo[0]) < epsx || fabs(x-hi[0]) < epsx) {
      if (y > lo[1]-epsy && y < hi[1]+epsy) {
        if (fabs(x-lo[0]) < epsx) x = lo[0];
        else x = hi[0];
      }
    }
    if (fabs(y-lo[1]) < epsy || fabs(y-hi[1]) < epsy) {
      if (x >= lo[0] && x <= hi[0]) {
        if (fabs(y-lo[1]) < epsy) y = lo[1];
        else y = hi[1];
      }
    }

  } else if (pushflag == 2) {
    if ((x >= lo[0] && x-lo[0] < epsx) || (x <= hi[0] && hi[0]-x < epsx)) {
      if (y >= lo[1]-epsy && y <= hi[1]+epsy) {
        if (x-lo[0] < epsx) x = lo[0];
        else x = hi[0];
      }
    }
    if ((y >= lo[1] && y-lo[1] < epsy) || (y <= hi[1] && hi[1]-y < epsy)) {
      if (x >= lo[0] && x <= hi[0]) {
        if (y-lo[1] < epsy) y = lo[1];
        else y = hi[1];
      }
    }

  } else if (pushflag == 3) {
    if ((x > lo[0]-epsx && x <= lo[0]) || (x < hi[0]+epsx && x >= hi[0])) {
      if (y >= lo[1]-epsy && y <= hi[1]+epsy) {
        if (x <= lo[0]) x = lo[0];
        else x = hi[0];
      }
    }
    if ((y > lo[1]-epsy && y <= lo[1]) || (y < hi[1]+epsy && y >= hi[1])) {
      if (x >= lo[0] && x <= hi[0]) {
        if (y <= lo[1]) y = lo[1];
        else y = hi[1];
      }
    }
  }

  int flag = 0;
  if (x != pt[0] || y != pt[1]) {
    flag = 1;
    pt[0] = x;
    pt[1] = y;
  }

  return flag;
}

/* ----------------------------------------------------------------------
   check if pts A,B are on same edge of cell
   both assumed to be on cell border
   return 1,2,3,4 = left,right,lower,upper if they are, 0 if not
------------------------------------------------------------------------- */

int Cut2d::sameedge(double *a, double *b)
{
  if (a[0] == lo[0] and b[0] == lo[0]) return 1;
  if (a[0] == hi[0] and b[0] == hi[0]) return 2;
  if (a[1] == lo[1] and b[1] == lo[1]) return 3;
  if (a[1] == hi[1] and b[1] == hi[1]) return 4;
  return 0;
}

/* ----------------------------------------------------------------------
   check which side of cell pt is on, assumed to be on border
   not called for corner pt, but return either side
   return 0,1,2,3 = lower,right,upper,left if on border, -1 if not
------------------------------------------------------------------------- */

int Cut2d::whichside(double *pt)
{
  if (pt[0] == lo[0]) return 3;
  if (pt[0] == hi[0]) return 1;
  if (pt[1] == lo[1]) return 0;
  if (pt[1] == hi[1]) return 2;
  return -1;
}

/* ----------------------------------------------------------------------
   debug: print out Clines list
------------------------------------------------------------------------- */

void Cut2d::print_clines()
{
  printf("ICELL %d\n",id);
  printf("  clines %d\n",clines.n);

  for (int i = 0; i < clines.n; i++) {
    printf("  line %d\n",i);
    printf("    xpoint: %g %g\n",clines[i].x[0],clines[i].x[1]);
    printf("    ypoint: %g %g\n",clines[i].y[0],clines[i].y[1]);
    printf("    line %d\n",clines[i].line);
  }
}

/* ----------------------------------------------------------------------
   debug: print out Points list
------------------------------------------------------------------------- */

void Cut2d::print_points()
{
  printf("ICELL %d\n",id);
  printf("  npoints %d\n",points.n);

  for (int i = 0; i < points.n; i++) {
    printf("  point %d\n",i);
    printf("    coord: %g %g\n",points[i].x[0],points[i].x[1]);
    printf("    type %d, next %d\n",points[i].type,points[i].next);
    if (points[i].type == ENTRY || points[i].type == TWO)
      printf("    line %d\n",points[i].line);
    else printf("    line -1\n");
    printf("    corner %d\n",points[i].corner);
    if (points[i].type != TWO) {
      printf("    cprev %d, cnext %d\n",points[i].cprev,points[i].cnext);
      printf("    side/value: %d %g\n",points[i].side,points[i].value);
    }
  }
}

/* ----------------------------------------------------------------------
   debug: print out Loops list
------------------------------------------------------------------------- */

void Cut2d::print_loops()
{
  printf("ICELL %d\n",id);
  printf("  loops %d\n",loops.n);

  for (int i = 0; i < loops.n; i++) {
    printf("  loop %d\n",i);
    printf("    active %d\n",loops[i].active);
    printf("    flag %d\n",loops[i].flag);
    printf("    area %g\n",loops[i].area);
    printf("    npoints %d\n",loops[i].n);
    printf("    points:\n");
    int ipt = loops[i].first;
    for (int j = 0; j < loops[i].n; j++) {
      printf("      %d %g %g\n",ipt,points[ipt].x[0],points[ipt].x[1]);
      ipt = points[ipt].next;
    }
  }
}

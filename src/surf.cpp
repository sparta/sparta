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

#include "ctype.h"
#include "surf.h"
#include "style_surf_collide.h"
#include "style_surf_react.h"
#include "surf_collide.h"
#include "surf_react.h"
#include "domain.h"
#include "region.h"
#include "grid.h"
#include "comm.h"
#include "geometry.h"
#include "input.h"
#include "math_extra.h"
#include "math_const.h"
#include "hash3.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{TALLYAUTO,TALLYREDUCE,TALLYRVOUS};         // same as Update
enum{REGION_ALL,REGION_ONE,REGION_CENTER};      // same as Grid
enum{TYPE,MOLECULE,ID};
enum{LT,LE,GT,GE,EQ,NEQ,BETWEEN};

#define DELTA 1024
#define DELTAMODEL 4
#define EPSSQ 1.0e-12
#define EPSILON_GRID 1.0e-3
#define BIG 1.0e20
#define MAXGROUP 32

/* ---------------------------------------------------------------------- */

Surf::Surf(SPARTA *sparta) : Pointers(sparta)
{
  exist = 0;
  implicit = 0;
  distributed = 0;
  surf_collision_check = 1;

  gnames = (char **) memory->smalloc(MAXGROUP*sizeof(char *),"surf:gnames");
  bitmask = (int *) memory->smalloc(MAXGROUP*sizeof(int),"surf:bitmask");
  inversemask = (int *) memory->smalloc(MAXGROUP*sizeof(int),
                                        "surf:inversemask");

  for (int i = 0; i < MAXGROUP; i++) bitmask[i] = 1 << i;
  for (int i = 0; i < MAXGROUP; i++) inversemask[i] = bitmask[i] ^ ~0;

  ngroup = 1;
  int n = strlen("all") + 1;
  gnames[0] = new char[n];
  strcpy(gnames[0],"all");

  nsurf = 0;

  nlocal = nghost = nmax = 0;
  lines = NULL;
  tris = NULL;
  pushflag = 1;

  nown = 0;
  mylines = NULL;
  mytris = NULL;

  nsc = maxsc = 0;
  sc = NULL;
  nsr = maxsr = 0;
  sr = NULL;

  tally_comm = TALLYAUTO;

  // allocate hash for surf IDs

  hash = new MySurfHash();
  hashfilled = 0;
}

/* ---------------------------------------------------------------------- */

Surf::~Surf()
{
  for (int i = 0; i < ngroup; i++) delete [] gnames[i];
  memory->sfree(gnames);
  memory->sfree(bitmask);
  memory->sfree(inversemask);

  memory->sfree(lines);
  memory->sfree(tris);
  memory->sfree(mylines);
  memory->sfree(mytris);

  for (int i = 0; i < nsc; i++) delete sc[i];
  memory->sfree(sc);
  for (int i = 0; i < nsr; i++) delete sr[i];
  memory->sfree(sr);

  hash->clear();
  delete hash;
}

/* ---------------------------------------------------------------------- */

void Surf::global(char *arg)
{
  if (exist) 
    error->all(FLERR,"Cannot set global surfs when surfaces already exist");

  if (strcmp(arg,"explicit") == 0) {
    implicit = 0;
    distributed = 0;
  } else if (strcmp(arg,"explicit/distributed") == 0) {
    implicit = 0;
    distributed = 1;
  } else if (strcmp(arg,"implicit") == 0) {
    implicit = 1;
    distributed = 1;
  } else error->all(FLERR,"Illegal global command");
}

/* ---------------------------------------------------------------------- */

void Surf::modify_params(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal surf_modify command");
  int igroup = find_group(arg[0]);
  if (igroup < 0) error->all(FLERR,"Surf_modify surface group is not defined");
  int groupbit = bitmask[igroup];

  int dim = domain->dimension;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"collide") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal surf_modify command");

      int isc = find_collide(arg[iarg+1]);
      if (isc < 0) error->all(FLERR,"Could not find surf_modify sc-ID");

      // NOTE: is this also needed for mylines and mytris?
      // set surf collision model for each surf in surface group

      if (dim == 2) {
        for (int i = 0; i < nlocal+nghost; i++)
          if (lines[i].mask & groupbit) lines[i].isc = isc;
      }
      if (dim == 3) {
        for (int i = 0; i < nlocal+nghost; i++)
          if (tris[i].mask & groupbit) tris[i].isc = isc;
      }

      iarg += 2;

    } else if (strcmp(arg[iarg],"react") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal surf_modify command");
      
      int isr;
      if (strcmp(arg[iarg+1],"none") == 0) isr = -1;
      else {
        isr = find_react(arg[iarg+1]);
        if (isr < 0) error->all(FLERR,"Could not find surf_modify sr-ID");
      }

      // set surf reaction model for each surf in surface group

      if (dim == 2) {
        for (int i = 0; i < nlocal+nghost; i++)
          if (lines[i].mask & groupbit) lines[i].isr = isr;
      }
      if (dim == 3) {
        for (int i = 0; i < nlocal+nghost; i++)
          if (tris[i].mask & groupbit) tris[i].isr = isr;
      }

      iarg += 2;

    } else error->all(FLERR,"Illegal surf_modify command");
  }
}

/* ---------------------------------------------------------------------- */

void Surf::init()
{
  // warn if surfs are distributed (explicit or implicit)
  //   and grid->cutoff < 0.0, since each proc will have copy of all cells

  if (exist && distributed && grid->cutoff < 0.0)
    if (comm->me == 0) 
      error->warning(FLERR,"Surfs are distributed with infinite grid cutoff");

  // check that every element is assigned to a surf collision model
  // skip if caller turned off the check, e.g. BalanceGrid
  // NOTE: add distributed logic

  int dim = domain->dimension;
  bigint flag,allflag;

  if (surf_collision_check) {
    flag = 0;
    if (dim == 2) {
      for (int i = 0; i < nlocal+nghost; i++)
        if (lines[i].isc < 0) flag++;
    } 
    if (dim == 3) {
      for (int i = 0; i < nlocal+nghost; i++)
        if (tris[i].isc < 0) flag++;
    }

    if (distributed)
      MPI_Allreduce(&flag,&allflag,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
    else allflag = flag;

    if (allflag) {
      char str[64];
      sprintf(str,BIGINT_FORMAT 
              " surface elements not assigned to a collision model",allflag);
      error->all(FLERR,str);
    }
  }

  // if a surf element is assigned a reaction model
  // must have a collision model that allows reactions

  if (surf_collision_check) {
    flag = 0;
    if (dim == 2) {
      for (int i = 0; i < nlocal+nghost; i++)
        if (lines[i].isr >= 0 && sc[lines[i].isc]->allowreact == 0) flag++;
    } 
    if (dim == 3) {
      for (int i = 0; i < nlocal+nghost; i++)
        if (tris[i].isr >= 0 && sc[tris[i].isc]->allowreact == 0) flag++;
    } 

    if (distributed)
      MPI_Allreduce(&flag,&allflag,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
    else allflag = flag;

    if (allflag) {
      char str[64];
      sprintf(str,BIGINT_FORMAT " surface elements with reaction model, "
              "but invalid collision model",allflag);
      error->all(FLERR,str);
    }
  }    

  // initialize surf collision and reaction models

  for (int i = 0; i < nsc; i++) sc[i]->init();
  for (int i = 0; i < nsr; i++) sr[i]->init();
}

/* ----------------------------------------------------------------------
   remove ghost surfs
------------------------------------------------------------------------- */

void Surf::remove_ghosts()
{
  nghost = 0;
}

/* ----------------------------------------------------------------------
   add a line to owned list
   called by ReadISurf
------------------------------------------------------------------------- */

void Surf::add_line(int itype, double *p1, double *p2)
{
  if (nlocal == nmax) {
    if ((bigint) nmax + DELTA > MAXSMALLINT)
      error->one(FLERR,"Surf add_line overflowed");
    nmax += DELTA;
    grow();
  }
  
  lines[nlocal].type = itype;
  lines[nlocal].mask = 1;
  lines[nlocal].isc = lines[nlocal].isr = -1;
  lines[nlocal].p1[0] = p1[0];
  lines[nlocal].p1[1] = p1[1];
  lines[nlocal].p1[2] = 0.0;
  lines[nlocal].p2[0] = p2[0];
  lines[nlocal].p2[1] = p2[1];
  lines[nlocal].p2[2] = 0.0;
  nlocal++;
}

/* ----------------------------------------------------------------------
   add a line to owned or ghost list, depending on ownflag
   called by Grid::unpack_one
------------------------------------------------------------------------- */

void Surf::add_line_copy(int ownflag, Line *line)
{
  int index;
  
  if (ownflag) {
    if (nlocal == nmax) {
      if ((bigint) nmax + DELTA > MAXSMALLINT)
        error->one(FLERR,"Surf add_line_copy overflowed");
      nmax += DELTA;
      grow();
    }
    index = nlocal;
    nlocal++;
    
  } else {
    if (nlocal+nghost == nmax) {
      if ((bigint) nmax + DELTA > MAXSMALLINT)
        error->one(FLERR,"Surf add_line_copy overflowed");
      nmax += DELTA;
      grow();
    }
    index = nlocal+nghost;
    nghost++;
  }
  
  memcpy(&lines[index],line,sizeof(Line));
}

/* ----------------------------------------------------------------------
   add a line to mylines list
   called by ReadSurf for distributed surfs
------------------------------------------------------------------------- */

void Surf::add_line_own(int itype, double *p1, double *p2)
{
  if (nown == maxown) {
    if ((bigint) maxown + DELTA > MAXSMALLINT)
      error->one(FLERR,"Surf add_line_own overflowed");
    maxown += DELTA;
    grow_own();
  }
  
  mylines[nown].type = itype;
  mylines[nown].mask = 1;
  mylines[nown].isc = mylines[nown].isr = -1;
  mylines[nown].p1[0] = p1[0];
  mylines[nown].p1[1] = p1[1];
  mylines[nown].p1[2] = 0.0;
  mylines[nown].p2[0] = p2[0];
  mylines[nown].p2[1] = p2[1];
  mylines[nown].p2[2] = 0.0;
  nown++;
}

/* ----------------------------------------------------------------------
   add a triangle to owned list
   called by ReadISurf
------------------------------------------------------------------------- */

void Surf::add_tri(int itype, double *p1, double *p2, double *p3)
{
  if (nlocal == nmax) {
    if ((bigint) nmax + DELTA > MAXSMALLINT)
      error->one(FLERR,"Surf add_tri overflowed");
    nmax += DELTA;
    grow();
  }

  tris[nlocal].type = itype;
  tris[nlocal].mask = 1;
  tris[nlocal].isc = tris[nlocal].isr = -1;
  tris[nlocal].p1[0] = p1[0];
  tris[nlocal].p1[1] = p1[1];
  tris[nlocal].p1[2] = p1[2];
  tris[nlocal].p2[0] = p2[0];
  tris[nlocal].p2[1] = p2[1];
  tris[nlocal].p2[2] = p2[2];
  tris[nlocal].p3[0] = p3[0];
  tris[nlocal].p3[1] = p3[1];
  tris[nlocal].p3[2] = p3[2];
  nlocal++;
}

/* ----------------------------------------------------------------------
   add a triangle to owned or ghost list, depending on ownflag
   called by Grid::unpack_one
------------------------------------------------------------------------- */

void Surf::add_tri_copy(int ownflag, Tri *tri)
{
  int index;
  
  if (ownflag) {
    if (nlocal == nmax) {
      if ((bigint) nmax + DELTA > MAXSMALLINT)
        error->one(FLERR,"Surf add_tri_copy overflowed");
      nmax += DELTA;
      grow();
    }
    index = nlocal;
    nlocal++;
    
  } else {
    if (nlocal+nghost == nmax) {
      if ((bigint) nmax + DELTA > MAXSMALLINT)
        error->one(FLERR,"Surf add_tri_copy overflowed");
      nmax += DELTA;
      grow();
    }
    index = nlocal+nghost;
    nghost++;
  }
  
  memcpy(&tris[index],tri,sizeof(Tri));
}

/* ----------------------------------------------------------------------
   add a tri to mytris list
   called by ReadSurf for distributed surfs
------------------------------------------------------------------------- */

void Surf::add_tri_own(int itype, double *p1, double *p2, double *p3)
{
  if (nown == maxown) {
    if ((bigint) maxown + DELTA > MAXSMALLINT)
      error->one(FLERR,"Surf add_tri_own overflowed");
    maxown += DELTA;
    grow_own();
  }
  
  mytris[nown].type = itype;
  mytris[nown].mask = 1;
  mytris[nown].isc = mytris[nown].isr = -1;
  mytris[nown].p1[0] = p1[0];
  mytris[nown].p1[1] = p1[1];
  mytris[nown].p1[2] = p1[2];
  mytris[nown].p2[0] = p2[0];
  mytris[nown].p2[1] = p2[1];
  mytris[nown].p2[2] = p2[2];
  mytris[nown].p3[0] = p3[0];
  mytris[nown].p3[1] = p3[1];
  mytris[nown].p3[2] = p3[2];
  nown++;
}

/* ----------------------------------------------------------------------
   hash all my nlocal surfs with key = ID, value = index
   called only for distributed explicit surfs
------------------------------------------------------------------------- */

void Surf::rehash()
{
  if (implicit) return;

  // hash all nlocal surfs
  // key = ID, value = index into lines or tris

  hash->clear();

  if (domain->dimension == 2) {
    for (int isurf = 0; isurf < nlocal; isurf++)
      (*hash)[lines[isurf].id] = isurf;
  } else {
    for (int isurf = 0; isurf < nlocal; isurf++)
      (*hash)[tris[isurf].id] = isurf;
  }

  hashfilled = 1;
}

/* ----------------------------------------------------------------------
   count owned surf elements = every Pth surf from global Nsurf list
   only called when surfs = explict (all or distributed)
   nothing to do when distributed b/c mylines/mytris already setup
------------------------------------------------------------------------ */

void Surf::setup_owned()
{
  if (distributed) return;

  nown = nsurf/comm->nprocs;
  if (comm->me < nsurf % comm->nprocs) nown++;
}

/* ----------------------------------------------------------------------
   set bounding box around all surfs based on their pts
   for 2d, set zlo,zhi to box bounds
   only called when surfs = explict (all or distributed)
------------------------------------------------------------------------- */

void Surf::setup_bbox()
{
  int i,j;
  double bblo_one[3],bbhi_one[3];

  for (j = 0; j < 3; j++) {
    bblo_one[j] = BIG;
    bbhi_one[j] = -BIG;
  }
  
  int dim = domain->dimension;
  double *x;

  if (dim == 2) {
    for (i = 0; i < nlocal; i++) {
      x = lines[i].p1;
      for (j = 0; j < 2; j++) {
        bblo_one[j] = MIN(bblo_one[j],x[j]);
        bbhi_one[j] = MAX(bbhi_one[j],x[j]);
      }
      x = lines[i].p2;
      for (j = 0; j < 2; j++) {
        bblo_one[j] = MIN(bblo_one[j],x[j]);
        bbhi_one[j] = MAX(bbhi_one[j],x[j]);
      }
    }
    bblo_one[2] = domain->boxlo[2];
    bbhi_one[2] = domain->boxhi[2];

  } else if (dim == 3) {
    for (i = 0; i < nlocal; i++) {
      x = tris[i].p1;
      for (j = 0; j < 3; j++) {
        bblo_one[j] = MIN(bblo_one[j],x[j]);
        bbhi_one[j] = MAX(bbhi_one[j],x[j]);
      }
      x = tris[i].p2;
      for (j = 0; j < 3; j++) {
        bblo_one[j] = MIN(bblo_one[j],x[j]);
        bbhi_one[j] = MAX(bbhi_one[j],x[j]);
      }
      x = tris[i].p3;
      for (j = 0; j < 3; j++) {
        bblo_one[j] = MIN(bblo_one[j],x[j]);
        bbhi_one[j] = MAX(bbhi_one[j],x[j]);
      }
    }
  }

  MPI_Allreduce(bblo_one,bblo,3,MPI_DOUBLE,MPI_MIN,world);
  MPI_Allreduce(bbhi_one,bbhi,3,MPI_DOUBLE,MPI_MAX,world);
}

/* ----------------------------------------------------------------------
   compute unit outward normal vectors of all lines starting at Nold
   outward normal = +z axis x (p2-p1)
------------------------------------------------------------------------- */

void Surf::compute_line_normal(int old)
{
  double z[3],delta[3];

  z[0] = 0.0; z[1] = 0.0; z[2] = 1.0;

  int n;
  Line *newlines;

  if (!implicit && distributed) {
    newlines = mylines;
    n = nown;
  } else {
    newlines = lines;
    n = nlocal;
  }

  for (int i = old; i < n; i++) {
    MathExtra::sub3(newlines[i].p2,newlines[i].p1,delta);
    MathExtra::cross3(z,delta,newlines[i].norm);
    MathExtra::norm3(newlines[i].norm);
    newlines[i].norm[2] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   compute unit outward normal vectors of all lines starting at Nold
   outward normal = (p2-p1) x (p3-p1)
------------------------------------------------------------------------- */

void Surf::compute_tri_normal(int old)
{
  int p1,p2,p3;
  double delta12[3],delta13[3];

  int n;
  Tri *newtris;

  if (!implicit && distributed) {
    newtris = mytris;
    n = nown;
  } else {
    newtris = tris;
    n = nlocal;
  }

  for (int i = old; i < n; i++) {
    MathExtra::sub3(newtris[i].p2,newtris[i].p1,delta12);
    MathExtra::sub3(newtris[i].p3,newtris[i].p1,delta13);
    MathExtra::cross3(delta12,delta13,newtris[i].norm);
    MathExtra::norm3(newtris[i].norm);
  }
}

/* ----------------------------------------------------------------------
   return coords of a corner point in a 2d quad
   icorner pts 1 to 4 are ordered by x, then by y
------------------------------------------------------------------------- */

void Surf::quad_corner_point(int icorner, double *lo, double *hi, double *pt)
{
  if (icorner % 2) pt[0] = hi[0];
  else pt[0] = lo[0];
  if (icorner / 2) pt[1] = hi[1];
  else pt[1] = lo[1];
  pt[2] = 0.0;
}

/* ----------------------------------------------------------------------
   return coords of a corner point in a 3d hex
   icorner pts 1 to 8 are ordered by x, then by y, then by z
------------------------------------------------------------------------- */

void Surf::hex_corner_point(int icorner, double *lo, double *hi, double *pt)
{
  if (icorner % 2) pt[0] = hi[0];
  else pt[0] = lo[0];
  if ((icorner/2) % 2) pt[1] = hi[1];
  else pt[1] = lo[1];
  if (icorner / 4) pt[2] = hi[2];
  else pt[2] = lo[2];
}

/* ----------------------------------------------------------------------
   return length of line M from lines list (not myline)
------------------------------------------------------------------------- */

double Surf::line_size(int m)
{
  return line_size(lines[m].p1,lines[m].p2);
}

/* ----------------------------------------------------------------------
   return length of line
------------------------------------------------------------------------- */

double Surf::line_size(Line *line)
{
  return line_size(line->p1,line->p2);
}

/* ----------------------------------------------------------------------
   return length of line bewteen 2 points
------------------------------------------------------------------------- */

double Surf::line_size(double *p1, double *p2)
{
  double delta[3];
  MathExtra::sub3(p2,p1,delta);
  return MathExtra::len3(delta);
}

/* ----------------------------------------------------------------------
   return area associated with rotating axisymmetric line around y=0 axis
------------------------------------------------------------------------- */

double Surf::axi_line_size(int m)
{
  double *x1 = lines[m].p1;
  double *x2 = lines[m].p2;
  double h = x2[0]-x1[0];
  double r = x2[1]-x1[1];
  double area = MY_PI*(x1[1]+x2[1])*sqrt(r*r+h*h);
  return area;
}

/* ----------------------------------------------------------------------
   return area associated with rotating axisymmetric line around y=0 axis
------------------------------------------------------------------------- */

double Surf::axi_line_size(Line *line)
{
  double *x1 = line->p1;
  double *x2 = line->p2;
  double h = x2[0]-x1[0];
  double r = x2[1]-x1[1];
  double area = MY_PI*(x1[1]+x2[1])*sqrt(r*r+h*h);
  return area;
}

/* ----------------------------------------------------------------------
   compute side length and area of triangle M from tri list (not mytri)
   return len = length of shortest edge of triangle M
   return area = area of triangle M
------------------------------------------------------------------------- */

double Surf::tri_size(int m, double &len)
{
  return tri_size(tris[m].p1,tris[m].p2,tris[m].p3,len);
}

/* ----------------------------------------------------------------------
   compute side length and area of triangle tri
   return len = length of shortest edge of triangle M
   return area = area of triangle M
------------------------------------------------------------------------- */

double Surf::tri_size(Tri *tri, double &len)
{
  return tri_size(tri->p1,tri->p2,tri->p3,len);
}

/* ----------------------------------------------------------------------
   compute side length and area of a triangle
   return len = length of shortest edge of triangle M
   return area = area of triangle M
------------------------------------------------------------------------- */

double Surf::tri_size(double *p1, double *p2, double *p3, double &len)
{
  double delta12[3],delta13[3],delta23[3],cross[3];

  MathExtra::sub3(p2,p1,delta12);
  MathExtra::sub3(p3,p1,delta13);
  MathExtra::sub3(p3,p2,delta23);
  len = MIN(MathExtra::len3(delta12),MathExtra::len3(delta13));
  len = MIN(len,MathExtra::len3(delta23));

  MathExtra::cross3(delta12,delta13,cross);
  double area = 0.5 * MathExtra::len3(cross);
  return area;
}

/* ----------------------------------------------------------------------
   check end points of lines
   each end point should appear exactly once as different ends of 2 lines
   exception: not required of end point on simulation box surface
   only check lines newer than old ones
------------------------------------------------------------------------- */

void Surf::check_watertight_2d(int nline_old)
{
  // NOTE: need to enable this for explicit and implicit distributed

  if (distributed) return;
  
  // hash end points of all lines
  // key = end point
  // value = 1 if first point, 2 if second point, 3 if both points
  // NOTE: could prealloc hash to correct size here

  MyHashPoint phash;
  MyPointIt it;

  // insert each end point into hash
  // should appear once at each end
  // error if any duplicate points

  double *p1,*p2;
  OnePoint2d key;
  int value;

  int nline_new = nsurf - nline_old;

  int ndup = 0;
  int m = nline_old;
  for (int i = 0; i < nline_new; i++) {
    p1 = lines[m].p1;
    key.pt[0] = p1[0]; key.pt[1] = p1[1];
    if (phash.find(key) == phash.end()) phash[key] = 1;
    else {
      value = phash[key];
      if (value == 2) phash[key] = 3;
      else ndup++;
    }

    p2 = lines[m].p2;
    key.pt[0] = p2[0]; key.pt[1] = p2[1];
    if (phash.find(key) == phash.end()) phash[key] = 2;
    else {
      value = phash[key];
      if (value == 1) phash[key] = 3;
      else ndup++;
    }

    m++;
  }
  
  if (ndup) {
    char str[128];
    sprintf(str,"Watertight check failed with %d duplicate points",ndup);
    error->all(FLERR,str);
  }

  // check that each end point has a match
  // allow for exception if end point on box surface

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  double *kpt;
  double pt[3];

  int nbad = 0;
  for (it = phash.begin(); it != phash.end(); ++it) {
    if (it->second != 3) {
      kpt = (double *) it->first.pt;
      pt[0] = kpt[0]; pt[1] = kpt[1]; pt[2] = 0.0;
      if (!Geometry::point_on_hex(pt,boxlo,boxhi)) nbad++;
    }
  }

  if (nbad) {
    char str[128];
    sprintf(str,"Watertight check failed with %d unmatched points",nbad);
    error->all(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   check directed triangle edges
   each edge should appear exactly once in each direction
   exception: not required of triangle edge on simulation box surface
   only check triangles newer than old ones
------------------------------------------------------------------------- */

void Surf::check_watertight_3d(int ntri_old)
{
  // NOTE: need to enable this for explicit and implicit
  if (distributed) return;

  // hash directed edges of all triangles
  // key = directed edge
  // value = 1 if appears once, 2 if reverse also appears once
  // NOTE: could prealloc hash to correct size here

  MyHash2Point phash;
  My2PointIt it;

  // insert each edge into hash
  // should appear once in each direction
  // error if any duplicate edges
  
  double *p1,*p2,*p3;
  TwoPoint3d key,keyinv;
  int value;

  int ntri_new = nsurf - ntri_old;

  int ndup = 0;
  int m = ntri_old;
  for (int i = 0; i < ntri_new; i++) {
    p1 = tris[m].p1;
    p2 = tris[m].p2;
    p3 = tris[m].p3;

    key.pts[0] = p1[0]; key.pts[1] = p1[1]; key.pts[2] = p1[2];
    key.pts[3] = p2[0]; key.pts[4] = p2[1]; key.pts[5] = p2[2];
    if (phash.find(key) == phash.end()) {
      keyinv.pts[0] = p2[0]; keyinv.pts[1] = p2[1]; keyinv.pts[2] = p2[2];
      keyinv.pts[3] = p1[0]; keyinv.pts[4] = p1[1]; keyinv.pts[5] = p1[2];
      if (phash.find(keyinv) == phash.end()) phash[key] = 1;
      else {
	value = phash[keyinv];
	if (value == 1) phash[keyinv] = 2;
	else ndup++;
      }
    } else ndup++;
    
    key.pts[0] = p2[0]; key.pts[1] = p2[1]; key.pts[2] = p2[2];
    key.pts[3] = p3[0]; key.pts[4] = p3[1]; key.pts[5] = p3[2];
    if (phash.find(key) == phash.end()) {
      keyinv.pts[0] = p3[0]; keyinv.pts[1] = p3[1]; keyinv.pts[2] = p3[2];
      keyinv.pts[3] = p2[0]; keyinv.pts[4] = p2[1]; keyinv.pts[5] = p2[2];
      if (phash.find(keyinv) == phash.end()) phash[key] = 1;
      else {
	value = phash[keyinv];
	if (value == 1) phash[keyinv] = 2;
	else ndup++;
      }
    } else ndup++;

    key.pts[0] = p3[0]; key.pts[1] = p3[1]; key.pts[2] = p3[2];
    key.pts[3] = p1[0]; key.pts[4] = p1[1]; key.pts[5] = p1[2];
    if (phash.find(key) == phash.end()) {
      keyinv.pts[0] = p1[0]; keyinv.pts[1] = p1[1]; keyinv.pts[2] = p1[2];
      keyinv.pts[3] = p3[0]; keyinv.pts[4] = p3[1]; keyinv.pts[5] = p3[2];
      if (phash.find(keyinv) == phash.end()) phash[key] = 1;
      else {
	value = phash[keyinv];
	if (value == 1) phash[keyinv] = 2;
	else ndup++;
      } 
    } else ndup++;
      
    m++;
  }

  if (ndup) {
    char str[128];
    sprintf(str,"Watertight check failed with %d duplicate edges",ndup);
    error->all(FLERR,str);
  }

  // check that each edge has an inverted match
  // allow for exception if edge on box surface
  // NOTE: this does not check if 2 pts of edge are on same box face

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  double *pts;

  int nbad = 0;
  for (it = phash.begin(); it != phash.end(); ++it) {
    if (it->second != 2) {
      pts = (double *) it->first.pts;
      if (!Geometry::point_on_hex(&pts[0],boxlo,boxhi)) nbad++;
      if (!Geometry::point_on_hex(&pts[3],boxlo,boxhi)) nbad++;
    }
  }

  if (nbad) {
    char str[128];
    sprintf(str,"Watertight check failed with %d unmatched edges",nbad);
    error->all(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   check if all points are inside or on surface of global simulation box
   called by ReadSurf for lines or triangles
   old = previous # of elements
------------------------------------------------------------------------- */

void Surf::check_point_inside(int old)
{
  int nbad;
  double *x;

  int dim = domain->dimension;
  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  if (dim == 2) {
    Line *newlines;
    int n;
    if (distributed) {
      newlines = mylines;
      n = nown;
    } else {
      newlines = lines;
      n = nlocal;
    }

    nbad = 0;
    for (int i = old; i < n; i++) {
      x = newlines[i].p1;
      if (x[0] < boxlo[0] || x[0] > boxhi[0] ||
	  x[1] < boxlo[1] || x[1] > boxhi[1] ||
	  x[2] < boxlo[2] || x[2] > boxhi[2]) nbad++;
      x = newlines[i].p2;
      if (x[0] < boxlo[0] || x[0] > boxhi[0] ||
	  x[1] < boxlo[1] || x[1] > boxhi[1] ||
	  x[2] < boxlo[2] || x[2] > boxhi[2]) nbad++;
    }

  } else if (dim == 3) {
    Tri *newtris;
    int n;
    if (distributed) {
      newtris = mytris;
      n = nown;
    } else {
      newtris = tris;
      n = nlocal;
    }

    nbad = 0;
    for (int i = old; i < n; i++) {
      x = newtris[i].p1;
      if (x[0] < boxlo[0] || x[0] > boxhi[0] ||
	  x[1] < boxlo[1] || x[1] > boxhi[1] ||
	  x[2] < boxlo[2] || x[2] > boxhi[2]) nbad++;
      x = newtris[i].p2;
      if (x[0] < boxlo[0] || x[0] > boxhi[0] ||
	  x[1] < boxlo[1] || x[1] > boxhi[1] ||
	  x[2] < boxlo[2] || x[2] > boxhi[2]) nbad++;
      x = newtris[i].p3;
      if (x[0] < boxlo[0] || x[0] > boxhi[0] ||
	  x[1] < boxlo[1] || x[1] > boxhi[1] ||
	  x[2] < boxlo[2] || x[2] > boxhi[2]) nbad++;
    }
  }

  int nbadall;
  if (distributed) MPI_Allreduce(&nbad,&nbadall,1,MPI_INT,MPI_SUM,world);
  else nbadall = nbad;

  if (nbadall) {
    char str[128];
    sprintf(str,"%d surface points are not inside simulation box",
	    nbadall);
    error->all(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   check nearness of all points to other lines in same cell
   error if point is on line, including duplicate point
   warn if closer than EPSILON_GRID = fraction of grid cell size
   NOTE: this can miss a close point/line pair in 2 different grid cells
------------------------------------------------------------------------- */

void Surf::check_point_near_surf_2d()
{
  int i,j,n;
  surfint *csurfs;
  double side,epssq;
  double *p1,*p2,*lo,*hi;
  Surf::Line *line;

  Surf::Line *lines = surf->lines;
  Grid::ChildCell *cells = grid->cells;
  int nglocal = grid->nlocal;

  int nerror = 0;
  int nwarn = 0;

  for (int icell = 0; icell < nglocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    n = cells[icell].nsurf;
    if (n == 0) continue;

    lo = cells[icell].lo;
    hi = cells[icell].hi;
    side = MIN(hi[0]-lo[0],hi[1]-lo[1]);
    epssq = (EPSILON_GRID*side) * (EPSILON_GRID*side);

    csurfs = cells[icell].csurfs;
    for (i = 0; i < n; i++) {
      line = &lines[csurfs[i]];
      for (j = 0; j < n; j++) {
        if (i == j) continue;
        p1 = lines[csurfs[j]].p1;
        p2 = lines[csurfs[j]].p2;
        point_line_compare(p1,line->p1,line->p2,epssq,nerror,nwarn);
        point_line_compare(p2,line->p1,line->p2,epssq,nerror,nwarn);
      }
    }
  }

  int all;
  MPI_Allreduce(&nerror,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,"Surface check failed with %d points on lines",all);
    error->all(FLERR,str);
  }

  MPI_Allreduce(&nwarn,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,"Surface check found %d points nearly on lines",all);
    if (comm->me == 0) error->warning(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   check nearness of all points to other triangles in same cell
   error if point is on triangle, including duplicate point
   warn if closer than EPSILON_GRID = fraction of grid cell size
   NOTE: this can miss a close point/triangle pair in 2 different grid cells
------------------------------------------------------------------------- */

void Surf::check_point_near_surf_3d()
{
  int i,j,n;
  surfint *csurfs;
  double side,epssq;
  double *p1,*p2,*p3,*lo,*hi;
  Surf::Tri *tri;

  Surf::Tri *tris = surf->tris;
  Grid::ChildCell *cells = grid->cells;
  int nglocal = grid->nlocal;

  int nerror = 0;
  int nwarn = 0;

  for (int icell = 0; icell < nglocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    n = cells[icell].nsurf;
    if (n == 0) continue;

    lo = cells[icell].lo;
    hi = cells[icell].hi;
    side = MIN(hi[0]-lo[0],hi[1]-lo[1]);
    side = MIN(side,hi[2]-lo[2]);
    epssq = (EPSILON_GRID*side) * (EPSILON_GRID*side);

    csurfs = cells[icell].csurfs;
    for (i = 0; i < n; i++) {
      tri = &tris[csurfs[i]];
      for (j = 0; j < n; j++) {
        if (i == j) continue;
        p1 = tris[csurfs[j]].p1;
        p2 = tris[csurfs[j]].p2;
        p3 = tris[csurfs[j]].p3;
        point_tri_compare(p1,tri->p1,tri->p2,tri->p3,tri->norm,
                          epssq,nerror,nwarn,icell,csurfs[i],csurfs[j]);
        point_tri_compare(p2,tri->p1,tri->p2,tri->p3,tri->norm,
                          epssq,nerror,nwarn,icell,csurfs[i],csurfs[j]);
        point_tri_compare(p3,tri->p1,tri->p2,tri->p3,tri->norm,
                          epssq,nerror,nwarn,icell,csurfs[i],csurfs[j]);
      }
    }
  }

  int all;
  MPI_Allreduce(&nerror,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,"Surface check failed with %d points on triangles",all);
    error->all(FLERR,str);
  }

  MPI_Allreduce(&nwarn,&all,1,MPI_INT,MPI_SUM,world);
  if (all) {
    char str[128];
    sprintf(str,"Surface check found %d points nearly on triangles",all);
    if (comm->me == 0) error->warning(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   compute extent of read-in surfs, including geometric transformations
------------------------------------------------------------------------- */

void Surf::output_extent(int old)
{
  // extent of surfs after geometric transformations
  // compute sizes of smallest surface elements

  double extent[3][2],extentall[3][2];
  extent[0][0] = extent[1][0] = extent[2][0] = BIG;
  extent[0][1] = extent[1][1] = extent[2][1] = -BIG;

  int dim = domain->dimension;

  if (dim == 2) {
    Line *newlines;
    int n;
    if (!implicit && distributed) {
      newlines = mylines;
      n = nown;
    } else {
      newlines = lines;
      n = nlocal;
    }

    for (int i = old; i < n; i++) {
      for (int j = 0; j < 3; j++) {
        extent[j][0] = MIN(extent[j][0],newlines[i].p1[j]);
        extent[j][0] = MIN(extent[j][0],newlines[i].p2[j]);
        extent[j][1] = MAX(extent[j][1],newlines[i].p1[j]);
        extent[j][1] = MAX(extent[j][1],newlines[i].p2[j]);
      }
    }

  } else {
    Tri *newtris;
    int n;
    if (!implicit && distributed) {
      newtris = mytris;
      n = nown;
    } else {
      newtris = tris;
      n = nlocal;
    }

    for (int i = old; i < n; i++) {
      for (int j = 0; j < 3; j++) {
        extent[j][0] = MIN(extent[j][0],newtris[i].p1[j]);
        extent[j][0] = MIN(extent[j][0],newtris[i].p2[j]);
        extent[j][0] = MIN(extent[j][0],newtris[i].p3[j]);
        extent[j][1] = MAX(extent[j][1],newtris[i].p1[j]);
        extent[j][1] = MAX(extent[j][1],newtris[i].p2[j]);
        extent[j][1] = MAX(extent[j][1],newtris[i].p3[j]);
      }
    }
  }

  extent[0][0] = -extent[0][0];
  extent[1][0] = -extent[1][0];
  extent[2][0] = -extent[2][0];
  MPI_Allreduce(extent,extentall,6,MPI_DOUBLE,MPI_MAX,world);
  extentall[0][0] = -extentall[0][0];
  extentall[1][0] = -extentall[1][0];
  extentall[2][0] = -extentall[2][0];

  double minlen,minarea; 
  if (dim == 2) minlen = shortest_line(old);
  if (dim == 3) smallest_tri(old,minlen,minarea);

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"  %g %g xlo xhi\n",extentall[0][0],extentall[0][1]);
      fprintf(screen,"  %g %g ylo yhi\n",extentall[1][0],extentall[1][1]);
      fprintf(screen,"  %g %g zlo zhi\n",extentall[2][0],extentall[2][1]);
      if (dim == 2)
	fprintf(screen,"  %g min line length\n",minlen);
      if (dim == 3) {
	fprintf(screen,"  %g min triangle edge length\n",minlen);
	fprintf(screen,"  %g min triangle area\n",minarea);
      }
    }
    if (logfile) {
      fprintf(logfile,"  %g %g xlo xhi\n",extentall[0][0],extentall[0][1]);
      fprintf(logfile,"  %g %g ylo yhi\n",extentall[1][0],extentall[1][1]);
      fprintf(logfile,"  %g %g zlo zhi\n",extentall[2][0],extentall[2][1]);
      if (dim == 2)
	fprintf(logfile,"  %g min line length\n",minlen);
      if (dim == 3) {
	fprintf(logfile,"  %g min triangle edge length\n",minlen);
	fprintf(logfile,"  %g min triangle area\n",minarea);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   return shortest line length
------------------------------------------------------------------------- */

double Surf::shortest_line(int old)
{
  double len = BIG;

  if (!implicit && distributed) {
    for (int i = old; i < nown; i++)
      len = MIN(len,line_size(&mylines[i]));
  } else {
    for (int i = old; i < nlocal; i++)
      len = MIN(len,line_size(&lines[i]));
  }

  double lenall;
  MPI_Allreduce(&len,&lenall,1,MPI_DOUBLE,MPI_MIN,world);

  return lenall;
}

/* ----------------------------------------------------------------------
   return shortest tri edge and smallest tri area
------------------------------------------------------------------------- */

void Surf::smallest_tri(int old, double &lenall, double &areaall)
{
  double lenone,areaone;
  double len = BIG;
  double area = BIG;

  if (!implicit && distributed) {
    for (int i = old; i < nown; i++) {
      areaone = tri_size(&mytris[i],lenone);
      len = MIN(len,lenone);
      area = MIN(area,areaone);
    }
  } else {
    for (int i = old; i < nlocal; i++) {
      areaone = tri_size(&tris[i],lenone);
      len = MIN(len,lenone);
      area = MIN(area,areaone);
    }
  }

  MPI_Allreduce(&len,&lenall,1,MPI_DOUBLE,MPI_MIN,world);
  MPI_Allreduce(&area,&areaall,1,MPI_DOUBLE,MPI_MIN,world);
}

/* ----------------------------------------------------------------------
   compute distance bewteen a point and line
   just return if point is an endpoint of line
   increment nerror if point on line
   increment nwarn if point is within epssq distance of line
------------------------------------------------------------------------- */

void Surf::point_line_compare(double *pt, double *p1, double *p2,
                              double epssq, int &nerror, int &nwarn)
{
  if (pt[0] == p1[0] && pt[1] == p1[1]) return;
  if (pt[0] == p2[0] && pt[1] == p2[1]) return;
  double rsq = Geometry::distsq_point_line(pt,p1,p2);
  if (rsq == 0.0) nerror++;
  else if (rsq < epssq) nwarn++;
}

/* ----------------------------------------------------------------------
   compute distance bewteen a point and triangle
   just return if point is an endpoint of triangle
   increment nerror if point on triangle
   increment nwarn if point is within epssq distance of triangle
------------------------------------------------------------------------- */

void Surf::point_tri_compare(double *pt, double *p1, double *p2, double *p3,
                             double *norm, double epssq, int &nerror, int &nwarn,
                             int, int, int) 
{
  if (pt[0] == p1[0] && pt[1] == p1[1] && pt[2] == p1[2]) return;
  if (pt[0] == p2[0] && pt[1] == p2[1] && pt[2] == p2[2]) return;
  if (pt[0] == p3[0] && pt[1] == p3[1] && pt[2] == p3[2]) return;
  double rsq = Geometry::distsq_point_tri(pt,p1,p2,p3,norm);
  if (rsq == 0.0) nerror++;
  else if (rsq < epssq) nwarn++;
}


/* ----------------------------------------------------------------------
   add a surface collision model
------------------------------------------------------------------------- */

void Surf::add_collide(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal surf_collide command");

  // error check

  for (int i = 0; i < nsc; i++)
    if (strcmp(arg[0],sc[i]->id) == 0)
      error->all(FLERR,"Reuse of surf_collide ID");

  // extend SurfCollide list if necessary

  if (nsc == maxsc) {
    maxsc += DELTAMODEL;
    sc = (SurfCollide **)
      memory->srealloc(sc,maxsc*sizeof(SurfCollide *),"surf:sc");
  }

  // create new SurfCollide class

  if (sparta->suffix_enable) {
    if (sparta->suffix) {
      char estyle[256];
      sprintf(estyle,"%s/%s",arg[1],sparta->suffix);

      if (0) return;

#define SURF_COLLIDE_CLASS
#define SurfCollideStyle(key,Class) \
      else if (strcmp(estyle,#key) == 0) { \
        sc[nsc] = new Class(sparta,narg,arg); \
        nsc++; \
        return; \
      }
#include "style_surf_collide.h"
#undef SurfCollideStyle
#undef SURF_COLLIDE_CLASS
    }
  }

  if (0) return;

#define SURF_COLLIDE_CLASS
#define SurfCollideStyle(key,Class) \
  else if (strcmp(arg[1],#key) == 0) \
    sc[nsc] = new Class(sparta,narg,arg);
#include "style_surf_collide.h"
#undef SurfCollideStyle
#undef SURF_COLLIDE_CLASS

  else error->all(FLERR,"Unrecognized surf_collide style");

  nsc++;
}

/* ----------------------------------------------------------------------
   find a surface collide model by ID
   return index of surf collide model or -1 if not found
------------------------------------------------------------------------- */

int Surf::find_collide(const char *id)
{
  int isc;
  for (isc = 0; isc < nsc; isc++)
    if (strcmp(id,sc[isc]->id) == 0) break;
  if (isc == nsc) return -1;
  return isc;
}

/* ----------------------------------------------------------------------
   add a surface reaction model
------------------------------------------------------------------------- */

void Surf::add_react(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal surf_react command");

  // error check

  for (int i = 0; i < nsr; i++)
    if (strcmp(arg[0],sr[i]->id) == 0)
      error->all(FLERR,"Reuse of surf_react ID");

  // extend SurfReact list if necessary

  if (nsr == maxsr) {
    maxsr += DELTAMODEL;
    sr = (SurfReact **)
      memory->srealloc(sr,maxsr*sizeof(SurfReact *),"surf:sr");
  }

  // create new SurfReact class

  if (0) return;

#define SURF_REACT_CLASS
#define SurfReactStyle(key,Class) \
  else if (strcmp(arg[1],#key) == 0) \
    sr[nsr] = new Class(sparta,narg,arg);
#include "style_surf_react.h"
#undef SurfReactStyle
#undef SURF_REACT_CLASS

  else error->all(FLERR,"Unrecognized surf_react style");

  nsr++;
}

/* ----------------------------------------------------------------------
   find a surface reaction model by ID
   return index of surf reaction model or -1 if not found
------------------------------------------------------------------------- */

int Surf::find_react(const char *id)
{
  int isr;
  for (isr = 0; isr < nsr; isr++)
    if (strcmp(id,sr[isr]->id) == 0) break;
  if (isr == nsr) return -1;
  return isr;
}

/* ----------------------------------------------------------------------
   group surf command called via input script
   NOTE: need to apply this also to mylines and mytris ??
------------------------------------------------------------------------- */

void Surf::group(int narg, char **arg)
{
  int i,flag;
  double x[3];

  if (narg < 3) error->all(FLERR,"Illegal group command");

  int dim = domain->dimension;

  int igroup = find_group(arg[0]);
  if (igroup < 0) igroup = add_group(arg[0]);
  int bit = bitmask[igroup];

  // style = type or id
  // add surf to group if matches types/ids or condition

  if (strcmp(arg[2],"type") == 0 || strcmp(arg[2],"id") == 0) {
    if (narg < 4) error->all(FLERR,"Illegal group command");

    int category;
    if (strcmp(arg[2],"type") == 0) category = TYPE;
    else if (strcmp(arg[2],"id") == 0) category = ID;

    // args = logical condition
    
    if (narg > 4 &&
        (strcmp(arg[3],"<") == 0 || strcmp(arg[3],">") == 0 ||
         strcmp(arg[3],"<=") == 0 || strcmp(arg[3],">=") == 0 ||
         strcmp(arg[3],"==") == 0 || strcmp(arg[3],"!=") == 0 ||
         strcmp(arg[3],"<>") == 0)) {

      int condition = -1;
      if (strcmp(arg[3],"<") == 0) condition = LT;
      else if (strcmp(arg[3],"<=") == 0) condition = LE;
      else if (strcmp(arg[3],">") == 0) condition = GT;
      else if (strcmp(arg[3],">=") == 0) condition = GE;
      else if (strcmp(arg[3],"==") == 0) condition = EQ;
      else if (strcmp(arg[3],"!=") == 0) condition = NEQ;
      else if (strcmp(arg[3],"<>") == 0) condition = BETWEEN;
      else error->all(FLERR,"Illegal group command");

      int bound1,bound2;
      bound1 = input->inumeric(FLERR,arg[4]);
      bound2 = -1;

      if (condition == BETWEEN) {
        if (narg != 6) error->all(FLERR,"Illegal group command");
        bound2 = input->inumeric(FLERR,arg[5]);
      } else if (narg != 5) error->all(FLERR,"Illegal group command");

      // add surf to group if meets condition

      if (category == ID) {
        if (condition == LT) {
          if (dim == 2) {
            for (i = 0; i < nlocal+nghost; i++) 
              if (lines[i].id < bound1) lines[i].mask |= bit;
          } else {
            for (i = 0; i < nlocal+nghost; i++) 
              if (tris[i].id < bound1) tris[i].mask |= bit;
          }
        } else if (condition == LE) {
          if (dim == 2) {
            for (i = 0; i < nlocal+nghost; i++) 
              if (lines[i].id <= bound1) lines[i].mask |= bit;
          } else {
            for (i = 0; i < nlocal+nghost; i++) 
              if (tris[i].id <= bound1) tris[i].mask |= bit;
          }
        } else if (condition == GT) {
          if (dim == 2) {
            for (i = 0; i < nlocal+nghost; i++) 
              if (lines[i].id > bound1) lines[i].mask |= bit;
          } else {
            for (i = 0; i < nlocal+nghost; i++) 
              if (tris[i].id > bound1) tris[i].mask |= bit;
          }
        } else if (condition == GE) {
          if (dim == 2) {
            for (i = 0; i < nlocal+nghost; i++) 
              if (lines[i].id >= bound1) lines[i].mask |= bit;
          } else {
            for (i = 0; i < nlocal+nghost; i++) 
              if (tris[i].id >= bound1) tris[i].mask |= bit;
          }
        } else if (condition == EQ) {
          if (dim == 2) {
            for (i = 0; i < nlocal+nghost; i++) 
              if (lines[i].id == bound1) lines[i].mask |= bit;
          } else {
            for (i = 0; i < nlocal+nghost; i++) 
              if (tris[i].id == bound1) tris[i].mask |= bit;
          }
        } else if (condition == NEQ) {
          if (dim == 2) {
            for (i = 0; i < nlocal+nghost; i++) 
              if (lines[i].id != bound1) lines[i].mask |= bit;
          } else {
            for (i = 0; i < nlocal+nghost; i++) 
              if (tris[i].id != bound1) tris[i].mask |= bit;
          }
        } else if (condition == BETWEEN) {
          if (dim == 2) {
            for (i = 0; i < nlocal+nghost; i++)
              if (lines[i].id >= bound1 && lines[i].id <= bound2) 
                lines[i].mask |= bit;
          } else {
            for (i = 0; i < nlocal+nghost; i++)
              if (tris[i].id >= bound1 && tris[i].id <= bound2) 
                tris[i].mask |= bit;
          }
        }
      } else if (category == TYPE) {
        if (condition == LT) {
          if (dim == 2) {
            for (i = 0; i < nlocal+nghost; i++) 
              if (lines[i].type < bound1) lines[i].mask |= bit;
          } else {
            for (i = 0; i < nlocal+nghost; i++) 
              if (tris[i].type < bound1) lines[i].mask |= bit;
          }
        } else if (condition == LE) {
          if (dim == 2) {
            for (i = 0; i < nlocal+nghost; i++) 
              if (lines[i].type <= bound1) lines[i].mask |= bit;
          } else {
            for (i = 0; i < nlocal+nghost; i++) 
              if (tris[i].type <= bound1) lines[i].mask |= bit;
          }
        } else if (condition == GT) {
          if (dim == 2) {
            for (i = 0; i < nlocal+nghost; i++) 
              if (lines[i].type > bound1) lines[i].mask |= bit;
          } else {
            for (i = 0; i < nlocal+nghost; i++) 
              if (tris[i].type > bound1) lines[i].mask |= bit;
          }
        } else if (condition == GE) {
          if (dim == 2) {
            for (i = 0; i < nlocal+nghost; i++) 
              if (lines[i].type >= bound1) lines[i].mask |= bit;
          } else {
            for (i = 0; i < nlocal+nghost; i++) 
              if (tris[i].type >= bound1) lines[i].mask |= bit;
          }
        } else if (condition == EQ) {
          if (dim == 2) {
            for (i = 0; i < nlocal+nghost; i++) 
              if (lines[i].type == bound1) lines[i].mask |= bit;
          } else {
            for (i = 0; i < nlocal+nghost; i++) 
              if (tris[i].type == bound1) lines[i].mask |= bit;
          }
        } else if (condition == NEQ) {
          if (dim == 2) {
            for (i = 0; i < nlocal+nghost; i++) 
              if (lines[i].type != bound1) lines[i].mask |= bit;
          } else {
            for (i = 0; i < nlocal+nghost; i++) 
              if (tris[i].type != bound1) lines[i].mask |= bit;
          }
        } else if (condition == BETWEEN) {
          if (dim == 2) {
            for (i = 0; i < nlocal+nghost; i++)
              if (lines[i].type >= bound1 && lines[i].type <= bound2) 
                lines[i].mask |= bit;
          } else {
            for (i = 0; i < nlocal+nghost; i++)
              if (tris[i].type >= bound1 && tris[i].type <= bound2) 
                tris[i].mask |= bit;
          }
        }
      }

    // args = list of values
      
    } else {
      char *ptr;
      int start,stop;

      for (int iarg = 3; iarg < narg; iarg++) {
        if (strchr(arg[iarg],':')) {
          ptr = strchr(arg[iarg],':'); 
          *ptr = '\0';
          start = input->inumeric(FLERR,arg[iarg]); 
          *ptr = ':';
          stop = input->inumeric(FLERR,ptr+1); 
        } else {
          start = stop = input->inumeric(FLERR,arg[iarg]);
        }

        // add surf to group if type/id matches value or sequence
      
        if (category == ID) {
          if (dim == 2) {
            for (i = 0; i < nlocal+nghost; i++)
              if (lines[i].id >= start && lines[i].id <= stop) 
                lines[i].mask |= bit;
          } else {
            for (i = 0; i < nlocal+nghost; i++)
              if (tris[i].id >= start && tris[i].id <= stop) 
                tris[i].mask |= bit;
          }
        } else if (category == TYPE) {
          if (dim == 2) {
            for (i = 0; i < nlocal+nghost; i++)
              if (lines[i].type >= start && lines[i].type <= stop) 
                lines[i].mask |= bit;
          } else {
            for (i = 0; i < nlocal+nghost; i++)
              if (tris[i].type >= start && tris[i].type <= stop) 
                tris[i].mask |= bit;
          }
        }
      }
    }

  // style = region
  // add surf to group if in region

  } else if (strcmp(arg[2],"region") == 0) {
    if (narg != 5) error->all(FLERR,"Illegal group command");
    int iregion = domain->find_region(arg[3]);
    if (iregion == -1) error->all(FLERR,"Group region ID does not exist");
    Region *region = domain->regions[iregion];

    int rstyle;
    if (strcmp(arg[4],"all") == 0) rstyle = REGION_ALL;
    else if (strcmp(arg[4],"one") == 0) rstyle = REGION_ONE;
    else if (strcmp(arg[4],"center") == 0) rstyle = REGION_CENTER;
    else error->all(FLERR,"Illegal group command");

    if (dim == 2) {
      if (rstyle == REGION_ALL) {
        for (i = 0; i < nlocal+nghost; i++) {
          flag = 1;
          if (!region->match(lines[i].p1)) flag = 0;
          if (!region->match(lines[i].p2)) flag = 0;
          if (flag) lines[i].mask |= bit;
        }
      } else if (rstyle == REGION_ONE) {
        for (i = 0; i < nlocal+nghost; i++) {
          flag = 0;
          if (region->match(lines[i].p1)) flag = 1;
          if (region->match(lines[i].p2)) flag = 1;
          if (flag) lines[i].mask |= bit;
        }
      } else if (rstyle == REGION_CENTER) {
        for (i = 0; i < nlocal+nghost; i++) {
          x[0] = 0.5 * (lines[i].p1[0] + lines[i].p2[0]);
          x[1] = 0.5 * (lines[i].p1[1] + lines[i].p2[1]);
          x[2] = 0.0;
          if (region->match(x)) lines[i].mask |= bit;
        }
      }

    } else if (dim == 3) {
      if (rstyle == REGION_ALL) {
        for (i = 0; i < nlocal+nghost; i++) {
          flag = 1;
          if (!region->match(tris[i].p1)) flag = 0;
          if (!region->match(tris[i].p2)) flag = 0;
          if (!region->match(tris[i].p3)) flag = 0;
          if (flag) tris[i].mask |= bit;
        }
      } else if (rstyle == REGION_ONE) {
        for (i = 0; i < nlocal+nghost; i++) {
          flag = 0;
          if (region->match(tris[i].p1)) flag = 1;
          if (region->match(tris[i].p2)) flag = 1;
          if (region->match(tris[i].p3)) flag = 1;
          if (flag) tris[i].mask |= bit;
        }
      } else if (rstyle == REGION_CENTER) {
        for (i = 0; i < nlocal+nghost; i++) {
          x[0] = (tris[i].p1[0] + tris[i].p2[0] + tris[i].p3[0]) / 3.0;
          x[1] = (tris[i].p1[1] + tris[i].p2[1] + tris[i].p3[1]) / 3.0;
          x[2] = (tris[i].p1[2] + tris[i].p2[2] + tris[i].p3[2]) / 3.0;
          if (region->match(x)) tris[i].mask |= bit;
        }
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

    if (dim == 2) {
      for (i = 0; i < nlocal+nghost; i++)
        if (lines[i].mask & otherbit) lines[i].mask |= bit;
    } else {
      for (i = 0; i < nlocal+nghost; i++)
        if (tris[i].mask & otherbit) tris[i].mask |= bit;
    }

    // remove surfs if they are in any of the other groups
    // AND with inverse mask removes the surf from group

    int inverse = inversemask[igroup];

    for (int ilist = 1; ilist < length; ilist++) {
      otherbit = bitmask[list[ilist]];
      if (dim == 2) {
        for (i = 0; i < nlocal+nghost; i++)
          if (lines[i].mask & otherbit) lines[i].mask &= inverse;
      } else {
        for (i = 0; i < nlocal+nghost; i++)
          if (tris[i].mask & otherbit) tris[i].mask &= inverse;
      }
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
      if (dim == 2) {
        for (i = 0; i < nlocal+nghost; i++)
          if (lines[i].mask & otherbit) lines[i].mask |= bit;
      } else {
        for (i = 0; i < nlocal+nghost; i++)
          if (tris[i].mask & otherbit) tris[i].mask |= bit;
      }
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

    if (dim == 2) {
      for (i = 0; i < nlocal+nghost; i++) {
        ok = 1;
        for (ilist = 0; ilist < length; ilist++) {
          otherbit = bitmask[list[ilist]];
          if ((lines[i].mask & otherbit) == 0) ok = 0;
        }
        if (ok) lines[i].mask |= bit;
      }
    } else {
      for (i = 0; i < nlocal+nghost; i++) {
        ok = 1;
        for (ilist = 0; ilist < length; ilist++) {
          otherbit = bitmask[list[ilist]];
          if ((tris[i].mask & otherbit) == 0) ok = 0;
        }
        if (ok) tris[i].mask |= bit;
      }
    }

    delete [] list;

  // style = clear
  // remove all surfs from group

  } else if (strcmp(arg[2],"clear") == 0) {
    if (igroup == 0) error->all(FLERR,"Cannot clear group all");
    int inversebits = inversemask[igroup];

    if (dim == 2) {
      for (i = 0; i < nlocal+nghost; i++) lines[i].mask &= inversebits;
    } else {
      for (i = 0; i < nlocal+nghost; i++) tris[i].mask &= inversebits;
    }
  }

  // print stats for changed group

  bigint n = 0;
  if (dim == 2) {
    for (i = 0; i < nlocal; i++) 
      if (lines[i].mask & bit) n++;
  } else {
    for (i = 0; i < nlocal; i++)
      if (tris[i].mask & bit) n++;
  }

  bigint nall;
  if (distributed) MPI_Allreduce(&n,&nall,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  else nall = n;

  if (comm->me == 0) {
    if (screen) 
      fprintf(screen,BIGINT_FORMAT " surfaces in group %s\n",
              nall,gnames[igroup]);
    if (logfile)
      fprintf(logfile,BIGINT_FORMAT " surfaces in group %s\n",
              nall,gnames[igroup]);
  }
}

/* ----------------------------------------------------------------------
   add a new surface group ID, assumed to be unique
------------------------------------------------------------------------- */

int Surf::add_group(const char *id)
{
  if (ngroup == MAXGROUP) 
    error->all(FLERR,"Cannot have more than 32 surface groups");

  int n = strlen(id) + 1;
  gnames[ngroup] = new char[n];
  strcpy(gnames[ngroup],id);

  for (int i = 0; i < n-1; i++)
    if (!isalnum(id[i]) && id[i] != '_')
      error->all(FLERR,"Group ID must be alphanumeric or "
                 "underscore characters");

  ngroup++;
  return ngroup-1;
}

/* ----------------------------------------------------------------------
   find a surface group ID
   return index of group or -1 if not found
------------------------------------------------------------------------- */

int Surf::find_group(const char *id)
{
  int igroup;
  for (igroup = 0; igroup < ngroup; igroup++)
    if (strcmp(id,gnames[igroup]) == 0) break;
  if (igroup == ngroup) return -1;
  return igroup;
}

/* ----------------------------------------------------------------------
   compress owned implicit surfs to account for migrating grid cells
   migrating grid cells are ones with proc != me
   store info on reordered nlocal surfs in shash
   only called for implicit surfs
------------------------------------------------------------------------- */

void Surf::compress_rebalance()
{
  int icell;
  cellint cellID;

  if (hashfilled) hash->clear();
  hashfilled = 1;

  Grid::ChildCell *cells = grid->cells;
  Grid::MyHash *ghash = grid->hash;
  int me = comm->me;
  int n = 0;

  if (domain->dimension == 2) {
    for (int i = 0; i < nlocal; i++) {
      icell = (*ghash)[lines[i].id] - 1;
      if (cells[icell].proc != me) continue;
      if (i != n) memcpy(&lines[n],&lines[i],sizeof(Line));
      if (hash->find(lines[n].id) == hash->end())
        (*hash)[lines[n].id] = n;
      n++;
    }

  } else {
    for (int i = 0; i < nlocal; i++) {
      icell = (*ghash)[tris[i].id];
      if (cells[icell].proc != me) continue;
      if (i != n) memcpy(&tris[n],&tris[i],sizeof(Tri));
      if (hash->find(lines[n].id) == hash->end())
        (*hash)[lines[n].id] = n;
      n++;
    }
  }

  nlocal = n;
}

/* ----------------------------------------------------------------------
   reset grid csurfs values for all my owned grid cells with surfs
   called for implicit surfs after grid cell compression
   b/c local surf list was also compressed
------------------------------------------------------------------------- */

void Surf::reset_csurfs_implicit()
{
  int m,isurf,nsurf;

  Grid::ChildCell *cells = grid->cells;
  int nslocal = grid->nlocal;

  for (int icell = 0; icell < nslocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    if (cells[icell].nsurf == 0) continue;
    isurf = (*hash)[cells[icell].id];
    nsurf = cells[icell].nsurf;
    for (m = 0; m < nsurf; m++)
      cells[icell].csurfs[m] = isurf++;
  }

  // can now clear surf hash created by compress_rebalance()

  hash->clear();
  hashfilled = 0;
}

/* ----------------------------------------------------------------------
   comm of tallies across all procs
   nrow = # of tally entries in input vector
   tally2surf = surf index of each entry in input vector
   in = input vector of tallies
   instride = stride between entries in input vector
   return out = summed tallies for explicit surfs I own
------------------------------------------------------------------------- */

void Surf::collate_vector(int nrow, int *tally2surf, 
                          double *in, int instride, double *out)
{
  // collate version depends on tally_comm setting

  if (tally_comm == TALLYAUTO) {
    if (comm->nprocs > nsurf) 
      collate_vector_reduce(nrow,tally2surf,in,instride,out);
    else collate_vector_rendezvous(nrow,tally2surf,in,instride,out);
  } else if (tally_comm == TALLYREDUCE) {
    collate_vector_reduce(nrow,tally2surf,in,instride,out);
  } else if (tally_comm == TALLYRVOUS) {
    collate_vector_rendezvous(nrow,tally2surf,in,instride,out);
  }
}

/* ----------------------------------------------------------------------
   allreduce version of collate
------------------------------------------------------------------------- */

void Surf::collate_vector_reduce(int nrow, int *tally2surf, 
                                 double *in, int instride, double *out)
{
  int i,j,m;

  if (nsurf > MAXSMALLINT)
    error->all(FLERR,"Two many surfs to tally reduce - "
               "use global surf/comm auto or rvous");
  
  int nglobal = nsurf;

  double *one,*all;
  memory->create(one,nglobal,"surf:one");
  memory->create(all,nglobal,"surf:all");

  // zero all values and add in values I accumulated
  
  for (i = 0; i < nglobal; i++) one[i] = 0.0;

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  int dim = domain->dimension;
  surfint id;

  j = 0;
  for (i = 0; i < nrow; i++) {
    if (dim == 2) m = (int) lines[tally2surf[i]].id - 1;
    else m = (int) tris[tally2surf[i]].id - 1;
    one[m] = in[j];
    j += instride;
  }

  // global allreduce
  
  MPI_Allreduce(one,all,nglobal,MPI_DOUBLE,MPI_SUM,world);

  // out = only surfs I own

  int me = comm->me;
  int nprocs = comm->nprocs;

  m = 0;
  for (i = me; i < nglobal; i += nprocs)
    out[m++] = all[i];

  // NOTE: could persist these for multiple invocations
  
  memory->destroy(one);
  memory->destroy(all);
}

/* ----------------------------------------------------------------------
   rendezvous version of collate
------------------------------------------------------------------------- */

void Surf::collate_vector_rendezvous(int nrow, int *tally2surf, 
                                     double *in, int instride, double *out)
{
  // allocate memory for rvous input

  int *proclist;
  memory->create(proclist,nrow,"surf:proclist");
  InRvousVec *in_rvous = 
    (InRvousVec *) memory->smalloc((bigint) nrow*sizeof(InRvousVec),
                                   "surf:in_rvous");

  // create rvous inputs
  // proclist = owner of each surf
  // logic of (id-1) % nprocs sends 
  //   surf IDs 1,11,21,etc on 10 procs to proc 0

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  int dim = domain->dimension;
  int nprocs = comm->nprocs;
  
  surfint id;
  
  int m = 0;
  for (int i = 0; i < nrow; i++) {
    if (dim == 2) id = lines[tally2surf[i]].id;
    else id = tris[tally2surf[i]].id;
    proclist[i] = (id-1) % nprocs;
    in_rvous[i].id = id;
    in_rvous[i].value = in[m];
    m += instride;
  }

  // perform rendezvous operation
  // each proc owns subset of surfs
  // receives all tally contributions to surfs it owns

  out_rvous = out;
  
  char *buf;
  int nout = comm->rendezvous(1,nrow,(char *) in_rvous,sizeof(InRvousVec),
			      0,proclist,rendezvous_vector,
			      0,buf,0,(void *) this);

  memory->destroy(proclist);
  memory->destroy(in_rvous);
}

/* ----------------------------------------------------------------------
   callback from rendezvous operation
   process tallies for surfs assigned to me
   inbuf = list of N Inbuf datums
   no outbuf
------------------------------------------------------------------------- */

int Surf::rendezvous_vector(int n, char *inbuf, int &flag, int *&proclist,
                            char *&outbuf, void *ptr)
{
  Surf *sptr = (Surf *) ptr;
  Memory *memory = sptr->memory;
  int nown = sptr->nown;
  double *out = sptr->out_rvous;
  int nprocs = sptr->comm->nprocs;
  int me = sptr->comm->me;
  
  // zero my owned surf values

  for (int i = 0; i < nown; i++) out[i] = 0.0;
  
  // accumulate per-surf values from different procs to my owned surfs
  // logic of (id-1-me) / nprocs maps 
  //   surf IDs [1,11,21,...] on 10 procs to [0,1,2,...] on proc 0
  
  Surf::InRvousVec *in_rvous = (Surf::InRvousVec *) inbuf;

  int m;
  for (int i = 0; i < n; i++) {
    m = (in_rvous[i].id-1-me) / nprocs;
    out[m] += in_rvous[i].value;
  }

  // flag = 0: no second comm needed in rendezvous

  flag = 0;
  return 0;
}

/* ----------------------------------------------------------------------
   comm of tallies across all procs
   nrow,ncol = # of entries and columns in input array
   tally2surf = global surf index of each entry in input array
   in = input array of tallies
   instride = stride between entries in input array
   return out = summed tallies for explicit surfs I own
------------------------------------------------------------------------- */

void Surf::collate_array(int nrow, int ncol, int *tally2surf, 
                         double **in, double **out)
{
  // collate version depends on tally_comm setting

  if (tally_comm == TALLYAUTO) {
    if (comm->nprocs > nsurf) 
      collate_array_reduce(nrow,ncol,tally2surf,in,out);
    else collate_array_rendezvous(nrow,ncol,tally2surf,in,out);
  } else if (tally_comm == TALLYREDUCE) {
    collate_array_reduce(nrow,ncol,tally2surf,in,out);
  } else if (tally_comm == TALLYRVOUS) {
    collate_array_rendezvous(nrow,ncol,tally2surf,in,out);
  }
}

/* ----------------------------------------------------------------------
   allreduce version of collate
------------------------------------------------------------------------- */

void Surf::collate_array_reduce(int nrow, int ncol, int *tally2surf, 
                                double **in, double **out)
{
  int i,j,m;

  bigint ntotal = (bigint) nsurf * ncol;

  if (ntotal > MAXSMALLINT)
    error->all(FLERR,"Two many surfs to tally reduce - "
               "use global surf/comm auto or rvous");
  
  int nglobal = nsurf;

  double **one,**all;
  memory->create(one,nglobal,ncol,"surf:one");
  memory->create(all,nglobal,ncol,"surf:all");

  // zero all values and add in values I accumulated
  
  for (i = 0; i < nglobal; i++)
    for (j = 0; j < ncol; j++)
      one[i][j] = 0.0;

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  int dim = domain->dimension;

  for (i = 0; i < nrow; i++) {
    if (dim == 2) m = (int) lines[tally2surf[i]].id - 1;
    else m = (int) tris[tally2surf[i]].id - 1;
    for (j = 0; j < ncol; j++) 
      one[m][j] = in[i][j];
  }

  // global allreduce

  MPI_Allreduce(&one[0][0],&all[0][0],ntotal,MPI_DOUBLE,MPI_SUM,world);

  // out = only surfs I own

  int me = comm->me;
  int nprocs = comm->nprocs;

  m = 0;
  for (i = me; i < nglobal; i += nprocs) {
    for (j = 0; j < ncol; j++) out[m][j] = all[i][j];
    m++;
  }
  
  // NOTE: could persist these for multiple invocations
  
  memory->destroy(one);
  memory->destroy(all);
}

/* ----------------------------------------------------------------------
   rendezvous version of collate
------------------------------------------------------------------------- */

void Surf::collate_array_rendezvous(int nrow, int ncol, int *tally2surf, 
                                    double **in, double **out)
{
  int i,j,m;

  // allocate memory for rvous input

  int *proclist;
  memory->create(proclist,nrow,"surf:proclist");
  double *in_rvous = (double *)     // worry about overflow
    memory->smalloc(nrow*(ncol+1)*sizeof(double*),"surf:in_rvous");

  // create rvous inputs
  // proclist = owner of each surf
  // logic of (id-1) % nprocs sends 
  //   surf IDs 1,11,21,etc on 10 procs to proc 0

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  int dim = domain->dimension;
  int nprocs = comm->nprocs;
  surfint id;
  
  m = 0;
  for (int i = 0; i < nrow; i++) {
    if (dim == 2) id = lines[tally2surf[i]].id;
    else id = tris[tally2surf[i]].id;
    proclist[i] = (id-1) % nprocs;
    in_rvous[m++] = ubuf(id).d;
    for (j = 0; j < ncol; j++)
      in_rvous[m++] = in[i][j];
  }

  // perform rendezvous operation
  // each proc owns subset of surfs
  // receives all tally contributions to surfs it owns

  ncol_rvous = ncol;
  if (out == NULL) out_rvous = NULL;
  else out_rvous = &out[0][0];
  int size = (ncol+1) * sizeof(double);

  char *buf;
  int nout = comm->rendezvous(1,nrow,(char *) in_rvous,size,
			      0,proclist,rendezvous_array,
			      0,buf,0,(void *) this);

  memory->destroy(proclist);
  memory->sfree(in_rvous);
}

/* ----------------------------------------------------------------------
   callback from rendezvous operation
   process tallies for surfs assigned to me
   inbuf = list of N Inbuf datums
   no outbuf
------------------------------------------------------------------------- */

int Surf::rendezvous_array(int n, char *inbuf,
                           int &flag, int *&proclist, char *&outbuf,
                           void *ptr)
{
  int i,j,k,m;

  Surf *sptr = (Surf *) ptr;
  Memory *memory = sptr->memory;
  int nown = sptr->nown;
  int ncol = sptr->ncol_rvous;
  double *out = sptr->out_rvous;
  int nprocs = sptr->comm->nprocs;
  int me = sptr->comm->me;
  
  // zero my owned surf values

  int ntotal = nown*ncol;
  for (m = 0; m < ntotal; m++) out[m] = 0.0;
  
  // accumulate per-surf values from different procs to my owned surfs
  // logic of (id-1-me) / nprocs maps 
  //   surf IDs [1,11,21,...] on 10 procs to [0,1,2,...] on proc 0
  
  double *in_rvous = (double *) inbuf;
  surfint id;

  m = 0;
  for (int i = 0; i < n; i++) {
    id = (surfint) ubuf(in_rvous[m++]).i;
    k = (id-1-me) / nprocs * ncol;
    for (j = 0; j < ncol; j++)
      out[k++] += in_rvous[m++];
  }

  // flag = 0: no second comm needed in rendezvous

  flag = 0;
  return 0;
}

/* ----------------------------------------------------------------------
   proc 0 writes surf geometry to restart file
------------------------------------------------------------------------- */

void Surf::write_restart(FILE *fp)
{
  fwrite(&ngroup,sizeof(int),1,fp);

  int n;
  for (int i = 0; i < ngroup; i++) {
    n = strlen(gnames[i]) + 1;
    fwrite(&n,sizeof(int),1,fp);
    fwrite(gnames[i],sizeof(char),n,fp);
  }

  if (domain->dimension == 2) {
    fwrite(&nsurf,sizeof(int),1,fp);
    for (int i = 0; i < nsurf; i++) {
      fwrite(&lines[i].id,sizeof(surfint),1,fp);
      fwrite(&lines[i].type,sizeof(int),1,fp);
      fwrite(&lines[i].mask,sizeof(int),1,fp);
      fwrite(&lines[i].p1,sizeof(double),3,fp);
      fwrite(&lines[i].p2,sizeof(double),3,fp);
    }
  }

  if (domain->dimension == 3) {
    fwrite(&nsurf,sizeof(int),1,fp);
    for (int i = 0; i < nsurf; i++) {
      fwrite(&tris[i].id,sizeof(surfint),1,fp);
      fwrite(&tris[i].type,sizeof(int),1,fp);
      fwrite(&tris[i].mask,sizeof(int),1,fp);
      fwrite(&tris[i].p1,sizeof(double),3,fp);
      fwrite(&tris[i].p2,sizeof(double),3,fp);
      fwrite(&tris[i].p3,sizeof(double),3,fp);
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads surf geometry from restart file
   bcast to other procs
   NOTE: needs to be generalized for different surf styles
------------------------------------------------------------------------- */

void Surf::read_restart(FILE *fp)
{
  int me = comm->me;

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

  if (domain->dimension == 2) {
    if (me == 0) fread(&nsurf,sizeof(int),1,fp);
    MPI_Bcast(&nsurf,1,MPI_INT,0,world);
    lines = (Line *) memory->smalloc(nsurf*sizeof(Line),"surf:lines");

    if (me == 0) {
      for (int i = 0; i < nsurf; i++) {
        fread(&lines[i].id,sizeof(surfint),1,fp);
        fread(&lines[i].type,sizeof(int),1,fp);
        fread(&lines[i].mask,sizeof(int),1,fp);
        lines[i].isc = lines[i].isr = -1;
        fread(lines[i].p1,sizeof(double),3,fp);
        fread(lines[i].p2,sizeof(double),3,fp);
        lines[i].norm[0] = lines[i].norm[1] = lines[i].norm[2] = 0.0;
      }
    }
    MPI_Bcast(lines,nsurf*sizeof(Line),MPI_CHAR,0,world);
  }

  if (domain->dimension == 3) {
    if (me == 0) fread(&nsurf,sizeof(int),1,fp);
    MPI_Bcast(&nsurf,1,MPI_INT,0,world);
    tris = (Tri *) memory->smalloc(nsurf*sizeof(Tri),"surf:tris");

    if (me == 0) {
      for (int i = 0; i < nsurf; i++) {
        fread(&tris[i].id,sizeof(surfint),1,fp);
        fread(&tris[i].type,sizeof(int),1,fp);
        fread(&tris[i].mask,sizeof(int),1,fp);
        tris[i].isc = tris[i].isr = -1;
        fread(tris[i].p1,sizeof(double),3,fp);
        fread(tris[i].p2,sizeof(double),3,fp);
        fread(tris[i].p3,sizeof(double),3,fp);
        tris[i].norm[0] = tris[i].norm[1] = tris[i].norm[2] = 0.0;
      }
    }
    MPI_Bcast(tris,nsurf*sizeof(Tri),MPI_CHAR,0,world);
  }
}

/* ---------------------------------------------------------------------- */

void Surf::grow()
{
  if (domain->dimension == 2)
    lines = (Surf::Line *) 
      memory->srealloc(lines,nmax*sizeof(Surf::Line),"surf:lines");
  else 
    tris = (Surf::Tri *) 
      memory->srealloc(tris,nmax*sizeof(Surf::Tri),"surf:tris");
}


/* ---------------------------------------------------------------------- */

void Surf::grow_own()
{
  if (domain->dimension == 2)
    mylines = (Surf::Line *) 
      memory->srealloc(mylines,maxown*sizeof(Surf::Line),"surf:lines");
  else 
    mytris = (Surf::Tri *) 
      memory->srealloc(mytris,maxown*sizeof(Surf::Tri),"surf:tris");
}

/* ---------------------------------------------------------------------- */

bigint Surf::memory_usage()
{
  bigint bytes = 0;

  if (implicit) {
    if (domain->dimension == 2) bytes += nlocal * sizeof(Line);
    else bytes += nlocal * sizeof(Tri);
  } else if (distributed) {
    if (domain->dimension == 2) bytes += (nlocal+nghost) * sizeof(Line);
    else bytes += (nlocal+nghost) * sizeof(Tri);
    if (domain->dimension == 2) bytes += nown * sizeof(Line);
    else bytes += nown * sizeof(Tri);
  } else {
    if (domain->dimension == 2) bytes += nsurf * sizeof(Line);
    else bytes += nsurf * sizeof(Tri);
    bytes += nlocal * sizeof(int);
  }

  return bytes;
}

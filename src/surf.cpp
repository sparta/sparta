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
#include "comm.h"
#include "geometry.h"
#include "input.h"
#include "math_extra.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

#ifdef SPARTA_MAP
#include <map>
#elif defined SPARTA_UNORDERED_MAP
#include <unordered_map>
#else
#include <tr1/unordered_map>
#endif

using namespace SPARTA_NS;
using namespace MathConst;

enum{TALLYAUTO,TALLYREDUCE,TALLYLOCAL};         // same as Update
enum{REGION_ALL,REGION_ONE,REGION_CENTER};      // same as Grid
enum{TYPE,MOLECULE,ID};
enum{LT,LE,GT,GE,EQ,NEQ,BETWEEN};

#define DELTA 4
#define EPSSQ 1.0e-12
#define BIG 1.0e20
#define MAXGROUP 32

/* ---------------------------------------------------------------------- */

Surf::Surf(SPARTA *sparta) : Pointers(sparta)
{
  exist = 0;
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

  nline = ntri = 0;
  lines = NULL;
  tris = NULL;
  pushflag = 1;

  nlocal = 0;
  mysurfs = NULL;

  nsc = maxsc = 0;
  sc = NULL;
  nsr = maxsr = 0;
  sr = NULL;

  tally_comm = TALLYAUTO;
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
  memory->sfree(mysurfs);

  for (int i = 0; i < nsc; i++) delete sc[i];
  memory->sfree(sc);
  for (int i = 0; i < nsr; i++) delete sr[i];
  memory->sfree(sr);
}

/* ---------------------------------------------------------------------- */

void Surf::modify_params(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal surf_modify command");
  int igroup = find_group(arg[0]);
  if (igroup < 0) error->all(FLERR,"Surf_modify surface group is not defined");
  int groupbit = bitmask[igroup];

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"collide") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal surf_modify command");

      int isc = find_collide(arg[iarg+1]);
      if (isc < 0) error->all(FLERR,"Could not find surf_modify sc-ID");

      // set surf collision model for each surf in surface group

      if (domain->dimension == 2) {
        for (int i = 0; i < nline; i++)
          if (lines[i].mask & groupbit) lines[i].isc = isc;
      }
      if (domain->dimension == 3) {
        for (int i = 0; i < ntri; i++)
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

      if (domain->dimension == 2) {
        for (int i = 0; i < nline; i++)
          if (lines[i].mask & groupbit) lines[i].isr = isr;
      }
      if (domain->dimension == 3) {
        for (int i = 0; i < ntri; i++)
          if (tris[i].mask & groupbit) tris[i].isr = isr;
      }

      iarg += 2;

    } else error->all(FLERR,"Illegal surf_modify command");
  }
}

/* ---------------------------------------------------------------------- */

void Surf::init()
{
  // check that every element is assigned to a surf collision model
  // skip if caller turned off the check, e.g. BalanceGrid

  if (surf_collision_check) {
    int flag = 0;
    if (domain->dimension == 2) {
      for (int i = 0; i < nline; i++)
        if (lines[i].isc < 0) flag++;
    } 
    if (domain->dimension == 3) {
      for (int i = 0; i < ntri; i++)
        if (tris[i].isc < 0) flag++;
    }
    if (flag) {
      char str[64];
      sprintf(str,"%d surface elements not assigned to a collision model",flag);
      error->all(FLERR,str);
    }
  }

  // if a surf element is assigned a reaction model
  // must have a collision model that allows reactions

  if (surf_collision_check) {
    int flag = 0;
    if (domain->dimension == 2) {
      for (int i = 0; i < nline; i++)
        if (lines[i].isr >= 0 && sc[lines[i].isc]->allowreact == 0) flag++;
    } 
    if (domain->dimension == 3) {
      for (int i = 0; i < ntri; i++)
        if (tris[i].isr >= 0 && sc[tris[i].isc]->allowreact == 0) flag++;
    } 
    if (flag) {
      char str[64];
      sprintf(str,"%d surface elements with reaction model, "
              "but invalid collision model",flag);
      error->all(FLERR,str);
    }
  }    

  // initialize surf collision and reaction models

  for (int i = 0; i < nsc; i++) sc[i]->init();
  for (int i = 0; i < nsr; i++) sr[i]->init();
}

/* ----------------------------------------------------------------------
   return # of lines or triangles
------------------------------------------------------------------------- */

int Surf::nelement()
{
  if (domain->dimension == 2) return nline;
  return ntri;
}

/* ----------------------------------------------------------------------
   setup owned surf elements
   create mysurfs list of owned surfs
   compute local index for owned cells
------------------------------------------------------------------------- */

void Surf::setup_surf()
{
  int i,j;

  int me = comm->me;
  int nprocs = comm->nprocs;

  int n = nelement();

  // assign every Pth surf element to this proc

  nlocal = n/nprocs;
  if (me < n % nprocs) nlocal++;

  memory->destroy(mysurfs);
  memory->create(mysurfs,nlocal,"surf:mysurfs");

  nlocal = 0;
  for (int m = me; m < n; m += nprocs)
    mysurfs[nlocal++] = m;

  // set bounding box of all surfs based on their pts
  // for 2d, set zlo,zhi to box bounds

  for (j = 0; j < 3; j++) {
    bblo[j] = BIG;
    bbhi[j] = -BIG;
  }
  
  double *x;

  if (domain->dimension == 2) {
    for (i = 0; i < nline; i++) {
      x = lines[i].p1;
      for (j = 0; j < 2; j++) {
        bblo[j] = MIN(bblo[j],x[j]);
        bbhi[j] = MAX(bbhi[j],x[j]);
      }
      x = lines[i].p2;
      for (j = 0; j < 2; j++) {
        bblo[j] = MIN(bblo[j],x[j]);
        bbhi[j] = MAX(bbhi[j],x[j]);
      }
    }
    bblo[2] = domain->boxlo[2];
    bbhi[2] = domain->boxhi[2];

  } else if (domain->dimension == 3) {
    for (i = 0; i < ntri; i++) {
      x = tris[i].p1;
      for (j = 0; j < 3; j++) {
        bblo[j] = MIN(bblo[j],x[j]);
        bbhi[j] = MAX(bbhi[j],x[j]);
      }
      x = tris[i].p2;
      for (j = 0; j < 3; j++) {
        bblo[j] = MIN(bblo[j],x[j]);
        bbhi[j] = MAX(bbhi[j],x[j]);
      }
      x = tris[i].p3;
      for (j = 0; j < 3; j++) {
        bblo[j] = MIN(bblo[j],x[j]);
        bbhi[j] = MAX(bbhi[j],x[j]);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   compute unit outward normal vectors of all lines starting at Nold
   outward normal = +z axis x (p2-p1)
------------------------------------------------------------------------- */

void Surf::compute_line_normal(int nold)
{
  double z[3],delta[3];

  z[0] = 0.0; z[1] = 0.0; z[2] = 1.0;

  for (int i = nold; i < nline; i++) {
    MathExtra::sub3(lines[i].p2,lines[i].p1,delta);
    MathExtra::cross3(z,delta,lines[i].norm);
    MathExtra::norm3(lines[i].norm);
    lines[i].norm[2] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   compute unit outward normal vectors of all lines starting at Nold
   outward normal = (p2-p1) x (p3-p1)
------------------------------------------------------------------------- */

void Surf::compute_tri_normal(int nold)
{
  int p1,p2,p3;
  double delta12[3],delta13[3];

  for (int i = nold; i < ntri; i++) {
    MathExtra::sub3(tris[i].p2,tris[i].p1,delta12);
    MathExtra::sub3(tris[i].p3,tris[i].p1,delta13);
    MathExtra::cross3(delta12,delta13,tris[i].norm);
    MathExtra::norm3(tris[i].norm);
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
   return length of line M
------------------------------------------------------------------------- */

double Surf::line_size(int m)
{
  return line_size(lines[m].p1,lines[m].p2);
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
   return area associated with rotating axisymmetric line M around y=0 axis
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
   compute side length and area of triangle M
   return len = length of shortest edge of triangle M
   return area = area of triangle M
------------------------------------------------------------------------- */

double Surf::tri_size(int m, double &len)
{
  return tri_size(tris[m].p1,tris[m].p2,tris[m].p3,len);
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
  // hash end points of all lines
  // key = end point
  // value = 1 if first point, 2 if second point, 3 if both points
  // NOTE: could prealloc hash to correct size here

#ifdef SPARTA_MAP
  std::map<OnePoint2d,int> hash;
  std::map<OnePoint2d,int>::iterator it;
#elif defined SPARTA_UNORDERED_MAP
  std::unordered_map<OnePoint2d,int,Hasher2d> hash;
  std::unordered_map<OnePoint2d,int,Hasher2d>::iterator it;
#else
  std::tr1::unordered_map<OnePoint2d,int,Hasher2d> hash;
  std::tr1::unordered_map<OnePoint2d,int,Hasher2d>::iterator it;
#endif

  // insert each end point into hash
  // should appear once at each end
  // error if any duplicate points

  double *p1,*p2;
  OnePoint2d key;
  int value;

  int nline_new = nline - nline_old;

  int ndup = 0;
  int m = nline_old;
  for (int i = 0; i < nline_new; i++) {
    p1 = lines[m].p1;
    key.pt[0] = p1[0]; key.pt[1] = p1[1];
    if (hash.find(key) == hash.end()) hash[key] = 1;
    else {
      value = hash[key];
      if (value == 2) hash[key] = 3;
      else ndup++;
    }

    p2 = lines[m].p2;
    key.pt[0] = p2[0]; key.pt[1] = p2[1];
    if (hash.find(key) == hash.end()) hash[key] = 2;
    else {
      value = hash[key];
      if (value == 1) hash[key] = 3;
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
  for (it = hash.begin(); it != hash.end(); ++it) {
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
  // hash directed edges of all triangles
  // key = directed edge
  // value = 1 if appears once, 2 if reverse also appears once
  // NOTE: could prealloc hash to correct size here

#ifdef SPARTA_MAP
  std::map<TwoPoint3d,int> hash;
  std::map<TwoPoint3d,int>::iterator it;
#elif defined SPARTA_UNORDERED_MAP
  std::unordered_map<TwoPoint3d,int,Hasher3d> hash;
  std::unordered_map<TwoPoint3d,int,Hasher3d>::iterator it;
#else
  std::tr1::unordered_map<TwoPoint3d,int,Hasher3d> hash;
  std::tr1::unordered_map<TwoPoint3d,int,Hasher3d>::iterator it;
#endif

  // insert each edge into hash
  // should appear once in each direction
  // error if any duplicate edges
  
  double *p1,*p2,*p3;
  TwoPoint3d key,keyinv;
  int value;

  int ntri_new = ntri - ntri_old;

  int ndup = 0;
  int m = ntri_old;
  for (int i = 0; i < ntri_new; i++) {
    p1 = tris[m].p1;
    p2 = tris[m].p2;
    p3 = tris[m].p3;

    key.pts[0] = p1[0]; key.pts[1] = p1[1]; key.pts[2] = p1[2];
    key.pts[3] = p2[0]; key.pts[4] = p2[1]; key.pts[5] = p2[2];
    if (hash.find(key) == hash.end()) {
      keyinv.pts[0] = p2[0]; keyinv.pts[1] = p2[1]; keyinv.pts[2] = p2[2];
      keyinv.pts[3] = p1[0]; keyinv.pts[4] = p1[1]; keyinv.pts[5] = p1[2];
      if (hash.find(keyinv) == hash.end()) hash[key] = 1;
      else {
	value = hash[keyinv];
	if (value == 1) hash[keyinv] = 2;
	else ndup++;
      }
    } else ndup++;
    
    key.pts[0] = p2[0]; key.pts[1] = p2[1]; key.pts[2] = p2[2];
    key.pts[3] = p3[0]; key.pts[4] = p3[1]; key.pts[5] = p3[2];
    if (hash.find(key) == hash.end()) {
      keyinv.pts[0] = p3[0]; keyinv.pts[1] = p3[1]; keyinv.pts[2] = p3[2];
      keyinv.pts[3] = p2[0]; keyinv.pts[4] = p2[1]; keyinv.pts[5] = p2[2];
      if (hash.find(keyinv) == hash.end()) hash[key] = 1;
      else {
	value = hash[keyinv];
	if (value == 1) hash[keyinv] = 2;
	else ndup++;
      }
    } else ndup++;

    key.pts[0] = p3[0]; key.pts[1] = p3[1]; key.pts[2] = p3[2];
    key.pts[3] = p1[0]; key.pts[4] = p1[1]; key.pts[5] = p1[2];
    if (hash.find(key) == hash.end()) {
      keyinv.pts[0] = p1[0]; keyinv.pts[1] = p1[1]; keyinv.pts[2] = p1[2];
      keyinv.pts[3] = p3[0]; keyinv.pts[4] = p3[1]; keyinv.pts[5] = p3[2];
      if (hash.find(keyinv) == hash.end()) hash[key] = 1;
      else {
	value = hash[keyinv];
	if (value == 1) hash[keyinv] = 2;
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
  for (it = hash.begin(); it != hash.end(); ++it) {
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
   nold = previous # of elements
   nnew = current # of elements
------------------------------------------------------------------------- */

void Surf::check_point_inside(int nold)
{
  int nbad;
  double *x;

  int dim = domain->dimension;
  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  if (dim == 2) {
    int nline_new = nline - nold;
    int m = nold;
    nbad = 0;
    for (int i = 0; i < nline_new; i++) {
      x = lines[m].p1;
      if (x[0] < boxlo[0] || x[0] > boxhi[0] ||
	  x[1] < boxlo[1] || x[1] > boxhi[1] ||
	  x[2] < boxlo[2] || x[2] > boxhi[2]) nbad++;
      x = lines[m].p2;
      if (x[0] < boxlo[0] || x[0] > boxhi[0] ||
	  x[1] < boxlo[1] || x[1] > boxhi[1] ||
	  x[2] < boxlo[2] || x[2] > boxhi[2]) nbad++;
      m++;
    }
  } else if (dim == 3) {
    int ntri_new = ntri - nold;
    int m = nold;
    nbad = 0;
    for (int i = 0; i < ntri_new; i++) {
      x = tris[m].p1;
      if (x[0] < boxlo[0] || x[0] > boxhi[0] ||
	  x[1] < boxlo[1] || x[1] > boxhi[1] ||
	  x[2] < boxlo[2] || x[2] > boxhi[2]) nbad++;
      x = tris[m].p2;
      if (x[0] < boxlo[0] || x[0] > boxhi[0] ||
	  x[1] < boxlo[1] || x[1] > boxhi[1] ||
	  x[2] < boxlo[2] || x[2] > boxhi[2]) nbad++;
      x = tris[m].p3;
      if (x[0] < boxlo[0] || x[0] > boxhi[0] ||
	  x[1] < boxlo[1] || x[1] > boxhi[1] ||
	  x[2] < boxlo[2] || x[2] > boxhi[2]) nbad++;
      m++;
    }
  }

  if (nbad) {
    char str[128];
    sprintf(str,"%d surface points are not inside simulation box",
	    nbad);
    error->all(FLERR,str);
  }
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
    maxsc += DELTA;
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
    maxsr += DELTA;
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
------------------------------------------------------------------------- */

void Surf::group(int narg, char **arg)
{
  int i,flag;
  double x[3];

  if (narg < 3) error->all(FLERR,"Illegal group command");

  int dimension = domain->dimension;

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
        printf("COND %d %d\n",condition,BETWEEN);
        if (condition == LT) {
          if (dimension == 2) {
            for (i = 0; i < nline; i++) 
              if (i+1 < bound1) lines[i].mask |= bit;
          } else {
            for (i = 0; i < ntri; i++) 
              if (i+1 < bound1) tris[i].mask |= bit;
          }
        } else if (condition == LE) {
          if (dimension == 2) {
            for (i = 0; i < nline; i++) 
              if (i+1 <= bound1) lines[i].mask |= bit;
          } else {
            for (i = 0; i < ntri; i++) 
              if (i+1 <= bound1) tris[i].mask |= bit;
          }
        } else if (condition == GT) {
          if (dimension == 2) {
            for (i = 0; i < nline; i++) 
              if (i+1 > bound1) lines[i].mask |= bit;
          } else {
            for (i = 0; i < ntri; i++) 
              if (i+1 > bound1) tris[i].mask |= bit;
          }
        } else if (condition == GE) {
          if (dimension == 2) {
            for (i = 0; i < nline; i++) 
              if (i+1 >= bound1) lines[i].mask |= bit;
          } else {
            for (i = 0; i < ntri; i++) 
              if (i+1 >= bound1) tris[i].mask |= bit;
          }
        } else if (condition == EQ) {
          if (dimension == 2) {
            for (i = 0; i < nline; i++) 
              if (i+1 == bound1) lines[i].mask |= bit;
          } else {
            for (i = 0; i < ntri; i++) 
              if (i+1 == bound1) tris[i].mask |= bit;
          }
        } else if (condition == NEQ) {
          if (dimension == 2) {
            for (i = 0; i < nline; i++) 
              if (i+1 != bound1) lines[i].mask |= bit;
          } else {
            for (i = 0; i < ntri; i++) 
              if (i+1 != bound1) tris[i].mask |= bit;
          }
        } else if (condition == BETWEEN) {
          if (dimension == 2) {
            for (i = 0; i < nline; i++)
              if (i+1 >= bound1 && i+1 <= bound2) lines[i].mask |= bit;
          } else {
            for (i = 0; i < ntri; i++)
              if (i+1 >= bound1 && i+1 <= bound2) tris[i].mask |= bit;
          }
        }
      } else if (category == TYPE) {
        if (condition == LT) {
          if (dimension == 2) {
            for (i = 0; i < nline; i++) 
              if (lines[i].type < bound1) lines[i].mask |= bit;
          } else {
            for (i = 0; i < ntri; i++) 
              if (tris[i].type < bound1) lines[i].mask |= bit;
          }
        } else if (condition == LE) {
          if (dimension == 2) {
            for (i = 0; i < nline; i++) 
              if (lines[i].type <= bound1) lines[i].mask |= bit;
          } else {
            for (i = 0; i < ntri; i++) 
              if (tris[i].type <= bound1) lines[i].mask |= bit;
          }
        } else if (condition == GT) {
          if (dimension == 2) {
            for (i = 0; i < nline; i++) 
              if (lines[i].type > bound1) lines[i].mask |= bit;
          } else {
            for (i = 0; i < ntri; i++) 
              if (tris[i].type > bound1) lines[i].mask |= bit;
          }
        } else if (condition == GE) {
          if (dimension == 2) {
            for (i = 0; i < nline; i++) 
              if (lines[i].type >= bound1) lines[i].mask |= bit;
          } else {
            for (i = 0; i < ntri; i++) 
              if (tris[i].type >= bound1) lines[i].mask |= bit;
          }
        } else if (condition == EQ) {
          if (dimension == 2) {
            for (i = 0; i < nline; i++) 
              if (lines[i].type == bound1) lines[i].mask |= bit;
          } else {
            for (i = 0; i < ntri; i++) 
              if (tris[i].type == bound1) lines[i].mask |= bit;
          }
        } else if (condition == NEQ) {
          if (dimension == 2) {
            for (i = 0; i < nline; i++) 
              if (lines[i].type != bound1) lines[i].mask |= bit;
          } else {
            for (i = 0; i < ntri; i++) 
              if (tris[i].type != bound1) lines[i].mask |= bit;
          }
        } else if (condition == BETWEEN) {
          if (dimension == 2) {
            for (i = 0; i < nline; i++)
              if (lines[i].type >= bound1 && lines[i].type <= bound2) 
                lines[i].mask |= bit;
          } else {
            for (i = 0; i < ntri; i++)
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
          if (dimension == 2) {
            for (i = 0; i < nline; i++)
              if (i+1 >= start && i+1 <= stop) lines[i].mask |= bit;
          } else {
            for (i = 0; i < ntri; i++)
              if (i+1 >= start && i+1 <= stop) tris[i].mask |= bit;
          }
        } else if (category == TYPE) {
          if (dimension == 2) {
            for (i = 0; i < nline; i++)
              if (lines[i].type >= start && lines[i].type <= stop) 
                lines[i].mask |= bit;
          } else {
            for (i = 0; i < ntri; i++)
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

    if (dimension == 2) {
      if (rstyle == REGION_ALL) {
        for (i = 0; i < nline; i++) {
          flag = 1;
          if (!region->match(lines[i].p1)) flag = 0;
          if (!region->match(lines[i].p2)) flag = 0;
          if (flag) lines[i].mask |= bit;
        }
      } else if (rstyle == REGION_ONE) {
        for (i = 0; i < nline; i++) {
          flag = 0;
          if (region->match(lines[i].p1)) flag = 1;
          if (region->match(lines[i].p2)) flag = 1;
          if (flag) lines[i].mask |= bit;
        }
      } else if (rstyle == REGION_CENTER) {
        for (i = 0; i < nline; i++) {
          x[0] = 0.5 * (lines[i].p1[0] + lines[i].p2[0]);
          x[1] = 0.5 * (lines[i].p1[1] + lines[i].p2[1]);
          x[2] = 0.0;
          if (region->match(x)) lines[i].mask |= bit;
        }
      }

    } else if (dimension == 3) {
      if (rstyle == REGION_ALL) {
        for (i = 0; i < ntri; i++) {
          flag = 1;
          if (!region->match(tris[i].p1)) flag = 0;
          if (!region->match(tris[i].p2)) flag = 0;
          if (!region->match(tris[i].p3)) flag = 0;
          if (flag) tris[i].mask |= bit;
        }
      } else if (rstyle == REGION_ONE) {
        for (i = 0; i < ntri; i++) {
          flag = 0;
          if (region->match(tris[i].p1)) flag = 1;
          if (region->match(tris[i].p2)) flag = 1;
          if (region->match(tris[i].p3)) flag = 1;
          if (flag) tris[i].mask |= bit;
        }
      } else if (rstyle == REGION_CENTER) {
        for (i = 0; i < ntri; i++) {
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

    if (dimension == 2) {
      for (i = 0; i < nline; i++)
        if (lines[i].mask & otherbit) lines[i].mask |= bit;
    } else {
      for (i = 0; i < ntri; i++)
        if (tris[i].mask & otherbit) tris[i].mask |= bit;
    }

    // remove surfs if they are in any of the other groups
    // AND with inverse mask removes the surf from group

    int inverse = inversemask[igroup];

    for (int ilist = 1; ilist < length; ilist++) {
      otherbit = bitmask[list[ilist]];
      if (dimension == 2) {
        for (i = 0; i < nline; i++)
          if (lines[i].mask & otherbit) lines[i].mask &= inverse;
      } else {
        for (i = 0; i < ntri; i++)
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
      if (dimension == 2) {
        for (i = 0; i < nline; i++)
          if (lines[i].mask & otherbit) lines[i].mask |= bit;
      } else {
        for (i = 0; i < ntri; i++)
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

    if (dimension == 2) {
      for (i = 0; i < nline; i++) {
        ok = 1;
        for (ilist = 0; ilist < length; ilist++) {
          otherbit = bitmask[list[ilist]];
          if ((lines[i].mask & otherbit) == 0) ok = 0;
        }
        if (ok) lines[i].mask |= bit;
      }
    } else {
      for (i = 0; i < ntri; i++) {
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

    if (dimension == 2) {
      for (i = 0; i < nline; i++) lines[i].mask &= inversebits;
    } else {
      for (i = 0; i < nline; i++) tris[i].mask &= inversebits;
    }
  }

  // print stats for changed group

  int n = 0;
  if (dimension == 2) {
    for (i = 0; i < nline; i++) 
      if (lines[i].mask & bit) n++;
  } else {
    for (i = 0; i < ntri; i++)
      if (tris[i].mask & bit) n++;
  }

  if (comm->me == 0) {
    if (screen) 
      fprintf(screen,"%d surfaces in group %s\n",n,gnames[igroup]);
    if (logfile)
      fprintf(logfile,"%d surfaces in group %s\n",n,gnames[igroup]);
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
   comm of vector of local tallies across all procs
   nrow = # of entries in input vector
   l2g = global surf index of each entry in input vector
   in = input vector
   instride = stride between entries in input vector
   return out = summed tallies for nlocal surfs I own
------------------------------------------------------------------------- */

void Surf::collate_vector(int nrow, int *l2g, 
                          double *in, int instride, double *out)
{
  collate_vector_allreduce(nrow,l2g,in,instride,out);

  //if (tally_comm == TALLYAUTO) 
  //  collate_vector_allreduce(nrow,l2g,in,instride,out);
  //else
  //  collate_vector_irregular(nrow,l2g,in,instride,out);
}

void Surf::collate_vector_allreduce(int nrow, int *l2g, 
                                    double *in, int instride, double *out)
{
  int i,m;

  int nglobal;
  if (domain->dimension == 2) nglobal = nline;
  else nglobal = ntri;
  if (nglobal == 0) return;

  double *one,*all;
  memory->create(one,nglobal,"surf:one");
  memory->create(all,nglobal,"surf:all");

  for (i = 0; i < nglobal; i++) one[i] = 0.0;

  m = 0;
  for (i = 0; i < nrow; i++) {
    one[l2g[i]] = in[m];
    m += instride;
  }

  MPI_Allreduce(one,all,nglobal,MPI_DOUBLE,MPI_SUM,world);

  for (i = 0; i < nlocal; i++) out[i] = all[mysurfs[i]];

  // NOTE: don't need to destroy them ?

  memory->destroy(one);
  memory->destroy(all);
}

void Surf::collate_vector_irregular(int, int *, 
                                    double *, int, double *)
{
}

/* ----------------------------------------------------------------------
   comm of array of local tallies across all procs
   nrow,ncol = # of entries and columns in input array
   l2g = global surf index of each entry in input vector
   in = input vector
   instride = stride between entries in input vector
   return out = summed tallies for nlocal surfs I own
------------------------------------------------------------------------- */

void Surf::collate_array(int nrow, int ncol, int *l2g, 
                         double **in, double **out)
{
  //if (tally_comm == TALLYAUTO) 
    collate_array_allreduce(nrow,ncol,l2g,in,out);
    //else
    // collate_array_irregular(nrow,ncol,l2g,in,out);
}

void Surf::collate_array_allreduce(int nrow, int ncol, int *l2g, 
                                   double **in, double **out)
{
  int i,j,m;

  int nglobal;
  if (domain->dimension == 2) nglobal = nline;
  else nglobal = ntri;
  if (nglobal == 0) return;

  double **one,**all;
  memory->create(one,nglobal,ncol,"surf:one");
  memory->create(all,nglobal,ncol,"surf:all");

  for (i = 0; i < nglobal; i++)
    for (j = 0; j < ncol; j++)
      one[i][j] = 0.0;

  for (i = 0; i < nrow; i++) {
    m = l2g[i];
    for (j = 0; j < ncol; j++) 
      one[m][j] = in[i][j];
  }

  MPI_Allreduce(&one[0][0],&all[0][0],nglobal*ncol,MPI_DOUBLE,MPI_SUM,world);

  for (i = 0; i < nlocal; i++) {
    m = mysurfs[i];
    for (j = 0; j < ncol; j++) 
      out[i][j] += all[m][j];
  }
  
  // NOTE: don't need to destroy them

  memory->destroy(one);
  memory->destroy(all);
}

void Surf::collate_array_irregular(int, int, int *, 
                                   double **, double **)
{
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
    fwrite(&nline,sizeof(int),1,fp);
    for (int i = 0; i < nline; i++) {
      fwrite(&lines[i].id,sizeof(surfint),1,fp);
      fwrite(&lines[i].type,sizeof(int),1,fp);
      fwrite(&lines[i].mask,sizeof(int),1,fp);
      fwrite(&lines[i].p1,sizeof(double),3,fp);
      fwrite(&lines[i].p2,sizeof(double),3,fp);
    }
  }
  if (domain->dimension == 3) {
    fwrite(&ntri,sizeof(int),1,fp);
    for (int i = 0; i < ntri; i++) {
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
    if (me == 0) fread(&nline,sizeof(int),1,fp);
    MPI_Bcast(&nline,1,MPI_INT,0,world);
    lines = (Line *) memory->smalloc(nline*sizeof(Line),"surf:lines");

    if (me == 0) {
      for (int i = 0; i < nline; i++) {
        fread(&lines[i].id,sizeof(surfint),1,fp);
        fread(&lines[i].type,sizeof(int),1,fp);
        fread(&lines[i].mask,sizeof(int),1,fp);
        lines[i].isc = lines[i].isr = -1;
        fread(lines[i].p1,sizeof(double),3,fp);
        fread(lines[i].p2,sizeof(double),3,fp);
        lines[i].norm[0] = lines[i].norm[1] = lines[i].norm[2] = 0.0;
      }
    }
    MPI_Bcast(lines,nline*sizeof(Line),MPI_CHAR,0,world);
  }

  if (domain->dimension == 3) {
    if (me == 0) fread(&ntri,sizeof(int),1,fp);
    MPI_Bcast(&ntri,1,MPI_INT,0,world);
    tris = (Tri *) memory->smalloc(ntri*sizeof(Tri),"surf:tris");

    if (me == 0) {
      for (int i = 0; i < ntri; i++) {
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
    MPI_Bcast(tris,ntri*sizeof(Tri),MPI_CHAR,0,world);
  }
}

/* ---------------------------------------------------------------------- */

bigint Surf::memory_usage()
{
  bigint bytes = 0;
  bytes += (bigint) nline * sizeof(Line);
  bytes += (bigint) ntri * sizeof(Tri);
  bytes += nlocal * sizeof(int);
  return bytes;
}

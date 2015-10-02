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
#include "cut3d.h"
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
enum{REGION_ALL,REGION_ONE,REGION_CENTER};
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

  npoint = nline = ntri = 0;
  pts = NULL;
  lines = NULL;
  tris = NULL;
  pushflag = 0;

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

  memory->sfree(pts);
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
  // skip if caller turned off the check, e.g. ReadRestart

  if (surf_collision_assign_check) {
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

  // set bounding box of all surfs based on pts
  // for 2d, set zlo,zhi to box bounds

  int i;
  for (i = 0; i < 3; i++) {
    bblo[i] = BIG;
    bbhi[i] = -BIG;
  }
  
  double *x;
  for (int ipt = 0; ipt < npoint; ipt++) {
    x = pts[ipt].x;
    for (i = 0; i < 3; i++) {
      bblo[i] = MIN(bblo[i],x[i]);
      bbhi[i] = MAX(bbhi[i],x[i]);
    }
  }

  if (domain->dimension == 2) {
    bblo[2] = domain->boxlo[2];
    bbhi[2] = domain->boxhi[2];
  }
}

/* ----------------------------------------------------------------------
   compute unit outward normal vectors of N lines starting at Nstart
   outward normal = +z axis x (p2-p1)
------------------------------------------------------------------------- */

void Surf::compute_line_normal(int nstart, int n)
{
  int p1,p2;
  double z[3],delta[3],norm[3];

  z[0] = 0.0; z[1] = 0.0; z[2] = 1.0;

  int m = nstart;
  for (int i = 0; i < n; i++) {
    p1 = lines[m].p1;
    p2 = lines[m].p2;
    MathExtra::sub3(pts[p2].x,pts[p1].x,delta);
    MathExtra::cross3(z,delta,norm);
    MathExtra::norm3(norm);
    lines[m].norm[0] = norm[0];
    lines[m].norm[1] = norm[1];
    lines[m].norm[2] = 0.0;
    m++;
  }
}

/* ----------------------------------------------------------------------
   compute unit outward normal vectors of N triangles starting at Nstart
   outward normal = (p2-p1) x (p3-p1)
------------------------------------------------------------------------- */

void Surf::compute_tri_normal(int nstart, int n)
{
  int p1,p2,p3;
  double delta12[3],delta13[3];

  int m = nstart;
  for (int i = 0; i < n; i++) {
    p1 = tris[m].p1;
    p2 = tris[m].p2;
    p3 = tris[m].p3;
    MathExtra::sub3(pts[p2].x,pts[p1].x,delta12);
    MathExtra::sub3(pts[p3].x,pts[p1].x,delta13);
    MathExtra::cross3(delta12,delta13,tris[m].norm);
    MathExtra::norm3(tris[m].norm);
    m++;
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
  double delta[3];
  MathExtra::sub3(pts[lines[m].p2].x,pts[lines[m].p1].x,delta);
  return MathExtra::len3(delta);
}

/* ----------------------------------------------------------------------
   return area associated with rotating axisymmetric line M around y=0 axis
------------------------------------------------------------------------- */

double Surf::axi_line_size(int m)
{
  double *x1 = pts[lines[m].p1].x;
  double *x2 = pts[lines[m].p2].x;
  double h = x2[0]-x1[0];
  double r = x2[1]-x1[1];
  double area = MY_PI*(x1[1]+x2[1])*sqrt(r*r+h*h);
  return area;
}

/* ----------------------------------------------------------------------
   compute side length and area of a triangle
   return len = length of shortest edge of triangle M
   return area = area of triangle M
------------------------------------------------------------------------- */

double Surf::tri_size(int m, double &len)
{
  double delta12[3],delta13[3],delta23[3],cross[3];

  MathExtra::sub3(pts[tris[m].p2].x,pts[tris[m].p1].x,delta12);
  MathExtra::sub3(pts[tris[m].p3].x,pts[tris[m].p1].x,delta13);
  MathExtra::sub3(pts[tris[m].p3].x,pts[tris[m].p2].x,delta23);
  len = MIN(MathExtra::len3(delta12),MathExtra::len3(delta13));
  len = MIN(len,MathExtra::len3(delta23));

  MathExtra::cross3(delta12,delta13,cross);
  double area = 0.5 * MathExtra::len3(cross);
  return area;
}

/* ----------------------------------------------------------------------
   check that points are each end point of exactly 2 new lines
   exception: not required of point on simulation box surface
   only check points and lines newer than old indices
------------------------------------------------------------------------- */

void Surf::check_watertight_2d(int npoint_old, int nline_old)
{
  int p1,p2;

  int npoint_new = npoint - npoint_old;
  int nline_new = nline - nline_old;

  // count[I] = # of lines that vertex I is part of

  int *count;
  memory->create(count,npoint_new,"readsurf:count");
  for (int i = 0; i < npoint_new; i++) count[i] = 0;

  int ndup = 0;
  int m = nline_old;
  for (int i = 0; i < nline_new; i++) {
    p1 = lines[m].p1 - npoint_old;
    p2 = lines[m].p2 - npoint_old;
    count[p1]++;
    count[p2]++;
    if (count[p1] > 2) ndup++;
    if (count[p2] > 2) ndup++;
    m++;
  }
  
  if (ndup) {
    char str[128];
    sprintf(str,"Surface check failed with %d duplicate points",ndup);
    error->all(FLERR,str);
  }

  // check that all counts are 2
  // allow for exception if point on box surface

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  int nbad = 0;
  for (int i = 0; i < npoint_new; i++) {
    if (count[i] == 0) nbad++;
    else if (count[i] == 1) {
      if (!Geometry::point_on_hex(pts[i+npoint_old].x,boxlo,boxhi)) nbad++;
    }
  }

  if (nbad) {
    char str[128];
    sprintf(str,"Surface check failed with %d unmatched points",nbad);
    error->all(FLERR,str);
  }

  // clean up

  memory->destroy(count);
}

/* ----------------------------------------------------------------------
   check directed triangle edges
   must be unique and match exactly one inverted edge
   exception: not required of triangle edge on simulation box surface
   only check points and lines newer than old indices
------------------------------------------------------------------------- */

void Surf::check_watertight_3d(int npoint_old, int ntri_old)
{
  int ntri_new = ntri - ntri_old;

  // hash directed edges of all triangles
  // key = directed edge, value = # of times it appears in any triangle
  // NOTE: could prealloc hash to correct size here

#ifdef SPARTA_MAP
  std::map<bigint,int> hash;
  std::map<bigint,int>::iterator it;
#elif defined SPARTA_UNORDERED_MAP
  std::unordered_map<bigint,int> hash;
  std::unordered_map<bigint,int>::iterator it;
#else
  std::tr1::unordered_map<bigint,int> hash;
  std::tr1::unordered_map<bigint,int>::iterator it;
#endif

  // insert each edge into hash
  // should appear once in each direction
  // error if any duplicate edges

  bigint p1,p2,p3,key;

  int ndup = 0;
  int m = ntri_old;
  for (int i = 0; i < ntri_new; i++) {
    p1 = tris[m].p1;
    p2 = tris[m].p2;
    p3 = tris[m].p3;
    key = (p1 << 32) | p2;
    if (hash.find(key) != hash.end()) ndup++;
    else hash[key] = 0;
    key = (p2 << 32) | p3;
    if (hash.find(key) != hash.end()) ndup++;
    else hash[key] = 0;
    key = (p3 << 32) | p1;
    if (hash.find(key) != hash.end()) ndup++;
    else hash[key] = 0;
    m++;
  }

  if (ndup) {
    char str[128];
    sprintf(str,"Surface check failed with %d duplicate edges",ndup);
    error->all(FLERR,str);
  }

  // check that each edge has an inverted match
  // allow for exception if edge on box surface

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  int nbad = 0;
  for (it = hash.begin(); it != hash.end(); ++it) {
    p1 = it->first >> 32;
    p2 = it->first & MAXSMALLINT;
    key = (p2 << 32) | p1;
    if (hash.find(key) == hash.end()) {
      if (!Geometry::point_on_hex(pts[p1].x,boxlo,boxhi) ||
          !Geometry::point_on_hex(pts[p2].x,boxlo,boxhi)) nbad++;
    }
  }

  if (nbad) {
    char str[128];
    sprintf(str,"Surface check failed with %d unmatched edges",nbad);
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
   check if all points are inside or on surface of global simulation box
------------------------------------------------------------------------- */

void Surf::check_point_inside(int npoint_old, int npoint_new)
{
  double *x;

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  int m = npoint_old;
  int nbad = 0;
  for (int i = 0; i < npoint_new; i++) {
    x = pts[m].x;
    if (x[0] < boxlo[0] || x[0] > boxhi[0] ||
	x[1] < boxlo[1] || x[1] > boxhi[1] ||
	x[2] < boxlo[2] || x[2] > boxhi[2]) nbad++;
    m++;
  }

  if (nbad) {
    char str[128];
    sprintf(str,"%d surface points are not inside simulation box",
	    nbad);
    error->all(FLERR,str);
  }
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
    else if (strcmp(arg[1],"id") == 0) category = ID;

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

      if (category == TYPE) {
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
      } else if (category == ID) {
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
          start = input->inumeric(FLERR,ptr); 
          *ptr = ':';
          stop = input->inumeric(FLERR,ptr+1); 
        } else {
          start = stop = input->inumeric(FLERR,arg[iarg]);
        }

        // add surf to group if type/id matches value or sequence
      
        if (category == TYPE) {
          if (dimension == 2) {
            for (i = 0; i < nline; i++)
              if (i+1 >= start && i+1 <= stop) lines[i].mask |= bit;
          } else {
            for (i = 0; i < ntri; i++)
              if (i+1 >= start && i+1 <= stop) tris[i].mask |= bit;
          }
        } else if (category == ID) {
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
          if (!region->match(pts[lines[i].p1].x)) flag = 0;
          if (!region->match(pts[lines[i].p2].x)) flag = 0;
          if (flag) lines[i].mask |= bit;
        }
      } else if (rstyle == REGION_ONE) {
        for (i = 0; i < nline; i++) {
          flag = 0;
          if (region->match(pts[lines[i].p1].x)) flag = 1;
          if (region->match(pts[lines[i].p2].x)) flag = 1;
          if (flag) lines[i].mask |= bit;
        }
      } else if (rstyle == REGION_CENTER) {
        for (i = 0; i < nline; i++) {
          x[0] = 0.5 * (pts[lines[i].p1].x[0] + pts[lines[i].p2].x[0]);
          x[1] = 0.5 * (pts[lines[i].p1].x[1] + pts[lines[i].p2].x[1]);
          x[2] = 0.0;
          if (region->match(x)) lines[i].mask |= bit;
        }
      }

    } else if (dimension == 3) {
      if (rstyle == REGION_ALL) {
        for (i = 0; i < ntri; i++) {
          flag = 1;
          if (!region->match(pts[tris[i].p1].x)) flag = 0;
          if (!region->match(pts[tris[i].p2].x)) flag = 0;
          if (!region->match(pts[tris[i].p3].x)) flag = 0;
          if (flag) tris[i].mask |= bit;
        }
      } else if (rstyle == REGION_ONE) {
        for (i = 0; i < ntri; i++) {
          flag = 0;
          if (region->match(pts[tris[i].p1].x)) flag = 1;
          if (region->match(pts[tris[i].p2].x)) flag = 1;
          if (region->match(pts[tris[i].p3].x)) flag = 1;
          if (flag) tris[i].mask |= bit;
        }
      } else if (rstyle == REGION_CENTER) {
        for (i = 0; i < ntri; i++) {
          x[0] = (pts[tris[i].p1].x[0] + pts[tris[i].p2].x[0] + 
                  pts[tris[i].p3].x[0]) / 3.0;
          x[1] = (pts[tris[i].p1].x[1] + pts[tris[i].p2].x[1] + 
                  pts[tris[i].p3].x[1]) / 3.0;
          x[2] = (pts[tris[i].p1].x[2] + pts[tris[i].p2].x[2] + 
                  pts[tris[i].p3].x[2]) / 3.0;
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

  } else if (strcmp(arg[1],"union") == 0) {
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

  } else if (strcmp(arg[1],"intersect") == 0) {
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

void Surf::collate_vector_irregular(int nrow, int *l2g, 
                                    double *in, int instride, double *out)
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

void Surf::collate_array_irregular(int nrow, int ncol, int *l2g, 
                                   double **in, double **out)
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
    fwrite(&bitmask[i],sizeof(int),1,fp);
  }

  fwrite(&npoint,sizeof(int),1,fp);
  fwrite(pts,sizeof(Point),npoint,fp);

  if (domain->dimension == 2) {
    fwrite(&nline,sizeof(int),1,fp);
    for (int i = 0; i < nline; i++) {
      fwrite(&lines[i].type,sizeof(int),1,fp);
      fwrite(&lines[i].mask,sizeof(int),1,fp);
      fwrite(&lines[i].p1,sizeof(int),2,fp);
    }
  }
  if (domain->dimension == 3) {
    fwrite(&ntri,sizeof(int),1,fp);
    for (int i = 0; i < ntri; i++) {
      fwrite(&tris[i].type,sizeof(int),1,fp);
      fwrite(&tris[i].mask,sizeof(int),1,fp);
      fwrite(&tris[i].p1,sizeof(int),3,fp);
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
    if (me == 0) fread(&bitmask[i],sizeof(int),1,fp);
    MPI_Bcast(&bitmask[i],1,MPI_INT,0,world);
  }

  if (me == 0) fread(&npoint,sizeof(int),1,fp);
  MPI_Bcast(&npoint,1,MPI_INT,0,world);
  pts = (Point *) memory->smalloc(npoint*sizeof(Point),"surf:pts");
  if (me == 0) fread(pts,sizeof(Point),npoint,fp);
  MPI_Bcast(pts,npoint*sizeof(Point),MPI_CHAR,0,world);

  if (domain->dimension == 2) {
    if (me == 0) fread(&nline,sizeof(int),1,fp);
    MPI_Bcast(&nline,1,MPI_INT,0,world);
    lines = (Line *) memory->smalloc(nline*sizeof(Line),"surf:lines");

    if (me == 0) {
      for (int i = 0; i < nline; i++) {
        fread(&lines[i].type,sizeof(int),1,fp);
        fread(&lines[i].mask,sizeof(int),1,fp);
        lines[i].isc = lines[i].isr = -1;
        fread(&lines[i].p1,sizeof(int),2,fp);
        lines[i].norm[0] = lines[i].norm[2] = lines[i].norm[2] = 0.0;
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
        fread(&tris[i].type,sizeof(int),1,fp);
        fread(&tris[i].mask,sizeof(int),1,fp);
        tris[i].isc = tris[i].isr = -1;
        fread(&tris[i].p1,sizeof(int),3,fp);
        tris[i].norm[0] = tris[i].norm[2] = tris[i].norm[2] = 0.0;
      }
    }
    MPI_Bcast(tris,ntri*sizeof(Tri),MPI_CHAR,0,world);
  }
}

/* ---------------------------------------------------------------------- */

bigint Surf::memory_usage()
{
  bigint bytes = 0;
  bytes += (bigint) npoint * sizeof(Point);
  bytes += (bigint) nline * sizeof(Line);
  bytes += (bigint) ntri * sizeof(Tri);
  bytes += nlocal * sizeof(int);
  return bytes;
}


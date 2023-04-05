/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "domain.h"
#include "style_region.h"
#include "update.h"
#include "grid.h"
#include "comm.h"
#include "particle.h"
#include "region.h"
#include "surf.h"
#include "surf_collide.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};         // several files
enum{PERIODIC,OUTFLOW,REFLECT,SURFACE,AXISYM};  // several files

#define DELTAREGION 4

/* ---------------------------------------------------------------------- */

Domain::Domain(SPARTA *sparta) : Pointers(sparta)
{
  box_exist = 0;
  dimension = 3;
  axisymmetric = 0;
  boundary_collision_check = 1;

  for (int i = 0; i < 6; i++) bflag[i] = PERIODIC;
  for (int i = 0; i < 6; i++) surf_collide[i] = surf_react[i] = -1;

  // surface normals of 6 box faces pointed inward towards particles

  norm[XLO][0] =  1.0; norm[XLO][1] =  0.0; norm[XLO][2] =  0.0;
  norm[XHI][0] = -1.0; norm[XHI][1] =  0.0; norm[XHI][2] =  0.0;
  norm[YLO][0] =  0.0; norm[YLO][1] =  1.0; norm[YLO][2] =  0.0;
  norm[YHI][0] =  0.0; norm[YHI][1] = -1.0; norm[YHI][2] =  0.0;
  norm[ZLO][0] =  0.0; norm[ZLO][1] =  0.0; norm[ZLO][2] =  1.0;
  norm[ZHI][0] =  0.0; norm[ZHI][1] =  0.0; norm[ZHI][2] = -1.0;

  nregion = maxregion = 0;
  regions = NULL;
  copy = copymode = 0;
}

/* ---------------------------------------------------------------------- */

Domain::~Domain()
{
  if (copy || copymode) return;

  for (int i = 0; i < nregion; i++) delete regions[i];
  memory->sfree(regions);
}

/* ---------------------------------------------------------------------- */

void Domain::init()
{
  if (axisymmetric && dimension != 2)
    error->all(FLERR,"Axi-symmetry only allowed for 2d simulation");
  if (dimension == 2 && (bflag[ZLO] != PERIODIC || bflag[ZHI] != PERIODIC))
    error->all(FLERR,"Z dimension must be periodic for 2d simulation");

  // check grid cutoff versus box size

  int cutflag = 0;
  if (bflag[0] == PERIODIC && grid->cutoff > xprd) cutflag = 1;
  if (bflag[2] == PERIODIC && grid->cutoff > yprd) cutflag = 1;
  if (dimension == 3 && bflag[4] == PERIODIC && grid->cutoff > zprd)
    cutflag = 1;
  if (cutflag) error->all(FLERR,"Grid cutoff is longer than "
                          "box length in a periodic dimension");

  // check that every SURFACE boundary is assigned to a surf collision model
  // skip if caller turned off the check, e.g. BalanceGrid

  if (boundary_collision_check) {
    for (int i = 0; i < 6; i++)
      if (bflag[i] == SURFACE && surf_collide[i] < 0)
        error->all(FLERR,"Box boundary not assigned a surf_collide ID");
  }

  // if a SURFACE boundary is assigned a reaction model
  // then its collision model must allow reactions

  for (int i = 0; i < 2*dimension; i++)
    if (surf_react[i] >= 0 && surf->sc[surf_collide[i]]->allowreact == 0)
      error->all(FLERR,"Box face with reaction model, "
                 "but collision model does not allow reactions");

  // surfreactany = 1 if any face has surface reactions assigned to it

  surfreactany = 0;
  for (int i = 0; i < 6; i++)
    if (surf_react[i] >= 0) surfreactany = 1;
}

/* ----------------------------------------------------------------------
   set initial global box
   assumes boxlo/hi already set
------------------------------------------------------------------------- */

void Domain::set_initial_box()
{
  if (boxlo[0] >= boxhi[0] || boxlo[1] >= boxhi[1] || boxlo[2] >= boxhi[2])
    error->one(FLERR,"Box bounds are invalid");
}

/* ----------------------------------------------------------------------
   set global box params
   assumes boxlo/hi already set
------------------------------------------------------------------------- */

void Domain::set_global_box()
{
  prd[0] = xprd = boxhi[0] - boxlo[0];
  prd[1] = yprd = boxhi[1] - boxlo[1];
  prd[2] = zprd = boxhi[2] - boxlo[2];
}

/* ----------------------------------------------------------------------
   boundary settings from input script
------------------------------------------------------------------------- */

void Domain::set_boundary(int narg, char **arg)
{
  if (domain->box_exist)
    error->all(FLERR,"Boundary command after simulation box is defined");

  if (narg != 3) error->all(FLERR,"Illegal boundary command");

  char c;
  int m = 0;
  for (int idim = 0; idim < 3; idim++)
    for (int iside = 0; iside < 2; iside++) {
      if (iside == 0) c = arg[idim][0];
      else if (iside == 1 && strlen(arg[idim]) == 1) c = arg[idim][0];
      else c = arg[idim][1];

      if (c == 'o') bflag[m] = OUTFLOW;
      else if (c == 'p') bflag[m] = PERIODIC;
      else if (c == 'r') bflag[m] = REFLECT;
      else if (c == 's') bflag[m] = SURFACE;
      else if (c == 'a') bflag[m] = AXISYM;
      else error->all(FLERR,"Illegal boundary command");

      surf_collide[m] = surf_react[m] = -1;
      m++;
    }

  if (dimension == 2 && (bflag[ZLO] != PERIODIC || bflag[ZHI] != PERIODIC))
    error->all(FLERR,"Z dimension must be periodic for 2d simulation");

  if (bflag[XLO] == AXISYM || bflag[XHI] == AXISYM ||
      bflag[YHI] == AXISYM || bflag[ZLO] == AXISYM || bflag[ZHI] == AXISYM)
    error->all(FLERR,"Only ylo boundary can be axi-symmetric");

  if (bflag[YLO] == AXISYM) {
    axisymmetric = 1;
    if (bflag[YHI] == PERIODIC)
      error->all(FLERR,"Y cannot be periodic for axi-symmetric");
  }

  for (m = 0; m < 6; m += 2)
    if (bflag[m] == PERIODIC || bflag[m+1] == PERIODIC) {
      if (bflag[m] != PERIODIC || bflag[m+1] != PERIODIC)
        error->all(FLERR,"Both sides of boundary must be periodic");
    }
}

/* ----------------------------------------------------------------------
   check whether simulation box is periodic in each dimension
   return dims[DIM] = 1/0 for periodic or not for DIM = 0,1,2
     for 2d simulation ZDIM is periodic
   return 1 if all dims are periodic, 0 if not
------------------------------------------------------------------------- */

int Domain::periodic(int *dims)
{
  dims[0] = dims[1] = dims[2] = 0;
  if (bflag[0] == PERIODIC) dims[0] = 1;
  if (bflag[2] == PERIODIC) dims[1] = 1;
  if (bflag[4] == PERIODIC) dims[2] = 1;
  if (dims[0] && dims[1] && dims[2]) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   boundary modifications from input script
   apply mods to list of faces
------------------------------------------------------------------------- */

void Domain::boundary_modify(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal bound_modify command");

  int faces[6];
  int nface = 0;

  int iarg = 0;
  while (iarg < 6) {
    if (strcmp(arg[iarg],"xlo") == 0) faces[nface++] = XLO;
    else if (strcmp(arg[iarg],"xhi") == 0) faces[nface++] = XHI;
    else if (strcmp(arg[iarg],"ylo") == 0) faces[nface++] = YLO;
    else if (strcmp(arg[iarg],"yhi") == 0) faces[nface++] = YHI;
    else if (strcmp(arg[iarg],"zlo") == 0) faces[nface++] = ZLO;
    else if (strcmp(arg[iarg],"zhi") == 0) faces[nface++] = ZHI;
    else break;
    iarg++;
  }

  if (nface == 0 || iarg == narg)
    error->all(FLERR,"Illegal bound_modify command");

  if (dimension != 3) {
    for (int i = 0; i < nface; i++)
      if (faces[i] == ZLO || faces[i] == ZHI)
        error->all(FLERR,
                   "Bound_modify cannot alter zlo/zhi faces for 2d system");
  }

  // setting surf_collide[] and surf_react[] indices here vs init()
  // assumes SurfCollide and SurfReact styles are static (never deleted)

  while (iarg < narg) {
    if (strcmp(arg[iarg],"collide") == 0) {
      if (iarg + 2 > narg) error->all(FLERR,"Illegal bound_modify command");
      int index = surf->find_collide(arg[iarg+1]);
      if (index < 0)
        error->all(FLERR,"Bound_modify surf_collide ID is unknown");
      for (int i = 0; i < nface; i++) {
        if (bflag[faces[i]] != SURFACE)
          error->all(FLERR,"Bound_modify surf requires boundary be a surface");
        surf_collide[faces[i]] = index;
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"react") == 0) {
      if (iarg + 2 > narg) error->all(FLERR,"Illegal bound_modify command");
      int index;
      if (strcmp(arg[iarg+1],"none") == 0) index = -1;
      else {
        index = surf->find_react(arg[iarg+1]);
        if (index < 0)
          error->all(FLERR,"Bound_modify surf_react ID is unknown");
      }
      for (int i = 0; i < nface; i++) {
        if (bflag[faces[i]] != SURFACE)
          error->all(FLERR,"Bound_modify react requires boundary be a surface");
        surf_react[faces[i]] = index;
      }
      iarg += 2;
    } else error->all(FLERR,"Illegal bound_modify command");
  }
}

/* ----------------------------------------------------------------------
   particle ip hits global boundary on face in icell
   called by Update::move()
   xnew = final position of particle at end of move
   return boundary type of global boundary
   return reaction = index of reaction (1 to N) that took place, 0 = no reaction
   if needed, update particle x,v,xnew due to collision
------------------------------------------------------------------------- */

int Domain::collide(Particle::OnePart *&ip, int face, int icell, double *xnew,
                    double &dtremain, Particle::OnePart *&jp, int &reaction)
{
  jp = NULL;
  reaction = 0;

  switch (bflag[face]) {

  // outflow boundary, particle deleted by caller

  case OUTFLOW:
    return OUTFLOW;

  // periodic boundary
  // set x to be on periodic box face
  // adjust xnew by periodic box length

  case PERIODIC:
    {
      double *x = ip->x;

      switch (face) {
      case XLO:
        x[0] = boxhi[0];
        xnew[0] += xprd;
        break;
      case XHI:
        x[0] = boxlo[0];
        xnew[0] -= xprd;
        break;
      case YLO:
        x[1] = boxhi[1];
        xnew[1] += yprd;
        break;
      case YHI:
        x[1] = boxlo[1];
        xnew[1] -= yprd;
        break;
      case ZLO:
        x[2] = boxhi[2];
        xnew[2] += zprd;
        break;
      case ZHI:
        x[2] = boxlo[2];
        xnew[2] -= zprd;
        break;
      }

      return PERIODIC;
    }

  // specular reflection boundary
  // adjust xnew and velocity

  case REFLECT:
    {
      double *v = ip->v;
      double *lo = grid->cells[icell].lo;
      double *hi = grid->cells[icell].hi;
      int dim = face / 2;

      if (face % 2 == 0) {
        xnew[dim] = lo[dim] + (lo[dim]-xnew[dim]);
        v[dim] = -v[dim];
      } else {
        xnew[dim] = hi[dim] - (xnew[dim]-hi[dim]);
        v[dim] = -v[dim];
      }

      return REFLECT;
    }

  // treat global boundary as a surface
  // particle velocity is changed by surface collision model
  // dtremain may be changed by collision model
  // reset all components of xnew, in case dtremain changed
  // if axisymmetric, caller will reset again, including xnew[2]
  // pass -face to collide() to distinguish from surf element collision

  case SURFACE:
    {
      jp = surf->sc[surf_collide[face]]->
        collide(ip,dtremain,-(face+1),norm[face],surf_react[face],reaction);

      if (ip) {
        double *x = ip->x;
        double *v = ip->v;
        xnew[0] = x[0] + dtremain*v[0];
        xnew[1] = x[1] + dtremain*v[1];
        if (dimension == 3) xnew[2] = x[2] + dtremain*v[2];
      }

      return SURFACE;
    }

  }

  return 0;
}

/* ----------------------------------------------------------------------
   undo the periodic remapping of particle coords
     performed by a previous call to collide()
   called by Update::move() if unable to move particle to a new cell
     that contains the remapped coords
------------------------------------------------------------------------- */

void Domain::uncollide(int face, double *x)
{
  switch (face) {
  case XLO:
    x[0] = boxlo[0];
    break;
  case XHI:
    x[0] = boxhi[0];
    break;
  case YLO:
    x[1] = boxlo[1];
    break;
  case YHI:
    x[1] = boxhi[1];
    break;
  case ZLO:
    x[2] = boxlo[2];
    break;
  case ZHI:
    x[2] = boxhi[2];
    break;
  }
}

/* ----------------------------------------------------------------------
   create a new region
------------------------------------------------------------------------- */

void Domain::add_region(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal region command");

  if (strcmp(arg[1],"delete") == 0) {
    delete_region(narg,arg);
    return;
  }

  if (find_region(arg[0]) >= 0) error->all(FLERR,"Reuse of region ID");

  // extend Region list if necessary

  if (nregion == maxregion) {
    maxregion += DELTAREGION;
    regions = (Region **)
      memory->srealloc(regions,maxregion*sizeof(Region *),"domain:regions");
  }

  // create the Region

  if (strcmp(arg[1],"none") == 0) error->all(FLERR,"Invalid region style");

#define REGION_CLASS
#define RegionStyle(key,Class) \
  else if (strcmp(arg[1],#key) == 0) \
    regions[nregion] = new Class(sparta,narg,arg);
#include "style_region.h"
#undef REGION_CLASS

  else error->all(FLERR,"Unrecognized region style");

  nregion++;
}

/* ----------------------------------------------------------------------
   delete a region
------------------------------------------------------------------------- */

void Domain::delete_region(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal region command");

  int iregion = find_region(arg[0]);
  if (iregion == -1) error->all(FLERR,"Delete region ID does not exist");

  delete regions[iregion];
  regions[iregion] = regions[nregion-1];
  nregion--;
}

/* ----------------------------------------------------------------------
   return region index if name matches existing region ID
   return -1 if no such region
------------------------------------------------------------------------- */

int Domain::find_region(char *name)
{
  for (int iregion = 0; iregion < nregion; iregion++)
    if (strcmp(name,regions[iregion]->id) == 0) return iregion;
  return -1;
}

/* ----------------------------------------------------------------------
   print box info
------------------------------------------------------------------------- */

void Domain::print_box(const char *str)
{
  if (comm->me == 0) {
    if (screen)
      fprintf(screen,"%sorthogonal box = (%g %g %g) to (%g %g %g)\n",
              str,boxlo[0],boxlo[1],boxlo[2],boxhi[0],boxhi[1],boxhi[2]);
    if (logfile)
      fprintf(logfile,"%sorthogonal box = (%g %g %g) to (%g %g %g)\n",
              str,boxlo[0],boxlo[1],boxlo[2],boxhi[0],boxhi[1],boxhi[2]);
  }
}

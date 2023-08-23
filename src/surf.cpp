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
#include "random_knuth.h"
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
  me = comm->me;
  nprocs = comm->nprocs;

  exist = 0;
  implicit = 0;
  distributed = 0;
  surf_collision_check = 1;
  localghost_changed_step = -1;
  
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

  nown = maxown = 0;
  mylines = NULL;
  mytris = NULL;

  nunique = 0;
  unique = NULL;
  uniqueID = NULL;
  
  nsc = maxsc = 0;
  sc = NULL;
  nsr = maxsr = 0;
  sr = NULL;

  tally_comm = TALLYAUTO;

  // RNG for unique surf assignment

  urandom = NULL;

  // custom per-surf vectors/arrays

  ncustom = 0;
  ename = NULL;
  etype = esize = estatus = ewhich = NULL;
  size_custom = 0;
  size_custom_local = NULL;
  
  ncustom_ivec = ncustom_iarray = 0;
  icustom_ivec = icustom_iarray = NULL;
  eivec = eivec_local = NULL;
  eiarray = eiarray_local = NULL;
  eicol = NULL;

  ncustom_dvec = ncustom_darray = 0;
  icustom_dvec = icustom_darray = NULL;
  edvec = edvec_local = NULL;
  edarray = edarray_local = NULL;
  edcol = NULL;

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

  memory->destroy(unique);
  memory->destroy(uniqueID);
  
  for (int i = 0; i < nsc; i++) delete sc[i];
  memory->sfree(sc);
  for (int i = 0; i < nsr; i++) delete sr[i];
  memory->sfree(sr);

  hash->clear();
  delete hash;

  for (int i = 0; i < ncustom; i++) delete [] ename[i];
  memory->sfree(ename);
  memory->destroy(etype);
  memory->destroy(esize);
  memory->destroy(estatus);
  memory->destroy(ewhich);
  memory->destroy(size_custom_local);

  for (int i = 0; i < ncustom_ivec; i++) {
    memory->destroy(eivec[i]);
    memory->destroy(eivec_local[i]);
  }
  for (int i = 0; i < ncustom_iarray; i++) {
    memory->destroy(eiarray[i]);
    memory->destroy(eiarray_local[i]);
  }
  for (int i = 0; i < ncustom_dvec; i++) {
    memory->destroy(edvec[i]);
    memory->destroy(edvec_local[i]);
  }
  for (int i = 0; i < ncustom_darray; i++) {
    memory->destroy(edarray[i]);
    memory->destroy(edarray_local[i]);
  }

  memory->destroy(icustom_ivec);
  memory->destroy(icustom_iarray);
  memory->sfree(eivec);
  memory->sfree(eivec_local);
  memory->sfree(eiarray);
  memory->sfree(eiarray_local);
  memory->destroy(eicol);
  
  memory->destroy(icustom_dvec);
  memory->destroy(icustom_darray);
  memory->sfree(edvec);
  memory->sfree(edvec_local);
  memory->sfree(edarray);
  memory->sfree(edarray_local);
  memory->destroy(edcol);

  delete urandom;
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

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"collide") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal surf_modify command");
      if (!exist) error->all(FLERR,"Surf_modify when surfs do not yet exist");

      int isc = find_collide(arg[iarg+1]);
      if (isc < 0) error->all(FLERR,"Could not find surf_modify sc-ID");

      // NOTE: is this also needed for mylines and mytris?
      // set surf collision model for each surf in surface group

      if (domain->dimension == 2) {
        for (int i = 0; i < nlocal+nghost; i++)
          if (lines[i].mask & groupbit) lines[i].isc = isc;
      } else {
        for (int i = 0; i < nlocal+nghost; i++)
          if (tris[i].mask & groupbit) tris[i].isc = isc;
      }

      iarg += 2;

    } else if (strcmp(arg[iarg],"react") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal surf_modify command");
      if (!exist) error->all(FLERR,"Surf_modify when surfs do not yet exist");

      int isr;
      if (strcmp(arg[iarg+1],"none") == 0) isr = -1;
      else {
        isr = find_react(arg[iarg+1]);
        if (isr < 0) error->all(FLERR,"Could not find surf_modify sr-ID");
      }

      // set surf reaction model for each surf in surface group

      if (domain->dimension == 2) {
        for (int i = 0; i < nlocal+nghost; i++)
          if (lines[i].mask & groupbit) lines[i].isr = isr;
      } else {
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

  // check that surf element types are all values >= 1

  bigint flag,allflag;

  flag = 0;
  if (domain->dimension == 2) {
    for (int i = 0; i < nlocal; i++)
      if (lines[i].type <= 0) flag++;
  } else {
    for (int i = 0; i < nlocal; i++)
      if (tris[i].type <= 0) flag++;
  }

  if (distributed)
    MPI_Allreduce(&flag,&allflag,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  else allflag = flag;

  if (allflag) {
    char str[64];
    sprintf(str,BIGINT_FORMAT
            " surface elements with invalid type <= 0",allflag);
    error->all(FLERR,str);
  }

  // check that every element is assigned to a surf collision model
  // skip if caller turned off the check, e.g. BalanceGrid, b/c too early

  if (surf_collision_check) {
    flag = 0;
    if (domain->dimension == 2) {
      for (int i = 0; i < nlocal+nghost; i++)
        if (lines[i].isc < 0) flag++;
    } else {
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
    if (domain->dimension == 2) {
      for (int i = 0; i < nlocal+nghost; i++)
        if (lines[i].isr >= 0 && sc[lines[i].isc]->allowreact == 0) flag++;
    } else {
      for (int i = 0; i < nlocal+nghost; i++)
        if (tris[i].isr >= 0 && sc[tris[i].isc]->allowreact == 0) flag++;
    }

    if (distributed)
      MPI_Allreduce(&flag,&allflag,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
    else allflag = flag;

    if (allflag) {
      char str[128];
      sprintf(str,BIGINT_FORMAT " surface elements with reaction model, "
              "but invalid collision model",allflag);
      error->all(FLERR,str);
    }
  }

  // checks on transparent surfaces
  // must be assigned to transparent surf collision model

  if (surf_collision_check) {
    flag = 0;
    if (domain->dimension == 2) {
      for (int i = 0; i < nlocal+nghost; i++) {
        if (!lines[i].transparent) continue;
        if (!sc[lines[i].isc]->transparent) flag++;
      }
    } else {
      for (int i = 0; i < nlocal+nghost; i++) {
        if (!tris[i].transparent) continue;
        if (!sc[tris[i].isc]->transparent) flag++;
      }
    }

    if (distributed)
      MPI_Allreduce(&flag,&allflag,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
    else allflag = flag;

    if (allflag) {
      char str[128];
      sprintf(str,BIGINT_FORMAT " transparent surface elements "
              "with invalid collision model or reaction model",allflag);
      error->all(FLERR,str);
    }
  }

  // initialize surf collision and reaction models

  for (int i = 0; i < nsc; i++) sc[i]->init();
  for (int i = 0; i < nsr; i++) sr[i]->init();
}

/* ----------------------------------------------------------------------
   remove all explicit surfs and deallocate memory
   called by add_surfs() via ReadSurf and RemoveSurf
------------------------------------------------------------------------- */

void Surf::clear_explicit()
{
  nsurf = 0;
  nlocal = nghost = nmax = 0;
  nown = maxown = 0;
  
  memory->sfree(lines);
  memory->sfree(tris);
  memory->sfree(mylines);
  memory->sfree(mytris);
  
  lines = NULL;
  tris = NULL;
  mylines = NULL;
  mytris = NULL;

  hash->clear();
  hashfilled = 0;
}

/* ----------------------------------------------------------------------
   remove all implicit surfs but do not deallocate memory
   called by FixAblate
------------------------------------------------------------------------- */

void Surf::clear_implicit()
{
  nsurf = 0;
  nlocal = nghost = 0;
  nown = 0;
  hash->clear();
  hashfilled = 0;
}

/* ----------------------------------------------------------------------
   remove ghost surfs
------------------------------------------------------------------------- */

void Surf::remove_ghosts()
{
  nghost = 0;
}

/* ----------------------------------------------------------------------
   add a line to lines list
   called by ReadISurf or FixAblate via Marching Squares
------------------------------------------------------------------------- */

void Surf::add_line(surfint id, int itype, double *p1, double *p2)
{
  if (nlocal == nmax) {
    if ((bigint) nmax + DELTA > MAXSMALLINT)
      error->one(FLERR,"Surf add_line overflowed");
    nmax += DELTA;
    grow(nmax-DELTA);
  }

  lines[nlocal].id = id;
  lines[nlocal].type = itype;
  lines[nlocal].mask = 1;
  lines[nlocal].isc = lines[nlocal].isr = -1;
  lines[nlocal].p1[0] = p1[0];
  lines[nlocal].p1[1] = p1[1];
  lines[nlocal].p1[2] = 0.0;
  lines[nlocal].p2[0] = p2[0];
  lines[nlocal].p2[1] = p2[1];
  lines[nlocal].p2[2] = 0.0;
  lines[nlocal].transparent = 0;
  nlocal++;
}

/* ----------------------------------------------------------------------
   add a line to owned or ghost list, depending on ownflag
   called by Grid::unpack_one() or Grid::coarsen_cell()
------------------------------------------------------------------------- */

void Surf::add_line_copy(int ownflag, Line *line)
{
  int index;

  if (ownflag) {
    if (nlocal == nmax) {
      if ((bigint) nmax + DELTA > MAXSMALLINT)
        error->one(FLERR,"Surf add_line_copy overflowed");
      nmax += DELTA;
      grow(nmax-DELTA);
    }
    index = nlocal;
    nlocal++;

  } else {
    if (nlocal+nghost == nmax) {
      if ((bigint) nmax + DELTA > MAXSMALLINT)
        error->one(FLERR,"Surf add_line_copy overflowed");
      nmax += DELTA;
      grow(nmax-DELTA);
    }
    index = nlocal+nghost;
    nghost++;
  }

  memcpy(&lines[index],line,sizeof(Line));
}

/* ----------------------------------------------------------------------
   add a triangle to tris list
   called by ReadISurf or FixAblate via Marching Cubes
------------------------------------------------------------------------- */

void Surf::add_tri(surfint id, int itype, double *p1, double *p2, double *p3)
{
  if (nlocal == nmax) {
    if ((bigint) nmax + DELTA > MAXSMALLINT)
      error->one(FLERR,"Surf add_tri overflowed");
    nmax += DELTA;
    grow(nmax-DELTA);
  }

  tris[nlocal].id = id;
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
  tris[nlocal].transparent = 0;
  nlocal++;
}

/* ----------------------------------------------------------------------
   add a triangle to owned or ghost list, depending on ownflag
   called by Grid::unpack_one or Grid::coarsen_cell()
------------------------------------------------------------------------- */

void Surf::add_tri_copy(int ownflag, Tri *tri)
{
  int index;

  if (ownflag) {
    if (nlocal == nmax) {
      if ((bigint) nmax + DELTA > MAXSMALLINT)
        error->one(FLERR,"Surf add_tri_copy overflowed");
      nmax += DELTA;
      grow(nmax-DELTA);
    }
    index = nlocal;
    nlocal++;

  } else {
    if (nlocal+nghost == nmax) {
      if ((bigint) nmax + DELTA > MAXSMALLINT)
        error->one(FLERR,"Surf add_tri_copy overflowed");
      nmax += DELTA;
      grow(nmax-DELTA);
    }
    index = nlocal+nghost;
    nghost++;
  }

  memcpy(&tris[index],tri,sizeof(Tri));
}

/* ----------------------------------------------------------------------
   add list of surfs created by caller to Surf data structs
     for explicit surfs only, all or distributed
   replace = 0: add newlines/newtris to existing ones
   replace = 1: wipe out existing surfs, replace with newlines/newtris
   ncount = # of new surfs contributed by this proc
     new surfs are distributed across procs
     can be in any order, but IDs need to be unique and contiguous
   newlines/newtris = list of new surfs, other ptr is NULL
   nc = # of custom attributes (vecs/arrays) for each surf
   index_custom = index for each custom vec or array in Surf custom lists
   cvalues = custom values for each surf in same order as newlines/newtris
     1st value is surf ID, remaining values are for nc vecs/arrays
   called by ReadSurf, RemoveSurf, ReadRestart
------------------------------------------------------------------------- */

void Surf::add_surfs(int replace, int ncount,
		     Line *newlines, Tri *newtris,
		     int nc, int *index_custom, double **cvalues)
{
  // if replace: remove all existing surfs and their memory
  // remove ghost surfs for replace or add
  
  if (replace) clear_explicit();
  remove_ghosts();
  
  // (re)allocate data structs for adding new surfs

  bigint bncount = ncount;
  bigint nsurf_new;
  MPI_Allreduce(&bncount,&nsurf_new,1,MPI_SPARTA_BIGINT,MPI_SUM,world);

  bigint nsurf_old = nsurf;
  nsurf += nsurf_new;

  int nlocal_old = nlocal;

  int nown_old = nown;
  nown = nsurf / nprocs;
  if (me < nsurf % nprocs) nown++;

  if (!distributed) {
    nlocal += nsurf_new;
    while (nmax < nlocal) nmax += DELTA;
    grow(nlocal_old);
  } else {
    maxown = nown;
    grow_own(nown_old);
  }

  // reallocate data structs for adding new custom values
  // wait unto now b/c nown has now been reset

  reallocate_custom();

  // offset IDs of new surfs by pre-existing nsurf_old
  // for both new lines/tris and custom data in cvalues

  if (domain->dimension == 2)
    for (int i = 0; i < ncount; i++)
      newlines[i].id += nsurf_old;
  else
    for (int i = 0; i < ncount; i++)
      newtris[i].id += nsurf_old;
 
  if (nc) 
    for (int i = 0; i < ncount; i++)
      cvalues[i][0] = ubuf(((surfint) ubuf(cvalues[i][0]).i) + nsurf_old).d;

  // redistribute surfs to correct layout in Surf data structs
  // also checks that new surfs have IDs contiguous from 1 to N
  
  redistribute_surfs(ncount,newlines,newtris,nc,index_custom,cvalues,
		     nsurf_new,nsurf_old);

  // check if new surf IDs are contiguous from 1 to Nsurf_new
  // if any ID = 0 in rendezvous output, new surf IDs were NOT contiguous

  int surfperproc = nsurf_new / nprocs;
  int ncontig = surfperproc;
  if (me == nprocs-1) ncontig = nsurf_new - (nprocs-1)*surfperproc;

  int flag = 0;
  if (domain->dimension == 2) {
    if (!distributed) {
      for (int i = nlocal_old; i < nlocal; i++)
	if (lines[i].id == 0) flag++;
    } else {
      for (int i = nown_old; i < nown; i++)
	if (surf->mylines[i].id == 0) flag++;
    }
  } else {
    if (!distributed) {
      for (int i = nlocal_old; i < nlocal; i++)
	if (tris[i].id == 0) flag++;
    } else {
      for (int i = nown_old; i < nown; i++)
	if (surf->mytris[i].id == 0) flag++;
    }
  }

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) {
    char str[128];
    sprintf(str,"Missing read_surf IDs = %d",flagall);
    error->all(FLERR,str);
  }
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
  hashfilled = 1;

  if (domain->dimension == 2) {
    for (int isurf = 0; isurf < nlocal; isurf++)
      (*hash)[lines[isurf].id] = isurf;
  } else {
    for (int isurf = 0; isurf < nlocal; isurf++)
      (*hash)[tris[isurf].id] = isurf;
  }
}

/* ----------------------------------------------------------------------
   return 1 if all surfs are transparent, else return 0
   called by set_inout()
------------------------------------------------------------------------- */

int Surf::all_transparent()
{
  // implicit surfs cannot be transparent

  if (implicit) return 0;

  // explicit surfs may be transparent

  int flag = 0;
  if (domain->dimension == 2) {
    for (int i = 0; i < nlocal; i++)
      if (!lines[i].transparent) flag = 1;
  } else {
    for (int i = 0; i < nlocal; i++)
      if (!tris[i].transparent) flag = 1;
  }

  int allflag;
  if (distributed)
    MPI_Allreduce(&flag,&allflag,1,MPI_INT,MPI_SUM,world);
  else allflag = flag;

  if (allflag) return 0;
  return 1;
}

/* ----------------------------------------------------------------------
   count owned surf elements = every Pth surf from global Nsurf list
   only called when surfs = explict (all or distributed)
   nothing to do when distributed b/c mylines/mytris already setup
------------------------------------------------------------------------ */

void Surf::setup_owned()
{
  if (distributed) return;

  nown = nsurf/nprocs;
  if (comm->me < nsurf % nprocs) nown++;
}

/* ----------------------------------------------------------------------
   set bounding box around all surfs based on their pts
   sets surf->bblo and surf->bbhi
   for 2d, set zlo,zhi to box bounds
   only called when surfs = explict (all or distributed)
------------------------------------------------------------------------- */

void Surf::bbox_all()
{
  int i,j;
  double bblo_one[3],bbhi_one[3];
  double *x;

  int istart,istop,idelta;
  Line *linelist;
  Tri *trilist;

  if (!distributed) {
    istart = me;
    istop = nlocal;
    idelta = nprocs;
    linelist = lines;
    trilist = tris;
  } else {
    istart = 0;
    istop = nown;
    idelta = 1;
    linelist = mylines;
    trilist = mytris;
  }

  for (j = 0; j < 3; j++) {
    bblo_one[j] = BIG;
    bbhi_one[j] = -BIG;
  }

  if (domain->dimension == 2) {
    for (i = istart; i < istop; i += idelta) {
      x = linelist[i].p1;
      for (j = 0; j < 2; j++) {
        bblo_one[j] = MIN(bblo_one[j],x[j]);
        bbhi_one[j] = MAX(bbhi_one[j],x[j]);
      }
      x = linelist[i].p2;
      for (j = 0; j < 2; j++) {
        bblo_one[j] = MIN(bblo_one[j],x[j]);
        bbhi_one[j] = MAX(bbhi_one[j],x[j]);
      }
    }
    bblo_one[2] = domain->boxlo[2];
    bbhi_one[2] = domain->boxhi[2];

  } else {
    for (i = istart; i < istop; i += idelta) {
      x = trilist[i].p1;
      for (j = 0; j < 3; j++) {
        bblo_one[j] = MIN(bblo_one[j],x[j]);
        bbhi_one[j] = MAX(bbhi_one[j],x[j]);
      }
      x = trilist[i].p2;
      for (j = 0; j < 3; j++) {
        bblo_one[j] = MIN(bblo_one[j],x[j]);
        bbhi_one[j] = MAX(bbhi_one[j],x[j]);
      }
      x = trilist[i].p3;
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
   set bounding box around one surf based on their pts
   caller passes in the Line or Tri, can be from lines/tris or mylines/mytris
   returns lo,hi which are allocated by caller
   for 2d, set zlo,zhi to box bounds
   only called when surfs = explict (all or distributed)
------------------------------------------------------------------------- */

void Surf::bbox_one(void *ptr, double *lo, double *hi)
{
  double *p1,*p2,*p3;

  if (domain->dimension == 2) {
    Line *line = (Line *) ptr;
    p1 = line->p1; p2 = line->p2;
    lo[0] = MIN(p1[0],p2[0]);
    lo[1] = MIN(p1[1],p2[1]);
    lo[2] = 0.0;
    hi[0] = MAX(p1[0],p2[0]);
    hi[1] = MAX(p1[1],p2[1]);
    hi[2] = 0.0;

  } else {
    Tri *tri = (Tri *) ptr;
    p1 = tri->p1; p2 = tri->p2; p3 = tri->p3;
    lo[0] = MIN(p1[0],p2[0]);
    lo[0] = MIN(lo[0],p3[0]);
    lo[1] = MIN(p1[1],p2[1]);
    lo[1] = MIN(lo[1],p3[1]);
    lo[2] = MIN(p1[2],p2[2]);
    lo[2] = MIN(lo[2],p3[2]);
    hi[0] = MAX(p1[0],p2[0]);
    hi[0] = MAX(hi[0],p3[0]);
    hi[1] = MAX(p1[1],p2[1]);
    hi[1] = MAX(hi[1],p3[1]);
    hi[2] = MAX(p1[2],p2[2]);
    hi[2] = MAX(hi[2],p3[2]);
  }
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
   return vector of masks for my owned surfs
   for distributed or non-distributed explicit surfs
   called by other classes which need to select owned surfs in surf group
------------------------------------------------------------------------- */

void Surf::extract_masks(int *masks)
{
  if (distributed) {
    if (domain->dimension == 2) {
      for (int i = 0; i < nown; i++)
        masks[i] = mylines[i].mask;
    } else {
      for (int i = 0; i < nown; i++)
	masks[i] = mytris[i].mask;
    }
  } else {
    int m = 0;
    if (domain->dimension == 2) {
      for (int i = me; i < nsurf; i += nprocs)
        masks[m++] = lines[i].mask;
    } else {
      for (int i = me; i < nsurf; i += nprocs)
        masks[m++] = tris[i].mask;
    }
  }
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
   check if 2d surf elements are watertight
   each end point should appear exactly once as different ends of 2 lines
   exception: not required of end points on simulation box surface
------------------------------------------------------------------------- */

void Surf::check_watertight_2d()
{
  if (distributed) check_watertight_2d_distributed();
  else check_watertight_2d_all();
}

/* ----------------------------------------------------------------------
   check if 2d surf elements are watertight
   this is for explicit non-distributed surfs where each proc has copy of all
   each proc tests the entire surface, no communication needed
------------------------------------------------------------------------- */

void Surf::check_watertight_2d_all()
{
  // hash end points of all lines
  // key = end point
  // value = 1 if first point, 2 if second point, 3 if both points

  MyHashPoint phash;
  MyPointIt it;

  // insert each end point into hash
  // should appear once at each end
  // error if any duplicate points

  double *p1,*p2;
  OnePoint2d key;
  int value;

  int ndup = 0;
  for (int i = 0; i < nsurf; i++) {
    if (lines[i].transparent) continue;
    p1 = lines[i].p1;
    key.pt[0] = p1[0]; key.pt[1] = p1[1];
    if (phash.find(key) == phash.end()) phash[key] = 1;
    else {
      value = phash[key];
      if (value == 2) phash[key] = 3;
      else ndup++;
    }

    p2 = lines[i].p2;
    key.pt[0] = p2[0]; key.pt[1] = p2[1];
    if (phash.find(key) == phash.end()) phash[key] = 2;
    else {
      value = phash[key];
      if (value == 1) phash[key] = 3;
      else ndup++;
    }
  }

  if (ndup) {
    char str[128];
    sprintf(str,"Watertight check failed with %d duplicate points",ndup);
    error->all(FLERR,str);
  }

  // check that each end point has a match (value = 3)
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
   check if 2d surf elements are watertight
   this is for explicit distributed surfs
   rendezvous communication used to check that each point appears twice
------------------------------------------------------------------------- */

void Surf::check_watertight_2d_distributed()
{
  int n;
  Line *lines_rvous;

  if (implicit) {
    n = nlocal;
    lines_rvous = lines;
  } else {
    n = nown;
    lines_rvous = mylines;
  }

  // allocate memory for rvous input

  int *proclist;
  memory->create(proclist,n*2,"surf:proclist");
  InRvousPoint *inpoint =
    (InRvousPoint *) memory->smalloc((bigint) n*2*sizeof(InRvousPoint),
                                     "surf:inpoint");

  // create rvous inputs
  // proclist = owner of each point
  // each line end point is sent with flag indicating first/second
  // hash of point coord (xy) determines which proc to send to

  int nrvous = 0;
  for (int i = 0; i < n; i++) {
    proclist[nrvous] = hashlittle(lines_rvous[i].p1,2*sizeof(double),0) % nprocs;
    inpoint[nrvous].x[0] = lines_rvous[i].p1[0];
    inpoint[nrvous].x[1] = lines_rvous[i].p1[1];
    inpoint[nrvous].which = 1;
    nrvous++;
    proclist[nrvous] = hashlittle(lines_rvous[i].p2,2*sizeof(double),0) % nprocs;
    inpoint[nrvous].x[0] = lines_rvous[i].p2[0];
    inpoint[nrvous].x[1] = lines_rvous[i].p2[1];
    inpoint[nrvous].which = 2;
    nrvous++;
  }

  // perform rendezvous operation
  // each proc assigned subset of points
  // receives all copies of points, checks if count of each point is valid

  char *buf;
  int nout = comm->rendezvous(1,nrvous,(char *) inpoint,sizeof(InRvousPoint),
                              0,proclist,rendezvous_watertight_2d,
                              0,buf,0,(void *) this);

  memory->destroy(proclist);
  memory->destroy(inpoint);
}

/* ----------------------------------------------------------------------
   callback from rendezvous operation
   process points assigned to me
   inbuf = list of N Inbuf datums
   no outbuf
------------------------------------------------------------------------- */

int Surf::rendezvous_watertight_2d(int n, char *inbuf, int &flag, int *&proclist,
                                   char *&outbuf, void *ptr)
{
  Surf *sptr = (Surf *) ptr;
  Domain *domain = sptr->domain;
  Error *error = sptr->error;
  MPI_Comm world = sptr->world;

  Surf::InRvousPoint *inpoint = (Surf::InRvousPoint *) inbuf;

  // hash all received end points
  // key = end point
  // value = 1 if first point, 2 if second point, 3 if both points

  Surf::MyHashPoint phash;
  Surf::MyPointIt it;

  // insert each point into hash
  // should appear once at each end of a line
  // error if any duplicate points

  Surf::OnePoint2d key;
  int which,value;

  int ndup = 0;
  for (int i = 0; i < n; i++) {
    key.pt[0] = inpoint[i].x[0]; key.pt[1] = inpoint[i].x[1];
    which = inpoint[i].which;
    if (phash.find(key) == phash.end()) phash[key] = which;
    else {
      value = phash[key];
      if (value == 3) ndup++;    // point already seen twice
      else if (value != which) phash[key] = 3;   // this is other point
      else ndup++;               // value = which, this is duplicate point
    }
  }

  int alldup;
  MPI_Allreduce(&ndup,&alldup,1,MPI_INT,MPI_SUM,world);
  if (alldup) {
    char str[128];
    sprintf(str,"Watertight check failed with %d duplicate points",alldup);
    error->all(FLERR,str);
  }

  // check that each end point has a match (value = 3)
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

  int allbad;
  MPI_Allreduce(&nbad,&allbad,1,MPI_INT,MPI_SUM,world);
  if (nbad) {
    char str[128];
    sprintf(str,"Watertight check failed with %d unmatched points",allbad);
    error->all(FLERR,str);
  }

  // flag = 0: no second comm needed in rendezvous

  flag = 0;
  return 0;
}

/* ----------------------------------------------------------------------
   check if 3d surf elements are watertight
   each edge should appear exactly once in each direction
   exception: not required of triangle edge on simulation box surface
------------------------------------------------------------------------- */

void Surf::check_watertight_3d()
{
  if (distributed) check_watertight_3d_distributed();
  else check_watertight_3d_all();
}

/* ----------------------------------------------------------------------
   check if 3d surf elements are watertight
   this is for explicit non-distributed surfs where each proc has copy of all
   each proc tests the entire surface, no communication needed
------------------------------------------------------------------------- */

void Surf::check_watertight_3d_all()
{
  // hash directed edges of all triangles
  // key = directed edge
  // value = 1 if appears once, 2 if reverse also appears once

  MyHash2Point phash;
  My2PointIt it;

  // insert each edge into hash
  // should appear once in each direction
  // error if any duplicate edges

  double *p1,*p2,*p3;
  TwoPoint3d key,keyinv;
  int value;

  int ndup = 0;
  for (int i = 0; i < nsurf; i++) {
    if (tris[i].transparent) continue;
    p1 = tris[i].p1;
    p2 = tris[i].p2;
    p3 = tris[i].p3;

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
  }

  if (ndup) {
    char str[128];
    sprintf(str,"Watertight check failed with %d duplicate edges",ndup);
    error->all(FLERR,str);
  }

  // check that each edge has an inverted match
  // allow for exception if edge is on box face

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  double *pts;

  int nbad = 0;
  for (it = phash.begin(); it != phash.end(); ++it) {
    if (it->second != 2) {
      pts = (double *) it->first.pts;
      if (Geometry::edge_on_hex_face(&pts[0],&pts[3],boxlo,boxhi) < 0) nbad++;
    }
  }

  if (nbad) {
    char str[128];
    sprintf(str,"Watertight check failed with %d unmatched edges",nbad);
    error->all(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   check if 3d surf elements are watertight
   this is for explicit distributed surfs
   rendezvous communication used to check that each edge appears twice
------------------------------------------------------------------------- */

void Surf::check_watertight_3d_distributed()
{
  int n;
  Tri *tris_rvous;

  if (implicit) {
    n = nlocal;
    tris_rvous = tris;
  } else {
    n = nown;
    tris_rvous = mytris;
  }

  // allocate memory for rvous input

  int *proclist;
  memory->create(proclist,n*6,"surf:proclist");
  InRvousEdge *inedge =
    (InRvousEdge *) memory->smalloc((bigint) n*6*sizeof(InRvousEdge),
                                     "surf:inedge");

  // create rvous inputs
  // proclist = owner of each point
  // each triangle edge is sent twice with flag indicating
  //   forward or reverse order
  // hash of edge coords (xyz for 2 pts) determines which proc to send to

  double edge[6];
  double *p1,*p2,*p3;

  int nbytes = 3*sizeof(double);

  int nrvous = 0;
  for (int i = 0; i < n; i++) {
    p1 = tris_rvous[i].p1;
    p2 = tris_rvous[i].p2;
    p3 = tris_rvous[i].p3;

    memcpy(&edge[0],p1,nbytes);
    memcpy(&edge[3],p2,nbytes);
    proclist[nrvous] = hashlittle(edge,2*nbytes,0) % nprocs;
    memcpy(inedge[nrvous].x1,p1,nbytes);
    memcpy(inedge[nrvous].x2,p2,nbytes);
    inedge[nrvous].which = 1;
    nrvous++;

    memcpy(&edge[0],p2,nbytes);
    memcpy(&edge[3],p1,nbytes);
    proclist[nrvous] = hashlittle(edge,2*nbytes,0) % nprocs;
    memcpy(inedge[nrvous].x1,p2,nbytes);
    memcpy(inedge[nrvous].x2,p1,nbytes);
    inedge[nrvous].which = 2;
    nrvous++;

    memcpy(&edge[0],p2,nbytes);
    memcpy(&edge[3],p3,nbytes);
    proclist[nrvous] = hashlittle(edge,2*nbytes,0) % nprocs;
    memcpy(inedge[nrvous].x1,p2,nbytes);
    memcpy(inedge[nrvous].x2,p3,nbytes);
    inedge[nrvous].which = 1;
    nrvous++;

    memcpy(&edge[0],p3,nbytes);
    memcpy(&edge[3],p2,nbytes);
    proclist[nrvous] = hashlittle(edge,2*nbytes,0) % nprocs;
    memcpy(inedge[nrvous].x1,p3,nbytes);
    memcpy(inedge[nrvous].x2,p2,nbytes);
    inedge[nrvous].which = 2;
    nrvous++;

    memcpy(&edge[0],p3,nbytes);
    memcpy(&edge[3],p1,nbytes);
    proclist[nrvous] = hashlittle(edge,2*nbytes,0) % nprocs;
    memcpy(inedge[nrvous].x1,p3,nbytes);
    memcpy(inedge[nrvous].x2,p1,nbytes);
    inedge[nrvous].which = 1;
    nrvous++;

    memcpy(&edge[0],p1,nbytes);
    memcpy(&edge[3],p3,nbytes);
    proclist[nrvous] = hashlittle(edge,2*nbytes,0) % nprocs;
    memcpy(inedge[nrvous].x1,p1,nbytes);
    memcpy(inedge[nrvous].x2,p3,nbytes);
    inedge[nrvous].which = 2;
    nrvous++;
  }

  // perform rendezvous operation
  // each proc assigned subset of edges
  // receives all copies of edges, checks if count of each edge is valid

  char *buf;
  int nout = comm->rendezvous(1,nrvous,(char *) inedge,sizeof(InRvousEdge),
                              0,proclist,rendezvous_watertight_3d,
                              0,buf,0,(void *) this);

  memory->destroy(proclist);
  memory->destroy(inedge);
}

/* ----------------------------------------------------------------------
   callback from rendezvous operation
   process points assigned to me
   inbuf = list of N Inbuf datums
   no outbuf
------------------------------------------------------------------------- */

int Surf::rendezvous_watertight_3d(int n, char *inbuf, int &flag, int *&proclist,
                                   char *&outbuf, void *ptr)
{
  Surf *sptr = (Surf *) ptr;
  Domain *domain = sptr->domain;
  Error *error = sptr->error;
  MPI_Comm world = sptr->world;

  Surf::InRvousEdge *inedge = (Surf::InRvousEdge *) inbuf;

  // hash all received end points
  // key = end point
  // value = 1 if first point, 2 if second point, 3 if both points

  Surf::MyHash2Point phash;
  Surf::My2PointIt it;

  // insert each edge into hash
  // should appear once in each direction
  // error if any duplicate edges

  Surf::TwoPoint3d key;
  double *x1,*x2;
  int which,value;

  int ndup = 0;
  for (int i = 0; i < n; i++) {
    x1 = inedge[i].x1; x2 = inedge[i].x2;
    key.pts[0] = x1[0]; key.pts[1] = x1[1]; key.pts[2] = x1[2];
    key.pts[3] = x2[0]; key.pts[4] = x2[1]; key.pts[5] = x2[2];
    which = inedge[i].which;
    if (phash.find(key) == phash.end()) phash[key] = which;
    else {
      value = phash[key];
      if (value == 3) ndup++;    // edge already seen twice
      else if (value != which) phash[key] = 3;   // this is flipped edge
      else ndup++;               // value = which, this is duplicate edge
    }
  }

  int alldup;
  MPI_Allreduce(&ndup,&alldup,1,MPI_INT,MPI_SUM,world);
  alldup /= 2;              // avoid double counting
  if (alldup) {
    char str[128];
    sprintf(str,"Watertight check failed with %d duplicate edges",alldup);
    error->all(FLERR,str);
  }

  // check that each edge has an inverted match(value = 3)
  // allow for exception if edge is on box face

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  double *pts;

  int nbad = 0;
  for (it = phash.begin(); it != phash.end(); ++it) {
    if (it->second != 3) {
      pts = (double *) it->first.pts;
      if (Geometry::edge_on_hex_face(&pts[0],&pts[3],boxlo,boxhi) < 0) nbad++;
    }
  }

  int allbad;
  MPI_Allreduce(&nbad,&allbad,1,MPI_INT,MPI_SUM,world);
  allbad /= 2;              // avoid double counting
  if (nbad) {
    char str[128];
    sprintf(str,"Watertight check failed with %d unmatched edges",allbad);
    error->all(FLERR,str);
  }

  // flag = 0: no second comm needed in rendezvous

  flag = 0;
  return 0;
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

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  if (domain->dimension == 2) {
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

  } else {
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
      // skip transparent surf elements
      if (line->transparent) continue;
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
      // skip transparent surf elements
      if (tri->transparent) continue;
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

  if (domain->dimension == 2) {
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
  if (domain->dimension == 2) minlen = shortest_line(old);
  else smallest_tri(old,minlen,minarea);

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"  %g %g xlo xhi\n",extentall[0][0],extentall[0][1]);
      fprintf(screen,"  %g %g ylo yhi\n",extentall[1][0],extentall[1][1]);
      fprintf(screen,"  %g %g zlo zhi\n",extentall[2][0],extentall[2][1]);
      if (domain->dimension == 2)
        fprintf(screen,"  %g min line length\n",minlen);
      else {
        fprintf(screen,"  %g min triangle edge length\n",minlen);
        fprintf(screen,"  %g min triangle area\n",minarea);
      }
    }
    if (logfile) {
      fprintf(logfile,"  %g %g xlo xhi\n",extentall[0][0],extentall[0][1]);
      fprintf(logfile,"  %g %g ylo yhi\n",extentall[1][0],extentall[1][1]);
      fprintf(logfile,"  %g %g zlo zhi\n",extentall[2][0],extentall[2][1]);
      if (domain->dimension == 2)
        fprintf(logfile,"  %g min line length\n",minlen);
      else {
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

  if (sparta->suffix_enable) {
    if (sparta->suffix) {
      char estyle[256];
      sprintf(estyle,"%s/%s",arg[1],sparta->suffix);

      if (0) return;

#define SURF_REACT_CLASS
#define SurfReactStyle(key,Class) \
      else if (strcmp(estyle,#key) == 0) { \
        sr[nsr] = new Class(sparta,narg,arg); \
        nsr++; \
        return; \
      }
#include "style_surf_react.h"
#undef SurfReactStyle
#undef SURF_REACT_CLASS
    }
  }

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

  int igroup = find_group(arg[0]);
  if (igroup < 0) igroup = add_group(arg[0]);
  int bit = bitmask[igroup];

  int dim = domain->dimension;

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

// ----------------------------------------------------------------------
// methods for write/read restart info
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   proc 0 writes surf groups to restart file
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
}

/* ----------------------------------------------------------------------
   proc 0 reads surf groups from restart file
   bcast to other procs
------------------------------------------------------------------------- */

void Surf::read_restart(FILE *fp)
{
  int tmp;

  int me = comm->me;

  // if any groups exist, clear them before reading new ones

  for (int i = 0; i < ngroup; i++) delete [] gnames[i];

  if (me == 0) tmp = fread(&ngroup,sizeof(int),1,fp);
  MPI_Bcast(&ngroup,1,MPI_INT,0,world);

  int n;
  for (int i = 0; i < ngroup; i++) {
    if (me == 0) tmp = fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    gnames[i] = new char[n];
    if (me == 0) tmp = fread(gnames[i],sizeof(char),n,fp);
    MPI_Bcast(gnames[i],n,MPI_CHAR,0,world);
  }
}

/* ----------------------------------------------------------------------
   return size of owned surf info for this proc
   using count of all owned surfs
   NOTE: worry about N overflowing int, and in IROUNDUP ???
------------------------------------------------------------------------- */

int Surf::size_restart_one()
{
  int n = sizeof(surfint);
  n = IROUNDUP(n);
  n += 3*sizeof(int);
  n = IROUNDUP(n);
  n += domain->dimension * 3*sizeof(double);
  n = IROUNDUP(n);
  n += sizeof_custom();
  return n;
}

int Surf::size_restart()
{
  int n = sizeof(int);
  n = IROUNDUP(n);
  n += nown * size_restart_one();
  return n;
}

/* ----------------------------------------------------------------------
   pack my owned surf info buf
   include persistent attributes of lines and tris
   include custom data as ints and doubles
   works for explicit/all or explict/distributed
   return n = # of packed bytes
   called by WriteRestart
     implicit surfs are not written to restart files
   NOTE: worry about N overflowing int, and in IROUNDUP ???
------------------------------------------------------------------------- */

int Surf::pack_restart(char *buf)
{
  int start,stop,stride;
  Line *plines;
  Tri *ptris;

  // all versus distributed
  
  if (!distributed) {
    start = me;
    stop = nlocal;
    stride = nprocs;
    plines = lines;
    ptris = tris;
  } else {
    start = 0;
    stop = nown;
    stride = 1;
    plines = mylines;
    ptris = mytris;
  }

  // pack nown = # of my owned surfs
  // pack data for each surf, including custom data
  
  char *ptr = buf;

  int *ibuf = (int *) ptr;
  ibuf[0] = nown;
  ptr += sizeof(int);
  ptr = ROUNDUP(ptr);

  if (domain->dimension == 2) {
    int count = 0;
    for (int m = start; m < stop; m += stride) {
      surfint *sbuf = (surfint *) ptr;
      sbuf[0] = plines[m].id;
      ptr += sizeof(surfint);
      ptr = ROUNDUP(ptr);
      
      int *ibuf = (int *) ptr;
      ibuf[0] = plines[m].type;
      ibuf[1] = plines[m].mask;
      ibuf[2] = plines[m].transparent;
      ptr += 3*sizeof(int);
      ptr = ROUNDUP(ptr);

      double *dbuf = (double *) ptr;
      memcpy(&dbuf[0],plines[m].p1,3*sizeof(double));
      memcpy(&dbuf[3],plines[m].p2,3*sizeof(double));
      ptr += 2 * 3*sizeof(double);
      ptr = ROUNDUP(ptr);

      ptr += pack_custom(count,ptr);
      count++;
    }
    
  } else if (domain->dimension == 3) {
    int count = 0;
    for (int m = start; m < stop; m += stride) {
      surfint *sbuf = (surfint *) ptr;
      sbuf[0] = ptris[m].id;
      ptr += sizeof(surfint);
      ptr = ROUNDUP(ptr);
      
      int *ibuf = (int *) ptr;
      ibuf[0] = ptris[m].type;
      ibuf[1] = ptris[m].mask;
      ibuf[2] = ptris[m].transparent;
      ptr += 3*sizeof(int);
      ptr = ROUNDUP(ptr);

      double *dbuf = (double *) ptr;
      memcpy(&dbuf[0],ptris[m].p1,3*sizeof(double));
      memcpy(&dbuf[3],ptris[m].p2,3*sizeof(double));
      memcpy(&dbuf[6],ptris[m].p3,3*sizeof(double));
      ptr += 3 * 3*sizeof(double);
      ptr = ROUNDUP(ptr);

      ptr += pack_custom(count,ptr);
      count++;
    }
  }
  
  return ptr - buf;
}

/* ----------------------------------------------------------------------
   grow lines or tris data struct
   zero added lines/tris beyond old
---------------------------------------------------------------------- */

void Surf::grow(int old)
{
  if (nmax <= old) return;

  if (domain->dimension == 2) {
    lines = (Surf::Line *)
      memory->srealloc(lines,nmax*sizeof(Line),"surf:lines");
    memset(&lines[old],0,(nmax-old)*sizeof(Line));
  } else {
    tris = (Surf::Tri *)
      memory->srealloc(tris,nmax*sizeof(Tri),"surf:tris");
    memset(&tris[old],0,(nmax-old)*sizeof(Tri));
  }
}

/* ----------------------------------------------------------------------
   grow (or shrink) mylines or mytris data struct
   zero added lines/tris beyond old
---------------------------------------------------------------------- */

void Surf::grow_own(int old)
{
  if (domain->dimension == 2) {
    mylines = (Surf::Line *)
      memory->srealloc(mylines,maxown*sizeof(Line),"surf:mylines");
    if (maxown > old) memset(&mylines[old],0,(maxown-old)*sizeof(Line));
  } else {
    mytris = (Surf::Tri *)
      memory->srealloc(mytris,maxown*sizeof(Tri),"surf:mytris");
    if (maxown > old) memset(&mytris[old],0,(maxown-old)*sizeof(Tri));
  }
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


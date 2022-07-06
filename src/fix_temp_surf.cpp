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
#include "stdlib.h"
#include "string.h"
#include "fix_temp_surf.h"
#include "domain.h"
#include "comm.h"
#include "surf.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "input.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{INT,DOUBLE};                      // several files
enum{COMPUTE,FIX};

// NOTE: should this fix output a temp per surf, b/c name of fix is temp/surf
//       maybe it should access the custom vector when that output is requested
//       output would be useful for debugging, or maybe dump surf is enough ?

/* ---------------------------------------------------------------------- */

FixTempSurf::FixTempSurf(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal fix temp/surf command");

  if (surf->implicit)
    error->all(FLERR,"Cannot use fix temp/surf with implicit surfs");

  distributed = surf->distributed;

  // disable implicit or distributed surfs with fix surf temp for now

  int igroup = surf->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Fix temp/surf group ID does not exist");
  groupbit = surf->bitmask[igroup];

  nevery = atoi(arg[3]);

  id_qw = NULL;

  if (strncmp(arg[4],"c_",2) == 0) {
    source = COMPUTE;
    int n = strlen(arg[4]);
    id_qw = new char[n];
    strcpy(id_qw,&arg[4][2]);

    char *ptr = strchr(id_qw,'[');
    if (ptr) {
      if (id_qw[strlen(id_qw)-1] != ']')
        error->all(FLERR,"Invalid qw in fix temp/surf command");
      qwindex = atoi(ptr+1);
      *ptr = '\0';
    } else qwindex = 0;

    // error checks

    icompute = modify->find_compute(id_qw);
    cqw = modify->compute[icompute];
    if (icompute < 0) error->all(FLERR,"Could not find fix temp/surf compute ID");
    if (cqw->per_surf_flag == 0)
      error->all(FLERR,"Fix temp/surf compute does not compute per-surf info");
    if (qwindex == 0 && cqw->size_per_surf_cols > 0)
      error->all(FLERR,"Fix temp/surf compute does not compute per-surf vector");
    if (qwindex > 0 && cqw->size_per_surf_cols == 0)
      error->all(FLERR,"Fix temp/surf compute does not compute per-surf array");
    if (qwindex > 0 && qwindex > cqw->size_per_surf_cols)
      error->all(FLERR,"Fix temp/surf compute array is accessed out-of-range");
    
  } else if (strncmp(arg[4],"f_",2) == 0) {
    source = FIX;
    int n = strlen(arg[4]);
    id_qw = new char[n];
    strcpy(id_qw,&arg[4][2]);

    char *ptr = strchr(id_qw,'[');
    if (ptr) {
      if (id_qw[strlen(id_qw)-1] != ']')
        error->all(FLERR,"Invalid qw in fix temp/surf command");
      qwindex = atoi(ptr+1);
      *ptr = '\0';
    } else qwindex = 0;

    // error checks

    ifix = modify->find_fix(id_qw);
    fqw = modify->fix[ifix];
    if (ifix < 0) error->all(FLERR,"Could not find fix temp/surf fix ID");
    if (fqw->per_surf_flag == 0)
      error->all(FLERR,"Fix temp/surf fix does not compute per-surf info");
    if (qwindex == 0 && fqw->size_per_surf_cols > 0)
      error->all(FLERR,"Fix temp/surf fix does not compute per-surf vector");
    if (qwindex > 0 && fqw->size_per_surf_cols == 0)
      error->all(FLERR,"Fix temp/surf fix does not compute per-surf array");
    if (qwindex > 0 && qwindex > fqw->size_per_surf_cols)
      error->all(FLERR,"Fix temp/surf fix array is accessed out-of-range");
    if (nevery % fqw->per_surf_freq)
      error->all(FLERR,"Fix temp/surf source not computed at compatible times");
  }

  twall = input->numeric(FLERR,arg[5]);
  if (twall <= 0.0) error->all(FLERR,"Fix temp/surf initial temp <= 0.0");

  emi = input->numeric(FLERR,arg[6]);
  if (emi <= 0.0 || emi > 1.0)
    error->all(FLERR,"Fix temp/surf emissivity must be > 0.0 and <= 1");
  prefactor = 1.0 / (emi *  5.670374419e-8);

  // NOTE: what units are prefactor ?
  // NOTE: what if surfs are setup after this fix is constructed ?
  //       maybe just check that count doesn't change ?
  
  // trigger setup of list of owned surf elements belonging to surf group

  firstflag = 1;

  // initialize data structures

  qw = qw_all = NULL;
  cglobal = clocal = NULL;

  // create per-surf temperature vector

  // NOTE: what if custom vec already exists
  //       could be OK if group for this fix does not overlap other group(s)

  tindex = surf->add_custom((char *) "temperature",DOUBLE,0);
}

/* ---------------------------------------------------------------------- */

FixTempSurf::~FixTempSurf()
{
  delete [] id_qw;
  memory->destroy(qw);
  memory->destroy(qw_all);
  memory->destroy(cglobal);
  memory->destroy(clocal);
  surf->remove_custom(tindex);
}

/* ---------------------------------------------------------------------- */

int FixTempSurf::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTempSurf::init()
{
  // one-time setup of lists of owned elements contributing to fix surf temp
  // NOTE: will need to recalculate, if allow addition of surf elements
  // nown = # of surf elements I own
  // nchoose = # of nown surf elements in surface group
  // cglobal[] = global indices for nchoose elements
  //             used to access lines/tris in Surf
  // clocal[] = local indices for nchoose elements
  //            used to access nown data from per-surf computes,fixes,variables

  // NOTE: should this be done once per run in case surfs change ?
  //       or group assignement changes ?

  if (!firstflag) return;
  firstflag = 0;
  
  int m;
  int me = comm->me;
  int nprocs = comm->nprocs;

  nchoose = 0;
  for (int i = 0; i < nown; i++) {
    if (dimension == 2) {
      if (!distributed) m = me + i*nprocs;
      else m = i;
      if (lines[m].mask & groupbit) nchoose++;
    } else {
      if (!distributed) m = me + i*nprocs;
      else m = i;
      if (tris[m].mask & groupbit) nchoose++;
    }
  }

  memory->create(cglobal,nchoose,"temp/surf:cglobal");
  memory->create(clocal,nchoose,"temp/surf:clocal");

  nchoose = 0;
  for (int i = 0; i < nown; i++) {
    if (dimension == 2) {
      if (!distributed) m = me + i*nprocs;
      else m = i;
      if (lines[m].mask & groupbit) {
        cglobal[nchoose] = m;
        clocal[nchoose++] = i;
      }
    } else {
      if (!distributed) m = me + i*nprocs;
      else m = i;
      if (tris[m].mask & groupbit) {
        cglobal[nchoose] = m;
        clocal[nchoose++] = i;
      }
    }
  }

  // one-time initialize of temperature of surfs in the group

  if (distributed) lines = surf->mylines;
  else lines = surf->lines;
  if (distributed) tris = surf->mytris;
  else tris = surf->tris;

  nown = surf->nown;
  dimension = domain->dimension;

  double *tvector = surf->edvec[surf->ewhich[tindex]];

  // NOTE: should loop bounds be nsurf or nown ??
  //       if nown and distributed, then need some comm ?
  // NOTE: custom vecs are allocated to Nlocal in length
  // NOTE: custom vecs are dumped with nchoose and clocal
  //       that might be for distributed and not

  if (dimension == 2) {
    for (int i = 0; i < nsurf; i++)
      if (lines[i].mask & groupbit) tvector[i] = twall;
      else tvector[i] = 0.0;
  } else {
    for (int i = 0; i < nsurf; i++)
      if (tris[i].mask & groupbit) tvector[i] = twall;
      else tvector[i] = 0.0;
  }


  memory->create(qw,nsurf,"temp/surf:init");
  memory->create(qw_all,nsurf,"temp/surf:init");
  for (int i = 0; i < nsurf; i++) qw[i] = qw_all[i] = 0.0;
}

/* ----------------------------------------------------------------------
   compute new surface element temperatures based on heat flux
   only invoked once every Nevery steps
------------------------------------------------------------------------- */

void FixTempSurf::end_of_step()
{
  // DEBUG
  return;

  // set qw vector from compute or source fix
  // NOTE: if compute, do values need to be collated ?
  
  if (source == COMPUTE) {
    if (qwindex == 0) {
      double *vector = fqw->vector_surf;
      for (int i = 0; i < nchoose; i++)
	qw[cglobal[i]] = vector[clocal[i]];
    } else if (source == FIX) {
      double **array = fqw->array_surf;
      int index = qwindex-1;
      for (int i = 0; i < nchoose; i++)
	qw[cglobal[i]] = array[clocal[i]][index];
    }
  } else if (source == FIX) {
  }

  // NOTE: is this the best way to do this comm ?
  //       what about collate in compute or fix ?
  // NOTE: is this really a distributed quantity coming
  //       from source = compute or fix ?
  
  MPI_Allreduce(qw,qw_all,nsurf,MPI_DOUBLE,MPI_SUM,world);

  // calculation of new temperature for each surf element
  // NOTE: does each proc do this for all surf elements ?
  
  double *tvector = surf->edvec[surf->ewhich[tindex]];

  for (int i = 0; i < nchoose; i++) {
    if (qw_all[cglobal[i]] > 1.0e-6) 
      tvector[cglobal[i]] = pow((prefactor * qw_all[cglobal[i]]),0.25);
    else tvector[cglobal[i]] = twall;
  }
}

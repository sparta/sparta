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

/* ---------------------------------------------------------------------- */

FixTempSurf::FixTempSurf(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal fix temp/surf command");

  // disable implicit or distributed surfs with fix surf temp for now

  implicit = surf->implicit;
  distributed = surf->distributed;
  if (implicit || distributed)
    error->all(FLERR,"Cannot currently use fix temp/surf with "
               "implicit or distributed surfs");

  int igroup = surf->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Fix temp/surf group ID does not exist");
  groupbit = surf->bitmask[igroup];

  nevery = atoi(arg[3]);

  id_qw = NULL;

  if (strncmp(arg[4],"f_",2) == 0) {
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

  // error check

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
    error->all(FLERR,"Fix temp/surf and fix ave/surf not "
               "computed at compatible times");
  }

  twall = input->numeric(FLERR,arg[5]);
  if (twall <= 0.0) error->all(FLERR,"Fix temp/surf initial temp <= 0.0");

  emi = input->numeric(FLERR,arg[6]);
  if (emi <= 0.0 || emi > 1.0)
    error->all(FLERR,"Emissivity of surface must be between 0 and 1");

  prefactor = 1.0 / (emi *  5.670374419e-8);

  // initialize data structures
  // trigger setup of list of owned surf elements belonging to surf group

  if (distributed && !implicit) {
    lines = surf->mylines;
    tris = surf->mytris;
    nsurf = surf->nown;
  }
  else {
    lines = surf->lines;
    tris = surf->tris;
    nsurf = surf->nlocal;
  }

  nown = surf->nown;

  dimension = domain->dimension;
  firstflag = 1;
  qw = qw_all = NULL;
  cglobal = clocal = NULL;

  // create per-surf temperature vector
  // initialize temperature of surfs in the group

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
  if (!firstflag) return;
  firstflag = 0;

  // one-time setup of lists of owned elements contributing to fix surf temp
  // NOTE: will need to recalculate, if allow addition of surf elements
  // nown = # of surf elements I own
  // nchoose = # of nown surf elements in surface group
  // cglobal[] = global indices for nchoose elements
  //             used to access lines/tris in Surf
  // clocal[] = local indices for nchoose elements
  //            used to access nown data from per-surf computes,fixes,variables

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

  memory->create(qw,nsurf,"temp/surf:init");
  memory->create(qw_all,nsurf,"temp/surf:init");
  for (int i = 0; i < nsurf; i++) qw[i] = qw_all[i] = 0.0;
}

/* ----------------------------------------------------------------------
   compute new surface temperature based on heat flux
------------------------------------------------------------------------- */

void FixTempSurf::end_of_step()
{
  if (qwindex == 0) {
    double *vector = fqw->vector_surf;
    for (int i = 0; i < nchoose; i++)
      qw[cglobal[i]] = vector[clocal[i]];
  } else {
    double **array = fqw->array_surf;
    int index = qwindex-1;
    for (int i = 0; i < nchoose; i++)
      qw[cglobal[i]] = array[clocal[i]][index];
  }

  // NOTE: is this the right way to do this ?
  // what about collate in compute or fix ?

  MPI_Allreduce(qw,qw_all,nsurf,MPI_DOUBLE,MPI_SUM,world);

  double *tvector = surf->edvec[surf->ewhich[tindex]];

  for (int i = 0; i < nchoose; i++) {
    if (qw_all[cglobal[i]] > 1.0e-6) 
      tvector[cglobal[i]] = pow((prefactor * qw_all[cglobal[i]]),0.25);
    else tvector[cglobal[i]] = twall;
  }
}

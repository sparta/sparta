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
#include "fix_surf_temp.h"
#include "surf.h"
#include "domain.h"
#include "comm.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "input.h"
#include "output.h"
#include "dump.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

FixSurfTemp::FixSurfTemp(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
//  MPI_Comm_rank(world,&me);

  if (narg < 5) error->all(FLERR,"Illegal fix surf temp command");

  implicit = surf->implicit;
  if (implicit)
    error->all(FLERR,"Cannot use fix surf temp with implicit surfs");

  int igroup = surf->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Fix surf temp group ID does not exist");
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
        error->all(FLERR,"Invalid qw in fix surf temp command");
      qwindex = atoi(ptr+1);
      *ptr = '\0';
    } else qwindex = 0;

  // error check

  ifix = modify->find_fix(id_qw);
  fqw = modify->fix[ifix];
  if (ifix < 0) error->all(FLERR,"Could not find fix surf temp fix ID");
  if (fqw->per_surf_flag == 0)
    error->all(FLERR,"Fix surf temp fix does not "
               "compute per-surf info");
  if (qwindex == 0 && fqw->size_per_surf_cols > 0)
    error->all(FLERR,"Fix surf temp fix does not "
               "compute per-surf vector");
  if (qwindex > 0 && fqw->size_per_surf_cols == 0)
    error->all(FLERR,"Fix surf temp fix does not "
               "compute per-surf array");
  if (qwindex > 0 && qwindex > fqw->size_per_surf_cols)
    error->all(FLERR,"Fix surf temp fix array is "
               "accessed out-of-range");
  if (nevery % fqw->per_surf_freq)
    error->all(FLERR,"Fix surf temp and fix ave/surf not computed at compatible times");
  }

  twall = input->numeric(FLERR,arg[5]);
  if (twall <= 0.0) error->all(FLERR,"Fix surf temp <= 0.0");

  emi = input->numeric(FLERR,arg[6]);
  if (emi <= 0.0 || emi > 1.0)
    error->all(FLERR,"Emissivity of surface has to be between 0 and 1");

  prefactor = 1.0 / (emi *  5.670374419e-8);

  // initialize data structures
  // trigger setup of list of owned surf elements belonging to surf group

  dimension = domain->dimension;
  nown = surf->nown;
  nsurf = surf->nsurf;
  firstflag = 1;
  qw = qw_all = NULL;
  cglobal = NULL;

}

/* ---------------------------------------------------------------------- */

FixSurfTemp::~FixSurfTemp()
{
  delete [] id_qw;
  memory->destroy(qw);
  memory->destroy(qw_all);
  memory->destroy(cglobal);
}

/* ---------------------------------------------------------------------- */

int FixSurfTemp::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSurfTemp::init()
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

    Surf::Line *lines;
    Surf::Tri *tris;
    distributed = surf->distributed;

    if (distributed && !implicit) lines = surf->mylines;
    else lines = surf->lines;
    if (distributed && !implicit) tris = surf->mytris;
    else tris = surf->tris;

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

    memory->create(cglobal,nchoose,"fix surf temp:cglobal");

    nchoose = 0;
    for (int i = 0; i < nown; i++) {
      if (dimension == 2) {
        if (!distributed) m = me + i*nprocs;
        else m = i;
        if (lines[m].mask & groupbit) {
          cglobal[nchoose++] = m;
        }
      } else {
        if (!distributed) m = me + i*nprocs;
        else m = i;
        if (tris[m].mask & groupbit) {
          cglobal[nchoose++] = m;
        }
      }
    }

    memory->create(qw,surf->nlocal,"fix surf temp:init");
    memory->create(qw_all,surf->nlocal,"fix surf temp:init");
    for (int i = 0; i < surf->nlocal; i++) qw[i] = qw_all[i] = 0.0;

}

/* ----------------------------------------------------------------------
   perform surface temperature computation based on heat flux
------------------------------------------------------------------------- */

void FixSurfTemp::end_of_step()
{
//  if (update->ntimestep >= nvalid) {
//    nvalid = (update->ntimestep/fqw->per_surf_freq)*fqw->per_surf_freq +
//      fqw->per_surf_freq;
//    nvalid -= (fqw->per_surf_freq-1)*fqw->nevery;
//    if (nvalid <= update->ntimestep) nvalid += fqw->per_surf_freq;

    if (qwindex == 0) {
      double *vector = fqw->vector_surf;
      for (int i = 0; i < nchoose; i++)
        qw[cglobal[i]] = vector[i];
    }
    else {
      double **array = fqw->array_surf;
      int index = qwindex-1;
      for (int i = 0; i < nchoose; i++)
        qw[cglobal[i]] = array[i][index];
    }

    MPI_Allreduce(qw,qw_all,surf->nlocal,MPI_DOUBLE,MPI_SUM,world);

    for (int i = 0; i < surf->nlocal; i++) {
        if (dimension == 2) {
            Surf::Line *lines = surf->lines;
            if (qw_all[i] > 1.0) lines[i].temp = pow((prefactor * qw_all[i]),0.25);
            else lines[i].temp = twall;
        } else {
            Surf::Tri *tris = surf->tris;
            if (qw_all[i] > 1.0) tris[i].temp = pow((prefactor * qw_all[i]),0.25);
            else tris[i].temp = twall;
        }
    }
//  }
}

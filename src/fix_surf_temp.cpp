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

/* ----------------------------------------------------------------------
   Contributing author: Arnaud Borner (NASA Ames)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_surf_temp.h"
#include "domain.h"
#include "comm.h"
#include "surf.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "input.h"
#include "update.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

// Stefan-Boltzmann constants for different units and dimensions

#define SB_SI 5.670374419e-8
#define SB_CGS 5.670374419e-5

enum{INT,DOUBLE};                      // several files
enum{COMPUTE,FIX};

/* ---------------------------------------------------------------------- */

FixSurfTemp::FixSurfTemp(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg != 8) error->all(FLERR,"Illegal fix surf/temp command");

  if (surf->implicit)
    error->all(FLERR,"Cannot use fix surf/temp with implicit surfs");
  if (surf->distributed)
    error->all(FLERR,"Cannot use fix surf/temp with distributed surfs");

  int igroup = surf->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Fix surf/temp group ID does not exist");
  groupbit = surf->bitmask[igroup];

  nevery = atoi(arg[3]);

  if (strncmp(arg[4],"c_",2) == 0) {
    source = COMPUTE;
    int n = strlen(arg[4]);
    id_qw = new char[n];
    strcpy(id_qw,&arg[4][2]);

    char *ptr = strchr(id_qw,'[');
    if (ptr) {
      if (id_qw[strlen(id_qw)-1] != ']')
        error->all(FLERR,"Invalid source in fix surf/temp command");
      qwindex = atoi(ptr+1);
      *ptr = '\0';
    } else qwindex = 0;

    // error checks

    icompute = modify->find_compute(id_qw);
    cqw = modify->compute[icompute];
    if (icompute < 0) error->all(FLERR,"Could not find fix surf/temp compute ID");
    if (cqw->per_surf_flag == 0)
      error->all(FLERR,"Fix surf/temp compute does not compute per-surf info");
    if (qwindex == 0 && cqw->size_per_surf_cols > 0)
      error->all(FLERR,"Fix surf/temp compute does not compute per-surf vector");
    if (qwindex > 0 && cqw->size_per_surf_cols == 0)
      error->all(FLERR,"Fix surf/temp compute does not compute per-surf array");
    if (qwindex > 0 && qwindex > cqw->size_per_surf_cols)
      error->all(FLERR,"Fix surf/temp compute array is accessed out-of-range");

  } else if (strncmp(arg[4],"f_",2) == 0) {
    source = FIX;
    int n = strlen(arg[4]);
    id_qw = new char[n];
    strcpy(id_qw,&arg[4][2]);

    char *ptr = strchr(id_qw,'[');
    if (ptr) {
      if (id_qw[strlen(id_qw)-1] != ']')
        error->all(FLERR,"Invalid source in fix surf/temp command");
      qwindex = atoi(ptr+1);
      *ptr = '\0';
    } else qwindex = 0;

    // error checks

    ifix = modify->find_fix(id_qw);
    fqw = modify->fix[ifix];
    if (ifix < 0) error->all(FLERR,"Could not find fix surf/temp fix ID");
    if (fqw->per_surf_flag == 0)
      error->all(FLERR,"Fix surf/temp fix does not compute per-surf info");
    if (qwindex == 0 && fqw->size_per_surf_cols > 0)
      error->all(FLERR,"Fix surf/temp fix does not compute per-surf vector");
    if (qwindex > 0 && fqw->size_per_surf_cols == 0)
      error->all(FLERR,"Fix surf/temp fix does not compute per-surf array");
    if (qwindex > 0 && qwindex > fqw->size_per_surf_cols)
      error->all(FLERR,"Fix surf/temp fix array is accessed out-of-range");
    if (nevery % fqw->per_surf_freq)
      error->all(FLERR,"Fix surf/temp source not computed at compatible times");

  } else error->all(FLERR,"Invalid source in fix surf/temp command");

  twall = input->numeric(FLERR,arg[5]);
  if (twall <= 0.0) error->all(FLERR,"Fix surf/temp initial temp <= 0.0");

  emi = input->numeric(FLERR,arg[6]);
  if (emi <= 0.0 || emi > 1.0)
    error->all(FLERR,"Fix surf/temp emissivity must be > 0.0 and <= 1");

  int n = strlen(arg[7]) + 1;
  char *id_custom = new char[n];
  strcpy(id_custom,arg[7]);

  // create per-surf temperature vector

  tindex = surf->add_custom(id_custom,DOUBLE,0);
  delete [] id_custom;

  // prefactor and threshold in Stefan/Boltzmann equation
  // units of prefactor (SI) is K^4 / (watt - m^2)
  // same in 3d vs 2d, since SPARTA treats 2d cell volume as 1 m in z

  int dimension = domain->dimension;

  if (strcmp(update->unit_style,"si") == 0) {
    prefactor = 1.0 / (emi * SB_SI);
    threshold = 1.0e-6;
  } else if (strcmp(update->unit_style,"cgs") == 0) {
    prefactor = 1.0 / (emi * SB_CGS);
    threshold = 1.0e-3;
  }

  // trigger setup of list of owned surf elements belonging to surf group

  firstflag = 1;

  // initialize data structure

  tvector_me = NULL;
}

/* ---------------------------------------------------------------------- */

FixSurfTemp::~FixSurfTemp()
{
  delete [] id_qw;
  memory->destroy(tvector_me);
  surf->remove_custom(tindex);
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
  // one-time initialization of temperature for all surfs in custom vector

  if (!firstflag) return;
  firstflag = 0;

  double *tvector = surf->edvec[surf->ewhich[tindex]];
  int nlocal = surf->nlocal;

  for (int i = 0; i < nlocal; i++)
    tvector[i] = twall;

  // allocate per-surf vector for explicit all surfs

  memory->create(tvector_me,nlocal,"surf/temp:tvector_me");
}

/* ----------------------------------------------------------------------
   compute new surface element temperatures based on heat flux
   only invoked once every Nevery steps
------------------------------------------------------------------------- */

void FixSurfTemp::end_of_step()
{
  int i,m,mask;
  double qw;

  int me = comm->me;
  int nprocs = comm->nprocs;
  int dimension = domain->dimension;

  // access source compute or fix
  // set new temperature via Stefan-Boltzmann eq for nown surfs I own
  // use Twall if surf is not in surf group or eng flux is too small
  // compute/fix output is just my nown surfs, indexed by M
  // store in tvector_me = all nlocal surfs, indexed by I

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;

  int nlocal = surf->nlocal;

  memset(tvector_me,0,nlocal*sizeof(double));

  if (qwindex == 0) {
    double *vector;
    if (source == COMPUTE) {
      cqw->post_process_surf();
      vector = cqw->vector_surf;
    } else vector = fqw->vector_surf;

    m = 0;
    for (i = me; i < nlocal; i += nprocs) {
      if (dimension == 3) mask = tris[i].mask;
      else mask = lines[i].mask;
      if (!(mask & groupbit)) tvector_me[i] = twall;
      else {
        qw = vector[m];
        if (qw > threshold) tvector_me[i] = pow(prefactor*qw,0.25);
        else tvector_me[i] = twall;
      }
      m++;
    }

  } else {
    double **array;
    if (source == COMPUTE) {
      cqw->post_process_surf();
      array = cqw->array_surf;
    } else array = fqw->array_surf;

    int icol = qwindex-1;

    m = 0;
    for (i = me; i < nlocal; i += nprocs) {
      if (dimension == 3) mask = tris[i].mask;
      else mask = lines[i].mask;
      if (!(mask & groupbit)) tvector_me[i] = twall;
      else {
        qw = array[m][icol];
        if (qw > threshold) tvector_me[i] = pow(prefactor*qw,0.25);
        else tvector_me[i] = twall;
      }
      m++;
    }
  }

  // Allreduce tvector_me with my owned surfs to tvector custom variable
  // so that all procs know new temperature of all surfs
  // NOTE: could possibly just Allreduce a vector size of surface group

  double *tvector = surf->edvec[surf->ewhich[tindex]];
  MPI_Allreduce(tvector_me,tvector,nlocal,MPI_DOUBLE,MPI_SUM,world);
}

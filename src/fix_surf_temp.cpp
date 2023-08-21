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

  // trigger one-time initialization of custom per-surf temperatures

  firstflag = 1;
}

/* ---------------------------------------------------------------------- */

FixSurfTemp::~FixSurfTemp()
{
  delete [] id_qw;
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
  if (!firstflag) return;
  firstflag = 0;

  // one-time initialization of temperature for all surfs in custom vector

  double *tvector = surf->edvec[surf->ewhich[tindex]];
  int nsown = surf->nown;

  for (int i = 0; i < nsown; i++)
    tvector[i] = twall;
}

/* ----------------------------------------------------------------------
   compute new surface element temperatures based on heat flux
   only invoked once every Nevery steps
------------------------------------------------------------------------- */

void FixSurfTemp::end_of_step()
{
  int m,mask;
  double qw;

  int me = comm->me;
  int nprocs = comm->nprocs;
  int dimension = domain->dimension;
  int distributed = surf->distributed;
  
  // access source compute or fix which is only surfs I own
  // set new temperature via Stefan-Boltzmann eq for nown surfs I own
  // NOTE: which of these 2 options (set doc page accordingly):
  //   set to Twall if eng flux to surf is too small
  //   no reset if eng flux < threshold

  Surf::Line *lines;
  Surf::Tri *tris;

  if (distributed) {
    lines = surf->mylines;
    tris = surf->mytris;
  } else {
    lines = surf->lines;
    tris = surf->tris;
  }
  
  double *tcustom = surf->edvec[surf->ewhich[tindex]];
  int nsown = surf->nown;
  
  if (qwindex == 0) {
    double *qwvector;
    if (source == COMPUTE) {
      cqw->post_process_surf();
      qwvector = cqw->vector_surf;
    } else qwvector = fqw->vector_surf;

    for (int i = 0; i < nsown; i++) {
      if (!distributed) m = me + i*nprocs;
      else m = i;
      if (dimension == 2) mask = lines[m].mask;
      else mask = tris[m].mask;
      if (mask & groupbit) {
        qw = qwvector[i];
	if (qw > threshold) tcustom[i] = pow(prefactor*qw,0.25);
        else tcustom[i] = twall;
      }
    }

  } else {
    double **qwarray;
    if (source == COMPUTE) {
      cqw->post_process_surf();
      qwarray = cqw->array_surf;
    } else qwarray = fqw->array_surf;

    int icol = qwindex-1;

    for (int i = 0; i < nsown; i++) {
      if (!distributed) m = me + i*nprocs;
      else m = i;
      if (dimension == 2) mask = lines[m].mask;
      else mask = tris[m].mask;
      if (mask & groupbit) {
        qw = qwarray[i][icol];
        if (qw > threshold) tcustom[i] = pow(prefactor*qw,0.25);
        else tcustom[i] = twall;
      }
    }
  }

  // flag custom attribute as updated

  surf->estatus[tindex] = 0;
}

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

#include "mpi.h"
#include "ctype.h"
#include "string.h"
#include "surf_collide.h"
#include "update.h"
#include "surf.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{INT,DOUBLE};                        // several files
enum{NUMERIC,CUSTOM,VARIABLE,VAREQUAL,VARSURF};   // surf_collide classes

/* ---------------------------------------------------------------------- */

SurfCollide::SurfCollide(SPARTA *sparta, int, char **arg) :
  Pointers(sparta)
{
  // ID and style
  // ID must be all alphanumeric chars or underscores

  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  for (int i = 0; i < n-1; i++)
    if (!isalnum(id[i]) && id[i] != '_')
      error->all(FLERR,"Surf_collide ID must be alphanumeric or "
                 "underscore characters");

  n = strlen(arg[1]) + 1;
  style = new char[n];
  strcpy(style,arg[1]);

  dynamicflag = 0;
  allowreact = 1;
  transparent = 0;
  vector_flag = 1;
  size_vector = 2;

  nsingle = ntotal = 0;
  tname = NULL;

  n_owned = n_localghost = 0;
  t_owned = t_localghost = NULL;

  kokkosable = copy = uncopy = copymode = 0;
}

/* ---------------------------------------------------------------------- */

SurfCollide::~SurfCollide()
{
  if (copy) return;

  delete [] id;
  delete [] style;
  delete [] tname;

  memory->destroy(t_owned);
  memory->destroy(t_localghost);
}

/* ---------------------------------------------------------------------- */

void SurfCollide::init()
{
  nsingle = ntotal = 0;
}

/* ---------------------------------------------------------------------- */

void SurfCollide::tally_reset()
{
  nsingle = 0;
}

/* ---------------------------------------------------------------------- */

void SurfCollide::tally_update()
{
  ntotal += nsingle;
}

/* ---------------------------------------------------------------------- */

void SurfCollide::parse_tsurf(char *str)
{
  tname = NULL;
  tfreq = 1;

  if (strstr(str,"v_") == str) {
    tmode = VARIABLE;
    dynamicflag = 1;
    int n = strlen(&str[2]) + 1;
    tname = new char[n];
    strcpy(tname,&str[2]);
  } else if (strstr(str,"s_") == str) {
    tmode = CUSTOM;
    dynamicflag = 1;
    int n = strlen(&str[2]) + 1;
    tname = new char[n];
    strcpy(tname,&str[2]);
  } else {
    tmode = NUMERIC;
    tsurf = input->numeric(FLERR,str);
    if (tsurf <= 0.0) error->all(FLERR,"Surf_collide tsurf <= 0.0");
  }

  // error checks

  if (tmode == VARIABLE) {
    tindex_var = input->variable->find(tname);
    if (tindex_var < 0)
      error->all(FLERR,"Surf_collide tsurf variable name does not exist");
    if (input->variable->equal_style(tindex_var)) tmode = VAREQUAL;
    else if (input->variable->surf_style(tindex_var)) tmode = VARSURF;
    else error->all(FLERR,"Surf_collide tsurf variable is invalid style");

  } else if (tmode == CUSTOM) {
    tindex_custom = surf->find_custom(tname);
    if (tindex_custom < 0)
      error->all(FLERR,"Surf_collide tsurf could not find custom attribute");
    if (surf->etype[tindex_custom] != DOUBLE)
      error->all(FLERR,"Surf_collide tsurf custom attribute is not a float");
    if (surf->esize[tindex_custom] > 0)
      error->all(FLERR,"Surf_collide tsurf custom attribute is not a vector");
  }
}

/* ---------------------------------------------------------------------- */

void SurfCollide::check_tsurf()
{
  if (tmode == VAREQUAL || tmode == VARSURF) {
    tindex_var = input->variable->find(tname);
    if (tindex_var < 0)
      error->all(FLERR,"Surf_collide tsurf variable name does not exist");
  } else if (tmode == CUSTOM) {
    int tindex_custom = surf->find_custom(tname);
    if (tindex_custom < 0)
      error->all(FLERR,"Surf_collide tsurf could not find custom attribute");
  }

  persurf_temperature = 0;
  if (tmode == VARSURF || tmode == CUSTOM) persurf_temperature = 1;
}

/* ----------------------------------------------------------------------
   recalculate Tsurf values which are dynamic
   called by Update::setup() and Update::run()
---------------------------------------------------------------------- */

void SurfCollide::dynamic()
{
  // VAREQUAL mode
  // equal-style variable sets single tsurf value for all surfs

  if (tmode == VAREQUAL) {

    // only evaluate variable if timestep is multiple of tfreq

    if (update->ntimestep % tfreq) return;
    tsurf = input->variable->compute_equal(tindex_var);
    if (tsurf <= 0.0) error->all(FLERR,"Surf_collide tsurf <= 0.0");

  // VARSURF mode
  // surf-style variable sets new tsurf values for all surfs
  // particle/surf collisions access t_persurf for local+ghost values

  } else if (tmode == VARSURF) {

    // only evaluate variable if timestep is multiple of tfreq

    int spreadflag = 0;
    if (update->ntimestep % tfreq == 0) {
      if (n_owned != surf->nown) {
	memory->destroy(t_owned);
	n_owned = surf->nown;
	memory->create(t_owned,n_owned,"surfcollide:t_owned");
      }

      input->variable->compute_surf(tindex_var,t_owned,1,0);
      spreadflag = 1;
    }

    // spread t_owned values to t_localghost values via spread_own2local()
    // if just re-computed variable OR surfs are
    //   distributed and load balance/adaptation took place on previous step

    if (spreadflag ||
	(surf->distributed && surf->localghost_changed_step == update->ntimestep-1)) {
      if (n_localghost != surf->nlocal + surf->nghost) {
	memory->destroy(t_localghost);
	n_localghost = surf->nlocal + surf->nghost;
	memory->create(t_localghost,n_localghost,"surfcollide:t_localghost");
      }

      surf->spread_own2local(1,DOUBLE,t_owned,t_localghost);
      t_persurf = t_localghost;
    }

  // CUSTOM mode
  // ensure access to custom per-surf vec for tsurf values for all surfs
  // particle/surf collisions access t_persurf for local+ghost values

  } else if (tmode == CUSTOM) {

    // spread owned values to local+ghost values via spread_custom()
    // estatus == 1 means owned values already spread to local+ghost values
    // if estatus == 0: owned values are new OR
    //   surfs are distributed and load balance/adaptation took place

    if (surf->estatus[tindex_custom] == 0) surf->spread_custom(tindex_custom);
    t_persurf = surf->edvec_local[tindex_custom];
  }
}

/* ---------------------------------------------------------------------- */

double SurfCollide::compute_vector(int i)
{
  one[0] = nsingle;
  one[1] = ntotal;
  MPI_Allreduce(one,all,2,MPI_DOUBLE,MPI_SUM,world);

  return all[i];
}

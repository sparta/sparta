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
  
  copy = copymode = 0;
}

/* ---------------------------------------------------------------------- */

SurfCollide::~SurfCollide()
{
  if (copy) return;

  delete [] id;
  delete [] style;
  delete [] tname;
}

/* ---------------------------------------------------------------------- */

void SurfCollide::init()
{
  nsingle = ntotal = 0;

  // NOTE: need to always invoke var/custom on this step?
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

  // allocate memory for t_owned and t_local
  
  if (tmode == VARSURF) {
  }
}

/* ----------------------------------------------------------------------
   recalculate Tsurf values which are dynamics
   called by Update::setup() and Update::run()
---------------------------------------------------------------------- */

void SurfCollide::dynamic()
{
  // VAREQUAL mode
  // evaulate equal-style variable to set new tsurf for all surfs
  // only evaluate if timestep is multiple of tfreq
 
  if (tmode == VAREQUAL) {
    if (update->ntimestep % tfreq) return;
    tsurf = input->variable->compute_equal(tindex_var);
    if (tsurf <= 0.0) error->all(FLERR,"Surf_collide tsurf <= 0.0");

  // VARSURF mode
  // evaulate surf-style variable to set new tsurf values for all surfs
  // only evaluate if timestep is multiple of tfreq
  // surf-style variable sets values for owned explicit surfs only
  // comm owned values to Surf::lines/tris via spread_vector()
  // particle/surf collisions access t_persurf for local+ghost values 
    
  } else if (tmode == VARSURF) {
    if (update->ntimestep % tfreq) return;
    input->variable->compute_surf(tindex_var,t_owned,1,0);

    int flag = 0;
    int nsown = surf->nown;
    for (int i = 0; i < nsown; i++)
      if (t_owned[i] <= 0.0) flag = 1;
    int flagall;
    MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
    if (flagall) error->all(FLERR,"Surf_collide tsurf <= 0.0 for one or more surfs");
    
    surf->spread_vector(t_owned,t_local);
    t_persurf = t_local;

  // CUSTOM mode
  // access custom per-surf vec for tsurf values for all surfs
  // only necessary if status of custom vector has changed
  // custom vector stores values for owned explicit surfs only
  // comm owned values to Surf::lines/tris via spread_custom_vector()
  // particle/surf collisions access t_persurf for local+ghost values 
    
  } else if (tmode == CUSTOM) {
    if (surf->estatus[tindex_custom] == 0) return;

    double *tcustom = surf->edvec[tindex_custom];
    int nsown = surf->nown;
    int flag = 0;
    for (int i = 0; i < nsown; i++)
      if (tcustom[i] <= 0.0) flag = 1;
    int flagall;
    MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
    if (flagall) error->all(FLERR,"Surf_collide tsurf <= 0.0 for one or more surfs");

    surf->spread_custom(tindex_custom);
    surf->estatus[tindex_custom] = 1;
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


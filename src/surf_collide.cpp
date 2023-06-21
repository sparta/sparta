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

  copy = copymode = 0;
}

/* ---------------------------------------------------------------------- */

SurfCollide::~SurfCollide()
{
  if (copy) return;

  delete [] id;
  delete [] style;
  delete [] tstr;
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
  tstr = NULL;
  tfreq = 1;

  if (strstr(str,"v_") == str) {
    tmode = VARIABLE;
    dynamicflag = 1;
    int n = strlen(&str[2]) + 1;
    tstr = new char[n];
    strcpy(tstr,&str[2]);
  } else if (strstr(str,"s_") == str) {
    tmode = CUSTOM;
    dynamicflag = 1;
    int n = strlen(&str[2]) + 1;
    tstr = new char[n];
    strcpy(tstr,&str[2]);
  } else {
    tmode = NUMERIC;
    tsurf = input->numeric(FLERR,str);
    if (tsurf <= 0.0) error->all(FLERR,"Surf_collide temp <= 0.0");
  }

  // error checks
  
  if (tmode == VARIABLE) {
    tindex_var = input->variable->find(tstr);
    if (tindex_var < 0)
      error->all(FLERR,"Surf_collide tsurf variable name does not exist");
    if (input->variable->equal_style(tindex_var)) tmode = VAREQUAL;
    else if (input->variable->surf_style(tindex_var)) tmode = VARSURF;
    else error->all(FLERR,"Surf_collide tsurf variable is invalid style");
  } else if (tmode == CUSTOM) {
    tindex_custom = surf->find_custom(tstr);
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
    tindex_var = input->variable->find(tstr);
    if (tindex_var < 0)
      error->all(FLERR,"Surf_collide tsurf variable name does not exist");
  } else if (tmode == CUSTOM) {
    int tindex_custom = surf->find_custom(tstr);
    if (tindex_custom < 0)
      error->all(FLERR,"Surf_collide tsurf could not find custom attribute");
  }

  // NOTE: allocate per-surf variable memory ?
}

/* ---------------------------------------------------------------------- */

void SurfCollide::dynamic()
{
  // for VAREQUAL mode
  // evaulate equal-style variable to set new tsurf
  // only if timestep is multiple of tfreq
 
  if (tmode == VAREQUAL) {
    if (update->ntimestep % tfreq) return;
    tsurf = input->variable->compute_equal(tindex_var);
    if (tsurf <= 0.0) error->all(FLERR,"Surf_collide tsurf <= 0.0");

  // for VARSURF mode
  // evaulate surf-style variable to set new t_persurf_own
  // only if timestep is multiple of tfreq
  // communicate t_persurf_own -> t_persurf_nlocal internal to Surf class
  // set t_persurf_nlocal pointer
    
  } else if (tmode == VARSURF) {
    if (update->ntimestep % tfreq) return;
    // when to allocate this vector?
    input->variable->compute_surf(tindex_var,t_persurf_own,1,0);
    // check if any value <= 0.0 ?
    // if (tsurf <= 0.0) error->all(FLERR,"Surf_collide tsurf <= 0.0");
    // trigger own -> nlocal comm ?

  // for CUSTOM mode
  // if status of custom per-surf vec has not changed, just return
  // set t_persurf_nlocal pointer
  // set t_persurf_own to point to custom per-surf vector
  // communicate t_persurf_own -> t_persurf_nlocal internal to Surf class
  // set t_persurf_nlocal pointer
    
  } else if (tmode == CUSTOM) {
    //if (surf->estatus[tindex_custom] == 0) continue;
    t_persurf_own = surf->edvec[surf->ewhich[tindex_custom]];
    //surf->estatus[tindex_custom] = 1;
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


/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
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
#include "fix_custom.h"
#include "custom.h"
#include "domain.h"
#include "grid.h"
#include "input.h"
#include "mixture.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{SET,FILESTYLE,REMOVE};
enum{EQUAL,PARTICLE,GRID,SURF};

/* ---------------------------------------------------------------------- */

FixCustom::FixCustom(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg != 8) error->all(FLERR,"Illegal fix custom command");

  nevery = atoi(arg[2]);

  int n = strlen(arg[3]) + 1;
  fname = new char[n];
  strcpy(fname,arg[3]);

  // NOTE: should file mode allow multiple custom attributes?
  //       so that one file can reset multiple attributes
  // NOTE: but set mode does not allow for multiple
  // NOTE: do we need optional keywords like custom command ?
  //       don't think so b/c custom attribute has to already exist ?
  //         or can this command create the attribute?
  //         this has to do with whether fix custom is invoked on step 0 ?

  if (strcmp(arg[4],"particle") == 0) mode = PARTICLE;
  else if (strcmp(arg[4],"grid") == 0) mode = GRID;
  else if (strcmp(arg[4],"surf") == 0) mode = SURF;
  else error->all(FLERR,"Illegal custom command");

  if (mode == PARTICLE && !particle->exist)
    error->all(FLERR,"Cannot use fix custom particle before particles are defined");
  if (mode == GRID && !grid->exist)
    error->all(FLERR,"Cannot use fix custom grid before a grid is defined");
  if (mode == SURF && !surf->exist)
    error->all(FLERR,"Cannot use fix custom surf before surfaces are defined");

  // attribute name
  // NOTE: should allow vec or array or single array col ?
  
  n = strlen(arg[5]) + 1;
  aname = new char[n];
  strcpy(aname,arg[5]);

  char *ptr = strchr(aname,'[');
  if (ptr) {
    if (aname[strlen(aname)-1] != ']')
      error->all(FLERR,"Fix custom command attribute name is invalid");
    ccol = atoi(ptr+1);
    *ptr = '\0';
  } else ccol = 0;

  // action

  if (strcmp(arg[2],"set") == 0) action = SET;
  else if (strcmp(arg[2],"file") == 0) action = FILESTYLE;
  else error->all(FLERR,"Illegal fix custom command action");

  // read args to reset a custom attribute via a variable
  // NOTE: should store these args and make the invocation here be in a function
  // NOTE: likewise make invocation by custom command be in the same function

  // cindex = index of existing custom attribute
  // otherwise create custom attribute if it does not exist
  // add_custom() initializes all values to zero

  if (mode == PARTICLE) {
    cindex = particle->find_custom(aname);
    if (cindex >= 0) {
      ctype = particle->etype[cindex];
      csize = particle->esize[cindex];
    } else {
      cindex = particle->add_custom(aname,ctype,csize);
    }

  } else if (mode == GRID) {
    cindex = grid->find_custom(aname);
    if (cindex >= 0) {
      ctype = grid->etype[cindex];
      csize = grid->esize[cindex];
    } else {
      cindex = grid->add_custom(aname,ctype,csize);
    }

  } else if (mode == SURF) {
    cindex = surf->find_custom(aname);
    if (cindex >= 0) {
      ctype = surf->etype[cindex];
      csize = surf->esize[cindex];
    } else {
      cindex = surf->add_custom(aname,ctype,csize);
    }
  }

  int iarg;
  vname = fname = NULL;

  if (action == SET) {

    if (narg < 6) error->all(FLERR,"Illegal custom command");
    
    // variable name

    variable = input->variable;

    if (strncmp(arg[3],"v_",2) == 0) {
      int n = strlen(arg[3]);
      vname = new char[n];
      strcpy(vname,&arg[3][2]);

      vindex = variable->find(vname);
      if (vindex < 0) error->all(FLERR,"Fix custom variable name does not exist");
      if (variable->equal_style(vindex)) vstyle = EQUAL;
      else if (variable->particle_style(vindex)) vstyle = PARTICLE;
      else if (variable->grid_style(vindex)) vstyle = GRID;
      else if (variable->surf_style(vindex)) vstyle = SURF;
      else error->all(FLERR,"Fix custom variable style is invalid");
      if (vstyle != EQUAL && vstyle != mode)
        error->all(FLERR,"Fix custom variable style is invalid");
      
    } else error->all(FLERR,"Fix custom variable name is invalid");

    // mixture or group ID

    if (mode == PARTICLE) {
      int imix = particle->find_mixture(arg[4]);
      if (imix < 0) error->all(FLERR,"Fix custom mixture ID does not exist");
      mixture = particle->mixture[imix];
      mixture->init();  // NOTE: do this here ?
    } else if (mode == GRID) {
      int igroup = grid->find_group(arg[4]);
      if (igroup < 0) error->all(FLERR,"Fix custom grid group ID does not exist");
      groupbit = grid->bitmask[igroup];
    } else if (mode == SURF) {
      int igroup = surf->find_group(arg[4]);
      if (igroup < 0) error->all(FLERR,"Fix custom surf group ID does not exist");
      groupbit = surf->bitmask[igroup];
    }

    // region ID

    if (strcmp(arg[5],"NULL") == 0) region = NULL;
    else {
      int iregion = domain->find_region(arg[5]);
      if (iregion < 0) error->all(FLERR,"Fix custom region ID does not exist");
      region = domain->regions[iregion];
    }

    iarg = 6;
  }

  // read args to reset a custom attribute via a file

  if (action == FILESTYLE) {

    if (narg < 4) error->all(FLERR,"Illegal custom command");
    if (mode == PARTICLE)
      error->all(FLERR,"Custom command cannot use action file with style particle");

    // file name
    
    int n = strlen(arg[3]) + 1;
    fname = new char[n];
    strcpy(fname,arg[3]);

    iarg = 4;
  }

  custom = new Custom(sparta);
}

/* ---------------------------------------------------------------------- */

FixCustom::~FixCustom()
{
}

/* ---------------------------------------------------------------------- */

int FixCustom::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCustom::init()
{
}

/* ----------------------------------------------------------------------
   invoked once every Nevery steps
------------------------------------------------------------------------- */

void FixCustom::end_of_step()
{
  int count;

  // for action SET, invoke variable and set attributes

  if (action == SET) {
    double scalar = 0.0;
    double *vector = NULL;

    if (vstyle == EQUAL) {
      scalar = variable->compute_equal(vindex);
    } else if (vstyle == PARTICLE) {
      memory->create(vector,particle->nlocal,"custom:vector");
      variable->compute_particle(vindex,vector,1,0);
    } else if (vstyle == GRID) {
      memory->create(vector,grid->nlocal,"custom:vector");
      variable->compute_grid(vindex,vector,1,0);
    } else if (vstyle == SURF) {
      memory->create(vector,surf->nown,"custom:vector");
      variable->compute_surf(vindex,vector,1,0);
    }

    // assign value(s) to custom attribute
    // convert to integer if necessary
    // no assignment if particle/grid/surf not in mixture or group or region

    /*
    if (mode == PARTICLE)
      count = custom->set_particle(cindex,ctype,csize,ccol,scalar,vector);
    else if (mode == GRID)
      count = custom->set_grid(cindex,ctype,csize,ccol,scalar,vector);
    else if (mode == SURF)
      count = custom->set_surf(cindex,ctype,csize,ccol,scalar,vector);
    */
    memory->destroy(vector);
  }
  
  // for action FILESTYLE, read file and set attributes
  // count = # of changed attributes

  /*
  if (action == FILESTYLE) {
    count = custom->read_file(mode,cindex,ctype,csize,ccol,fname,0);
  }
  */
  
  // for mode = SURF, set estatus of custom vec/array to 0

  if (mode == SURF) surf->estatus[cindex] = 0;
}

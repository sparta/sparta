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

#include "stdio.h"
#include "string.h"
#include "modify.h"
#include "domain.h"
#include "update.h"
#include "compute.h"
#include "fix.h"
#include "style_compute.h"
#include "style_fix.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

#define DELTA 4

// mask settings - same as in fix.cpp

#define START_OF_STEP  1
#define END_OF_STEP    2

/* ---------------------------------------------------------------------- */

Modify::Modify(SPARTA *sparta) : Pointers(sparta)
{
  nfix = maxfix = 0;
  n_start_of_step = n_end_of_step = 0;

  fix = NULL;
  fmask = NULL;
  list_start_of_step = list_end_of_step = NULL;

  end_of_step_every = NULL;
  list_pergrid = NULL;
  list_update_custom = NULL;
  list_gas_react = NULL;
  list_surf_react = NULL;
  list_timeflag = NULL;

  ncompute = maxcompute = 0;
  compute = NULL;

  // n_pergrid needs to be initialized here because ReadSurf calls
  //  Modify::reset_grid_count without calling Modify::init

  n_pergrid = 0;
}

/* ---------------------------------------------------------------------- */

Modify::~Modify()
{
  // delete all fixes
  // do it via delete_fix() so callbacks are also updated correctly

  while (nfix) delete_fix(fix[0]->id);
  memory->sfree(fix);
  memory->destroy(fmask);

  // delete all computes

  for (int i = 0; i < ncompute; i++) delete compute[i];
  memory->sfree(compute);

  delete [] list_start_of_step;
  delete [] list_end_of_step;

  delete [] end_of_step_every;
  delete [] list_pergrid;
  delete [] list_update_custom;
  delete [] list_gas_react;
  delete [] list_surf_react;
  delete [] list_timeflag;
}

/* ----------------------------------------------------------------------
   initialize all fixes and computes
------------------------------------------------------------------------- */

void Modify::init()
{
  int i;

  // create lists of fixes with masks for calling at each stage of run

  list_init(START_OF_STEP,n_start_of_step,list_start_of_step);
  list_init_end_of_step(END_OF_STEP,n_end_of_step,list_end_of_step);

  // create other lists of fixes and computes

  list_init_fixes();
  list_init_computes();

  // init each fix

  for (i = 0; i < nfix; i++) fix[i]->init();

  // init each compute
  // set invoked_scalar,vector,etc to -1 to force new run to re-compute them
  // add initial timestep to all computes that store invocation times
  //   since any of them may be invoked by initial thermo
  // do not clear out invocation times stored within a compute,
  //   b/c some may be holdovers from previous run, like for ave fixes

  for (i = 0; i < ncompute; i++) {
    compute[i]->init();
    compute[i]->invoked_scalar = -1;
    compute[i]->invoked_vector = -1;
    compute[i]->invoked_array = -1;
    compute[i]->invoked_per_particle = -1;
    compute[i]->invoked_per_grid = -1;
    compute[i]->invoked_per_surf = -1;
  }
  addstep_compute_all(update->ntimestep);
}

/* ---------------------------------------------------------------------- */

void Modify::setup()
{
  // setup each fix

  for (int i = 0; i < nfix; i++) fix[i]->setup();
}

/* ----------------------------------------------------------------------
   start-of-timestep call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::start_of_step()
{
  for (int i = 0; i < n_start_of_step; i++)
    fix[list_start_of_step[i]]->start_of_step();
}

/* ----------------------------------------------------------------------
   end-of-timestep call, only for relevant fixes
   only call fix->end_of_step() on timesteps that are multiples of nevery
------------------------------------------------------------------------- */

void Modify::end_of_step()
{
  for (int i = 0; i < n_end_of_step; i++)
    if (update->ntimestep % end_of_step_every[i] == 0)
      fix[list_end_of_step[i]]->end_of_step();
}

/* ----------------------------------------------------------------------
   pack_grid_one call, only for relevant fixes
   invoked by load balancer when grid cells migrate
------------------------------------------------------------------------- */

int Modify::pack_grid_one(int icell, char *buf, int memflag)
{
  char *ptr = buf;
  for (int i = 0; i < n_pergrid; i++)
    ptr += fix[list_pergrid[i]]->pack_grid_one(icell,ptr,memflag);
  return ptr-buf;
}

/* ----------------------------------------------------------------------
   unpack_grid_one call, only for relevant fixes
   invoked by load balancer when grid cells migrate
------------------------------------------------------------------------- */

int Modify::unpack_grid_one(int icell, char *buf)
{
  char *ptr = buf;
  for (int i = 0; i < n_pergrid; i++)
    ptr += fix[list_pergrid[i]]->unpack_grid_one(icell,ptr);
  return ptr-buf;
}

/* ----------------------------------------------------------------------
   copy_grid call, only for relevant fixes
   invoked when a grod cell is removed
------------------------------------------------------------------------- */

void Modify::copy_grid_one(int icell, int jcell)
{
  for (int i = 0; i < n_pergrid; i++)
    fix[list_pergrid[i]]->copy_grid_one(icell,jcell);
}

/* ----------------------------------------------------------------------
   add_grid_one call, only for relevant fixes
   invoked by adapt_grid and fix adapt when new child cells are created
------------------------------------------------------------------------- */

void Modify::add_grid_one()
{
  for (int i = 0; i < n_pergrid; i++)
    fix[list_pergrid[i]]->add_grid_one();
}

/* ----------------------------------------------------------------------
   reset_grid call, only for relevant fixes
   invoked after all grid cell removals
------------------------------------------------------------------------- */

void Modify::reset_grid_count(int nlocal)
{
  for (int i = 0; i < n_pergrid; i++)
    fix[list_pergrid[i]]->reset_grid_count(nlocal);
}

/* ----------------------------------------------------------------------
   grid_changed call, only for relevant fixes
   invoked after per-processor list of grid cells has changed
------------------------------------------------------------------------- */

void Modify::grid_changed()
{
  for (int i = 0; i < n_pergrid; i++)
    fix[list_pergrid[i]]->grid_changed();
}

/* ----------------------------------------------------------------------
   invoke update_custom() method, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::update_custom(int index, double temp_thermal,
                          double temp_rot, double temp_vib, double *vstream)
{
  for (int i = 0; i < n_update_custom; i++)
    fix[list_update_custom[i]]->update_custom(index,temp_thermal,temp_rot,
                                              temp_vib,vstream);
}

/* ----------------------------------------------------------------------
   invoke gas_react() method, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::gas_react(int index)
{
  for (int i = 0; i < n_gas_react; i++)
    fix[list_gas_react[i]]->gas_react(index);
}

/* ----------------------------------------------------------------------
   invoke surf_react() method, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::surf_react(Particle::OnePart *iorig, int &i, int &j)
{
  for (int m = 0; m < n_surf_react; m++) {
    fix[list_surf_react[m]]->surf_react(iorig,i,j);
  }
}

/* ----------------------------------------------------------------------
   add a new fix or replace one with same ID
------------------------------------------------------------------------- */

void Modify::add_fix(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Fix command before simulation box is defined");
  if (narg < 2) error->all(FLERR,"Illegal fix command");

  // if fix ID exists:
  //   set newflag = 0 so create new fix in same location in fix list
  //   error if new style does not match old style
  //     since can't replace it (all when-to-invoke ptrs would be invalid)
  //   delete old fix, but do not call update_callback(),
  //     since will replace this fix and thus other fix locs will not change
  //   set ptr to NULL in case new fix scans list of fixes,
  //     e.g. scan will occur in add_callback() if called by new fix
  // if fix ID does not exist:
  //   set newflag = 1 so create new fix
  //   extend fix and fmask lists as necessary

  int ifix,newflag;
  for (ifix = 0; ifix < nfix; ifix++)
    if (strcmp(arg[0],fix[ifix]->id) == 0) break;

  if (ifix < nfix) {
    newflag = 0;
    if (strcmp(arg[1],fix[ifix]->style) != 0)
      error->all(FLERR,"Replacing a fix, but new style != old style");
    delete fix[ifix];
    fix[ifix] = NULL;
  } else {
    newflag = 1;
    if (nfix == maxfix) {
      maxfix += DELTA;
      fix = (Fix **) memory->srealloc(fix,maxfix*sizeof(Fix *),"modify:fix");
      memory->grow(fmask,maxfix,"modify:fmask");
    }
  }

  // create the Fix

  int found = 0;
  if (sparta->suffix_enable) {
    if (sparta->suffix) {
      char estyle[256];
      sprintf(estyle,"%s/%s",arg[1],sparta->suffix);

      if (0) return;

#define FIX_CLASS
#define FixStyle(key,Class) \
      else if (strcmp(estyle,#key) == 0) { \
        fix[ifix] = new Class(sparta,narg,arg); \
        found = 1; \
      }
#include "style_fix.h"
#undef FixStyle
#undef FIX_CLASS
    }
  }

  if (!found) {

    if (0) return;

#define FIX_CLASS
#define FixStyle(key,Class) \
    else if (strcmp(arg[1],#key) == 0) fix[ifix] = new Class(sparta,narg,arg);
#include "style_fix.h"
#undef FixStyle
#undef FIX_CLASS

    else error->all(FLERR,"Unrecognized fix style");
  }

  // set fix mask values and increment nfix (if new)

  fmask[ifix] = fix[ifix]->setmask();
  if (newflag) nfix++;
}

/* ----------------------------------------------------------------------
   delete a Fix from list of Fixes
------------------------------------------------------------------------- */

void Modify::delete_fix(const char *id)
{
  int ifix = find_fix(id);
  if (ifix < 0) error->all(FLERR,"Could not find fix ID to delete");
  delete fix[ifix];
  //atom->update_callback(ifix);

  // move other Fixes and fmask down in list one slot

  for (int i = ifix+1; i < nfix; i++) fix[i-1] = fix[i];
  for (int i = ifix+1; i < nfix; i++) fmask[i-1] = fmask[i];
  nfix--;
}

/* ----------------------------------------------------------------------
   find a fix by ID
   return index of fix or -1 if not found
------------------------------------------------------------------------- */

int Modify::find_fix(const char *id)
{
  int ifix;
  for (ifix = 0; ifix < nfix; ifix++)
    if (strcmp(id,fix[ifix]->id) == 0) break;
  if (ifix == nfix) return -1;
  return ifix;
}

/* ----------------------------------------------------------------------
   add a new compute
------------------------------------------------------------------------- */

void Modify::add_compute(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal compute command");

  // error check

  for (int icompute = 0; icompute < ncompute; icompute++)
    if (strcmp(arg[0],compute[icompute]->id) == 0)
      error->all(FLERR,"Reuse of compute ID");

  // extend Compute list if necessary

  if (ncompute == maxcompute) {
    maxcompute += DELTA;
    compute = (Compute **)
      memory->srealloc(compute,maxcompute*sizeof(Compute *),"modify:compute");
  }

  // create the Compute

  int found = 0;
  if (sparta->suffix_enable) {
    if (sparta->suffix) {
      char estyle[256];
      sprintf(estyle,"%s/%s",arg[1],sparta->suffix);

      if (0) return;

#define COMPUTE_CLASS
#define ComputeStyle(key,Class) \
      else if (strcmp(estyle,#key) == 0) { \
        compute[ncompute] = new Class(sparta,narg,arg); \
        found = 1; \
      }
#include "style_compute.h"
#undef ComputeStyle
#undef COMPUTE_CLASS
    }
  }

  if (!found) {

    if (0) return;

#define COMPUTE_CLASS
#define ComputeStyle(key,Class) \
    else if (strcmp(arg[1],#key) == 0) \
      compute[ncompute] = new Class(sparta,narg,arg);
#include "style_compute.h"
#undef ComputeStyle
#undef COMPUTE_CLASS

    else error->all(FLERR,"Unrecognized compute style");
  }

  ncompute++;
}

/* ----------------------------------------------------------------------
   delete a Compute from list of Computes
------------------------------------------------------------------------- */

void Modify::delete_compute(const char *id)
{
  int icompute = find_compute(id);
  if (icompute < 0) error->all(FLERR,"Could not find compute ID to delete");
  delete compute[icompute];

  // move other Computes down in list one slot

  for (int i = icompute+1; i < ncompute; i++) compute[i-1] = compute[i];
  ncompute--;
}

/* ----------------------------------------------------------------------
   find a compute by ID
   return index of compute or -1 if not found
------------------------------------------------------------------------- */

int Modify::find_compute(const char *id)
{
  int icompute;
  for (icompute = 0; icompute < ncompute; icompute++)
    if (strcmp(id,compute[icompute]->id) == 0) break;
  if (icompute == ncompute) return -1;
  return icompute;
}

/* ----------------------------------------------------------------------
   clear invoked flag of all computes
   called everywhere that computes are used, before computes are invoked
   invoked flag used to avoid re-invoking same compute multiple times
   and to flag computes that store invocation times as having been invoked
------------------------------------------------------------------------- */

void Modify::clearstep_compute()
{
  for (int icompute = 0; icompute < ncompute; icompute++)
    compute[icompute]->invoked_flag = 0;
}

/* ----------------------------------------------------------------------
   loop over computes that store invocation times
   if its invoked flag set on this timestep, schedule next invocation
   called everywhere that computes are used, after computes are invoked
------------------------------------------------------------------------- */

void Modify::addstep_compute(bigint newstep)
{
  for (int icompute = 0; icompute < n_timeflag; icompute++)
    if (compute[list_timeflag[icompute]]->invoked_flag)
      compute[list_timeflag[icompute]]->addstep(newstep);
}

/* ----------------------------------------------------------------------
   loop over all computes
   schedule next invocation for those that store invocation times
   called when not sure what computes will be needed on newstep
   do not loop only over n_timeflag, since may not be set yet
------------------------------------------------------------------------- */

void Modify::addstep_compute_all(bigint newstep)
{
  for (int icompute = 0; icompute < ncompute; icompute++)
    if (compute[icompute]->timeflag) compute[icompute]->addstep(newstep);
}

/* ----------------------------------------------------------------------
   create list of fix indices for fixes which match mask
------------------------------------------------------------------------- */

void Modify::list_init(int mask, int &n, int *&list)
{
  delete [] list;

  n = 0;
  for (int i = 0; i < nfix; i++) if (fmask[i] & mask) n++;
  list = new int[n];

  n = 0;
  for (int i = 0; i < nfix; i++) if (fmask[i] & mask) list[n++] = i;
}

/* ----------------------------------------------------------------------
   create list of fix indices for end_of_step fixes
   also create end_of_step_every[]
------------------------------------------------------------------------- */

void Modify::list_init_end_of_step(int mask, int &n, int *&list)
{
  delete [] list;
  delete [] end_of_step_every;

  n = 0;
  for (int i = 0; i < nfix; i++) if (fmask[i] & mask) n++;
  list = new int[n];
  end_of_step_every = new int[n];

  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fmask[i] & mask) {
      list[n] = i;
      end_of_step_every[n++] = fix[i]->nevery;
    }
}

/* ----------------------------------------------------------------------
   create list of indices for fixes with various attributes
------------------------------------------------------------------------- */

void Modify::list_init_fixes()
{
  delete [] list_pergrid;
  delete [] list_update_custom;
  delete [] list_gas_react;
  delete [] list_surf_react;

  n_pergrid = n_update_custom = n_gas_react = n_surf_react = 0;
  for (int i = 0; i < nfix; i++) {
    if (fix[i]->gridmigrate) n_pergrid++;
    if (fix[i]->flag_update_custom) n_update_custom++;
    if (fix[i]->flag_gas_react) n_gas_react++;
    if (fix[i]->flag_surf_react) n_surf_react++;
  }

  list_pergrid = new int[n_pergrid];
  list_update_custom = new int[n_update_custom];
  list_gas_react = new int[n_gas_react];
  list_surf_react = new int[n_surf_react];

  n_pergrid = n_update_custom = n_gas_react = n_surf_react = 0;
  for (int i = 0; i < nfix; i++) {
    if (fix[i]->gridmigrate) list_pergrid[n_pergrid++] = i;
    if (fix[i]->flag_update_custom) list_update_custom[n_update_custom++] = i;
    if (fix[i]->flag_gas_react) list_gas_react[n_gas_react++] = i;
    if (fix[i]->flag_surf_react) list_surf_react[n_surf_react++] = i;
  }
}

/* ----------------------------------------------------------------------
   create list of indices for computes which various attributes
------------------------------------------------------------------------- */

void Modify::list_init_computes()
{
  delete [] list_timeflag;

  n_timeflag = 0;
  for (int i = 0; i < ncompute; i++)
    if (compute[i]->timeflag) n_timeflag++;
  list_timeflag = new int[n_timeflag];

  n_timeflag = 0;
  for (int i = 0; i < ncompute; i++)
    if (compute[i]->timeflag) list_timeflag[n_timeflag++] = i;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory from all fixes
------------------------------------------------------------------------- */

bigint Modify::memory_usage()
{
  bigint bytes = 0;
  for (int i = 0; i < nfix; i++) bytes += fix[i]->memory_usage();
  for (int i = 0; i < ncompute; i++) bytes += compute[i]->memory_usage();
  return bytes;
}

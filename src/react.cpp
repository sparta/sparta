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
#include "string.h"
#include "react.h"
#include "update.h"
#include "grid.h"
#include "comm.h"
#include "input.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "error.h"
#include "modify.h"
#include "memory.h"
#include "compute.h"
#include "fix.h"

using namespace SPARTA_NS;
enum{NONE,COMPUTE,FIX};

#define INVOKED_PER_GRID 16

/* ---------------------------------------------------------------------- */

React::React(SPARTA *sparta, int, char **arg) : Pointers(sparta)
{
  int n = strlen(arg[0]) + 1;
  style = new char[n];
  strcpy(style,arg[0]);

  recombflag_user = 1;
  recomb_boost = 1000.0;
  recomb_boost_inverse = 0.001;
  computeChemRates = 0;
  partialEnergy = 1;
  id_temp = NULL;
  temp = NULL;
  nglocal = 0;

  random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,comm->me,100);

  copy = copymode = 0;
}

/* ---------------------------------------------------------------------- */

React::~React()
{
  if (copy) return;

  delete [] style;
  delete random;

  delete [] id_temp;
  memory->destroy(temp);
}

/* ---------------------------------------------------------------------- */

void React::modify_params(int narg, char **arg)
{
  if (narg == 0) error->all(FLERR,"Illegal react_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"recomb") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal react_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) recombflag_user = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) recombflag_user = 0;
      else error->all(FLERR,"Illegal react_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"rboost") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal react_modify command");
      recomb_boost = input->numeric(FLERR,arg[iarg+1]);
      if (recomb_boost < 1.0) error->all(FLERR,"Illegal react_modify command");
      recomb_boost_inverse = 1.0 / recomb_boost;
      iarg += 2;
    } else if (strcmp(arg[iarg],"compute_chem_rates") == 0) {
        if (iarg+2 > narg) error->all(FLERR,"Illegal react_modify command");
        if (strcmp(arg[iarg+1],"yes") == 0) computeChemRates = 1;
        else if (strcmp(arg[iarg+1],"no") == 0) computeChemRates = 0;
        else error->all(FLERR,"Illegal react_modify command");
        iarg += 2;
    } else if (strcmp(arg[iarg],"partial_energy") == 0) {
      if (strcmp(arg[iarg+1],"yes") == 0) {
          if (iarg+2 > narg) error->all(FLERR,"Illegal react_modify command");
          partialEnergy = 1;
          iarg += 2;
      } else if (strcmp(arg[iarg+1],"no") == 0) {
          if (iarg+3 > narg) error->all(FLERR,"Illegal react_modify command");
          partialEnergy = 0;

          if (strncmp(arg[iarg+2],"c_",2) == 0 || strncmp(arg[iarg+2],"f_",2) == 0) {
            int n = strlen(arg[iarg+2]);
            id_temp = new char[n];
            strcpy(id_temp,&arg[iarg+2][2]);

            char *ptr = strchr(id_temp,'[');
            if (ptr) {
              if (id_temp[strlen(id_temp)-1] != ']')
                error->all(FLERR,"Invalid temp in react_modify partialEnergy yes command");
              tempindex = atoi(ptr+1);
              *ptr = '\0';
            } else tempindex = 0;

            if (strncmp(arg[iarg+2],"c_",2) == 0) tempwhich = COMPUTE;
            else tempwhich = FIX;

            if (tempwhich == COMPUTE) {
              int n = modify->find_compute(id_temp);
              if (n < 0)
                error->all(FLERR,"Could not find react_modify partialEnergy yes compute ID");
              if (modify->compute[n]->per_grid_flag == 0)
                error->all(FLERR,"React_modify partialEnergy yes compute does not "
                           "compute per-grid info");
              if (tempindex == 0 && modify->compute[n]->size_per_grid_cols > 0)
                error->all(FLERR,
                           "React_modify partialEnergy yes compute does not "
                           "compute per-grid vector");
              if (tempindex > 0 && modify->compute[n]->size_per_grid_cols == 0)
                error->all(FLERR,
                           "React_modify partialEnergy yes compute does not "
                           "compute per-grid array");
              if (tempindex > 0 && tempindex > modify->compute[n]->size_per_grid_cols)
                error->all(FLERR,"React_modify partialEnergy yes compute vector is "
                           "accessed out-of-range");
              ctemp = modify->compute[n];
            } else {
              int n = modify->find_fix(id_temp);
              if (n < 0) error->all(FLERR,"Could not find react_modofy partialEnergy yes fix ID");
              if (modify->fix[n]->per_grid_flag == 0)
                error->all(FLERR,"React_modify partialEnergy yes fix does not "
                           "compute per-grid info");
              if (tempindex == 0 && modify->fix[n]->size_per_grid_cols > 0)
                error->all(FLERR,"React_modify partialEnergy yes fix does not "
                           "compute per-grid vector");
              if (tempindex > 0 && modify->fix[n]->size_per_grid_cols == 0)
                error->all(FLERR,"React_modify partialEnergy yes fix does not "
                           "compute per-grid array");
              if (tempindex > 0 && tempindex > modify->fix[n]->size_per_grid_cols)
                error->all(FLERR,"React_modify partialEnergy yes fix array is "
                           "accessed out-of-range");
              ftemp = modify->fix[n];
            }
          } else if (strcmp(arg[iarg+2],"NULL") == 0) tempwhich = NONE;
          else error->all(FLERR,"Illegal react_modofy partialEnergy yes command");
          iarg += 3;
          nglocal = grid->nlocal;
          memory->create(temp,nglocal,"lambda/grid:temp");
      }

    } else error->all(FLERR,"Illegal react_modify command");

  }
}

/* ---------------------------------------------------------------------- */

void React::compute_per_grid()
{
  invoked_per_grid = update->ntimestep;

  if (tempwhich == FIX && update->ntimestep % ftemp->per_grid_freq)
    error->all(FLERR,"React_modify fix not computed at compatible time");

  // grab temp value from compute or fix
  // invoke temp compute as needed

  if (tempwhich == COMPUTE) {
    if (!(ctemp->invoked_flag & INVOKED_PER_GRID)) {
      ctemp->compute_per_grid();
      ctemp->invoked_flag |= INVOKED_PER_GRID;
    }
    if (ctemp->post_process_grid_flag)
      ctemp->post_process_grid(tempindex,1,NULL,NULL,NULL,1);

    if (tempindex == 0 || ctemp->post_process_grid_flag)
      memcpy(temp,ctemp->vector_grid,nglocal*sizeof(double));
    else {
      int index = tempindex-1;
      double **array = ctemp->array_grid;
      for (int i = 0; i < nglocal; i++)
        temp[i] = array[i][index];
    }

  } else if (tempwhich == FIX) {
    if (tempindex == 0) {
      memcpy(temp,ftemp->vector_grid,nglocal*sizeof(double));
    }
    else {
      double **array = ftemp->array_grid;
      int index = tempindex-1;
      for (int i = 0; i < nglocal; i++)
        temp[i] = array[i][index];
    }
  }
}

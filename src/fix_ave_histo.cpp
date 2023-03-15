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

#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "fix_ave_histo.h"
#include "update.h"
#include "particle.h"
#include "mixture.h"
#include "grid.h"
#include "domain.h"
#include "region.h"
#include "modify.h"
#include "compute.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{X,V,F,COMPUTE,FIX,VARIABLE};
enum{ONE,RUNNING};
enum{SCALAR,VECTOR,WINDOW};
enum{GLOBAL,PERPARTICLE,PERGRID};
enum{IGNORE,END,EXTRA};

#define INVOKED_SCALAR 1
#define INVOKED_VECTOR 2
#define INVOKED_ARRAY 4
#define INVOKED_PER_PARTICLE 8
#define INVOKED_PER_GRID 16

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

FixAveHisto::FixAveHisto(SPARTA *spa, int narg, char **arg) :
  Fix(spa, narg, arg)
{
  if (narg < 9) error->all(FLERR,"Illegal fix ave/histo command");

  MPI_Comm_rank(world,&me);
  weightflag = 0;

  nevery = input->inumeric(FLERR,arg[2]);
  nrepeat = input->inumeric(FLERR,arg[3]);
  nfreq = input->inumeric(FLERR,arg[4]);

  time_depend = 1;
  global_freq = nfreq;
  vector_flag = 1;
  size_vector = 4;
  array_flag = 1;
  size_array_cols = 3;

  lo = input->numeric(FLERR,arg[5]);
  hi = input->numeric(FLERR,arg[6]);
  nbins = input->inumeric(FLERR,arg[7]);

  // scan values to count them
  // then read options so know mode = SCALAR/VECTOR before re-reading values

  nvalues = 0;

  int iarg = 8;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"x") == 0 ||
        strcmp(arg[iarg],"y") == 0 ||
        strcmp(arg[iarg],"z") == 0 ||
        strcmp(arg[iarg],"vx") == 0 ||
        strcmp(arg[iarg],"vy") == 0 ||
        strcmp(arg[iarg],"vz") == 0 ||
        strncmp(arg[iarg],"c_",2) == 0 ||
        strncmp(arg[iarg],"f_",2) == 0 ||
        strncmp(arg[iarg],"v_",2) == 0) {
      nvalues++;
      iarg++;
    } else break;
  }

  if (nvalues == 0) error->all(FLERR,"No values in fix ave/histo command");

  options(iarg,narg,arg);

  // expand args if any have wildcard character "*"
  // this can reset nvalues

  int expand = 0;
  char **earg;
  nvalues = input->expand_args(nvalues,&arg[8],mode,earg);

  if (earg != &arg[8]) expand = 1;
  arg = earg;

  // parse values

  which = new int[nvalues];
  argindex = new int[nvalues];
  value2index = new int[nvalues];
  ids = new char*[nvalues];

  for (int i = 0; i < nvalues; i++) {
    if (strcmp(arg[i],"x") == 0) {
      which[i] = X;
      argindex[i] = 0;
      ids[i] = NULL;
    } else if (strcmp(arg[i],"y") == 0) {
      which[i] = X;
      argindex[i] = 1;
      ids[i] = NULL;
    } else if (strcmp(arg[i],"z") == 0) {
      which[i] = X;
      argindex[i] = 2;
      ids[i] = NULL;

    } else if (strcmp(arg[i],"vx") == 0) {
      which[i] = V;
      argindex[i] = 0;
      ids[i] = NULL;
    } else if (strcmp(arg[i],"vy") == 0) {
      which[i] = V;
      argindex[i] = 1;
      ids[i] = NULL;
    } else if (strcmp(arg[i],"vz") == 0) {
      which[i] = V;
      argindex[i] = 2;
      ids[i] = NULL;

    } else if ((strncmp(arg[i],"c_",2) == 0) ||
        (strncmp(arg[i],"f_",2) == 0) ||
        (strncmp(arg[i],"v_",2) == 0)) {
      if (arg[i][0] == 'c') which[i] = COMPUTE;
      else if (arg[i][0] == 'f') which[i] = FIX;
      else if (arg[i][0] == 'v') which[i] = VARIABLE;

      int n = strlen(arg[i]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[i][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Illegal fix ave/histo command");
        argindex[i] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[i] = 0;

      n = strlen(suffix) + 1;
      ids[i] = new char[n];
      strcpy(ids[i],suffix);
      delete [] suffix;
    }
  }

  // if wildcard expansion occurred, free earg memory from expand_args()

  if (expand) {
    for (int i = 0; i < nvalues; i++) delete [] earg[i];
    memory->sfree(earg);
  }

  // setup and error check
  // kind = inputs are all global, or all per-atom, or all local
  // for fix inputs, check that fix frequency is acceptable

  if (nevery <= 0 || nrepeat <= 0 || nfreq <= 0)
    error->all(FLERR,"Illegal fix ave/histo command");
  if (nfreq % nevery || nrepeat*nevery > nfreq)
    error->all(FLERR,"Illegal fix ave/histo command");
  if (lo >= hi) error->all(FLERR,"Illegal fix ave/histo command");
  if (nbins <= 0) error->all(FLERR,"Illegal fix ave/histo command");
  if (ave != RUNNING && overwrite)
    error->all(FLERR,"Illegal fix ave/histo command");

  int kindflag;
  for (int i = 0; i < nvalues; i++) {
    if (which[i] == X || which[i] == V) kindflag = PERPARTICLE;
    else if (which[i] == COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix ave/histo does not exist");
      Compute *compute = modify->compute[icompute];
      if (compute->scalar_flag || compute->vector_flag || compute->array_flag)
        kindflag = GLOBAL;
      else if (compute->per_particle_flag) kindflag = PERPARTICLE;
      else if (compute->per_grid_flag) kindflag = PERGRID;
      else error->all(FLERR,"Fix ave/histo input is invalid compute");
    } else if (which[i] == FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for fix ave/histo does not exist");
      Fix *fix = modify->fix[ifix];
      if (fix->scalar_flag || fix->vector_flag || fix->array_flag)
        kindflag = GLOBAL;
      else if (fix->per_particle_flag) kindflag = PERPARTICLE;
      else if (fix->per_grid_flag) kindflag = PERGRID;
      else error->all(FLERR,"Fix ave/histo input is invalid fix");
    } else if (which[i] == VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for fix ave/histo does not exist");
      if (input->variable->equal_style(ivariable)) kindflag = GLOBAL;
      else if (input->variable->particle_style(ivariable))
        kindflag = PERPARTICLE;
      else if (input->variable->grid_style(ivariable)) kindflag = PERGRID;
      else error->all(FLERR,"Fix ave/histo input is invalid variable");
    }
    if (i == 0) kind = kindflag;
    else if (kindflag != kind)
      error->all(FLERR,"Fix ave/histo inputs are not all "
                 "global, per-particle, or per-grid");
  }

  if (kind == PERPARTICLE && mode == SCALAR)
    error->all(FLERR,
               "Fix ave/histo cannot input per-particle values in scalar mode");
  if (kind == PERGRID && mode == SCALAR)
    error->all(FLERR,
               "Fix ave/histo cannot input per-grid values in scalar mode");

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE && kind == GLOBAL && mode == SCALAR) {
      int icompute = modify->find_compute(ids[i]);
      if (argindex[i] == 0 && modify->compute[icompute]->scalar_flag == 0)
        error->all(FLERR,
                   "Fix ave/histo compute does not calculate a global scalar");
      if (argindex[i] && modify->compute[icompute]->vector_flag == 0)
        error->all(FLERR,
                   "Fix ave/histo compute does not calculate a global vector");
      if (argindex[i] && argindex[i] > modify->compute[icompute]->size_vector)
        error->all(FLERR,
                   "Fix ave/histo compute vector is accessed out-of-range");

    } else if (which[i] == COMPUTE && kind == GLOBAL && mode == VECTOR) {
      int icompute = modify->find_compute(ids[i]);
      if (argindex[i] == 0 && modify->compute[icompute]->vector_flag == 0)
        error->all(FLERR,
                   "Fix ave/histo compute does not calculate a global vector");
      if (argindex[i] && modify->compute[icompute]->array_flag == 0)
        error->all(FLERR,
                   "Fix ave/histo compute does not calculate a global array");
      if (argindex[i] &&
          argindex[i] > modify->compute[icompute]->size_array_cols)
        error->all(FLERR,
                   "Fix ave/histo compute array is accessed out-of-range");

    } else if (which[i] == COMPUTE && kind == PERPARTICLE) {
      int icompute = modify->find_compute(ids[i]);
      if (modify->compute[icompute]->per_particle_flag == 0)
        error->all(FLERR,"Fix ave/histo compute does not "
                   "calculate per-particle values");
      if (argindex[i] == 0 &&
          modify->compute[icompute]->size_per_particle_cols != 0)
        error->all(FLERR,"Fix ave/histo compute does not "
                   "calculate a per-particle vector");
      if (argindex[i] && modify->compute[icompute]->size_per_particle_cols == 0)
        error->all(FLERR,"Fix ave/histo compute does not "
                   "calculate a per-particle array");
      if (argindex[i] &&
          argindex[i] > modify->compute[icompute]->size_per_particle_cols)
        error->all(FLERR,
                   "Fix ave/histo compute array is accessed out-of-range");

    } else if (which[i] == COMPUTE && kind == PERGRID) {
      int icompute = modify->find_compute(ids[i]);
      if (modify->compute[icompute]->per_grid_flag == 0)
        error->all(FLERR,"Fix ave/histo compute does not "
                   "calculate per-grid values");
      if (argindex[i] == 0 &&
          modify->compute[icompute]->size_per_grid_cols != 0)
        error->all(FLERR,"Fix ave/histo compute does not "
                   "calculate a per-grid vector");
      if (argindex[i] && modify->compute[icompute]->size_per_grid_cols == 0)
        error->all(FLERR,"Fix ave/histo compute does not "
                   "calculate a per-grid array");
      if (argindex[i] &&
          argindex[i] > modify->compute[icompute]->size_per_grid_cols)
        error->all(FLERR,
                   "Fix ave/histo compute array is accessed out-of-range");

    } else if (which[i] == FIX && kind == GLOBAL && mode == SCALAR) {
      int ifix = modify->find_fix(ids[i]);
      if (argindex[i] == 0 && modify->fix[ifix]->scalar_flag == 0)
        error->all(FLERR,
                   "Fix ave/histo fix does not calculate a global scalar");
      if (argindex[i] && modify->fix[ifix]->vector_flag == 0)
        error->all(FLERR,
                   "Fix ave/histo fix does not calculate a global vector");
      if (argindex[i] && argindex[i] > modify->fix[ifix]->size_vector)
        error->all(FLERR,"Fix ave/histo fix vector is accessed out-of-range");
      if (nevery % modify->fix[ifix]->global_freq)
        error->all(FLERR,
                   "Fix for fix ave/histo not computed at compatible time");

    } else if (which[i] == FIX && kind == GLOBAL && mode == VECTOR) {
      int ifix = modify->find_fix(ids[i]);
      if (argindex[i] == 0 && modify->fix[ifix]->vector_flag == 0)
        error->all(FLERR,
                   "Fix ave/histo fix does not calculate a global vector");
      if (argindex[i] && modify->fix[ifix]->array_flag == 0)
        error->all(FLERR,"Fix ave/histo fix does not calculate a global array");
      if (argindex[i] && argindex[i] > modify->fix[ifix]->size_array_cols)
        error->all(FLERR,"Fix ave/histo fix array is accessed out-of-range");
      if (nevery % modify->fix[ifix]->global_freq)
        error->all(FLERR,
                   "Fix for fix ave/histo not computed at compatible time");

    } else if (which[i] == FIX && kind == PERPARTICLE) {
      int ifix = modify->find_fix(ids[i]);
      if (modify->fix[ifix]->per_particle_flag == 0)
        error->all(FLERR,
                   "Fix ave/histo fix does not calculate per-particle values");
      if (argindex[i] == 0 &&
          modify->fix[ifix]->size_per_particle_cols != 0)
        error->all(FLERR,"Fix ave/histo fix does not "
                   "calculate a per-particle vector");
      if (argindex[i] && modify->fix[ifix]->size_per_particle_cols == 0)
        error->all(FLERR,"Fix ave/histo fix does not "
                   "calculate a per-particle array");
      if (argindex[i] &&
          argindex[i] > modify->fix[ifix]->size_per_particle_cols)
        error->all(FLERR,"Fix ave/histo fix array is accessed out-of-range");
      if (nevery % modify->fix[ifix]->per_particle_freq)
        error->all(FLERR,
                   "Fix for fix ave/histo not computed at compatible time");

    } else if (which[i] == FIX && kind == PERGRID) {
      int ifix = modify->find_fix(ids[i]);
      if (modify->fix[ifix]->per_grid_flag == 0)
        error->all(FLERR,
                   "Fix ave/histo fix does not calculate per-grid values");
      if (argindex[i] == 0 &&
          modify->fix[ifix]->size_per_grid_cols != 0)
        error->all(FLERR,"Fix ave/histo fix does not "
                   "calculate a per-grid vector");
      if (argindex[i] && modify->fix[ifix]->size_per_grid_cols == 0)
        error->all(FLERR,"Fix ave/histo fix does not "
                   "calculate a per-grid array");
      if (argindex[i] &&
          argindex[i] > modify->fix[ifix]->size_per_grid_cols)
        error->all(FLERR,"Fix ave/histo fix array is accessed out-of-range");
      if (nevery % modify->fix[ifix]->per_grid_freq)
        error->all(FLERR,
                   "Fix for fix ave/histo not computed at compatible time");

    } else if (which[i] == VARIABLE && kind == GLOBAL && mode == SCALAR) {
      int ivariable = input->variable->find(ids[i]);
      if (input->variable->equal_style(ivariable) == 0)
        error->all(FLERR,"Fix ave/histo variable is not equal-style variable");
      if (argindex[i])
        error->all(FLERR,"Fix ave/histo variable cannot have an index");

    } else if (which[i] == VARIABLE && kind == PERPARTICLE) {
      int ivariable = input->variable->find(ids[i]);
      if (argindex[i] == 0 && input->variable->particle_style(ivariable) == 0)
        error->all(FLERR,
                   "Fix ave/histo variable is not particle-style variable");
      if (argindex[i])
        error->all(FLERR,"Fix ave/histo variable cannot be indexed");

    } else if (which[i] == VARIABLE && kind == PERGRID) {
      int ivariable = input->variable->find(ids[i]);
      if (argindex[i] == 0 && input->variable->grid_style(ivariable) == 0)
        error->all(FLERR,
                   "Fix ave/histo variable is not grid-style variable");
      if (argindex[i])
        error->all(FLERR,"Fix ave/histo variable cannot be indexed");
    }
  }

  // print file comment lines

  if (fp && me == 0) {
    clearerr(fp);
    if (title1) fprintf(fp,"%s\n",title1);
    else fprintf(fp,"# Histogrammed data for fix %s\n",id);
    if (title2) fprintf(fp,"%s\n",title2);
    else fprintf(fp,"# TimeStep Number-of-bins "
                 "Total-counts Missing-counts Min-value Max-value\n");
    if (title3) fprintf(fp,"%s\n",title3);
    else fprintf(fp,"# Bin Coord Count Count/Total\n");

    //else if (strcmp(style,"ave/histo") == 0)
    //  fprintf(fp,"# Bin Coord Count Count/Total\n");
    //else fprintf(fp,"# Bin Coord Count WtCount WtCount/Total\n");

    if (ferror(fp))
      error->one(FLERR,"Error writing file header");

    filepos = ftell(fp);
  }

  delete [] title1;
  delete [] title2;
  delete [] title3;

  // allocate and initialize memory for averaging

  if (beyond == EXTRA) nbins += 2;
  size_array_rows = nbins;

  bin = new double[nbins];
  bin_total = new double[nbins];
  bin_all = new double[nbins];
  coord = new double[nbins];

  stats_list = NULL;
  bin_list = NULL;
  vector = NULL;
  maxvector = 0;

  if (ave == WINDOW) {
    memory->create(stats_list,nwindow,4,"ave/histo:stats_list");
    memory->create(bin_list,nwindow,nbins,"ave/histo:bin_list");
  }

  // initializations
  // set coord to bin centers

  if (beyond == EXTRA) {
    binsize = (hi-lo)/(nbins-2);
    bininv = 1.0/binsize;
  } else {
    binsize = (hi-lo)/nbins;
    bininv = 1.0/binsize;
  }

  if (beyond == EXTRA) {
    coord[0] = lo;
    coord[nbins-1] = hi;
    for (int i = 1; i < nbins-1; i++)
      coord[i] = lo + (i-1+0.5)*binsize;
  } else {
    for (int i = 0; i < nbins; i++)
      coord[i] = lo + (i+0.5)*binsize;
  }

  irepeat = 0;
  iwindow = window_limit = 0;

  stats_total[0] = stats_total[1] = stats_total[2] = stats_total[3] = 0.0;
  for (int i = 0; i < nbins; i++) bin_total[i] = 0.0;

  // nvalid = next step on which end_of_step does something
  // add nvalid to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  nvalid = nextvalid();
  modify->addstep_compute_all(nvalid);
}

/* ---------------------------------------------------------------------- */

FixAveHisto::~FixAveHisto()
{
  if (copymode) return;

  delete [] which;
  delete [] argindex;
  delete [] value2index;
  for (int i = 0; i < nvalues; i++) delete [] ids[i];
  delete [] ids;

  if (fp && me == 0) fclose(fp);

  delete [] bin;
  delete [] bin_total;
  delete [] bin_all;
  delete [] coord;
  memory->destroy(stats_list);
  memory->destroy(bin_list);
  memory->destroy(vector);
}

/* ---------------------------------------------------------------------- */

int FixAveHisto::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAveHisto::init()
{
  // set current indices for all computes,fixes,variables

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix ave/histo does not exist");
      value2index[i] = icompute;

    } else if (which[i] == FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for fix ave/histo does not exist");
      value2index[i] = ifix;

    } else if (which[i] == VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for fix ave/histo does not exist");
      value2index[i] = ivariable;
    }
  }
}

/* ----------------------------------------------------------------------
   only does something if nvalid = current timestep
------------------------------------------------------------------------- */

void FixAveHisto::setup()
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveHisto::end_of_step()
{
  int i,j,m;

  // skip if not step which requires doing something

  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;

  // zero if first step

  if (irepeat == 0) {
    stats[0] = stats[1] = 0.0;
    stats[2] = BIG;
    stats[3] = -BIG;
    for (i = 0; i < nbins; i++) bin[i] = 0.0;
  }

  // accumulate results of computes,fixes,variables to local copy
  // compute/fix/variable may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  // for fix ave/histo/weight, nvalues will be 2
  // first calculate weight factors, then histogram single value

  int ncount = nvalues;
  if (weightflag) {
    calculate_weights();
    ncount = 1;
  }

  for (i = 0; i < ncount; i++) {
    m = value2index[i];
    j = argindex[i];

    // invoke compute if not previously invoked

    if (which[i] == COMPUTE) {
      Compute *compute = modify->compute[m];

      if (kind == GLOBAL && mode == SCALAR) {
        if (j == 0) {
          if (!(compute->invoked_flag & INVOKED_SCALAR)) {
            compute->compute_scalar();
            compute->invoked_flag |= INVOKED_SCALAR;
          }
          bin_one(compute->scalar);
        } else {
          if (!(compute->invoked_flag & INVOKED_VECTOR)) {
            compute->compute_vector();
            compute->invoked_flag |= INVOKED_VECTOR;
          }
          bin_one(compute->vector[j-1]);
        }
      } else if (kind == GLOBAL && mode == VECTOR) {
        if (j == 0) {
          if (!(compute->invoked_flag & INVOKED_VECTOR)) {
            compute->compute_vector();
            compute->invoked_flag |= INVOKED_VECTOR;
          }
          bin_vector(compute->size_vector,compute->vector,1);
        } else {
          if (!(compute->invoked_flag & INVOKED_ARRAY)) {
            compute->compute_array();
            compute->invoked_flag |= INVOKED_ARRAY;
          }
          if (compute->array)
            bin_vector(compute->size_array_rows,&compute->array[0][j-1],
                       compute->size_array_cols);
        }

      } else if (kind == PERPARTICLE) {
        if (!(compute->invoked_flag & INVOKED_PER_PARTICLE)) {
          compute->compute_per_particle();
          compute->invoked_flag |= INVOKED_PER_PARTICLE;
        }
        if (j == 0)
          bin_particles(compute->vector_particle,1);
        else if (compute->array_particle)
          bin_particles(&compute->array_particle[0][j-1],
                        compute->size_per_particle_cols);

      } else if (kind == PERGRID) {
        if (!(compute->invoked_flag & INVOKED_PER_GRID)) {
          compute->compute_per_grid();
          compute->invoked_flag |= INVOKED_PER_GRID;
        }

        if (compute->post_process_grid_flag)
          compute->post_process_grid(j,1,NULL,NULL,NULL,1);
        else if (compute->post_process_isurf_grid_flag)
          compute->post_process_isurf_grid();

        if (j == 0 || compute->post_process_grid_flag)
          bin_grid_cells(compute->vector_grid,1);
        else if (compute->array_grid)
          bin_grid_cells(&compute->array_grid[0][j-1],
                         compute->size_per_grid_cols);
      }

    // access fix fields, guaranteed to be ready

    } else if (which[i] == FIX) {

      Fix *fix = modify->fix[m];

      if (kind == GLOBAL && mode == SCALAR) {
        if (j == 0) bin_one(fix->compute_scalar());
        else bin_one(fix->compute_vector(j-1));

      } else if (kind == GLOBAL && mode == VECTOR) {
        if (j == 0) {
          int n = fix->size_vector;
          for (i = 0; i < n; i++) bin_one(fix->compute_vector(i));
        } else {
          int n = fix->size_vector;
          for (i = 0; i < n; i++) bin_one(fix->compute_array(i,j-1));
        }

      } else if (kind == PERPARTICLE) {
        if (j == 0) bin_particles(fix->vector_particle,1);
        else if (fix->array_particle)
          bin_particles(&fix->array_particle[0][j-1],fix->size_per_particle_cols);

      } else if (kind == PERGRID) {
        if (j == 0) bin_grid_cells(fix->vector_grid,1);
        else if (fix->array_grid)
          bin_grid_cells(&fix->array_grid[0][j-1],fix->size_per_grid_cols);
      }

    // evaluate equal-style or particle-style or grid-style variable

    } else if (which[i] == VARIABLE) {
      if (kind == GLOBAL && mode == SCALAR) {
        bin_one(input->variable->compute_equal(m));

      } else if (which[i] == VARIABLE && kind == PERPARTICLE) {
        if (particle->maxlocal > maxvector) {
          memory->destroy(vector);
          maxvector = particle->maxlocal;
          memory->create(vector,maxvector,"ave/histo:vector");
        }
        input->variable->compute_particle(m,vector,1,0);
        bin_particles(vector,1);

      } else if (which[i] == VARIABLE && kind == PERGRID) {
        if (grid->maxlocal > maxvector) {
          memory->destroy(vector);
          maxvector = grid->maxlocal;
          memory->create(vector,maxvector,"ave/histo:vector");
        }
        input->variable->compute_grid(m,vector,1,0);
        bin_grid_cells(vector,1);
      }

    // explicit per-particle attributes

    } else bin_particles(which[i],j);
  }

  // done if irepeat < nrepeat
  // else reset irepeat and nvalid

  irepeat++;
  if (irepeat < nrepeat) {
    nvalid += nevery;
    modify->addstep_compute(nvalid);
    return;
  }

  irepeat = 0;
  nvalid = ntimestep + nfreq - (nrepeat-1)*nevery;
  modify->addstep_compute(nvalid);

  // merge histogram stats across procs if necessary

  if (kind == PERPARTICLE || kind == PERGRID) {
    MPI_Allreduce(stats,stats_all,2,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&stats[2],&stats_all[2],1,MPI_DOUBLE,MPI_MIN,world);
    MPI_Allreduce(&stats[3],&stats_all[3],1,MPI_DOUBLE,MPI_MAX,world);
    MPI_Allreduce(bin,bin_all,nbins,MPI_DOUBLE,MPI_SUM,world);

    stats[0] = stats_all[0];
    stats[1] = stats_all[1];
    stats[2] = stats_all[2];
    stats[3] = stats_all[3];
    for (i = 0; i < nbins; i++) bin[i] = bin_all[i];
  }

  // if ave = ONE, only single Nfreq timestep value is needed
  // if ave = RUNNING, combine with all previous Nfreq timestep values
  // if ave = WINDOW, combine with nwindow most recent Nfreq timestep values

  if (ave == ONE) {
    stats_total[0] = stats[0];
    stats_total[1] = stats[1];
    stats_total[2] = stats[2];
    stats_total[3] = stats[3];
    for (i = 0; i < nbins; i++) bin_total[i] = bin[i];

  } else if (ave == RUNNING) {
    stats_total[0] += stats[0];
    stats_total[1] += stats[1];
    stats_total[2] = MIN(stats_total[2],stats[2]);
    stats_total[3] = MAX(stats_total[3],stats[3]);
    for (i = 0; i < nbins; i++) bin_total[i] += bin[i];

  } else if (ave == WINDOW) {
    stats_total[0] += stats[0];
    if (window_limit) stats_total[0] -= stats_list[iwindow][0];
    stats_list[iwindow][0] = stats[0];
    stats_total[1] += stats[1];
    if (window_limit) stats_total[1] -= stats_list[iwindow][1];
    stats_list[iwindow][1] = stats[1];

    if (window_limit) m = nwindow;
    else m = iwindow+1;

    stats_list[iwindow][2] = stats[2];
    stats_total[2] = stats_list[0][2];
    for (i = 1; i < m; i++)
      stats_total[2] = MIN(stats_total[2],stats_list[i][2]);
    stats_list[iwindow][3] = stats[3];
    stats_total[3] = stats_list[0][3];
    for (i = 1; i < m; i++)
      stats_total[3] = MAX(stats_total[3],stats_list[i][3]);

    for (i = 0; i < nbins; i++) {
      bin_total[i] += bin[i];
      if (window_limit) bin_total[i] -= bin_list[iwindow][i];
      bin_list[iwindow][i] = bin[i];
    }

    iwindow++;
    if (iwindow == nwindow) {
      iwindow = 0;
      window_limit = 1;
    }
  }

  // output result to file

  if (fp && me == 0) {
    clearerr(fp);
    if (overwrite) fseek(fp,filepos,SEEK_SET);
    fprintf(fp,BIGINT_FORMAT " %d %g %g %g %g\n",ntimestep,nbins,
            stats_total[0],stats_total[1],stats_total[2],stats_total[3]);
    if (stats_total[0] != 0.0)
      for (i = 0; i < nbins; i++)
        fprintf(fp,"%d %g %g %g\n",
                i+1,coord[i],bin_total[i],bin_total[i]/stats_total[0]);
    else
      for (i = 0; i < nbins; i++)
        fprintf(fp,"%d %g %g %g\n",i+1,coord[i],0.0,0.0);

    if (ferror(fp))
      error->one(FLERR,"Error writing out histogram data");

    fflush(fp);
    if (overwrite) {
      long fileend = ftell(fp);
      if (fileend > 0) int tmp = ftruncate(fileno(fp),fileend);
    }
  }
}

/* ----------------------------------------------------------------------
   return Ith vector value
------------------------------------------------------------------------- */

double FixAveHisto::compute_vector(int i)
{
  return stats_total[i];
}

/* ----------------------------------------------------------------------
   return I,J array value
------------------------------------------------------------------------- */

double FixAveHisto::compute_array(int i, int j)
{
  if (j == 0) return coord[i];
  else if (j == 1) return bin_total[i];
  else if (stats_total[0] != 0.0) return bin_total[i]/stats_total[0];
  return 0.0;
}

/* ----------------------------------------------------------------------
   bin a single value
------------------------------------------------------------------------- */

void FixAveHisto::bin_one(double value)
{
  stats[2] = MIN(stats[2],value);
  stats[3] = MAX(stats[3],value);

  if (value < lo) {
    if (beyond == IGNORE) {
      stats[1] += 1.0;
      return;
    } else bin[0] += 1.0;
  } else if (value > hi) {
    if (beyond == IGNORE) {
      stats[1] += 1.0;
      return;
    } else bin[nbins-1] += 1.0;
  } else {
    int ibin = static_cast<int> ((value-lo)*bininv);
    ibin = MIN(ibin,nbins-1);
    if (beyond == EXTRA) ibin++;
    bin[ibin] += 1.0;
  }

  stats[0] += 1.0;
}

/* ----------------------------------------------------------------------
   bin a vector of values with stride
------------------------------------------------------------------------- */

void FixAveHisto::bin_vector(int n, double *values, int stride)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    bin_one(values[m]);
    m += stride;
  }
}

/* ----------------------------------------------------------------------
   bin a per-particle attribute
   index is 0,1,2 if attribute is X or V
------------------------------------------------------------------------- */

void FixAveHisto::bin_particles(int attribute, int index)
{
  Particle::OnePart *particles = particle->particles;
  int nlocal = particle->nlocal;

  Region *region;
  if (regionflag) region = domain->regions[iregion];

  if (attribute == X) {
    if (regionflag && mixflag) {
      int *s2g = particle->mixture[imix]->species2group;
      for (int i = 0; i < nlocal; i++) {
        if (region->match(particles[i].x) &&
            s2g[particles[i].ispecies] >= 0) bin_one(particles[i].x[index]);
      }
    } else if (regionflag) {
      for (int i = 0; i < nlocal; i++) {
        if (region->match(particles[i].x)) bin_one(particles[i].x[index]);
      }
    } else if (mixflag) {
      int *s2g = particle->mixture[imix]->species2group;
      for (int i = 0; i < nlocal; i++) {
        if (s2g[particles[i].ispecies] >= 0) bin_one(particles[i].x[index]);
      }
    } else {
      for (int i = 0; i < nlocal; i++)
        bin_one(particles[i].x[index]);
    }

  } else if (attribute == V) {
    if (regionflag && mixflag) {
      int *s2g = particle->mixture[imix]->species2group;
      for (int i = 0; i < nlocal; i++) {
        if (region->match(particles[i].x) &&
            s2g[particles[i].ispecies] < 0) bin_one(particles[i].v[index]);
      }
    } else if (regionflag) {
      for (int i = 0; i < nlocal; i++) {
        if (region->match(particles[i].x)) bin_one(particles[i].v[index]);
      }
    } else if (mixflag) {
      int *s2g = particle->mixture[imix]->species2group;
      for (int i = 0; i < nlocal; i++) {
        if (s2g[particles[i].ispecies] >= 0) bin_one(particles[i].v[index]);
      }
    } else {
      for (int i = 0; i < nlocal; i++)
        bin_one(particles[i].v[index]);
    }
  }
}

/* ----------------------------------------------------------------------
   bin a per-particle vector of values with stride
------------------------------------------------------------------------- */

void FixAveHisto::bin_particles(double *values, int stride)
{
  Particle::OnePart *particles = particle->particles;
  int nlocal = particle->nlocal;

  Region *region;
  if (regionflag) region = domain->regions[iregion];

  int m = 0;

  if (regionflag && mixflag) {
    int *s2g = particle->mixture[imix]->species2group;
    for (int i = 0; i < nlocal; i++) {
      if (region->match(particles[i].x) &&
          s2g[particles[i].ispecies] >= 0) bin_one(values[m]);
      m += stride;
    }
  } else if (regionflag) {
    for (int i = 0; i < nlocal; i++) {
      if (region->match(particles[i].x)) bin_one(values[m]);
      m += stride;
    }
  } else if (mixflag) {
    int *s2g = particle->mixture[imix]->species2group;
    for (int i = 0; i < nlocal; i++) {
      if (s2g[particles[i].ispecies] < 0) bin_one(values[m]);
      m += stride;
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      bin_one(values[m]);
      m += stride;
    }
  }
}

/* ----------------------------------------------------------------------
   bin a per-grid vector of values with stride
------------------------------------------------------------------------- */

void FixAveHisto::bin_grid_cells(double *values, int stride)
{
  Grid::ChildInfo *cinfo = grid->cinfo;
  int nglocal = grid->nlocal;

  int m = 0;

  if (groupflag) {
    for (int i = 0; i < nglocal; i++) {
      if (cinfo[i].mask & groupbit) bin_one(values[m]);
      m += stride;
    }
  } else {
    for (int i = 0; i < nglocal; i++) {
      bin_one(values[m]);
      m += stride;
    }
  }
}

/* ----------------------------------------------------------------------
   parse optional args
------------------------------------------------------------------------- */

void FixAveHisto::options(int iarg, int narg, char **arg)
{
  // option defaults

  fp = NULL;
  ave = ONE;
  startstep = 0;
  mode = SCALAR;
  beyond = IGNORE;
  overwrite = 0;
  regionflag = 0;
  mixflag = 0;
  groupflag = 0;
  title1 = NULL;
  title2 = NULL;
  title3 = NULL;

  // optional args

  while (iarg < narg) {
    if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/histo command");
      if (me == 0) {
        fp = fopen(arg[iarg+1],"w");
        if (fp == NULL) {
          char str[128];
          sprintf(str,"Cannot open fix ave/histo file %s",arg[iarg+1]);
          error->one(FLERR,str);
        }
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"ave") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/histo command");
      if (strcmp(arg[iarg+1],"one") == 0) ave = ONE;
      else if (strcmp(arg[iarg+1],"running") == 0) ave = RUNNING;
      else if (strcmp(arg[iarg+1],"window") == 0) ave = WINDOW;
      else error->all(FLERR,"Illegal fix ave/histo command");
      if (ave == WINDOW) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal fix ave/histo command");
        nwindow = input->inumeric(FLERR,arg[iarg+2]);
        if (nwindow <= 0) error->all(FLERR,"Illegal fix ave/histo command");
      }
      iarg += 2;
      if (ave == WINDOW) iarg++;
    } else if (strcmp(arg[iarg],"start") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/histo command");
      startstep = input->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"mode") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/histo command");
      if (strcmp(arg[iarg+1],"scalar") == 0) mode = SCALAR;
      else if (strcmp(arg[iarg+1],"vector") == 0) mode = VECTOR;
      else error->all(FLERR,"Illegal fix ave/histo command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"beyond") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/histo command");
      if (strcmp(arg[iarg+1],"ignore") == 0) beyond = IGNORE;
      else if (strcmp(arg[iarg+1],"end") == 0) beyond = END;
      else if (strcmp(arg[iarg+1],"extra") == 0) beyond = EXTRA;
      else error->all(FLERR,"Illegal fix ave/histo command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"overwrite") == 0) {
      overwrite = 1;
      iarg += 1;

    } else if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/histo command");
      regionflag = 1;
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Fix ave/histo region ID does not exist");
      iarg += 2;
    } else if (strcmp(arg[iarg],"mix") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/histo command");
      mixflag = 1;
      imix = particle->find_mixture(arg[iarg+1]);
      if (imix == -1)
        error->all(FLERR,"Fix ave/histo mixture ID does not exist");
      iarg += 2;
    } else if (strcmp(arg[iarg],"group") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/histo command");
      groupflag = 1;
      int igroup = grid->find_group(arg[iarg+1]);
      if (igroup == -1)
        error->all(FLERR,"Fix ave/histo group ID does not exist");
      groupbit = grid->bitmask[igroup];
      iarg += 2;

    } else if (strcmp(arg[iarg],"title1") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/histo command");
      delete [] title1;
      int n = strlen(arg[iarg+1]) + 1;
      title1 = new char[n];
      strcpy(title1,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title2") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/histo command");
      delete [] title2;
      int n = strlen(arg[iarg+1]) + 1;
      title2 = new char[n];
      strcpy(title2,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title3") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/histo command");
      delete [] title3;
      int n = strlen(arg[iarg+1]) + 1;
      title3 = new char[n];
      strcpy(title3,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix ave/histo command");
  }
}

/* ----------------------------------------------------------------------
   calculate nvalid = next step on which end_of_step does something
   can be this timestep if multiple of nfreq and nrepeat = 1
   else backup from next multiple of nfreq
   startstep is lower bound on nfreq multiple
------------------------------------------------------------------------- */

bigint FixAveHisto::nextvalid()
{
  bigint nvalid = (update->ntimestep/nfreq)*nfreq + nfreq;
  while (nvalid < startstep) nvalid += nfreq;
  if (nvalid-nfreq == update->ntimestep && nrepeat == 1)
    nvalid = update->ntimestep;
  else
    nvalid -= (nrepeat-1)*nevery;
  if (nvalid < update->ntimestep) nvalid += nfreq;
  return nvalid;
}

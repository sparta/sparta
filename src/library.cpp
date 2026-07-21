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

// C or Fortran style library interface to SPARTA
// customize by adding new SPARTA-specific functions
//
// the function names and semantics mirror the LAMMPS C library
// interface (lammps_* -> sparta_*) where a LAMMPS equivalent exists

#include "spatype.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "stdio.h"
#include "library.h"
#include "sparta.h"
#include "spaexception.h"
#include "input.h"
#include "update.h"
#include "particle.h"
#include "mixture.h"
#include "domain.h"
#include "region.h"
#include "grid.h"
#include "surf.h"
#include "surf_collide.h"
#include "surf_react.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "output.h"
#include "dump.h"
#include "stats.h"
#include "timer.h"
#include "variable.h"
#include "error.h"
#include "version.h"

#include <map>
#include <string>
#include <vector>

using namespace SPARTA_NS;

// ----------------------------------------------------------------------
// utility macros and functions
// ----------------------------------------------------------------------

// wrap library function bodies so that errors thrown by the Error
// class are caught and recorded instead of aborting the process.
// for an error thrown by Error::one() in a parallel run the other
// MPI ranks cannot be unwound, so abort as before

#define BEGIN_CAPTURE try

#define END_CAPTURE(sparta) \
  catch (SpartaAbortException &e) { \
    SPARTA *_spa = (SPARTA *) sparta; \
    _spa->error->set_last_error(e.what(),Error::ERROR_ABORT); \
    int _nprocs; \
    MPI_Comm_size(_spa->world,&_nprocs); \
    if (_nprocs > 1) MPI_Abort(e.get_universe(),1); \
  } \
  catch (SpartaException &e) { \
    ((SPARTA *) sparta)->error->set_last_error(e.what(),Error::ERROR_NORMAL); \
  } \
  catch (std::exception &e) { \
    ((SPARTA *) sparta)->error->set_last_error(e.what(),Error::ERROR_NORMAL); \
  }

/* copy a C string into a fixed-size buffer, always NULL terminated
   return 1 if successful, 0 if not */

static int copy_string(const char *src, char *buffer, int buf_size)
{
  if (!src || !buffer || buf_size <= 0) return 0;
  strncpy(buffer,src,buf_size-1);
  buffer[buf_size-1] = '\0';
  return 1;
}

/* convert version string "24 Sep 2025" to numeric 20250924 */

static int version2num(const char *version)
{
  static const char *months[] = {"Jan","Feb","Mar","Apr","May","Jun",
                                 "Jul","Aug","Sep","Oct","Nov","Dec"};
  int day,year;
  char month[16];
  if (sscanf(version,"%d %15s %d",&day,month,&year) != 3) return 0;
  int imonth = 0;
  for (int i = 0; i < 12; i++)
    if (strncmp(month,months[i],3) == 0) imonth = i+1;
  if (imonth == 0) return 0;
  return 10000*year + 100*imonth + day;
}

/* build table of styles for each style category,
   extracted from the generated style_*.h header files.
   defining the *_CLASS macros makes each per-style header expand to
   just its Style(key,Class) macro line, and the Style macros below
   record the style name, so no class code is instantiated here */

static const std::map<std::string,std::vector<std::string> > &style_map()
{
  static std::map<std::string,std::vector<std::string> > smap;

  if (smap.empty()) {
    {
      std::vector<std::string> &v = smap["command"];
#define COMMAND_CLASS
#define CommandStyle(key,Class) v.push_back(#key);
#include "style_command.h"
#undef CommandStyle
#undef COMMAND_CLASS
    }
    {
      std::vector<std::string> &v = smap["compute"];
#define COMPUTE_CLASS
#define ComputeStyle(key,Class) v.push_back(#key);
#include "style_compute.h"
#undef ComputeStyle
#undef COMPUTE_CLASS
    }
    {
      std::vector<std::string> &v = smap["fix"];
#define FIX_CLASS
#define FixStyle(key,Class) v.push_back(#key);
#include "style_fix.h"
#undef FixStyle
#undef FIX_CLASS
    }
    {
      std::vector<std::string> &v = smap["dump"];
#define DUMP_CLASS
#define DumpStyle(key,Class) v.push_back(#key);
#include "style_dump.h"
#undef DumpStyle
#undef DUMP_CLASS
    }
    {
      std::vector<std::string> &v = smap["region"];
#define REGION_CLASS
#define RegionStyle(key,Class) v.push_back(#key);
#include "style_region.h"
#undef RegionStyle
#undef REGION_CLASS
    }
    {
      std::vector<std::string> &v = smap["collide"];
#define COLLIDE_CLASS
#define CollideStyle(key,Class) v.push_back(#key);
#include "style_collide.h"
#undef CollideStyle
#undef COLLIDE_CLASS
    }
    {
      std::vector<std::string> &v = smap["react"];
#define REACT_CLASS
#define ReactStyle(key,Class) v.push_back(#key);
#include "style_react.h"
#undef ReactStyle
#undef REACT_CLASS
    }
    {
      std::vector<std::string> &v = smap["surf_collide"];
#define SURF_COLLIDE_CLASS
#define SurfCollideStyle(key,Class) v.push_back(#key);
#include "style_surf_collide.h"
#undef SurfCollideStyle
#undef SURF_COLLIDE_CLASS
    }
    {
      std::vector<std::string> &v = smap["surf_react"];
#define SURF_REACT_CLASS
#define SurfReactStyle(key,Class) v.push_back(#key);
#include "style_surf_react.h"
#undef SurfReactStyle
#undef SURF_REACT_CLASS
    }
  }
  return smap;
}

// ----------------------------------------------------------------------
// library initialization and finalization
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   create an instance of SPARTA and return pointer to it
   pass in command-line args and MPI communicator to run on
   also returns the pointer through the last argument if it is not NULL
------------------------------------------------------------------------- */

void *sparta_open(int argc, char **argv, MPI_Comm communicator, void **ptr)
{
  SPARTA *sparta = NULL;

  try {
    sparta = new SPARTA(argc,argv,communicator);
  } catch (SpartaException &) {
    // error message was already printed by the Error class
    sparta = NULL;
  }

  if (ptr) *ptr = (void *) sparta;
  return (void *) sparta;
}

/* ----------------------------------------------------------------------
   create an instance of SPARTA and return pointer to it
   caller doesn't know MPI communicator, so use MPI_COMM_WORLD
   intialize MPI if needed
------------------------------------------------------------------------- */

void *sparta_open_no_mpi(int argc, char **argv, void **ptr)
{
  int flag;
  MPI_Initialized(&flag);

  if (!flag) {
    int argc = 0;
    char **argv = NULL;
    MPI_Init(&argc,&argv);
  }

  MPI_Comm communicator = MPI_COMM_WORLD;
  return sparta_open(argc,argv,communicator,ptr);
}

/* ----------------------------------------------------------------------
   destruct an instance of SPARTA
------------------------------------------------------------------------- */

void sparta_close(void *ptr)
{
  SPARTA *sparta = (SPARTA *) ptr;
  delete sparta;
}

// ----------------------------------------------------------------------
// executing commands
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   process an input script in filename str
------------------------------------------------------------------------- */

void sparta_file(void *ptr, const char *str)
{
  SPARTA *sparta = (SPARTA *) ptr;

  BEGIN_CAPTURE {
    sparta->input->file(str);
  }
  END_CAPTURE(sparta)
}

/* ----------------------------------------------------------------------
   process a single input command in str
------------------------------------------------------------------------- */

char *sparta_command(void *ptr, const char *str)
{
  SPARTA *sparta = (SPARTA *) ptr;
  char *result = NULL;

  BEGIN_CAPTURE {
    result = sparta->input->one(str);
  }
  END_CAPTURE(sparta)

  return result;
}

/* ----------------------------------------------------------------------
   process multiple input commands in str, separated by newlines
   a line ending in "&" is concatenated with the following line
   Input::line_num tracks the line number of the executing command
------------------------------------------------------------------------- */

void sparta_commands_string(void *ptr, const char *str)
{
  SPARTA *sparta = (SPARTA *) ptr;
  if (!str) return;

  BEGIN_CAPTURE {
    sparta->input->line_num = 0;

    int nline = 0;
    std::string buffer;
    const char *p = str;

    while (1) {
      const char *nl = strchr(p,'\n');
      std::string oneline;
      if (nl) oneline.assign(p,nl-p);
      else oneline.assign(p);
      nline++;

      // strip trailing carriage return and whitespace

      while (!oneline.empty() &&
             (oneline.back() == '\r' || oneline.back() == ' ' ||
              oneline.back() == '\t')) oneline.pop_back();

      // line continuation: strip '&' and append the next line

      if (!oneline.empty() && oneline.back() == '&') {
        oneline.pop_back();
        if (buffer.empty()) sparta->input->line_num = nline;
        buffer += oneline;
      } else {
        if (buffer.empty()) sparta->input->line_num = nline;
        buffer += oneline;
        if (!buffer.empty()) sparta->input->one(buffer.c_str());
        buffer.clear();
      }

      if (!nl) break;
      p = nl + 1;
      if (*p == '\0') break;
    }

    // process any trailing continuation line

    if (!buffer.empty()) sparta->input->one(buffer.c_str());
  }
  END_CAPTURE(sparta)
}

/* ----------------------------------------------------------------------
   clean-up function to free memory allocated by lib and returned to caller
------------------------------------------------------------------------- */

void sparta_free(void *ptr)
{
  free(ptr);
}

// ----------------------------------------------------------------------
// data extraction
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   extract a pointer to an internal SPARTA global entity
   name = desired quantity, e.g. dt or fnum, etc
   returns a void pointer to the entity
     which the caller can cast to the proper data type
   returns a NULL if name not listed below
   customize by adding names
------------------------------------------------------------------------- */

void *sparta_extract_global(void *ptr, const char *name)
{
  SPARTA *sparta = (SPARTA *) ptr;

  // these entries do not require an instance,
  // e.g. so a version check can be done before calling sparta_open()

  // optional build-time git provenance (guarded; absence is harmless so a
  // plain non-git build still returns an empty/"(unknown)" string)
#ifndef SPARTA_GIT_COMMIT
#define SPARTA_GIT_COMMIT ""
#endif
#ifndef SPARTA_GIT_BRANCH
#define SPARTA_GIT_BRANCH ""
#endif

  if (strcmp(name,"sparta_version") == 0) return (void *) SPARTA_VERSION;
  if (strcmp(name,"git_branch") == 0) return (void *) SPARTA_GIT_BRANCH;
  if (strcmp(name,"git_commit") == 0) return (void *) SPARTA_GIT_COMMIT;

  if (!sparta) return NULL;

  if (strcmp(name,"dt") == 0) return (void *) &sparta->update->dt;
  if (strcmp(name,"fnum") == 0) return (void *) &sparta->update->fnum;
  if (strcmp(name,"nrho") == 0) return (void *) &sparta->update->nrho;
  if (strcmp(name,"nplocal") == 0) return (void *) &sparta->particle->nlocal;
  if (strcmp(name,"ntimestep") == 0)
    return (void *) &sparta->update->ntimestep;
  if (strcmp(name,"boxlo") == 0) return (void *) sparta->domain->boxlo;
  if (strcmp(name,"boxhi") == 0) return (void *) sparta->domain->boxhi;
  if (strcmp(name,"units") == 0) return (void *) sparta->update->unit_style;
  return NULL;
}

/* ----------------------------------------------------------------------
   extract a pointer to an internal SPARTA compute-based entity
   id = compute ID
   style = 0 for global data, 1 for per particle data,
     2 for per grid data, 3 for per surf data, 4 for per tally data
     NOTE: for per tally data, some columns of double **array
           may be int/bigint ubuf data, up to caller to decode correctly
   type
     for style global: 0 for scalar, 1 for vector, 2 for array
     for style = particle, grid, surf:
       type = 0 if compute produces a vector
       type = 1 to N if compute produces an array, type = which column of array
   for global data, returns a pointer to the
     compute's internal data structure for the entity
     caller should cast it to (double *) for a scalar or vector
     caller should cast it to (double **) for an array
   for per particle or grid or surf or tally data, returns a pointer to the
     compute's internal data structure for the entity
     caller should cast it to (double *) for a vector
     caller should cast it to (double **) for an array
   returns a void pointer to the compute's internal data structure
     for the entity which the caller can cast to the proper data type
   returns a NULL if id is not recognized or style/type not supported
   IMPORTANT: if the compute is not current it will be invoked
     SPARTA cannot easily check here if it is valid to invoke the compute,
     so caller must insure that it is OK
------------------------------------------------------------------------- */

void *sparta_extract_compute(void *ptr, const char *id, int style, int type)
{
  SPARTA *sparta = (SPARTA *) ptr;
  void *result = NULL;

  BEGIN_CAPTURE {

  int icompute = sparta->modify->find_compute(id);
  if (icompute < 0) return NULL;
  Compute *compute = sparta->modify->compute[icompute];

  if (style == 0) {
    if (type == 0) {
      if (!compute->scalar_flag) return NULL;
      if (compute->invoked_scalar != sparta->update->ntimestep)
        compute->compute_scalar();
      result = (void *) &compute->scalar;
    }
    if (type == 1) {
      if (!compute->vector_flag) return NULL;
      if (compute->invoked_vector != sparta->update->ntimestep)
        compute->compute_vector();
      result = (void *) compute->vector;
    }
    if (type == 2) {
      if (!compute->array_flag) return NULL;
      if (compute->invoked_array != sparta->update->ntimestep)
        compute->compute_array();
      result = (void *) compute->array;
    }
  }

  if (style == 1) {
    if (!compute->per_particle_flag) return NULL;
    if (type == 1) {
      if (compute->invoked_per_particle != sparta->update->ntimestep)
        compute->compute_per_particle();
      result = (void *) compute->vector_particle;
    }
    if (type > 1) {
      if (compute->invoked_per_particle != sparta->update->ntimestep)
        compute->compute_per_particle();
      result = (void *) compute->array_particle;
    }
  }

  if (style == 2) {
    if (!compute->per_grid_flag) return NULL;
    if (type == 1) {
      if (compute->invoked_per_grid != sparta->update->ntimestep)
        compute->compute_per_grid();
      if (compute->post_process_grid_flag)
        compute->post_process_grid(0,1,NULL,NULL,NULL,1);
      else if (compute->post_process_isurf_grid_flag)
        compute->post_process_isurf_grid();
      result = (void *) compute->vector_grid;
    }
    if (type > 1) {
      if (compute->invoked_per_grid != sparta->update->ntimestep)
        compute->compute_per_grid();
      if (compute->post_process_grid_flag) {
        compute->post_process_grid(type,1,NULL,NULL,NULL,1);
        result = (void *) compute->vector_grid;
      } else if (compute->post_process_isurf_grid_flag) {
        compute->post_process_isurf_grid();
        result = (void *) compute->array_grid;
      } else result = (void *) compute->array_grid;
    }
  }

  if (style == 3) {
    if (!compute->per_surf_flag) return NULL;
    if (type == 1) {
      if (compute->invoked_per_surf != sparta->update->ntimestep)
        compute->compute_per_surf();
      compute->post_process_surf();
      result = (void *) compute->vector_surf;
    }
    if (type > 1) {
      if (compute->invoked_per_surf != sparta->update->ntimestep)
        compute->compute_per_surf();
      compute->post_process_surf();
      result = (void *) compute->array_surf;
    }
  }

  if (style == 4) {
    if (!compute->per_tally_flag) return NULL;
    if (type == 1) {
      if (compute->invoked_per_tally != sparta->update->ntimestep)
        compute->compute_per_tally();
      result = (void *) compute->vector_tally;
    }
    if (type > 1) {
      if (compute->invoked_per_tally != sparta->update->ntimestep)
        compute->compute_per_tally();
      result = (void *) compute->array_tally;
    }
  }

  }
  END_CAPTURE(sparta)

  return result;
}

/* ----------------------------------------------------------------------
   extract a pointer to an internal SPARTA fix-based entity
   id = fix ID
   style = 0 for global data, 1 for per particle data,
     2 for per grid data, 3 for per surf data
   type
     for style global: 0 for scalar, 1 for vector
       returns a pointer to allocated memory holding the value(s)
       which the caller must free via sparta_free() to avoid a leak
       for a vector the length is the size_vector member of the fix
     for style = particle, grid, surf:
       type = 0 if fix produces a vector
       type = 1 to N if fix produces an array, type = which column of array
       returns a pointer to the fix's internal data structure
   returns a NULL if id is not recognized or style/type not supported
   IMPORTANT: SPARTA cannot easily check here if the fix data is current,
     so caller must insure that it is OK to access it
------------------------------------------------------------------------- */

void *sparta_extract_fix(void *ptr, const char *id, int style, int type)
{
  SPARTA *sparta = (SPARTA *) ptr;
  void *result = NULL;

  BEGIN_CAPTURE {

  int ifix = sparta->modify->find_fix(id);
  if (ifix < 0) return NULL;
  Fix *fix = sparta->modify->fix[ifix];

  if (style == 0) {
    if (type == 0) {
      if (!fix->scalar_flag) return NULL;
      double *dptr = (double *) malloc(sizeof(double));
      *dptr = fix->compute_scalar();
      result = (void *) dptr;
    }
    if (type == 1) {
      if (!fix->vector_flag) return NULL;
      double *vector = (double *) malloc(fix->size_vector*sizeof(double));
      for (int i = 0; i < fix->size_vector; i++)
        vector[i] = fix->compute_vector(i);
      result = (void *) vector;
    }
  }

  if (style == 1) {
    if (!fix->per_particle_flag) return NULL;
    if (type == 0 || type == 1) result = (void *) fix->vector_particle;
    if (type > 1) result = (void *) fix->array_particle;
  }

  if (style == 2) {
    if (!fix->per_grid_flag) return NULL;
    if (type == 0 || type == 1) result = (void *) fix->vector_grid;
    if (type > 1) result = (void *) fix->array_grid;
  }

  if (style == 3) {
    if (!fix->per_surf_flag) return NULL;
    if (type == 0 || type == 1) result = (void *) fix->vector_surf;
    if (type > 1) result = (void *) fix->array_surf;
  }

  }
  END_CAPTURE(sparta)

  return result;
}

/* ----------------------------------------------------------------------
   extract a pointer to an internal SPARTA evaluated variable
   name = variable name, must be equal-style or particle-style variable
   for equal-style variable, returns a pointer to a memory location
     which is allocated by this function
     which the caller can cast to a (double *) which points to the value
   for particle-style variable, returns a pointer to the
     vector of per-atom values on each processor,
     which the caller can cast to a (double *) which points to the values
   returns a NULL if name is not recognized or not equal-style or particle-style
   IMPORTANT: for both equal-style and particle-style variables,
     this function allocates memory to store the variable data in
     so the caller must free this memory to avoid a leak
     e.g. for equal-style variables
       double *dptr = (double *) sparta_extract_variable();
       double value = *dptr;
       sparta_free(dptr);
     e.g. for particle-style variables
       double *vector = (double *) sparta_extract_variable();
       use the vector values
       sparta_free(vector);
   IMPORTANT: SPARTA cannot easily check here when it is valid to evaluate
     the variable or any fixes or computes or stats info it references,
     so caller must insure that it is OK
------------------------------------------------------------------------- */

void *sparta_extract_variable(void *ptr, const char *name)
{
  SPARTA *sparta = (SPARTA *) ptr;
  void *result = NULL;

  BEGIN_CAPTURE {

  int ivar = sparta->input->variable->find((char *) name);
  if (ivar < 0) return NULL;

  if (sparta->input->variable->equal_style(ivar)) {
    double *dptr = (double *) malloc(sizeof(double));
    *dptr = sparta->input->variable->compute_equal(ivar);
    result = (void *) dptr;
  } else if (sparta->input->variable->particle_style(ivar)) {
    int nlocal = sparta->particle->nlocal;
    double *vector = (double *) malloc(nlocal*sizeof(double));
    sparta->input->variable->compute_particle(ivar,vector,1,0);
    result = (void *) vector;
  }

  }
  END_CAPTURE(sparta)

  return result;
}

/* ----------------------------------------------------------------------
   return the style of a variable as a SPARTA_VAR_ constant
   returns -1 if name is not a defined variable
------------------------------------------------------------------------- */

int sparta_extract_variable_datatype(void *ptr, const char *name)
{
  SPARTA *sparta = (SPARTA *) ptr;
  Variable *variable = sparta->input->variable;

  int ivar = variable->find((char *) name);
  if (ivar < 0) return -1;

  if (variable->internal_style(ivar)) return SPARTA_VAR_INTERNAL;
  if (variable->equal_style(ivar)) return SPARTA_VAR_EQUAL;
  if (variable->particle_style(ivar)) return SPARTA_VAR_PARTICLE;
  if (variable->grid_style(ivar)) return SPARTA_VAR_GRID;
  if (variable->surf_style(ivar)) return SPARTA_VAR_SURF;
  return SPARTA_VAR_STRING;
}

/* ----------------------------------------------------------------------
   extract simulation settings and configuration as an integer
   returns -1 if name is not recognized
------------------------------------------------------------------------- */

int sparta_extract_setting(void *ptr, const char *name)
{
  SPARTA *sparta = (SPARTA *) ptr;

  if (strcmp(name,"bigint") == 0) return sizeof(bigint);

  if (strcmp(name,"dimension") == 0) return sparta->domain->dimension;
  if (strcmp(name,"box_exist") == 0) return sparta->domain->box_exist;
  if (strcmp(name,"grid_exist") == 0) return sparta->grid->exist;
  if (strcmp(name,"surf_exist") == 0) return sparta->surf->exist;

  if (strcmp(name,"world_size") == 0) {
    int nprocs;
    MPI_Comm_size(sparta->world,&nprocs);
    return nprocs;
  }
  if (strcmp(name,"world_rank") == 0) {
    int me;
    MPI_Comm_rank(sparta->world,&me);
    return me;
  }

  if (strcmp(name,"nplocal") == 0) return sparta->particle->nlocal;
  if (strcmp(name,"nspecies") == 0) return sparta->particle->nspecies;
  if (strcmp(name,"nmixture") == 0) return sparta->particle->nmixture;
  if (strcmp(name,"ncompute") == 0) return sparta->modify->ncompute;
  if (strcmp(name,"nfix") == 0) return sparta->modify->nfix;
  if (strcmp(name,"ndump") == 0) return sparta->output->ndump;
  if (strcmp(name,"nregion") == 0) return sparta->domain->nregion;
  if (strcmp(name,"nvariable") == 0)
    return sparta->input->variable->nvar_active();
  if (strcmp(name,"ngroup_grid") == 0) return sparta->grid->ngroup;
  if (strcmp(name,"ngroup_surf") == 0) return sparta->surf->ngroup;
  if (strcmp(name,"nsurf") == 0) return (int) sparta->surf->nsurf;
  if (strcmp(name,"nlocal_surf") == 0) return sparta->surf->nlocal;

  if (strcmp(name,"stats_every") == 0) return sparta->output->stats_every;

  return -1;
}

/* ----------------------------------------------------------------------
   evaluate a stats keyword and return its value as a double
   step, np and cpu-related keywords are always available,
   other stats keywords only while a run is in progress
   returns 0.0 if keyword cannot be evaluated
------------------------------------------------------------------------- */

double sparta_get_thermo(void *ptr, const char *name)
{
  SPARTA *sparta = (SPARTA *) ptr;
  double dvalue = 0.0;

  BEGIN_CAPTURE {

  if (strcmp(name,"step") == 0) {
    dvalue = (double) sparta->update->ntimestep;
  } else if (strcmp(name,"np") == 0) {
    dvalue = (double) sparta->particle->nglobal;
  } else if (strcmp(name,"dt") == 0) {
    dvalue = sparta->update->dt;
  } else if (strcmp(name,"cpu") == 0) {
    if (sparta->update->runflag)
      dvalue = sparta->timer->elapsed(TIME_LOOP);
  } else if (strcmp(name,"cpuremain") == 0) {
    if (sparta->update->runflag &&
        sparta->update->ntimestep > sparta->update->firststep) {
      double cpu = sparta->timer->elapsed(TIME_LOOP);
      bigint done = sparta->update->ntimestep - sparta->update->firststep;
      bigint todo = sparta->update->laststep - sparta->update->ntimestep;
      dvalue = cpu*todo/done;
    }
  } else if (strcmp(name,"cpuuse") == 0) {
    dvalue = sparta->timer->cpu_usage();
  } else if (sparta->update->runflag) {
    sparta->output->stats->evaluate_keyword((char *) name,&dvalue);
  }

  }
  END_CAPTURE(sparta)

  return dvalue;
}

/* ----------------------------------------------------------------------
   access the cached data of the last computed stats output
   this cache can be read safely, e.g. from a GUI thread,
   while a simulation is running on another thread
   what = "lock"    acquire the cache mutex, returns NULL
   what = "unlock"  release the cache mutex, returns NULL
   what = "step"    pointer to bigint with timestep of the cached data
   what = "num"     pointer to int with # of cached columns
   what = "setup"   pointer to int, 1 if in setup phase (data invalid)
   what = "line"    pointer to int with input script line of run command
   what = "imagename"  most recently written dump image filename (char *)
   what = "keyword" name of stats column with given index (char *)
   what = "type"    pointer to int with data type of column (SPARTA_INT,
                    SPARTA_DOUBLE, or SPARTA_INT64)
   what = "data"    pointer to value of column with given index
   returns NULL if what is not recognized or index out of range
------------------------------------------------------------------------- */

void *sparta_last_thermo(void *ptr, const char *what, int index)
{
  SPARTA *sparta = (SPARTA *) ptr;
  return sparta->output->stats->last_thermo(what,index);
}

// ----------------------------------------------------------------------
// introspection
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   return SPARTA version as an integer in YYYYMMDD format
------------------------------------------------------------------------- */

int sparta_version(void *)
{
  return version2num(SPARTA_VERSION);
}

/* ----------------------------------------------------------------------
   return the number of styles available in a style category
   category = command, compute, fix, dump, region, collide, react,
     surf_collide, or surf_react
   returns 0 for an unknown category
------------------------------------------------------------------------- */

int sparta_style_count(void *, const char *category)
{
  const std::map<std::string,std::vector<std::string> > &smap = style_map();
  std::map<std::string,std::vector<std::string> >::const_iterator entry =
    smap.find(category);
  if (entry == smap.end()) return 0;
  return (int) entry->second.size();
}

/* ----------------------------------------------------------------------
   copy the name of style number idx in a style category into buffer
   returns 1 if successful, 0 if not
------------------------------------------------------------------------- */

int sparta_style_name(void *, const char *category, int idx,
                      char *buffer, int buf_size)
{
  const std::map<std::string,std::vector<std::string> > &smap = style_map();
  std::map<std::string,std::vector<std::string> >::const_iterator entry =
    smap.find(category);
  if (entry == smap.end()) return 0;
  if (idx < 0 || idx >= (int) entry->second.size()) return 0;
  return copy_string(entry->second[idx].c_str(),buffer,buf_size);
}

/* ----------------------------------------------------------------------
   return the number of defined instances in an ID category
   category = compute, dump, fix, mixture, region, species,
     surf_collide, surf_react, variable, group_grid, or group_surf
   returns -1 for an unknown category
------------------------------------------------------------------------- */

int sparta_id_count(void *ptr, const char *category)
{
  SPARTA *sparta = (SPARTA *) ptr;

  if (strcmp(category,"compute") == 0) return sparta->modify->ncompute;
  if (strcmp(category,"dump") == 0) return sparta->output->ndump;
  if (strcmp(category,"fix") == 0) return sparta->modify->nfix;
  if (strcmp(category,"mixture") == 0) return sparta->particle->nmixture;
  if (strcmp(category,"region") == 0) return sparta->domain->nregion;
  if (strcmp(category,"species") == 0) return sparta->particle->nspecies;
  if (strcmp(category,"surf_collide") == 0) return sparta->surf->nsc;
  if (strcmp(category,"surf_react") == 0) return sparta->surf->nsr;
  if (strcmp(category,"variable") == 0)
    return sparta->input->variable->nvar_active();
  if (strcmp(category,"group_grid") == 0) return sparta->grid->ngroup;
  if (strcmp(category,"group_surf") == 0) return sparta->surf->ngroup;
  return -1;
}

/* ----------------------------------------------------------------------
   copy the name of instance idx in an ID category into buffer
   returns 1 if successful, 0 if not
------------------------------------------------------------------------- */

int sparta_id_name(void *ptr, const char *category, int idx,
                   char *buffer, int buf_size)
{
  SPARTA *sparta = (SPARTA *) ptr;

  if (idx < 0) return 0;

  if (strcmp(category,"compute") == 0) {
    if (idx >= sparta->modify->ncompute) return 0;
    return copy_string(sparta->modify->compute[idx]->id,buffer,buf_size);
  }
  if (strcmp(category,"dump") == 0) {
    if (idx >= sparta->output->ndump) return 0;
    return copy_string(sparta->output->dump[idx]->id,buffer,buf_size);
  }
  if (strcmp(category,"fix") == 0) {
    if (idx >= sparta->modify->nfix) return 0;
    return copy_string(sparta->modify->fix[idx]->id,buffer,buf_size);
  }
  if (strcmp(category,"mixture") == 0) {
    if (idx >= sparta->particle->nmixture) return 0;
    return copy_string(sparta->particle->mixture[idx]->id,buffer,buf_size);
  }
  if (strcmp(category,"region") == 0) {
    if (idx >= sparta->domain->nregion) return 0;
    return copy_string(sparta->domain->regions[idx]->id,buffer,buf_size);
  }
  if (strcmp(category,"species") == 0) {
    if (idx >= sparta->particle->nspecies) return 0;
    return copy_string(sparta->particle->species[idx].id,buffer,buf_size);
  }
  if (strcmp(category,"surf_collide") == 0) {
    if (idx >= sparta->surf->nsc) return 0;
    return copy_string(sparta->surf->sc[idx]->id,buffer,buf_size);
  }
  if (strcmp(category,"surf_react") == 0) {
    if (idx >= sparta->surf->nsr) return 0;
    return copy_string(sparta->surf->sr[idx]->id,buffer,buf_size);
  }
  if (strcmp(category,"variable") == 0) {
    if (idx >= sparta->input->variable->nvar_active()) return 0;
    return copy_string(sparta->input->variable->name(idx),buffer,buf_size);
  }
  if (strcmp(category,"group_grid") == 0) {
    if (idx >= sparta->grid->ngroup) return 0;
    return copy_string(sparta->grid->gnames[idx],buffer,buf_size);
  }
  if (strcmp(category,"group_surf") == 0) {
    if (idx >= sparta->surf->ngroup) return 0;
    return copy_string(sparta->surf->gnames[idx],buffer,buf_size);
  }
  return 0;
}

/* ----------------------------------------------------------------------
   return 1 if a specific ID exists in an ID category, 0 if not
   returns -1 for an unknown category
------------------------------------------------------------------------- */

int sparta_has_id(void *ptr, const char *category, const char *name)
{
  int count = sparta_id_count(ptr,category);
  if (count < 0) return -1;

  char buffer[256];
  for (int i = 0; i < count; i++) {
    if (sparta_id_name(ptr,category,i,buffer,256) &&
        strcmp(buffer,name) == 0) return 1;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   copy a human-readable description of variable number idx into buffer
   (name, style, and definition), matching the LAMMPS lammps_variable_info
   contract so the GUI can display one variable per line
   returns 1 if successful, 0 if not
------------------------------------------------------------------------- */

int sparta_variable_info(void *ptr, int idx, char *buffer, int buf_size)
{
  SPARTA *sparta = (SPARTA *) ptr;

  if (idx < 0 || idx >= sparta->input->variable->nvar_active()) return 0;
  std::string info = sparta->input->variable->get_info(idx);
  return copy_string(info.c_str(),buffer,buf_size);
}

// ----------------------------------------------------------------------
// run control
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   return 1 if a run is in progress, 0 if not
------------------------------------------------------------------------- */

int sparta_is_running(void *ptr)
{
  SPARTA *sparta = (SPARTA *) ptr;
  return sparta->update->runflag;
}

/* ----------------------------------------------------------------------
   trigger a timeout, so that a run stops cleanly at the next timestep
   the timeout setting is restored at the beginning of the next run
------------------------------------------------------------------------- */

void sparta_force_timeout(void *ptr)
{
  SPARTA *sparta = (SPARTA *) ptr;
  sparta->timer->force_timeout();
}

// ----------------------------------------------------------------------
// error handling
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   return 1 if an error occurred in a previous library call, 0 if not
------------------------------------------------------------------------- */

int sparta_has_error(void *ptr)
{
  SPARTA *sparta = (SPARTA *) ptr;
  return (sparta->error->get_last_error_type() != Error::ERROR_NONE) ? 1 : 0;
}

/* ----------------------------------------------------------------------
   copy the last error message into buffer and clear the error status
   returns 1 if the error was recoverable (Error::all),
   2 if it was not (Error::one in a parallel run), 0 if there was no error
------------------------------------------------------------------------- */

int sparta_get_last_error_message(void *ptr, char *buffer, int buf_size)
{
  SPARTA *sparta = (SPARTA *) ptr;
  Error *error = sparta->error;

  if (error->get_last_error_type() == Error::ERROR_NONE) {
    if (buffer && buf_size > 0) buffer[0] = '\0';
    return 0;
  }

  int type = error->get_last_error_type();
  copy_string(error->get_last_error(),buffer,buf_size);
  error->set_last_error(NULL,Error::ERROR_NONE);
  return type;
}

// ----------------------------------------------------------------------
// compile time configuration
// ----------------------------------------------------------------------

/* return 1 if SPARTA was compiled with a real MPI library, 0 if not */

int sparta_config_has_mpi_support()
{
#ifdef MPI_STUBS
  return 0;
#else
  return 1;
#endif
}

/* return 1 if SPARTA was compiled with PNG support, 0 if not */

int sparta_config_has_png_support()
{
#ifdef SPARTA_PNG
  return 1;
#else
  return 0;
#endif
}

/* return 1 if SPARTA was compiled with JPEG support, 0 if not */

int sparta_config_has_jpeg_support()
{
#ifdef SPARTA_JPEG
  return 1;
#else
  return 0;
#endif
}

/* return 1 if SPARTA was compiled with FFmpeg support, 0 if not */

int sparta_config_has_ffmpeg_support()
{
#ifdef SPARTA_FFMPEG
  return 1;
#else
  return 0;
#endif
}

/* return 1 if SPARTA was compiled with gzip support, 0 if not */

int sparta_config_has_gzip_support()
{
#ifdef SPARTA_GZIP
  return 1;
#else
  return 0;
#endif
}

/* return 1 if errors are recoverable via C++ exceptions
   always true since exceptions are unconditionally enabled */

int sparta_config_has_exceptions()
{
  return 1;
}

/* return 1 if SPARTA was compiled with the named package, 0 if not */

int sparta_config_has_package(const char *name)
{
  if (strcmp(name,"KOKKOS") == 0) {
#ifdef SPARTA_KOKKOS
    return 1;
#else
    return 0;
#endif
  }
  if (strcmp(name,"FFT") == 0) {
    // FFT package presence is detected via its compute style
    const std::vector<std::string> &v = style_map().at("compute");
    for (size_t i = 0; i < v.size(); i++)
      if (v[i] == "fft/grid") return 1;
    return 0;
  }
  return 0;
}

/* return availability of an accelerator package setting
   package = accelerator package name (only KOKKOS)
   category/setting = api:serial/openmp/cuda etc.
   simplified relative to LAMMPS: returns 1 if the package is available
   and the category/setting combination is plausible, 0 if not */

int sparta_config_accelerator(const char *package, const char *category,
                              const char *setting)
{
#ifdef SPARTA_KOKKOS
  if (strcmp(package,"KOKKOS") == 0) {
    if (strcmp(category,"api") == 0) {
      if (strcmp(setting,"serial") == 0) return 1;
#ifdef KOKKOS_ENABLE_OPENMP
      if (strcmp(setting,"openmp") == 0) return 1;
#endif
#ifdef KOKKOS_ENABLE_CUDA
      if (strcmp(setting,"cuda") == 0) return 1;
#endif
#ifdef KOKKOS_ENABLE_HIP
      if (strcmp(setting,"hip") == 0) return 1;
#endif
      return 0;
    }
    if (strcmp(category,"precision") == 0) {
      if (strcmp(setting,"double") == 0) return 1;
      return 0;
    }
    return 0;
  }
#endif
  (void) package;
  (void) category;
  (void) setting;
  return 0;
}

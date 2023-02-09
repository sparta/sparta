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

// C or Fortran style library interface to SPARTA
// customize by adding new SPARTA-specific functions

#include "spatype.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "library.h"
#include "sparta.h"
#include "input.h"
#include "update.h"
#include "particle.h"
#include "modify.h"
#include "compute.h"
#include "variable.h"

using namespace SPARTA_NS;

/* ----------------------------------------------------------------------
   create an instance of SPARTA and return pointer to it
   pass in command-line args and MPI communicator to run on
------------------------------------------------------------------------- */

void sparta_open(int argc, char **argv, MPI_Comm communicator, void **ptr)
{
  SPARTA *sparta = new SPARTA(argc,argv,communicator);
  *ptr = (void *) sparta;
}

/* ----------------------------------------------------------------------
   create an instance of SPARTA and return pointer to it
   caller doesn't know MPI communicator, so use MPI_COMM_WORLD
   intialize MPI if needed
------------------------------------------------------------------------- */

void sparta_open_no_mpi(int argc, char **argv, void **ptr)
{
  int flag;
  MPI_Initialized(&flag);

  if (!flag) {
    int argc = 0;
    char **argv = NULL;
    MPI_Init(&argc,&argv);
  }

  MPI_Comm communicator = MPI_COMM_WORLD;

  SPARTA *sparta = new SPARTA(argc,argv,communicator);
  *ptr = (void *) sparta;
}

/* ----------------------------------------------------------------------
   destruct an instance of SPARTA
------------------------------------------------------------------------- */

void sparta_close(void *ptr)
{
  SPARTA *sparta = (SPARTA *) ptr;
  delete sparta;
}

/* ----------------------------------------------------------------------
   process an input script in filename str
------------------------------------------------------------------------- */

void sparta_file(void *ptr, char *str)
{
  SPARTA *sparta = (SPARTA *) ptr;
  sparta->input->file(str);
}

/* ----------------------------------------------------------------------
   process a single input command in str
------------------------------------------------------------------------- */

char *sparta_command(void *ptr, char *str)
{
  SPARTA *sparta = (SPARTA *) ptr;
  return sparta->input->one(str);
}

/* ----------------------------------------------------------------------
   clean-up function to free memory allocated by lib and returned to caller
------------------------------------------------------------------------- */

void sparta_free(void *ptr)
{
  free(ptr);
}

/* ----------------------------------------------------------------------
   add SPARTA-specific library functions
   all must receive SPARTA pointer as argument
   customize by adding a function here and in library.h header file
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   extract a pointer to an internal SPARTA global entity
   name = desired quantity, e.g. dt or fnum, etc
   returns a void pointer to the entity
     which the caller can cast to the proper data type
   returns a NULL if name not listed below
   customize by adding names
------------------------------------------------------------------------- */

void *sparta_extract_global(void *ptr, char *name)
{
  SPARTA *sparta = (SPARTA *) ptr;

  if (strcmp(name,"dt") == 0) return (void *) &sparta->update->dt;
  if (strcmp(name,"fnum") == 0) return (void *) &sparta->update->fnum;
  if (strcmp(name,"nrho") == 0) return (void *) &sparta->update->nrho;
  if (strcmp(name,"nplocal") == 0) return (void *) &sparta->particle->nlocal;
  return NULL;
}

/* ----------------------------------------------------------------------
   extract a pointer to an internal SPARTA compute-based entity
   id = compute ID
   style = 0 for global data, 1 for per particle data,
     2 for per grid data, 3 for per surf data
   type
     for style global: 0 for scalar, 1 for vector, 2 for array
     for style = particle, grid, surf:
       type = 0 if compute produces a vector
       type = 1 to N if compute produces an array, type = which column of array
   for global data, returns a pointer to the
     compute's internal data structure for the entity
     caller should cast it to (double *) for a scalar or vector
     caller should cast it to (double **) for an array
   for per particle or grid or surf data, returns a pointer to the
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

void *sparta_extract_compute(void *ptr, char *id, int style, int type)
{
  SPARTA *sparta = (SPARTA *) ptr;

  int icompute = sparta->modify->find_compute(id);
  if (icompute < 0) return NULL;
  Compute *compute = sparta->modify->compute[icompute];

  if (style == 0) {
    if (type == 0) {
      if (!compute->scalar_flag) return NULL;
      if (compute->invoked_scalar != sparta->update->ntimestep)
        compute->compute_scalar();
      return (void *) &compute->scalar;
    }
    if (type == 1) {
      if (!compute->vector_flag) return NULL;
      if (compute->invoked_vector != sparta->update->ntimestep)
        compute->compute_vector();
      return (void *) compute->vector;
    }
    if (type == 2) {
      if (!compute->array_flag) return NULL;
      if (compute->invoked_array != sparta->update->ntimestep)
        compute->compute_array();
      return (void *) compute->array;
    }
  }

  if (style == 1) {
    if (!compute->per_particle_flag) return NULL;
    if (type == 1) {
      if (compute->invoked_per_particle != sparta->update->ntimestep)
        compute->compute_per_particle();
      return (void *) compute->vector_particle;
    }
    if (type > 1) {
      if (compute->invoked_per_particle != sparta->update->ntimestep)
        compute->compute_per_particle();
      return (void *) compute->array_particle;
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
      return (void *) compute->vector_grid;
    }
    if (type > 1) {
      if (compute->invoked_per_grid != sparta->update->ntimestep)
        compute->compute_per_grid();
      if (compute->post_process_grid_flag) {
        compute->post_process_grid(type,1,NULL,NULL,NULL,1);
        return (void *) compute->vector_grid;
      }
      if (compute->post_process_isurf_grid_flag) {
        compute->post_process_isurf_grid();
        return (void *) compute->array_grid;
      }
      return (void *) compute->array_grid;
    }
  }

  if (style == 3) {
    if (!compute->per_surf_flag) return NULL;
    if (type == 1) {
      if (compute->invoked_per_surf != sparta->update->ntimestep)
        compute->compute_per_surf();
      compute->post_process_surf();
      return (void *) compute->vector_surf;
    }
    if (type > 1) {
      if (compute->invoked_per_surf != sparta->update->ntimestep)
        compute->compute_per_surf();
      compute->post_process_surf();
      return (void *) compute->array_surf;
    }
  }

  return NULL;
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

void *sparta_extract_variable(void *ptr, char *name)
{
  SPARTA *sparta = (SPARTA *) ptr;

  int ivar = sparta->input->variable->find(name);
  if (ivar < 0) return NULL;

  if (sparta->input->variable->equal_style(ivar)) {
    double *dptr = (double *) malloc(sizeof(double));
    *dptr = sparta->input->variable->compute_equal(ivar);
    return (void *) dptr;
  }

  if (sparta->input->variable->particle_style(ivar)) {
    int nlocal = sparta->particle->nlocal;
    double *vector = (double *) malloc(nlocal*sizeof(double));
    sparta->input->variable->compute_particle(ivar,vector,1,0);
    return (void *) vector;
  }

  return NULL;
}

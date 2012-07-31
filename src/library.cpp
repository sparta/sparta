/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   www.sandia.gov/sparta.html
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2012) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

// C or Fortran style library interface to SPARTA
// customize by adding new SPARTA-specific functions

#include "sptype.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "library.h"
#include "sparta.h"
#include "input.h"
#include "input.h"
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
   name = desired quantity, e.g. ??? or ??? or ???
   returns a void pointer to the entity
     which the caller can cast to the proper data type
   returns a NULL if name not listed below
   customize by adding names
------------------------------------------------------------------------- */

void *sparta_extract_global(void *ptr, char *name)
{
  SPARTA *sparta = (SPARTA *) ptr;

  //if (strcmp(name,"dt") == 0) return (void *) &sparta->update->dt;
  return NULL;
}

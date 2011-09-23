/* ----------------------------------------------------------------------
   DSMC - Sandia parallel DSMC code
   www.sandia.gov/~sjplimp/dsmc.html
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2011) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level DSMC directory.
------------------------------------------------------------------------- */

// C or Fortran style library interface to DSMC
// customize by adding new DSMC-specific functions

#include "dsmctype.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "library.h"
#include "dsmc.h"
#include "input.h"
#include "input.h"
#include "variable.h"

using namespace DSMC_NS;

/* ----------------------------------------------------------------------
   create an instance of DSMC and return pointer to it
   pass in command-line args and MPI communicator to run on
------------------------------------------------------------------------- */

void dsmc_open(int argc, char **argv, MPI_Comm communicator, void **ptr)
{
  DSMC *dsmc = new DSMC(argc,argv,communicator);
  *ptr = (void *) dsmc;
}

/* ----------------------------------------------------------------------
   create an instance of DSMC and return pointer to it
   caller doesn't know MPI communicator, so use MPI_COMM_WORLD
   intialize MPI if needed
------------------------------------------------------------------------- */

void dsmc_open_no_mpi(int argc, char **argv, void **ptr)
{
  int flag;
  MPI_Initialized(&flag);

  if (!flag) {
    int argc = 0;
    char **argv = NULL;
    MPI_Init(&argc,&argv);
  }

  MPI_Comm communicator = MPI_COMM_WORLD;

  DSMC *dsmc = new DSMC(argc,argv,communicator);
  *ptr = (void *) dsmc;
}

/* ----------------------------------------------------------------------
   destruct an instance of DSMC
------------------------------------------------------------------------- */

void dsmc_close(void *ptr)
{
  DSMC *dsmc = (DSMC *) ptr;
  delete dsmc;
}

/* ----------------------------------------------------------------------
   process an input script in filename str
------------------------------------------------------------------------- */

void dsmc_file(void *ptr, char *str)
{
  DSMC *dsmc = (DSMC *) ptr;
  dsmc->input->file(str);
}

/* ----------------------------------------------------------------------
   process a single input command in str
------------------------------------------------------------------------- */

char *dsmc_command(void *ptr, char *str)
{
  DSMC *dsmc = (DSMC *) ptr;
  return dsmc->input->one(str);
}

/* ----------------------------------------------------------------------
   clean-up function to free memory allocated by lib and returned to caller
------------------------------------------------------------------------- */

void dsmc_free(void *ptr)
{
  free(ptr);
}

/* ----------------------------------------------------------------------
   add DSMC-specific library functions
   all must receive DSMC pointer as argument
   customize by adding a function here and in library.h header file
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   extract a pointer to an internal DSMC global entity
   name = desired quantity, e.g. ??? or ??? or ???
   returns a void pointer to the entity
     which the caller can cast to the proper data type
   returns a NULL if name not listed below
   customize by adding names
------------------------------------------------------------------------- */

void *dsmc_extract_global(void *ptr, char *name)
{
  DSMC *dsmc = (DSMC *) ptr;

  //if (strcmp(name,"dt") == 0) return (void *) &dsmc->update->dt;
  return NULL;
}

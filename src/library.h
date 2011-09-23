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

/* 
   C or Fortran style library interface to DSMC
   new DSMC-specific functions can be added
*/

#include "mpi.h"

/* ifdefs allow this file to be included in a C program */

#ifdef __cplusplus
extern "C" {
#endif

void dsmc_open(int, char **, MPI_Comm, void **);
void dsmc_open_no_mpi(int, char **, void **);
void dsmc_close(void *);
void dsmc_file(void *, char *);
char *dsmc_command(void *, char *);
void dsmc_free(void *);

void *dsmc_extract_global(void *, char *);

#ifdef __cplusplus
}
#endif

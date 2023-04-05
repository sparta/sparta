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

/*
   C or Fortran style library interface to SPARTA
   new SPARTA-specific functions can be added
*/

#include "mpi.h"

/* ifdefs allow this file to be included in a C program */

#ifdef __cplusplus
extern "C" {
#endif

void sparta_open(int, char **, MPI_Comm, void **);
void sparta_open_no_mpi(int, char **, void **);
void sparta_close(void *);
void sparta_file(void *, char *);
char *sparta_command(void *, char *);
void sparta_free(void *);

void *sparta_extract_global(void *, char *);
void *sparta_extract_compute(void *, char *, int, int);
void *sparta_extract_variable(void *, char *);

#ifdef __cplusplus
}
#endif

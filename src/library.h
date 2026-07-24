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

/*
   C or Fortran style library interface to SPARTA
   new SPARTA-specific functions can be added

   the function names and semantics mirror the LAMMPS C library
   interface (lammps_* -> sparta_*) where a LAMMPS equivalent exists,
   so that tools written for liblammps (e.g. LAMMPS-GUI) can be
   adapted to SPARTA with minimal changes
*/

#include "mpi.h"

/* data type constants for extracted or returned data */

enum _SPARTA_DATATYPE {
  SPARTA_NONE = -1,     /* no data type assigned (yet) */
  SPARTA_INT = 0,       /* 32-bit integer (array) */
  SPARTA_INT_2D = 1,    /* two-dimensional 32-bit integer array */
  SPARTA_DOUBLE = 2,    /* 64-bit double (array) */
  SPARTA_DOUBLE_2D = 3, /* two-dimensional 64-bit double array */
  SPARTA_INT64 = 4,     /* 64-bit integer (array) */
  SPARTA_INT64_2D = 5,  /* two-dimensional 64-bit integer array */
  SPARTA_STRING = 6     /* C-String */
};

/* style constants for variable data */

enum _SPARTA_VAR_STYLE {
  SPARTA_VAR_EQUAL = 0,    /* equal-style (and compatible) variables */
  SPARTA_VAR_PARTICLE = 1, /* particle-style variables */
  SPARTA_VAR_GRID = 2,     /* grid-style variables */
  SPARTA_VAR_SURF = 3,     /* surf-style variables */
  SPARTA_VAR_STRING = 4,   /* string-valued variables, e.g. index or loop */
  SPARTA_VAR_INTERNAL = 5  /* internal variables */
};

/* ifdefs allow this file to be included in a C program */

#ifdef __cplusplus
extern "C" {
#endif

/* ----------------------------------------------------------------------
 * library initialization and finalization
 * ---------------------------------------------------------------------- */

void *sparta_open(int, char **, MPI_Comm, void **);
void *sparta_open_no_mpi(int, char **, void **);
void sparta_close(void *);

/* ----------------------------------------------------------------------
 * executing commands
 * ---------------------------------------------------------------------- */

void sparta_file(void *, const char *);
char *sparta_command(void *, const char *);
void sparta_commands_string(void *, const char *);

/* ----------------------------------------------------------------------
 * data extraction
 * ---------------------------------------------------------------------- */

void *sparta_extract_global(void *, const char *);
void *sparta_extract_compute(void *, const char *, int, int);
void *sparta_extract_fix(void *, const char *, int, int);
void *sparta_extract_variable(void *, const char *);
int sparta_extract_variable_datatype(void *, const char *);
int sparta_extract_setting(void *, const char *);
double sparta_get_thermo(void *, const char *);
void *sparta_last_thermo(void *, const char *, int);
void sparta_free(void *);

/* ----------------------------------------------------------------------
 * introspection
 * ---------------------------------------------------------------------- */

int sparta_version(void *);
int sparta_style_count(void *, const char *);
int sparta_style_name(void *, const char *, int, char *, int);
int sparta_id_count(void *, const char *);
int sparta_id_name(void *, const char *, int, char *, int);
int sparta_has_id(void *, const char *, const char *);
int sparta_variable_info(void *, int, char *, int);

/* ----------------------------------------------------------------------
 * run control
 * ---------------------------------------------------------------------- */

int sparta_is_running(void *);
void sparta_force_timeout(void *);

/* ----------------------------------------------------------------------
 * error handling
 * ---------------------------------------------------------------------- */

int sparta_has_error(void *);
int sparta_get_last_error_message(void *, char *, int);

/* ----------------------------------------------------------------------
 * compile time configuration
 * ---------------------------------------------------------------------- */

int sparta_config_has_mpi_support();
int sparta_config_has_png_support();
int sparta_config_has_jpeg_support();
int sparta_config_has_ffmpeg_support();
int sparta_config_has_gzip_support();
int sparta_config_has_exceptions();
int sparta_config_has_package(const char *);
int sparta_config_accelerator(const char *, const char *, const char *);

#ifdef __cplusplus
}
#endif

/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2012) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

// define integer data types used by SPARTA and associated size limits

// smallint = variables for on-procesor system (nlocal, nmax, etc)
// bigint = variables for total system (natoms, ntimestep, etc)

// smallint must be an int, as defined by C compiler
// bigint can be 32-bit or 64-bit int, must be >= smallint

// MPI_SPARTA_BIGINT = MPI data type corresponding to a bigint

#ifndef SPARTA_SPTYPE_H
#define SPARTA_SPTYPE_H

#define __STDC_LIMIT_MACROS
#define __STDC_FORMAT_MACROS

#include "limits.h"
#include "stdint.h"
#include "inttypes.h"

// grrr - IBM Power6 does not provide this def in their system header files

#ifndef PRId64
#define PRId64 "ld"
#endif

namespace SPARTA_NS {

// default to 32-bit smallint, 64-bit bigint

#if !defined(SPARTA_SMALL)
#define SPARTA_BIG
#endif

// allow user override of LONGLONG to LONG, necessary for some machines/MPI

#ifdef SPARTA_LONGLONG_TO_LONG
#define MPI_LL MPI_LONG
#define ATOLL atoll
#else
#define MPI_LL MPI_LONG_LONG
#define ATOLL atol
#endif

// for problems that exceed 2 billion (2^31) particles
// 32-bit smallint, 64-bit bigint

#ifdef SPARTA_BIG

typedef int smallint;
typedef int64_t bigint;

#define MAXSMALLINT INT_MAX
#define MAXBIGINT INT64_MAX
#define MPI_SPARTA_BIGINT MPI_LL
#define BIGINT_FORMAT "%" PRId64
#define ATOBIGINT ATOLL

#endif

// for machines that do not support 64-bit ints
// 32-bit smallint and bigint

#ifdef SPARTA_SMALL

typedef int smallint;
typedef int bigint;

#define MAXSMALLINT INT_MAX
#define MAXBIGINT INT_MAX
#define MPI_LMP_BIGINT MPI_INT
#define BIGINT_FORMAT "%d"
#define ATOBIGINT atoi

#endif

}

// settings to enable SPARTA to build under Windows

#ifdef _WIN32
#include "spartawindows.h"
#endif

#endif

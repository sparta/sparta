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

// define integer data types used by SPARTA and associated size limits

// smallint = variables for on-procesor system (nlocal, nmax, etc)
// bigint = variables for total system (natoms, ntimestep, etc)
// cellint = variables for cell IDs

// smallint must be an int, as defined by C compiler
// bigint can be 32-bit or 64-bit int, must be >= smallint
// cellint can be unsigned 32-bit or 64-bit int

// MPI_SPARTA_BIGINT = MPI data type corresponding to a bigint

#ifndef SPARTA_SPTYPE_H
#define SPARTA_SPTYPE_H

#define __STDC_LIMIT_MACROS
#define __STDC_FORMAT_MACROS

#include "limits.h"
#include "stdint.h"
#include "inttypes.h"
#include "accelerator_kokkos_defs.h"

// grrr - IBM Power6 does not provide this def in their system header files

#ifndef PRId64
#define PRId64 "ld"
#endif

namespace SPARTA_NS {

// enum used for KOKKOS host/device flags

enum ExecutionSpace{Host,Device};

// struct alignment for GPUs

#if defined(SPARTA_KOKKOS_GPU)
#define SPARTA_ALIGN(n) alignas(n)
#define SPARTA_GET_ALIGN(type) alignof(type)
#else
#define SPARTA_ALIGN(n)
#define SPARTA_GET_ALIGN(type) 0
#endif

// default settings: 32-bit smallint, 64-bit bigint, 32-bit cellint

#if !defined(SPARTA_SMALL) && !defined(SPARTA_BIG) && !defined(SPARTA_BIGBIG)
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

// default, sufficient for problems with up to 2B grid cells
// 32-bit smallint, 64-bit bigint, 32-bit cellint

#ifdef SPARTA_BIG

typedef int smallint;
typedef uint32_t cellint;
typedef int surfint;
typedef int64_t bigint;

#define MAXSMALLINT INT_MAX
#define MAXBIGINT INT64_MAX
#define MPI_SPARTA_BIGINT MPI_LL
#define CELLINT_FORMAT "%u"
#define SURFINT_FORMAT "%d"
#define BIGINT_FORMAT "%" PRId64
#define ATOCELLINT atoi
#define ATOSURFINT atoi
#define ATOBIGINT ATOLL

#endif

// for problems with more than 2B grid cells
// 32-bit smallint, 64-bit bigint, 64-bit cellint

#ifdef SPARTA_BIGBIG

typedef int smallint;
typedef uint64_t cellint;
typedef int64_t surfint;
typedef int64_t bigint;

#define MAXSMALLINT INT_MAX
#define MAXBIGINT INT64_MAX
#define MPI_SPARTA_BIGINT MPI_LL
#define CELLINT_FORMAT "%" PRIu64
#define SURFINT_FORMAT "%" PRId64
#define BIGINT_FORMAT "%" PRId64
#define ATOCELLINT ATOLL
#define ATOSURFINT ATOLL
#define ATOBIGINT ATOLL

#endif

// for machines that do not support 64-bit ints
// 32-bit smallint and bigint and cellint

#ifdef SPARTA_SMALL

typedef int smallint;
typedef uint32_t cellint;
typedef int surfint;
typedef int bigint;

#define MAXSMALLINT INT_MAX
#define MAXBIGINT INT_MAX
#define MPI_SPARTA_BIGINT MPI_INT
#define CELLINT_FORMAT "%u"
#define SURFINT_FORMAT "%d"
#define BIGINT_FORMAT "%d"
#define ATOCELLINT atoi
#define ATOSURFINT atoi
#define ATOBIGINT atoi

#endif

}

// settings to enable SPARTA to build under Windows

#ifdef _WIN32
#include "spawindows.h"
#endif

#endif

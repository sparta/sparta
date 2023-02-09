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

#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "ctype.h"
#include "read_surf_kokkos.h"
#include "math_extra.h"
#include "surf_kokkos.h"
#include "domain.h"
#include "grid.h"
#include "comm.h"
#include "geometry.h"
#include "input.h"
#include "write_surf.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "error.h"
#include "kokkos_type.h"
#include "sparta_masks.h"
#include "kokkos.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{NEITHER,BAD,GOOD};
enum{NONE,CHECK,KEEP};
enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};           // several files

#define MAXLINE 256
#define CHUNK 1024
#define EPSILON_NORM 1.0e-12
#define EPSILON_GRID 1.0e-3
#define BIG 1.0e20
#define DELTA 128           // must be 2 or greater

/* ---------------------------------------------------------------------- */

ReadSurfKokkos::ReadSurfKokkos(SPARTA *sparta) : ReadSurf(sparta)
{

}

/* ---------------------------------------------------------------------- */

void ReadSurfKokkos::command(int narg, char **arg)
{
  ReadSurf::command(narg,arg);

  SurfKokkos* surf_kk = (SurfKokkos*) surf;
  surf_kk->modify(Host,ALL_MASK);
}

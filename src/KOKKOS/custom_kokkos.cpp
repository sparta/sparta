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

#include "stdlib.h"
#include "string.h"
#include "custom.h"
#include "domain.h"
#include "comm.h"
#include "particle.h"
#include "grid.h"
#include "surf.h"
#include "mixture.h"
#include "region.h"
#include "input.h"
#include "variable.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{SET,REMOVE};
enum{EQUAL,PARTICLE,GRID,SURF};
enum{INT,DOUBLE};                       // several files

/* ---------------------------------------------------------------------- */

CustomKokkos::CustomKokkos(SPARTA *sparta) : Custom(sparta) {}

/* ---------------------------------------------------------------------- */

int CustomKokkos::set_particle(double scalar, double *vector)
{

}

/* ---------------------------------------------------------------------- */

int CustomKokkos::set_grid(double scalar, double *vector)
{

}

/* ---------------------------------------------------------------------- */

int CustomKokkos::set_surf(double scalar, double *vector)
{

}

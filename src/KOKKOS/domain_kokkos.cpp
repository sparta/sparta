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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "domain_kokkos.h"
#include "style_region.h"
#include "update.h"
#include "grid.h"
#include "comm.h"
#include "particle.h"
#include "region.h"
#include "surf.h"
#include "surf_collide.h"
#include "comm.h"
#include "memory_kokkos.h"
#include "error.h"

using namespace SPARTA_NS;

enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};         // several files
enum{PERIODIC,OUTFLOW,REFLECT,SURFACE,AXISYM};  // several files

#define DELTAREGION 4

/* ---------------------------------------------------------------------- */

DomainKokkos::DomainKokkos(SPARTA *sparta) : Domain(sparta)
{

}

/* ---------------------------------------------------------------------- */

DomainKokkos::~DomainKokkos()
{

}

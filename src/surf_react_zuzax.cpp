/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "ctype.h"
#include "stdlib.h"
#include "string.h"
#include "surf_react_zuzax.h"
#include "input.h"
#include "update.h"
#include "comm.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"


enum{DISSOCIATION,EXCHANGE,RECOMBINATION};        // other surf react files
enum{SIMPLE};                                     // other surf react files

#define MAXREACTANT 1
#define MAXPRODUCT 2
#define MAXCOEFF 2

#define MAXLINE 1024
#define DELTALIST 16

#ifdef USE_ZSURF

namespace SPARTA_NS {


//=================================================================================================
SurfReactZuzax::SurfReactZuzax(SPARTA *sparta, int narg, char **arg) : 
    SurfReact(sparta, narg, arg)
{

}

//=================================================================================================

SurfReactZuzax::~SurfReactZuzax()
{
}

/* ---------------------------------------------------------------------- */

void SurfReactZuzax::init()
{
  nsingle = ntotal = 0;
}

/* ---------------------------------------------------------------------- */

void SurfReactZuzax::init_reactions() 
{
}

/* ---------------------------------------------------------------------- */

int SurfReactZuzax::react(Particle::OnePart *&ip, double *tmpp, Particle::OnePart *&jp)
{
  return 0;
}
//=================================================================================================
}

#endif

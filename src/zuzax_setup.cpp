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
#include "zuzax_setup.h"
#include "zuzax/thermo/IdealGasPhase.h"

namespace SPARTA_NS {


//=================================================================================================
ZuzaxSetup::ZuzaxSetup(SPARTA *sparta) : 
    Pointers(sparta) 
{
}
//=================================================================================================
ZuzaxSetup::~ZuzaxSetup()
{
}
//=================================================================================================
void ZuzaxSetup::init()
{
}
//=================================================================================================
void ZuzaxSetup::initGasSetup(int nargs, char** args)
{
 gasThermo =  new Zuzax::IdealGasPhase(args[1], "");
}
//=================================================================================================
}


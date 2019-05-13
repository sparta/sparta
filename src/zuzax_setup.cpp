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
#include <cstring>
#include "zuzax_setup.h"
#include "zuzax/thermo/IdealGasPhase.h"

#include "particle.h"
#include "memory.h"
#include "error.h"
#include <string>

#ifdef USE_ZUZAX

namespace SPARTA_NS {


//=================================================================================================
ZuzaxSetup::ZuzaxSetup(SPARTA *sparta) : 
    Pointers(sparta) ,
    gasThermo(nullptr),
    SptoZu_speciesMap(nullptr)
{
}
//=================================================================================================
ZuzaxSetup::~ZuzaxSetup()
{
   delete (gasThermo);
   memory->destroy(SptoZu_speciesMap);
}
//=================================================================================================
void ZuzaxSetup::init()
{
}
//=================================================================================================
void ZuzaxSetup::initGasSetup(int nargs, char** args)
{
  bool found = false;
  gasThermo =  new Zuzax::IdealGasPhase(args[0], "");
  Particle::Species *species = particle->species;
  memory->create(SptoZu_speciesMap, particle->nspecies, "ZuzaxSetup::SptoZu_speciesMap");

  for (int k = 0; k < particle->nspecies; ++k) {
      Particle::Species& spk = species[k]; 
      found = false;
      std::string sspName(spk.id);
      for (size_t kz = 0; kz < gasThermo->nSpecies(); ++kz) {
          if (sspName == gasThermo->speciesName(kz)) {
              found = true;
              SptoZu_speciesMap[k] = (int) kz;
              break;
          }
      }
      if (!found) {
        char estring[128];
        sprintf(estring, "Can't find a corresponding Zuzax species for the Sparta species, %s",
                        sspName.c_str());
        error->all(FLERR, estring);
      } 
      spk.zuzax_indexGasPhase = SptoZu_speciesMap[k];
  }

   

}
//=================================================================================================
}
#endif


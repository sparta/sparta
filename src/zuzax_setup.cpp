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
#include "zuzax/thermo/SpeciesThermo.h"
#include "zuzax/thermo/StatMech.h"

#include "particle.h"
#include "memory.h"
#include "error.h"
#include <string>

#ifdef USE_ZUZAX

//----------------------------------------------------------------------------------------------------------------------------------
namespace SPARTA_NS {
//==================================================================================================================================
ZuzaxSetup::ZuzaxSetup(SPARTA *sparta) : 
    Pointers(sparta) ,
    gasThermo_(nullptr),
    SptoZu_speciesMap(nullptr)
{
}
//==================================================================================================================================
ZuzaxSetup::~ZuzaxSetup()
{
   delete (gasThermo_);
   memory->destroy(SptoZu_speciesMap);
   memory->destroy(ZutoSp_speciesMap);
}
//=================================================================================================
void ZuzaxSetup::init()
{
}
//=================================================================================================
void ZuzaxSetup::initGasSetup(int nargs, char** args)
{
  int k;
  bool found = false;
  gasThermo_ =  new Zuzax::IdealGasPhase(args[0], "");
  Particle::Species *species = particle->species;
  size_t nsp = gasThermo_->nSpecies();
  memory->create(SptoZu_speciesMap, particle->nspecies, "ZuzaxSetup::SptoZu_speciesMap");
  memory->create(ZutoSp_speciesMap, nsp, "ZuzaxSetup::ZutoSp_speciesMap");
  for (k = 0; k < nsp; ++k) {
    ZutoSp_speciesMap[k] = -1;
  }

  // Loop over all of the gas species that are currently defined within sparta
  for (k = 0; k < particle->nspecies; ++k) {
      SptoZu_speciesMap[k] = -1; 
      Particle::Species& spk = species[k]; 
      std::string sspName(spk.id);
      // Establish the mapping between the Sparta species into the Zuzax ThermoPhase
      // -> The mapping is carried out via a string comparison on the species name
      found = false;
      for (size_t kz = 0; kz < nsp; ++kz) {
          if (sspName == gasThermo_->speciesName(kz)) {
              if (ZutoSp_speciesMap[kz] == -1) { 
                  ZutoSp_speciesMap[kz] = k;
                  found = true;
                  SptoZu_speciesMap[k] = (int) kz;
                  break;
              }
          }
      }
      if (!found) {
        char estring[128];
        sprintf(estring, "Can't find a corresponding Zuzax species for the Sparta species, %s\n",
                        sspName.c_str());
        error->all(FLERR, estring); 
      } 
      // insert the Zuzax species index into the particle class
      spk.zuzax_indexGasPhase = SptoZu_speciesMap[k];

      // Calculate the ezero value with the Sparta's Species struct.
      spk.ezero = calcEzero(spk, spk.zuzax_indexGasPhase);
  }

}
//=================================================================================================
// Calculations follow the notes in SurfaceAblationNotes.docx by hkm
double ZuzaxSetup::calcEzero(Particle::Species& spk, int kgas, int doTDep, double temp_thermal)
{
    double T = 298.15;
    if (doTDep) {
        T = temp_thermal;
    }
    const Zuzax::SpeciesThermoInterpType* stit = gasThermo_->speciesThermo().provideTempDepSTIT(kgas, T);
    const Zuzax::StatMech* sm = dynamic_cast<const Zuzax::StatMech*>(stit);
    const Zuzax::StatMech::speciesStatMechInput& sI = sm->statMechInput();

    if (!sm) {
        Zuzax::writelogf("ZuzaxSetup::calcEzero Warning(): Not StatMech type for species %s, Setting to zero\n",
                         gasThermo_->speciesName(kgas).c_str());
        return 0.0;
    }

    double Hf298 = gasThermo_->Hf298SS(kgas) / Zuzax::Avogadro;
    double H_pv298 = Zuzax::Boltzmann * 298.15;
    double UtranT = 3./2. * Zuzax::Boltzmann * T;
    double Utran298 = 3./2. * Zuzax::Boltzmann * 298.15;
    int nrotdof = spk.rotdof;
    if (nrotdof != 0 && nrotdof != 2 && nrotdof != 3) {
        throw Zuzax::ZuzaxError("ZuzaxSetup::calcEzero", " unknown rotation option\n");
    }
    double UrotT = nrotdof / 2.0 * Zuzax::Boltzmann * T;
    double Urot298 = nrotdof / 2.0 * Zuzax::Boltzmann * 298.15;

    // Note, until Sparta has electron partition function coverage, Uelectron not worth it.
    double UelectronT = 0.0;
    double Uelectron298 = 0.0;
    double theta0 = sI.theta[0];
    double theta0_sp = spk.vibtemp[0];
    if (theta0 != theta0_sp) {
        throw Zuzax::ZuzaxError("ZuzaxSetup::calcEzero", "%s: zuzax and Sparta different theta0: %g %g\n", 
                                gasThermo_->speciesName(kgas).c_str(), theta0 , theta0_sp); 
    }
    double tstar = theta0 / T;
    double UvibT = Zuzax::Boltzmann *  theta0 / 2.0 + Zuzax::Boltzmann *  theta0 / (exp(tstar) + 1.0);
    double tstar298 = theta0 / 298.15;
    double Uvib298 = Zuzax::Boltzmann *  theta0 / 2.0 + Zuzax::Boltzmann *  theta0 / (exp(tstar298) + 1.0);
    double Etran_sp = 3./2. * Zuzax::Boltzmann * T; 
    // Assume smooth rotational modes;
    // get rotdofs
    double nrotdof_sp = spk.rotdof;
    if (nrotdof_sp != nrotdof) {
        throw Zuzax::ZuzaxError("ZuzaxSetup::calcEzero", "zuzax and Sparta different nrotdofs: %d %d\n", 
                                nrotdof, nrotdof_sp);
    }
    double Erot_sp = nrotdof / 2.0 * Zuzax::Boltzmann * T;

    double Evib_sp = Zuzax::Boltzmann * theta0 / (exp(tstar) + 1.0);

    double U_z = Hf298 - H_pv298 + UtranT - Utran298 + UrotT - Urot298
               + UvibT - Uvib298 + UelectronT - Uelectron298;
    double Ezero = U_z - Etran_sp - Erot_sp - Evib_sp;
    return Ezero;
}
//==================================================================================================================================

//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif


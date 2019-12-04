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

#ifndef SPARTA_ZUZAX_SETUP_H
#define SPARTA_ZUZAX_SETUP_H

#include "cstdio"
#include "pointers.h"
#include "particle.h"

#ifdef USE_ZUZAX
#include "zuzax/base/ct_defs.h"
#include "zuzax/thermo/ThermoPhase.h"

//----------------------------------------------------------------------------------------------------------------------------------
namespace SPARTA_NS {

//==================================================================================================================================
//! class to help set up Zuzax
/*!
 *  Note, this class must be set up on all processors. How is this done?
 *
 */
class ZuzaxSetup : protected Pointers {
public:

    //! Constructor 
    /*!
     *  @param(in)           sparta              Pointer to the SPARTA class containing a list
     *                                           of common pointers
     */
    ZuzaxSetup(class SPARTA * sparta);

    //! Virtual destructor
    virtual ~ZuzaxSetup();

    //! initializer
    void init();

    //! Read the Zuzax gas xml file, setting up NASA polynomials 
    /*!
     *  Nasa polynomials allow for a consistent energy level between gas, surf, and solid species.
     *
     *  This will contain one argument, the name of the gas, ThermoPhase file.
     *  There must be zuzax species for every sparta species. Extra zuzax species are treated as
     *  having zero concentrations.
     *  If there isn't a zuzax species for a sparta species, then an error is thrown.
     *
     *  @param[in]           nargs               number of args
     *  @param[in]           args                vector of args
     */
    void initGasSetup(int nargs, char** args);

    //! Calculate the ezero value
    /*!
     *  @param[in]           spk                 Species structure
     *  @param[in]           kgas                index within the gasThermo_ corresponding to the species
     *
     *  @param[in]                               Returns the internal energy offset needed for comparison to NASA polynomials.
     *                                             (units Joules/kmol)
     *  @param[in]           doTDep              calculate the Ezero value as a function of the thermal temperature
     *  @param[in]           temp_thermal        Current Thermal temperature (used if doTDep is true)
     */
    double calcEzero(Particle::Species& spk, int kgas, int doTDep = false, double temp_thermal = 298.15);


 protected:

    //! Pointer to the Gas ThermoPhase object
    Zuzax::thermo_t_double* gasThermo_;

    //! Sparta to Zuzax species map
    /*!
     *  Maps the sparta species index into the Zuzax species index within the gasThermo ThermoPhase routine
     *  Length:   particle->nspecies 
     *  Index:    index of the species within Sparta's list of species
     *  Value:    index of the species within the Zuzax Gas phase ThermoPhase class
     */
    int* SptoZu_speciesMap;

};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif

#endif


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

#ifdef USE_ZSURF
#include "zuzax/base/ct_defs.h"
#include "zuzax/thermo/ThermoPhase.h"

namespace SPARTA_NS {

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
     *  This will contain one argument, the name of the gas ThermoPhase file.
     *  There must be zuzax species for every sparta species. Extra zuzax species are treated as
     *  having zero concentrations.
     *
     *  @param[in]           nargs               number of args
     *  @param[in]           args                vector of args
     */
    void initGasSetup(int nargs, char** args);


 protected:

    Zuzax::thermo_t_double* gasThermo;

};

}
#endif

#endif


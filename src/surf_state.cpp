/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov
   Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "surf_state.h"

#include "memory.h"

// Thing about adding pointers in as a pointer itseslf

//-------------------------------------------------------------------------------------------------
namespace SPARTA_NS
{
/* ---------------------------------------------------------------------- */

SurfState::SurfState()
{
}

/* ---------------------------------------------------------------------- */

void SurfState::init()
{
#ifdef USE_ZSURF

  int nun = net->neq();
  yUnknowns = new double[nun];

  Temp = net->reactor(0).temperature();
  Press = net->reactor(0).pressure();

  int nsp =       net->reservoir(0).thermo(0).nSpecies();
  gasState  = new double[nsp + 2];
  net->reservoir(0).thermo(0).saveState(nsp +2, gasState);

  net->getInitialConditions(0.0, nun, yUnknowns);

  saveState();
  
#endif

}

/* ---------------------------------------------------------------------- */

SurfState::~SurfState()
{
  
#ifdef USE_ZSURF
  delete [] yUnknowns;
  delete [] gasState;
  if (net) {
      delete net;
  }
#endif

}

/* ---------------------------------------------------------------------- */

// do all of the setup for the surface here
void SurfState::setupNewTimeStep(int ntimestep, double dt, double fnum)
{
  // Restore the state within the net object to that of this surface
  setState(ntimestep, dt); 

  // Zero various arrays and counters and update counters that need to be updated at the start of a time step
  net->initTimeStepArrays();

  // A) get gas concentrations, temperature, and pressure

  // m_gas->

  // B) Get surface states


  // C) AdvanceThetaBDF

    // For a time step of size deltaT, this advances the thetas from t_n to t_n+1 using linear stabilization
    // advanceThetaLinearBDF(deltaT);

  // D) Get probability map

     //Create a probability of reaction table for each gas phase species wrt to surface reactions
     /*!
      *  The end result is a fully formed probability map for all gas phase species wrt the current surface, 
      *       m_probMapGasSpecies[k] 
      */
  net->createProbabilityTable();

  // Save probabilty map
  m_probMapGasSpecies = net->m_probMapGasSpecies;


  // E)  Get expected frequency of surface events.  Push that to sparta.

  net->getSurfaceInitiatedEventTables(dt, fnum);

  // Save the surface initiated events table to state_surf

   //! Vector of surface initiated events
    /*!
     *  This is unexpanded list of events where each type of event is listed with a number of times per time step counter.
     *  Thus, its size is small.
     *  This list is created in getSurfaceInitiatedEventTables().
     */
  m_surfInitPSTaskList = net->m_surfInitPSTaskList;


}

/* ---------------------------------------------------------------------- */
// Write the state from the net object to the SurfState object
void SurfState::saveState()
{
#ifdef USE_ZSURF

  int nun = net->neq();

  Temp = net->reactor(0).temperature();
  Press = net->reactor(0).pressure();

  int nsp = net->reservoir(0).thermo(0).nSpecies();
  net->reservoir(0).thermo(0).saveState(nsp +2, gasState);

  double t0 = 0.0;
  net->getInitialConditions(t0, nun, yUnknowns);
  
#endif

}
/* ---------------------------------------------------------------------- */
// Read the state into the net object
void SurfState::setState(int ntimestep, double dt) const
{
#ifdef USE_ZSURF

  int nun = net->neq();

  net->reactor(0).setState_TP(Temp, Press);

  int nsp = net->reservoir(0).thermo(0).nSpecies();
  net->reservoir(0).thermo(0).restoreState(nsp +2, gasState);

  
  double tcurr = ntimestep * dt;
  net->updateState(tcurr,  yUnknowns);
  
#endif

}


}
//-------------------------------------------------------------------------------------------------

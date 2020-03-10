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

#ifndef SPARTA_SURF_STATE_H
#define SPARTA_SURF_STATE_H

#include "stdio.h"
#include "pointers.h"
#ifdef USE_ZSURF
#include "zuzax/zeroD/SurfPropagationSparta.h"
#endif

namespace SPARTA_NS {

// Could add the sparta* pointer  to SurfState 

class SurfState {
  public:

  SurfState(double area);

  virtual ~SurfState();

#ifdef USE_ZSURF
  mutable Zuzax::SurfPropagationSparta* net {nullptr};
#endif   

  virtual void init();

  //! Hook into setting up the Surface calculations for the new time step
  /*!
   *
   */
  virtual void setupNewTimeStep(int ntimestep, double dt, double fnum);

  //! Write the state from the net object to the SurfState object
  /*!
   *  We need to store the state of the surface after every operation, because
   *  the implementation object is reused. So we transfer the state of the surface
   *  from the net object to the SurfState object
   */
  virtual void saveState();

  //! Read in the state of the surface into the net object
  /*!
   *  Transfer the state information from the SurfState object into the net object
   *
   *  @param[in]             ntimestep           time step number
   *  @param[in]             dt                  Delta time step
   */
  virtual void setState(int ntimestep, double dt) const;

  //! Special purpose writing routine
  /*!
   *  @param[in]             ntimestep           time step number
   *  @param[in]             dt                  Delta time step
   */
  void write_step_results(int ntimestep, double dt);

  // ----------------------------------------------- DATA -----------------------------------------
  //! Surface temperature
  double Temp {300.0};

  double Press {1.0E5};

  //! Putting a copy of the area here, because surface evolver ALWAYS needs to know this
  //! to initialize itself.
  double Area {1.0};

  //! Unknowns in the substrate element problem
  double* yUnknowns {nullptr};

  //! State of the gas
  double* gasState {nullptr};

#ifdef USE_ZSURF
  //! Probability map for all the gas phase species
  /*!
   *  This is a map from 0 to one containing the probabilities of what happens when a species collides with this surface
   *
   *    Length:  nsp in the gas phase
   *    index:   gas phase species index
   *    value:   ProbMap struct for each gas phase species. Determines the probability of what happens to the gas species
   *             as it hits the surface.
   */
  std::vector<Zuzax::ProbMap> m_probMapGasSpecies;

  //! Vector of tasks that must be carried out
  /*!
   *  Note that currently this is not used for anything!
   */
  std::vector< struct Zuzax::ZTask > surfInitTaskList;

  //! Vector of surface initiated events to be handled by this surface in this time step
  /*!
   *  This is unexpanded list of events where each type of event is listed with a number of times
   *  per time step counter. Thus, its size is small.
   *  This list is created in getSurfaceInitiatedEventTables().
   *
   *    Length:  Number of event types
   *    index:   event index
   *    value:   ProbMap struct for each gas phase species. Determines the probability of what happens to the gas species
   *             as it hits the surface.
   */
  std::vector< struct Zuzax::PartToSurf > m_surfInitPSTaskList;

#endif
  

};

}

#endif


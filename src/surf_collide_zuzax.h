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

#ifdef SURF_COLLIDE_CLASS

// Turn off recognition of zuzax reactions if this is not defined. This is the
// desired behavior because it will produce an error when reading the input file.
#ifdef USE_ZSURF
// This section is used surf::add_react() to create a 
//    sr[nsr] = new SurfCollideZuzax(sparta,narg,arg);
// statement
SurfCollideStyle(zuzax,SurfCollideZuzax)
#endif

#else

#ifndef SPARTA_SURF_COLLIDE_ZUZAX_H
#define SPARTA_SURF_COLLIDE_ZUZAX_H

#include "surf_collide.h"
#ifdef USE_ZSURF
#include "zuzax/zeroD/SurfPropagationSparta.h"
#include "zuzax/zeroD/ReactorNetDAE.h"
#endif

namespace SPARTA_NS {

#ifdef USE_ZSURF
//! Surface Collision class for zuzax integration
/*!
 *  This class sets up a surface which will be modeled with the continuum code Zuzax.
 *
 *  Mass balances of all substances will be created.
 *
 *  
 *
 */
class SurfCollideZuzax : public SurfCollide {
 public:
  SurfCollideZuzax(class SPARTA *, int, char **);
  SurfCollideZuzax(class SPARTA *sparta) : SurfCollide(sparta) {}
  virtual ~SurfCollideZuzax();
  void init();
  virtual Particle::OnePart *collide(Particle::OnePart *&, double *, double &, int,
                                     SurfState* surfaceState) override;

  virtual void dynamic();

  //! Provide a state object that will be assigned to each surface that will hold the state of 
  //! of the surface
  /*!
   *  (virtual from surf_collide)
   *
   *  @param[in]             area                Area of the surface or fase
   *
   *  @return                                    Returns a pointer to void that will be used
   */
  virtual SurfState* provideStateObject(double area) const override;


  //! Initialize the Network model with all of the ThermoPhase classes
  void initNetwork();

  //! Setup a new time step
  /*!
   *  Zuzax uses this hook to calculate the probability table for reactions on a surface
   *  and the event table for surface initiated events that will take place during
   *  the next time step.
   */
  virtual void setupNewTimeStep() override;

  void rollDiceOnParticleSurfInteraction(Particle::OnePart *&ip, SurfState* surfState,
                                         int& irxn, int& idir);


 protected:
  double twall;              // surface temperature (Kelvin)
  double pwall;              // surface pressure ( Pascal )
  double acc;                // surface accomodation coeff
  double vx,vy,vz;           // translational velocity of surface
  double wx,wy,wz;           // angular velocity of surface
  double px,py,pz;           // point to rotate surface around
  int tflag,rflag;           // flags for translation and rotation
  int trflag;                // 1 if either tflag or rflag is set

  char *tstr {nullptr};      // temperature variable name (NULL if constant)
  int tvar;                  // index of equal-style variable
  int m_initNetwork {0};     // Flag indicating network has been initialized
  char *inputConfigFile {nullptr};       // input file

  double vstream[3];
  class RanPark *random {nullptr};     // RNG for particle reflection

  void diffuse(Particle::OnePart *, double *);

  //! Pointer to a malloced net that will be used as a base in a copy constructor
  //! to malloc state objects for all surfaces
  Zuzax::SurfPropagationSparta* baseNet {nullptr};

  //! If this is true, the surface site fractions are adjusted to a pseudo-steady
  //! state condition at the start of the calculation.
  int initAsPseudoSteadyState {0};

};

#endif
}

#endif
#endif


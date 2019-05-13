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

#ifndef SPARTA_SURF_REACT_ZUZAX_H
#define SPARTA_SURF_REACT_ZUZAX_H

#include "pointers.h"
#include "particle.h"
#include "surf_react.h"
#ifdef USE_ZSURF
#include "zuzax/base/ct_defs.h"
#endif

namespace SPARTA_NS {

#ifdef USE_ZSURF
//! Class which carries out complex surface reactions with etching or growth
/*!
 *
 */
class SurfReactZuzax : public SurfReact
{
public:
    //! Constructor
    SurfReactZuzax(class SPARTA *, int, char **);

    //! Destructor
    virtual ~SurfReactZuzax();

    //! Initializer
    virtual void init() override;

    void init_reactions();

    //! Main routine which handles collisions
    virtual int react(Particle::OnePart *&, double *, Particle::OnePart *&)  override;

 private:

     class RanPark *random;     // RNG for reaction probabilities
};

#endif

}
#endif

/* ERROR/WARNING messages:

E: Surf_react ID must be alphanumeric or underscore characters

Self-explanatory.

*/

/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS

CommandStyle(scale_particles,ScaleParticles)

#else

#ifndef SPARTA_SCALE_PARTICLES_H
#define SPARTA_SCALE_PARTICLES_H

#include "pointers.h"

namespace SPARTA_NS {

class ScaleParticles : protected Pointers {

 public:
  ScaleParticles(class SPARTA *);
  void command(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/

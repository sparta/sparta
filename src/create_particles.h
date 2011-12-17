/* ----------------------------------------------------------------------
   DSMC - Sandia parallel DSMC code
   www.sandia.gov/~sjplimp/dsmc.html
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2011) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level DSMC directory.
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS

CommandStyle(create_particles,CreateParticles)

#else

#ifndef DSMC_CREATE_PARTICLES_H
#define DSMC_CREATE_PARTICLES_H

#include "pointers.h"

namespace DSMC_NS {

class CreateParticles : protected Pointers {
 public:
  CreateParticles(class DSMC *);
  void command(int, char **);

 private:
  int seed;

  void create_all(bigint);
  void create_local(bigint);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot create_particles before simulation box is defined

UNDOCUMENTED

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running DSMC to see the offending line.

E: Created incorrect # of particles = %ld

UNDOCUMENTED

*/

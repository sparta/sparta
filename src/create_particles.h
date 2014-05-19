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

#ifdef COMMAND_CLASS

CommandStyle(create_particles,CreateParticles)

#else

#ifndef SPARTA_CREATE_PARTICLES_H
#define SPARTA_CREATE_PARTICLES_H

#include "pointers.h"

namespace SPARTA_NS {

class CreateParticles : protected Pointers {

 public:
  CreateParticles(class SPARTA *);
  void command(int, char **);
  int evib(int);
  double erot(int);

 private:
  int imix,single,mspecies;
  double xp,yp,zp,vx,vy,vz;

  void create_single();
  void create_local(bigint);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot create particles before simulation box is defined

UNDOCUMENTED

E: Cannot create particles  before grid is defined

UNDOCUMENTED

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Create_particles mixture ID does not exist

UNDOCUMENTED

E: Create_particles species ID does not exist

UNDOCUMENTED

E: Cannot use n and single in create_particles command

UNDOCUMENTED

E: Created incorrect # of particles: %ld versus %ld

UNDOCUMENTED

E: Create_particles single requires z = 0 for 2d simulation

UNDOCUMENTED

E: Could not create a single particle

UNDOCUMENTED

*/

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

CommandStyle(velocity,Velocity)

#else

#ifndef DSMC_VELOCITY_H
#define DSMC_VELOCITY_H

#include "pointers.h"

namespace DSMC_NS {

class Velocity : protected Pointers {
 public:
  Velocity(class DSMC *);
  void command(int, char **);

 private:
  double t_desired,vx,vy,vz;
  int seed;
  int sum_flag;

  void create();
  void set();
  void options(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running DSMC to see the offending line.

E: Velocity command before simulation box is defined

The velocity command cannot be used before a read_data, read_restart,
or create_box command.

E: Velocity command with no particles existing

UNDOCUMENTED

E: Cannot set velocity to non-zero z value for 2d simulation

UNDOCUMENTED

*/
